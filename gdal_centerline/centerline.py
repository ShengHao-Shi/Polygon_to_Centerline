# -*- coding: utf-8 -*-
"""
centerline.py
=============
Open-source Polygon-to-Centerline algorithm.

No ArcPy required.  Dependencies (all available in Anaconda / conda-forge):
    numpy, scipy, shapely, geopandas, scikit-image, networkx, matplotlib

Algorithm overview
------------------
Two complementary strategies are implemented and can be selected via the
`method` parameter of :func:`polygon_to_centerline`:

Method A – ``"voronoi"`` (default)
    A vector-based Voronoi / Thiessen skeleton approach:
    1. Densify the polygon boundary so that vertices are spaced at most
       *densify_distance* apart.
    2. Compute the Voronoi tessellation of those boundary points using
       ``scipy.spatial.Voronoi``.
    3. Retain only those Voronoi ridge segments whose **both endpoints** lie
       **inside** the original polygon (using ``shapely.contains``).
    4. Build a graph from the retained segments (NetworkX).
    5. Iteratively prune dead-end branches (degree-1 nodes) until only the
       main skeleton trunk remains.
    6. Extract the single longest path through the graph (two-pass Dijkstra
       diameter algorithm) and return it as a non-branching LineString.

Method B – ``"skeleton"``
    A raster-based morphological skeletonisation approach:
    1. Rasterise the polygon to a binary image at a configurable resolution.
    2. Apply ``skimage.morphology.skeletonize`` (Zhang–Suen thinning).
    3. Build a weighted pixel graph (8-connected, edge weights = Euclidean
       distance between pixel centres).
    4. Extract the single longest path using the two-pass Dijkstra diameter
       algorithm and return it as a non-branching LineString.

Both methods accept a ``single_line=True`` parameter (the default).  Set it
to ``False`` to recover the original full-skeleton / multi-branch behaviour.

Usage
-----
    from centerline import polygon_to_centerline
    import geopandas as gpd

    gdf = gpd.read_file("my_polygons.geojson")
    centerlines = polygon_to_centerline(gdf, method="voronoi")
    centerlines.to_file("my_centerlines.geojson", driver="GeoJSON")
"""

from __future__ import annotations

import math
import warnings
from typing import Literal, Sequence

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
)
from shapely.ops import linemerge, unary_union

# scikit-image import – optional for method="skeleton"
try:
    from skimage.morphology import skeletonize as _skeletonize

    _SKIMAGE_AVAILABLE = True
except ImportError:  # pragma: no cover
    _SKIMAGE_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def polygon_to_centerline(
    source: gpd.GeoDataFrame | gpd.GeoSeries,
    method: Literal["voronoi", "skeleton"] = "voronoi",
    densify_distance: float = 1.0,
    prune_threshold: float = 0.0,
    smooth_sigma: float = 0.0,
    raster_resolution: float | None = None,
    single_line: bool = True,
) -> gpd.GeoDataFrame:
    """
    Convert polygon features to centerline polylines.

    Parameters
    ----------
    source : GeoDataFrame or GeoSeries
        Input polygon layer.  Any CRS is accepted; coordinates are used as-is.
    method : {"voronoi", "skeleton"}
        Algorithm to use.  "voronoi" (default) is vector-based and generally
        produces cleaner results.  "skeleton" is raster-based and works well
        for very complex or jagged shapes.
    densify_distance : float
        Maximum spacing (in CRS units) between vertices inserted along the
        polygon boundary before Voronoi tessellation.  Smaller = finer
        skeleton (method="voronoi" only).
    prune_threshold : float
        Remove skeleton branches shorter than this value (in CRS units).
        0 = no pruning (default).  Only used by method="voronoi".
    smooth_sigma : float
        Standard-deviation (in CRS units) of a Gaussian blur applied to the
        raster before skeletonisation.  0 = no blur (default).
        Only used by method="skeleton".
    raster_resolution : float or None
        Cell size for rasterisation (method="skeleton").  If None, it is set
        automatically to densify_distance.
    single_line : bool
        When True (default), the output is always a single, non-branching
        LineString — the longest trunk path through the skeleton graph,
        found by the two-pass Dijkstra diameter algorithm.  Set to False to
        return the full skeleton as a MultiLineString (may contain forks).

    Returns
    -------
    GeoDataFrame
        Polyline GeoDataFrame in the same CRS as the input, with one row per
        input polygon.  Non-polygon geometries are silently skipped.
    """
    if isinstance(source, gpd.GeoSeries):
        gdf = gpd.GeoDataFrame(geometry=source)
    else:
        gdf = source.copy()

    results = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            results.append(None)
            continue
        # Explode multi-polygons so each part is processed independently
        if isinstance(geom, MultiPolygon):
            parts = list(geom.geoms)
        elif isinstance(geom, Polygon):
            parts = [geom]
        else:
            results.append(None)
            continue

        part_lines = []
        for poly in parts:
            if poly.is_empty or not poly.is_valid:
                poly = poly.buffer(0)
            try:
                if method == "voronoi":
                    cl = _centerline_voronoi(poly, densify_distance, prune_threshold, single_line)
                elif method == "skeleton":
                    res = raster_resolution if raster_resolution else densify_distance
                    cl = _centerline_skeleton(poly, res, smooth_sigma, single_line)
                else:
                    raise ValueError(f"Unknown method: {method!r}")
            except (RuntimeError, ValueError) as exc:
                warnings.warn(f"Skipping polygon: {exc}", stacklevel=2)
                cl = None
            if cl is not None and not cl.is_empty:
                part_lines.append(cl)

        if part_lines:
            merged = linemerge(unary_union(part_lines)) if len(part_lines) > 1 else part_lines[0]
            results.append(merged)
        else:
            results.append(None)

    out_gdf = gdf.copy()
    out_gdf["geometry"] = results
    out_gdf = out_gdf[out_gdf.geometry.notna() & ~out_gdf.geometry.is_empty]
    out_gdf = out_gdf.set_geometry("geometry")
    return out_gdf


# ---------------------------------------------------------------------------
# Method A – Voronoi skeleton
# ---------------------------------------------------------------------------


def _densify(polygon: Polygon, max_distance: float) -> np.ndarray:
    """
    Return an (N, 2) array of the densified exterior-ring coordinates.
    Interior rings (holes) are also included so the Voronoi diagram
    correctly represents polygons with holes.
    """
    rings = [polygon.exterior] + list(polygon.interiors)
    points = []
    for ring in rings:
        coords = np.array(ring.coords)
        for i in range(len(coords) - 1):
            p0, p1 = coords[i], coords[i + 1]
            seg_len = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
            n_insert = max(1, math.ceil(seg_len / max_distance))
            for k in range(n_insert):
                t = k / n_insert
                points.append(p0 + t * (p1 - p0))
    return np.array(points)


def _centerline_voronoi(
    polygon: Polygon,
    densify_distance: float,
    prune_threshold: float,
    single_line: bool = True,
) -> LineString | MultiLineString | None:
    """
    Compute the centerline of *polygon* using a Voronoi skeleton.

    Steps
    -----
    1. Densify the polygon boundary.
    2. Compute the Voronoi diagram of the boundary points.
    3. Keep only ridge segments fully inside the polygon, filtering out
       same-side artefacts via a generator-distance threshold.
    4. Build a graph; iteratively prune degree-1 branches shorter than
       *prune_threshold*.
    5. If *single_line* is True (default), extract the single longest trunk
       path using the two-pass Dijkstra diameter algorithm and return it as
       a LineString.  Otherwise return all surviving edges (may contain forks).
    """
    pts = _densify(polygon, densify_distance)
    if len(pts) < 4:
        return None

    # Voronoi diagram
    vor = Voronoi(pts)

    # Minimum distance between the two generating boundary points.
    # Same-side adjacent points are ~densify_distance apart and produce
    # parallel-to-boundary artefact ridges (the "fan" spikes visible at
    # tight bends).  True centerline ridges bridge OPPOSITE sides of the
    # polygon and have generators separated by at least ~1 polygon-width.
    # Requiring gen_dist > 3 × densify_distance reliably removes same-side
    # artefacts while preserving all cross-polygon (centerline) ridges for
    # polygons wider than ~3 × densify_distance (the typical use case).
    min_gen_dist = 3.0 * densify_distance

    # Collect ridge segments that are fully contained inside the polygon
    segments = []
    for ridge_idx, ridge in enumerate(vor.ridge_vertices):
        i, j = ridge
        if i < 0 or j < 0:
            # Infinite ridge — skip
            continue

        # --- Generator-distance filter (same-side artefact removal) ---
        p_idx, q_idx = vor.ridge_points[ridge_idx]
        gp0 = pts[p_idx]
        gp1 = pts[q_idx]
        gen_dist = math.hypot(gp1[0] - gp0[0], gp1[1] - gp0[1])
        if gen_dist < min_gen_dist:
            continue

        v0 = vor.vertices[i]
        v1 = vor.vertices[j]
        seg = LineString([v0, v1])

        # Full-segment containment: the entire ridge must lie inside the
        # polygon, not just the midpoint.  This removes ridges that merely
        # have a midpoint inside a concave "bay" of the polygon.
        if polygon.contains(seg):
            segments.append((tuple(v0), tuple(v1)))

    if not segments:
        return None

    # Build NetworkX graph
    G = nx.Graph()
    for p, q in segments:
        length = math.hypot(q[0] - p[0], q[1] - p[1])
        G.add_edge(p, q, weight=length)

    # Prune dead-end branches shorter than prune_threshold
    if prune_threshold > 0:
        G = _prune_branches(G, prune_threshold)

    if G.number_of_edges() == 0:
        return None

    if single_line:
        # Extract the single longest trunk (no forks / branches)
        path_nodes = _extract_longest_path(G)
        if len(path_nodes) < 2:
            return None
        return LineString(list(path_nodes))

    # single_line=False: return the full skeleton (may contain forks)
    lines = [LineString([u, v]) for u, v in G.edges()]
    return linemerge(lines) if lines else None


def _prune_branches(G: nx.Graph, threshold: float) -> nx.Graph:
    """
    Iteratively remove skeleton branches (dead-end paths) whose **total
    path length** is shorter than *threshold*.

    A branch is defined as a chain of nodes that starts at a leaf (degree 1)
    and runs until either:
      - the first junction node (degree ≥ 3), or
      - another leaf (isolated edge between two leaves).

    The entire branch is removed only when its accumulated length < threshold.
    Junction nodes are kept to preserve connectivity of the remaining skeleton.
    """
    changed = True
    while changed:
        changed = False
        leaves = [n for n in list(G.nodes()) if G.degree(n) == 1]
        for leaf in leaves:
            if leaf not in G or G.degree(leaf) != 1:
                continue
            # Walk the branch from the leaf toward the first junction
            path = [leaf]
            total_len = 0.0
            current = leaf
            while True:
                neighbors = list(G.neighbors(current))
                # Pick the next unvisited node; stop if all neighbours are
                # already in the path (shouldn't happen for a leaf-started
                # walk, but acts as a safe guard).
                nxt = None
                for nb in neighbors:
                    if nb not in path:
                        nxt = nb
                        break
                if nxt is None:
                    break
                edge_len = G[current][nxt].get("weight", 0.0)
                total_len += edge_len
                path.append(nxt)
                current = nxt
                # Stop walking if we've reached a junction or another leaf
                if G.degree(current) != 2:
                    break
                # Stop early if already longer than threshold (optimisation)
                if total_len >= threshold:
                    break

            if total_len < threshold:
                # Remove all branch nodes except the terminal junction/leaf
                nodes_to_remove = path[:-1]  # keep the junction end
                if len(path) >= 2 and G.degree(path[-1]) == 1:
                    # Isolated edge: both endpoints are leaves — remove all
                    nodes_to_remove = path
                for node in nodes_to_remove:
                    if node in G:
                        G.remove_node(node)
                changed = True
    return G


def _extract_longest_path(G: nx.Graph) -> list:
    """
    Find the single longest (by sum of edge weights) path between two leaf
    nodes (degree = 1) in *G*, using the two-pass Dijkstra diameter algorithm.

    This is the standard O(V + E log V) tree-diameter algorithm:
      Pass 1 – from an arbitrary leaf, find the farthest leaf *u*.
      Pass 2 – from *u*, find the farthest leaf *v*.
    The path u → v is the weighted diameter of the tree and corresponds to
    the main trunk of the skeleton (no branches / forks).

    For graphs with cycles (rare in Voronoi skeletons) the algorithm still
    returns the pair of leaves with the maximum shortest-path distance between
    them, which is a robust proxy for the longest trunk.

    Returns
    -------
    list
        Ordered list of graph nodes from one end of the trunk to the other.
    """
    if G.number_of_nodes() < 2:
        return list(G.nodes())

    leaves = [n for n in G.nodes() if G.degree(n) == 1]
    if not leaves:
        # No leaves (all junctions) — pick any two nodes as endpoints
        leaves = list(G.nodes())

    # Pass 1: distances from an arbitrary starting leaf
    start = leaves[0]
    dist1, path1 = nx.single_source_dijkstra(G, start, weight="weight")
    # Farthest leaf from 'start'
    u = max(leaves, key=lambda n: dist1.get(n, 0.0))

    # Pass 2: distances from u — finds the true diameter endpoint
    dist2, path2 = nx.single_source_dijkstra(G, u, weight="weight")
    v = max(leaves, key=lambda n: dist2.get(n, 0.0))

    return path2.get(v, [u, v])


# ---------------------------------------------------------------------------
# Method B – Raster skeleton
# ---------------------------------------------------------------------------


def _centerline_skeleton(
    polygon: Polygon,
    resolution: float,
    smooth_sigma: float = 0.0,
    single_line: bool = True,
) -> LineString | MultiLineString | None:
    """
    Compute the centerline of *polygon* using raster skeletonisation.

    Steps
    -----
    1. Rasterise the polygon to a binary image at *resolution* cell size.
    2. Optionally apply Gaussian blurring before thinning.
    3. Apply Zhang–Suen thinning (scikit-image ``skeletonize``).
    4. Build a NetworkX graph of the skeleton pixels (8-connected neighbours,
       edge weights = Euclidean pixel distance).
    5. If *single_line* is True (default), extract the longest trunk path
       using ``_extract_longest_path`` and return it as a single LineString.
       Otherwise return all pixel-pair edges as a MultiLineString.
    """
    if not _SKIMAGE_AVAILABLE:
        raise ImportError(
            "scikit-image is required for method='skeleton'. "
            "Install it with: pip install scikit-image"
        )

    minx, miny, maxx, maxy = polygon.bounds

    # Add a 1-cell border so the skeleton never touches the image edge
    border = resolution
    minx -= border
    miny -= border
    maxx += border
    maxy += border

    cols = max(2, math.ceil((maxx - minx) / resolution))
    rows = max(2, math.ceil((maxy - miny) / resolution))

    # Rasterise the polygon to a binary numpy array
    binary = _rasterize_polygon(polygon, minx, miny, resolution, rows, cols)

    # Optional Gaussian smoothing (helps with jagged boundaries)
    if smooth_sigma > 0:
        from scipy.ndimage import gaussian_filter

        sigma_px = smooth_sigma / resolution
        blurred = gaussian_filter(binary.astype(float), sigma=sigma_px)
        binary = blurred > 0.5

    # Morphological skeletonisation (Zhang–Suen thinning)
    skeleton = _skeletonize(binary)

    # Build a weighted graph of skeleton pixels
    G = _build_skeleton_graph(skeleton, minx, miny, resolution)

    if G.number_of_edges() == 0:
        return None

    if single_line:
        # Extract the single longest trunk (no forks / branches)
        path_pixels = _extract_longest_path(G)
        if len(path_pixels) < 2:
            return None
        coords = [G.nodes[p]["coord"] for p in path_pixels]
        return LineString(coords)

    # single_line=False: return all skeleton edges as a MultiLineString
    lines = _graph_to_lines(G)
    return linemerge(lines) if lines else None


def _rasterize_polygon(
    polygon: Polygon,
    minx: float,
    miny: float,
    resolution: float,
    rows: int,
    cols: int,
) -> np.ndarray:
    """
    Rasterise *polygon* into a boolean numpy array of shape (rows, cols).
    Pixels whose **centre** falls inside the polygon are True.
    """
    binary = np.zeros((rows, cols), dtype=bool)
    # Use shapely's vectorised contains check via a grid of centroids
    xs = minx + (np.arange(cols) + 0.5) * resolution
    ys = miny + (np.arange(rows) + 0.5) * resolution
    # Build meshgrid and check containment row by row to limit memory usage
    for r in range(rows):
        y = ys[r]
        # Quick bounding-box check: skip rows entirely outside polygon
        if y < polygon.bounds[1] or y > polygon.bounds[3]:
            continue
        row_points = [Point(x, y) for x in xs]
        for c, pt in enumerate(row_points):
            if polygon.contains(pt):
                binary[r, c] = True
    return binary


def _build_skeleton_graph(
    skeleton: np.ndarray,
    minx: float,
    miny: float,
    resolution: float,
) -> nx.Graph:
    """
    Build a weighted NetworkX graph from a binary skeleton image.

    Nodes are (row, col) pixel indices with a ``coord`` attribute giving the
    world-coordinate centre of each pixel.  Edge weights are the Euclidean
    pixel distances (1 for cardinal neighbours, √2 for diagonal neighbours).
    """
    # Local helper: pixel (r, c) → world-coordinate centre.
    # Defined as a closure over minx/miny/resolution to avoid passing
    # three extra arguments through every call site.
    def px_to_coord(r: int, c: int) -> tuple[float, float]:
        return (minx + (c + 0.5) * resolution, miny + (r + 0.5) * resolution)

    yx = np.argwhere(skeleton)
    if len(yx) == 0:
        return nx.Graph()

    pixel_set = set(map(tuple, yx))
    G = nx.Graph()
    for r, c in yx:
        G.add_node((r, c), coord=px_to_coord(r, c))
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                nb = (r + dr, c + dc)
                if nb in pixel_set:
                    # Weight = Euclidean distance between pixel centres
                    G.add_edge((r, c), nb, weight=math.hypot(dr, dc) * resolution)

    return G


def _graph_to_lines(G: nx.Graph) -> list[LineString]:
    """Convert graph edges to a list of two-point LineString segments."""
    lines = []
    for u, v in G.edges():
        cu = G.nodes[u]["coord"]
        cv = G.nodes[v]["coord"]
        lines.append(LineString([cu, cv]))
    return lines


# ---------------------------------------------------------------------------
# File I/O helpers
# ---------------------------------------------------------------------------


def read_polygons(path: str) -> gpd.GeoDataFrame:
    """Read a polygon vector file (GeoJSON, Shapefile, GPKG, …)."""
    return gpd.read_file(path)


def write_centerlines(gdf: gpd.GeoDataFrame, path: str, driver: str | None = None) -> None:
    """Write the centerline GeoDataFrame to a file."""
    if driver is None:
        ext = path.rsplit(".", 1)[-1].lower()
        driver_map = {"geojson": "GeoJSON", "shp": "ESRI Shapefile", "gpkg": "GPKG"}
        driver = driver_map.get(ext, "GeoJSON")
    gdf.to_file(path, driver=driver)
