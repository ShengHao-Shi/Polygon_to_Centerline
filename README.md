# Polygon to Centerline

An open-source and ArcGIS toolbox repository that converts elongated polygon
features into polyline centerlines.

---

## Repository Structure

```
Polygon_to_Centerline/
├── arcpy_toolbox/          ArcGIS Python Toolbox (.pyt) – requires ArcGIS Pro / ArcMap
│   ├── Polygon_to_Centerline.pyt
│   └── README.md
├── gdal_centerline/        Open-source implementation – no ArcPy required
│   ├── centerline.py       Core algorithm (importable module)
│   ├── cli.py              Command-line interface
│   ├── requirements.txt    pip / conda dependencies
│   └── README.md           Full algorithm explanation
└── README.md               This file
```

---

## Two Implementations

### 1. `arcpy_toolbox/` — ArcGIS Python Toolbox

Uses ArcPy (bundled with ArcGIS Pro / ArcMap) and implements the
Voronoi / Thiessen-polygon skeleton approach via built-in ArcGIS tools.

| Requirement | |
|---|---|
| ArcGIS Pro 2.x+ or ArcMap 10.x | **Standard / Advanced** license |
| arcpy | Bundled with ArcGIS |

See [`arcpy_toolbox/README.md`](arcpy_toolbox/README.md) for installation
and usage.

---

### 2. `gdal_centerline/` — Open-Source (no ArcPy)

Uses only standard scientific Python libraries available in **Anaconda /
conda-forge**.  Two algorithms are provided:

| Method | Libraries | Best for |
|---|---|---|
| `"voronoi"` (default) | `scipy`, `shapely`, `networkx` | Road / river / building polygons |
| `"skeleton"` | `scikit-image`, `shapely`, `networkx` | Organic / complex shapes |

**Quick start:**

```bash
# Install dependencies
conda install -c conda-forge numpy scipy shapely geopandas scikit-image networkx matplotlib
# or
pip install -r gdal_centerline/requirements.txt

# Run on your data
python gdal_centerline/cli.py input.geojson output.geojson

# With options
python gdal_centerline/cli.py rivers.gpkg centerlines.gpkg \
    --method voronoi --densify 0.5 --prune 5.0
```

**Use as a library:**

```python
import geopandas as gpd
from gdal_centerline.centerline import polygon_to_centerline

gdf = gpd.read_file("my_polygons.geojson")
centerlines = polygon_to_centerline(gdf, method="voronoi", densify_distance=1.0)
centerlines.to_file("my_centerlines.geojson", driver="GeoJSON")
```

See [`gdal_centerline/README.md`](gdal_centerline/README.md) for the full
algorithm description and parameter reference.

---

## Algorithm Overview (both implementations)

Both implementations are based on the **Voronoi / Thiessen skeleton**:

1. **Densify** the polygon boundary (insert extra vertices every N metres).
2. **Voronoi tessellation** — partition space so every point belongs to
   its nearest boundary vertex.
3. **Filter interior edges** — keep only Voronoi ridges that lie fully
   inside the polygon and whose two generating vertices come from
   *opposite* sides of the polygon (not adjacent same-side vertices).
4. **Prune dead-end branches** shorter than a configurable threshold.
5. **Output** the surviving edges as a polyline feature class / GeoJSON.

The open-source implementation additionally offers a **raster skeleton**
method (morphological thinning via `scikit-image`).

---

## Example Output

![Centerline demo](https://github.com/user-attachments/assets/b360be4f-c1d0-4add-b5a0-dd552580c379)

*Left: Voronoi (no pruning) · Centre: Voronoi (pruned) · Right: Raster skeleton*
*Top row: simple rectangle · Bottom row: sinuous strip polygon*

