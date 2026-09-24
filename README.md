# cross_section_plugin
Cross section plugin for QGIS, generates cross-sections from a line, then extracts elevations from TIFs, saves them and compare them.
=======
# Cross-section profiles from QGIS rasters

Generate river or coastal cross-sections in QGIS, sample one or two elevation GeoTIFFs, and plot each profile with comparison statistics.

## Requirements

- **Generation:** QGIS Python environment with GDAL (`osgeo`). `cross_section_tool.py` contains the `CrossSectionPlugin` class for integration into an existing QGIS plugin; it is not a standalone command-line program or a complete installable plugin package.
- **Plotting:** Python 3 and Matplotlib (`pip install matplotlib`). QGIS is not needed for plotting.

## Workflow

1. Load a line or multipart line and up to two local `.tif` / `.tiff` rasters in QGIS. Open the cross-section tool.
2. Select a **projected CRS in metres** (EPSG:32633 is the initial choice for geographic input). Choose the interval between sections, distance **on each side** of the line, raster sampling spacing, and nearest or bilinear interpolation.
3. Select an existing output directory. With at least one raster selected, the tool writes one `cs00001_<raster>.csv` per cross-section. It also adds section lines to QGIS. Optionally save section lines and sampled points, including raster values, to a GeoPackage. Without rasters it can generate the section lines alone.
4. Put `plot_cross_sections.py` next to the `cs*.csv` files (or edit `INPUT_DIR` near the top of that script) and run:

   ```bash
   python3 plot_cross_sections.py
   ```

The plotting script creates `plot_cross-sections/` with one PNG per section and `cross_section_stats.csv`. PNGs show both raster profiles when present, section chainage, centre coordinates, average depths, and paired comparison statistics.

## Data and statistics

Profile CSV columns are `distance_m,x,y,z_tif1[,z_tif2]`, with `#` comments identifying the input line, rasters, working CRS, source feature, part, and chainage. The row at `distance_m=0` gives the section centre. Negative distances are to the **right** and positive distances to the **left** of the directed input line. Missing raster values are `nan`.

The statistics CSV includes `center_x`, `center_y`, `working_crs`, source feature ID, part number, chainage, valid-sample counts, mean raster elevation (`mean_z_*`), mean water depth (`mean_depth_*`), and, where both rasters have data at the same positions, their paired mean differences and elevation RMSE. Use the recorded CRS when mapping the centre coordinates.

By default, depth is measured below `DATUM = 0` from **negative** raster elevations; positions above the datum count as zero depth. Edit `DATUM` and `DEPTH_IS_NEGATIVE_Z` in `plot_cross_sections.py` to match your raster convention. The paired depth difference is TIF2 minus TIF1: a positive value means TIF2 is deeper at the shared valid samples.
