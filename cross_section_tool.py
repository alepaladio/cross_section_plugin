"""QGIS cross sections and bathymetry/topography profiles.

The working CRS must be projected in metres. Signed profile distance is
negative on the right and positive on the left of the directed input line.
"""

import csv
import math
import os
import re
import tempfile

from osgeo import gdal
from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import (
    QAction, QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QDoubleSpinBox, QVBoxLayout,
)
from qgis.core import (
    QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsFeature,
    QgsField, QgsGeometry, QgsMapLayerProxyModel, QgsPointXY, QgsProject,
    QgsUnitTypes, QgsVectorFileWriter, QgsVectorLayer, QgsWkbTypes,
    QgsMessageLog,
)
from qgis.gui import QgsMapLayerComboBox, QgsProjectionSelectionWidget


def safe_stem(name):
    """Create a portable filename component."""
    stem = os.path.splitext(os.path.basename(name))[0]
    return re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_") or "raster"


def sample_distances(half_length, spacing):
    """Include the centre and both endpoints, even with uneven spacing."""
    positions = [0.0, -half_length, half_length]
    for i in range(1, int(math.floor(half_length / spacing)) + 1):
        value = i * spacing
        if value < half_length - 1e-9:
            positions.extend((-value, value))
    return sorted(set(positions))


def section_chainages(length, interval):
    positions = [i * interval for i in range(int(math.floor(length / interval)) + 1)]
    if length - positions[-1] > 1e-8:
        positions.append(length)
    return positions


def tangent_at(line, chainage):
    """At an interior vertex choose the forward nonzero segment."""
    cumulative, last = 0.0, None
    for start, end in zip(line, line[1:]):
        dx, dy = end.x() - start.x(), end.y() - start.y()
        segment_length = math.hypot(dx, dy)
        if segment_length <= 1e-12:
            continue
        last = (dx / segment_length, dy / segment_length)
        if chainage < cumulative + segment_length - 1e-8:
            return last
        cumulative += segment_length
    return last


def point_on_line(line, chainage):
    cumulative = 0.0
    for start, end in zip(line, line[1:]):
        segment_length = math.hypot(end.x() - start.x(), end.y() - start.y())
        if segment_length <= 1e-12:
            continue
        if chainage <= cumulative + segment_length:
            fraction = min(1.0, max(0.0, (chainage - cumulative) / segment_length))
            return QgsPointXY(start.x() + fraction * (end.x() - start.x()),
                              start.y() + fraction * (end.y() - start.y()))
        cumulative += segment_length
    return QgsPointXY(line[-1])


class RasterSampler:
    """Sample band one in the working CRS using nearest or bilinear."""

    def __init__(self, path, target_crs, interpolation, temporary_dir, index):
        dataset = gdal.Open(path, gdal.GA_ReadOnly)
        if dataset is None or dataset.RasterCount < 1:
            raise ValueError("Cannot read raster: " + path)
        source_crs = QgsCoordinateReferenceSystem.fromWkt(dataset.GetProjection())
        if not source_crs.isValid():
            raise ValueError("Raster has no valid CRS: " + path)
        if source_crs != target_crs:
            warped_path = os.path.join(temporary_dir, "raster_{}.tif".format(index))
            warped = gdal.Warp(
                warped_path, dataset, format="GTiff", dstSRS=target_crs.toWkt(),
                resampleAlg="bilinear" if interpolation == "Bilinear" else "near",
                outputType=gdal.GDT_Float32, multithread=True,
            )
            if warped is None:
                raise ValueError("Could not reproject raster: " + path)
            warped = None
            dataset = None
            dataset = gdal.Open(warped_path, gdal.GA_ReadOnly)
            if dataset is None:
                raise ValueError("Could not reopen reprojected raster: " + path)
        self.dataset = dataset
        self.band = dataset.GetRasterBand(1)
        self.nodata = self.band.GetNoDataValue()
        self.inverse = gdal.InvGeoTransform(dataset.GetGeoTransform())
        if self.inverse is None:
            raise ValueError("Raster has an invalid geotransform: " + path)
        self.interpolation = interpolation

    def _value(self, col, row):
        if not (0 <= col < self.dataset.RasterXSize and 0 <= row < self.dataset.RasterYSize):
            return math.nan
        data = self.band.ReadAsArray(int(col), int(row), 1, 1)
        if data is None:
            return math.nan
        value = float(data[0, 0])
        if not math.isfinite(value) or (self.nodata is not None and value == self.nodata):
            return math.nan
        return value

    def sample(self, x, y):
        pixel_x = self.inverse[0] + self.inverse[1] * x + self.inverse[2] * y
        pixel_y = self.inverse[3] + self.inverse[4] * x + self.inverse[5] * y
        if not (0 <= pixel_x < self.dataset.RasterXSize and
                0 <= pixel_y < self.dataset.RasterYSize):
            return math.nan
        if self.interpolation == "Nearest":
            return self._value(math.floor(pixel_x), math.floor(pixel_y))

        # Pixel coordinates from GDAL refer to corners. Interpolate at centres.
        cx, cy = pixel_x - 0.5, pixel_y - 0.5
        c0, r0 = math.floor(cx), math.floor(cy)
        fx, fy = cx - c0, cy - r0
        corners = (
            (c0, r0, (1 - fx) * (1 - fy)),
            (c0 + 1, r0, fx * (1 - fy)),
            (c0, r0 + 1, (1 - fx) * fy),
            (c0 + 1, r0 + 1, fx * fy),
        )
        total = 0.0
        for col, row, weight in corners:
            if weight <= 1e-12:
                continue
            value = self._value(col, row)
            if math.isnan(value):
                return math.nan
            total += value * weight
        return total


class CrossSectionDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Cross Section Generator")
        self.setMinimumWidth(540)
        layout = QVBoxLayout(self)

        line_group = QGroupBox("Main line")
        line_form = QFormLayout(line_group)
        self.main_line_combo = QgsMapLayerComboBox()
        self.main_line_combo.setFilters(QgsMapLayerProxyModel.LineLayer)
        line_form.addRow("Line layer:", self.main_line_combo)
        layout.addWidget(line_group)

        crs_group = QGroupBox("Working and output CRS")
        crs_form = QFormLayout(crs_group)
        self.crs_widget = QgsProjectionSelectionWidget()
        crs_form.addRow("Projected CRS (metres):", self.crs_widget)
        self.crs_warning = QLabel()
        self.crs_warning.setWordWrap(True)
        crs_form.addRow(self.crs_warning)
        layout.addWidget(crs_group)

        parameters = QGroupBox("Distances")
        parameters_form = QFormLayout(parameters)
        self.interval_spin = self.distance_spin(10.0)
        parameters_form.addRow("Interval between sections:", self.interval_spin)
        self.half_length_spin = self.distance_spin(25.0)
        parameters_form.addRow("Distance each side of centre:", self.half_length_spin)
        self.spacing_spin = self.distance_spin(2.0)
        parameters_form.addRow("Profile sample spacing:", self.spacing_spin)
        layout.addWidget(parameters)

        raster_group = QGroupBox("Elevation rasters (optional)")
        raster_form = QFormLayout(raster_group)
        self.raster1, self.raster2 = self.raster_combo(), self.raster_combo()
        raster_form.addRow("TIF 1:", self.raster1)
        raster_form.addRow("TIF 2:", self.raster2)
        self.interpolation = QComboBox()
        self.interpolation.addItems(["Nearest", "Bilinear"])
        raster_form.addRow("Sampling:", self.interpolation)
        layout.addWidget(raster_group)

        output_group = QGroupBox("Output")
        output_form = QFormLayout(output_group)
        folder_row = QHBoxLayout()
        self.output_dir = QLineEdit()
        browse = QPushButton("Browse…")
        browse.clicked.connect(self.choose_folder)
        folder_row.addWidget(self.output_dir)
        folder_row.addWidget(browse)
        output_form.addRow("Profile CSV directory:", folder_row)
        self.output_name = QLineEdit("cross_sections")
        output_form.addRow("Section layer / file name:", self.output_name)
        self.save_gpkg = QCheckBox("Save section lines and sampled points as GeoPackage")
        output_form.addRow(self.save_gpkg)
        layout.addWidget(output_group)

        buttons = QHBoxLayout()
        self.generate_btn = QPushButton("Generate Cross Sections")
        cancel = QPushButton("Cancel")
        buttons.addWidget(self.generate_btn)
        buttons.addStretch()
        buttons.addWidget(cancel)
        layout.addLayout(buttons)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        self.generate_btn.clicked.connect(self.accept)
        cancel.clicked.connect(self.reject)
        self.main_line_combo.layerChanged.connect(self.set_default_crs)
        self.set_default_crs(self.main_line_combo.currentLayer())

    @staticmethod
    def distance_spin(default):
        spin = QDoubleSpinBox()
        spin.setRange(0.01, 10000000.0)
        spin.setDecimals(2)
        spin.setValue(default)
        spin.setSuffix(" m")
        return spin

    @staticmethod
    def raster_combo():
        combo = QgsMapLayerComboBox()
        combo.setFilters(QgsMapLayerProxyModel.RasterLayer)
        combo.setAllowEmptyLayer(True)
        combo.setLayer(None)
        return combo

    def set_default_crs(self, layer):
        if layer and layer.crs().isValid() and not layer.crs().isGeographic() and \
                layer.crs().mapUnits() == QgsUnitTypes.DistanceMeters:
            self.crs_widget.setCrs(layer.crs())
            self.crs_warning.setText("")
        else:
            self.crs_widget.setCrs(QgsCoordinateReferenceSystem("EPSG:32633"))
            if layer and layer.crs().isGeographic():
                self.crs_warning.setText(
                    "⚠ The line uses geographic coordinates. Choose a suitable "
                    "projected CRS in metres; EPSG:32633 is the initial choice."
                )
            elif layer:
                self.crs_warning.setText(
                    "⚠ The line CRS is missing or does not use metres. Choose a "
                    "suitable projected CRS in metres."
                )
            else:
                self.crs_warning.setText("")

    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Choose output directory", self.output_dir.text() or os.path.expanduser("~")
        )
        if folder:
            self.output_dir.setText(folder)

    def get_parameters(self):
        return {
            "line": self.main_line_combo.currentLayer(),
            "crs": self.crs_widget.crs(),
            "interval": self.interval_spin.value(),
            "half_length": self.half_length_spin.value(),
            "spacing": self.spacing_spin.value(),
            "rasters": [r for r in (self.raster1.currentLayer(), self.raster2.currentLayer()) if r],
            "interpolation": self.interpolation.currentText(),
            "directory": self.output_dir.text().strip(),
            "name": self.output_name.text().strip(),
            "gpkg": self.save_gpkg.isChecked(),
        }


class CrossSectionPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.dlg = None

    def initGui(self):
        icon_path = os.path.join(os.path.dirname(__file__), "icon_cross-section.png")
        self.action = QAction(
            QIcon(icon_path), "Generate Cross Sections", self.iface.mainWindow()
        )
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&Cross Section", self.action)

    def unload(self):
        self.iface.removePluginMenu("&Cross Section", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        if self.dlg is None:
            self.dlg = CrossSectionDialog(self.iface.mainWindow())
        if self.dlg.exec_():
            self.generate_cross_sections()

    @staticmethod
    def make_memory_layer(geometry_type, crs, name, fields):
        layer = QgsVectorLayer(
            "{}?crs={}".format(geometry_type, crs.authid()), name, "memory"
        )
        if not layer.isValid():
            raise ValueError("Unable to create output layer in " + crs.authid())
        layer.dataProvider().addAttributes(fields)
        layer.updateFields()
        return layer

    @staticmethod
    def save_layer(layer, path, layer_name, overwrite_file):
        options = QgsVectorFileWriter.SaveVectorOptions()
        options.driverName = "GPKG"
        options.layerName = layer_name
        options.actionOnExistingFile = (
            QgsVectorFileWriter.CreateOrOverwriteFile if overwrite_file
            else QgsVectorFileWriter.CreateOrOverwriteLayer
        )
        result = QgsVectorFileWriter.writeAsVectorFormatV3(
            layer, path, QgsProject.instance().transformContext(), options
        )
        if result[0] != QgsVectorFileWriter.NoError:
            raise RuntimeError("Cannot save {}: {}".format(layer_name, result[1]))

    def generate_cross_sections(self):
        dialog = self.dlg
        try:
            params = dialog.get_parameters()
            line_layer, crs = params["line"], params["crs"]
            if line_layer is None or not line_layer.isValid():
                raise ValueError("Select a valid line layer.")
            if not line_layer.crs().isValid():
                raise ValueError("The line layer has no valid CRS. Assign its actual CRS first.")
            if not crs.isValid() or crs.isGeographic() or \
                    crs.mapUnits() != QgsUnitTypes.DistanceMeters:
                raise ValueError("Choose a projected working CRS with metre units.")
            if not params["name"] or safe_stem(params["name"]) != params["name"]:
                raise ValueError("Use only letters, digits, _ and - in the output name.")
            rasters = params["rasters"]
            if len(rasters) == 2 and rasters[0].id() == rasters[1].id():
                raise ValueError("Select two different raster layers, or leave one empty.")
            raster_paths = []
            for raster in rasters:
                path = raster.source().split("|")[0]
                if not os.path.isfile(path) or not path.lower().endswith((".tif", ".tiff")):
                    raise ValueError("Choose local TIFF raster layers: " + raster.name())
                raster_paths.append(path)
            directory = params["directory"]
            if rasters or params["gpkg"]:
                if not os.path.isdir(directory):
                    raise ValueError("Choose an existing output directory.")
            gpkg_path = os.path.join(directory, params["name"] + ".gpkg") if params["gpkg"] else None
            if gpkg_path and os.path.exists(gpkg_path):
                raise ValueError("GeoPackage already exists: " + gpkg_path)

            section_fields = [
                QgsField("section_id", QVariant.Int),
                QgsField("source_fid", QVariant.LongLong),
                QgsField("part_no", QVariant.Int),
                QgsField("chainage_m", QVariant.Double),
                QgsField("length_m", QVariant.Double),
                QgsField("x_center", QVariant.Double),
                QgsField("y_center", QVariant.Double),
                QgsField("angle_deg", QVariant.Double),
            ]
            sections = self.make_memory_layer("LineString", crs, params["name"], section_fields)
            points = None
            if rasters and params["gpkg"]:
                points = self.make_memory_layer("Point", crs, "profile_samples", [
                    QgsField("section_id", QVariant.Int),
                    QgsField("distance_m", QVariant.Double),
                    QgsField("x", QVariant.Double),
                    QgsField("y", QVariant.Double),
                    QgsField("z_tif1", QVariant.Double),
                    QgsField("z_tif2", QVariant.Double),
                ])

            transform = None
            if line_layer.crs() != crs:
                transform = QgsCoordinateTransform(line_layer.crs(), crs, QgsProject.instance())
            parts = []
            for feature in line_layer.getFeatures():
                geometry = QgsGeometry(feature.geometry())
                if geometry.isNull() or geometry.isEmpty():
                    continue
                if QgsWkbTypes.geometryType(geometry.wkbType()) != QgsWkbTypes.LineGeometry:
                    continue
                if transform:
                    geometry.transform(transform)
                polylines = geometry.asMultiPolyline() if geometry.isMultipart() else [geometry.asPolyline()]
                for part_no, line in enumerate(polylines, start=1):
                    if len(line) < 2:
                        continue
                    length = sum(math.hypot(b.x() - a.x(), b.y() - a.y())
                                 for a, b in zip(line, line[1:]))
                    if length > 1e-8:
                        parts.append((feature.id(), part_no, line, section_chainages(length, params["interval"])))
            if not parts:
                raise ValueError("No nonempty line parts found in the selected layer.")
            total = sum(len(item[3]) for item in parts)
            raster_stems = [safe_stem(path) for path in raster_paths]
            if len(raster_stems) == 2 and raster_stems[0] == raster_stems[1]:
                raster_stems[1] += "_2"
            suffix = "_".join(raster_stems)
            names = ["cs{:05d}_{}.csv".format(i, suffix) for i in range(1, total + 1)] if rasters else []
            existing = [name for name in names if os.path.exists(os.path.join(directory, name))]
            if existing:
                raise ValueError("Profile already exists (existing files are protected): " + existing[0])
            line_name = os.path.basename(line_layer.source().split("|")[0]) or line_layer.name()

            dialog.progress_bar.setVisible(True)
            dialog.progress_bar.setMaximum(total)
            dialog.progress_bar.setValue(0)
            section_id = 0
            # Write into a staging directory so failed runs leave no incomplete profiles.
            with tempfile.TemporaryDirectory(prefix="cross_section_run_", dir=directory or None) as stage_dir:
                with tempfile.TemporaryDirectory(prefix="cross_section_rasters_") as temp_dir:
                    samplers = [RasterSampler(path, crs, params["interpolation"], temp_dir, i)
                                for i, path in enumerate(raster_paths, 1)]
                    for source_fid, part_no, line, chainages in parts:
                        for chainage in chainages:
                            tangent = tangent_at(line, chainage)
                            if tangent is None:
                                continue
                            section_id += 1
                            centre = point_on_line(line, chainage)
                            nx, ny = -tangent[1], tangent[0]
                            half = params["half_length"]
                            start = QgsPointXY(centre.x() - half * nx, centre.y() - half * ny)
                            end = QgsPointXY(centre.x() + half * nx, centre.y() + half * ny)
                            feature = QgsFeature(sections.fields())
                            feature.setGeometry(QgsGeometry.fromPolylineXY([start, end]))
                            feature.setAttributes([
                                section_id, source_fid, part_no, chainage, 2 * half,
                                centre.x(), centre.y(), math.degrees(math.atan2(ny, nx)) % 360,
                            ])
                            if not sections.dataProvider().addFeature(feature):
                                raise RuntimeError("Could not add section " + str(section_id))

                            if samplers:
                                csv_path = os.path.join(stage_dir, names[section_id - 1])
                                # Exclusive creation prevents replacement if a file appeared mid-run.
                                with open(csv_path, "x", newline="", encoding="utf-8") as output:
                                    output.write("# vector line: {}\n".format(line_name))
                                    for i, path in enumerate(raster_paths, 1):
                                        output.write("# tif{}: {}\n".format(i, os.path.basename(path)))
                                    output.write("# working CRS: {}\n".format(crs.authid()))
                                    output.write("# source feature id: {}; part: {}; chainage_m: {:.3f}\n".format(
                                        source_fid, part_no, chainage))
                                    output.write("# columns: distance_m,x,y,{}\n".format(
                                        ",".join("z_tif{}".format(i) for i in range(1, len(samplers) + 1))))
                                    writer = csv.writer(output)
                                    for distance in sample_distances(half, params["spacing"]):
                                        x, y = centre.x() + distance * nx, centre.y() + distance * ny
                                        z = [sampler.sample(x, y) for sampler in samplers]
                                        writer.writerow(["{:.3f}".format(distance), "{:.3f}".format(x),
                                                         "{:.3f}".format(y)] + [
                                                             "{:.6f}".format(v) if math.isfinite(v) else "nan"
                                                             for v in z])
                                        if points is not None:
                                            sample = QgsFeature(points.fields())
                                            sample.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
                                            sample.setAttributes([section_id, distance, x, y] + [
                                                None if not math.isfinite(v) else v
                                                for v in (z + [math.nan])[:2]])
                                            if not points.dataProvider().addFeature(sample):
                                                raise RuntimeError("Could not add sample point")
                            dialog.progress_bar.setValue(section_id)
                            QCoreApplication.processEvents()

                sections.updateExtents()
                if points is not None:
                    points.updateExtents()
                if gpkg_path:
                    self.save_layer(sections, os.path.join(stage_dir, os.path.basename(gpkg_path)), "cross_sections", True)
                    if points is not None:
                        self.save_layer(points, os.path.join(stage_dir, os.path.basename(gpkg_path)), "profile_samples", False)
                for name in names:
                    source = os.path.join(stage_dir, name)
                    target = os.path.join(directory, name)
                    if os.path.exists(target):
                        raise FileExistsError("Profile appeared during the run: " + target)
                    os.rename(source, target)
                if gpkg_path:
                    if os.path.exists(gpkg_path):
                        raise FileExistsError("GeoPackage appeared during the run: " + gpkg_path)
                    os.rename(os.path.join(stage_dir, os.path.basename(gpkg_path)), gpkg_path)
            QgsProject.instance().addMapLayer(sections)
            if points is not None:
                QgsProject.instance().addMapLayer(points)
            QMessageBox.information(self.iface.mainWindow(), "Cross sections",
                                    "Created {} cross sections{}{}.".format(
                                        section_id,
                                        " and {} CSV profiles".format(section_id) if samplers else "",
                                        " in " + directory if directory and (samplers or gpkg_path) else ""))
        except Exception as error:
            QgsMessageLog.logMessage(str(error), "Cross Section")
            QMessageBox.critical(self.iface.mainWindow(), "Cross sections", str(error))
        finally:
            dialog.progress_bar.setVisible(False)
