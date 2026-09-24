"""Plot CSV profiles exported by cross_section_tool.py.

Install: pip install matplotlib
Put this script beside the cs*.csv files, or change INPUT_DIR below.
Run: python plot_cross_sections.py
"""

import csv
import math
from pathlib import Path
from statistics import mean

import matplotlib
matplotlib.use("Agg")  # Save PNGs without needing a graphical desktop.
import matplotlib.pyplot as plt


# ---- Change these settings if needed --------------------------------------
INPUT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = INPUT_DIR / "plot_cross-sections"
DATUM = 0.0  # Same vertical units/datum as the rasters.
DEPTH_IS_NEGATIVE_Z = True  # True: negative z is underwater. False: positive z is depth.
# ---------------------------------------------------------------------------


def depth(z):
    """Water depth below DATUM; dry points count as zero in section averages."""
    return max(0.0, DATUM - z if DEPTH_IS_NEGATIVE_Z else z - DATUM)


def average(values):
    return mean(values) if values else math.nan


def label_number(value):
    return "{:.3f}".format(value) if math.isfinite(value) else "n/a"


def read_profile(path):
    metadata = {}
    rows = []
    with path.open(newline="", encoding="utf-8-sig") as file:
        for line in file:
            if line.startswith("#"):
                key, _, value = line[1:].strip().partition(":")
                metadata[key.strip().lower()] = value.strip()
            elif line.strip():
                rows.append(next(csv.reader([line])))
    columns = metadata.get("columns", "").split(",")
    columns = [column.strip() for column in columns]
    if "distance_m" not in columns or "z_tif1" not in columns:
        raise ValueError("Missing '# columns: distance_m,x,y,z_tif1...' header")
    if any(len(row) != len(columns) for row in rows):
        raise ValueError("A data row has a different number of columns")
    if not rows:
        raise ValueError("No samples")
    data = {key: [float(row[i]) for row in rows] for i, key in enumerate(columns)}
    return metadata, data


def plot_profile(path):
    meta, data = read_profile(path)
    x = data["distance_m"]
    # The generator always exports a sample at signed distance 0: the centre.
    centre_rows = [i for i, distance in enumerate(x) if abs(distance) < 1e-8]
    if len(centre_rows) != 1 or "x" not in data or "y" not in data:
        raise ValueError("Expected one centre sample (distance_m=0) with x,y coordinates")
    centre_x = data["x"][centre_rows[0]]
    centre_y = data["y"][centre_rows[0]]
    if not math.isfinite(centre_x) or not math.isfinite(centre_y):
        raise ValueError("Centre coordinates are not finite")
    z1 = data["z_tif1"]
    z2 = data.get("z_tif2")
    tif1 = meta.get("tif1", "TIF 1")
    tif2 = meta.get("tif2", "TIF 2") if z2 is not None else ""
    # The exporter stores source feature, part and chainage in one comment.
    source_comment = meta.get("source feature id", "")
    source_fid = source_comment.split(";", 1)[0].strip()
    source_info = {}
    for piece in source_comment.split(";")[1:]:
        if ":" in piece:
            key, _, value = piece.partition(":")
            source_info[key.strip()] = value.strip()
    part_no = source_info.get("part", "")
    chainage = source_info.get("chainage_m", "unknown")
    section = path.name.split("_", 1)[0].upper()

    valid1 = [v for v in z1 if math.isfinite(v)]
    valid2 = [v for v in z2 if math.isfinite(v)] if z2 is not None else []
    mean_depth1 = average([depth(v) for v in valid1])
    mean_depth2 = average([depth(v) for v in valid2])
    paired = [(a, b) for a, b in zip(z1, z2)
              if math.isfinite(a) and math.isfinite(b)] if z2 is not None else []

    stats = {
        "section": section, "source_csv": path.name,
        "vector_line": meta.get("vector line", ""),
        "source_feature_id": source_fid, "part_no": part_no,
        "chainage_m": chainage, "center_x": centre_x, "center_y": centre_y,
        "working_crs": meta.get("working crs", ""),
        "tif1": tif1, "tif2": tif2,
        "samples_tif1": len(valid1), "samples_tif2": len(valid2),
        "mean_z_tif1": average(valid1), "mean_z_tif2": average(valid2),
        "mean_depth_tif1": mean_depth1, "mean_depth_tif2": mean_depth2,
        "paired_samples": len(paired),
        "mean_z_difference_tif2_minus_tif1": average([b - a for a, b in paired]),
        "rmse_z": math.sqrt(mean([(b - a) ** 2 for a, b in paired])) if paired else math.nan,
        "mean_depth_difference_tif2_minus_tif1": average(
            [depth(b) - depth(a) for a, b in paired]
        ),
    }

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.plot(x, z1, color="tab:blue", linewidth=1.8, label=tif1)
    if z2 is not None:
        ax.plot(x, z2, color="tab:orange", linewidth=1.8, label=tif2)
    ax.axhline(DATUM, color="gray", linestyle=":", linewidth=1, label="Depth datum")
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_title("{}  |  Chainage: {} m  |  Centre: ({:.3f}, {:.3f}) {}\n{}".format(
        section, chainage, centre_x, centre_y, meta.get("working crs", ""),
        tif1 + ("  vs  " + tif2 if z2 is not None else "")
    ))
    ax.set_xlabel("Signed distance from centre (m)")
    ax.set_ylabel("Elevation / bathymetry z (raster units)")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)

    note = "{}: mean depth {} (n={})".format(tif1, label_number(mean_depth1), len(valid1))
    if z2 is not None:
        note += "\n{}: mean depth {} (n={})".format(
            tif2, label_number(mean_depth2), len(valid2)
        )
        note += "\nPaired Δdepth (TIF2 − TIF1): {}  |  RMSE(z): {}  |  n={}".format(
            label_number(stats["mean_depth_difference_tif2_minus_tif1"]),
            label_number(stats["rmse_z"]), len(paired)
        )
    ax.text(0.02, 0.8, note, transform=ax.transAxes, va="bottom", fontsize=8.5,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / (path.stem + ".png"), dpi=170)
    plt.close(fig)
    return stats


def main():
    paths = sorted(INPUT_DIR.glob("cs*.csv"))
    if not paths:
        raise SystemExit("No cs*.csv files found in {}".format(INPUT_DIR))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    for path in paths:
        try:
            results.append(plot_profile(path))
            print("Saved:", OUTPUT_DIR / (path.stem + ".png"))
        except (ValueError, OSError, KeyError) as error:
            print("Skipped {}: {}".format(path.name, error))
    if results:
        with (OUTPUT_DIR / "cross_section_stats.csv").open(
            "w", newline="", encoding="utf-8"
        ) as file:
            writer = csv.DictWriter(file, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print("Saved:", OUTPUT_DIR / "cross_section_stats.csv")


if __name__ == "__main__":
    main()
