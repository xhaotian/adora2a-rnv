#!/usr/bin/env python3
"""Shared style, measured layout QA, and PLOS export helpers."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.text import Text
import pandas as pd
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
CONTRACT = json.loads((HERE / "figure_contract.json").read_text(encoding="utf-8"))
STYLE = CONTRACT["style"]
QA = CONTRACT["qa"]

BLUE, ORANGE, GREEN, MAGENTA = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
GREY, LIGHT, DARK = "#6F6F6F", "#D4D4D4", "#222222"


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [STYLE["font_family"], STYLE["font_fallback"]],
        "font.size": STYLE["tick_pt"],
        "axes.titlesize": STYLE["panel_heading_pt"],
        "axes.titleweight": "semibold",
        "axes.labelsize": STYLE["axis_label_pt"],
        "xtick.labelsize": STYLE["tick_pt"],
        "ytick.labelsize": STYLE["tick_pt"],
        "legend.fontsize": STYLE["legend_pt"],
        "axes.linewidth": STYLE["line_width_pt"],
        "lines.linewidth": STYLE["line_width_pt"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def panel_label(ax, label: str) -> Text:
    text = ax.text(
        QA["panel_label_x"], QA["panel_label_y"], label,
        transform=ax.transAxes, fontsize=STYLE["panel_label_pt"],
        fontweight="bold", va="bottom", ha="left", clip_on=False,
    )
    text.set_gid(f"panel_label:{label}")
    return text


def panel_title(ax, label: str, title: str, pad: float = 8) -> Text:
    text = ax.set_title(title, loc="left", pad=pad, fontsize=STYLE["panel_heading_pt"], fontweight="semibold")
    text.set_gid(f"panel_title:{label}")
    return text


def reference_line(ax, x: float = 0, *, linestyle: str = "-", linewidth: float | None = None):
    return ax.axvline(
        x, color=STYLE["reference_line_color"], linestyle=linestyle,
        linewidth=linewidth or STYLE["reference_line_width_pt"], zorder=0,
    )


def register_box_text(fig, patch, text: Text, name: str) -> None:
    fig.__dict__.setdefault("_qa_box_text", []).append((patch, text, name))


def register_legend_data(fig, legend, ax, xs, ys, name: str, radius_pt: float = 5) -> None:
    fig.__dict__.setdefault("_qa_legend_data", []).append((legend, ax, list(xs), list(ys), name, radius_pt))


def _bbox_string(bbox) -> str:
    return f"{bbox.x0:.1f},{bbox.y0:.1f},{bbox.x1:.1f},{bbox.y1:.1f}"


def _intersects(a, b) -> bool:
    return a.x0 < b.x1 and a.x1 > b.x0 and a.y0 < b.y1 and a.y1 > b.y0


def audit_figure(fig, figure_name: str) -> pd.DataFrame:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    tol = float(QA["canvas_tolerance_px"])
    rows: list[dict[str, str]] = []

    def add(element: str, bbox, violation: str, passed: bool) -> None:
        rows.append({
            "figure": figure_name,
            "element": element,
            "bbox": _bbox_string(bbox),
            "violation_type": "NONE" if passed else violation,
            "status": "PASS" if passed else "FAIL",
        })

    visible_text = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for index, text in enumerate(visible_text):
        bbox = text.get_window_extent(renderer)
        inside = (bbox.x0 >= canvas.x0 - tol and bbox.y0 >= canvas.y0 - tol and
                  bbox.x1 <= canvas.x1 + tol and bbox.y1 <= canvas.y1 + tol)
        add(f"text:{text.get_gid() or index}:{text.get_text()[:45]}", bbox, "TEXT_OUTSIDE_CANVAS", inside)
        add(f"font:{text.get_gid() or index}:{text.get_text()[:45]}", bbox, "FONT_BELOW_8_PT",
            float(text.get_fontsize()) >= float(QA["minimum_font_pt"]))

    by_gid = {t.get_gid(): t for t in visible_text if t.get_gid()}
    labels = {k.split(":", 1)[1]: v for k, v in by_gid.items() if k.startswith("panel_label:")}
    titles = {k.split(":", 1)[1]: v for k, v in by_gid.items() if k.startswith("panel_title:")}
    for label in sorted(set(labels) & set(titles)):
        lb = labels[label].get_window_extent(renderer)
        tb = titles[label].get_window_extent(renderer)
        add(f"panel_label_title:{label}", lb, "PANEL_LABEL_TITLE_INTERSECTION", not _intersects(lb, tb))

    axes = [ax for ax in fig.axes if ax.get_visible()]
    for index, ax in enumerate(axes):
        tick_text = [*ax.get_xticklabels(), *ax.get_yticklabels()]
        for tick_index, text in enumerate(tick_text):
            if not text.get_visible() or not text.get_text().strip():
                continue
            bbox = text.get_window_extent(renderer)
            intrudes = any(_intersects(bbox, other.bbox) for other in axes if other is not ax)
            add(f"tick_neighbor:axes{index}:{tick_index}:{text.get_text()}", bbox,
                "TICK_LABEL_IN_NEIGHBOR_AXES", not intrudes)

    required_padding = float(QA["box_text_padding_mm"]) * fig.dpi / 25.4
    for patch, text, name in fig.__dict__.get("_qa_box_text", []):
        pb = patch.get_window_extent(renderer)
        tb = text.get_window_extent(renderer)
        padding = min(tb.x0 - pb.x0, pb.x1 - tb.x1, tb.y0 - pb.y0, pb.y1 - tb.y1)
        add(f"box_padding:{name}", tb, "BOX_TEXT_PADDING_BELOW_1.5_MM", padding >= required_padding)

    registered_legends = fig.__dict__.get("_qa_legend_data", [])
    for legend, ax, xs, ys, name, radius_pt in registered_legends:
        lb = legend.get_window_extent(renderer)
        radius_px = radius_pt * fig.dpi / 72.0
        overlaps = False
        for x, y in zip(xs, ys):
            px, py = ax.transData.transform((x, y))
            if lb.x0 - radius_px <= px <= lb.x1 + radius_px and lb.y0 - radius_px <= py <= lb.y1 + radius_px:
                overlaps = True
                break
        add(f"legend_data:{name}", lb, "LEGEND_OVERLAPS_DATA_POINT", not overlaps)

    registered_ids = {id(item[0]) for item in registered_legends}
    for index, legend in enumerate(fig.findobj(Legend)):
        if not legend.get_visible() or id(legend) in registered_ids:
            continue
        lb = legend.get_window_extent(renderer)
        intersects_plotting_area = any(_intersects(lb, ax.bbox) for ax in axes)
        add(f"legend_axes:{index}", lb, "LEGEND_OVERLAPS_PLOTTING_AREA", not intersects_plotting_area)

    if not rows:
        rows.append({"figure": figure_name, "element": "figure", "bbox": _bbox_string(canvas),
                     "violation_type": "NONE", "status": "PASS"})
    return pd.DataFrame(rows)


def _qa_overlay(fig, qa: pd.DataFrame, path: Path) -> None:
    fig.canvas.draw()
    image = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).convert("RGB")
    draw = ImageDraw.Draw(image)
    height = image.height
    for _, row in qa.iterrows():
        if row["element"].startswith(("panel_label_title", "box_padding", "legend_data")) or row["status"] == "FAIL":
            x0, y0, x1, y1 = [float(v) for v in row["bbox"].split(",")]
            color = "#D62728" if row["status"] == "FAIL" else "#2CA02C"
            draw.rectangle((x0, height - y1, x1, height - y0), outline=color, width=2)
    image.save(path)


def save_main_figure(fig, stem: str, out: Path) -> pd.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    (out / "review_png").mkdir(exist_ok=True)
    (out / "qa_overlays").mkdir(exist_ok=True)
    qa = audit_figure(fig, stem)
    qa.to_csv(out / f"{stem}_layout_qa.tsv", sep="\t", index=False)
    _qa_overlay(fig, qa, out / "qa_overlays" / f"{stem}_QA_overlay.png")
    if qa["status"].eq("FAIL").any():
        failures = qa.loc[qa.status.eq("FAIL"), ["element", "violation_type"]]
        raise RuntimeError(f"{stem} layout QA failed:\n{failures.to_string(index=False)}")
    fig.savefig(out / f"{stem}.pdf", facecolor="white")
    preview = out / "review_png" / f"{stem}.png"
    fig.savefig(preview, dpi=CONTRACT["output"]["review_dpi"], facecolor="white")
    png = out / f".{stem}.png"
    fig.savefig(png, dpi=CONTRACT["output"]["dpi"], facecolor="white")
    with Image.open(png) as im:
        im.convert("RGB").save(out / f"{stem}.tif", compression="tiff_lzw",
                               dpi=(CONTRACT["output"]["dpi"], CONTRACT["output"]["dpi"]))
    png.unlink()
    plt.close(fig)
    return qa


def save_supplement_page(fig, pdf, stem: str, out: Path) -> pd.DataFrame:
    qa = audit_figure(fig, stem)
    qa.to_csv(out / f"{stem}_layout_qa.tsv", sep="\t", index=False)
    (out / "review_png").mkdir(exist_ok=True)
    (out / "qa_overlays").mkdir(exist_ok=True)
    _qa_overlay(fig, qa, out / "qa_overlays" / f"{stem}_QA_overlay.png")
    if qa["status"].eq("FAIL").any():
        failures = qa.loc[qa.status.eq("FAIL"), ["element", "violation_type"]]
        raise RuntimeError(f"{stem} layout QA failed:\n{failures.to_string(index=False)}")
    fig.savefig(out / "review_png" / f"{stem}.png", dpi=CONTRACT["output"]["review_dpi"], facecolor="white")
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)
    return qa


apply_style()
