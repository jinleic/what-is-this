"""Read Fig11 vector regions exactly; test a conditional local-order reading."""
from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from fractions import Fraction as F
from pathlib import Path


NUM = r"[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?"
SPOKES = ("N", "E", "S", "W", "NE", "SW")


def points(element):
    tokens = re.findall(r"[A-Za-z]|" + NUM, element.attrib["d"])
    result = []
    index = 0
    while index < len(tokens):
        command = tokens[index]
        index += 1
        assert command in ("M", "L", "C", "Z"), command
        count = {"M": 2, "L": 2, "C": 6, "Z": 0}[command]
        values = list(map(F, tokens[index:index + count]))
        index += count
        result.extend(zip(values[::2], values[1::2]))
    transform = element.get("transform")
    if transform:
        assert transform.startswith("matrix(")
        a, b, c, d, e, f = map(F, re.findall(NUM, transform))
        result = [(a * x + c * y + e, b * x + d * y + f) for x, y in result]
    return result


def centre(vertices):
    return tuple((min(p[i] for p in vertices) + max(p[i] for p in vertices)) / 2 for i in (0, 1))


def inside(point, polygon):
    x, y = point
    answer = False
    for (ax, ay), (bx, by) in zip(polygon, polygon[1:] + polygon[:1]):
        if (ax, ay) == (bx, by):
            continue
        cross = (x - ax) * (by - ay) - (y - ay) * (bx - ax)
        assert cross or not (min(ax, bx) <= x <= max(ax, bx) and min(ay, by) <= y <= max(ay, by)), "Node on residual boundary"
        if (ay > y) != (by > y):
            crossing_x = ax + (y - ay) * (bx - ax) / (by - ay)
            if x < crossing_x:
                answer = not answer
    return answer


def colour(value):
    if not value or not value.startswith("rgb("):
        return None
    red, green, blue = map(F, re.findall(NUM, value))
    if red > 90 and green < 10 and blue < 10:
        return "red"
    if green > 50 and green > red and green > blue:
        return "green"
    if blue > 60 and blue > green and blue > red:
        return "blue"
    return None


def normalize(support):
    mask = sum(1 << SPOKES.index(s) for s in support)
    return min(mask, 63 ^ mask)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    root = ET.parse(args.run / "figure11.svg").getroot()
    paths = root.findall("{http://www.w3.org/2000/svg}path")
    centres = []
    nodes = []
    polygons = []
    for element in paths:
        opacity = F(element.get("fill-opacity", "1"))
        fill = element.get("fill")
        if opacity < 1 and fill not in (None, "none"):
            assert "C" not in element.get("d"), "Residual is not polygonal"
            polygons.append(points(element))
        elif opacity == 1:
            marker_colour = colour(fill) or colour(element.get("stroke"))
            if marker_colour and "C" not in element.get("d"):
                centres.append({"colour": marker_colour, "point": centre(points(element))})
            elif fill in ("rgb(0%, 0%, 0%)", "rgb(100%, 100%, 100%)") and "C" in element.get("d"):
                vertices = points(element)
                width = max(p[0] for p in vertices) - min(p[0] for p in vertices)
                height = max(p[1] for p in vertices) - min(p[1] for p in vertices)
                assert 6 < width < 12 and 6 < height < 12
                nodes.append(centre(vertices))
    assert len(centres) == 6 and len(nodes) == 36 and len(polygons) == 18, (len(centres), len(nodes), len(polygons))
    sorted_y = sorted(c["point"][1] for c in centres)
    middle_y = (sorted_y[2] + sorted_y[3]) / 2
    for c in centres:
        c["panel"] = ("Z" if c["point"][1] < middle_y else "X") + "_" + c["colour"]
    node_labels = []
    for node in nodes:
        c = min(centres, key=lambda p: sum((a - b) ** 2 for a, b in zip(p["point"], node)))
        dx, dy = (node[i] - c["point"][i] for i in (0, 1))
        if abs(dx) < 1:
            spoke = "N" if dy < 0 else "S"
        elif abs(dy) < 1:
            spoke = "W" if dx < 0 else "E"
        else:
            assert dx * dy < 0
            spoke = "NE" if dx > 0 else "SW"
        node_labels.append({"panel": c["panel"], "spoke": spoke, "point": node})
    panel_regions = {c["panel"]: [] for c in centres}
    raw_regions = []
    for polygon in polygons:
        members = [n for n in node_labels if inside(n["point"], polygon)]
        assert len(members) in (2, 3), members
        panel = members[0]["panel"]
        assert all(n["panel"] == panel for n in members)
        support = sorted(n["spoke"] for n in members)
        panel_regions[panel].append(support)
        raw_regions.append({"panel": panel, "support": support, "vertices": [[str(x), str(y)] for x, y in polygon]})
    annotations = json.loads((args.run / "inputs/legacy_transcription.json").read_text())
    comparisons = []
    for panel, regions in sorted(panel_regions.items()):
        assert sorted(map(len, regions)) == [2, 2, 3]
        kind, shade = panel.split("_")
        style = "filled" if kind == "Z" else "hollow"
        values = dict(zip(annotations["edge_order"], annotations["figure12_times"][style][shade]))
        order = sorted(values, key=values.get)
        expected = [order[:2], order[:3], order[-2:]]
        comparisons.append({"panel": panel, "figure11_regions": regions, "conditional_ascending_annotation_order": order, "conditional_residuals": expected, "match_modulo_check": sorted(map(normalize, regions)) == sorted(map(normalize, expected))})
    report = {"source": "Pinned Figure11 PDF exported with pdftocairo -svg", "method": "Exact rational transforms and polygon membership of circle centers", "centres": [{**c, "point": [str(v) for v in c["point"]]} for c in centres], "nodes": [{**n, "point": [str(v) for v in n["point"]]} for n in node_labels], "raw_regions": raw_regions, "comparisons": comparisons, "all_six_match": all(c["match_modulo_check"] for c in comparisons), "scope": "Conditional local-order correspondence only; no absolute timing or physical circuit admission"}
    (args.run / "residual_vector_correspondence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"comparisons": comparisons, "all_six_match": report["all_six_match"], "scope": report["scope"]}, indent=2))


if __name__ == "__main__":
    main()
