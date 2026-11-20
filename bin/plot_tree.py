#!/usr/bin/env python3
"""Render a compact rectangular PNG from a Newick tree without Biopython."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0


def parse_newick(text: str) -> Node:
    index = 0

    def token() -> str:
        nonlocal index
        start = index
        while index < len(text) and text[index] not in ",():;":
            index += 1
        return text[start:index].strip()

    def node() -> Node:
        nonlocal index
        children: list[Node] = []
        if index < len(text) and text[index] == "(":
            index += 1
            while True:
                children.append(node())
                if text[index] == ",":
                    index += 1
                    continue
                if text[index] == ")":
                    index += 1
                    break
        name = token()
        length = 0.0
        if index < len(text) and text[index] == ":":
            index += 1
            raw = token()
            try:
                length = float(raw)
            except ValueError:
                length = 0.0
        return Node(name=name, length=length, children=children)

    return node()


def layout(root: Node) -> list[Node]:
    leaves: list[Node] = []

    def walk(current: Node, parent_x: float) -> None:
        current.x = parent_x + max(current.length, 0.0)
        if not current.children:
            current.y = float(len(leaves))
            leaves.append(current)
            return
        for child in current.children:
            walk(child, current.x)
        current.y = sum(child.y for child in current.children) / len(current.children)

    walk(root, 0.0)
    if max((leaf.x for leaf in leaves), default=0.0) == 0:
        def depth(current: Node, level: int) -> None:
            current.x = float(level)
            for child in current.children:
                depth(child, level + 1)
        depth(root, 0)
    return leaves


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = parse_newick(args.tree.read_text().strip())
    leaves = layout(root)
    fig, ax = plt.subplots(figsize=(11, max(4, len(leaves) * 0.22)))

    def draw(current: Node) -> None:
        if current.children:
            ys = [child.y for child in current.children]
            ax.plot([current.x, current.x], [min(ys), max(ys)], color="black", linewidth=0.8)
            for child in current.children:
                ax.plot([current.x, child.x], [child.y, child.y], color="black", linewidth=0.8)
                draw(child)
    draw(root)
    offset = max((leaf.x for leaf in leaves), default=1.0) * 0.01 or 0.05
    for leaf in leaves:
        ax.text(leaf.x + offset, leaf.y, leaf.name, va="center", fontsize=7)
    ax.set_title(args.title)
    ax.set_yticks([])
    ax.set_xlabel("Branch length")
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
