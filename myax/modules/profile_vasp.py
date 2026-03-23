#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional


@dataclass
class Node:
    name: str
    time: float
    calls: int
    level: int
    parent: Optional["Node"] = None
    children: list["Node"] = field(default_factory=list)

    def path(self) -> str:
        parts = []
        cur: Optional["Node"] = self
        while cur is not None:
            parts.append(cur.name.strip())
            cur = cur.parent
        return " > ".join(reversed(parts))


TREE_ROW = re.compile(
    r"^\s*(\d+)\s+(.*?)\s+([0-9]+(?:\.[0-9]+)?)\s+(\d+)\s+(\d+)\s*$"
)

# OUTCAR: "Iteration      1(  25)"  -> take the number inside parentheses
ITERATION_LINE = re.compile(r"Iteration\s+\d+\(\s*(\d+)\s*\)", re.IGNORECASE)


def parse_tree(lines: Iterable[str]) -> Node:
    lines = list(lines)
    start = None
    for i, s in enumerate(lines):
        if s.strip().startswith("index") and "routine" in s and "level" in s:
            start = i + 2
            break
    if start is None:
        raise ValueError("Tree header (index/routine/time/calls/level) not found in OUTCAR.")

    rows: list[tuple[str, float, int, int]] = []
    for s in lines[start:]:
        if not s.strip():
            break
        if s.lstrip().startswith("Flat profile"):
            break
        m = TREE_ROW.match(s)
        if not m:
            continue
        _, name, t, calls, level = m.groups()
        rows.append((name, float(t), int(calls), int(level)))

    if not rows:
        raise ValueError("No parsable tree rows found.")

    root: Optional[Node] = None
    stack: list[Optional[Node]] = [None]

    for name, t, calls, level in rows:
        node = Node(name=name, time=t, calls=calls, level=level)
        if level == 1:
            root = node
            stack = [None, root]
            continue

        while len(stack) <= level:
            stack.append(None)

        parent = stack[level - 1]
        if parent is None:
            raise ValueError(f"Broken hierarchy at '{name}' (level={level}).")

        node.parent = parent
        parent.children.append(node)

        stack[level] = node
        for k in range(level + 1, len(stack)):
            stack[k] = None

    if root is None:
        raise ValueError("Root node not found.")
    return root

def total_iterations_from_outcar(lines: Iterable[str]) -> int:
    total = 0
    for s in lines:
        m = ITERATION_LINE.search(s)   # match ではなく search
        if m:
            total += int(m.group(1))
    return total

def fmt_seconds_aligned(x: float) -> str:
    return f"{x:12.6f} s"


def is_under_elmin(node: Node) -> bool:
    cur: Optional[Node] = node
    while cur is not None:
        if cur.name.strip() == "elmin":
            return True
        cur = cur.parent
    return False


def show_children(node: Node, top: int = 20, total_iters: int = 0) -> list[Node]:
    children = sorted(node.children, key=lambda n: n.time, reverse=True)
    if top > 0:
        children = children[:top]

    print()
    print(node.path())
    print(f"time: {fmt_seconds_aligned(node.time)}   calls: {node.calls}   level: {node.level}")

    if total_iters > 0 and is_under_elmin(node):
        avg = node.time / total_iters
        print(f"iterations(total): {total_iters}   avg per iteration: {fmt_seconds_aligned(avg)}")

    print()

    if not children:
        print("(no children)")
        return []

    denom = node.time if node.time > 0 else 1.0
    for i, c in enumerate(children, 1):
        pct = 100.0 * c.time / denom
        print(
            f"{i:>2}. {c.name:<35} {fmt_seconds_aligned(c.time):>14}  {pct:6.2f}%   calls={c.calls}"
        )
    return children


def repl(root: Node, top: int, total_iters: int, shown_initial: list[Node]) -> None:
    cur = root
    shown = shown_initial

    cmd = "Commands: number=enter, b=back, r=root, a=all children, t=top N, p=path, q=quit"

    while True:
        print()
        print(cmd)
        s = input("> ").strip()
        if not s:
            continue

        if s in {"q", "quit", "exit"}:
            return

        if s in {"b", "back"}:
            if cur.parent is None:
                print("(already at root)")
            else:
                cur = cur.parent
                shown = show_children(cur, top=top, total_iters=total_iters)
            continue

        if s in {"r", "root"}:
            cur = root
            shown = show_children(cur, top=top, total_iters=total_iters)
            continue

        if s in {"p", "path"}:
            print(cur.path())
            continue

        if s in {"a", "all"}:
            shown = show_children(cur, top=0, total_iters=total_iters)
            continue

        if s.startswith("t"):
            parts = s.split()
            if len(parts) == 2 and parts[1].isdigit():
                top = int(parts[1])
                shown = show_children(cur, top=top, total_iters=total_iters)
            else:
                print("Usage: t 20")
            continue

        if s.isdigit():
            k = int(s)
            if 1 <= k <= len(shown):
                cur = shown[k - 1]
                shown = show_children(cur, top=top, total_iters=total_iters)
            else:
                print("out of range")
            continue

        print("Unknown command")


def main() -> None:
    ap = argparse.ArgumentParser(description="Interactive analyzer for VASP PROFILING tree in OUTCAR.")
    ap.add_argument("file", nargs="?", default="OUTCAR", help="VASP OUTCAR (default: OUTCAR)")
    ap.add_argument("--top", type=int, default=20, help="How many children to show by default (0 = all)")
    args = ap.parse_args()

    path = Path(args.file)
    lines = path.read_text(errors="replace").splitlines()

    root = parse_tree(lines)
    total_iters = total_iterations_from_outcar(lines)

    print()
    print("Summary (level-2 under total_time):")
    shown_initial = show_children(root, top=args.top, total_iters=total_iters)

    repl(root, top=args.top, total_iters=total_iters, shown_initial=shown_initial)


if __name__ == "__main__":
    main()
