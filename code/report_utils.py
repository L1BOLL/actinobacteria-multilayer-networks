from __future__ import annotations

from pathlib import Path
import pandas as pd


def replace_section(report_path: Path, sentinel: str, heading: str, body_md: str) -> None:
    """Replace content from ``sentinel`` up to the next the next section marker."""
    text = report_path.read_text(encoding="utf-8")
    section = f"{sentinel}\n\n## {heading}\n\n{body_md.strip()}\n"
    if sentinel in text:
        head, tail = text.split(sentinel, 1)
        next_marker = tail.find("\n<!-- ")
        if next_marker != -1:
            remainder = tail[next_marker + 1 :]
        else:
            remainder = ""
        text = head + section + remainder
    else:
        if not text.endswith("\n"):
            text += "\n"
        text += "\n" + section
    report_path.write_text(text, encoding="utf-8")


def append_summary(report_path: Path, summary_md: str) -> None:
    text = report_path.read_text(encoding="utf-8")
    marker = "## 11. What this changes in the paper\n\n"
    if marker in text:
        head, _ = text.split(marker, 1)
        text = head + marker + summary_md.strip() + "\n"
    else:
        text += "\n" + marker + summary_md.strip() + "\n"
    report_path.write_text(text, encoding="utf-8")


def dataframe_to_md(df: pd.DataFrame, index: bool = False) -> str:
    if df.empty:
        return "_No rows._"

    if index:
        cols = [df.index.name or "index", *df.columns.tolist()]
        rows = [[idx, *row.tolist()] for idx, row in df.iterrows()]
    else:
        cols = df.columns.tolist()
        rows = [row.tolist() for _, row in df.iterrows()]

    def fmt(x: object) -> str:
        if isinstance(x, float):
            return f"{x:.6g}"
        return str(x)

    header = "| " + " | ".join(map(str, cols)) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in rows]
    return "\n".join([header, sep, *body])
