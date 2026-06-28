import argparse
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
}


def _clean(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def _collect_text_from_run(node):
    parts = []
    for t in node.findall(".//w:t", NS):
        if t.text:
            parts.append(t.text)
    return "".join(parts)


def _collect_text_from_math(node):
    parts = []
    for t in node.findall(".//m:t", NS):
        if t.text:
            parts.append(t.text)
    return "".join(parts)


def extract_docx(docx_path: Path) -> str:
    with zipfile.ZipFile(docx_path, "r") as zf:
        xml_bytes = zf.read("word/document.xml")
    root = ET.fromstring(xml_bytes)

    lines = []
    eq_count = 0
    para_count = 0

    for p in root.findall(".//w:body/w:p", NS):
        para_parts = []
        math_parts = []

        for child in list(p):
            tag = child.tag
            if tag.endswith("}r"):
                txt = _collect_text_from_run(child)
                if txt:
                    para_parts.append(txt)
            elif tag.endswith("}oMath") or tag.endswith("}oMathPara"):
                mtxt = _collect_text_from_math(child)
                if mtxt:
                    math_parts.append(mtxt)

        plain = _clean("".join(para_parts))
        if plain:
            para_count += 1
            lines.append(plain)

        for mtxt in math_parts:
            mtxt = _clean(mtxt)
            if mtxt:
                eq_count += 1
                lines.append(f"[EQ{eq_count:04d}] {mtxt}")

    header = [
        f"# Extracted from: {docx_path}",
        f"# Paragraphs: {para_count}",
        f"# Equations detected: {eq_count}",
        "",
    ]
    return "\n".join(header + lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Extract text and equation text from .docx")
    parser.add_argument("docx_files", nargs="+", help="Path(s) to .docx files")
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for extracted .txt files",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for item in args.docx_files:
        src = Path(item).resolve()
        if not src.exists():
            print(f"SKIP (missing): {src}")
            continue
        if src.suffix.lower() != ".docx":
            print(f"SKIP (not docx): {src}")
            continue
        text = extract_docx(src)
        dst = out_dir / f"{src.stem}.extracted.txt"
        dst.write_text(text, encoding="utf-8")
        print(f"WROTE: {dst}")


if __name__ == "__main__":
    main()

