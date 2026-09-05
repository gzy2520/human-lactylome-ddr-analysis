#!/usr/bin/env python3
"""Apply audited Methods corrections as Word tracked changes.

This deliberately edits only text nodes in the Methods paragraphs of the
manuscript.  The source DOCX is copied to a dated output first; existing
revisions elsewhere in the document are left untouched.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "DDR_Kla_manuscript_V3.1.docx"
OUTPUT = ROOT / "manuscript" / "DDR_Kla_manuscript_V3.1_methods_revised_gzy_20260905.docx"

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_NS = "http://www.w3.org/XML/1998/namespace"
NS = {"w": W_NS}


def w(local: str) -> str:
    return f"{{{W_NS}}}{local}"


def local_name(tag: str) -> str:
    return etree.QName(tag).localname


def xml_space(text: str) -> str | None:
    return "preserve" if text and (text[0].isspace() or text[-1].isspace()) else None


def make_run(template: etree._Element, text: str, text_tag: str = "t") -> etree._Element:
    """Make a visually equivalent run containing text or deleted text."""

    run = etree.Element(w("r"))
    for key, value in template.attrib.items():
        run.set(key, value)
    rpr = template.find(w("rPr"))
    if rpr is not None:
        run.append(deepcopy(rpr))
    text_node = etree.SubElement(run, w(text_tag))
    space = xml_space(text)
    if space:
        text_node.set(f"{{{XML_NS}}}space", space)
    text_node.text = text
    return run


def make_revision(
    template: etree._Element,
    kind: str,
    revision_id: int,
    text: str,
    author: str,
    date: str,
) -> etree._Element:
    if kind not in {"ins", "del"}:
        raise ValueError(kind)
    revision = etree.Element(
        w(kind),
        {
            w("id"): str(revision_id),
            w("author"): author,
            w("date"): date,
        },
    )
    revision.append(make_run(template, text, "delText" if kind == "del" else "t"))
    return revision


def next_revision_id(root: etree._Element) -> int:
    values: list[int] = []
    for element in root.iter():
        for key, value in element.attrib.items():
            if etree.QName(key).localname != "id":
                continue
            try:
                values.append(int(value))
            except (TypeError, ValueError):
                continue
    return max(values, default=0) + 1


def text_nodes(paragraph: etree._Element) -> list[etree._Element]:
    # Methods paragraphs being corrected contain ordinary runs.  Restricting
    # this to direct runs prevents accidentally changing text in an existing
    # revision or an equation object.
    return paragraph.xpath("./w:r/w:t", namespaces=NS)


def replace_in_paragraph(
    paragraph: etree._Element,
    old: str,
    new: str,
    revision_id: int,
    author: str,
    date: str,
    occurrence: str = "first",
) -> int:
    """Replace one substring in one text node with a tracked delete/insert."""

    matches: list[tuple[etree._Element, int]] = []
    for node in text_nodes(paragraph):
        value = node.text or ""
        start = value.find(old)
        if start >= 0:
            matches.append((node, start))
    if not matches:
        paragraph_text = "".join(paragraph.itertext())
        raise RuntimeError(f"Could not find {old!r} in paragraph: {paragraph_text!r}")
    node, start = matches[-1] if occurrence == "last" else matches[0]
    value = node.text or ""
    end = start + len(old)
    template_run = node.getparent()
    if template_run is None or local_name(template_run.tag) != "r":
        raise RuntimeError("Text node is not inside a run")
    parent = template_run.getparent()
    if parent is None:
        raise RuntimeError("Run has no parent")
    index = parent.index(template_run)

    replacement_nodes: list[etree._Element] = []
    prefix = value[:start]
    suffix = value[end:]
    if prefix:
        replacement_nodes.append(make_run(template_run, prefix))
    replacement_nodes.append(make_revision(template_run, "del", revision_id, old, author, date))
    replacement_nodes.append(make_revision(template_run, "ins", revision_id + 1, new, author, date))
    if suffix:
        replacement_nodes.append(make_run(template_run, suffix))

    parent.remove(template_run)
    for offset, child in enumerate(replacement_nodes):
        parent.insert(index + offset, child)
    return revision_id + 2


def paragraph_map(root: etree._Element) -> list[etree._Element]:
    return root.xpath(".//w:body/w:p", namespaces=NS)


def main() -> None:
    if not SOURCE.exists():
        raise SystemExit(f"Missing source manuscript: {SOURCE}")
    if OUTPUT.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {OUTPUT}")

    with ZipFile(SOURCE) as archive:
        document_xml = archive.read("word/document.xml")
        settings_xml = archive.read("word/settings.xml")
        entries = {name: archive.read(name) for name in archive.namelist()}

    parser = etree.XMLParser(remove_blank_text=False)
    document_root = etree.fromstring(document_xml, parser)
    settings_root = etree.fromstring(settings_xml, parser)
    paragraphs = paragraph_map(document_root)
    if len(paragraphs) < 103 or "Methods" not in "".join(paragraphs[49].itertext()):
        raise SystemExit("Unexpected manuscript paragraph layout; refusing to patch")

    author = "gzy"
    date = "2026-09-05T00:00:00Z"
    revision_id = next_revision_id(document_root)
    operations: list[tuple[int, str, str, str]] = [
        (
            51,
            "30 Kla groups: nine non-tumor tissues, two tumor tissues",
            "31 Kla groups: nine non-tumor tissues, three tumor tissues",
            "Update the expanded 31-group category scope.",
        ),
        (
            52,
            "sample-level membership records",
            "group-level protein membership records",
            "Table S1 contains group-level membership, not the plotted sample rows.",
        ),
        (
            52,
            ".",
            ". Source-level sample observations used for the sample-level figures were retained in the analysis inputs and provenance records.",
            "Document the source-level figure inputs retained for the sample-only plots.",
        ),
        (
            69,
            "prioritized  the following",
            "prioritized in the following",
            "Correct the missing preposition in the reference-selection rule.",
        ),
        (
            71,
            "The final 30 Kla groups mapped to 28 unique whole-proteome display rows because two pairs of Kla studies used identical reference matrices and sample subsets; the analytical linkage remained at 30 rows.",
            "The final 31 Kla groups mapped to 29 unique whole-proteome display rows because two pairs of Kla studies used identical reference matrices and sample subsets; the analytical linkage remained at 31 rows. Reference reuse was retained for per-group linkage, flagged in provenance and not silently deduplicated in the descriptive plots.",
            "Update the expanded scope, deduplicated reference-display count and retained-reference provenance rule.",
        ),
        (
            73,
            "DNA DamageResponse",
            "DNA Damage Response",
            "Correct the DDR heading spacing.",
        ),
        (
            76,
            "For each sample group, the Kla and whole-proteome protein sets were intersected independently with the DDR set. The DDR fraction was calculated as",
            "For each publication group and each retained source-level sample observation, the Kla and whole-proteome protein sets were intersected independently with the DDR set. For the sample-level boxplots and pathway barplots, only observations labelled `ObservationType = \"sample\"` were retained; pooled, condition-level, model-level, dataset-union and other aggregate records were excluded from the plotted data and associated statistics. The DDR fraction was calculated as",
            "Specify the sample-only contract used for boxplot/barplot points and statistics.",
        ),
        (
            82,
            " summed and transformed as\u00a0",
            " summed; the summed values were transformed as\u00a0",
            "Make the within-sample log transformation explicit after feature summation.",
        ),
        (
            82,
            " ranked with average ranks",
            ", then ranked with average ranks",
            "Make the order of transformation and ranking explicit.",
        ),
        (
            83,
            "Regulator percentiles were summarized by the median across conditions or replicates in the same sample group. An undetected regulator was assigned 0 only as a visualization convention; this value does not denote measured biological zero abundance.",
            "Regulator percentiles were summarized by the median across conditions or replicates in the same sample group. For any sample in which a regulator was not detected, its percentile was set to 0 before the median as a visualization convention; this zero-coding affects the summary but does not denote measured biological zero abundance.",
            "Describe the per-sample zero-coding used before the group-level median.",
        ),
        (
            88,
            "Exact Four-Set Venn Analyses",
            "Exact Four-Set Intersection (UpSet) Analyses",
            "Rename the analysis to match the current UpSet figures.",
        ),
        (
            90,
            "The diagrams used fixed schematic geometry. Region area was not used to encode protein abundance or intersection size; printed counts were the quantitative encoding. Zero-count regions were retained, and all displayed counts were reconstructed from the underlying membership table rather than inferred from diagram area.",
            "The four-set intersections were visualized with UpSet plots: bars encode the exact protein count for each intersection, the dot matrix encodes set membership, and set-size bars show the marginal size of each category. All 15 mutually exclusive regions, including zero-count regions, were retained and reconstructed from the underlying membership table. For direct comparison, the whole-proteome DDR plot used the intersection order of the Kla-DDR plot.",
            "Describe the actual UpSet encoding rather than the former schematic Venn geometry.",
        ),
        (
            92,
            "399 unique",
            "401 unique",
            "Update the current Kla-DDR union size.",
        ),
        (
            92,
            "excluded from",
            "excluded from ",
            "Restore the missing space before the styled Table S5 reference.",
        ),
        (
            95,
            "The pathway rows were displayed from the lowest to the highest coefficient: BER, NER, MMR, FA, HR, AEJ and NHEJ.",
            "The score weights were assigned in the order BER, NER, MMR, FA, HR, AEJ and NHEJ; the rendered matrix rows were displayed in the figure-defined order BER, NER, MMR, FA, HR, NHEJ and AEJ.",
            "Distinguish score-weight order from the rendered matrix row order.",
        ),
        (
            96,
            "178",
            "192",
            "Update the tumor-tissue pathway matrix size.",
        ),
        (
            96,
            "292) .",
            "292).",
            "Remove the stray space before the sentence-ending period.",
        ),
        (
            98,
            "Data processing, statistical summaries, and visualization were performed in R 4.4.3 using the following packages: data.table 1.18.2.1, dplyr 1.2.0, readr 2.2.0, tidyr 1.3.2, ggplot2 4.0.2, ggVennDiagram 1.5.7, readxl 1.4.5, ragg 1.5.1 and digest 0.6.39. Dataset-specific evidence extraction also used Python 3.12.10. Figure and spreadsheet labels used Arial Unicode MS.",
            "Data processing, statistical summaries, and visualization were performed in R 4.4.3 using the following packages: data.table 1.18.2.1, dplyr 1.2.0, readr 2.2.0, tidyr 1.3.2, ggplot2 4.0.2, patchwork 1.3.2, readxl 1.4.5, openxlsx2 1.25, ragg 1.5.1, systemfonts 1.3.2 and digest 0.6.39. Dataset-specific evidence extraction also used Python 3.12.10. Figure and spreadsheet labels used Arial Unicode MS.",
            "Update the software list for the current UpSet/patchwork and workbook workflow.",
        ),
        (
            99,
            "Automated checks verified the 30-group scope, category counts of 9/2/12/7, 399 current Kla-DDR proteins, four pathway-matrix set sizes of 183/178/381/292, the 399-row filtered score table and exact reconstruction of all Venn regions.",
            "Automated checks verified the 31-group scope, category counts of 9/3/12/7, 401 current Kla-DDR proteins, four pathway-matrix set sizes of 183/192/381/292, the 401-row filtered score table and exact reconstruction of all 15 UpSet intersection regions.",
            "Update the automated QC statement to the current data contract and figure type.",
        ),
    ]

    applied: list[dict[str, object]] = []
    for index, old, new, reason in operations:
        occurrence = "last" if index == 52 and old == "." else "first"
        revision_id = replace_in_paragraph(
            paragraphs[index], old, new, revision_id, author, date, occurrence=occurrence
        )
        applied.append({"paragraph": f"P{index:04d}", "reason": reason, "old": old, "new": new})

    settings_changed = False
    if settings_root.find(w("trackRevisions")) is None:
        settings_root.insert(0, etree.Element(w("trackRevisions")))
        settings_changed = True

    entries["word/document.xml"] = etree.tostring(
        document_root, xml_declaration=True, encoding="UTF-8", standalone="yes"
    )
    # Preserve the source settings bytes when change tracking is already
    # enabled.  This keeps the patch limited to document.xml rather than
    # introducing a serializer-only diff in an unrelated package part.
    if settings_changed:
        entries["word/settings.xml"] = etree.tostring(
            settings_root, xml_declaration=True, encoding="UTF-8", standalone="yes"
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(OUTPUT, "w", compression=ZIP_DEFLATED) as archive:
        for name, data in entries.items():
            archive.writestr(name, data)

    print(f"Wrote {OUTPUT}")
    print(f"Applied {len(applied)} tracked replacements; next unused revision id: {revision_id}")
    for item in applied:
        print(f"{item['paragraph']}: {item['reason']}")


if __name__ == "__main__":
    main()
