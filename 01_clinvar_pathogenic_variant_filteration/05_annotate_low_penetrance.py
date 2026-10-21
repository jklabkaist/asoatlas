#!/usr/bin/env python3
"""
Annotate SCV-level ClinVar comments with a 'low_penetrance' TRUE/FALSE column
based on a curated set of keyword/regex patterns covering:
  - Direct incomplete/reduced/low/variable/partial penetrance phrasing
  - Non-penetrance / nonpenetrant
  - "penetrance is/was low|reduced|incomplete" and "penetrance of X is low|..."
  - "varying/variable clinical consequence(s)" (CFTR2 terminology)
  - "variable severity dependent/depending on..." (severity gated by other variants)
  - cis-modifier phrasing (effect influenced/modulated by an adjacent sequence element)
"""

import re
import sys

# Each entry: (label, regex). Labels are what we store in the keyword set file.
PATTERNS = [
    # Direct penetrance qualifiers (adjective/adverb + penetrance/penetrant)
    ("incomplete penetrance / incompletely penetrant",
     re.compile(r"\bincomplete(?:ly)?[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("reduced penetrance",
     re.compile(r"\breduced[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("decreased penetrance",
     re.compile(r"\bdecreased[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("low penetrance / low penetrant",
     re.compile(r"\blow[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("lower penetrance",
     re.compile(r"\blower[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("lowered penetrance",
     re.compile(r"\blowered[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("variable penetrance / variably penetrant",
     re.compile(r"\bvariabl[ye][-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("partial penetrance / partially penetrant",
     re.compile(r"\bpartial(?:ly)?[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("imperfect penetrance",
     re.compile(r"\bimperfect[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("diminished penetrance",
     re.compile(r"\bdiminished[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("less penetrant",
     re.compile(r"\bless[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("weakly penetrant",
     re.compile(r"\bweakly[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("age-related penetrance",
     re.compile(r"\bage[-\s]?related[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("age-dependent penetrance",
     re.compile(r"\bage[-\s]?dependent[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),
    ("incomplete clinical penetrance",
     re.compile(r"\bincomplete[-\s]+clinical[-\s]+penetra(?:nce|nt)\w*", re.IGNORECASE)),

    # Non-penetrance / nonpenetrant
    ("non-penetrance / nonpenetrant",
     re.compile(r"\bnon[-\s]?penetra(?:nce|nt)\w*", re.IGNORECASE)),

    # 'penetrance is/was low|reduced|incomplete|variable|partial|decreased|not complete|not 100%'
    ("penetrance is/was low|reduced|incomplete|variable|partial|decreased|not 100%",
     re.compile(
         r"\bpenetrance\s+(?:is|was)\s+(?:low|reduced|incomplete|variable|partial|decreased|not\s+complete|not\s+100\s*%?|<\s*\d+\s*%?)",
         re.IGNORECASE,
     )),

    # 'penetrance of X is/was low|reduced|incomplete' (X up to 100 chars, no period)
    ("penetrance of X is/was low|reduced|incomplete",
     re.compile(
         r"\bpenetrance\s+of\s+[^.]{1,100}?\s+(?:is|was)\s+(?:low|reduced|incomplete)\b",
         re.IGNORECASE,
     )),

    # CFTR2-style 'varying/variable clinical consequence(s)'
    ("varying clinical consequence(s)",
     re.compile(r"\bvarying\s+clinical\s+consequences?\b", re.IGNORECASE)),
    ("variable clinical consequence(s)",
     re.compile(r"\bvariable\s+clinical\s+consequences?\b", re.IGNORECASE)),

    # Variable severity gated by other variants
    ("variable severity dependent/depending on (presence of other variants)",
     re.compile(r"\bvariable\s+severity\s+depend(?:ing|ent)\b", re.IGNORECASE)),
    ("severity dependent/depending on the presence",
     re.compile(r"\bseverity\s+depend(?:ing|ent)\s+on\s+the\s+presence\b", re.IGNORECASE)),

    # General cis-modifier phrasing: a variant's effect is influenced/modulated/modified
    # by an adjacent sequence element. The concept (pathogenicity depends on co-occurring
    # cis variants → not 100% penetrant in isolation) applies to any gene, not just CFTR.
    ("(effect) influenced/modulated/modified by the adjacent (cis modifier)",
     re.compile(r"\b(?:influenced|modulated|modified)\s+by\s+the\s+adjacent\b", re.IGNORECASE)),
]


# Exclusion 1 — homozygous-state context (sentence-level check).
# Rationale: homozygous-state non-penetrance is already handled upstream via nhomalt;
# this annotation should only flag generic / in-trans penetrance descriptions.
HOMOZYG_RE = re.compile(r"\bhomozyg\w*", re.IGNORECASE)
SENT_BOUNDARY_RE = re.compile(r"(?<=[.!?])\s+")  # Split on sentence boundaries


def _get_sentence_containing(text: str, m_start: int, m_end: int) -> str:
    """Extract the sentence containing the match."""
    # Find the start of the sentence: look back for . ! ? followed by space
    sent_start = 0
    for m in SENT_BOUNDARY_RE.finditer(text[:m_start]):
        sent_start = m.end()

    # Find the end of the sentence: look forward for . ! ?
    sent_end = len(text)
    m = SENT_BOUNDARY_RE.search(text[m_end:])
    if m:
        sent_end = m_end + m.start()

    return text[sent_start:sent_end]


def _is_homo_specific(text: str, m_start: int, m_end: int) -> bool:
    """Return True if a penetrance-keyword match is in a sentence containing 'homozyg'."""
    sentence = _get_sentence_containing(text, m_start, m_end)
    return bool(HOMOZYG_RE.search(sentence))

# Exclusion 2 — tentative/uncertain phrasing.
# E.g. "albeit with potentially reduced penetrance", "may result in reduced penetrance"
# These express possibility rather than confirmed incomplete penetrance.
TENTATIVE_RE = re.compile(
    r"\b(?:potentially|may|might|possible|possibly|suggesting|suggest|could|suggest\s+a|apparent(?:ly)?|possible(?:ly)?|albeit)\b",
    re.IGNORECASE
)


def _is_tentative(text: str, m_start: int, m_end: int) -> bool:
    """Return True if penetrance phrase is preceded by tentative language in same sentence."""
    sentence = _get_sentence_containing(text, m_start, m_end)
    # Find the part of sentence before the match (within the same sentence)
    match_pos_in_sent = m_start - (text.find(sentence))
    before_match = sentence[:match_pos_in_sent]
    return bool(TENTATIVE_RE.search(before_match))

# Exclusion 3 — VCGS classification criteria format.
# E.g. "0112 - The condition associated with this gene has incomplete penetrance (PMID: 19372260). (I)"
# VCGS criteria describe general gene properties, not specific variant phenotypes.
VCGS_CRITERION_RE = re.compile(r"\b\d{4}\s*-\s*", re.IGNORECASE)  # "0112 - ..." anywhere in sentence


def _is_vcgs_criterion(text: str, m_start: int, m_end: int) -> bool:
    """Return True if penetrance phrase is part of a VCGS criterion (0XXX format)
    that describes gene/condition properties, not variant-specific phenotype.

    Check the immediate context (before and after the match) rather than whole sentence,
    since one sentence can mix gene-level and variant-level descriptions."""

    # Extract immediate context around the match (±100 chars, respecting sentence boundaries)
    back_start = max(0, m_start - 100)
    fwd_end = min(len(text), m_end + 100)

    # Respect sentence boundaries
    sent_start = text.rfind('.', 0, m_start)
    if sent_start == -1:
        sent_start = 0
    else:
        sent_start += 1  # After the period

    sent_end = text.find('.', m_end)
    if sent_end == -1:
        sent_end = len(text)

    # Use tighter context: from last period to next period
    context_start = sent_start
    context_end = sent_end
    context = text[context_start:context_end]

    # Check if this specific context (around the penetrance keyword) is VCGS criterion format
    # and describes gene/condition properties

    # Is there a VCGS criterion marker (0XXX - ) in this context?
    if not VCGS_CRITERION_RE.search(context):
        return False

    # If immediate context before penetrance mentions "this variant/allele", it's variant-level
    # Look back from the match to find what the subject is
    before_penetrance = text[context_start:m_start]
    if re.search(r"\bthis\s+(?:variant|allele|mutation|change|substitution)\b", before_penetrance, re.IGNORECASE):
        return False

    # If immediate context before penetrance mentions "observed in" or "reported in"
    # referring to patients/families, it's variant-level description
    if re.search(r"\b(?:observed|reported|identified|found)\s+(?:in|among)\s+(?:patients?|families|individuals)", before_penetrance, re.IGNORECASE):
        return False

    return True

# Exclusion 4 — molecular/cellular usage of "penetrant".
# E.g. "incompletely penetrant splicing defect" describes a partial splicing
# defect at the transcript level, not clinical incomplete penetrance.
MOLECULAR_NEXT_RE = re.compile(
    r"^\s*(?:splic(?:e|ing)|mRNA|exon|intron|transcript|isoform)\b",
    re.IGNORECASE,
)


def _is_molecular_specific(text: str, m_end: int) -> bool:
    """Return True if the match is immediately followed by a molecular/cellular noun
    (e.g. 'incompletely penetrant splicing defect' rather than clinical penetrance)."""
    fwd_hi = min(len(text), m_end + 30)
    return bool(MOLECULAR_NEXT_RE.match(text[m_end:fwd_hi]))


def match_keywords(text: str):
    """Return a sorted list of pattern labels that hit in `text`, excluding matches
    whose only context is the homozygous state, tentative phrasing, VCGS criteria,
    or molecular/cellular usage."""
    hits = []
    for label, pat in PATTERNS:
        for m in pat.finditer(text):
            if _is_homo_specific(text, m.start(), m.end()):
                continue
            if _is_tentative(text, m.start(), m.end()):
                continue
            if _is_vcgs_criterion(text, m.start(), m.end()):
                continue
            if _is_molecular_specific(text, m.end()):
                continue
            hits.append(label)
            break  # at least one valid match for this label is enough
    return hits


def main(infile, outfile, keywords_file):
    n_total = 0
    n_true = 0
    matched_labels = set()

    with open(infile, "r", encoding="utf-8", errors="replace") as fin, \
         open(outfile, "w", encoding="utf-8") as fout:

        header = fin.readline().rstrip("\n")
        fout.write(header + "\tlow_penetrance\n")

        for line in fin:
            # Don't strip - preserve embedded whitespace in last field; just drop the trailing newline.
            line_stripped = line.rstrip("\n")
            fields = line_stripped.split("\t")
            if len(fields) < 14:
                # Malformed line - emit as FALSE to keep alignment
                fout.write(line_stripped + "\tFALSE\n")
                n_total += 1
                continue

            comment = fields[13]
            hits = match_keywords(comment)
            flag = "TRUE" if hits else "FALSE"
            if hits:
                n_true += 1
                matched_labels.update(hits)
            fout.write(line_stripped + "\t" + flag + "\n")
            n_total += 1

            if n_total % 500000 == 0:
                print(f"  processed {n_total:,} rows, TRUE so far: {n_true:,}", file=sys.stderr)

    if keywords_file:
        with open(keywords_file, "w", encoding="utf-8") as kf:
            kf.write("# Keyword/regex patterns used to flag 'low_penetrance'\n")
            kf.write("# label, regex source (case-insensitive), and whether the pattern matched any submission.\n")
            kf.write("#\n")
            for label, pat in PATTERNS:
                used = "yes" if label in matched_labels else "no"
                kf.write(f"[{used}]\t{label}\t::\t{pat.pattern}\n")

    print(f"Done. Total rows: {n_total:,}; flagged TRUE: {n_true:,}", file=sys.stderr)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", action="store", required=True)
    parser.add_argument("-o", "--output", dest="output", action="store", default="")
    parser.add_argument("-k", "--keywords", dest="keywords", action="store", default="")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace(".tsv", ".lowP.tsv")

    main(infile, outfile, args.keywords)
