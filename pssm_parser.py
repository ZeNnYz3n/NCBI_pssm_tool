"""
pssm_parser.py

Parses NCBI PSSM "Scoremat" files in ASN.1 *value notation* (text) form,
i.e. the human-readable .asn output of psiblast (-out_pssm_asn_ascii /
-save_pssm_after_last_round with text output), NOT the binary ASN.1
encoding. This is the format produced when psiblast is asked for an
ASCII/text ASN.1 PSSM rather than the binary .asn or .smp.

Only the fields needed for downstream analysis are extracted:
  - numRows / numColumns
  - query sequence + id/title
  - intermediateData.freqRatios  (weighted observed / background ratio)
  - finalData.scores             (integer PSSM log-odds scores)
  - finalData.lambda / kappa / h (Karlin-Altschul statistics)

Design note on layout: NCBI stores freqRatios/scores as a flat list that
is COLUMN-major over a fixed 28-letter alphabet (one block of 28 values
per sequence position, in the ncbistdaa order below). This was verified
against the uploaded file: position 0's scores block is
[-32768(Gap), -1(A), -32768(B), -2(C), ...] which only makes sense if
X (index 21) gets a small penalty, Z/U/O/J (rare/ambiguous) get -32768,
and '*' (stop) gets a fixed -4 -- exactly matching known PSI-BLAST
scoring conventions.
"""

import re
import json
from dataclasses import dataclass, field, asdict

# ncbistdaa alphabet order used by NCBI C++ toolkit PSSM matrices (28 letters)
NCBISTDAA_ORDER = [
    "Gap", "A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
    "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y", "Z", "U", "*", "O", "J",
]

# The 20 standard amino acids we care about for analysis
STANDARD_AA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


@dataclass
class ParsedPSSM:
    num_rows: int
    num_columns: int
    query_id: str
    query_title: str
    sequence: str
    scores: list          # list of dict{residue: score}, one per column, standard AA only
    freq_ratios: list      # list of dict{residue: ratio}, one per column, standard AA only
    scores_raw: list       # list of dict over full 28-letter alphabet (for edge cases)
    lambda_: float = None
    kappa: float = None
    h: float = None


def _to_float_triplet(mantissa, base, exponent):
    """NCBI encodes floats as {mantissa, base, exponent} => mantissa * base**exponent."""
    return mantissa * (base ** exponent)


def _extract_block(text, key):
    """Extract the balanced-brace contents of `key { ... }` starting after `key`."""
    m = re.search(re.escape(key) + r"\s*{", text)
    if not m:
        return None
    start = m.end() - 1  # position of the opening '{'
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    raise ValueError(f"Unbalanced braces while extracting block '{key}'")


def _parse_triplet_list(block_text):
    """Parse a sequence of `{ m, b, e }` triplets into floats."""
    triplets = re.findall(r"\{\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\}", block_text)
    return [_to_float_triplet(int(m_), int(b_), int(e_)) for m_, b_, e_ in triplets]


def _parse_int_list(block_text):
    """Parse a flat comma-separated list of integers (skips any stray braces)."""
    nums = re.findall(r"-?\d+", block_text)
    return [int(n) for n in nums]


def parse_pssm_asn_text(path: str) -> ParsedPSSM:
    with open(path, "r", errors="replace") as f:
        text = f.read()

    num_rows = int(re.search(r"numRows\s+(\d+)", text).group(1))
    num_columns = int(re.search(r"numColumns\s+(\d+)", text).group(1))

    # --- query id / title ---
    id_match = re.search(r'local str\s+"([^"]+)"', text)
    query_id = id_match.group(1) if id_match else "unknown"

    title_match = re.search(r'title\s+"([^"]*(?:\n[^"]*)*?)"', text)
    query_title = re.sub(r"\s+", " ", title_match.group(1)).strip() if title_match else ""

    # --- sequence (ncbieaa "....") ---
    seq_match = re.search(r'seq-data ncbieaa\s+"([^"]*(?:\n[^"]*)*?)"', text)
    if not seq_match:
        raise ValueError("Could not find seq-data ncbieaa sequence block")
    sequence = re.sub(r"\s+", "", seq_match.group(1))

    if len(sequence) != num_columns:
        raise ValueError(
            f"Sequence length ({len(sequence)}) != numColumns ({num_columns}); "
            "file may use a different layout than expected."
        )

    # --- freqRatios (intermediateData) ---
    inter_block = _extract_block(text, "intermediateData")
    freq_block = _extract_block(inter_block, "freqRatios")
    freq_flat = _parse_triplet_list(freq_block)

    # --- scores (finalData) ---
    final_block = _extract_block(text, "finalData")
    scores_block = _extract_block(final_block, "scores")
    scores_flat = _parse_int_list(scores_block)

    expected = num_rows * num_columns
    if len(scores_flat) != expected:
        raise ValueError(f"Expected {expected} scores, got {len(scores_flat)}")
    if len(freq_flat) != expected:
        raise ValueError(f"Expected {expected} freqRatios, got {len(freq_flat)}")

    # --- lambda / kappa / h ---
    def _scalar_triplet(name):
        blk = _extract_block(final_block, name)
        vals = _parse_triplet_list("{" + blk + "}") if blk else []
        return vals[0] if vals else None

    lambda_ = _scalar_triplet("lambda")
    kappa = _scalar_triplet("kappa")
    h = _scalar_triplet("h")

    # --- reshape column-major: block of num_rows values per sequence position ---
    scores_per_col = []
    freq_per_col = []
    scores_raw_per_col = []
    for col in range(num_columns):
        block_scores = scores_flat[col * num_rows:(col + 1) * num_rows]
        block_freqs = freq_flat[col * num_rows:(col + 1) * num_rows]
        row_scores = dict(zip(NCBISTDAA_ORDER, block_scores))
        row_freqs = dict(zip(NCBISTDAA_ORDER, block_freqs))
        scores_raw_per_col.append(row_scores)
        scores_per_col.append({aa: row_scores[aa] for aa in STANDARD_AA})
        freq_per_col.append({aa: row_freqs[aa] for aa in STANDARD_AA})

    return ParsedPSSM(
        num_rows=num_rows,
        num_columns=num_columns,
        query_id=query_id,
        query_title=query_title,
        sequence=sequence,
        scores=scores_per_col,
        freq_ratios=freq_per_col,
        scores_raw=scores_raw_per_col,
        lambda_=lambda_,
        kappa=kappa,
        h=h,
    )


if __name__ == "__main__":
    import sys
    p = parse_pssm_asn_text(sys.argv[1])
    print(f"Query: {p.query_id}")
    print(f"Title: {p.query_title}")
    print(f"Length: {p.num_columns}, Alphabet size: {p.num_rows}")
    print(f"Lambda={p.lambda_}, Kappa={p.kappa}, H={p.h}")
    print("First 3 positions (scores):")
    for i in range(3):
        print(p.sequence[i], p.scores[i])
