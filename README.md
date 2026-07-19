# PSSM Interpreter (v1)

Turns an NCBI PSSM ASN.1 text file (psiblast's ASCII ASN.1 Scoremat output)
into a conservation analysis: per-position information content, functional
interpretation, a mutation-effect predictor, and plots -- instead of a raw
matrix.

## Install

```bash
git clone https://github.com/<your-username>/pssm-interpreter.git
cd pssm-interpreter
pip install -r requirements.txt
```

## Quickstart

```bash
python3 run_analysis.py path/to/your_pssm.asn output_folder/
```

This writes `pssm_analysis.csv`, `pssm_analysis.json`, `conservation_profile.png`,
`score_heatmap.png`, and `report.md` into `output_folder/`. See
`example_output/` in this repo for what that looks like on a real file
(human PCNA).

## What's here

- `pssm_parser.py` -- parses the ASN.1 text format directly (no external
  ASN.1 library needed; this format is text, not the binary encoding).
  Extracts the query sequence, the score matrix, and the frequency data.
- `pssm_analyzer.py` -- computes information content (bits) per position
  from the profile's own data, classifies conservation, and can answer
  "what if position X mutates to Y?"
- `run_analysis.py` -- CLI that ties it together and writes a CSV, JSON,
  two plots, and a markdown report.
- `example_output/` -- a full run against your uploaded PCNA PSSM.

## Usage

```bash
python3 run_analysis.py path/to/your_pssm.asn output_folder/
```

Or use it programmatically:

```python
from pssm_parser import parse_pssm_asn_text
from pssm_analyzer import analyze_pssm, predict_mutation

parsed = parse_pssm_asn_text("your_pssm.asn")
results = analyze_pssm(parsed)

# "what if position 87 becomes Asp?"
print(predict_mutation(parsed, 87, "D"))
```

## Format notes (things that took real debugging to pin down)

1. **This only handles the text/ASCII ASN.1 form**, i.e. what you get from
   psiblast when it writes PSSM output in ASN.1 *value notation* (human
   readable, like your uploaded file). It does **not** parse the binary
   `.asn`/`.smp` encoding -- that needs NCBI's C++ toolkit or a real ASN.1
   BER/DER decoder, which is a separate (and much more annoying) piece of
   work if you need it later.
2. **Row order** is the 28-letter `ncbistdaa` alphabet: Gap, A, B, C, D, E,
   F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, X, Y, Z, U, `*`, O, J. This
   isn't documented anywhere obvious -- it was reverse-engineered by
   checking that self-match scores (wild-type residue vs. itself) came out
   highest, and that ambiguous codes (B/Z/U/O/J) carried the `-32768`
   sentinel for "disallowed."
3. **`freqRatios` is misleadingly named.** In this file the values sum to
   ~1 across the 20 standard residues, i.e. they're the weighted observed
   frequencies themselves, not `observed/background`. Information content
   is computed as the standard relative-entropy formula
   `IC = sum p_i * log2(p_i / background_i)` using standard BLOSUM62
   background frequencies. If you run this against a PSSM from a different
   psiblast version and the numbers stop summing to ~1, that assumption
   needs re-checking against the new file before trusting the IC values.


