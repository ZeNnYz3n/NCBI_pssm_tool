"""
pssm_analyzer.py

Turns a ParsedPSSM into per-position conservation statistics and a
plain-language interpretation. All numbers here are computed directly
from the file's own freqRatios / scores -- nothing is fetched from
external databases (no Pfam/CDD/literature yet; see README for why
those are deliberately out of scope for v1).

Conservation math:
  freqRatio_i (from the PSSM) = observed_weighted_freq_i / background_freq_i
  => observed_weighted_freq_i = freqRatio_i * background_freq_i
  => information content at a position = sum_i p_i * log2(freqRatio_i)
     (this is exactly the relative-entropy/IC formula PSI-BLAST itself
     uses internally, computed here from the stored ratios rather than
     re-derived from scratch)

Background frequencies are the standard Robinson & Robinson amino-acid
frequencies BLAST uses as its default (BLOSUM62 background).
"""

import math
from pssm_parser import STANDARD_AA

BACKGROUND_FREQ = {
    "A": 0.0776, "R": 0.0510, "N": 0.0405, "D": 0.0546, "C": 0.0157,
    "Q": 0.0393, "E": 0.0672, "G": 0.0691, "H": 0.0227, "I": 0.0591,
    "L": 0.0943, "K": 0.0487, "M": 0.0212, "F": 0.0400, "P": 0.0523,
    "S": 0.0715, "T": 0.0596, "W": 0.0119, "Y": 0.0288, "V": 0.0632,
}


def information_content(freq_ratios_at_pos):
    """
    IC in bits for one position.

    Empirically (verified against the uploaded file), the values NCBI
    stores in `freqRatios` for this PSSM sum to ~1 across the 20 standard
    residues -- i.e. they are the weighted OBSERVED frequencies p_i
    themselves, not p_i/background_i despite the field's name. IC is
    therefore the standard relative-entropy formula:
        IC = sum_i p_i * log2(p_i / background_i)
    """
    ic = 0.0
    for aa in STANDARD_AA:
        p = freq_ratios_at_pos[aa]
        if p <= 0:
            continue
        ic += p * math.log2(p / BACKGROUND_FREQ[aa])
    return max(ic, 0.0)


def classify_conservation(ic_bits):
    if ic_bits >= 2.0:
        return "highly conserved"
    elif ic_bits >= 1.0:
        return "moderately conserved"
    elif ic_bits >= 0.5:
        return "variable"
    else:
        return "hypervariable"


FUNCTIONAL_HINTS = (
    "Positions this conserved are commonly associated with active sites, "
    "ligand/DNA binding surfaces, protein-protein interfaces, or residues "
    "required to maintain the structural fold -- not necessarily any one "
    "of these, but conservation this strong rarely happens by chance."
)


def analyze_position(pos_index, wild_residue, scores_at_pos, freq_at_pos,
                      favored_threshold=3, rejected_threshold=-3):
    """
    favored_threshold / rejected_threshold are in raw PSSM score units
    (roughly half-bits for PSI-BLAST-style matrices). Defaults of +3/-3
    are a reasonable starting point -- tune per-matrix if your scores
    run hotter or colder than typical BLOSUM62-derived profiles.
    """
    ic = information_content(freq_at_pos)
    conservation = classify_conservation(ic)

    ranked = sorted(scores_at_pos.items(), key=lambda kv: kv[1], reverse=True)
    best_aa, best_score = ranked[0]
    second_aa, second_score = ranked[1]
    worst_aa, worst_score = ranked[-1]
    wild_score = scores_at_pos[wild_residue]

    favored = [(aa, s) for aa, s in ranked if s >= favored_threshold]
    rejected = [(aa, s) for aa, s in ranked if s <= rejected_threshold]

    return {
        "position": pos_index + 1,  # 1-indexed for humans
        "wild_residue": wild_residue,
        "wild_score": wild_score,
        "information_content_bits": round(ic, 3),
        "conservation": conservation,
        "best_residue": best_aa,
        "best_score": best_score,
        "second_best_residue": second_aa,
        "second_best_score": second_score,
        "worst_residue": worst_aa,
        "worst_score": worst_score,
        "favored_residues": ";".join(f"{aa}({s})" for aa, s in favored),
        "rejected_residues": ";".join(f"{aa}({s})" for aa, s in rejected),
    }


def analyze_pssm(parsed):
    """Returns a list of per-position analysis dicts."""
    results = []
    for i, wild in enumerate(parsed.sequence):
        results.append(
            analyze_position(i, wild, parsed.scores[i], parsed.freq_ratios[i])
        )
    return results


def predict_mutation(parsed, position_1indexed, mutant_residue):
    """Module 7: 'what if position X becomes residue Y' predictor."""
    i = position_1indexed - 1
    if i < 0 or i >= parsed.num_columns:
        raise ValueError(f"Position {position_1indexed} out of range (1-{parsed.num_columns})")
    if mutant_residue not in STANDARD_AA:
        raise ValueError(f"'{mutant_residue}' is not a standard amino acid code")

    wild = parsed.sequence[i]
    wild_score = parsed.scores[i][wild]
    mut_score = parsed.scores[i][mutant_residue]
    ic = information_content(parsed.freq_ratios[i])
    conservation = classify_conservation(ic)

    delta = mut_score - wild_score
    if delta >= 0:
        verdict = "This substitution is favored or neutral by the profile -- likely tolerated."
    elif delta > -4:
        verdict = "Mildly disfavored by the profile; plausibly tolerated depending on context."
    elif delta > -8:
        verdict = "Disfavored by the profile; evolution rarely accepts this change here."
    else:
        verdict = "Strongly disfavored; this substitution is essentially never seen in homologs and is likely deleterious."

    return {
        "position": position_1indexed,
        "wild_residue": wild,
        "wild_score": wild_score,
        "mutant_residue": mutant_residue,
        "mutant_score": mut_score,
        "score_delta": delta,
        "position_conservation": conservation,
        "position_ic_bits": round(ic, 3),
        "verdict": verdict,
    }


def summary_narrative(parsed, results):
    n = len(results)
    high = [r for r in results if r["conservation"] == "highly conserved"]
    mod = [r for r in results if r["conservation"] == "moderately conserved"]
    pct_high = 100 * len(high) / n
    pct_mod = 100 * len(mod) / n

    top_positions = sorted(high, key=lambda r: r["information_content_bits"], reverse=True)[:8]
    top_str = ", ".join(f"{r['wild_residue']}{r['position']}" for r in top_positions)

    text = (
        f"The PSSM for {parsed.query_title or parsed.query_id} ({n} positions) shows "
        f"{pct_high:.0f}% of residues are highly conserved and a further {pct_mod:.0f}% "
        f"are moderately conserved, based on information content computed from the "
        f"profile's own frequency ratios. "
    )
    if top_positions:
        text += (
            f"The strongest conservation signals are at positions {top_str}, which are "
            f"good candidates for functionally or structurally important residues. "
        )
    text += FUNCTIONAL_HINTS
    return text
