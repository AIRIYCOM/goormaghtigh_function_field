# Changelog

## v1.0.2 — 2026-06-16

**Documentation / correctness patch. No change to any computed result; the
scripts produce the same numbers and the same `STATUS = ok` verdict as v1.0.1.**

This release brings the in-code documentation into agreement with the corrected,
published version of the paper, whose main theorem is now stated over
**characteristic 0** (an earlier draft over-claimed a "tame positive
characteristic" range that is in fact false — see below).

Changes:

- **Scope corrected to characteristic 0.** The header docstring of
  `src/genus_verification.py` and `README.md` previously described the theorem
  over a "tame field `k` (char 0, or char `p` with `p ∤ mn(m-1)(n-1)`)". The
  positive-characteristic statement is false: the critical-value disjointness
  underlying the proof is archimedean and can fail modulo `p`, even at tame
  primes. The smallest failure is `(m,n) = (3,4)` at `p = 19`, where
  `Res_c(Q_3, Q_4) = 19/16 ≡ 0 (mod 19)`, the curve acquires a node and becomes
  rational over `F̄_19`. The theorem is therefore a characteristic-0 statement.

- **Hypothesis corrected to `m ≠ n` (coprimality not needed).** Earlier text
  claimed `gcd(m,n) = 1` and that "gcd = 1 is necessary / gcd > 1 makes the curve
  reducible". This is false: for all `m ≠ n` with `min(m,n) ≥ 3` the curve is
  (absolutely) irreducible with genus ≥ 1, regardless of `gcd(m,n)`. The only
  reducible case is the diagonal `m = n` (the factor `X − Y`). TASK 2 is
  re-labelled a "gcd > 1 corroboration" accordingly.

- **TASK 5 re-framed.** The `F_q[t]` search is now described as an exploratory
  positive-characteristic probe that is *outside* the characteristic-0 theorem,
  and the false claim "tameness `p ∤ mn` is necessary" is removed. The genuine
  char-`p` obstruction is `p | Res_c(Q_m, Q_n)`.

- `is_absolutely_irreducible` docstring and the run's summary labels updated to
  match.

No change to `src/monodromy_verification.py`, `src/disjointness_resultant.py`,
or `src/enestrom_kakeya_check.py` (their outputs and claims were already
characteristic-0-correct).

## v1.0.1 — 2026-06-16

Initial archived release: four verification scripts (closed-form genus and
Riemann–Hurwitz cross-check, full symmetric monodromy of the repunit
polynomials, critical-value resultant disjointness, Eneström–Kakeya /
forced-modulus checks), reproducing every numerical claim of the paper.
