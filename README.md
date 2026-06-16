# Verification code — *On the function field analogue of the Goormaghtigh equation*

This repository reproduces every numerical claim in the paper **"On the
function field analogue of the Goormaghtigh equation"**. All results in the
paper are proved unconditionally; the scripts here are independent
**corroboration**, not a logical input to any proof.

The main theorem: for an algebraically closed field `k` of characteristic `0`,
`m ≠ n` and `min(m,n) ≥ 3`, the equation

```
(X^m − 1)/(X − 1) = (Y^n − 1)/(Y − 1)
```

has no non-constant solution `(X,Y) ∈ k(t)^2`. The proof computes the
geometric genus of the associated curve in closed form,

```
g(C_{m,n}) = 1 + ½ [ (m−1)(n−2) − (n−1) − gcd(m−1,n−1) ]  ≥ 1,
```

and excludes rational (ℙ¹) parametrisations by genus non-increase.
Coprimality of `(m,n)` is **not** required (only `m ≠ n`); the result is a
**characteristic-0** statement — in positive characteristic the underlying
critical-value disjointness can fail (even at tame primes), the smallest case
being `(m,n)=(3,4)` at `p=19`.

## Contents

| script | verifies | paper item |
|---|---|---|
| `src/genus_verification.py` | closed-form genus = Riemann–Hurwitz count (two independent projections agree); `g ≥ 1`; diagonal genera = triangular numbers (OEIS A000217); `min(m,n)<3` control and `gcd>1` corroboration (irreducible, `g≥1` — coprimality not needed); `disc_Y` squarefree; **and an exploratory positive-characteristic probe over `F_q[t]`** (`q∈{5,7,11}`, `m,n≤7`, `deg≤4`; char p is outside the char-0 theorem) | §Numerical (genus table, disc squarefree, char-p probe), Prop. (genus formula), Rem. (triangular) |
| `src/monodromy_verification.py` | `Mon(f_m) = S_{m−1}` for `3 ≤ m ≤ 30` (PARI `polgalois` for `m−1≤11` when cypari2 present, else Dedekind/Frobenius cycle types over ~600 primes: `d`-cycle + transposition + odd element); `h_m=G_m/(X−1)^2` separable; critical values pairwise distinct | §Numerical (monodromy), Lemma (monodromy), Lemma (critical values distinct) |
| `src/disjointness_resultant.py` | `Res_c(Q_m,Q_n) ≠ 0` (critical-value sets disjoint ⇒ affine smoothness) for all coprime and non-coprime pairs to `m ≤ 20`; `P_m(c)` form | Lemma (cross-disjointness), Lemma (smoothness) |
| `src/enestrom_kakeya_check.py` | Eneström–Kakeya bound `|ρ|<1`, `|c|<k`; forced modulus `L^{1/(m−n)} > n`; `G(m,n) > 0` | Lemma (Eneström–Kakeya), Lemma (cross-disjointness) |

`data/genus_table.csv` is the genus table reproduced by `genus_verification.py`.

## Requirements

- Python ≥ 3.10
- [`sympy`](https://www.sympy.org), [`mpmath`](https://mpmath.org)
- optional: [`cypari2`](https://github.com/sagemath/cypari2) (PARI/GP) for the
  exact Galois-group route in `monodromy_verification.py` (`m−1 ≤ 11`); the
  script falls back to Dedekind/Frobenius cycle types without it.

Each script also carries a [PEP 723](https://peps.python.org/pep-0723/)
inline dependency block, so with [`uv`](https://docs.astral.sh/uv/) no manual
install is needed.

## Reproduce

```bash
# with uv (recommended — resolves deps from the inline PEP 723 block)
bash reproduce_all.sh

# or run individually
uv run src/genus_verification.py
uv run src/monodromy_verification.py
uv run src/disjointness_resultant.py
uv run src/enestrom_kakeya_check.py

# or with a plain environment
pip install -r requirements.txt
python src/genus_verification.py
```

Total runtime is a few minutes on a single core. The ranges used match the
paper exactly: `enestrom_kakeya_check.py` checks the Eneström–Kakeya band for
`3 ≤ k ≤ 59` and the forced-modulus / `G(m,n)>0` inequalities for
`3 ≤ n < m ≤ 399`; `monodromy_verification.py` covers `3 ≤ m ≤ 30`;
`disjointness_resultant.py` the 101 coprime pairs (`3 ≤ n < m ≤ 20` plus
`(21,20),(25,24)`) and 11 non-coprime controls; `genus_verification.py` the
64 coprime pairs `3 ≤ n < m ≤ 14` plus samples `m ∈ {15,16,20}` and the
`F_q[t]` failure search (`q ∈ {5,7,11}`, `m,n ≤ 7`, `deg ≤ 4`).

## License

Code: MIT (see `LICENSE`). Please cite the paper and this archive if you use
the code; see `CITATION.cff`.
