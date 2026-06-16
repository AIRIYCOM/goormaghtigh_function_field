# /// script
# requires-python = ">=3.10"
# dependencies = ["mpmath"]
# ///
r"""
Eneström--Kakeya bound and forced-modulus check  --  function field Goormaghtigh
================================================================================

Numerical corroboration of the two load-bearing inequalities behind
Lemma (Eneström--Kakeya bound) and Lemma (cross-disjointness of critical
values) in the paper "On the function field analogue of the Goormaghtigh
equation".  These lemmas make the closed-form genus exact UNCONDITIONALLY,
so the numerics here are corroboration, not a logical input.

Two facts are checked, both with high-precision mpmath:

  (EK)   Every finite critical point rho of f_k(X) = (X^k-1)/(X-1) lies in
         the open unit disk |rho| < 1, hence every critical value satisfies
         |c| = k |rho|^{k-1} < k.
         [f_k'(X) = 1 + 2X + ... + (k-1)X^{k-2} has strictly increasing
          positive coefficients => Eneström--Kakeya => zeros in |z|<1.]

  (FM)   A common critical value of f_m and f_n would force
            |c| = L_{m,n}^{1/(m-n)}  with
            L_{m,n} = n^{n(m-1)} / m^{m(n-1)} * ((m-1)/(n-1))^{(m-1)(n-1)},
         and L_{m,n}^{1/(m-n)} > n for every m > n >= 3, contradicting (EK)
         applied with k = n (which gives |c| < n).  The contradiction is
         routed through the SMALLER index n.
         Equivalently log L - (m-n) log n = (n-1) G(m,n) with
            G(m,n) = -m log(m/n) + (m-1) log((m-1)/(n-1)),
         and dG/dn = (n-m)/(n(n-1)) < 0 with G(m,m^-) = 0, so G > 0.

Run:  uv run enestrom_kakeya_check.py
Expect: "ALL CHECKS PASSED" and zero failures.
"""
import mpmath as mp

mp.mp.dps = 40


def ek_radius(k):
    """Max modulus of a finite critical point of f_k (root of f_k')."""
    # f_k'(X) = sum_{j=1}^{k-1} j X^{j-1}; coefficients (k-1),...,2,1 (leading first)
    coeffs = [mp.mpf(j) for j in range(k - 1, 0, -1)]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    return max(abs(r) for r in roots), [k * (r ** (k - 1)) for r in roots]


def L(m, n):
    return (mp.mpf(n) ** (n * (m - 1)) / mp.mpf(m) ** (m * (n - 1))
            * (mp.mpf(m - 1) / mp.mpf(n - 1)) ** ((m - 1) * (n - 1)))


def G(m, n):
    return -m * mp.log(mp.mpf(m) / n) + (m - 1) * mp.log(mp.mpf(m - 1) / (n - 1))


def main():
    fails = 0

    # (EK) critical points in open unit disk; critical values |c| < k
    print("== Eneström--Kakeya: |rho| < 1 and max|c| < k ==")
    ek_kmax = 40   # paper checks to k=59; default kept modest for a ~1 min reproduce
    worst_rho = mp.mpf(0)
    for k in range(3, ek_kmax + 1):
        mr, cvs = ek_radius(k)
        worst_rho = max(worst_rho, mr)
        max_cv = max(abs(c) for c in cvs)
        if not (mr < 1):
            print(f"  FAIL k={k}: max|rho|={mr}")
            fails += 1
        if not (max_cv < k):
            print(f"  FAIL k={k}: max|c|={max_cv} !< {k}")
            fails += 1
    print(f"  k = 3..{ek_kmax}: max over all critical points |rho| = {float(worst_rho):.6f} (< 1)")

    # (FM) forced modulus L^{1/(m-n)} > n, and > max|critical value of f_n|
    print("== forced modulus: L^{1/(m-n)} > n  (contradiction via smaller index n) ==")
    mmax = 80   # paper checks to m=200 (separation) / m=399 (G>0); default modest for speed
    min_ratio = mp.mpf("1e9")
    argmin = None
    # precompute max|CV(f_n)| once per n
    cv_cache = {}
    for n in range(3, mmax):
        _, cvs = ek_radius(n)
        cv_cache[n] = max(abs(c) for c in cvs)
    for n in range(3, mmax):
        max_cv_n = cv_cache[n]
        for m in range(n + 1, mmax + 1):
            val = L(m, n) ** (mp.mpf(1) / (m - n))
            ratio = val / n
            if ratio < min_ratio:
                min_ratio = ratio
                argmin = (m, n)
            if not (val > n):
                print(f"  FAIL routing m={m},n={n}: L^(1/(m-n))={float(val):.6f} !> n={n}")
                fails += 1
            if not (val > max_cv_n):
                print(f"  FAIL sep m={m},n={n}: forced={float(val):.6f} <= max|CV(f_n)|={float(max_cv_n):.6f}")
                fails += 1
    print(f"  3<=n<m<={mmax}: min (L^(1/(m-n)) / n) = {float(min_ratio):.6f} at {argmin} (> 1)")

    # (G monotonicity) dG/dn < 0 and G(m, m^-) -> 0, hence G > 0
    print("== G(m,n) > 0  (monotone decreasing in n, limit 0 at n=m) ==")
    minG = mp.mpf("1e9")
    for m in range(4, mmax):
        for n in range(3, m):
            g = G(m, n)
            minG = min(minG, g)
            if not (g > 0):
                print(f"  FAIL G m={m},n={n}: G={float(g)}")
                fails += 1
    print(f"  3<=n<m<{mmax}: min G(m,n) = {float(minG):.6e} (> 0)")

    print()
    if fails == 0:
        print("ALL CHECKS PASSED (0 failures).")
    else:
        print(f"FAILURES: {fails}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
