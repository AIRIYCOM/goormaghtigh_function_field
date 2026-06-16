# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy"]
# ///
r"""
Monodromy verification  --  function field Goormaghtigh
========================================================================

Corroborates the "Monodromy of f_m" paragraph of the paper's numerical
section: for 3 <= m <= 30,

    Mon(f_m) = Gal( f_m(X) - t / \bar Q(t) ) = S_{m-1},     f_m = (X^m-1)/(X-1),

together with the two geometric inputs that the "critical values" Lemma
proves unconditionally:

    (S)  h_m = G_m/(X-1)^2  is separable          (critical points simple),
    (D)  the m-2 critical values are pairwise distinct,

where G_m = (m-1)X^m - m X^{m-1} + 1.

NOTE.  The paper's PROOF of Mon(f_m)=S_{m-1} is the elementary
transposition-tree argument (no computation).  The certification here is
INDEPENDENT corroborating evidence, exactly as the paper states.

Galois group of the generic fibre (Hilbert specialisation t -> t0):
  * exact:  PARI/GP `polgalois` for degree d = m-1 <= 11 (cypari2, optional),
            confirming order = d! (= S_d).
  * Frobenius cycle types (Dedekind), all m:  factor f_m(X)-t0 mod p over
    ~600 good primes; the factor-degree multiset is the Frobenius cycle type.
    We detect
       - a d-cycle      (one degree-d irreducible factor)  => transitive;
       - a transposition: a cycle type with exactly one even part, equal to 2,
         all other parts odd (raise such a Frobenius to lcm(odd parts) to leave
         a single 2-cycle = transposition);
       - an odd permutation                                  => group not in A_d.
    For PRIME d, a d-cycle gives primitivity (transitive of prime degree), and
    primitive + transposition => S_d by Jordan's theorem.  For composite d the
    same evidence is reported; the rigorous step is the transposition-tree Lemma.

Run (from the repository root):
  uv run src/monodromy_verification.py
"""
import sys
from math import factorial
from sympy import symbols, Poly, GF, resultant, gcd as sgcd

X, c = symbols("X c")

try:
    import cypari2
    _pari = cypari2.Pari()
    HAVE_PARI = True
except Exception:
    HAVE_PARI = False


def fm_expr(m):
    return sum(X**i for i in range(m))          # 1 + X + ... + X^{m-1}


def Gm(m):
    return Poly((m - 1) * X**m - m * X**(m - 1) + 1, X)


def hm_squarefree(m):
    """(S): h_m = G_m/(X-1)^2 is separable."""
    G = Gm(m)
    h, r = G.div(Poly((X - 1)**2, X))
    if not r.is_zero:
        return False
    return sgcd(h, h.diff(X)).degree() == 0


def Qm(m):
    """Critical-value polynomial Res_X(f_m(X)-c, f_m'(X)) with the spurious
    factor (c-m) stripped; its roots are the m-2 genuine critical values."""
    f = Poly(fm_expr(m), X)
    R = Poly(resultant(Poly(f.as_expr() - c, X), Poly(f.diff(X), X), X), c)
    while R.degree() >= 1:
        q, rem = R.div(Poly(c - m, c))
        if rem.is_zero:
            R = q
        else:
            break
    return R


def distinct_critical_values(m):
    """(D): critical values pairwise distinct <=> Q_m squarefree."""
    Q = Qm(m)
    return Q.degree() < 1 or sgcd(Q, Q.diff(c)).degree() == 0


def primes_upto(n):
    sieve = [True] * (n + 1)
    out = []
    for i in range(2, n + 1):
        if sieve[i]:
            out.append(i)
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return out


_PRIMES = primes_upto(9000)


def cycle_type_mod_p(m, t0, p):
    """Frobenius cycle type = factor-degree multiset of f_m(X)-t0 over GF(p),
    or None if ramified (non-separable) at p."""
    P = Poly(fm_expr(m) - t0, X, domain=GF(p))
    if P.degree() != m - 1 or sgcd(P, P.diff(X)).degree() != 0:
        return None
    typ = []
    for fac, e in P.factor_list()[1]:
        typ += [fac.degree()] * e
    return tuple(sorted(typ))


def yields_transposition(typ):
    """Literal transposition [2,1,...], or a type with exactly one even part
    (=2) and the rest odd (powers down to a single transposition)."""
    evens = [a for a in typ if a % 2 == 0]
    return len(evens) == 1 and evens[0] == 2


def is_odd_perm(typ):
    return sum(a - 1 for a in typ) % 2 == 1


def galois_evidence(m, n_primes=600):
    d = m - 1
    t0 = 100003
    primes = [p for p in _PRIMES if p > m][:n_primes]
    saw_dcycle = saw_transp = saw_odd = False
    n_types = set()
    for p in primes:
        ct = cycle_type_mod_p(m, t0, p)
        if ct is None:
            continue
        n_types.add(ct)
        if ct == (d,):
            saw_dcycle = True
        if yields_transposition(ct):
            saw_transp = True
        if is_odd_perm(ct):
            saw_odd = True
    d_prime = d >= 2 and all(d % k for k in range(2, int(d**0.5) + 1))
    certified = saw_dcycle and saw_transp and d_prime          # Jordan, prime d
    corroborated = saw_dcycle and saw_transp and saw_odd
    return saw_dcycle, saw_transp, saw_odd, d_prime, certified, corroborated, len(n_types)


def polgalois_order(m):
    """Order of Gal(f_m(X)-t0/Q) via PARI polgalois (degree d=m-1<=11), or None."""
    if not HAVE_PARI or m - 1 > 11:
        return None
    t0 = 100003
    coeffs = [1] * m
    coeffs[0] = 1 - t0                       # f_m - t0
    try:
        poly = _pari.Pol(list(reversed(coeffs)))
        return int(_pari.polgalois(poly)[0])
    except Exception:
        return None


def main():
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    print("=" * 78)
    print("MONODROMY VERIFICATION  --  Mon(f_m)=S_{m-1},  3 <= m <= 30")
    print("PARI polgalois: %s" % ("available" if HAVE_PARI
                                  else "unavailable (cycle-type evidence only)"))
    print("=" * 78)
    print("  m  d   (S)hsep (D)distinct  polgalois     dcyc transp odd   verdict")
    all_sd = 0
    all_sep_dist = True
    for m in range(3, 31):
        d = m - 1
        sep = hm_squarefree(m)
        dist = distinct_critical_values(m)
        order = polgalois_order(m)
        pg = "-" if order is None else ("%d/%d%s" % (order, factorial(d),
              "OK" if order == factorial(d) else "!!"))
        dc, tr, od, dp, cert, corr, _ = galois_evidence(m)
        exact_ok = (order == factorial(d)) if order is not None else False
        if cert or exact_ok:
            verdict = "S_%d (certified)" % d
            all_sd += 1
        elif corr:
            verdict = "S_%d (evidence)" % d
            all_sd += 1
        else:
            verdict = "?? CHECK"
        if not (sep and dist):
            all_sep_dist = False
        print("  %2d %2d   %-5s   %-5s      %-12s  %-4s %-5s %-4s  %s"
              % (m, d, sep, dist, pg, dc, tr, od, verdict))
    print("-" * 78)
    print("  (S)+(D) hold for all 3<=m<=30 : %s" % all_sep_dist)
    print("  Mon(f_m)=S_{m-1} certified/corroborated : %d/28" % all_sd)
    print("  (Composite-(m-1) rigour = transposition-tree Lemma, non-computational;")
    print("   cycle-type evidence d-cycle + transposition + odd element corroborates it.)")
    print("\nSUMMARY monodromy sep_distinct_all=%s S_d_all=%d/28 pari=%s"
          % (all_sep_dist, all_sd, HAVE_PARI))
    if not (all_sep_dist and all_sd == 28):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
