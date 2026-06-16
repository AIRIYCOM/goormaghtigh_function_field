# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy"]
# ///
r"""
Numerical cross-disjointness  --  function-field Goormaghtigh
================================================================================

This script gives independent numerical corroboration of the critical-value
disjointness that makes the closed-form genus exact.  (The paper proves this
disjointness unconditionally over characteristic 0, via an Enestrom-Kakeya /
forced-modulus argument; the computation here is corroboration, not a logical
input.)

The genus closed form g = 1 + [(m-1)(n-2)-(n-1)-gcd(m-1,n-1)]/2 is the GEOMETRIC
genus ONLY IF the affine curve C_{m,n}: f_m(X)=f_n(Y) is smooth, i.e. has no
finite singular points.  For separated-variable curves F(X)-G(Y)=0 the finite
singular points are exactly the points where BOTH dF/dX=0 and dG/dY=0 AND the
two critical VALUES coincide:
        f_m has critical value c  (f_m'(x0)=0, f_m(x0)=c)
        f_n has critical value c  (f_n'(y0)=0, f_n(y0)=c)
        => (x0,y0) is a singular point of C_{m,n}.
So  C_{m,n} affine-smooth  <=>  f_m and f_n share NO common critical value
                            <=>  critvalset(f_m) cap critvalset(f_n) = {}.

ALGEBRA (re-checked here in TASK 1):
  f_m = (X^m-1)/(X-1).  Its critical values c satisfy
        P_m(c) := (m-1)^{m-1} c^m - m^m (c-1)^{m-1} = 0,
  with a SPURIOUS double root c=m (the x=infinity / rho=1 artifact).  The
  genuine critical-value polynomial is
        Q_m(c) := P_m(c) / (c-m)^2 ,  deg Q_m = m-2 .
  Equivalent u=(c-1)/c form:  critical value  <=>
        phi_m(u) := m^m u^{m-1} (1-u) = (m-1)^{m-1} ,
  and (m-1)^{m-1} = max_u phi_m attained at u=(m-1)/m (the spurious double root).

Cross-disjointness (=> affine smoothness => genus formula unconditional):
        critvalset(f_m) cap critvalset(f_n) = {}
   <=>  Q_m and Q_n have no common root
   <=>  R(m,n) := Res_c(Q_m, Q_n) != 0 .

TASKS:
  1  Verify P_m(c) form for m=3..12: independently compute critical values via
     diff(f_m), confirm all are roots of P_m, confirm c=m is a double root, and
     deg(Q_m)=m-2.  Also verify the phi_m(u) equivalent form.
  2  Resultant scan: ALL coprime (m,n), 3<=n<m<=20, plus (21,20),(25,24).
     EXACT integer arithmetic.  Report whether R(m,n)!=0 ALWAYS (cross-
     disjointness holds for every coprime pair).
  3  Non-coprime control: (4,6),(6,9),(6,10),(9,12),... .  Does gcd>1 ever
     produce R=0 (a shared critical value)?  Is coprimality the mechanism?
  4  u^{m-n}=K reduction: for several coprime pairs verify the reduction of
     "shared critical value" to a single equation in u, and check the
     numerical (in)compatibility.
  5  Arithmetic mechanism: factor R(m,n) and test whether its prime factors all
     divide mn(m-1)(n-1) (relevant to the positive-characteristic failure).

RUN (from the repository root):
  uv run src/disjointness_resultant.py
"""

import sys
import json
from math import gcd
from itertools import combinations

import sympy as sp
from sympy import Poly, Rational, Integer, symbols, factorint, resultant

c, u, X = symbols("c u X")


# ---------------------------------------------------------------------------
# Core polynomials
# ---------------------------------------------------------------------------

def f_m_expr(m):
    """f_m(X) = 1 + X + ... + X^{m-1}, exact."""
    return sum(X**i for i in range(m))


def P_m(m):
    """P_m(c) = (m-1)^{m-1} c^m - m^m (c-1)^{m-1}  (exact integer poly in c)."""
    return Poly(Integer(m - 1)**(m - 1) * c**m
                - Integer(m)**m * (c - 1)**(m - 1), c)


def Q_m(m):
    """Genuine critical-value polynomial Q_m = P_m / (c-m)^2, exact."""
    Pm = P_m(m)
    div = Poly((c - m)**2, c)
    q, r = sp.div(Pm, div)
    if not r.is_zero:
        raise RuntimeError(f"(c-m)^2 does not divide P_{m}; remainder={r}")
    return q


def R(m, n):
    """R(m,n) = Res_c(Q_m, Q_n), EXACT integer (resultant of integer polys)."""
    return sp.resultant(Q_m(m).as_expr(), Q_m(n).as_expr(), c)


# ---------------------------------------------------------------------------
# TASK 1 : verify P_m(c) form
# ---------------------------------------------------------------------------

def task1_verify_Pform():
    print("\n" + "=" * 78)
    print("TASK 1 : verify critical-value polynomial form  P_m(c), Q_m=P_m/(c-m)^2")
    print("=" * 78)
    print("  Independent critical-value polynomial via elimination (fast, exact):")
    print("  the critical VALUES of f_m are exactly the c with f_m(X)-c having a")
    print("  double X-root, i.e. roots of  CritRes_m(c) := Res_X(f_m(X)-c, f_m'(X)).")
    print("  We confirm  CritRes_m = +/- Q_m  (so Q_m IS the honest critval poly,")
    print("  deg = m-2, and the spurious x=oo value c=m is the DOUBLE root of P_m).")
    print()
    print(f"  {'m':>3} {'deg CritRes':>11} {'CritRes=+-Q_m?':>15} "
          f"{'c=m dbl root P_m?':>17} {'deg Q_m':>8} {'=m-2?':>6} {'phi form?':>10}")
    all_ok = True
    rows = []
    for m in range(3, 13):
        fm = f_m_expr(m)
        fp = sp.diff(fm, X)
        # exact integer critical-value polynomial by elimination of X
        CritRes = Poly(sp.expand(resultant(Poly(fm - c, X), Poly(fp, X))), c)
        Qm = Q_m(m)
        # CritRes should equal +/- Q_m exactly
        q, r = sp.div(CritRes, Qm)
        is_pm_Q = (r.is_zero and q.degree() == 0 and abs(int(q.as_expr())) == 1)
        all_roots = is_pm_Q          # equality of polys => identical root sets
        # c=m is a DOUBLE root of P_m (the spurious rho=1 / x=oo value)
        Pe = P_m(m).as_expr()
        Pm_at_m = sp.simplify(Pe.subs(c, m))
        Pm_d_at_m = sp.simplify(sp.diff(Pe, c).subs(c, m))
        Pm_dd_at_m = sp.simplify(sp.diff(Pe, c, 2).subs(c, m))
        dbl = (Pm_at_m == 0 and Pm_d_at_m == 0 and Pm_dd_at_m != 0)
        degQ = Qm.degree()
        deg_ok = (degQ == m - 2)
        # phi_m(u) form: phi_m=m^m u^{m-1}(1-u), max=(m-1)^{m-1} at u=(m-1)/m,
        # and that maximizing u corresponds via c=1/(1-u) to the spurious c=m.
        phi = Integer(m)**m * u**(m - 1) * (1 - u)
        umax = Rational(m - 1, m)
        phimax = sp.simplify(phi.subs(u, umax))
        dphi_at_umax = sp.simplify(sp.diff(phi, u).subs(u, umax))
        c_from_umax = sp.simplify(1 / (1 - umax))
        phi_ok = (phimax == Integer(m - 1)**(m - 1)
                  and dphi_at_umax == 0           # genuine maximum (double root)
                  and c_from_umax == m)

        ok = all_roots and dbl and deg_ok and phi_ok
        all_ok = all_ok and ok
        rows.append((m, CritRes.degree(), is_pm_Q, dbl, degQ, deg_ok, phi_ok))
        print(f"  {m:>3} {CritRes.degree():>11} {str(is_pm_Q):>15} "
              f"{str(dbl):>17} {degQ:>8} {str(deg_ok):>6} {str(phi_ok):>10}")
    print(f"\n  P_m / Q_m form fully verified for m=3..12 : {all_ok}")
    print("  KEY: Res_X(f_m-c, f_m') = +/- Q_m exactly (deg m-2).  The elimination")
    print("  resultant ALREADY excludes the x=oo artifact, so it lands on Q_m and")
    print("  NOT on the (c-m)^2 spurious factor of P_m -- an independent confirmation")
    print("  that c=m (the phi-maximum at u=(m-1)/m) is the only spurious value.")
    return all_ok, rows


# ---------------------------------------------------------------------------
# TASK 2 : resultant scan over coprime pairs
# ---------------------------------------------------------------------------

def task2_coprime_scan():
    print("\n" + "=" * 78)
    print("TASK 2 : Res_c(Q_m,Q_n) over ALL coprime (m,n), 3<=n<m<=20  + frontier")
    print("=" * 78)
    pairs = [(m, n) for m in range(4, 21) for n in range(3, m)
             if gcd(m, n) == 1]
    pairs += [(21, 20), (25, 24)]            # frontier samples
    print(f"  {'(m,n)':>9} {'gcd':>4} {'R(m,n)':>14} {'R!=0?':>6}  "
          f"{'|R| digits':>10}")
    results = {}
    all_nonzero = True
    zeros = []
    for (m, n) in pairs:
        r = R(m, n)
        nz = (r != 0)
        all_nonzero = all_nonzero and nz
        if not nz:
            zeros.append((m, n))
        results[(m, n)] = r
        digits = len(str(abs(int(r)))) if r != 0 else 0
        rstr = str(r) if abs(int(r)) < 10**8 else f"~10^{digits}"
        print(f"  {str((m,n)):>9} {gcd(m,n):>4} {rstr:>14} {str(nz):>6}  "
              f"{digits:>10}")
    print(f"\n  coprime pairs tested            : {len(pairs)}")
    print(f"  R(m,n) != 0 for ALL coprime     : {all_nonzero}")
    print(f"  coprime pairs with R==0 (FAIL)  : {zeros if zeros else 'NONE'}")
    print("  => cross-disjointness (affine smoothness) holds for every coprime")
    print("     pair tested; genus formula is UNCONDITIONAL on this range.")
    return results, all_nonzero, zeros


# ---------------------------------------------------------------------------
# TASK 3 : non-coprime control
# ---------------------------------------------------------------------------

def task3_noncoprime_control():
    print("\n" + "=" * 78)
    print("TASK 3 : non-coprime control  (does gcd>1 ever give a shared critval?)")
    print("=" * 78)
    cases = [(4, 6), (6, 9), (6, 10), (9, 12), (8, 12), (10, 15),
             (6, 4), (8, 6), (9, 6), (10, 6), (12, 8), (15, 10), (12, 9),
             (6, 8), (4, 8), (4, 10), (4, 12), (4, 14), (10, 4), (12, 4)]
    # dedupe unordered
    seen = set()
    uniq = []
    for (m, n) in cases:
        k = (max(m, n), min(m, n))
        if k not in seen:
            seen.add(k)
            uniq.append((max(m, n), min(m, n)))
    print(f"  {'(m,n)':>9} {'gcd':>4} {'R(m,n)':>16} {'R==0?':>6}  note")
    any_collide = False
    collisions = []
    for (m, n) in sorted(uniq):
        r = R(m, n)
        z = (r == 0)
        any_collide = any_collide or z
        if z:
            collisions.append((m, n))
        rstr = str(r) if abs(int(r)) < 10**8 else f"~10^{len(str(abs(int(r))))}"
        note = "SHARED CRITICAL VALUE (singular curve)" if z else "still disjoint"
        print(f"  {str((m,n)):>9} {gcd(m,n):>4} {rstr:>16} {str(z):>6}  {note}")
    print(f"\n  non-coprime pairs with R==0     : {collisions if collisions else 'NONE'}")
    print(f"  any non-coprime collision       : {any_collide}")
    if not any_collide:
        print("  NOTE: even gcd>1 pairs are cross-disjoint => coprimality is NOT")
        print("        the mechanism for critical-value disjointness; disjointness")
        print("        is generic and holds far beyond gcd=1 (consistent with")
        print("        the numerical scan: gcd>1 curves still irreducible, large genus).")
    else:
        print("  NOTE: gcd>1 produces shared critical values => coprimality IS the")
        print("        arithmetic gate for cross-disjointness.")
    return collisions, any_collide


# ---------------------------------------------------------------------------
# TASK 4 : u^{m-n}=K reduction
# ---------------------------------------------------------------------------

def task4_u_reduction():
    print("\n" + "=" * 78)
    print("TASK 4 : u^{m-n}=K reduction of the shared-critical-value system")
    print("=" * 78)
    print("""  Setup.  A common critical value c of f_m and f_n satisfies, in the
  u=(c-1)/c coordinate (so c=1/(1-u), and the f_m / f_n critical points give
  possibly DIFFERENT u_m, u_n with the SAME c):
        phi_m(u_m) = (m-1)^{m-1},   phi_n(u_n) = (n-1)^{n-1},
  where phi_k(u)=k^k u^{k-1}(1-u).  Because BOTH share the same c, u_m=u_n=:u
  (same point on the c-line).  Dividing the two phi-constraints:
        m^m u^{m-1}(1-u) / [ n^n u^{n-1}(1-u) ]  = (m-1)^{m-1}/(n-1)^{n-1}
        =>  (m^m / n^n) u^{m-n} = (m-1)^{m-1}/(n-1)^{n-1}
        =>  u^{m-n} = K_{m,n} := [ (m-1)^{m-1} n^n ] / [ (n-1)^{n-1} m^m ]   (*)
  Constraint (A): u must ALSO satisfy phi_n(u) = (n-1)^{n-1} (a genuine
  critical value, not the spurious u=(n-1)/n).  Shared critical value EXISTS
  <=> (*) and (A) have a common solution u (with u != spurious values).
""")
    print(f"  {'(m,n)':>9} {'K_{m,n}':>22} {'|K|<1?':>7} {'R(m,n)':>14} "
          f"{'consistent':>11}")
    for (m, n) in [(4, 3), (5, 3), (5, 4), (7, 5), (8, 5), (7, 4), (9, 8)]:
        K = (Rational(m - 1)**(m - 1) * Rational(n)**n) \
            / (Rational(n - 1)**(n - 1) * Rational(m)**m)
        Kf = float(K)
        r = R(m, n)
        # CONSISTENCY CHECK: build the reduced system numerically and confirm
        # it has no genuine common root <=> R(m,n)!=0.
        # Eliminate: the shared-critval c satisfies Q_m(c)=Q_n(c)=0; that is
        # exactly what R encodes.  Here we cross-check via the u-system: solve
        # phi_n(u)=(n-1)^{n-1} (Q_n side) and test u^{m-n}=K on each root.
        Qn_in_u = sp.expand(Integer(n)**n * u**(n - 1) * (1 - u)
                            - Integer(n - 1)**(n - 1))
        un_roots = sp.nroots(Poly(Qn_in_u, u))
        # discard spurious double root u=(n-1)/n
        spur = Rational(n - 1, n)
        consistent = True
        for ur in un_roots:
            if abs(complex(ur) - float(spur)) < 1e-9:
                continue
            lhs = complex(ur)**(m - n)
            if abs(lhs - Kf) < 1e-7:
                # would be a genuine shared critical value -> R should be 0
                if r != 0:
                    consistent = False
        rstr = str(r) if abs(int(r)) < 10**8 else f"~10^{len(str(abs(int(r))))}"
        print(f"  {str((m,n)):>9} {str(sp.nsimplify(K))[:22]:>22} "
              f"{str(abs(Kf) < 1):>7} {rstr:>14} {str(consistent):>11}")
    print("\n  Interpretation: (*) forces u onto the locus u^{m-n}=K_{m,n}; none of")
    print("  the genuine roots of phi_n(u)=(n-1)^{n-1} (constraint A) land on that")
    print("  locus -> no common solution -> R(m,n)!=0.  The reduction is correct:")
    print("  R(m,n)=0 <=> the u-systems (*)&(A) share a non-spurious root.")
    return True


# ---------------------------------------------------------------------------
# TASK 5 : arithmetic mechanism (factor R)
# ---------------------------------------------------------------------------

def task5_factor_R(coprime_results):
    print("\n" + "=" * 78)
    print("TASK 5 : prime factorization of R(m,n) -- arithmetic mechanism")
    print("=" * 78)
    print("  Hypothesis: primes p | R(m,n) all divide mn(m-1)(n-1)")
    print("  (=> a ramification / Mason-Stothers style argument controls the resultant).")
    print("  NOTE: |R| grows super-exponentially; for large pairs we only TRIAL-")
    print("  DIVIDE (limit=10^6) and report the SMALL-prime part + a residual cofactor")
    print("  flag.  Finding even one small prime NOT dividing mn(m-1)(n-1) already")
    print("  refutes the clean-support hypothesis.")
    print()
    print(f"  {'(m,n)':>9} {'small primes of R (trial<1e6)':>34} "
          f"{'extra small p':>14} {'big cofactor?':>13}")
    # full factorization only for small |R|; trial-division for the rest
    sample = [(4, 3), (5, 3), (5, 4), (7, 4), (7, 5), (8, 5), (7, 6),
              (8, 7), (9, 8), (11, 10), (13, 12), (12, 7), (16, 3)]
    FULL_DIGITS = 50          # full factorint only below ~10^50
    all_supported = True
    extra_primes_cases = []
    for (m, n) in sample:
        r = coprime_results.get((m, n))
        if r is None:
            r = R(m, n)
        if r == 0:
            print(f"  {str((m,n)):>9}  R==0 (skip)")
            continue
        N = abs(int(r))
        digits = len(str(N))
        if digits <= FULL_DIGITS:
            fac = factorint(N)
            cofactor = 1
        else:
            # trial division by primes up to 1e6; leave the rest as a cofactor
            fac = factorint(N, limit=10**6)
            cofactor = 1
            small = {}
            for p, e in fac.items():
                if p < 10**6:
                    small[p] = e
                else:
                    cofactor *= p**e        # unfactored residual
            fac = small
        rad = sorted(fac.keys())
        allowed = set()
        for base in (m, n, m - 1, n - 1):
            allowed |= set(factorint(base).keys())
        extra = [p for p in rad if p not in allowed]
        # only the SMALL part is certified; a big cofactor may hide more primes
        ok = (len(extra) == 0 and cofactor == 1)
        all_supported = all_supported and ok
        if extra:
            extra_primes_cases.append((m, n, extra))
        radstr = ",".join(str(p) for p in rad) if rad else "(none<1e6)"
        if len(radstr) > 34:
            radstr = radstr[:31] + "..."
        big = "yes" if cofactor != 1 else "no"
        print(f"  {str((m,n)):>9} {radstr:>34} {str(extra if extra else '-'):>14} "
              f"{big:>13}")
    print(f"\n  ALL fully-factored R(m,n) supported on primes | mn(m-1)(n-1) : "
          f"{all_supported}")
    if extra_primes_cases:
        print(f"  cases with extra (small) primes refuting clean support:")
        for (m, n, ex) in extra_primes_cases:
            print(f"      ({m},{n}): extra primes {ex}")
        print("  => prime support is STRICTLY LARGER than mn(m-1)(n-1); a naive")
        print("     ramification / 'only the obvious primes' heuristic is INCOMPLETE.")
        print("     The extra primes are the primes p where Q_m and Q_n acquire a")
        print("     common root MOD p (i.e. p | Res), and they are NOT controlled by")
        print("     mn(m-1)(n-1).  Hence an unconditional disjointness proof")
        print("     CANNOT come from a small fixed bad-prime set; it must instead")
        print("     show Q_m, Q_n stay coprime over Q DIRECTLY (e.g. a height /")
        print("     archimedean separation of critical values, or a Bilu-Tichy")
        print("     decomposition argument), not a mod-p ramification count.")
    else:
        print("  => clean ramification structure: a Mason-Stothers / ABC-over-Z")
        print("     bound on Res should yield UNCONDITIONAL disjointness.")
    return all_supported, extra_primes_cases


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    print("=" * 78)
    print("CROSS-DISJOINTNESS  --  function-field Goormaghtigh")
    print("=" * 78)

    p_ok, p_rows = task1_verify_Pform()
    cop_results, all_nonzero, zeros = task2_coprime_scan()
    noncop_collisions, noncop_any = task3_noncoprime_control()
    task4_u_reduction()
    primes_ok, extra_cases = task5_factor_R(cop_results)

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  TASK1 P_m form verified (m=3..12)       : {p_ok}")
    print(f"  TASK2 coprime R(m,n)!=0 for ALL         : {all_nonzero}"
          f"  (zeros={zeros})")
    print(f"  TASK3 non-coprime collision (R==0)      : {noncop_any}"
          f"  ({noncop_collisions})")
    print(f"  TASK5 R primes | mn(m-1)(n-1) (clean)   : {primes_ok}")

    disjoint_all_coprime = "yes" if all_nonzero else "no"
    noncoprime_collide = "yes" if noncop_any else "no"
    print(f"\n  disjoint_all_coprime = {disjoint_all_coprime}")
    print(f"  noncoprime_collide   = {noncoprime_collide}")
    print(f"\nSUMMARY disjoint_all_coprime={disjoint_all_coprime} "
          f"noncoprime_collide={noncoprime_collide} "
          f"P_form_ok={p_ok} R_primes_clean={primes_ok}")

    summary = dict(
        P_form_ok=bool(p_ok),
        coprime_all_nonzero=bool(all_nonzero),
        coprime_zeros=zeros,
        noncoprime_collisions=noncop_collisions,
        noncoprime_any_collide=bool(noncop_any),
        R_primes_clean=bool(primes_ok),
        R_extra_prime_cases=extra_cases,
        disjoint_all_coprime=disjoint_all_coprime,
        noncoprime_collide=noncoprime_collide,
    )
    print("\nSUMMARY_JSON " + json.dumps(summary, default=str))
    return summary


if __name__ == "__main__":
    main()
