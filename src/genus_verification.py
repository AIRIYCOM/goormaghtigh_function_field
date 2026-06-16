# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy"]
# ///
r"""
L3 NUMERICAL GENUS — Proposition P1 (function-field Goormaghtigh theorem)
========================================================================

Proposition P1.  Over a tame field k (char 0, or char p with p \nmid mn),
gcd(m,n)=1 and min(m,n) >= 3, the equation

        (X^m - 1)/(X - 1) = (Y^n - 1)/(Y - 1)

has NO non-constant solution (X, Y) in \bar k(t).

Core lemma:  g(\widetilde C_{m,n}) >= 1, where \widetilde C_{m,n} is the
smooth projective model of the affine curve

        C_{m,n} :  f_m(X) = g_n(Y),
        f_m(X) = 1 + X + ... + X^{m-1}  (degree m-1),
        g_n(Y) = 1 + Y + ... + Y^{n-1}  (degree n-1).

If g >= 1 then P^1 admits NO non-constant map to \widetilde C_{m,n}
(Riemann-Hurwitz: a non-constant map P^1 -> C is dominant, forcing
0 = g(P^1) >= g(C) >= 1, contradiction).  Hence every \bar k(t)-point is
constant.  This is the clean, unconditional core of the proof.

This L3 script PUSHES the genus evidence to the computational frontier and
adds four cross-checks:

  TASK 1  Multi-scale genus scan over all gcd(m,n)=1, 3 <= n < m <= 14,
          plus samples m in {15,16,20}.  Verify  gcd=1 & min>=3 => g>=1 with
          NO exception.  CSV -> data/L3num_genus.csv.  Genus by Riemann-
          Hurwitz on pi_X : C -> P^1_X (deg n-1), cross-checked against the
          symmetric pi_Y projection (genus is symmetric in m,n).

  TASK 2  gcd>1 control: (4,6),(6,9),(4,8),(6,10) — genus + reducibility test
          (verifies gcd=1 is NECESSARY; non-coprime curves can be reducible /
          have rational components / lower the effective obstruction).

  TASK 3  min<3 control: (2,3),(2,5),(2,7) — should give g=0 (verifies
          min>=3 is NECESSARY: m=2 makes f_2(X)=1+X linear, so X is a
          polynomial in Y, the "curve" is rational and HAS non-constant
          solutions).

  TASK 4  OEIS lookup of the m=3 genus sequence and closed-form genus fits.

  TASK 5  failure search: brute-force F_q[t] (tame q in {5,7,11}, deg<=4,
          m,n<=7) for non-constant solutions of f_m(X(t))=g_n(Y(t)) (expect
          NONE), and a WILD-characteristic probe (p | mn) to see whether
          non-constant solutions appear once we leave the tame regime.

Genus-at-infinity formula used:
    R_inf = (n-1) - gcd(m-1, n-1)      [Newton polygon of f_m(X)=g_n(Y) at oo]
Derivation: near infinity f_m(X) ~ X^{m-1}, g_n(Y) ~ Y^{n-1}, so the places
over X=oo are governed by X^{m-1} = Y^{n-1}, which has gcd(m-1,n-1) places
each of ramification index (n-1)/gcd(m-1,n-1) for the pi_X cover; hence
R_inf = gcd*( (n-1)/gcd - 1 ) = (n-1) - gcd(m-1,n-1).
This is validated INTERNALLY by the two-projection genus agreement
(g via pi_X == g via pi_Y) for every case — a strong independent check, since
the two projections use different degrees and different R_inf values.
"""

import csv
import os
import sys
import json
import subprocess
from math import gcd
from itertools import product

from sympy import symbols, Poly, discriminant, factor_list
from sympy.abc import X, Y

OUT_CSV = None  # set in main from argv


# ---------------------------------------------------------------------------
# Genus via Riemann-Hurwitz
# ---------------------------------------------------------------------------

def f_poly(var, m):
    """1 + var + ... + var^{m-1} = (var^m-1)/(var-1), exact polynomial."""
    return Poly(sum(var**i for i in range(m)), var)


def genus_projection(m, n):
    """
    Geometric genus of the normalization of f_m(X)=g_n(Y) via
    pi_X : C -> P^1_X, (X,Y) |-> X, of degree d = n-1.

        2 g - 2 = -2 d + R_fin + R_inf
        R_fin   = deg_X Disc_Y( g_n(Y) - f_m(X) )      (exact for tame covers:
                  ord_x(Disc) = sum_i (e_i - 1) = local ramification)
        R_inf   = (n-1) - gcd(m-1, n-1)                (Newton polygon at oo)

    Returns (genus_or_None, detail_dict).
    """
    fX = f_poly(X, m)
    gY = f_poly(Y, n)
    d = n - 1
    PXY = Poly(gY.as_expr() - fX.as_expr(), Y)
    discX = Poly(discriminant(PXY, Y), X)
    R_fin = discX.degree()
    gg = gcd(m - 1, n - 1)
    R_inf = (n - 1) - gg
    two_g_minus_2 = -2 * d + R_fin + R_inf
    detail = dict(d=d, R_fin=R_fin, R_inf=R_inf, two_g_minus_2=two_g_minus_2,
                  disc_sqfree=_is_squarefree(discX))
    if two_g_minus_2 % 2 != 0:
        return None, detail
    return (two_g_minus_2 + 2) // 2, detail


def _is_squarefree(poly):
    """True iff a univariate sympy Poly is squarefree (rad deg == deg)."""
    if poly.degree() <= 0:
        return True
    g = poly.gcd(poly.diff())
    return g.degree() == 0


def genus_both(m, n):
    """Genus via pi_X and via pi_Y (swap roles).  Returns (gX, gY, dX, dY)."""
    gX, dX = genus_projection(m, n)
    gY, dY = genus_projection(n, m)   # symmetric projection to P^1_Y
    return gX, gY, dX, dY


def is_absolutely_irreducible(m, n):
    """
    Heuristic absolute-irreducibility check of F = g_n(Y) - f_m(X) over Q:
    return (irreducible_over_Q, n_factors).  Absolute irreducibility is
    stronger, but for separated-variable curves with gcd(m,n)=1 the Q-factor
    count of 1 together with the coprime branch profile at infinity certifies
    geometric connectedness (advocate angle 2).  We report the Q-factor count.
    """
    F = (f_poly(X, m).as_expr() - f_poly(Y, n).as_expr())
    c, facs = factor_list(F)
    nfac = sum(mult for _, mult in facs)
    return (nfac == 1), nfac


# ---------------------------------------------------------------------------
# TASK 1 : multi-scale genus scan
# ---------------------------------------------------------------------------

def task1_scan():
    print("\n" + "=" * 74)
    print("TASK 1 : multi-scale genus scan  (gcd=1, 3<=n<m<=14  +  m in 15,16,20)")
    print("=" * 74)

    pairs = []
    for m in range(3, 15):
        for n in range(3, m):
            pairs.append((m, n))
    # samples at the frontier
    for m in (15, 16, 20):
        for n in range(3, m):
            if gcd(m, n) == 1:           # sample only coprime at the frontier
                pairs.append((m, n))

    rows = []                 # full rows for CSV (both orders, for symmetry CSV)
    coprime_min3 = []         # (m,n,g) with gcd=1, min>=3
    sym_fail = []
    print(f"{'(m,n)':>8} {'gcd':>4} {'g':>5} {'>=1':>4} {'symOK':>6} "
          f"{'discSF':>7}  {'Rfin':>5} {'Rinf':>5}")
    for (m, n) in pairs:
        gX, gY, dX, dY = genus_both(m, n)
        g_val = gX if gX is not None else gY
        sym = (gX == gY)
        if not sym:
            sym_fail.append((m, n, gX, gY))
        gg = gcd(m, n)
        rows.append((m, n, gg, g_val, dX['R_fin'], dX['R_inf'],
                     dX['two_g_minus_2'], int(dX['disc_sqfree'])))
        if gg == 1 and min(m, n) >= 3:
            coprime_min3.append((m, n, g_val))
            ge1 = (g_val is not None and g_val >= 1)
            flag = "" if ge1 else "  <<< VIOLATION g<1 !!!"
            print(f"{str((m,n)):>8} {gg:>4} {str(g_val):>5} {str(ge1):>4} "
                  f"{str(sym):>6} {str(dX['disc_sqfree']):>7}  "
                  f"{dX['R_fin']:>5} {dX['R_inf']:>5}{flag}")

    # write CSV: canonical m,n,gcd(m-1,n-1),genus  (one row per unordered pair, m>n).
    # The gcd column is gcd(m-1,n-1) -- the quantity that enters the genus
    # formula and the paper's Table 1 -- NOT gcd(m,n).
    seen = set()
    csv_rows = []
    for (m, n, gg, g_val, rf, ri, t2, sf) in rows:
        key = (max(m, n), min(m, n))
        if key in seen:
            continue
        seen.add(key)
        M, N = max(m, n), min(m, n)
        csv_rows.append((M, N, gcd(M - 1, N - 1), g_val))
    csv_rows.sort()
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["m", "n", "gcd(m-1,n-1)", "genus"])
        for r in csv_rows:
            w.writerow(r)
    print(f"\n  CSV written: {OUT_CSV}  ({len(csv_rows)} unordered pairs)")

    # verdict
    viol = [(m, n, g) for (m, n, g) in coprime_min3
            if not (g is not None and g >= 1)]
    genus_min = min((g for (_, _, g) in coprime_min3 if g is not None),
                    default=None)
    print(f"\n  coprime & min>=3 cases tested : {len(coprime_min3)}")
    print(f"  genus >= 1 violations         : {len(viol)}  {viol if viol else ''}")
    print(f"  MIN genus over coprime&min>=3 : {genus_min}")
    print(f"  two-projection symmetry fails : {len(sym_fail)}  "
          f"{sym_fail if sym_fail else '(none -> R_inf formula validated)'}")
    # which cases are exactly g=1 (the elliptic boundary, g>=2 would be false)
    g1 = sorted((m, n) for (m, n, g) in coprime_min3 if g == 1)
    print(f"  cases with g == 1 (g>=2 FALSE): {g1}")
    return dict(coprime_min3=coprime_min3, genus_min=genus_min,
                violations=viol, sym_fail=sym_fail, all_rows=rows)


# ---------------------------------------------------------------------------
# TASK 2 : gcd>1 control
# ---------------------------------------------------------------------------

def task2_gcd_control():
    print("\n" + "=" * 74)
    print("TASK 2 : gcd>1 control  (verifies gcd=1 is necessary)")
    print("=" * 74)
    cases = [(4, 6), (6, 9), (4, 8), (6, 10)]
    print(f"{'(m,n)':>8} {'gcd':>4} {'g':>5} {'irred/Q':>8} {'#Qfac':>6}  note")
    results = []
    for (m, n) in cases:
        gX, gY, dX, dY = genus_both(m, n)
        g_val = gX if gX is not None else gY
        irr, nfac = is_absolutely_irreducible(m, n)
        note = "irreducible" if irr else f"REDUCIBLE ({nfac} Q-factors)"
        results.append((m, n, gcd(m, n), g_val, irr, nfac))
        print(f"{str((m,n)):>8} {gcd(m,n):>4} {str(g_val):>5} {str(irr):>8} "
              f"{nfac:>6}  {note}")
    any_reducible = any(not irr for (_, _, _, _, irr, _) in results)
    print(f"\n  any gcd>1 case reducible over Q : {any_reducible}")
    print("  (reducibility / shared factor (X-1)(Y-1) structure shows why the")
    print("   coprimality hypothesis is genuinely used.)")
    return results


# ---------------------------------------------------------------------------
# TASK 3 : min<3 control
# ---------------------------------------------------------------------------

def task3_min_control():
    print("\n" + "=" * 74)
    print("TASK 3 : min<3 control  (verifies min>=3 is necessary)")
    print("=" * 74)
    cases = [(2, 3), (2, 5), (2, 7)]
    print(f"{'(m,n)':>8} {'gcd':>4} {'g':>5} {'g==0?':>6}  note")
    results = []
    allzero = True
    for (m, n) in cases:
        # m=2: f_2(X)=1+X is linear -> X = g_n(Y)-1, a rational parametrization
        gX, gY, dX, dY = genus_both(m, n)
        g_val = gX if gX is not None else gY
        is0 = (g_val == 0)
        allzero = allzero and is0
        note = "rational (X is a poly in Y) -> non-constant solns EXIST" if is0 \
               else "unexpected"
        results.append((m, n, gcd(m, n), g_val))
        print(f"{str((m,n)):>8} {gcd(m,n):>4} {str(g_val):>5} {str(is0):>6}  {note}")
    print(f"\n  all min<3 cases have g==0 : {allzero}")
    print("  (m=2 makes f_2 linear so X is a polynomial in Y: the curve is")
    print("   rational and DOES carry non-constant k(t)-solutions; this is")
    print("   exactly why min>=3 is required.)")
    return results, allzero


# ---------------------------------------------------------------------------
# TASK 4 : OEIS lookup + genus formula fit
# ---------------------------------------------------------------------------

def oeis_query(seq):
    """Query OEIS in text mode (NOT arxiv -> no rate limit). Returns A-numbers
    + names found, or a status string. Uses curl per task spec."""
    q = ",".join(str(x) for x in seq)
    try:
        out = subprocess.run(
            ["curl", "-sG", "https://oeis.org/search?fmt=text",
             "--data-urlencode", f"q={q}"],
            capture_output=True, text=True, timeout=30).stdout
    except Exception as e:
        return f"(curl failed: {e})", []
    hits = []
    cur_id = None
    for line in out.splitlines():
        if line.startswith("%I "):
            parts = line.split()
            cur_id = parts[1] if len(parts) > 1 else "?"
        elif line.startswith("%N ") and cur_id:
            name = line[3:].split(" ", 1)[-1] if " " in line[3:] else line[3:]
            # strip leading A-number token
            toks = line.split(None, 2)
            name = toks[2] if len(toks) > 2 else line
            hits.append((cur_id, name))
            cur_id = None
    nres = "?"
    for line in out.splitlines():
        if "Search:" in line or "results" in line.lower():
            nres = line.strip()
            break
    return (out[:200] if not hits else f"{len(hits)} hit(s)"), hits


def task4_oeis(scan):
    print("\n" + "=" * 74)
    print("TASK 4 : OEIS lookup + genus formula fit")
    print("=" * 74)
    gmap = {}
    rfin_map = {}
    for (m, n, gg, g_val, rf, ri, t2, sf) in scan['all_rows']:
        gmap[(m, n)] = g_val
        gmap[(n, m)] = g_val          # genus is symmetric
        rfin_map[(m, n)] = rf         # R_fin for the pi_X projection of (m,n)

    # m=3 sequence: n=4,5,7,8,10,11,13,14 (coprime to 3, >=4) — genus values
    ns_m3 = [4, 5, 7, 8, 10, 11, 13, 14]
    seq_m3 = [gmap.get((3, n)) for n in ns_m3]
    print(f"  m=3, n in {ns_m3}:  genus = {seq_m3}")

    # also the dense m=3 coprime sequence n=4,5,7,8,10,11,13,14,16,17 (formula floor((n-1)/?)...)
    queries = {
        "m=3 (n=4,5,7,8,10,11,13,14)": [g for g in seq_m3 if g is not None],
    }
    # generic small coprime diagonal genus g(n,n+1): (3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10)
    diag = [gmap.get((n + 1, n)) for n in range(3, 10)]
    diag = [g for g in diag if g is not None]
    queries["diagonal g(n+1,n) n=3..9"] = diag

    oeis_results = {}
    for label, seq in queries.items():
        if len(seq) < 4:
            print(f"\n  [{label}] sequence too short, skip OEIS")
            continue
        print(f"\n  [{label}] querying OEIS with: {seq}")
        status, hits = oeis_query(seq)
        oeis_results[label] = (seq, status, hits)
        print(f"    -> {status}")
        for aid, name in hits[:5]:
            print(f"       {aid}: {name[:90]}")

    # closed-form genus fits on coprime, min>=3 (m<n)
    print("\n  ## closed-form genus fits (coprime, min>=3):")
    cop = sorted((min(m, n), max(m, n), g)
                 for (m, n, g) in scan['coprime_min3'] if g is not None)
    cop = sorted(set(cop))
    cands = {
        "(m-2)(n-2)/2  [floor]":      lambda a, b: ((a - 2) * (b - 2)) // 2,
        "ceil((m-2)(n-2)/2)":         lambda a, b: -((-(a - 2) * (b - 2)) // 2),
        "(m-2)(n-2)  [plane bidegree]": lambda a, b: (a - 2) * (b - 2),
        "RH: 1+((m-2)(n-2)-gcd(m-1,n-1)+1)//2":
            lambda a, b: 1 + ((a - 2) * (b - 2) - gcd(a - 1, b - 1) + 1) // 2,
        "RH-exact: ((m-2)(n-2)+1-gcd(m-1,n-1))//2 + ... ":
            lambda a, b: ((a - 1) * (b - 1) - (b - 1) + (b - 1)
                          - gcd(a - 1, b - 1) - (a - 1) + 2) // 2,
    }
    # The clean exact RH form derived from 2g-2 = -2(n-1) + (m-2)(n-1) + (n-1) - gcd(m-1,n-1)
    #   R_fin = deg_X Disc = (m-2)(n-1)   [verified pattern], R_inf=(n-1)-gcd
    #   2g-2 = -2(n-1) + (m-2)(n-1) + (n-1) - gcd(m-1,n-1)
    #        = (n-1)(m-3) - gcd(m-1,n-1) ... let's just give the closed form:
    def g_closed(a, b):
        # exact RH closed form; a=larger, b=smaller (genus symmetric so order
        # the inputs).  g = 1 + [ (M-1)(N-2) - (N-1) - gcd(M-1,N-1) ] / 2
        M, N = max(a, b), min(a, b)
        val = (M - 1) * (N - 2) - (N - 1) - gcd(M - 1, N - 1)
        return 1 + val // 2 if val % 2 == 0 else None
    cands["g_closed (exact RH closed form)"] = g_closed

    best = None
    for name, fn in cands.items():
        exact = 0
        lb_ok = 0
        for (a, b, g) in cop:
            try:
                v = fn(a, b)
            except Exception:
                v = None
            if v is not None and v == g:
                exact += 1
            if v is not None and v <= g:
                lb_ok += 1
        print(f"    {name:50s} exact={exact}/{len(cop)} lowerbound_ok={lb_ok}/{len(cop)}")
        if best is None or exact > best[1]:
            best = (name, exact)

    # verify R_fin pattern explicitly.  In the stored orientation (m,n) the
    # cover pi_X : C -> P^1_X has degree d=n-1 and
    #   R_fin = deg_X Disc_Y( g_n(Y) - f_m(X) ) = (m-1)(n-2).
    print("\n  ## R_fin pattern check  (pi_X: R_fin =?= (m-1)(n-2)):")
    rfin_ok = True
    bad = []
    for (m, n, gg, g_val, rf, ri, t2, sf) in scan['all_rows']:
        pred = (m - 1) * (n - 2)
        if rf != pred:
            rfin_ok = False
            bad.append((m, n, rf, pred))
    print(f"    R_fin == (m-1)(n-2) for ALL scanned pairs : {rfin_ok}"
          f"{'' if rfin_ok else '  mismatches: '+str(bad[:5])}")
    # exact RH closed form (m>n orientation), R_inf=(n-1)-gcd(m-1,n-1):
    #   2g-2 = -2(n-1) + (m-1)(n-2) + (n-1) - gcd(m-1,n-1)
    #        = (n-2)(m-1) - (n-1) - gcd(m-1,n-1)
    #   g = 1 + [ (m-1)(n-2) - (n-1) - gcd(m-1,n-1) ] / 2
    print("    => exact genus (m>n):")
    print("       g = 1 + [ (m-1)(n-2) - (n-1) - gcd(m-1,n-1) ] / 2")
    # verify this closed form reproduces the table
    cf_ok = True
    for (m, n, gg, g_val, rf, ri, t2, sf) in scan['all_rows']:
        num = (m - 1) * (n - 2) - (n - 1) - gcd(m - 1, n - 1)
        gcf = 1 + num // 2 if num % 2 == 0 else None
        if gcf != g_val:
            cf_ok = False
    print(f"       closed form reproduces ENTIRE genus table : {cf_ok}")
    return dict(seq_m3=seq_m3, oeis=oeis_results, best_fit=best,
                rfin_pattern_ok=rfin_ok, closed_form_ok=cf_ok)


# ---------------------------------------------------------------------------
# TASK 5 : failure search over F_q[t]  (tame + wild probe)
# ---------------------------------------------------------------------------

def _pmul(a, b, p):
    if not a or not b:
        return []
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    res[i + j] = (res[i + j] + ai * bj) % p
    return res


def _padd(a, b, p):
    L = max(len(a), len(b))
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
            for i in range(L)]


def _psub(a, b, p):
    L = max(len(a), len(b))
    return [((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p
            for i in range(L)]


def _deg(a):
    d = -1
    for i, c in enumerate(a):
        if c != 0:
            d = i
    return d


def fm_of(coeffs, m, p):
    """f_m(X(t)) = 1 + X + ... + X^{m-1} in F_p[t]."""
    out = [1]
    powX = [1]
    for _ in range(1, m):
        powX = _pmul(powX, coeffs, p)
        out = _padd(out, powX, p)
    return out


def eq_holds(cx, cy, m, n, p):
    """Honest equation f_m(X)=g_n(Y) as a polynomial identity in F_p[t]."""
    diff = _psub(fm_of(cx, m, p), fm_of(cy, n, p), p)
    return all(c == 0 for c in diff)


def search_fq(p, D, m, n, max_enum=200000):
    """
    Brute-force X(t),Y(t) in F_p[t], deg in [1..D], for non-constant solutions
    of f_m(X)=g_n(Y).  Both must be non-constant (f_m of a non-constant poly
    has positive degree, can't equal g_n of a constant) and satisfy the
    degree balance (m-1)deg X = (n-1)deg Y.  Returns list of solutions.

    FAST: for each balanced (dx,dy), hash all f_m(X) values (X of degree dx)
    into a dict keyed by the tuple of coefficients, then probe with each
    g_n(Y) value (Y of degree dy).  This is an O(N) hash-join instead of the
    O(N^2) double loop.
    """
    vecs = list(product(range(p), repeat=D + 1))
    by_deg = {}
    for v in vecs:
        by_deg.setdefault(_deg(v), []).append(v)
    sols = []
    for dy in range(1, D + 1):
        num = (n - 1) * dy
        if num % (m - 1) != 0:
            continue
        dx = num // (m - 1)
        if dx < 1 or dx > D:
            continue
        xs = by_deg.get(dx, [])
        ys = by_deg.get(dy, [])
        if not xs or not ys:
            continue
        # hash f_m(X) -> X candidate(s)
        table = {}
        for cx in xs:
            key = tuple(c % p for c in fm_of(cx, m, p))
            # trim trailing zeros for canonical key
            while key and key[-1] == 0:
                key = key[:-1]
            table.setdefault(key, []).append(cx)
        for cy in ys:
            key = tuple(c % p for c in fm_of(cy, n, p))
            while key and key[-1] == 0:
                key = key[:-1]
            if key in table:
                for cx in table[key]:
                    sols.append((cx, cy, dx, dy))
    return sols


def task5_failure_search():
    print("\n" + "=" * 74)
    print("TASK 5 : failure search over F_q[t]  (tame + wild probe)")
    print("=" * 74)

    # --- tame search: q in {5,7,11}, deg<=4, m,n<=7, p \nmid mn -------------
    print("\n  --- TAME search (q in {5,7,11}, deg<=4, m,n<=7, p|mn skipped) ---")
    print(f"  {'q':>3} {'(m,n)':>8} {'tame?':>6} {'deg<=':>6}  result")
    qs = [5, 7, 11]
    mn = sorted(set((a, b) for a in range(3, 8) for b in range(3, 8) if a != b))
    Dmax = 4
    tame_total = 0
    for p in qs:
        for (m, n) in mn:
            tame = (m * n) % p != 0
            if not tame:
                print(f"  {p:>3} {str((m,n)):>8} {'WILD':>6} {'-':>6}  skip (p|mn)")
                continue
            # hash-join search is O(N): deg<=4 affordable for q up to 11
            D = Dmax
            while (p ** (D + 1)) > 300000 and D > 2:
                D -= 1
            sols = search_fq(p, D, m, n)
            tame_total += len(sols)
            res = "none" if not sols else f"FOUND {len(sols)}: {sols[:2]}"
            print(f"  {p:>3} {str((m,n)):>8} {'tame':>6} {D:>6}  {res}")
    print(f"\n  TAME non-constant solutions found (total): {tame_total}")

    # --- wild probe: p | mn (leaving tame regime) --------------------------
    print("\n  --- WILD probe (p | mn): do non-constant solutions appear? ---")
    print(f"  {'q':>3} {'(m,n)':>8} {'p|mn':>6} {'deg<=':>6}  result")
    wild_total = 0
    wild_hits = []
    # pick pairs where a small prime divides mn
    wild_cases = [
        (5, 3, 5),   # p=5 | n=5
        (5, 5, 3),   # p=5 | m=5
        (5, 5, 6),   # p=5 | m=5
        (7, 7, 3),   # p=7 | m=7
        (7, 3, 7),   # p=7 | n=7
        (3, 3, 4),   # p=3 | m=3
        (3, 4, 3),   # p=3 | n=3
        (3, 6, 5),   # p=3 | m=6
        (2, 4, 3),   # p=2 | m=4
        (2, 3, 4),   # p=2 | n=4
    ]
    for (p, m, n) in wild_cases:
        D = 4
        while (p ** (D + 1)) > 300000 and D > 2:
            D -= 1
        sols = search_fq(p, D, m, n)
        wild_total += len(sols)
        if sols:
            wild_hits.append((p, m, n, sols[:3]))
        res = "none" if not sols else f"FOUND {len(sols)}: {sols[:3]}"
        print(f"  {p:>3} {str((m,n)):>8} {'yes':>6} {D:>6}  {res}")
    print(f"\n  WILD non-constant solutions found (total): {wild_total}")
    if wild_hits:
        print("  >>> wild characteristic DOES admit non-constant solutions:")
        for (p, m, n, s) in wild_hits:
            print(f"      q={p} (m,n)=({m},{n}): {s}")
        print("  (consistent with: tameness p\\nmid mn is necessary for P1.)")
    else:
        print("  (no wild non-constant solutions at this degree bound; the wild")
        print("   failure mode, if any, needs higher degree / Frobenius twists.)")

    return dict(tame_total=tame_total, wild_total=wild_total,
                wild_hits=wild_hits)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    global OUT_CSV
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    OUT_CSV = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "data", "genus_table.csv")
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

    print("=" * 74)
    print("L3 NUMERICAL GENUS  --  Proposition P1 (function-field Goormaghtigh)")
    print("=" * 74)

    scan = task1_scan()
    gcd_ctrl = task2_gcd_control()
    min_ctrl, min_allzero = task3_min_control()
    oeis = task4_oeis(scan)
    fail = task5_failure_search()

    # --------------------------- VERDICT -----------------------------------
    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    n_cop = len(scan['coprime_min3'])
    viol = scan['violations']
    genus_min = scan['genus_min']
    all_ge1 = (len(viol) == 0)
    nonconstant = fail['tame_total']
    print(f"  TASK1 coprime&min>=3 cases       : {n_cop}")
    print(f"  TASK1 g>=1 holds for ALL         : {all_ge1}  (violations={viol})")
    print(f"  TASK1 MIN genus                  : {genus_min}")
    print(f"  TASK1 two-projection symmetry    : {'OK' if not scan['sym_fail'] else 'FAIL'}")
    print(f"  TASK2 gcd>1 control              : "
          f"{[(m,n,g) for (m,n,_,g,_,_) in gcd_ctrl]} "
          f"reducible_any={any(not irr for *_ ,irr,_ in [(0,0,0,0,r[4],r[5]) for r in gcd_ctrl])}")
    print(f"  TASK3 min<3 control all g==0     : {min_allzero}")
    print(f"  TASK4 R_fin=(m-1)(n-2) pattern   : {oeis['rfin_pattern_ok']}")
    print(f"  TASK4 closed-form genus reproduces table : {oeis['closed_form_ok']}")
    print(f"  TASK4 best genus fit             : {oeis['best_fit']}")
    print(f"  TASK5 tame non-constant solns    : {nonconstant}  (expect 0)")
    print(f"  TASK5 wild non-constant solns    : {fail['wild_total']}  "
          f"hits={len(fail['wild_hits'])}")

    status = 'ok' if (all_ge1 and nonconstant == 0 and not scan['sym_fail']) \
        else 'FAIL'
    oeis_hit = any(hits for (_, (_, _, hits)) in oeis['oeis'].items())
    print(f"\n  STATUS = {status}  genus_min_gcd1={genus_min}  "
          f"nonconstant={nonconstant}  oeis_hit={oeis_hit}")
    print(f"\nDECISION_LINE status={status} genus_min_gcd1={genus_min} "
          f"nonconstant={nonconstant} oeis_hit={oeis_hit} "
          f"wild_nonconstant={fail['wild_total']}")

    # machine-readable summary
    summary = dict(status=status, genus_min_gcd1=genus_min,
                   nonconstant=nonconstant, n_coprime_min3=n_cop,
                   violations=viol, two_proj_symmetry=(not scan['sym_fail']),
                   wild_nonconstant=fail['wild_total'],
                   rfin_pattern_ok=oeis['rfin_pattern_ok'],
                   best_fit=oeis['best_fit'], oeis_hit=oeis_hit)
    print("\nSUMMARY_JSON " + json.dumps(summary, default=str))
    return summary


if __name__ == "__main__":
    main()
