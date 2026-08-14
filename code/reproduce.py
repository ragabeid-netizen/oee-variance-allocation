# -*- coding: utf-8 -*-
"""
Reproduces every quantity reported in

    "Why Overall Equipment Effectiveness Tracks Its Performance Term:
     Variance Allocation and Operating-Range Asymmetry in a Six-Year Plant Record"

Run from the repository root:   python code/reproduce.py
Requires: pandas, numpy.  Runtime: about two minutes.
"""
import pandas as pd, numpy as np, os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHIFTS = os.path.join(HERE, 'data', 'shift_records_2020_2025.csv')
MAINT = os.path.join(HERE, 'data', 'maintenance_daily_register.csv')
SEED = 20260815

d = pd.read_csv(SHIFTS).dropna(subset=['oee']).sort_values('seq').reset_index(drop=True)
m = pd.read_csv(MAINT)
m.columns = [c.strip() for c in m.columns]
MON = ["January", "February", "March", "April", "May", "June", "July", "August",
       "September", "October", "November", "December"]
m['mn'] = m['Month'].map({x: i + 1 for i, x in enumerate(MON)})
m = m[m.Year.between(2020, 2025)].copy()
m['A_day'] = (m['Operating (min)'] - m['Downtime (min)']) / m['Operating (min)']
mon = m.groupby(['Year', 'mn'], as_index=False).agg(D=('Downtime (min)', 'sum'),
                                                    O=('Operating (min)', 'sum'))
mon['A_month'] = (mon.O - mon.D) / mon.O
d = d.merge(mon[['Year', 'mn', 'A_month']], left_on=['year', 'month'],
            right_on=['Year', 'mn'], how='left')


def shares(A, P, Q):
    """s(X) = 100 * Cov(ln X, ln OEE) / Var(ln OEE); the three sum to 100 exactly."""
    A, P, Q = map(np.asarray, (A, P, Q))
    lo = np.log(A * P * Q); v = lo.var(ddof=1)
    return np.array([100 * np.cov(np.log(X), lo, ddof=1)[0, 1] / v for X in (A, P, Q)])


def hdr(t):
    print('\n' + '=' * 78); print(t); print('=' * 78)


A0, P0, Q0 = d.availability.values, d.performance.values, d.quality.values
n = len(d)
base = shares(A0, P0, Q0)

hdr('TABLE 2  -  variance allocation')
print(f'  n = {n}')
for k, nm in enumerate(['Performance', 'Availability', 'Quality']):
    i = [1, 0, 2][k]
    print(f'  {nm:<14} {base[i]:>8.2f}%   distinct values '
          f'{[len(np.unique(x)) for x in (P0, A0, Q0)][k]:>5}')
print(f'  r(OEE, Performance) = {np.corrcoef(A0*P0*Q0, P0)[0,1]:.4f}')

hdr('TABLE 3  -  five resampling schemes, 5,000 replicates each')
rng = np.random.default_rng(SEED)
B = 5000
schemes = {'Ordinary i.i.d. shifts': ('iid', None)}
d['cl_my'] = d.machine.astype(str) + '|' + d.year.astype(str)
d['cl_mm'] = d.cl_my + '|' + d.month.astype(str)
for lab, col in [('Cluster, machine-month', 'cl_mm'), ('Cluster, machine-year', 'cl_my'),
                 ('Cluster, calendar-year', 'year')]:
    schemes[lab] = ('cluster', [np.asarray(v) for v in d.groupby(col).indices.values()])
schemes['Moving block, L = 30'] = ('block', 30)

print(f"  {'scheme':<26}{'s(A) 95% CI':>22}{'s(P) 95% CI':>24}")
for lab, (kind, arg) in schemes.items():
    out = np.empty((B, 3))
    for b in range(B):
        if kind == 'iid':
            i = rng.integers(0, n, n)
        elif kind == 'cluster':
            G = len(arg); i = np.concatenate([arg[j] for j in rng.integers(0, G, G)])
        else:
            L = arg
            i = np.concatenate([np.arange(s, s + L) for s in
                                rng.integers(0, n - L, int(np.ceil(n / L)))])[:n]
        out[b] = shares(A0[i], P0[i], Q0[i])
    lo, hi = np.nanpercentile(out, [2.5, 97.5], axis=0)
    print(f'  {lab:<26}[{lo[0]:>8.2f},{hi[0]:>8.2f}]   [{lo[1]:>9.2f},{hi[1]:>9.2f}]')

hdr('TABLE 5  -  resolution ladder')
ann_s, mon_s = shares(A0, P0, Q0), shares(d.A_month, P0, Q0)
day_map = {k: g['A_day'].values for k, g in m.groupby(['Year', 'mn'])}
key = list(zip(d.year.astype(int), d.month.astype(int)))
R = 500
dr = np.array([shares(np.array([rng.choice(day_map[k]) for k in key]), P0, Q0)
               for _ in range(R)])
dm, dlo, dhi = dr.mean(0), *np.percentile(dr, [2.5, 97.5], axis=0)
print(f"  {'resolution':<24}{'distinct':>10}{'SD (pp)':>10}{'s(A)':>9}{'s(P)':>10}")
print(f"  {'annual (as recorded)':<24}{len(np.unique(A0)):>10}{A0.std(ddof=1)*100:>10.2f}"
      f"{ann_s[0]:>9.2f}{ann_s[1]:>10.2f}")
print(f"  {'monthly':<24}{d.A_month.nunique():>10}{d.A_month.std(ddof=1)*100:>10.2f}"
      f"{mon_s[0]:>9.2f}{mon_s[1]:>10.2f}")
print(f"  {'daily (mean of 500)':<24}{m.A_day.nunique():>10}{m.A_day.std(ddof=1)*100:>10.2f}"
      f"{dm[0]:>9.2f}{dm[1]:>10.2f}")
print(f"  daily 95% interval: s(A) [{dlo[0]:.2f}, {dhi[0]:.2f}]  "
      f"s(P) [{dlo[1]:.2f}, {dhi[1]:.2f}]")

hdr('SECTION 4.2  -  CV values, Eq. (4), and the covariance sweep')
cv = lambda x: np.std(x, ddof=1) / np.mean(x)
cA, cP, cQ = cv(A0), cv(P0), cv(Q0)
print(f'  CV(A) = {cA:.4f}   CV(P) = {cP:.4f}   CV(Q) = {cQ:.4f}')
print(f'  CV(P)/CV(A) = {cP/cA:.1f}     CV(P)/CV(Q) = {cP/cQ:.0f}')
print(f'  Eq.(4) approximation for s(P) = '
      f'{100*cP**2/(cA**2+cP**2+cQ**2):.2f}%   exact = {base[1]:.2f}%')
lA, lP, lQ = np.log(A0), np.log(P0), np.log(Q0)
C = np.corrcoef(np.vstack([lA, lP, lQ]))
sA, sP, sQ = lA.std(ddof=1), lP.std(ddof=1), lQ.std(ddof=1)
print(f'  Cov(lnA,lnP) = {np.cov(lA,lP,ddof=1)[0,1]:.3e}   Var(lnA) = {lA.var(ddof=1):.3e}')
print(f'\n  {"r(lnA,lnP)":>12}{"s(A)":>9}{"s(P)":>9}')
for r in [-0.9, -0.5, C[0, 1], 0.0, 0.5, 0.9]:
    cAP, cAQ, cPQ = r*sA*sP, C[0, 2]*sA*sQ, C[1, 2]*sP*sQ
    V = sA**2+sP**2+sQ**2+2*(cAP+cAQ+cPQ)
    tag = '   <- observed' if abs(r-C[0, 1]) < 1e-9 else ''
    print(f'  {r:>12.3f}{100*(sA**2+cAP+cAQ)/V:>9.2f}{100*(sP**2+cAP+cPQ)/V:>9.2f}{tag}')

hdr('SECTION 5.3  -  dispersion thresholds from Eq. (5)')
for tgt in [0.01, 0.05, 0.10]:
    need = np.sqrt(tgt/(1-tgt)*(cP**2+cQ**2))
    print(f'  s(A) = {tgt:>5.0%}  requires CV(A) = {need:.5f}  '
          f'i.e. SD = {need*np.mean(A0)*100:5.2f} pp  = {need/cA:5.1f}x observed')
print(f'  s(Q) = 10%  requires CV(Q) to rise {np.sqrt(0.10/0.90*cP**2)/cQ:.0f}x')

hdr('TABLE 6  -  observed operating range')
print(f"  {'rate':<14}{'min':>9}{'median':>9}{'max':>9}{'range pp':>10}{'SD pp':>8}")
for nm, x in [('Availability', A0), ('Quality', Q0), ('Performance', P0)]:
    print(f'  {nm:<14}{x.min():>9.4f}{np.median(x):>9.4f}{x.max():>9.4f}'
          f'{(x.max()-x.min())*100:>10.2f}{x.std(ddof=1)*100:>8.2f}')

hdr('TABLE 4  -  strata')
for col, lab in [('machine', 'Machine'), ('factory', 'Factory'), ('product', 'Product'),
                 ('speed_setting', 'Speed preset'), ('shift_regime', 'Shift regime'),
                 ('year', 'Calendar year')]:
    v = [shares(g.availability, g.performance, g.quality)[1] for _, g in d.groupby(col)
         if len(g) > 5]
    r = [np.corrcoef(g.availability*g.performance*g.quality, g.performance)[0, 1]
         for _, g in d.groupby(col) if len(g) > 5]
    print(f'  {lab:<15}{len(v):>3} strata   s(P) {min(v):7.2f} - {max(v):7.2f}%   '
          f'r {min(r):.4f} - {max(r):.4f}')

hdr('SERIAL DEPENDENCE')
lo_ = np.log(A0*P0*Q0)
print('  lag-1..8 autocorrelation of ln OEE: ' +
      '  '.join(f'{np.corrcoef(lo_[:-k], lo_[k:])[0,1]:+.3f}' for k in range(1, 9)))
r1 = np.corrcoef(lo_[:-1], lo_[1:])[0, 1]
print(f'  pooled rho = {r1:+.4f}  ->  effective n = {int(n*(1-r1)/(1+r1))} of {n}')
print('\nDone.')
