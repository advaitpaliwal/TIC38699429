# Figure 1 for the TIC 38699429 note. Usage: python -I fig1.py <npz> <shape_fit.json> <out.png>
import sys, json, numpy as np, batman
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

d = np.load(sys.argv[1]); t, f, s = d['t'], d['f'], d['s']
fit = json.load(open(sys.argv[2]))['density_TIC_prior']
P, T0 = 11.718154156303166, 1329.529459184984          # box-refined ephemeris, as in shape.py
tc, rp, b, lrho, q1, q2, f0 = fit['best']
G = 6.674e-8; a = (10**lrho * G * (P*86400)**2 / (3*np.pi))**(1/3)

def model(x):
    pr = batman.TransitParams(); pr.t0 = tc; pr.per = P; pr.rp = rp; pr.a = a
    pr.inc = np.degrees(np.arccos(b/a)); pr.ecc = 0; pr.w = 90; pr.limb_dark = 'quadratic'
    pr.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
    return batman.TransitModel(pr, x, supersample_factor=9, exp_time=30/1440).light_curve(pr) * f0

dt = (((t-T0)/P+0.5) % 1 - 0.5) * P
ep = np.round((t-T0)/P).astype(int)
ppt = lambda y: (y-1)*1e3

def binned(x, y, w):
    e = np.arange(x.min(), x.max()+w, w); i = np.digitize(x, e)
    c = np.array([(e[k-1]+e[k])/2 for k in range(1, len(e)) if (i == k).sum() > 2])
    m = np.array([np.mean(y[i == k]) for k in range(1, len(e)) if (i == k).sum() > 2])
    r = np.array([np.std(y[i == k])/np.sqrt((i == k).sum()) for k in range(1, len(e)) if (i == k).sum() > 2])
    return c, m, r

fig = plt.figure(figsize=(7.5, 9.6))
gs = fig.add_gridspec(4, 4, height_ratios=[1.25, 0.9, 0.8, 0.9], hspace=0.55, wspace=0.35)
ink, grey, accent = 'k', '0.75', '#c4552b'
xm = np.linspace(-0.25, 0.25, 800)

# (a) folded light curve
ax = fig.add_subplot(gs[0, :]); w = np.abs(dt) < 0.25
ax.plot(dt[w]*24, ppt(f[w]), '.', color=grey, ms=2, zorder=1)
c, m, r = binned(dt[w], f[w], 15/1440)
ax.errorbar(c*24, ppt(m), r*1e3, fmt='o', color=ink, ms=3, lw=0.8, zorder=3)
ax.plot(xm*24, ppt(model(xm)), '-', color=accent, lw=1.4, zorder=4)
ax.set_xlim(-6, 6); ax.set_ylim(-45, 25); ax.set_xlabel('Hours from mid-transit'); ax.set_ylabel('Relative flux (ppt)')
ax.set_title('(a) 31 sectors folded on P = 11.71815 d', fontsize=9, loc='left')

# (b) unfolded Sector 63
ax = fig.add_subplot(gs[1, :]); m63 = s == 63
ax.plot(t[m63], ppt(f[m63]), '.', color=ink, ms=2)
for e in sorted(set(ep[m63])):
    tt = T0+e*P
    if t[m63].min() < tt < t[m63].max(): ax.axvline(tt, color=accent, lw=0.8, ymin=0.9, ymax=1.0)
ax.set_ylim(-60, 40); ax.set_xlabel('Time (BJD$_{\\rm TDB}$ $-$ 2457000)'); ax.set_ylabel('Rel. flux (ppt)')
ax.set_title('(b) Sector 63, unfolded, 30-min bins; ticks mark predicted transits', fontsize=9, loc='left')

# (c) single transits
for k, e in enumerate([144, 145, 146, 148]):
    ax = fig.add_subplot(gs[2, k]); x = t-(T0+e*P); w = np.abs(x) < 0.25
    ax.plot(x[w]*24, ppt(f[w]), 'o', color=ink, ms=2.5)
    ax.plot(xm*24, ppt(model(xm)), '-', color=accent, lw=1)
    ax.set_xlim(-6, 6); ax.set_ylim(-60, 30)
    ax.set_title(f'(c) epoch {e}, S{int(np.median(s[w]))}', fontsize=8, loc='left'); ax.set_xlabel('Hours', fontsize=8)
    ax.tick_params(labelsize=7)
    if k == 0: ax.set_ylabel('Rel. flux (ppt)', fontsize=8)

# (d) odd and even transits; (e) phase 0.5
ax = fig.add_subplot(gs[3, :2]); w = np.abs(dt) < 0.25
for sel, col, lab, off in [(ep % 2 == 1, ink, 'odd', 0), (ep % 2 == 0, accent, 'even', 0.1/24)]:
    c, m, r = binned(dt[w & sel], f[w & sel], 30/1440)
    ax.errorbar(c*24+off*24, ppt(m), r*1e3, fmt='o', color=col, ms=3, lw=0.8, label=lab)
ax.set_xlim(-6, 6); ax.set_ylim(-45, 25); ax.legend(fontsize=7, frameon=False, loc='lower right')
ax.set_xlabel('Hours from mid-transit'); ax.set_ylabel('Rel. flux (ppt)')
ax.set_title('(d) odd and even transits', fontsize=9, loc='left')

ax = fig.add_subplot(gs[3, 2:]); d5 = (((t-T0-P/2)/P+0.5) % 1 - 0.5)*P; w = np.abs(d5) < 0.25
ax.plot(d5[w]*24, ppt(f[w]), '.', color=grey, ms=2)
c, m, r = binned(d5[w], f[w], 30/1440)
ax.errorbar(c*24, ppt(m), r*1e3, fmt='o', color=ink, ms=3, lw=0.8)
ax.set_xlim(-6, 6); ax.set_ylim(-45, 25); ax.set_xlabel('Hours from phase 0.5')
ax.set_title('(e) phase 0.5', fontsize=9, loc='left')

fig.savefig(sys.argv[3], dpi=200, bbox_inches='tight')
print('model depth ppm', (1-model(np.array([tc]))[0])*1e6, 'a/Rs', a)
