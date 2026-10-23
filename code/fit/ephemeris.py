# Linear ephemeris (T0, P) of TIC 38699429.01 and predicted mid-transit times, 2026-10-23 to 2027-12-31.
# The transit shape is fixed at the best density-prior fit (data/vet/01b_shape_fit.json, from code/vet/shape.py);
# only the reference mid-time and the period are free (emcee). Every posterior sample is propagated (T0 + n*P),
# so the T0-P correlation is kept. UTC times are geocentric: BJD_TDB minus the light-travel delay to the star, then TDB -> UTC.
# Usage: python -I code/fit/ephemeris.py <38699429.npz> data/vet/01b_shape_fit.json data/fit/ephemeris_fit.json ephemeris.csv
import sys, csv, json, numpy as np, batman, emcee
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
npz, shape, fitout, out = sys.argv[1:5]
d = np.load(npz); t, f = d['t'], d['f']
best = json.load(open(shape))['density_TIC_prior']['best']; sig = json.load(open(shape))['point_sigma_ppm']*1e-6
tc0, rp, b, lrho, q1, q2, f0 = best
P0 = 11.718154156303166; T00 = 1329.529459184984 + tc0          # box-refined ephemeris used by shape.py
nref = int(np.round((np.median(t) - T00)/P0)); Tref0 = T00 + nref*P0   # reference epoch near the middle of the data
G = 6.674e-8; a = (10**lrho*G*(P0*86400)**2/(3*np.pi))**(1/3)
dt0 = (((t-Tref0)/P0+0.5) % 1-0.5)*P0; m = np.abs(dt0) < 0.3; tt, ff = t[m], f[m]
pr = batman.TransitParams(); pr.t0 = 0; pr.per = P0; pr.rp = rp; pr.a = a; pr.inc = np.degrees(np.arccos(b/a))
pr.ecc = 0; pr.w = 90; pr.limb_dark = 'quadratic'; pr.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
def lp(th):
    T, P = th
    if abs(T-Tref0) > 0.05 or abs(P-P0) > 0.002: return -np.inf
    x = (((tt-T)/P+0.5) % 1-0.5)*P
    mod = batman.TransitModel(pr, x, supersample_factor=9, exp_time=30/1440).light_curve(pr)*f0
    return -0.5*np.sum(((ff-mod)/sig)**2)
nw = 16; p0 = np.array([Tref0, P0]) + np.array([1e-3, 1e-5])*np.random.default_rng(2).standard_normal((nw, 2))
S = emcee.EnsembleSampler(nw, 2, lp); S.run_mcmc(p0, 3000, progress=False)
ch = S.get_chain(discard=1000, flat=True)
Tref = ch[:, 0] + 2457000.0; P = ch[:, 1]
res = dict(T0_BJD_TDB=float(np.median(Tref)), T0_sigma_d=float(Tref.std()), P_d=float(np.median(P)), P_sigma_d=float(P.std()),
           corr_T0_P=float(np.corrcoef(Tref, P)[0, 1]), reference_epoch_from_TLS_T0=nref, n_points=int(m.sum()),
           point_sigma_ppm=sig*1e6, acceptance=float(np.mean(S.acceptance_fraction)), n_samples=int(len(P)))
json.dump(res, open(fitout, 'w'), indent=1); print(json.dumps(res, indent=1))
STAR = SkyCoord(66.1882793443072*u.deg, -62.0085242002372*u.deg)  # TIC v8.2, J2000
GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
T14_H = 1.47
t_start = Time('2026-10-23T00:00:00', scale='utc').tdb.jd; t_end = Time('2027-12-31T23:59:59', scale='utc').tdb.jd
t0m, pm = np.median(Tref), np.median(P)
n = np.arange(int(np.ceil((t_start-t0m)/pm))-1, int(np.floor((t_end-t0m)/pm))+2)
pred = Tref[None, :] + n[:, None]*P[None, :]
bjd = np.median(pred, 1); sg = pred.std(1)
def bjd_to_utc(bb):
    tq = Time(bb, format='jd', scale='tdb', location=GEO)
    for _ in range(3):
        tq = Time(bb, format='jd', scale='tdb', location=GEO) - tq.light_travel_time(STAR, 'barycentric')
        tq = Time(tq.jd, format='jd', scale='tdb', location=GEO)
    return tq.utc
utc = bjd_to_utc(bjd)
keep = (utc.jd >= Time('2026-10-23', scale='utc').jd) & (utc.jd <= Time('2028-01-01', scale='utc').jd)
half = T14_H/48
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['epoch', 'mid_BJD_TDB', 'sigma_d', 'sigma_min', 'mid_UTC', 'ingress_UTC', 'egress_UTC'])
    for i in np.where(keep)[0]:
        mm = utc[i]
        w.writerow([int(n[i]+nref), '%.5f' % bjd[i], '%.5f' % sg[i], '%.1f' % (sg[i]*1440),
                    mm.isot[:16], (mm-half*u.day).isot[:16], (mm+half*u.day).isot[:16]])
print('rows', keep.sum(), 'sigma_min range', sg[keep].min()*1440, sg[keep].max()*1440)
