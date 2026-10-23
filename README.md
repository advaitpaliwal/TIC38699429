# TIC 38699429.01: a 2.7 Earth-radius planet candidate in the habitable zone of an M dwarf

Discovered by Advait Paliwal on 2026-10-20 in public TESS full-frame-image photometry (search result file `data/search/38699429_search_result.json`, written 2026-10-20 11:49 UTC).

## Summary

TIC 38699429 (Gaia DR3 4676612253064169344) is a late M dwarf (Teff 2914 K, 0.166 R☉, 0.137 M☉) at 57.5 pc in the TESS southern continuous viewing zone. TESS-SPOC light curves from 31 sectors (2018 to 2023) show 60 transits with period 11.718167 ± 0.000013 d, depth 22,130 ± 750 ppm and duration 1.47 h. A transit fit with a stellar-density prior gives a radius of 2.72 (+0.13 / −0.12) Earth radii. The star's insolation at that orbit is 0.66 ± 0.15 times Earth's. Odd and even depths agree (0.41σ), there is no secondary eclipse, and the depth is steady across five years. Difference images rule out all 25 Gaia stars within 150″ brighter than G = 20.5 as the source.

**This is a planet candidate, not a validated or confirmed planet.** No high-resolution imaging exists for this star.

False-positive probability (TRICERATOPS v1.1.0, 10 runs of 10⁶ draws, no contrast curve):

- With no neighbours removed: **FPP 0.064 ± 0.007, NFPP 0.039 ± 0.005.** Most of the NFPP is an eclipsing binary on TIC 38696697 (41″), which the pixel fits rule out at Δχ² = 193.
- With the 23 neighbours that the TESS pixels rule out in every sector removed (Giacalone et al. 2021, Section 2.2): **FPP 0.026 ± 0.003, NFPP 0** in all ten runs. This exclusion rests on the TESS pixels, not on ground-based photometry. What remains is mainly a planet around an unresolved bound companion (0.023).

Both values are above the 0.015 threshold for statistical validation.

"Habitable zone" here means the 0.25 to 1.7 times Earth's insolation window used by the search. No habitable-zone boundary model is applied.

As of 2026-10-16 the star is not a TOI or a CTOI (ExoFOP), and it is not in the NASA Exoplanet Archive. TESS-SPOC data validation for the multi-sector run over Sectors 56 to 69 reported it as a threshold-crossing event (P = 11.7183 d, MES 11.9, not flagged as an eclipsing binary, radius 3.9 ± 4.7 R⊕), but it never became a TOI (`data/vet/02b_tess_spoc_tce_params.txt`, `data/vet/02_known.txt`).

| | |
|---|---|
| Star | TIC 38699429, RA 66.18828°, Dec −62.00852° (J2000), Tmag 15.27, Gaia G 16.80 |
| Period | 11.718167 ± 0.000013 d |
| Mid-transit (T0) | BJD_TDB 2459196.67148 ± 0.00072 |
| Depth | 22,130 ± 750 ppm |
| Duration | 1.47 ± 0.05 h |
| Rp/R★ | 0.150 ± 0.007 |
| Impact parameter | 0.27 (+0.16 / −0.18) |
| Radius | 2.72 (+0.13 / −0.12) R⊕ |
| Insolation | 0.66 ± 0.15 S⊕ |
| Equilibrium temperature | 230 ± 13 K |
| False-positive probability | 0.064 (no neighbours removed), 0.026 (pixel-excluded neighbours removed) |

Period and T0 are from the linear-ephemeris fit in `data/fit/ephemeris_fit.json`. Transit Least Squares on the same light curve gives P = 11.71810 ± 0.00046 d (`data/vet/01c_tls_refit_narrow.json`), the value quoted in the note. Depth is the box fit in `data/vet/01a_fold_secondary_oddeven.json`; shape and radius are from `data/vet/01b_shape_fit.json`; insolation and equilibrium temperature are from the SPOC data validation report.

## Contents

- `note/TIC38699429_note.md`, `note/TIC38699429_note.pdf`: the research note (abstract, star and search, transit fit, vetting, false-positive analysis).
- `note/figure1_TIC38699429.png`: all 31 sectors folded on the ephemeris with the transit model, an unfolded sector, single transits, odd and even transits, and phase 0.5.
- `ephemeris.csv`: predicted mid-transit times from 2026-10-23 to 2027-12-31.
- `code/search/`: the red-dwarf search (GPU box search with cuvarbase and the earlier Transit Least Squares search, light-curve preparation, and the LEO-vetter gate).
- `code/vet/`: the vetting scripts (fold, secondary and odd/even, transit-shape fit, depths by data product, difference images, TRICERATOPS with and without neighbour removal).
- `code/fit/ephemeris.py`: the ephemeris fit and transit predictions.
- `code/figure/fig1.py`: Figure 1.
- `data/search/`: the search result, the pool row and the light curve the search used (`38699429_lc.npz`, 31 TESS-SPOC sectors).
- `data/vet/`: small outputs of the vetting scripts, including all TRICERATOPS runs (`trc/` no removal, `trc_rm/` with removal), difference-image fits and plots.
- `data/fit/ephemeris_fit.json`: the ephemeris posterior summary.
- `data/exofop_ctoi_row.txt`: the prepared ExoFOP community-candidate parameter row (not yet uploaded).

## Next predicted transits

From `ephemeris.csv` (1σ timing uncertainty about 3.5 minutes now, growing to 4.1 minutes by end of 2027):

| Epoch | Mid-transit BJD_TDB | Mid-transit UTC | ±1σ |
|---|---|---|---|
| 257 | 2461341.09604 | 2026-10-27 14:15 | 3.5 min |
| 258 | 2461352.81421 | 2026-11-08 07:29 | 3.5 min |
| 259 | 2461364.53238 | 2026-11-20 00:44 | 3.5 min |
| 260 | 2461376.25054 | 2026-12-01 17:58 | 3.5 min |
| 261 | 2461387.96871 | 2026-12-13 11:13 | 3.5 min |
| 262 | 2461399.68688 | 2026-12-25 04:27 | 3.6 min |

The transit shape is held at the best density-prior fit and only T0 and P are fitted (emcee, 1,702 points within 0.3 d of transit). Times are propagated from every posterior sample, so the correlation between T0 and P is kept. UTC times are geocentric (BJD_TDB minus the light-travel delay to the star, then TDB to UTC); ingress and egress are mid-transit ± 0.735 h. Epochs count from the first Sector 1 transit (BTJD 1329.528); T0 is epoch 74. The transit is about 2% deep on a T = 15.3 star.

## Reproducing

Python 3.12 with numpy, scipy, astropy, astroquery, lightkurve, wotan, transitleastsquares, cuvarbase 0.3.0 (GPU search), LEO-vetter v1.2.0, batman-package, emcee, TESS_PRF and triceratops 1.1.0.

1. Light curves. The search builds each star's light curve from MAST (`code/search/run1.py` and `run1_body.py`: SPOC 2-minute PDCSAP, else TESS-SPOC PDCSAP, else QLP; wotan biweight 0.75 d; 30-minute bins). The light curve used for this star is included as `data/search/38699429_lc.npz` (arrays `t` in BTJD, `f`, `s` sector). FITS light curves and TESScut cutouts are not included; they are public at MAST.
2. Search. `code/search/run_star_gpu.py` (cuvarbase box search, up to four passes) driven by `drive_gpu.py`; `run_star.py` and `drive.py` are the Transit Least Squares version used for the first 1,133 stars. Light-curve preparation helpers are `v4_test.py`, `fix_test.py`, `trim_test.py` and `v3_test.py`. Peaks go through `gate.py` and `leo_glue.py` (LEO-vetter flux-level tests); `control.py` runs the known-planet and known-artifact controls.
3. Vetting. Run from a folder containing the light curve: `python -I code/vet/fold.py 38699429.npz`, `code/vet/shape.py`, `code/vet/tls_refit_narrow.py`, `code/vet/products.py` and `products_stats.py` (depths by data product), `code/vet/diff.py` (difference images and PRF fits at each Gaia position), `code/vet/trc_setup.py`, `trc_run.py` and `trc_all.sh` (TRICERATOPS, no removal), `trc_rm_match.py`, `trc_rm_run.py` and `trc_rm_all.sh` (with removal), and `trc_summary.py`. `gaia_patch.py` routes TRICERATOPS's Gaia field query through VizieR. Outputs are in `data/vet/`.
4. Ephemeris. `python -I code/fit/ephemeris.py data/search/38699429_lc.npz data/vet/01b_shape_fit.json data/fit/ephemeris_fit.json ephemeris.csv`.
5. Figure. `python -I code/figure/fig1.py data/search/38699429_lc.npz data/vet/01b_shape_fit.json figure1_TIC38699429.png`.

Some scripts use paths relative to the original working folders (for example the search lanes `hz-south/` and `edgefix/`); adjust them to your layout.

## Credits

This work uses data collected by the TESS mission, funded by the NASA Science Mission Directorate, obtained from the Mikulski Archive for Space Telescopes (MAST), including TESS-SPOC (Caldwell et al. 2020) light curves, TGLC (Han & Brandt 2023) photometry and TESScut cutouts. It uses data from the European Space Agency mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium, and the TESS Input Catalog (Stassun et al. 2019; Paegert et al. 2021). Full references are in the note.

## Licence and citation

Code is under the MIT licence (`LICENSE`). Text, figures and data are under CC BY 4.0 (`LICENSE-CC-BY-4.0`). Please cite with `CITATION.cff`.

## Contact

Advait Paliwal, advaitpaliwal@gmail.com
