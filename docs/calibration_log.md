# Calibration log: every case, the values ParFlow read, and the outcome

Built by `scripts/case_table.py` from each run's `.pfidb` and driver files; outcomes from `results/<run>_daily.csv`.
Observed WY2008: Pemberton 3.27 annual / 1.66 summer m3/s, Toms 4.50 / 2.77.  Geology K as shipped: g5 0.06 m/h = 1.4 m/d, vertical ratio 0.1.

| run | start state | flow barrier | geology K x (g5 m/d) | vertical/horizontal (g5) | vertical K g5 m/d | Manning's n channel / hillslope min | overland | interception alpha | dewmx | cropland LAI / vcmx25 | urban top classes | CLM start | Pemberton ann / summer m3/s | Pemberton KGE | Toms ann / summer m3/s | Toms KGE | ET bias US-Slt | ET bias US-Ced |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pb10_kin | Delaware end state (dry) | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 1.0 | 0.2 | 6.00 / 50. | no | cold | 0.00 / nan | nan | 0.00 / nan | nan | +0.63 | -0.92 |
| pb_spin_c1 | CONUS2 steady state | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 1.0 | 0.2 | 6.00 / 50. | no | cold | 3.33 / 0.01 | -1.48 | 3.52 / 0.01 | -2.00 | +0.59 | +0.96 |
| pb2_y1 | kinematic steady state, package | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 3.65 / 0.02 | -1.65 | 3.77 / 0.01 | -2.11 | +0.56 | +0.66 |
| pb2_rec_c2 | pb2 recursion | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.81 / 0.02 | -1.22 | 2.94 / 0.01 | -1.91 | +0.52 | +0.58 |
| pb2_rec_c4 | pb2 recursion, converged | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.80 / 0.02 | -1.22 | 2.93 / 0.01 | -1.91 | +0.52 | +0.58 |
| pb2_k01 | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 0.1 (0) | 0.1 | 0.01 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.81 / 0.02 | -1.25 | 2.93 / 0.01 | -1.92 | +0.52 | +0.59 |
| pb2_k10 | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 10 (14) | 0.1 | 1.44 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.89 / 0.04 | -1.04 | 3.12 / 0.05 | -1.74 | +0.56 | +0.46 |
| pb2_nofb | pb2 recursion cycle 1 end (= control start) | off | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.86 / 0.02 | -1.04 | 3.10 / 0.02 | -1.72 | +0.55 | +0.54 |
| pb2_k10_nofb | pb2 recursion cycle 1 end (= control start) | off | 10 (14) | 0.1 | 1.44 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 3.36 / 0.14 | -0.80 | 3.82 / 0.27 | -1.42 | +0.56 | +0.21 |
| pb2_an001 | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 1 (1) | 0.01 | 0.01 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.81 / 0.02 | -1.24 | 2.93 / 0.01 | -1.92 | +0.52 | +0.59 |
| pb2_an1 | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 1 (1) | 1.0 | 1.44 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.83 / 0.02 | -1.09 | 3.02 / 0.02 | -1.83 | +0.54 | +0.55 |
| pb2_mann | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.15 / 0.100 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.77 / 0.03 | -0.71 | 2.89 / 0.03 | -1.18 | +0.53 | +0.58 |
| pb2_diff | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | diffusive correction, Lagged | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.81 / 0.02 | -1.16 | 2.93 / 0.01 | -2.05 | +0.53 | +0.58 |
| pb2_dwe | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 1 (1) | 0.1 | 0.14 | 0.09 / 0.020 | OverlandDiffusive (DWE) | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.82 / 0.02 | -1.20 | 2.95 / 0.01 | -2.10 | +0.53 | +0.58 |
| pb2_kh100 | pb2 recursion cycle 1 end (= control start) | 0.001 at 7 m | 100 (144) | 0.001 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 2.97 / 0.04 | -1.06 | 3.27 / 0.11 | -1.55 | +0.52 | +0.44 |
| pb2_fb42_k30 | pb2 recursion cycle 1 end (= control start) | 0.001 at 42 m | 30 (43) | 0.1 | 4.32 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 3.68 / 0.49 | -0.85 | 4.11 / 0.61 | -1.38 | +0.54 | -0.07 |
| pb2_fb42_kh100 | pb2 recursion cycle 1 end (= control start) | 0.001 at 42 m | 100 (144) | 0.001 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 4.82 / 1.54 | -0.83 | 4.82 / 1.15 | -1.14 | +0.61 | -0.66 |
| pb3_y1 | pb2_fb42_kh100 end state | 0.001 at 42 m | 100 (144) | 0.001 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 4.20 / 1.38 | -0.30 | 4.11 / 1.04 | -0.69 | +0.62 | -0.84 |
| pb3_rec_c1 | pb3 kinematic steady state | 0.001 at 42 m | 100 (144) | 0.001 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 4.01 / 1.31 | -0.21 | 4.28 / 1.03 | -0.75 | +0.63 | -0.90 |
| pb3_rec_c4 | pb3 recursion, converged | 0.001 at 42 m | 100 (144) | 0.001 | 0.14 | 0.09 / 0.020 | kinematic | 0.25 | 0.1 | 4.5 / 24 | yes | restart | 3.83 / 1.29 | -0.09 | 3.93 / 1.01 | -0.59 | +0.60 | -0.90 |
