# MiniBA results

Scene: `desk.txt`: 150 images, 40137 points, 790261 observations, one SIMPLE_RADIAL camera 1280x720. Seed 7, 12 threads.

Pose errors are measured against COLMAP/GLOMAP's solution after a 7-DoF similarity alignment (Umeyama on camera centres), as a percentage of the scene size (median camera distance to the centroid = 3.5526 model units).

## 1. Check: same projection model as COLMAP

| metric | MiniBA (C++) | pycolmap / COLMAP |
|---|---|---|
| mean over all observations | **0.624341 px** | 0.624341 px (pycolmap's own projection, `tools/check_colmap.py`) |
| mean of per-point means | **0.516551 px** | 0.516551 px (`Reconstruction.compute_mean_reprojection_error()`) |

COLMAP's headline number averages each 3D point's track first, then over points, so long tracks (which have larger errors) count less. Both match to 6 digits, so the projection model is the same. 0 observations behind a camera. All errors below are means over observations.

## 2. Polish: BA started from COLMAP's solution

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| COLMAP solution | 0.624 | 0.404 | 0.918 | 0.000 / 0.000 | 0.0000 / 0.0000 | 1323.24 | -0.01082 | - | - |
| polish (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0706 / 0.0855 | 1322.22 | -0.01124 | 5 | 3.94 |

A squared loss minimises RMS, so RMS can only go down; the solution barely moves.

## 3. Noisy start, clean matches

Every camera rotated by a random 0.5 deg (std), camera centres and all points moved by 1 % of the scene size (std per axis), focal length 5 % too long, distortion reset to 0.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 22.336 | 21.372 | 25.096 | 1.555 / 3.349 | 0.4926 / 1.5608 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0713 / 0.0864 | 1322.17 | -0.01126 | 8 | 5.56 |

## 4. Noisy start + 5 % bad matches

39768 of 790261 observations moved 20-80 px in a random direction (wrong matches). Errors below are over the clean observations only, measured against their true pixels. Robust losses start down-weighting at 2 px.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 22.335 | 21.371 | 25.095 | 1.555 / 3.349 | 0.4926 / 1.5608 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 2.022 | 1.468 | 2.977 | 0.812 / 3.858 | 0.4206 / 0.7756 | 1344.84 | -0.00212 | 35 | 19.82 |
| BA (Huber) | 0.657 | 0.434 | 1.027 | 0.067 / 0.698 | 0.0230 / 0.0533 | 1323.46 | -0.01072 | 49 | 27.10 |
| BA (Cauchy) | 0.626 | 0.401 | 0.982 | 0.019 / 0.206 | 0.0115 / 0.0189 | 1323.17 | -0.01079 | 70 | 37.68 |

Flagging every observation above 4 px after the solve:

| run | precision | recall |
|---|---|---|
| BA (squared) | 0.325 | 0.999 |
| BA (Huber) | 0.945 | 0.998 |
| BA (Cauchy) | 0.911 | 0.998 |

## 5. Linear solvers (noisy start, clean matches, Huber)

Points are eliminated with the Schur complement; what differs is how the reduced camera system (902 unknowns) is solved.

| linear solver | iters | time s | s / iter | final cost | mean px |
|---|---|---|---|---|---|
| DENSE_SCHUR | 60 | 33.10 | 0.552 | 304484.4 | 0.631 |
| SPARSE_SCHUR | 57 | 30.84 | 0.541 | 304484.4 | 0.631 |
| ITERATIVE_SCHUR + SCHUR_JACOBI | 38 | 17.04 | 0.448 | 304484.4 | 0.631 |
