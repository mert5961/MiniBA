# MiniBA results

Scene: `desk.txt`: 150 images, 40137 points, 790261 observations, one SIMPLE_RADIAL camera 1280x720. Seed 42, 12 threads.

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
| polish (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0706 / 0.0855 | 1322.22 | -0.01124 | 5 | 3.96 |

A squared loss minimises RMS, so RMS can only go down; the solution barely moves.

## 3. Noisy start, clean matches

Every camera rotated by a random 0.5 deg (std), camera centres and all points moved by 1 % of the scene size (std per axis), focal length 5 % too long, distortion reset to 0.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 21.264 | 20.241 | 24.016 | 1.634 / 3.790 | 0.4501 / 1.6931 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0713 / 0.0864 | 1322.17 | -0.01126 | 8 | 5.53 |

## 4. Noisy start + 5 % bad matches

39581 of 790261 observations moved 20-80 px in a random direction (wrong matches). Errors below are over the clean observations only, measured against their true pixels. Robust losses start down-weighting at 2 px.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 21.267 | 20.245 | 24.019 | 1.634 / 3.790 | 0.4501 / 1.6931 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 2.016 | 1.472 | 2.956 | 0.730 / 4.583 | 0.0939 / 0.3964 | 1325.53 | -0.00934 | 40 | 22.45 |
| BA (Huber) | 0.657 | 0.434 | 1.033 | 0.075 / 0.591 | 0.0570 / 0.0872 | 1322.61 | -0.01102 | 99 | 52.74 |
| BA (Cauchy) | 0.626 | 0.401 | 0.978 | 0.026 / 0.260 | 0.0291 / 0.0364 | 1323.13 | -0.01078 | 119 | 62.15 |

Flagging every observation above 4 px after the solve:

| run | precision | recall |
|---|---|---|
| BA (squared) | 0.327 | 0.999 |
| BA (Huber) | 0.945 | 0.998 |
| BA (Cauchy) | 0.911 | 0.998 |

## 5. Linear solvers (noisy start, clean matches, Huber)

Points are eliminated with the Schur complement; what differs is how the reduced camera system (902 unknowns) is solved.

| linear solver | iters | time s | s / iter | final cost | mean px |
|---|---|---|---|---|---|
| DENSE_SCHUR | 69 | 37.84 | 0.548 | 304484.4 | 0.631 |
| SPARSE_SCHUR | 70 | 37.11 | 0.530 | 304484.4 | 0.631 |
| ITERATIVE_SCHUR + SCHUR_JACOBI | 45 | 20.40 | 0.453 | 304484.4 | 0.631 |
