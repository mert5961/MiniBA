# MiniBA results

Scene: `desk.txt`: 150 images, 40137 points, 790261 observations, one SIMPLE_RADIAL camera 1280x720. Seed 123, 12 threads.

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
| polish (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0706 / 0.0855 | 1322.22 | -0.01124 | 5 | 3.97 |

A squared loss minimises RMS, so RMS can only go down; the solution barely moves.

## 3. Noisy start, clean matches

Every camera rotated by a random 0.5 deg (std), camera centres and all points moved by 1 % of the scene size (std per axis), focal length 5 % too long, distortion reset to 0.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 21.687 | 20.968 | 24.180 | 1.703 / 3.481 | 0.4956 / 1.6896 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 0.635 | 0.417 | 0.902 | 0.063 / 0.360 | 0.0713 / 0.0864 | 1322.17 | -0.01126 | 8 | 5.60 |

## 4. Noisy start + 5 % bad matches

39811 of 790261 observations moved 20-80 px in a random direction (wrong matches). Errors below are over the clean observations only, measured against their true pixels. Robust losses start down-weighting at 2 px.

| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |
|---|---|---|---|---|---|---|---|---|---|
| perturbed start | 21.687 | 20.968 | 24.180 | 1.703 / 3.481 | 0.4956 / 1.6896 | 1389.40 | 0.00000 | - | - |
| BA (squared) | 2.020 | 1.471 | 2.969 | 0.811 / 5.596 | 0.3551 / 0.8232 | 1325.71 | -0.00787 | 37 | 21.01 |
| BA (Huber) | 0.657 | 0.434 | 1.029 | 0.071 / 0.389 | 0.0573 / 0.0703 | 1322.61 | -0.01106 | 59 | 32.31 |
| BA (Cauchy) | 0.626 | 0.401 | 0.987 | 0.023 / 0.277 | 0.0190 / 0.0254 | 1323.14 | -0.01083 | 81 | 43.49 |

Flagging every observation above 4 px after the solve:

| run | precision | recall |
|---|---|---|
| BA (squared) | 0.328 | 0.999 |
| BA (Huber) | 0.944 | 0.998 |
| BA (Cauchy) | 0.910 | 0.998 |

## 5. Linear solvers (noisy start, clean matches, Huber)

Points are eliminated with the Schur complement; what differs is how the reduced camera system (902 unknowns) is solved.

| linear solver | iters | time s | s / iter | final cost | mean px |
|---|---|---|---|---|---|
| DENSE_SCHUR | 72 | 39.64 | 0.551 | 304484.4 | 0.631 |
| SPARSE_SCHUR | 71 | 38.06 | 0.536 | 304484.4 | 0.631 |
| ITERATIVE_SCHUR + SCHUR_JACOBI | 48 | 21.78 | 0.454 | 304484.4 | 0.631 |
