"""Reprojection error of a COLMAP model computed with pycolmap's own projection, two ways.

Usage: python tools/check_colmap.py <colmap_model_dir>

COLMAP's compute_mean_reprojection_error() is the mean of per-point mean errors. MiniBA reports the mean over
all observations. This prints both, so MiniBA's projection model can be checked against COLMAP's.
"""
import sys

import numpy as np
import pycolmap


def main(model_dir: str) -> None:
    rec = pycolmap.Reconstruction(model_dir)
    per_obs, per_point = [], []
    for point in rec.points3D.values():
        errs = []
        for el in point.track.elements:
            image = rec.images[el.image_id]
            pose = image.cam_from_world
            pose = pose() if callable(pose) else pose
            uv = np.asarray(image.camera.img_from_cam(pose * point.xyz)).ravel()[:2]
            errs.append(np.linalg.norm(uv - image.points2D[el.point2D_idx].xy))
        per_obs += errs
        per_point.append(np.mean(errs))
    print(f"mean over observations:  {np.mean(per_obs):.6f} px  ({len(per_obs)} observations)")
    print(f"mean of per-point means: {np.mean(per_point):.6f} px  ({len(per_point)} points)")
    print(f"COLMAP's own number:     {rec.compute_mean_reprojection_error():.6f} px")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
