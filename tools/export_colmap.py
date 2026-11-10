"""Export a COLMAP sparse model (one SIMPLE_RADIAL camera) to MiniBA's plain-text format.

Usage: python tools/export_colmap.py <colmap_model_dir> <out.txt>

Format (all numbers in full precision):
    MINIBA 1
    CAMERA SIMPLE_RADIAL <width> <height> <f> <cx> <cy> <k>
    <n_images> <n_points> <n_observations>
    <n_images lines>        name rx ry rz tx ty tz      (cam_from_world, angle-axis + translation)
    <n_points lines>        X Y Z                       (world)
    <n_observations lines>  image_index point_index x y (pixels, COLMAP convention)
"""
import sys

import numpy as np
import pycolmap
from scipy.spatial.transform import Rotation


def cam_from_world(image):
    pose = image.cam_from_world
    return pose() if callable(pose) else pose  # method in pycolmap >= 3.12, property before


def main(model_dir: str, out_path: str) -> None:
    rec = pycolmap.Reconstruction(model_dir)
    if len(rec.cameras) != 1:
        sys.exit(f"expected one shared camera, found {len(rec.cameras)}")
    cam = next(iter(rec.cameras.values()))
    if cam.model.name != "SIMPLE_RADIAL":
        sys.exit(f"only SIMPLE_RADIAL is supported, found {cam.model.name}")

    images = sorted((im for im in rec.images.values() if im.has_pose), key=lambda im: im.name)
    image_index = {im.image_id: i for i, im in enumerate(images)}
    point_ids = sorted(rec.points3D.keys())
    point_index = {pid: j for j, pid in enumerate(point_ids)}

    obs = []
    for pid in point_ids:
        for el in rec.points3D[pid].track.elements:
            if el.image_id not in image_index:
                continue
            xy = rec.images[el.image_id].points2D[el.point2D_idx].xy
            obs.append((image_index[el.image_id], point_index[pid], xy[0], xy[1]))

    with open(out_path, "w") as f:
        f.write("MINIBA 1\n")
        f.write("CAMERA SIMPLE_RADIAL %d %d %s\n" % (cam.width, cam.height,
                                                    " ".join("%.17g" % p for p in cam.params)))
        f.write("%d %d %d\n" % (len(images), len(point_ids), len(obs)))
        for im in images:
            pose = cam_from_world(im)
            rvec = Rotation.from_matrix(np.asarray(pose.rotation.matrix())).as_rotvec()
            t = np.asarray(pose.translation)
            f.write("%s %s\n" % (im.name, " ".join("%.17g" % v for v in (*rvec, *t))))
        for pid in point_ids:
            f.write("%s\n" % " ".join("%.17g" % v for v in rec.points3D[pid].xyz))
        for i, j, x, y in obs:
            f.write("%d %d %.17g %.17g\n" % (i, j, x, y))

    print(f"{len(images)} images, {len(point_ids)} points, {len(obs)} observations -> {out_path}")
    print(f"COLMAP mean reprojection error: {rec.compute_mean_reprojection_error():.6f} px")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
