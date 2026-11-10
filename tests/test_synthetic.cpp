// Synthetic scenes with known answers. Run: ctest (or ./build/test_synthetic).
#include <cstdio>
#include <random>

#include "../src/ba.h"

using namespace miniba;

static int failures = 0;
#define EXPECT(cond, ...)                                  \
  do {                                                    \
    if (!(cond)) {                                        \
      ++failures;                                         \
      std::printf("FAIL %s:%d: ", __FILE__, __LINE__);    \
      std::printf(__VA_ARGS__);                           \
      std::printf("\n");                                  \
    }                                                     \
  } while (0)

// Cameras on a ring looking at the origin, points in a cube around it.
static Scene MakeScene(std::mt19937& rng, double pixel_noise) {
  Scene s;
  s.cam = {1280, 720, 1000.0, 640.0, 360.0, -0.02};
  const int n_cams = 12, n_pts = 400;
  for (int i = 0; i < n_cams; ++i) {
    const double a = 2 * M_PI * i / n_cams;
    const Eigen::Vector3d C(4 * std::cos(a), 0.5 * std::sin(3 * a), 4 * std::sin(a));
    const Eigen::Vector3d z = (-C).normalized();                           // look at origin
    const Eigen::Vector3d x = Eigen::Vector3d::UnitY().cross(z).normalized();
    const Eigen::Vector3d y = z.cross(x);
    Eigen::Matrix3d R;
    R.row(0) = x; R.row(1) = y; R.row(2) = z;                             // cam_from_world
    std::array<double, 6> pose{};
    SetRotation(pose, R);
    const Eigen::Vector3d t = -R * C;
    pose[3] = t.x(); pose[4] = t.y(); pose[5] = t.z();
    s.poses.push_back(pose);
    s.names.push_back("cam" + std::to_string(i));
  }
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::normal_distribution<double> n(0.0, pixel_noise);
  for (int j = 0; j < n_pts; ++j) s.points.push_back({u(rng), u(rng), u(rng)});
  const double fk[2] = {s.cam.f, s.cam.k};
  for (int i = 0; i < n_cams; ++i)
    for (int j = 0; j < n_pts; ++j) {
      double uv[2];
      if (!Project(s.poses[i].data(), fk, s.cam.cx, s.cam.cy, s.points[j].data(), uv)) continue;
      if (uv[0] < 0 || uv[1] < 0 || uv[0] >= s.cam.width || uv[1] >= s.cam.height) continue;
      s.obs.push_back({i, j, uv[0] + n(rng), uv[1] + n(rng)});
    }
  return s;
}

int main() {
  std::mt19937 rng(7);

  // A point on the optical axis projects to the principal point, whatever f and k are.
  {
    const std::array<double, 6> pose{0, 0, 0, 0, 0, 0};
    const double X[3] = {0, 0, 5}, fk[2] = {1234.0, 0.3};
    double uv[2];
    EXPECT(Project(pose.data(), fk, 640.0, 360.0, X, uv), "on-axis point should be in front");
    EXPECT(std::abs(uv[0] - 640) < 1e-12 && std::abs(uv[1] - 360) < 1e-12, "got %f %f", uv[0], uv[1]);
    const double behind[3] = {0, 0, -1};
    EXPECT(!Project(pose.data(), fk, 640.0, 360.0, behind, uv), "point behind the camera must be rejected");
  }

  // Hand-checked SIMPLE_RADIAL: u = 0.2, v = 0.1, r2 = 0.05, d = 1 + k r2.
  {
    const std::array<double, 6> pose{0, 0, 0, 0, 0, 0};
    const double X[3] = {0.4, 0.2, 2.0}, fk[2] = {1000.0, 0.1};
    double uv[2];
    Project(pose.data(), fk, 640.0, 360.0, X, uv);
    EXPECT(std::abs(uv[0] - (640 + 1000 * 1.005 * 0.2)) < 1e-9, "x %f", uv[0]);
    EXPECT(std::abs(uv[1] - (360 + 1000 * 1.005 * 0.1)) < 1e-9, "y %f", uv[1]);
  }

  // Noise-free scene: zero error at the truth; BA from a perturbed start returns to it.
  {
    const Scene truth = MakeScene(rng, 0.0);
    EXPECT(Stats(ReprojectionErrors(truth, truth.obs)).mean < 1e-9, "truth should reproject exactly");
    Scene s = truth;
    std::normal_distribution<double> n(0.0, 1.0);
    for (auto& p : s.poses) for (int c = 0; c < 6; ++c) p[c] += 0.01 * n(rng);
    for (auto& X : s.points) for (double& v : X) v += 0.02 * n(rng);
    s.cam.f *= 1.05;
    s.cam.k = 0.0;
    EXPECT(Stats(ReprojectionErrors(s, s.obs)).mean > 5.0, "perturbation should be visible");
    SolveOptions o;
    o.threads = 4;
    Solve(s, o);
    const double e = Stats(ReprojectionErrors(s, s.obs)).mean;
    const PoseError pe = ComparePoses(s, truth);
    EXPECT(e < 1e-4, "noise-free BA should reach ~0 px, got %g", e);
    EXPECT(pe.centre_max_pct < 0.01, "centres should match after alignment, max %g %%", pe.centre_max_pct);
    EXPECT(pe.rotation_max_deg < 1e-3, "rotations should match after alignment, max %g deg", pe.rotation_max_deg);
    EXPECT(std::abs(s.cam.f - truth.cam.f) / truth.cam.f < 1e-4, "focal %f vs %f", s.cam.f, truth.cam.f);
    EXPECT(std::abs(s.cam.k - truth.cam.k) < 1e-4, "k %f vs %f", s.cam.k, truth.cam.k);
  }

  // 0.5 px pixel noise: BA ends at the noise level (RMS ~ 0.5 * sqrt(2) per 2D residual), not below it.
  {
    const Scene truth = MakeScene(rng, 0.5);
    Scene s = truth;
    std::normal_distribution<double> n(0.0, 1.0);
    for (auto& X : s.points) for (double& v : X) v += 0.02 * n(rng);
    SolveOptions o;
    o.threads = 4;
    Solve(s, o);
    const ErrorStats st = Stats(ReprojectionErrors(s, s.obs));
    EXPECT(st.rms > 0.55 && st.rms < 0.75, "RMS should sit near the 0.71 px noise floor, got %g", st.rms);
  }

  // Robust loss: 10% gross outliers hurt the squared loss, Cauchy shrugs them off.
  {
    const Scene truth = MakeScene(rng, 0.3);
    Scene bad = truth;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<bool> inlier(bad.obs.size(), true);
    for (size_t i = 0; i < bad.obs.size(); ++i)
      if (u(rng) < 0.10) { bad.obs[i].x += 40 + 40 * u(rng); bad.obs[i].y -= 30; inlier[i] = false; }
    Scene sq = bad, ca = bad;
    SolveOptions o;
    o.threads = 4;
    Solve(sq, o);
    o.loss = Loss::kCauchy;
    Solve(ca, o);
    const double e_sq = Stats(ReprojectionErrors(sq, truth.obs), inlier).mean;
    const double e_ca = Stats(ReprojectionErrors(ca, truth.obs), inlier).mean;
    EXPECT(e_ca < 0.6, "Cauchy inlier error should stay near the noise, got %g", e_ca);
    EXPECT(e_sq > 2 * e_ca, "squared loss should be clearly worse with outliers (%g vs %g)", e_sq, e_ca);
  }

  std::printf(failures ? "%d check(s) failed\n" : "all checks passed\n", failures);
  return failures ? 1 : 0;
}
