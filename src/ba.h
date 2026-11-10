// MiniBA: a small bundle adjustment on a real COLMAP reconstruction, with Ceres.
// Header-only core shared by the CLI (miniba.cpp) and the synthetic test.
#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <memory>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace miniba {

// One shared SIMPLE_RADIAL camera, as in COLMAP: f, principal point (cx, cy), one radial term k.
// cx, cy stay fixed; f and k are optimised (a 2-double parameter block).
struct Intrinsics {
  int width = 0, height = 0;
  double f = 0, cx = 0, cy = 0, k = 0;
};

struct Observation {
  int image = 0, point = 0;
  double x = 0, y = 0;
};

struct Scene {
  Intrinsics cam;
  std::vector<std::string> names;
  std::vector<std::array<double, 6>> poses;   // cam_from_world: angle-axis (3) + translation (3)
  std::vector<std::array<double, 3>> points;  // world
  std::vector<Observation> obs;
};

inline Scene LoadScene(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::string magic, tag, model;
  int version = 0;
  in >> magic >> version;
  if (magic != "MINIBA" || version != 1) throw std::runtime_error("not a MINIBA v1 file");
  Scene s;
  in >> tag >> model >> s.cam.width >> s.cam.height >> s.cam.f >> s.cam.cx >> s.cam.cy >> s.cam.k;
  if (tag != "CAMERA" || model != "SIMPLE_RADIAL") throw std::runtime_error("expected CAMERA SIMPLE_RADIAL");
  size_t ni = 0, np = 0, no = 0;
  in >> ni >> np >> no;
  s.names.resize(ni);
  s.poses.resize(ni);
  for (size_t i = 0; i < ni; ++i) {
    in >> s.names[i];
    for (double& v : s.poses[i]) in >> v;
  }
  s.points.resize(np);
  for (auto& p : s.points) in >> p[0] >> p[1] >> p[2];
  s.obs.resize(no);
  for (auto& o : s.obs) in >> o.image >> o.point >> o.x >> o.y;
  if (!in) throw std::runtime_error("truncated file " + path);
  return s;
}

// World point -> pixel, exactly like COLMAP's SIMPLE_RADIAL model.
// Returns false when the point is behind the camera.
template <typename T>
bool Project(const T* pose, const T* fk, double cx, double cy, const T* X, T* uv) {
  T p[3];
  ceres::AngleAxisRotatePoint(pose, X, p);
  p[0] += pose[3];
  p[1] += pose[4];
  p[2] += pose[5];
  if (!(p[2] > T(1e-12))) return false;
  const T u = p[0] / p[2], v = p[1] / p[2];
  const T d = T(1.0) + fk[1] * (u * u + v * v);
  uv[0] = fk[0] * d * u + T(cx);
  uv[1] = fk[0] * d * v + T(cy);
  return true;
}

// Ceres residual: projected minus observed pixel (2 residuals).
// Parameter blocks: pose (6), intrinsics {f, k} (2), point (3).
struct ReprojectionError {
  double x, y, cx, cy;
  template <typename T>
  bool operator()(const T* pose, const T* fk, const T* X, T* r) const {
    T uv[2];
    if (!Project(pose, fk, cx, cy, X, uv)) {
      r[0] = r[1] = T(0);
      return false;  // tells Ceres this step is invalid, it shrinks the trust region
    }
    r[0] = uv[0] - T(x);
    r[1] = uv[1] - T(y);
    return true;
  }
  static ceres::CostFunction* Create(const Observation& o, double cx, double cy) {
    return new ceres::AutoDiffCostFunction<ReprojectionError, 2, 6, 2, 3>(
        new ReprojectionError{o.x, o.y, cx, cy});
  }
};

// Per-observation reprojection error in pixels (NaN if behind the camera).
inline std::vector<double> ReprojectionErrors(const Scene& s, const std::vector<Observation>& obs) {
  const double fk[2] = {s.cam.f, s.cam.k};
  std::vector<double> e(obs.size());
  for (size_t i = 0; i < obs.size(); ++i) {
    const auto& o = obs[i];
    double uv[2];
    if (!Project(s.poses[o.image].data(), fk, s.cam.cx, s.cam.cy, s.points[o.point].data(), uv)) {
      e[i] = std::nan("");
      continue;
    }
    e[i] = std::hypot(uv[0] - o.x, uv[1] - o.y);
  }
  return e;
}

struct ErrorStats {
  double mean = 0, median = 0, rms = 0, p95 = 0;
  size_t n = 0, behind = 0;
};

// Statistics over the observations where mask[i] is true (all when mask is empty).
inline ErrorStats Stats(const std::vector<double>& e, const std::vector<bool>& mask = {}) {
  std::vector<double> v;
  v.reserve(e.size());
  ErrorStats st;
  for (size_t i = 0; i < e.size(); ++i) {
    if (!mask.empty() && !mask[i]) continue;
    if (std::isnan(e[i])) { ++st.behind; continue; }
    v.push_back(e[i]);
  }
  st.n = v.size();
  if (v.empty()) return st;
  double sum = 0, sq = 0;
  for (double x : v) { sum += x; sq += x * x; }
  st.mean = sum / v.size();
  st.rms = std::sqrt(sq / v.size());
  std::sort(v.begin(), v.end());
  st.median = v[v.size() / 2];
  st.p95 = v[static_cast<size_t>(0.95 * (v.size() - 1))];
  return st;
}

inline Eigen::Matrix3d RotationOf(const std::array<double, 6>& pose) {
  const Eigen::Vector3d aa(pose[0], pose[1], pose[2]);
  const double angle = aa.norm();
  if (angle < 1e-15) return Eigen::Matrix3d::Identity();
  return Eigen::AngleAxisd(angle, aa / angle).toRotationMatrix();
}

inline void SetRotation(std::array<double, 6>& pose, const Eigen::Matrix3d& R) {
  const Eigen::AngleAxisd aa(R);
  const Eigen::Vector3d v = aa.axis() * aa.angle();
  pose[0] = v.x(); pose[1] = v.y(); pose[2] = v.z();
}

// Camera centre in world coordinates: C = -R^T t.
inline Eigen::Vector3d CentreOf(const std::array<double, 6>& pose) {
  return -RotationOf(pose).transpose() * Eigen::Vector3d(pose[3], pose[4], pose[5]);
}

// Pose accuracy against a reference, after removing the 7-DoF gauge (scale, rotation, translation)
// with an Umeyama similarity fit on the camera centres. BA alone cannot know the world frame.
struct PoseError {
  double centre_mean_pct = 0;  // mean centre error, % of the scene size (median camera distance to centroid)
  double centre_max_pct = 0;
  double rotation_mean_deg = 0;
  double rotation_max_deg = 0;
};

inline double SceneScale(const Scene& s) {
  Eigen::Vector3d c = Eigen::Vector3d::Zero();
  for (const auto& p : s.poses) c += CentreOf(p);
  c /= double(s.poses.size());
  std::vector<double> d;
  for (const auto& p : s.poses) d.push_back((CentreOf(p) - c).norm());
  std::sort(d.begin(), d.end());
  return d[d.size() / 2];
}

inline PoseError ComparePoses(const Scene& est, const Scene& ref) {
  const int n = static_cast<int>(est.poses.size());
  Eigen::Matrix3Xd A(3, n), B(3, n);
  for (int i = 0; i < n; ++i) {
    A.col(i) = CentreOf(est.poses[i]);
    B.col(i) = CentreOf(ref.poses[i]);
  }
  const Eigen::Matrix4d T = Eigen::umeyama(A, B, true);  // B ~ s R A + t
  const Eigen::Matrix3d sR = T.block<3, 3>(0, 0);
  const double scale = std::cbrt(sR.determinant());
  const Eigen::Matrix3d Ralign = sR / scale;
  const double ref_scale = SceneScale(ref);
  PoseError pe;
  for (int i = 0; i < n; ++i) {
    const Eigen::Vector3d a = sR * A.col(i) + T.block<3, 1>(0, 3);
    const double dc = 100.0 * (a - B.col(i)).norm() / ref_scale;
    // cam_from_world rotation expressed in the reference world: R_est * Ralign^T
    const Eigen::Matrix3d Rd = RotationOf(ref.poses[i]).transpose() * (RotationOf(est.poses[i]) * Ralign.transpose());
    const double dr = Eigen::AngleAxisd(Rd).angle() * 180.0 / M_PI;
    pe.centre_mean_pct += dc / n;
    pe.rotation_mean_deg += dr / n;
    pe.centre_max_pct = std::max(pe.centre_max_pct, dc);
    pe.rotation_max_deg = std::max(pe.rotation_max_deg, dr);
  }
  return pe;
}

enum class Loss { kSquared, kHuber, kCauchy };

inline const char* LossName(Loss l) {
  switch (l) {
    case Loss::kSquared: return "squared";
    case Loss::kHuber: return "Huber";
    case Loss::kCauchy: return "Cauchy";
  }
  return "?";
}

struct SolveOptions {
  Loss loss = Loss::kSquared;
  double loss_scale_px = 2.0;  // where the robust loss starts to down-weight a residual
  ceres::LinearSolverType linear_solver = ceres::SPARSE_SCHUR;
  ceres::PreconditionerType preconditioner = ceres::SCHUR_JACOBI;
  bool refine_intrinsics = true;
  int max_iterations = 200;
  int threads = 1;
  bool verbose = false;
};

struct SolveReport {
  double initial_cost = 0, final_cost = 0, seconds = 0;
  int iterations = 0;
  bool converged = false;
  std::string termination;
};

// Bundle adjustment in place: poses, points and (optionally) f and k.
inline SolveReport Solve(Scene& s, const SolveOptions& opt) {
  ceres::Problem::Options popt;
  popt.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  ceres::Problem problem(popt);
  std::unique_ptr<ceres::LossFunction> loss;
  if (opt.loss == Loss::kHuber) loss = std::make_unique<ceres::HuberLoss>(opt.loss_scale_px);
  if (opt.loss == Loss::kCauchy) loss = std::make_unique<ceres::CauchyLoss>(opt.loss_scale_px);

  double fk[2] = {s.cam.f, s.cam.k};
  for (const auto& o : s.obs) {
    problem.AddResidualBlock(ReprojectionError::Create(o, s.cam.cx, s.cam.cy), loss.get(),
                             s.poses[o.image].data(), fk, s.points[o.point].data());
  }
  if (!opt.refine_intrinsics) problem.SetParameterBlockConstant(fk);
  // No pose is held fixed: the 7-DoF gauge stays free and Levenberg-Marquardt damping keeps it
  // well behaved. Accuracy is measured after a similarity alignment (ComparePoses).

  ceres::Solver::Options so;
  so.linear_solver_type = opt.linear_solver;
  so.preconditioner_type = opt.preconditioner;
  so.minimizer_progress_to_stdout = opt.verbose;
  so.max_num_iterations = opt.max_iterations;
  so.num_threads = opt.threads;
  so.function_tolerance = 1e-8;
  ceres::Solver::Summary sum;
  ceres::Solve(so, &problem, &sum);
  s.cam.f = fk[0];
  s.cam.k = fk[1];

  SolveReport r;
  r.initial_cost = sum.initial_cost;
  r.final_cost = sum.final_cost;
  r.seconds = sum.total_time_in_seconds;
  r.iterations = sum.num_successful_steps + sum.num_unsuccessful_steps;
  r.converged = sum.termination_type == ceres::CONVERGENCE;
  r.termination = ceres::TerminationTypeToString(sum.termination_type);
  return r;
}

}  // namespace miniba
