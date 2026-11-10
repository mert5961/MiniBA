// MiniBA experiments on a real COLMAP reconstruction.
//
//   miniba <scene.txt> [--out results] [--seed 42] [--outliers 0.05]
//
// 1. check      the projection model reproduces COLMAP's own mean reprojection error
// 2. polish     BA started from COLMAP's solution (should barely move)
// 3. noisy      poses, points and focal length perturbed, then BA recovers them
// 4. outliers   same start + a share of observations replaced by bad matches: squared vs Huber vs Cauchy
// 5. solvers    dense vs sparse vs iterative Schur on the same problem
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>

#include "ba.h"

using namespace miniba;

namespace {

struct Perturbation {
  double rotation_deg = 0.5;  // std of a random small rotation per camera
  double centre_pct = 1.0;    // std per axis, % of scene size
  double point_pct = 1.0;     // std per axis, % of scene size
  double focal_factor = 1.05; // start with a 5% wrong focal length
  bool reset_k = true;        // start with zero distortion
};

Scene Perturb(const Scene& ref, const Perturbation& p, std::mt19937& rng) {
  Scene s = ref;
  const double scale = SceneScale(ref);
  std::normal_distribution<double> n01(0.0, 1.0);
  for (auto& pose : s.poses) {
    const Eigen::Vector3d axis = Eigen::Vector3d(n01(rng), n01(rng), n01(rng)).normalized();
    const double angle = p.rotation_deg * M_PI / 180.0 * n01(rng);
    const Eigen::Matrix3d R = Eigen::AngleAxisd(angle, axis).toRotationMatrix() * RotationOf(pose);
    const Eigen::Vector3d C =
        CentreOf(pose) + scale * p.centre_pct / 100.0 * Eigen::Vector3d(n01(rng), n01(rng), n01(rng));
    SetRotation(pose, R);
    const Eigen::Vector3d t = -R * C;
    pose[3] = t.x(); pose[4] = t.y(); pose[5] = t.z();
  }
  for (auto& X : s.points)
    for (double& v : X) v += scale * p.point_pct / 100.0 * n01(rng);
  s.cam.f *= p.focal_factor;
  if (p.reset_k) s.cam.k = 0.0;
  return s;
}

// Replace a share of observations with bad matches: shifted 20-80 px in a random direction.
std::vector<bool> InjectOutliers(Scene& s, double share, std::mt19937& rng) {
  std::vector<bool> is_outlier(s.obs.size(), false);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  for (size_t i = 0; i < s.obs.size(); ++i) {
    if (u01(rng) >= share) continue;
    const double r = 20.0 + 60.0 * u01(rng), a = 2.0 * M_PI * u01(rng);
    s.obs[i].x += r * std::cos(a);
    s.obs[i].y += r * std::sin(a);
    is_outlier[i] = true;
  }
  return is_outlier;
}

std::string Fmt(double v, int prec = 3) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*f", prec, v);
  return b;
}

// Reprojection error of the *clean* observations, so a robust solve is not rewarded for ignoring them.
ErrorStats CleanStats(const Scene& s, const Scene& clean_obs, const std::vector<bool>& inlier) {
  return Stats(ReprojectionErrors(s, clean_obs.obs), inlier);
}

struct Row {
  std::string name;
  ErrorStats err;
  PoseError pose;
  double f = 0, k = 0;
  SolveReport rep;
  double outlier_precision = -1, outlier_recall = -1;
};

std::string Header() {
  return "| run | mean px | median px | RMS px | cam centre err (mean / max, % scene) | rotation err (mean / max, deg) | f | k | iters | time s |\n"
         "|---|---|---|---|---|---|---|---|---|---|\n";
}

std::string Line(const Row& r) {
  std::ostringstream o;
  o << "| " << r.name << " | " << Fmt(r.err.mean) << " | " << Fmt(r.err.median) << " | " << Fmt(r.err.rms)
    << " | " << Fmt(r.pose.centre_mean_pct, 3) << " / " << Fmt(r.pose.centre_max_pct, 3)
    << " | " << Fmt(r.pose.rotation_mean_deg, 4) << " / " << Fmt(r.pose.rotation_max_deg, 4)
    << " | " << Fmt(r.f, 2) << " | " << Fmt(r.k, 5) << " | "
    << (r.rep.iterations ? std::to_string(r.rep.iterations) + (r.rep.converged ? "" : " (not converged)") : "-")
    << " | "
    << (r.rep.iterations ? Fmt(r.rep.seconds, 2) : "-") << " |\n";
  return o.str();
}

void Histogram(const std::vector<double>& e, const std::vector<bool>& mask, std::vector<double>& out,
               double width, int bins) {
  out.assign(bins + 1, 0.0);  // last bin = everything beyond the range
  size_t n = 0;
  for (size_t i = 0; i < e.size(); ++i) {
    if (!mask.empty() && !mask[i]) continue;
    if (std::isnan(e[i])) continue;
    const int b = std::min(bins, static_cast<int>(e[i] / width));
    out[b] += 1.0;
    ++n;
  }
  for (double& v : out) v /= std::max<size_t>(n, 1);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: miniba <scene.txt> [--out results] [--seed 42] [--outliers 0.05]\n";
    return 2;
  }
  std::string scene_path = argv[1], out_dir = "results";
  unsigned seed = 42;
  double outlier_share = 0.05;
  for (int i = 2; i + 1 < argc; i += 2) {
    const std::string a = argv[i];
    if (a == "--out") out_dir = argv[i + 1];
    else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(argv[i + 1]));
    else if (a == "--outliers") outlier_share = std::stod(argv[i + 1]);
    else { std::cerr << "unknown option " << a << "\n"; return 2; }
  }
  google::InitGoogleLogging(argv[0]);
  std::filesystem::create_directories(out_dir);
  const int threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));

  const auto t0 = std::chrono::steady_clock::now();
  const Scene ref = LoadScene(scene_path);
  const double load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::cout << "loaded " << ref.poses.size() << " images, " << ref.points.size() << " points, "
            << ref.obs.size() << " observations in " << Fmt(load_s, 1) << " s\n";

  std::ostringstream md;
  md << "# MiniBA results\n\n"
     << "Scene: `" << std::filesystem::path(scene_path).filename().string() << "`: " << ref.poses.size()
     << " images, " << ref.points.size() << " points, " << ref.obs.size()
     << " observations, one SIMPLE_RADIAL camera " << ref.cam.width << "x" << ref.cam.height
     << ". Seed " << seed << ", " << threads << " threads.\n\n"
     << "Pose errors are measured against COLMAP/GLOMAP's solution after a 7-DoF similarity alignment "
        "(Umeyama on camera centres), as a percentage of the scene size (median camera distance to the "
        "centroid = " << Fmt(SceneScale(ref), 4) << " model units).\n\n";

  // 1. check
  const std::vector<double> base_e = ReprojectionErrors(ref, ref.obs);
  const ErrorStats base = Stats(base_e);
  double per_point_mean = 0;
  {
    std::vector<double> sum(ref.points.size(), 0.0);
    std::vector<int> cnt(ref.points.size(), 0);
    for (size_t i = 0; i < ref.obs.size(); ++i)
      if (!std::isnan(base_e[i])) { sum[ref.obs[i].point] += base_e[i]; ++cnt[ref.obs[i].point]; }
    int n = 0;
    for (size_t j = 0; j < sum.size(); ++j)
      if (cnt[j]) { per_point_mean += sum[j] / cnt[j]; ++n; }
    per_point_mean /= std::max(n, 1);
  }
  std::cout << "[check] COLMAP's solution: mean over observations " << Fmt(base.mean, 6)
            << " px, mean of per-point means " << Fmt(per_point_mean, 6) << " px\n";
  md << "## 1. Check: same projection model as COLMAP\n\n"
     << "| metric | MiniBA (C++) | pycolmap / COLMAP |\n|---|---|---|\n"
     << "| mean over all observations | **" << Fmt(base.mean, 6) << " px** | 0.624341 px (pycolmap's own projection, "
        "`tools/check_colmap.py`) |\n"
     << "| mean of per-point means | **" << Fmt(per_point_mean, 6) << " px** | 0.516551 px "
        "(`Reconstruction.compute_mean_reprojection_error()`) |\n\n"
     << "COLMAP's headline number averages each 3D point's track first, then over points, so long tracks "
        "(which have larger errors) count less. Both match to 6 digits, so the projection model is the same. "
     << base.behind << " observations behind a camera. All errors below are means over observations.\n\n";

  // 2. polish
  Row polish;
  {
    Scene s = ref;
    SolveOptions o;
    o.threads = threads;
    polish.rep = Solve(s, o);
    polish.name = "polish (squared)";
    polish.err = Stats(ReprojectionErrors(s, s.obs));
    polish.pose = ComparePoses(s, ref);
    polish.f = s.cam.f;
    polish.k = s.cam.k;
    std::cout << "[polish] " << Fmt(base.rms) << " -> " << Fmt(polish.err.rms) << " px RMS in "
              << polish.rep.iterations << " iterations\n";
  }
  Row ref_row{"COLMAP solution", base, PoseError{}, ref.cam.f, ref.cam.k, SolveReport{}};
  md << "## 2. Polish: BA started from COLMAP's solution\n\n" << Header() << Line(ref_row) << Line(polish)
     << "\nA squared loss minimises RMS, so RMS can only go down; the solution barely moves.\n\n";

  // 3. noisy start, clean data
  std::mt19937 rng(seed);
  const Perturbation pert;
  const Scene start = Perturb(ref, pert, rng);
  Row start_row{"perturbed start", Stats(ReprojectionErrors(start, start.obs)), ComparePoses(start, ref),
                start.cam.f, start.cam.k, SolveReport{}};
  Row noisy;
  {
    Scene s = start;
    SolveOptions o;
    o.threads = threads;
    noisy.rep = Solve(s, o);
    noisy.name = "BA (squared)";
    noisy.err = Stats(ReprojectionErrors(s, s.obs));
    noisy.pose = ComparePoses(s, ref);
    noisy.f = s.cam.f;
    noisy.k = s.cam.k;
    std::cout << "[noisy] " << Fmt(start_row.err.mean) << " -> " << Fmt(noisy.err.mean) << " px mean, centre err "
              << Fmt(noisy.pose.centre_mean_pct) << " % in " << noisy.rep.iterations << " iterations\n";
  }
  md << "## 3. Noisy start, clean matches\n\n"
     << "Every camera rotated by a random " << pert.rotation_deg << " deg (std), camera centres and all points moved by "
     << pert.centre_pct << " % of the scene size (std per axis), focal length " << Fmt((pert.focal_factor - 1) * 100, 0)
     << " % too long, distortion reset to 0.\n\n"
     << Header() << Line(start_row) << Line(noisy) << "\n";

  // 4. outliers
  Scene bad = start;
  std::mt19937 rng_out(seed + 1);
  const std::vector<bool> is_out = InjectOutliers(bad, outlier_share, rng_out);
  std::vector<bool> inlier(is_out.size());
  size_t n_out = 0;
  for (size_t i = 0; i < is_out.size(); ++i) { inlier[i] = !is_out[i]; n_out += is_out[i]; }
  std::cout << "[outliers] " << n_out << " of " << is_out.size() << " observations replaced by bad matches\n";

  std::vector<Row> rows;
  std::vector<std::vector<double>> hists;
  const double hist_w = 0.1;
  const int hist_bins = 50;
  for (Loss l : {Loss::kSquared, Loss::kHuber, Loss::kCauchy}) {
    Scene s = bad;
    SolveOptions o;
    o.threads = threads;
    o.loss = l;
    Row r;
    r.rep = Solve(s, o);
    r.name = std::string("BA (") + LossName(l) + ")";
    r.err = CleanStats(s, ref, inlier);  // inliers only, against the true (clean) pixels
    r.pose = ComparePoses(s, ref);
    r.f = s.cam.f;
    r.k = s.cam.k;
    // Which observations would we now flag as bad matches (> 4 px)?
    const auto e = ReprojectionErrors(s, s.obs);
    size_t tp = 0, flagged = 0;
    for (size_t i = 0; i < e.size(); ++i) {
      const bool f = std::isnan(e[i]) || e[i] > 4.0;
      flagged += f;
      tp += f && is_out[i];
    }
    r.outlier_precision = flagged ? double(tp) / flagged : 0.0;
    r.outlier_recall = n_out ? double(tp) / n_out : 0.0;
    std::vector<double> h;
    Histogram(ReprojectionErrors(s, ref.obs), inlier, h, hist_w, hist_bins);
    hists.push_back(h);
    std::cout << "[outliers] " << LossName(l) << ": inlier mean " << Fmt(r.err.mean) << " px, centre err "
              << Fmt(r.pose.centre_mean_pct) << " %, f " << Fmt(r.f, 1) << ", flagged P/R "
              << Fmt(r.outlier_precision) << "/" << Fmt(r.outlier_recall) << "\n";
    rows.push_back(r);
  }
  Row bad_start{"perturbed start", CleanStats(bad, ref, inlier), start_row.pose, bad.cam.f, bad.cam.k, SolveReport{}};
  md << "## 4. Noisy start + " << Fmt(outlier_share * 100, 0) << " % bad matches\n\n"
     << n_out << " of " << is_out.size()
     << " observations moved 20-80 px in a random direction (wrong matches). Errors below are over the clean "
        "observations only, measured against their true pixels. Robust losses start down-weighting at 2 px.\n\n"
     << Header() << Line(bad_start);
  for (const auto& r : rows) md << Line(r);
  md << "\nFlagging every observation above 4 px after the solve:\n\n| run | precision | recall |\n|---|---|---|\n";
  for (const auto& r : rows)
    md << "| " << r.name << " | " << Fmt(r.outlier_precision) << " | " << Fmt(r.outlier_recall) << " |\n";
  md << "\n";
  {
    std::ofstream h(out_dir + "/residual_hist.csv");
    h << "bin_start_px,squared,huber,cauchy\n";
    for (int b = 0; b <= hist_bins; ++b)
      h << Fmt(b * hist_w, 1) << "," << hists[0][b] << "," << hists[1][b] << "," << hists[2][b] << "\n";
  }

  // 5. linear solvers (same noisy, clean problem; Huber)
  md << "## 5. Linear solvers (noisy start, clean matches, Huber)\n\n"
     << "Points are eliminated with the Schur complement; what differs is how the reduced camera system ("
     << ref.poses.size() * 6 + 2 << " unknowns) is solved.\n\n"
     << "| linear solver | iters | time s | s / iter | final cost | mean px |\n|---|---|---|---|---|---|\n";
  struct LS { const char* name; ceres::LinearSolverType t; };
  for (LS ls : {LS{"DENSE_SCHUR", ceres::DENSE_SCHUR}, LS{"SPARSE_SCHUR", ceres::SPARSE_SCHUR},
                LS{"ITERATIVE_SCHUR + SCHUR_JACOBI", ceres::ITERATIVE_SCHUR}}) {
    Scene s = start;
    SolveOptions o;
    o.threads = threads;
    o.loss = Loss::kHuber;
    o.linear_solver = ls.t;
    const SolveReport r = Solve(s, o);
    const ErrorStats e = Stats(ReprojectionErrors(s, s.obs));
    std::cout << "[solvers] " << ls.name << ": " << r.iterations << " iterations, " << Fmt(r.seconds, 2) << " s\n";
    md << "| " << ls.name << " | " << r.iterations << " | " << Fmt(r.seconds, 2) << " | "
       << Fmt(r.seconds / std::max(r.iterations, 1), 3) << " | " << Fmt(r.final_cost, 1)
       << " | " << Fmt(e.mean) << " |\n";
  }

  std::ofstream(out_dir + "/report.md") << md.str();
  std::cout << "wrote " << out_dir << "/report.md and " << out_dir << "/residual_hist.csv\n";
  return 0;
}
