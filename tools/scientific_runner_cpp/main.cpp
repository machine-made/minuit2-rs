#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using ROOT::Minuit2::FCNBase;
using ROOT::Minuit2::FunctionMinimum;
using ROOT::Minuit2::MnHesse;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnStrategy;
using ROOT::Minuit2::MnUserParameters;

namespace {

enum class RunMode
{
   Full,
   LoadOnly,
   SolveOnly
};

std::string trim(std::string const &s)
{
   std::size_t b = s.find_first_not_of(" \t\r\n");
   if (b == std::string::npos) {
      return {};
   }
   std::size_t e = s.find_last_not_of(" \t\r\n");
   return s.substr(b, e - b + 1);
}

std::vector<std::string> split_csv(std::string const &line)
{
   std::vector<std::string> out;
   std::stringstream ss(line);
   std::string tok;
   while (std::getline(ss, tok, ',')) {
      out.push_back(trim(tok));
   }
   return out;
}

bool parse_double(std::string s, double &out)
{
   if (s.empty()) {
      return false;
   }
   if (s.front() == '.') {
      s = "0" + s;
   } else if (s.rfind("-.", 0) == 0) {
      s = "-0" + s.substr(1);
   }
   try {
      std::size_t idx = 0;
      out = std::stod(s, &idx);
      return idx > 0;
   } catch (...) {
      return false;
   }
}

std::vector<double> parse_floats(std::string line)
{
   for (char &c : line) {
      if (c == ',' || c == '=' || c == ':' || c == '(' || c == ')' || c == '[' || c == ']' || c == '*') {
         c = ' ';
      }
   }
   std::stringstream ss(line);
   std::vector<double> vals;
   std::string tok;
   while (ss >> tok) {
      double v = 0.0;
      if (parse_double(tok, v)) {
         vals.push_back(v);
      }
   }
   return vals;
}

struct LeastSquaresFCN final : public FCNBase {
   using ModelFn = std::function<double(std::vector<double> const &, double)>;
   std::vector<double> x;
   std::vector<double> y;
   ModelFn model;

   double operator()(std::vector<double> const &p) const override
   {
      double rss = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i) {
         double pred = model(p, x[i]);
         if (!std::isfinite(pred)) {
            return 1e30;
         }
         double r = y[i] - pred;
         rss += r * r;
      }
      return rss;
   }

   double Up() const override { return 1.0; }
};

struct NoaaFCN final : public FCNBase {
   std::vector<double> t;
   std::vector<double> y;
   std::vector<double> sigma;

   static double model(std::vector<double> const &p, double x)
   {
      double const w1 = 2.0 * M_PI * x;
      double const w2 = 4.0 * M_PI * x;
      return p[0] + p[1] * x + p[2] * x * x + p[3] * std::sin(w1) + p[4] * std::cos(w1) + p[5] * std::sin(w2) +
             p[6] * std::cos(w2) + p[7] * x * std::sin(w1);
   }

   double operator()(std::vector<double> const &p) const override
   {
      double chi2 = 0.0;
      for (std::size_t i = 0; i < t.size(); ++i) {
         double pred = model(p, t[i]);
         if (!std::isfinite(pred)) {
            return 1e30;
         }
         double r = (y[i] - pred) / sigma[i];
         chi2 += r * r;
      }
      return chi2;
   }

   double Up() const override { return 1.0; }
};

struct HistFCN final : public FCNBase {
   std::vector<double> x;
   std::vector<double> y;
   std::vector<double> sigma;

   static double model(std::vector<double> const &p, double x)
   {
      double amp = p[0];
      double mu = p[1];
      double sig = p[2];
      double c0 = p[3];
      double c1 = p[4];
      if (sig <= 0.05) {
         return std::numeric_limits<double>::quiet_NaN();
      }
      double z = (x - mu) / sig;
      double peak = amp * std::exp(-0.5 * z * z);
      double bg = c0 + c1 * (x - 91.0);
      return std::max(peak + bg, 1e-9);
   }

   double operator()(std::vector<double> const &p) const override
   {
      double chi2 = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i) {
         double pred = model(p, x[i]);
         if (!std::isfinite(pred)) {
            return 1e30;
         }
         double r = (y[i] - pred) / sigma[i];
         chi2 += r * r;
      }
      return chi2;
   }

   double Up() const override { return 1.0; }
};

struct NistDataset {
   std::vector<double> x;
   std::vector<double> y;
   std::vector<double> start1;
   std::vector<double> start2;
   std::vector<double> certified;
};

bool parse_noaa(std::string const &path, NoaaFCN &fcn)
{
   std::ifstream in(path);
   if (!in) {
      return false;
   }
   std::string line;
   bool first = true;
   double t0 = 0.0;
   while (std::getline(in, line)) {
      line = trim(line);
      if (line.empty() || line[0] == '#') {
         continue;
      }
      if (line.rfind("year,month,decimal date", 0) == 0) {
         continue;
      }
      auto cols = split_csv(line);
      if (cols.size() < 8) {
         continue;
      }
      double td = 0.0, y = 0.0, unc = 0.0;
      if (!parse_double(cols[2], td) || !parse_double(cols[3], y) || !parse_double(cols[7], unc)) {
         continue;
      }
      if (unc <= 0.0) {
         continue;
      }
      if (first) {
         t0 = td;
         first = false;
      }
      fcn.t.push_back(td - t0);
      fcn.y.push_back(y);
      fcn.sigma.push_back(std::max(unc, 1e-6));
   }
   return !fcn.t.empty();
}

bool parse_usgs_magnitudes(std::string const &path, std::vector<double> &mags)
{
   std::ifstream in(path);
   if (!in) {
      return false;
   }
   std::string line;
   bool header = true;
   while (std::getline(in, line)) {
      if (header) {
         header = false;
         continue;
      }
      auto cols = split_csv(line);
      if (cols.size() < 5) {
         continue;
      }
      double m = 0.0;
      if (parse_double(cols[4], m) && std::isfinite(m)) {
         mags.push_back(m);
      }
   }
   return !mags.empty();
}

void build_cumulative(std::vector<double> const &mags, double mmin, double mmax, double dm, std::vector<double> &mvals,
                      std::vector<double> &counts)
{
   for (double m = mmin; m <= mmax + 1e-12; m += dm) {
      double n = 0.0;
      for (double v : mags) {
         if (v >= m) {
            n += 1.0;
         }
      }
      if (n > 0.0) {
         mvals.push_back(m);
         counts.push_back(n);
      }
   }
}

bool parse_nist_dat(std::string const &path, std::size_t nparam, NistDataset &out)
{
   std::ifstream in(path);
   if (!in) {
      return false;
   }
   std::string line;
   bool in_data = false;
   while (std::getline(in, line)) {
      std::string t = trim(line);
      if (t.empty()) {
         continue;
      }
      if (auto pos = t.find('='); pos != std::string::npos) {
         std::string lhs = trim(t.substr(0, pos));
         if (!lhs.empty() && lhs[0] == 'b') {
            auto nums = parse_floats(t.substr(pos + 1));
            if (nums.size() >= 4) {
               out.start1.push_back(nums[0]);
               out.start2.push_back(nums[1]);
               out.certified.push_back(nums[2]);
            }
         }
      }
      if (t.rfind("Data:", 0) == 0) {
         std::string tail = trim(t.substr(5));
         if (!tail.empty() && tail[0] == 'y') {
            in_data = true;
            continue;
         }
      }
      if (in_data) {
         auto nums = parse_floats(t);
         if (nums.size() >= 2) {
            out.y.push_back(nums[0]);
            out.x.push_back(nums[1]);
         }
      }
   }
   return out.x.size() > 0 && out.start1.size() == nparam && out.start2.size() == nparam &&
          out.certified.size() == nparam;
}

bool parse_mass_column(std::string const &path, std::string const &column, std::vector<double> &masses)
{
   std::ifstream in(path);
   if (!in) {
      return false;
   }
   std::string header;
   if (!std::getline(in, header)) {
      return false;
   }
   auto names = split_csv(header);
   std::size_t idx = names.size();
   for (std::size_t i = 0; i < names.size(); ++i) {
      if (names[i] == column) {
         idx = i;
         break;
      }
   }
   if (idx >= names.size()) {
      return false;
   }
   std::string line;
   while (std::getline(in, line)) {
      auto cols = split_csv(line);
      if (idx >= cols.size()) {
         continue;
      }
      double m = 0.0;
      if (parse_double(cols[idx], m) && std::isfinite(m)) {
         masses.push_back(m);
      }
   }
   return !masses.empty();
}

bool parse_zmumu_reco_mass(std::string const &path, std::vector<double> &masses)
{
   std::ifstream in(path);
   if (!in) {
      return false;
   }
   std::string header;
   if (!std::getline(in, header)) {
      return false;
   }
   auto names = split_csv(header);
   auto find_idx = [&](std::string const &name) -> int {
      for (std::size_t i = 0; i < names.size(); ++i) {
         if (names[i] == name) {
            return static_cast<int>(i);
         }
      }
      return -1;
   };

   int i_pt1 = find_idx("pt1");
   int i_eta1 = find_idx("eta1");
   int i_phi1 = find_idx("phi1");
   int i_pt2 = find_idx("pt2");
   int i_eta2 = find_idx("eta2");
   int i_phi2 = find_idx("phi2");
   if (i_pt1 < 0 || i_eta1 < 0 || i_phi1 < 0 || i_pt2 < 0 || i_eta2 < 0 || i_phi2 < 0) {
      return false;
   }

   std::string line;
   while (std::getline(in, line)) {
      auto cols = split_csv(line);
      if (cols.size() <= static_cast<std::size_t>(i_phi2)) {
         continue;
      }
      double pt1 = 0.0, eta1 = 0.0, phi1 = 0.0, pt2 = 0.0, eta2 = 0.0, phi2 = 0.0;
      if (!parse_double(cols[i_pt1], pt1) || !parse_double(cols[i_eta1], eta1) || !parse_double(cols[i_phi1], phi1) ||
          !parse_double(cols[i_pt2], pt2) || !parse_double(cols[i_eta2], eta2) ||
          !parse_double(cols[i_phi2], phi2)) {
         continue;
      }
      double m2 = 2.0 * pt1 * pt2 * (std::cosh(eta1 - eta2) - std::cos(phi1 - phi2));
      if (m2 > 0.0 && std::isfinite(m2)) {
         masses.push_back(std::sqrt(m2));
      }
   }

   return !masses.empty();
}

void histogram(std::vector<double> const &masses, double low, double high, std::size_t bins, std::vector<double> &x,
               std::vector<double> &y)
{
   y.assign(bins, 0.0);
   double const w = (high - low) / static_cast<double>(bins);
   for (double m : masses) {
      if (m < low || m >= high) {
         continue;
      }
      std::size_t idx = static_cast<std::size_t>(std::floor((m - low) / w));
      if (idx >= bins) {
         idx = bins - 1;
      }
      y[idx] += 1.0;
   }
   x.resize(bins);
   for (std::size_t i = 0; i < bins; ++i) {
      x[i] = low + (static_cast<double>(i) + 0.5) * w;
   }
}

FunctionMinimum run_migrad(FCNBase const &fcn, MnUserParameters const &upars, unsigned int maxfcn = 300000,
                           double tolerance = 0.01)
{
   MnMigrad migrad(fcn, upars, MnStrategy(2));
   return migrad(maxfcn, tolerance);
}

template <typename SolveFn> bool bench_solve_times(SolveFn &&solve_once, int repeats, int warmups)
{
   for (int i = 0; i < warmups; ++i) {
      if (!solve_once()) {
         return false;
      }
   }
   std::vector<double> times;
   times.reserve(static_cast<std::size_t>(std::max(0, repeats)));
   for (int i = 0; i < repeats; ++i) {
      auto const t0 = std::chrono::steady_clock::now();
      if (!solve_once()) {
         return false;
      }
      auto const t1 = std::chrono::steady_clock::now();
      std::chrono::duration<double> dt = t1 - t0;
      times.push_back(dt.count());
   }

   std::ostringstream oss;
   oss.setf(std::ios::fixed);
   oss << std::setprecision(9);
   for (std::size_t i = 0; i < times.size(); ++i) {
      if (i > 0) {
         oss << ",";
      }
      oss << times[i];
   }
   std::cout << "BENCH_TIMES_S:" << oss.str() << "\n";
   return true;
}

bool run_case_noaa(RunMode mode, int bench_repeats, int bench_warmups)
{
   NoaaFCN fcn;
   if (!parse_noaa("examples/data/noaa/co2_mm_mlo.csv", fcn)) {
      std::cerr << "failed to parse NOAA data\n";
      return false;
   }
   if (mode == RunMode::LoadOnly) {
      return true;
   }

   MnUserParameters u;
   u.Add("a0", fcn.y.front(), 0.5);
   u.Add("a1", 2.0, 0.2);
   u.Add("a2", 0.0, 0.01);
   u.Add("b1", 2.0, 0.2);
   u.Add("c1", 0.0, 0.2);
   u.Add("b2", 0.5, 0.1);
   u.Add("c2", 0.0, 0.1);
   u.Add("d1", 0.0, 0.01);

   auto solve_once = [&]() -> bool {
      auto min = run_migrad(fcn, u);
      MnHesse hesse;
      hesse(fcn, min);
      return min.IsValid();
   };

   if (bench_repeats > 0) {
      return bench_solve_times(solve_once, bench_repeats, bench_warmups);
   }
   return solve_once();
}

std::tuple<std::vector<double>, double>
fit_nist_dataset_with_starts(NistDataset const &ds, LeastSquaresFCN::ModelFn const &model, bool b4_positive,
                             std::vector<std::vector<double>> const &starts)
{
   LeastSquaresFCN fcn;
   fcn.x = ds.x;
   fcn.y = ds.y;
   fcn.model = model;

   auto fit_with_start = [&](std::vector<double> const &start) {
      MnUserParameters u;
      for (std::size_t i = 0; i < start.size(); ++i) {
         std::string name = "b" + std::to_string(i + 1);
         double step = std::max(std::abs(start[i]) * 0.05, 1e-6);
         u.Add(name, start[i], step);
      }
      if (b4_positive && start.size() > 3) {
         u.SetLowerLimit(3, 1e-6);
      }
      auto min = run_migrad(fcn, u, 600000, 0.001);
      return std::make_tuple(min.UserState().Params(), min.Fval());
   };

   std::vector<double> best_params;
   double best_f = std::numeric_limits<double>::infinity();
   for (auto const &start : starts) {
      auto [p, f] = fit_with_start(start);
      if (std::isfinite(f) && f < best_f) {
         best_f = f;
         best_params = std::move(p);
      }
   }
   if (!std::isfinite(best_f)) {
      return {std::vector<double>{}, std::numeric_limits<double>::infinity()};
   }
   return {best_params, best_f};
}

std::tuple<std::vector<double>, double> fit_nist_dataset(NistDataset const &ds, LeastSquaresFCN::ModelFn const &model,
                                                          bool b4_positive)
{
   std::vector<std::vector<double>> starts{ds.start1, ds.start2};
   return fit_nist_dataset_with_starts(ds, model, b4_positive, starts);
}

std::tuple<std::vector<double>, double> fit_hahn_dataset(NistDataset const &ds, LeastSquaresFCN::ModelFn const &model)
{
   std::vector<std::vector<double>> starts{ds.start1, ds.start2, ds.certified};
   std::vector<double> mid(ds.start1.size(), 0.0);
   for (std::size_t i = 0; i < ds.start1.size(); ++i) {
      mid[i] = 0.5 * (ds.start1[i] + ds.start2[i]);
   }
   starts.push_back(mid);
   for (double scale : {0.25, 0.5, 1.5, 2.0}) {
      std::vector<double> s1 = ds.start1;
      std::vector<double> s2 = ds.start2;
      for (double &v : s1) {
         v *= scale;
      }
      for (double &v : s2) {
         v *= scale;
      }
      starts.push_back(std::move(s1));
      starts.push_back(std::move(s2));
   }
   for (std::size_t idx = 0; idx < ds.start2.size(); ++idx) {
      std::vector<double> up = ds.start2;
      std::vector<double> down = ds.start2;
      up[idx] *= 1.35;
      down[idx] *= 0.65;
      starts.push_back(std::move(up));
      starts.push_back(std::move(down));
   }

   return fit_nist_dataset_with_starts(ds, model, false, starts);
}

bool run_case_nist(RunMode mode, int bench_repeats, int bench_warmups)
{
   NistDataset misra;
   NistDataset hahn;
   NistDataset rat;
   if (!parse_nist_dat("examples/data/nist/Misra1a.dat", 2, misra)) {
      std::cerr << "failed to parse Misra1a\n";
      return false;
   }
   if (!parse_nist_dat("examples/data/nist/Hahn1.dat", 7, hahn)) {
      std::cerr << "failed to parse Hahn1\n";
      return false;
   }
   if (!parse_nist_dat("examples/data/nist/Rat43.dat", 4, rat)) {
      std::cerr << "failed to parse Rat43\n";
      return false;
   }
   if (mode == RunMode::LoadOnly) {
      return true;
   }

   auto misra_model = [](std::vector<double> const &p, double x) {
      return p[0] * (1.0 - std::exp(-p[1] * x));
   };
   auto hahn_model = [](std::vector<double> const &p, double x) {
      double x2 = x * x;
      double x3 = x2 * x;
      double num = p[0] + p[1] * x + p[2] * x2 + p[3] * x3;
      double den = 1.0 + p[4] * x + p[5] * x2 + p[6] * x3;
      if (std::abs(den) < 1e-14) {
         return std::numeric_limits<double>::quiet_NaN();
      }
      return num / den;
   };
   auto rat_model = [](std::vector<double> const &p, double x) {
      if (p[3] <= 0.0) {
         return std::numeric_limits<double>::quiet_NaN();
      }
      double expo = std::clamp(p[1] - p[2] * x, -700.0, 700.0);
      double base = 1.0 + std::exp(expo);
      return p[0] / std::pow(base, 1.0 / p[3]);
   };

   auto solve_once = [&]() -> bool {
      auto [_, f_m] = fit_nist_dataset(misra, misra_model, false);
      auto [__, f_h] = fit_hahn_dataset(hahn, hahn_model);
      auto [___, f_r] = fit_nist_dataset(rat, rat_model, true);
      return std::isfinite(f_m) && std::isfinite(f_h) && std::isfinite(f_r);
   };
   if (bench_repeats > 0) {
      return bench_solve_times(solve_once, bench_repeats, bench_warmups);
   }
   return solve_once();
}

bool fit_hist(std::vector<double> const &x, std::vector<double> const &y, double mu0, double sig0)
{
   HistFCN fcn;
   fcn.x = x;
   fcn.y = y;
   fcn.sigma.resize(y.size());
   for (std::size_t i = 0; i < y.size(); ++i) {
      fcn.sigma[i] = std::sqrt(std::max(y[i], 1.0));
   }

   double max_count = *std::max_element(y.begin(), y.end());
   double mean_bg = 0.0;
   for (double v : y) {
      mean_bg += v;
   }
   mean_bg /= static_cast<double>(std::max<std::size_t>(1, y.size()));

   MnUserParameters u;
   u.Add("amp", max_count, std::max(1.0, max_count * 0.1));
   u.Add("mu", mu0, std::max(sig0 * 0.2, 0.01));
   u.Add("sigma", sig0, std::max(sig0 * 0.1, 0.01));
   u.Add("c0", std::max(1.0, mean_bg), 0.5);
   u.Add("c1", 0.0, 0.05);
   u.SetLowerLimit(0, 0.0);
   u.SetLowerLimit(2, 0.05);
   u.SetLowerLimit(3, 0.0);

   auto min = run_migrad(fcn, u);
   return min.IsValid();
}

bool run_case_cern(RunMode mode, int bench_repeats, int bench_warmups)
{
   std::vector<double> murun_masses;
   std::vector<double> zmumu_masses;
   if (!parse_mass_column("examples/data/cern/MuRun2010B_0.csv", "M", murun_masses)) {
      std::cerr << "failed to parse MuRun masses\n";
      return false;
   }
   if (!parse_zmumu_reco_mass("examples/data/cern/Zmumu.csv", zmumu_masses)) {
      std::cerr << "failed to parse Zmumu masses\n";
      return false;
   }

   std::vector<double> murun_jpsi;
   for (double m : murun_masses) {
      if (m >= 2.0 && m <= 5.0) {
         murun_jpsi.push_back(m);
      }
   }
   std::vector<double> zmumu_z;
   for (double m : zmumu_masses) {
      if (m >= 60.0 && m <= 120.0) {
         zmumu_z.push_back(m);
      }
   }

   std::vector<double> x1, y1, x2, y2;
   histogram(murun_jpsi, 2.0, 5.0, 60, x1, y1);
   histogram(zmumu_z, 60.0, 120.0, 60, x2, y2);
   if (mode == RunMode::LoadOnly) {
      return true;
   }

   auto solve_once = [&]() -> bool {
      bool ok1 = fit_hist(x1, y1, 3.10, 0.12);
      bool ok2 = fit_hist(x2, y2, 91.0, 2.5);
      return ok1 && ok2;
   };
   if (bench_repeats > 0) {
      return bench_solve_times(solve_once, bench_repeats, bench_warmups);
   }
   return solve_once();
}

bool run_case_usgs(RunMode mode, int bench_repeats, int bench_warmups)
{
   std::vector<double> mags;
   if (!parse_usgs_magnitudes("examples/data/usgs/earthquakes_2025_m4p5.csv", mags)) {
      std::cerr << "failed to parse USGS data\n";
      return false;
   }
   double mmax = std::floor(*std::max_element(mags.begin(), mags.end()));

   std::vector<double> mvals;
   std::vector<double> counts;
   build_cumulative(mags, 4.5, mmax, 0.1, mvals, counts);

   std::vector<double> logn;
   std::vector<double> sigma;
   logn.reserve(counts.size());
   sigma.reserve(counts.size());
   for (double c : counts) {
      logn.push_back(std::log10(c));
      sigma.push_back(1.0 / (std::log(10.0) * std::sqrt(c)));
   }
   if (mode == RunMode::LoadOnly) {
      return true;
   }

   struct USGSFCN final : public FCNBase {
      std::vector<double> m, logn, sigma;
      double operator()(std::vector<double> const &p) const override
      {
         double chi2 = 0.0;
         for (std::size_t i = 0; i < m.size(); ++i) {
            double pred = p[0] - p[1] * m[i];
            if (!std::isfinite(pred)) {
               return 1e30;
            }
            double r = (logn[i] - pred) / sigma[i];
            chi2 += r * r;
         }
         return chi2;
      }
      double Up() const override { return 1.0; }
   } fcn;
   fcn.m = std::move(mvals);
   fcn.logn = std::move(logn);
   fcn.sigma = std::move(sigma);

   MnUserParameters u;
   u.Add("a", 5.0, 0.1);
   u.Add("b", 1.0, 0.05);

   auto solve_once = [&]() -> bool {
      auto min = run_migrad(fcn, u);
      MnHesse hesse;
      hesse(fcn, min);
      return min.IsValid();
   };
   if (bench_repeats > 0) {
      return bench_solve_times(solve_once, bench_repeats, bench_warmups);
   }
   return solve_once();
}

} // namespace

int main(int argc, char **argv)
{
   std::string case_id;
   RunMode mode = RunMode::Full;
   int bench_repeats = 0;
   int bench_warmups = 0;
   for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--case" && i + 1 < argc) {
         case_id = argv[++i];
      } else if (arg.rfind("--case=", 0) == 0) {
         case_id = arg.substr(7);
      } else if (arg == "--bench-repeats" && i + 1 < argc) {
         try {
            bench_repeats = std::stoi(argv[++i]);
         } catch (...) {
            std::cerr << "invalid --bench-repeats value\n";
            return 2;
         }
         if (bench_repeats < 0) {
            std::cerr << "--bench-repeats must be non-negative\n";
            return 2;
         }
      } else if (arg.rfind("--bench-repeats=", 0) == 0) {
         try {
            bench_repeats = std::stoi(arg.substr(16));
         } catch (...) {
            std::cerr << "invalid --bench-repeats value\n";
            return 2;
         }
         if (bench_repeats < 0) {
            std::cerr << "--bench-repeats must be non-negative\n";
            return 2;
         }
      } else if (arg == "--bench-warmups" && i + 1 < argc) {
         try {
            bench_warmups = std::stoi(argv[++i]);
         } catch (...) {
            std::cerr << "invalid --bench-warmups value\n";
            return 2;
         }
         if (bench_warmups < 0) {
            std::cerr << "--bench-warmups must be non-negative\n";
            return 2;
         }
      } else if (arg.rfind("--bench-warmups=", 0) == 0) {
         try {
            bench_warmups = std::stoi(arg.substr(16));
         } catch (...) {
            std::cerr << "invalid --bench-warmups value\n";
            return 2;
         }
         if (bench_warmups < 0) {
            std::cerr << "--bench-warmups must be non-negative\n";
            return 2;
         }
      } else if (arg == "--mode" && i + 1 < argc) {
         std::string const value = argv[++i];
         if (value == "full") {
            mode = RunMode::Full;
         } else if (value == "load-only") {
            mode = RunMode::LoadOnly;
         } else if (value == "solve-only") {
            mode = RunMode::SolveOnly;
         } else {
            std::cerr << "invalid mode: " << value << " (expected full|load-only|solve-only)\n";
            return 2;
         }
      } else if (arg.rfind("--mode=", 0) == 0) {
         std::string const value = arg.substr(7);
         if (value == "full") {
            mode = RunMode::Full;
         } else if (value == "load-only") {
            mode = RunMode::LoadOnly;
         } else if (value == "solve-only") {
            mode = RunMode::SolveOnly;
         } else {
            std::cerr << "invalid mode: " << value << " (expected full|load-only|solve-only)\n";
            return 2;
         }
      }
   }

   if (case_id.empty()) {
      std::cerr << "usage: scientific_runner --case <noaa_co2|nist_strd|usgs_earthquakes|cern_dimuon> "
                   "[--mode full|load-only|solve-only] [--bench-repeats N --bench-warmups W]\n";
      return 2;
   }

   bool ok = false;
   if (case_id == "noaa_co2") {
      ok = run_case_noaa(mode, bench_repeats, bench_warmups);
   } else if (case_id == "nist_strd") {
      ok = run_case_nist(mode, bench_repeats, bench_warmups);
   } else if (case_id == "usgs_earthquakes") {
      ok = run_case_usgs(mode, bench_repeats, bench_warmups);
   } else if (case_id == "cern_dimuon") {
      ok = run_case_cern(mode, bench_repeats, bench_warmups);
   } else {
      std::cerr << "unknown case: " << case_id << "\n";
      return 3;
   }

   if (!ok) {
      std::cerr << "case failed: " << case_id << "\n";
      return 1;
   }
   std::cout << "ok " << case_id << "\n";
   return 0;
}
