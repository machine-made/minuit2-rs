#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnContours.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnTraceObject.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/ModularFunctionMinimizer.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserCovariance.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using ROOT::Minuit2::FCNBase;
using ROOT::Minuit2::FunctionMinimum;
using ROOT::Minuit2::MinosError;
using ROOT::Minuit2::MnContours;
using ROOT::Minuit2::MnHesse;
using ROOT::Minuit2::MnMinimize;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnMinos;
using ROOT::Minuit2::MnScan;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnStrategy;
using ROOT::Minuit2::MnUserCovariance;
using ROOT::Minuit2::MnUserParameterState;
using ROOT::Minuit2::MnTraceObject;
using ROOT::Minuit2::MnUserParameters;

namespace {

class JsonlTraceObject final : public MnTraceObject {
public:
   explicit JsonlTraceObject(char const *path) : fOut(path, std::ios::out | std::ios::trunc) {}

   void operator()(int iter, ROOT::Minuit2::MinimumState const &state) override
   {
      double grad_norm = 0.0;
      auto const &grad = state.Gradient().Grad();
      for (std::size_t i = 0; i < grad.size(); ++i) {
         double const g = grad(i);
         grad_norm += g * g;
      }
      grad_norm = std::sqrt(grad_norm);

      fOut << "{\"runner\":\"root-minuit2\",\"iter\":" << iter << ",\"nfcn\":" << state.NFcn()
           << ",\"fval\":" << std::setprecision(17) << state.Fval() << ",\"edm\":" << state.Edm()
           << ",\"lambda\":null,\"grad_norm\":" << grad_norm << ",\"dcovar\":" << state.Error().Dcovar()
           << "}\n";
   }

private:
   std::ofstream fOut;
};

JsonlTraceObject *maybe_attach_trace(MnMigrad &migrad)
{
   char const *path = std::getenv("MINUIT2_ROOT_TRACE_JSONL");
   if (path == nullptr || path[0] == '\0') {
      return nullptr;
   }
   auto *trace = new JsonlTraceObject(path);
   migrad.Minimizer().Builder().SetTraceObject(*trace);
   return trace;
}

class Quadratic3 final : public FCNBase {
public:
   double operator()(std::vector<double> const &p) const override
   {
      const double x = p.at(0);
      const double y = p.at(1);
      const double z = p.at(2);
      return x * x + 10.0 * y * y + 100.0 * z * z + 2.0 * x * y + 4.0 * x * z + 8.0 * y * z;
   }

   bool HasGradient() const override { return true; }

   std::vector<double> Gradient(std::vector<double> const &p) const override
   {
      const double x = p.at(0);
      const double y = p.at(1);
      const double z = p.at(2);
      return {
         2.0 * x + 2.0 * y + 4.0 * z,
         2.0 * x + 20.0 * y + 8.0 * z,
         4.0 * x + 8.0 * y + 200.0 * z,
      };
   }

   double Up() const override { return 1.0; }
};

class Rosenbrock2 final : public FCNBase {
public:
   double operator()(std::vector<double> const &p) const override
   {
      const double x = p.at(0);
      const double y = p.at(1);
      const double t1 = y - x * x;
      const double t2 = 1.0 - x;
      return 100.0 * t1 * t1 + t2 * t2;
   }

   double Up() const override { return 1.0; }
};

class Quadratic2 final : public FCNBase {
public:
   double operator()(std::vector<double> const &p) const override
   {
      const double x = p.at(0);
      const double y = p.at(1);
      const double dx = x - 1.0;
      const double dy = y + 2.0;
      return dx * dx + 4.0 * dy * dy + 0.3 * x * y;
   }

   double Up() const override { return 1.0; }
};

class QuadraticNoG2 final : public FCNBase {
public:
   double operator()(std::vector<double> const &p) const override
   {
      const double x = p.at(0);
      const double y = p.at(1);
      const double dx = x - 1.0;
      const double dy = y + 2.0;
      return dx * dx + dy * dy;
   }

   bool HasGradient() const override { return true; }

   std::vector<double> Gradient(std::vector<double> const &p) const override
   {
      return {2.0 * (p.at(0) - 1.0), 2.0 * (p.at(1) + 2.0)};
   }

   bool HasHessian() const override { return true; }

   std::vector<double> Hessian(std::vector<double> const &) const override
   {
      // Packed upper triangle of 2x2 identity-scaled Hessian.
      return {2.0, 0.0, 2.0};
   }

   bool HasG2() const override { return false; }

   std::vector<double> G2(std::vector<double> const &) const override { return {}; }

   double Up() const override { return 1.0; }
};

struct RunResult {
   std::string workload;
   std::string algorithm;
   bool valid = false;
   double fval = 0.0;
   double edm = 0.0;
   unsigned int nfcn = 0;
   std::vector<double> params;
   std::vector<double> errors;
   bool has_covariance = false;
   std::vector<std::vector<double>> covariance;
   bool has_minos = false;
   bool minos_valid = false;
   unsigned int minos_parameter = 0;
   double minos_lower = 0.0;
   double minos_upper = 0.0;
};

std::string to_json_array(std::vector<double> const &v)
{
   std::ostringstream oss;
   oss << "[";
   for (std::size_t i = 0; i < v.size(); ++i) {
      if (i != 0) {
         oss << ",";
      }
      oss << std::setprecision(17) << v[i];
   }
   oss << "]";
   return oss.str();
}

std::string to_json_matrix(std::vector<std::vector<double>> const &m)
{
   std::ostringstream oss;
   oss << "[";
   for (std::size_t i = 0; i < m.size(); ++i) {
      if (i != 0) {
         oss << ",";
      }
      oss << to_json_array(m[i]);
   }
   oss << "]";
   return oss.str();
}

std::string to_json(RunResult const &r)
{
   std::ostringstream oss;
   oss << "{";
   oss << "\"runner\":\"root-minuit2\",";
   oss << "\"workload\":\"" << r.workload << "\",";
   oss << "\"algorithm\":\"" << r.algorithm << "\",";
   oss << "\"valid\":" << (r.valid ? "true" : "false") << ",";
   oss << "\"fval\":" << std::setprecision(17) << r.fval << ",";
   oss << "\"edm\":" << std::setprecision(17) << r.edm << ",";
   oss << "\"nfcn\":" << r.nfcn << ",";
   oss << "\"params\":" << to_json_array(r.params) << ",";
   oss << "\"errors\":" << to_json_array(r.errors) << ",";
   oss << "\"has_covariance\":" << (r.has_covariance ? "true" : "false") << ",";
   if (r.has_covariance) {
      oss << "\"covariance\":" << to_json_matrix(r.covariance) << ",";
   } else {
      oss << "\"covariance\":null,";
   }
   oss << "\"has_minos\":" << (r.has_minos ? "true" : "false") << ",";
   if (r.has_minos) {
      oss << "\"minos\":{";
      oss << "\"valid\":" << (r.minos_valid ? "true" : "false") << ",";
      oss << "\"parameter\":" << r.minos_parameter << ",";
      oss << "\"lower\":" << std::setprecision(17) << r.minos_lower << ",";
      oss << "\"upper\":" << std::setprecision(17) << r.minos_upper;
      oss << "}";
   } else {
      oss << "\"minos\":null";
   }
   oss << "}";
   return oss.str();
}

std::vector<std::vector<double>> covariance_to_dense(MnUserCovariance const &cov)
{
   const unsigned int n = cov.Nrow();
   std::vector<std::vector<double>> out(n, std::vector<double>(n, 0.0));
   for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
         out[i][j] = cov(i, j);
      }
   }
   return out;
}

void fill_common(RunResult &result, FunctionMinimum const &minimum)
{
   MnUserParameterState const &state = minimum.UserState();
   result.valid = minimum.IsValid();
   result.fval = minimum.Fval();
   result.edm = minimum.Edm();
   result.nfcn = minimum.NFcn();
   result.params = state.Params();
   result.errors = state.Errors();

   if (state.HasCovariance()) {
      result.has_covariance = true;
      result.covariance = covariance_to_dense(state.Covariance());
   }
}

RunResult run_quadratic3_fixx_migrad()
{
   Quadratic3 fcn;
   MnUserParameters upar;
   upar.Add("x", 1.0, 0.1);
   upar.Add("y", 2.0, 0.1);
   upar.Add("z", 3.0, 0.1);
   upar.Fix(0);

   MnMigrad migrad(fcn, upar);
   JsonlTraceObject *trace = maybe_attach_trace(migrad);
   FunctionMinimum minimum = migrad();
   delete trace;

   RunResult result;
   result.workload = "quadratic3_fixx_migrad";
   result.algorithm = "migrad";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic3_fixx_hesse()
{
   Quadratic3 fcn;
   MnUserParameters upar;
   upar.Add("x", 1.0, 0.1);
   upar.Add("y", 2.0, 0.1);
   upar.Add("z", 3.0, 0.1);
   upar.Fix(0);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnHesse hesse;
   hesse(fcn, minimum);

   RunResult result;
   result.workload = "quadratic3_fixx_hesse";
   result.algorithm = "migrad+hesse";
   fill_common(result, minimum);
   return result;
}

RunResult run_rosenbrock2_migrad()
{
   Rosenbrock2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.0, 0.1);
   upar.Add("y", 0.0, 0.1);

   MnMigrad migrad(fcn, upar);
   JsonlTraceObject *trace = maybe_attach_trace(migrad);
   FunctionMinimum minimum = migrad();
   delete trace;

   RunResult result;
   result.workload = "rosenbrock2_migrad";
   result.algorithm = "migrad";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_minos_p0()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnHesse hesse;
   hesse(fcn, minimum);

   MnMinos minos(fcn, minimum);
   MinosError me = minos.Minos(0);

   RunResult result;
   result.workload = "quadratic2_minos_p0";
   result.algorithm = "migrad+hesse+minos";
   fill_common(result, minimum);
   result.has_minos = true;
   result.minos_valid = me.IsValid();
   result.minos_parameter = me.Parameter();
   result.minos_lower = me.Lower();
   result.minos_upper = me.Upper();
   return result;
}

RunResult run_quadratic2_minos_p1()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnHesse hesse;
   hesse(fcn, minimum);

   MnMinos minos(fcn, minimum);
   MinosError me = minos.Minos(1);

   RunResult result;
   result.workload = "quadratic2_minos_p1";
   result.algorithm = "migrad+hesse+minos";
   fill_common(result, minimum);
   result.has_minos = true;
   result.minos_valid = me.IsValid();
   result.minos_parameter = me.Parameter();
   result.minos_lower = me.Lower();
   result.minos_upper = me.Upper();
   return result;
}

RunResult run_quadratic2_simplex()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnSimplex simplex(fcn, upar);
   FunctionMinimum minimum = simplex();

   RunResult result;
   result.workload = "quadratic2_simplex";
   result.algorithm = "simplex";
   fill_common(result, minimum);
   return result;
}

RunResult run_rosenbrock2_minimize()
{
   Rosenbrock2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.0, 0.1);
   upar.Add("y", 0.0, 0.1);

   MnMinimize minimize(fcn, upar);
   FunctionMinimum minimum = minimize();

   RunResult result;
   result.workload = "rosenbrock2_minimize";
   result.algorithm = "minimize";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_limited_migrad()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1, 0.0, 2.0);
   upar.Add("y", -1.0, 0.1, -3.0, -1.0);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();

   RunResult result;
   result.workload = "quadratic2_limited_migrad";
   result.algorithm = "migrad";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_lower_limited_migrad()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);
   upar.SetLowerLimit(0, 0.0);
   upar.SetLowerLimit(1, -2.5);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();

   RunResult result;
   result.workload = "quadratic2_lower_limited_migrad";
   result.algorithm = "migrad";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_upper_limited_migrad()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);
   upar.SetUpperLimit(0, 1.8);
   upar.SetUpperLimit(1, -1.5);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();

   RunResult result;
   result.workload = "quadratic2_upper_limited_migrad";
   result.algorithm = "migrad";
   fill_common(result, minimum);
   return result;
}

RunResult run_rosenbrock2_migrad_strategy2()
{
   Rosenbrock2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.0, 0.1);
   upar.Add("y", 0.0, 0.1);

   MnMigrad migrad(fcn, upar, MnStrategy(2));
   FunctionMinimum minimum = migrad();

   RunResult result;
   result.workload = "rosenbrock2_migrad_strategy2";
   result.algorithm = "migrad_s2";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_scan_p0()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnScan scan(fcn, minimum.UserState());
   auto points = scan.Scan(0, 61, 0.0, 0.0);
   (void)points;

   RunResult result;
   result.workload = "quadratic2_scan_p0";
   result.algorithm = "migrad+scan";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_scan_p1_limited()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1, 0.0, 2.0);
   upar.Add("y", -1.0, 0.1, -3.0, -1.0);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnScan scan(fcn, minimum.UserState());
   auto points = scan.Scan(1, 61, 0.0, 0.0);
   (void)points;

   RunResult result;
   result.workload = "quadratic2_scan_p1_limited";
   result.algorithm = "migrad+scan";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_contours_01()
{
   Quadratic2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();
   MnHesse hesse;
   hesse(fcn, minimum);
   MnContours contours(fcn, minimum);
   auto points = contours(0, 1, 12);
   (void)points;

   RunResult result;
   result.workload = "quadratic2_contours_01";
   result.algorithm = "migrad+hesse+contours";
   fill_common(result, minimum);
   return result;
}

RunResult run_quadratic2_no_g2_migrad()
{
   QuadraticNoG2 fcn;
   MnUserParameters upar;
   upar.Add("x", 0.4, 0.1);
   upar.Add("y", -1.0, 0.1);

   MnMigrad migrad(fcn, upar);
   FunctionMinimum minimum = migrad();

   RunResult result;
   result.workload = "quadratic2_no_g2_migrad";
   result.algorithm = "migrad_no_g2";
   fill_common(result, minimum);
   return result;
}

std::string parse_workload_arg(int argc, char **argv)
{
   for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--workload" && i + 1 < argc) {
         return argv[i + 1];
      }
      if (arg.rfind("--workload=", 0) == 0) {
         return arg.substr(std::string("--workload=").size());
      }
   }
   return {};
}

} // namespace

int main(int argc, char **argv)
{
   const std::string workload = parse_workload_arg(argc, argv);
   if (workload.empty()) {
      std::cerr << "usage: ref_runner --workload <id>\n";
      return 2;
   }

   RunResult result;
   if (workload == "quadratic3_fixx_migrad") {
      result = run_quadratic3_fixx_migrad();
   } else if (workload == "quadratic3_fixx_hesse") {
      result = run_quadratic3_fixx_hesse();
   } else if (workload == "rosenbrock2_migrad") {
      result = run_rosenbrock2_migrad();
   } else if (workload == "quadratic2_minos_p0") {
      result = run_quadratic2_minos_p0();
   } else if (workload == "quadratic2_minos_p1") {
      result = run_quadratic2_minos_p1();
   } else if (workload == "quadratic2_simplex") {
      result = run_quadratic2_simplex();
   } else if (workload == "rosenbrock2_minimize") {
      result = run_rosenbrock2_minimize();
   } else if (workload == "quadratic2_limited_migrad") {
      result = run_quadratic2_limited_migrad();
   } else if (workload == "quadratic2_lower_limited_migrad") {
      result = run_quadratic2_lower_limited_migrad();
   } else if (workload == "quadratic2_upper_limited_migrad") {
      result = run_quadratic2_upper_limited_migrad();
   } else if (workload == "rosenbrock2_migrad_strategy2") {
      result = run_rosenbrock2_migrad_strategy2();
   } else if (workload == "quadratic2_scan_p0") {
      result = run_quadratic2_scan_p0();
   } else if (workload == "quadratic2_scan_p1_limited") {
      result = run_quadratic2_scan_p1_limited();
   } else if (workload == "quadratic2_contours_01") {
      result = run_quadratic2_contours_01();
   } else if (workload == "quadratic2_no_g2_migrad") {
      result = run_quadratic2_no_g2_migrad();
   } else {
      std::cerr << "unknown workload: " << workload << "\n";
      return 3;
   }

   std::cout << to_json(result) << "\n";
   return 0;
}
