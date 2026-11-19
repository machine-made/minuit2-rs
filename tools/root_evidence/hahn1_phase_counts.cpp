// ROOT Minuit2 v6-36-08: phase counts for Hahn1, NIST start 2.
// From the repository root, build with the existing reference libraries:
//   g++ -O2 -g -rdynamic -std=c++23 -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/hahn1_phase_counts.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o target/nfcn-accounting-4/hahn1_phase_counts_root
//   target/nfcn-accounting-4/hahn1_phase_counts_root --trace
// The dataset argument defaults to tools/hahn1_probe/hahn1_dataset.txt.
// Phase classification is in phase_counts.h; M2_TRACE_STACKS prints stacks.
#include "phase_counts.h"
#include <fstream>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <cmath>
using namespace ROOT::Minuit2;

// Same objective and arithmetic as tools/hahn1_probe/hahn1_probe.cpp.
struct CountedHahn1 : ROOT::Minuit2::FCNBase {
  mutable PhaseCounts counts;
  std::vector<double> x, y;
  double operator()(const std::vector<double>& p) const override {
    counts.record(p);
    double rss = 0;
    for (size_t i = 0; i < x.size(); i++) {
      double X = x[i], x2 = X*X, x3 = x2*X;
      double den = 1.0 + p[4]*X + p[5]*x2 + p[6]*x3;
      if (std::fabs(den) < 1e-300) return 1e30;
      double pred = (p[0] + p[1]*X + p[2]*x2 + p[3]*x3)/den;
      if (!std::isfinite(pred)) return 1e30;
      double r = y[i]-pred; rss += r*r;
    }
    return rss;
  }
  double Up() const override { return 1.0; }
};
int main(int argc, char **argv)
{
   CountedHahn1 fcn;
   std::string path = "tools/hahn1_probe/hahn1_dataset.txt";
   for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "--trace") fcn.counts.trace = true;
      else path = argv[i];
   }
   std::ifstream input(path);
   int n = 0;
   if (!(input >> n) || n <= 0) throw std::runtime_error("invalid Hahn1 dataset");
   fcn.x.resize(n); fcn.y.resize(n);
   for (int i = 0; i < n; ++i)
      if (!(input >> fcn.x[i] >> fcn.y[i])) throw std::runtime_error("invalid Hahn1 observation");
   MnUserParameters parameters;
   for (int i = 0; i < 7; ++i) {
      double start;
      if (!(input >> start)) throw std::runtime_error("invalid Hahn1 start");
      parameters.Add("b" + std::to_string(i + 1), start, 0.1);
   }
   auto minimum = MnMigrad(fcn, parameters, MnStrategy(1))(1000000, 0.1);
   std::cout << std::boolalpha << std::setprecision(17)
             << "{\"workload\":\"hahn1_start2\",\"reported_nfcn\":" << minimum.NFcn()
             << ",\"valid\":" << minimum.IsValid() << ",\"fval\":" << minimum.Fval();
   fcn.counts.print();
   std::cout << "}\n";
   return std::cout ? 0 : 1;
}
