// ROOT Minuit2 tag v6-36-08: user-side call accounting evidence.
// From the repository root, with an existing ROOT reference build
// (scripts/build_root_reference_runner.sh v6-36-08 creates the archives).
// Each run prints one JSON document, which is the stored file:
//   mkdir -p target/nfcn-accounting-4
//   g++ -O2 -g -rdynamic -std=c++23 -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/nfcn_accounting.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o target/nfcn-accounting-4/nfcn_accounting_root
//   target/nfcn-accounting-4/nfcn_accounting_root --trace quadratic2_no_g2_migrad \
//     > reports/verification/raw/trace/nfcn_accounting_quadratic2_no_g2_migrad.root.json
//   target/nfcn-accounting-4/nfcn_accounting_root --trace rosenbrock2_maxfcn140 \
//     > reports/verification/raw/trace/nfcn_accounting_rosenbrock2_maxfcn140.root.json
//   target/nfcn-accounting-4/nfcn_accounting_root --trace quadratic2_minos_before_minos \
//     > reports/verification/raw/trace/nfcn_accounting_quadratic2_minos_before_minos.root.json
//   target/nfcn-accounting-4/nfcn_accounting_root --trace quadratic2_upper_limited_migrad \
//     > reports/verification/raw/trace/nfcn_accounting_quadratic2_upper_limited_migrad.root.json
// Built in the C++23 mode of the house standard; this file uses no C++23 feature.
// The stack parsing targets glibc on Linux (see phase_counts.h).
// Reuse the runner's objectives without changing its stored reference results.
// The document holds "provenance" (the build and run commands above) and "result".
// The program exits with 1 if reported_nfcn and user_calls differ.
// Pass --trace to include every evaluation and its phase. Set M2_TRACE_STACKS to
// print call stacks on stderr for phase identification.
#include "phase_counts.h"

#define main reference_runner_main
#include "../ref_runner_cpp/main.cpp"
#undef main

// Commands that build this program; the run command is added to the provenance.
static std::string const program = "tools/root_evidence/nfcn_accounting.cpp";
static std::string const buildCommand =
   "mkdir -p target/nfcn-accounting-4 && g++ -O2 -g -rdynamic -std=c++23 -Ithird_party/root_ref/math/minuit2/inc "
   "tools/root_evidence/nfcn_accounting.cpp third_party/root_ref_build/minuit2/lib/libMinuit2.a "
   "third_party/root_ref_build/minuit2/lib/libMinuit2Math.a -lpthread "
   "-o target/nfcn-accounting-4/nfcn_accounting_root";

template <class Objective> class Counted final : public FCNBase {
public:
   mutable PhaseCounts counts;
   mutable unsigned int hessianCalls = 0;

   double operator()(std::vector<double> const &p) const override
   {
      counts.record(p);
      return objective(p);
   }
   double Up() const override { return objective.Up(); }
   bool HasGradient() const override { return objective.HasGradient(); }
   std::vector<double> Gradient(std::vector<double> const &p) const override { return objective.Gradient(p); }
   bool HasHessian() const override { return objective.HasHessian(); }
   std::vector<double> Hessian(std::vector<double> const &p) const override
   {
      ++hessianCalls;
      return objective.Hessian(p);
   }
   bool HasG2() const override { return objective.HasG2(); }
   std::vector<double> G2(std::vector<double> const &p) const override { return objective.G2(p); }

private:
   Objective objective;
};

int main(int argc, char **argv)
{
   char const *usage = "usage: nfcn_accounting_root [--trace] "
                       "quadratic2_no_g2_migrad|rosenbrock2_maxfcn140|quadratic2_minos_before_minos|quadratic2_upper_limited_migrad";
   bool trace = false;
   bool consistent = true;
   std::string workload;
   std::string command = "target/nfcn-accounting-4/nfcn_accounting_root";
   for (int i = 1; i < argc; ++i) {
      std::string const arg = argv[i];
      command += " " + arg;
      if (arg == "--trace") {
         trace = true;
      } else if (arg.rfind("--", 0) == 0 || !workload.empty()) {
         std::cerr << "rejected argument " << arg << '\n' << usage << '\n';
         return 1;
      } else {
         workload = arg;
      }
   }
   if (workload == "quadratic2_upper_limited_migrad") {
      Counted<Quadratic2> quadratic;
      quadratic.counts.trace = trace;
      MnUserParameters parameters;
      parameters.Add("x", 0.4, 0.1);
      parameters.Add("y", -1.0, 0.1);
      parameters.SetUpperLimit(0, 1.8);
      parameters.SetUpperLimit(1, -1.5);
      auto const minimum = MnMigrad(quadratic, parameters)(0, 0.1);
      beginRecord(provenanceJson(program, buildCommand, command), workload);
      field("reported_nfcn", minimum.NFcn());
      field("valid", minimum.IsValid());
      field("edm", minimum.Edm());
      quadratic.counts.print();
      endRecord();
      consistent = quadratic.counts.reportedMatches(minimum.NFcn());
   } else if (workload == "quadratic2_no_g2_migrad") {
      Counted<QuadraticNoG2> quadratic;
      quadratic.counts.trace = trace;
      MnUserParameters parameters;
      parameters.Add("x", 0.4, 0.1);
      parameters.Add("y", -1.0, 0.1);
      auto const minimum = MnMigrad(quadratic, parameters)();
      beginRecord(provenanceJson(program, buildCommand, command), workload);
      field("reported_nfcn", minimum.NFcn());
      field("hessian_calls", quadratic.hessianCalls);
      std::cout << ",\n    \"state_nfcn\": [";
      for (std::size_t i = 0; i < minimum.States().size(); ++i)
         std::cout << (i ? "," : "") << "\n      " << minimum.States()[i].NFcn();
      std::cout << "\n    ]";
      quadratic.counts.print();
      endRecord();
      consistent = quadratic.counts.reportedMatches(minimum.NFcn());
   } else if (workload == "rosenbrock2_maxfcn140") {
      Counted<Rosenbrock2> rosenbrock;
      rosenbrock.counts.trace = trace;
      MnUserParameters start;
      start.Add("x", 0.0, 0.1);
      start.Add("y", 0.0, 0.1);
      auto const continuation = MnMigrad(rosenbrock, start)(140, 0.1);
      beginRecord(provenanceJson(program, buildCommand, command), workload);
      field("reported_nfcn", continuation.NFcn());
      field("valid", continuation.IsValid());
      field("call_limit", continuation.HasReachedCallLimit());
      field("above_edm", continuation.IsAboveMaxEdm());
      field("edm", continuation.Edm());
      rosenbrock.counts.print();
      endRecord();
      consistent = rosenbrock.counts.reportedMatches(continuation.NFcn());
   } else if (workload == "quadratic2_minos_before_minos") {
      // Match the Minos runner up to the point whose NFcn it stores.
      Counted<Quadratic2> minos;
      minos.counts.trace = trace;
      MnUserParameters minosStart;
      minosStart.Add("x", 0.4, 0.1);
      minosStart.Add("y", -1.0, 0.1);
      auto minosMinimum = MnMigrad(minos, minosStart)();
      auto migradCalls = minosMinimum.NFcn();
      minos.counts.beginStandalone();
      MnHesse{}(minos, minosMinimum);
      beginRecord(provenanceJson(program, buildCommand, command), workload);
      field("migrad_nfcn", migradCalls);
      field("reported_nfcn", minosMinimum.NFcn());
      field("hessian_calls", minos.hessianCalls);
      minos.counts.print();
      endRecord();
      consistent = minos.counts.reportedMatches(minosMinimum.NFcn());
   } else {
      std::cerr << usage << '\n';
      return 1;
   }
   return consistent && std::cout ? 0 : 1;
}
