// ROOT Minuit2 tag v6-36-08: a second explicit Hesse pass after strategy-2 Migrad.
// From the repository root, with an existing ROOT reference build
// (scripts/build_root_reference_runner.sh v6-36-08 creates the archives).
// The run prints one JSON document, which is the stored file:
//   mkdir -p target/nfcn-accounting-4
//   g++ -O2 -g -rdynamic -std=c++23 -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/repeated_hesse.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o target/nfcn-accounting-4/repeated_hesse_root
//   target/nfcn-accounting-4/repeated_hesse_root \
//     > reports/verification/raw/trace/repeated_hesse_strategy2.root.json
// Built in the C++23 mode of the house standard; this file uses no C++23 feature.
// Reuse the runner's Rosenbrock objective without changing its stored reference results.
// The document holds "provenance" (the build and run commands above) and "result".
// The workload is strategy-2 Migrad on Rosenbrock from (-1.2, 1.0) with steps 0.2,
// followed by one explicit strategy-2 Hesse pass on the returned minimum.
// The matrices are Error().InvHessian(), the error matrix Minuit2 stores and the
// Migrad state returns. ROOT's user covariance is 2 * Up() times this matrix.
// ROOT stores the refined gradient steps of the Migrad-time Hesse pass.
// The explicit pass starts from them, so the second matrix differs from the first.
#include "phase_counts.h"

#define main reference_runner_main
#include "../ref_runner_cpp/main.cpp"
#undef main

// Commands that build this program; the run command is added to the provenance.
static std::string const program = "tools/root_evidence/repeated_hesse.cpp";
static std::string const buildCommand =
   "mkdir -p target/nfcn-accounting-4 && g++ -O2 -g -rdynamic -std=c++23 -Ithird_party/root_ref/math/minuit2/inc "
   "tools/root_evidence/repeated_hesse.cpp third_party/root_ref_build/minuit2/lib/libMinuit2.a "
   "third_party/root_ref_build/minuit2/lib/libMinuit2Math.a -lpthread "
   "-o target/nfcn-accounting-4/repeated_hesse_root";

// One matrix as rows of numbers, one row per line.
static void matrixField(char const *name, ROOT::Minuit2::MnAlgebraicSymMatrix const &matrix)
{
   std::cout << ",\n    " << jsonString(name) << ": [";
   for (unsigned int i = 0; i < matrix.Nrow(); ++i) {
      std::cout << (i ? "," : "") << "\n      [";
      for (unsigned int j = 0; j < matrix.Nrow(); ++j) std::cout << (j ? ", " : "") << jsonNumber(matrix(i, j));
      std::cout << "]";
   }
   std::cout << "\n    ]";
}

int main(int argc, char **)
{
   if (argc != 1) {
      std::cerr << "usage: repeated_hesse_root\n";
      return 1;
   }
   Rosenbrock2 rosenbrock;
   MnUserParameters start;
   start.Add("x", -1.2, 0.2);
   start.Add("y", 1.0, 0.2);
   auto minimum = MnMigrad(rosenbrock, start, MnStrategy{2})(100000, 0.1);
   beginRecord(provenanceJson(program, buildCommand, "target/nfcn-accounting-4/repeated_hesse_root"),
               "rosenbrock2_strategy2_repeated_hesse");
   field("migrad_valid", minimum.IsValid());
   field("migrad_nfcn", minimum.NFcn());
   field("migrad_fval", minimum.Fval());
   field("migrad_x", minimum.UserState().Value(0));
   field("migrad_y", minimum.UserState().Value(1));
   matrixField("migrad_inverse_hessian", minimum.Error().InvHessian());
   // The explicit pass updates the minimum in place and appends a state.
   MnHesse(2)(rosenbrock, minimum, 100000);
   field("hesse_valid", minimum.IsValid());
   field("hesse_nfcn", minimum.NFcn());
   matrixField("hesse_inverse_hessian", minimum.Error().InvHessian());
   endRecord();
   return std::cout ? 0 : 1;
}
