// ROOT Minuit2 tag v6-36-08: user-side call accounting evidence.
// From the repository root, with an existing ROOT reference build:
//   g++ -O2 -std=c++23 -g -rdynamic -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/nfcn_accounting.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o target/nfcn-accounting-4/nfcn_accounting_root
//   target/nfcn-accounting-4/nfcn_accounting_root --trace
// Reuse the runner's objectives without changing its stored reference results.
// Pass --trace to include every evaluation and its phase in the JSON.
// Set M2_TRACE_STACKS to include call stacks on stderr for phase identification.
#include "phase_counts.h"

#define main reference_runner_main
#include "../ref_runner_cpp/main.cpp"
#undef main

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
   std::cout << std::boolalpha;
   bool trace = argc > 1 && std::string(argv[1]) == "--trace";
   Counted<QuadraticNoG2> quadratic;
   quadratic.counts.trace = trace;
   MnUserParameters parameters;
   parameters.Add("x", 0.4, 0.1);
   parameters.Add("y", -1.0, 0.1);
   auto const minimum = MnMigrad(quadratic, parameters)();
   std::cout << "{\"workload\":\"quadratic2_no_g2_migrad\",\"reported_nfcn\":" << minimum.NFcn()
             << ",\"hessian_calls\":" << quadratic.hessianCalls
             << ",\"state_nfcn\":[";
   for (std::size_t i = 0; i < minimum.States().size(); ++i) {
      if (i != 0)
         std::cout << ',';
      std::cout << minimum.States()[i].NFcn();
   }
   std::cout << "]";
   quadratic.counts.print();
   std::cout << "}\n";

   Counted<Rosenbrock2> rosenbrock;
   rosenbrock.counts.trace = trace;
   MnUserParameters start;
   start.Add("x", 0.0, 0.1);
   start.Add("y", 0.0, 0.1);
   auto const continuation = MnMigrad(rosenbrock, start)(140, 0.1);
   std::cout << std::setprecision(17)
             << "{\"workload\":\"rosenbrock2_maxfcn140\",\"reported_nfcn\":" << continuation.NFcn()
             << ",\"valid\":" << continuation.IsValid()
             << ",\"call_limit\":" << continuation.HasReachedCallLimit()
             << ",\"above_edm\":" << continuation.IsAboveMaxEdm() << ",\"edm\":" << continuation.Edm();
   rosenbrock.counts.print();
   std::cout << "}\n";
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
   std::cout << "{\"workload\":\"quadratic2_minos_before_minos\",\"migrad_nfcn\":" << migradCalls
             << ",\"reported_nfcn\":" << minosMinimum.NFcn() << ",\"hessian_calls\":" << minos.hessianCalls;
   minos.counts.print();
   std::cout << "}\n";
   return std::cout ? 0 : 1;
}
