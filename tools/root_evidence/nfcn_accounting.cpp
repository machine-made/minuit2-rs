// ROOT Minuit2 tag v6-36-08: user-side call accounting evidence.
// From the repository root, with an existing ROOT reference build:
//   evidence=$(mktemp -d)
//   g++ -O2 -std=c++17 -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/nfcn_accounting.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o "$evidence/nfcn_accounting"
//   "$evidence/nfcn_accounting"
// Reuse the runner's objectives without changing its stored reference results.
// Pass --trace to print every Rosenbrock evaluation with 17 significant digits.
// Set M2_TRACE_STACKS to include call stacks on stderr for phase identification.
#include <execinfo.h>
#include <cstdlib>

#define main reference_runner_main
#include "../ref_runner_cpp/main.cpp"
#undef main

template <class Objective> class Counted final : public FCNBase {
public:
   mutable unsigned int calls = 0;
   mutable unsigned int hessianCalls = 0;
   bool trace = false;

   double operator()(std::vector<double> const &p) const override
   {
      ++calls;
      if (trace) {
         std::cout << std::setprecision(17) << "TRACE " << calls;
         for (double value : p)
            std::cout << ' ' << value;
         std::cout << '\n';
         if (std::getenv("M2_TRACE_STACKS")) {
            std::cerr << "CALL " << calls << '\n';
            void *frames[64];
            int size = backtrace(frames, 64);
            backtrace_symbols_fd(frames, size, 2);
         }
      }
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
   Counted<QuadraticNoG2> quadratic;
   MnUserParameters parameters;
   parameters.Add("x", 0.4, 0.1);
   parameters.Add("y", -1.0, 0.1);
   auto const minimum = MnMigrad(quadratic, parameters)();
   std::cout << "{\"workload\":\"quadratic2_no_g2_migrad\",\"reported_nfcn\":" << minimum.NFcn()
             << ",\"user_calls\":" << quadratic.calls << ",\"hessian_calls\":" << quadratic.hessianCalls
             << ",\"state_nfcn\":[";
   for (std::size_t i = 0; i < minimum.States().size(); ++i) {
      if (i != 0)
         std::cout << ',';
      std::cout << minimum.States()[i].NFcn();
   }
   std::cout << "]}\n";

   Counted<Rosenbrock2> rosenbrock;
   rosenbrock.trace = argc > 1 && std::string(argv[1]) == "--trace";
   MnUserParameters start;
   start.Add("x", 0.0, 0.1);
   start.Add("y", 0.0, 0.1);
   auto const continuation = MnMigrad(rosenbrock, start)(140, 0.1);
   std::cout << std::setprecision(17)
             << "{\"workload\":\"rosenbrock2_maxfcn140\",\"reported_nfcn\":" << continuation.NFcn()
             << ",\"user_calls\":" << rosenbrock.calls << ",\"valid\":" << continuation.IsValid()
             << ",\"call_limit\":" << continuation.HasReachedCallLimit()
             << ",\"above_edm\":" << continuation.IsAboveMaxEdm() << ",\"edm\":" << continuation.Edm() << "}\n";
   return std::cout ? 0 : 1;
}
