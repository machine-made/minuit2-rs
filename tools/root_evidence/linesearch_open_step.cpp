// ROOT Minuit2 tag v6-36-08: direct MnLineSearch evidence.
// From the repository root, with an existing ROOT reference build:
//   evidence=$(mktemp -d)
//   g++ -O2 -std=c++17 -Ithird_party/root_ref/math/minuit2/inc \
//     tools/root_evidence/linesearch_open_step.cpp \
//     third_party/root_ref_build/minuit2/lib/libMinuit2.a \
//     third_party/root_ref_build/minuit2/lib/libMinuit2Math.a \
//     -lpthread -o "$evidence/linesearch_open_step"
//   "$evidence/linesearch_open_step" > "$evidence/linesearch_open_step.root.json"
// This writes the evidence stored at
// reports/verification/raw/trace/linesearch_open_step.root.json.
// Compare with cmp before copying it there.

#include <Minuit2/FCNBase.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MnFcn.h>
#include <Minuit2/MnLineSearch.h>
#include <Minuit2/MnMachinePrecision.h>
#include <Minuit2/MnParabolaPoint.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace ROOT::Minuit2;

struct Quartic final : FCNBase {
   double operator()(std::vector<double> const &p) const override
   {
      double const x = p[0];
      return std::pow(x, 4) / 4 - 10 * std::pow(x, 3) / 3 - 23 * x * x / 2 - 12 * x;
   }
   double Up() const override { return 1.0; }
};

int main()
{
   Quartic objective;
   MnFcn fcn(objective);
   MnAlgebraicVector x(1), step(1);
   x(0) = 0;
   step(0) = 1;
   double const gdel = -12;
   auto const result = MnLineSearch()(fcn, MinimumParameters(x, 0), step, gdel, MnMachinePrecision());
   std::cout << std::setprecision(17)
             << "{\n  \"root_version\": \"v6-36-08\",\n"
                "  \"source\": \"math/minuit2/src/MnLineSearch.cxx\",\n"
                "  \"start\": [\n    "
             << x(0) << "\n  ],\n  \"step\": [\n    " << step(0) << "\n  ],\n  \"gdel\": " << gdel << ",\n"
             << "  \"lambda\": " << result.X() << ",\n  \"fval\": " << result.Y()
             << ",\n  \"nfcn\": " << fcn.NumOfCalls() << "\n}\n";
   return std::cout ? 0 : 1;
}
