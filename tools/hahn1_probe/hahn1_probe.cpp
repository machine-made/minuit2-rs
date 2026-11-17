// Hahn1 fit with ROOT Minuit2, tag v6-36-08.
//
// This program fits the NIST Hahn1 model from NIST start 2 with ROOT's Migrad.
// It uses strategy 1, an initial error of 0.1 per parameter, and a call limit of 1,000,000.
// The first argument is the EDM tolerance and defaults to 0.1.
// The output is the validity flag, the minimum, the EDM, the call count and the seven parameters.
// Each parameter line also shows the certified value and the relative difference.
//
// The program opens the file hahn1.txt in the current directory.
// In this repository that file is hahn1_dataset.txt: copy it to hahn1.txt before you run.
// Line 1 holds the number of points n.
// The next n lines hold x and y.
// The next line holds the seven start-2 values and the line after it holds the seven certified values.
//
// Build with g++ -O2 -std=c++17 -w.
// Add the include directories math/minuit2/inc, math/mathcore/inc and core/foundation/inc of a ROOT v6-36-08 source tree.
// Add the include directory inc and the build directory of a ROOT Minuit2 build of the same tag.
// Link the static libraries libMinuit2.a and libMinuit2Math.a from that build, then -lpthread.
// In this repository the source tree is third_party/root_ref and the build is third_party/root_ref_build/minuit2.
//
// The ROOT output at tolerance 0.1 is stored in
// reports/verification/raw/ref/hahn1_core_start2_migrad_root.txt.
// The test tests/hahn1_core_parity.rs pins the port to that output.
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/FunctionMinimum.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace ROOT::Minuit2;
struct F : FCNBase {
  std::vector<double> x, y;
  double operator()(const std::vector<double>& p) const override {
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
int main(int argc, char** argv) { double tol = argc > 1 ? atof(argv[1]) : 0.1;
  FILE* f = fopen("hahn1.txt","r"); int n; fscanf(f,"%d",&n);
  F fcn; fcn.x.resize(n); fcn.y.resize(n);
  for (int i=0;i<n;i++) fscanf(f,"%lf %lf",&fcn.x[i],&fcn.y[i]);
  std::vector<double> s(7), c(7);
  for (int i=0;i<7;i++) fscanf(f,"%lf",&s[i]);
  for (int i=0;i<7;i++) fscanf(f,"%lf",&c[i]);
  MnUserParameters up;
  for (int i=0;i<7;i++) { char nm[8]; snprintf(nm,8,"b%d",i+1); up.Add(nm, s[i], 0.1); }
  MnMigrad migrad(fcn, up, MnStrategy(1));
  FunctionMinimum min = migrad(1000000, tol);
  printf("valid=%d fval=%.10g edm=%.6g nfcn=%d\n", min.IsValid(), min.Fval(), min.Edm(), min.NFcn());
  for (int i=0;i<7;i++) {
    double v = min.UserParameters().Params()[i];
    printf("b%d %.10e cert %.10e rel %.3e\n", i+1, v, c[i], std::fabs(v-c[i])/std::fabs(c[i]));
  }
}
