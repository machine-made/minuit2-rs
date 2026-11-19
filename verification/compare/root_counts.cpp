// Count calls inside the ROOT objective, using the existing Minuit2 libraries.
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <iostream>

class Quadratic : public ROOT::Minuit2::FCNBase {
public:
    mutable unsigned calls = 0;
    double operator()(const std::vector<double> &p) const override {
        ++calls;
        return (p[0]-1)*(p[0]-1)+4*(p[1]+2)*(p[1]+2)+0.3*p[0]*p[1];
    }
    double Up() const override { return 1; }
};

int main() {
    Quadratic f;
    ROOT::Minuit2::MnUserParameters parameters;
    parameters.Add("x", 0.4, 0.1);
    parameters.Add("y", -1.0, 0.1);
    ROOT::Minuit2::MnMigrad migrad(f, parameters);
    auto minimum = migrad(0, 0.1);
    unsigned before = f.calls;
    ROOT::Minuit2::MnHesse hesse;
    hesse(f, minimum);
    std::cout << "{\"migrad\":" << before << ",\"total\":" << f.calls
              << ",\"hesse\":" << f.calls-before << "}\n";
}
