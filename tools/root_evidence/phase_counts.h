// Evidence-only stack classification for ROOT Minuit2 v6-36-08.
// Build with -rdynamic so backtrace_symbols exposes the Minuit2 entry points.
// The first numerical Hesse call is its centre. Later calls change one
// coordinate (diagonal) or two (mixed); the refinement has its own stack frame.
#pragma once
#include <execinfo.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct PhaseCounts {
   unsigned int calls = 0;
   unsigned int passes = 0;
   unsigned int beforeHesse = 0;
   bool inHesse = false;
   bool standalone = false;
   bool trace = false;
   std::vector<double> centre;
   std::map<std::string, unsigned int> phases;
   std::vector<std::string> evaluations;

   void beginStandalone() { standalone = true; inHesse = false; }

   void record(std::vector<double> const &p)
   {
      void *frames[64];
      int size = backtrace(frames, 64);
      char **symbols = backtrace_symbols(frames, size);
      if (!symbols) throw std::runtime_error("backtrace_symbols failed");
      std::string stack;
      for (int i = 0; i < size; ++i) stack += std::string(symbols[i]) + '\n';
      std::free(symbols);
      bool hesse = stack.find("MnHesse") != std::string::npos;
      std::string phase;
      if (hesse) {
         if (!inHesse) {
            if (passes == 0) beforeHesse = calls;
            ++passes;
            centre = p;
            phase = "centre";
         } else if (stack.find("HessianGradientCalculator") != std::string::npos) {
            phase = "refinement";
         } else {
            unsigned int changed = 0;
            for (std::size_t i = 0; i < p.size(); ++i) changed += p[i] != centre[i];
            if (changed != 1 && changed != 2) throw std::runtime_error("unclassified Hesse call");
            phase = changed == 1 ? "diagonal" : "mixed";
         }
         phase = (standalone ? "standalone_hesse_" : "hesse_") + phase;
      } else if (passes > 0) {
         phase = "continuation";
      } else {
         phase = stack.find("MnSeedGenerator") != std::string::npos ? "seed" : "migrad_iteration";
      }
      inHesse = hesse;
      ++calls;
      ++phases[phase];
      if (std::getenv("M2_TRACE_STACKS")) std::cerr << "CALL " << calls << '\n' << stack;
      if (trace) {
         std::ostringstream out;
         out << std::setprecision(17) << "{\"call\":" << calls << ",\"phase\":\"" << phase << "\",\"p\":[";
         for (std::size_t i = 0; i < p.size(); ++i) { if (i) out << ','; out << p[i]; }
         out << "]}";
         evaluations.push_back(out.str());
      }
   }

   void print() const
   {
      std::cout << ",\"user_calls\":" << calls << ",\"hesse_passes\":" << passes
                << ",\"before_hesse\":" << beforeHesse << ",\"phases\":{";
      bool comma = false;
      for (auto const &[phase, count] : phases) {
         if (comma) std::cout << ',';
         std::cout << '"' << phase << "\":" << count;
         comma = true;
      }
      std::cout << '}';
      if (trace) {
         std::cout << ",\"trace\":[";
         for (std::size_t i = 0; i < evaluations.size(); ++i) { if (i) std::cout << ','; std::cout << evaluations[i]; }
         std::cout << ']';
      }
   }
};
