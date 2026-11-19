// Evidence-only stack classification for ROOT Minuit2 v6-36-08.
// Build with -rdynamic so backtrace_symbols exposes the Minuit2 entry points.
// The first numerical Hesse call is its centre. Later calls change one
// coordinate (diagonal) or two (mixed); the refinement has its own stack frame.
// Frames are matched on their symbol names only, never on the executable path.
#pragma once
#include <execinfo.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// JSON has no NaN or infinity: write null for a non-finite number. Otherwise write
// the first of 15, 16 and 17 significant digits that reads back as the same double.
inline std::string jsonNumber(double value)
{
   if (!std::isfinite(value)) return "null";
   char text[32];
   for (int digits = 15; digits <= 17; ++digits) {
      std::snprintf(text, sizeof text, "%.*g", digits, value);
      if (std::strtod(text, nullptr) == value) break;
   }
   return text;
}

// A record is the whole output of one run, in the layout of the stored files:
// two-space indent, one key per line, a trailing newline. The caller prints the
// result fields, then PhaseCounts::print(), then endRecord().
inline void beginRecord(std::string const &provenance, std::string const &workload)
{
   std::cout << "{\n  \"provenance\": " << provenance << ",\n  \"result\": {\n    \"workload\": \"" << workload << '"';
}

template <class T> void field(char const *name, T const &value)
{
   std::cout << ",\n    \"" << name << "\": " << value;
}

inline void field(char const *name, bool value) { field(name, value ? "true" : "false"); }

inline void endRecord() { std::cout << "\n  }\n}\n"; }

// Symbol name of one backtrace_symbols frame, "file(symbol+0xoffset) [0xaddress]".
// The name is empty when the frame has no symbol; the file path is never part of it.
inline std::string frameSymbol(std::string const &frame)
{
   auto const close = frame.rfind(") [");
   if (close == std::string::npos) return {};
   auto const open = frame.rfind('(', close);
   if (open == std::string::npos) return {};
   std::string const inside = frame.substr(open + 1, close - open - 1);
   return inside.substr(0, inside.rfind('+'));
}

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
      std::string stack, names;
      for (int i = 0; i < size; ++i) {
         stack += std::string(symbols[i]) + '\n';
         names += frameSymbol(symbols[i]) + '\n';
      }
      std::free(symbols);
      bool hesse = names.find("MnHesse") != std::string::npos;
      std::string phase;
      if (hesse) {
         if (!inHesse) {
            if (passes == 0) beforeHesse = calls;
            ++passes;
            centre = p;
            phase = "centre";
         } else if (names.find("HessianGradientCalculator") != std::string::npos) {
            phase = "refinement";
         } else {
            unsigned int changed = 0;
            for (std::size_t i = 0; i < p.size(); ++i) changed += p[i] != centre[i];
            if (changed != 1 && changed != 2) {
               std::ostringstream message;
               message << std::setprecision(17) << "unclassified Hesse call " << calls + 1 << ": " << changed
                       << " changed coordinates at p = [";
               for (std::size_t i = 0; i < p.size(); ++i) message << (i ? "," : "") << p[i];
               throw std::runtime_error(message.str() + "]");
            }
            phase = changed == 1 ? "diagonal" : "mixed";
         }
         phase = (standalone ? "standalone_hesse_" : "hesse_") + phase;
      } else if (passes > 0) {
         phase = "continuation";
      } else {
         phase = names.find("MnSeedGenerator") != std::string::npos ? "seed" : "migrad_iteration";
      }
      inHesse = hesse;
      ++calls;
      ++phases[phase];
      if (std::getenv("M2_TRACE_STACKS")) std::cerr << "CALL " << calls << '\n' << stack;
      if (trace) {
         std::ostringstream out;
         out << "{\"call\":" << calls << ",\"phase\":\"" << phase << "\",\"p\":[";
         for (std::size_t i = 0; i < p.size(); ++i) { if (i) out << ','; out << jsonNumber(p[i]); }
         out << "]}";
         evaluations.push_back(out.str());
      }
   }

   void print() const
   {
      std::cout << ",\n    \"user_calls\": " << calls << ",\n    \"hesse_passes\": " << passes
                << ",\n    \"before_hesse\": " << beforeHesse << ",\n    \"phases\": {";
      bool comma = false;
      for (auto const &[phase, count] : phases) {
         std::cout << (comma ? "," : "") << "\n      \"" << phase << "\": " << count;
         comma = true;
      }
      std::cout << "\n    }";
      if (trace) {
         std::cout << ",\n    \"trace\": [";
         for (std::size_t i = 0; i < evaluations.size(); ++i) std::cout << (i ? "," : "") << "\n      " << evaluations[i];
         std::cout << "\n    ]";
      }
   }
};
