// Evidence-only stack classification for ROOT Minuit2 v6-36-08.
// This tool targets glibc on Linux: it parses the backtrace_symbols layout of glibc
// and stops with a message on any other layout.
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
#include <type_traits>
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

// JSON string literal for any text: quote, backslash and control characters escaped.
inline std::string jsonString(std::string const &text)
{
   std::string out = "\"";
   for (unsigned char c : text) {
      if (c == '"') out += "\\\"";
      else if (c == '\\') out += "\\\\";
      else if (c < 0x20) {
         char escape[8];
         std::snprintf(escape, sizeof escape, "\\u%04X", static_cast<unsigned int>(c));
         out += escape;
      } else out += static_cast<char>(c);
   }
   return out + '"';
}

// Provenance block of a stored file: the program, the ROOT tag and the script that
// built the archives, then the commands that build and run the program.
inline std::string provenanceJson(std::string const &program, std::string const &buildCommand,
                                  std::string const &runCommand)
{
   return "{\n    \"program\": " + jsonString(program) + ",\n    \"root_tag\": \"v6-36-08\""
          ",\n    \"library_build_command\": \"scripts/build_root_reference_runner.sh v6-36-08\""
          ",\n    \"build_command\": " + jsonString(buildCommand) +
          ",\n    \"command\": " + jsonString(runCommand) + "\n  }";
}

// A record is the whole output of one run, in the layout of the stored files:
// two-space indent, one key per line, a trailing newline. The caller prints the
// result fields, then PhaseCounts::print(), then endRecord().
inline void beginRecord(std::string const &provenance, std::string const &workload)
{
   std::cout << "{\n  \"provenance\": " << provenance << ",\n  \"result\": {\n    \"workload\": "
             << jsonString(workload);
}

// Fields take integers, bool or double only, so no number reaches the output unchecked.
template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0> void field(char const *name, T value)
{
   std::cout << ",\n    " << jsonString(name) << ": " << value;
}

inline void field(char const *name, bool value)
{
   std::cout << ",\n    " << jsonString(name) << ": " << (value ? "true" : "false");
}

inline void field(char const *name, double value)
{
   std::cout << ",\n    " << jsonString(name) << ": " << jsonNumber(value);
}

inline void endRecord() { std::cout << "\n  }\n}\n"; }

// True for the glibc layouts "file(symbol+0xoffset) [0xaddress]", "file(+0xoffset) [0xaddress]"
// and "[0xaddress]". Frames of libc itself often have no symbol name ("libc.so.6(+0x2a1ca)").
inline bool glibcFrame(std::string const &frame)
{
   if (frame.size() < 6 || frame.back() != ']') return false;
   if (frame.compare(0, 3, "[0x") == 0) return true;
   auto const close = frame.rfind(") [0x");
   return close != std::string::npos && frame.rfind('(', close) != std::string::npos;
}

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
      std::vector<std::string> frameTexts(symbols, symbols + size);
      std::free(symbols);
      std::string stack, names;
      for (std::string const &frame : frameTexts) {
         if (!glibcFrame(frame))
            throw std::runtime_error("unsupported backtrace_symbols layout; this tool needs glibc: " + frame);
         stack += frame + '\n';
         names += frameSymbol(frame) + '\n';
      }
      // The first frame is this harness, in the executable: unnamed means no -rdynamic.
      if (frameTexts.empty() || frameSymbol(frameTexts[0]).empty())
         throw std::runtime_error("the first frame of the call stack has no symbol name; build with -rdynamic");
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
         out << "{\"call\":" << calls << ",\"phase\":" << jsonString(phase) << ",\"p\":[";
         for (std::size_t i = 0; i < p.size(); ++i) { if (i) out << ','; out << jsonNumber(p[i]); }
         out << "]}";
         evaluations.push_back(out.str());
      }
   }

   // The invariant of the call-count change: Minuit2 reports exactly the calls the user function saw.
   bool reportedMatches(int reportedNfcn) const
   {
      if (reportedNfcn >= 0 && static_cast<unsigned int>(reportedNfcn) == calls) return true;
      std::cerr << "reported_nfcn " << reportedNfcn << " differs from user_calls " << calls << '\n';
      return false;
   }

   void print() const
   {
      std::cout << ",\n    \"user_calls\": " << calls << ",\n    \"hesse_passes\": " << passes
                << ",\n    \"before_hesse\": " << beforeHesse << ",\n    \"phases\": {";
      bool comma = false;
      for (auto const &[phase, count] : phases) {
         std::cout << (comma ? "," : "") << "\n      " << jsonString(phase) << ": " << count;
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
