#pragma once

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace FROST {

// Very small command line helper. Supported forms:
//   --key value
//   --flag
//   --key=value
//
// Unknown/extra positional arguments are rejected to prevent silent mistakes.
class CLI {
 public:
  CLI(int argc, char** argv) { Parse(argc, argv); }

  bool Has(const std::string& key) const {
    return kv_.count(Norm(key)) > 0 || flags_.count(Norm(key)) > 0;
  }

  std::string Get(const std::string& key, const std::string& def = "") const {
    auto it = kv_.find(Norm(key));
    if (it == kv_.end()) return def;
    return it->second;
  }

  int GetInt(const std::string& key, int def) const {
    auto s = Get(key, "");
    if (s.empty()) return def;
    return std::stoi(s);
  }

  double GetDouble(const std::string& key, double def) const {
    auto s = Get(key, "");
    if (s.empty()) return def;
    return std::stod(s);
  }

  void Require(const std::string& key) const {
    if (!Has(key)) throw std::runtime_error("Missing required option: " + key);
  }

 private:
  static std::string Norm(const std::string& key) {
    if (key.rfind("--", 0) == 0) return key.substr(2);
    return key;
  }

  void Parse(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      std::string a(argv[i]);
      if (a.rfind("--", 0) != 0) {
        throw std::runtime_error("Unexpected positional argument: " + a);
      }

      // --key=value
      auto eq = a.find('=');
      if (eq != std::string::npos) {
        auto k = Norm(a.substr(0, eq));
        auto v = a.substr(eq + 1);
        kv_[k] = v;
        continue;
      }

      // --flag OR --key value
      auto k = Norm(a);
      if (i + 1 < argc) {
        std::string next(argv[i + 1]);
        if (next.rfind("--", 0) != 0) {
          kv_[k] = next;
          ++i;
          continue;
        }
      }
      flags_.insert(k);
    }
  }

  std::unordered_map<std::string, std::string> kv_;
  std::unordered_set<std::string> flags_;
};

}  // namespace FROST
