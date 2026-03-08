#pragma once

#include <string>
#include <unordered_map>

namespace FROST {

inline void ReplaceAll(std::string& s, const std::string& from, const std::string& to) {
  if (from.empty()) return;
  size_t pos = 0;
  while ((pos = s.find(from, pos)) != std::string::npos) {
    s.replace(pos, from.size(), to);
    pos += to.size();
  }
}

// Simple template substitution using {name} placeholders.
inline std::string ApplyTemplate(std::string pattern,
                                 const std::unordered_map<std::string, std::string>& vars) {
  for (const auto& kv : vars) {
    ReplaceAll(pattern, "{" + kv.first + "}", kv.second);
  }
  return pattern;
}

}  // namespace FROST
