// merge_mcdata.C
// Merge ROOT files in a directory by grouping files that share the same base name
// before the last "_<number>.root" suffix.
//
// Example group:
// cosmicmuon..._20260205_0.root ... cosmicmuon..._20260205_20.root
// -> cosmicmuon..._20260205.root
//
// How to run:
// root -l
// .L merge_mcdata.C+
// merge_mcdata();

#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TFileMerger.h>

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>

// Return true if s is all digits
static bool isAllDigits(const std::string& s) {
  if (s.empty()) return false;
  for (char c : s) if (!std::isdigit((unsigned char)c)) return false;
  return true;
}

// If filename matches "..._<digits>.root", return base without "_<digits>.root" and the index.
// Otherwise return false.
static bool parseBaseAndIndex(const std::string& fname, std::string& base, int& idx) {
  const std::string ext = ".root";
  if (fname.size() <= ext.size()) return false;
  if (fname.substr(fname.size() - ext.size()) != ext) return false;

  const std::string noext = fname.substr(0, fname.size() - ext.size());
  const std::size_t pos = noext.rfind('_');
  if (pos == std::string::npos) return false;

  const std::string tail = noext.substr(pos + 1); // supposed digits
  if (!isAllDigits(tail)) return false;

  base = noext.substr(0, pos); // without last _digits
  idx = std::stoi(tail);
  return true;
}

// Main entry
void merge_mcdata(std::string inDirArg = "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/rawdata/x", std::string outDirArg = "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/mergedata/x") {
  const std::string inDir  = inDirArg;
  const std::string outDir = outDirArg;

  gSystem->mkdir(outDir.c_str(), true);

  // Scan directory
  TSystemDirectory dir("indir", inDir.c_str());
  TList* files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "[ERROR] Cannot list directory: " << inDir << "\n";
    return;
  }

  // base -> list of (idx, fullpath)
  std::map<std::string, std::vector<std::pair<int,std::string>>> groups;

  TIter next(files);
  while (TSystemFile* f = (TSystemFile*)next()) {
    const char* nameC = f->GetName();
    if (!nameC) continue;

    std::string fname = nameC;
    if (f->IsDirectory()) continue;
    if (fname == "." || fname == "..") continue;

    std::string base;
    int idx = -1;
    if (!parseBaseAndIndex(fname, base, idx)) continue; // only *_digits.root

    const std::string fullpath = inDir + "/" + fname;
    groups[base].push_back({idx, fullpath});
  }

  if (groups.empty()) {
    std::cerr << "[WARN] No matching files (*_digits.root) found in: " << inDir << "\n";
    return;
  }

  // Merge each group
  int nMerged = 0;
  for (auto& kv : groups) {
    const std::string& base = kv.first;
    auto& v = kv.second;

    // Sort by index
    std::sort(v.begin(), v.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    const std::string outPath = outDir + "/" + base + ".root";

    std::cout << "[INFO] Merging group: " << base
              << "  (files=" << v.size() << ") -> " << outPath << "\n";

    TFileMerger merger(/*isLocal=*/false, /*verbose=*/false);
    merger.OutputFile(outPath.c_str(), "RECREATE");

    for (const auto& it : v) {
      // it.first is idx, it.second is full path
      merger.AddFile(it.second.c_str());
    }

    const bool ok = merger.Merge();
    if (!ok) {
      std::cerr << "[ERROR] Merge failed for: " << outPath << "\n";
      continue;
    }

    ++nMerged;
  }

  std::cout << "[INFO] Done. Merged outputs: " << nMerged
            << "  -> " << outDir << "\n";
}
