// mppccorrection.C
// ROOT macro to apply MPPC effects correction to MC merged data.
// - Reads all ROOT files in:
//   eg: /group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/mergedata/x
// - For each file, reads TTree "wls" which contains photonx1..132 and photony1..140 (Double_t)
// - Computes corrected light yields and stores them into new branches:
//   lightyieldx[132], lightyieldy[140]
// - Drops photonx*/photony* branches from the output tree (keeps all other branches as-is)
// - Writes output files into:
//   /group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/aftermppccorrection/x
// - Output name rule:
//   if input ends with _YYYYMMDD.root, remove that suffix and append _aftermppccorrection.root
//   else just replace .root -> _aftermppccorrection.root
//
// How to run:
//   root -l
//   .L mppc_correction_all.C+
//   make_mppc_correction_all();

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <dirent.h>
#include <sys/stat.h>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <regex>

// ---------------- Geometry ----------------
static const int NFIBX = 132;
static const int NFIBY = 140;

// ---------------- MPPC / digitizer settings (same as your old code) ----------------
static const double INTEGRAL_RANGE = 1.0e-6;   // integration range [s]
static const int    MPPC_NPIXEL    = 559;    // MPPC number of pixels

// Noise / crosstalk parameters (same as your old code)
static const double NOISE_RATE     = 5.0e5; // typical noise rate [Hz]
static const double CROSSTALK_PROB = 0.07;  // typical crosstalk prob

// ---------------- Utilities ----------------
static bool isRegularFile(const std::string& path) {
  struct stat st;
  if (stat(path.c_str(), &st) != 0) return false;
  return S_ISREG(st.st_mode);
}

static std::vector<std::string> listRootFiles(const std::string& dirPath) {
  std::vector<std::string> files;

  DIR* dp = opendir(dirPath.c_str());
  if (!dp) {
    std::cerr << "[ERROR] Cannot open directory: " << dirPath << "\n";
    return files;
  }

  while (true) {
    struct dirent* ent = readdir(dp);
    if (!ent) break;

    const std::string name = ent->d_name;
    if (name == "." || name == "..") continue;
    if (name.size() < 6) continue;
    if (name.substr(name.size() - 5) != ".root") continue;

    const std::string full = dirPath + "/" + name;
    if (!isRegularFile(full)) continue;

    files.push_back(full);
  }
  closedir(dp);

  std::sort(files.begin(), files.end());
  return files;
}

static std::string basenameOnly(const std::string& path) {
  const auto pos = path.find_last_of('/');
  if (pos == std::string::npos) return path;
  return path.substr(pos + 1);
}

static std::string makeOutputName(const std::string& inBase) {
  // If input matches: (anything)_(8digits).root -> (anything)_aftermppccorrection.root
  // Example: cosmic..._20260205.root -> cosmic..._aftermppccorrection.root
  static const std::regex reDateSuffix(R"(^(.*)_([0-9]{8})\.root$)");

  std::smatch m;
  if (std::regex_match(inBase, m, reDateSuffix)) {
    return m[1].str() + "_aftermppccorrection.root";
  }

  // Otherwise: replace .root -> _aftermppccorrection.root
  if (inBase.size() >= 5 && inBase.substr(inBase.size() - 5) == ".root") {
    return inBase.substr(0, inBase.size() - 5) + "_aftermppccorrection.root";
  }
  return inBase + "_aftermppccorrection.root";
}

// Convert a (non-negative) double count into an integer hit count
static inline int toNonNegativeInt(double x) {
  if (!std::isfinite(x) || x <= 0.0) return 0;
  // photon branches are effectively integer-like; rounding is reasonable here
  return (int)std::llround(x);
}

// Apply: add noise (Poisson) + crosstalk (Binomial)
static inline int applyNoiseAndCrosstalk(TRandom3& r, int nIn) {
  const int nNoise = r.Poisson(NOISE_RATE * INTEGRAL_RANGE);
  int n = nIn + nNoise;
  if (n < 0) n = 0;

  const int nXT = r.Binomial(n, CROSSTALK_PROB);
  n += nXT;
  return n;
}

// Apply MPPC pixel saturation effect by counting unique fired pixels
static int applyPixelSaturation(TRandom3& r, int nPhotoElectrons) {
  if (nPhotoElectrons <= 0) return 0;

  // Generate pixel IDs and count unique
  std::vector<int> pix;
  pix.reserve((size_t)nPhotoElectrons);

  for (int i = 0; i < nPhotoElectrons; ++i) {
    pix.push_back((int)r.Integer(MPPC_NPIXEL)); // 0..MPPC_NPIXEL-1
  }

  std::sort(pix.begin(), pix.end());
  const auto it = std::unique(pix.begin(), pix.end());
  return (int)std::distance(pix.begin(), it);
}

// ---------------- Core per-file processing ----------------
static bool processOneFile(const std::string& inPath, const std::string& outDir) {
  TFile fin(inPath.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "[ERROR] Failed to open input: " << inPath << "\n";
    return false;
  }

  TTree* tin = (TTree*)fin.Get("wls");
  if (!tin) {
    std::cerr << "[ERROR] TTree 'wls' not found in: " << inPath << "\n";
    return false;
  }

  const std::string inBase = basenameOnly(inPath);
  const std::string outBase = makeOutputName(inBase);
  const std::string outPath = outDir + "/" + outBase;

  // Make sure output directory exists
  gSystem->mkdir(outDir.c_str(), true);

  // ---- Prepare input reading ----
  // Enable all branches for safety, then set addresses for photon branches
  tin->SetBranchStatus("*", 1);

  double photonx_in[NFIBX];
  double photony_in[NFIBY];

  for (int i = 0; i < NFIBX; ++i) {
    photonx_in[i] = 0.0;
    tin->SetBranchAddress(Form("photonx%d", i + 1), &photonx_in[i]);
  }
  for (int i = 0; i < NFIBY; ++i) {
    photony_in[i] = 0.0;
    tin->SetBranchAddress(Form("photony%d", i + 1), &photony_in[i]);
  }

  // ---- Create output and clone tree structure WITHOUT photon branches ----
  // IMPORTANT:
  //   - We disable photonx*/photony* BEFORE CloneTree(0) so output does not contain them.
  //   - After cloning, we re-enable them for reading events (GetEntry).
  for (int i = 0; i < NFIBX; ++i) tin->SetBranchStatus(Form("photonx%d", i + 1), 0);
  for (int i = 0; i < NFIBY; ++i) tin->SetBranchStatus(Form("photony%d", i + 1), 0);

  TFile fout(outPath.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Failed to create output: " << outPath << "\n";
    return false;
  }

  TTree* tout = tin->CloneTree(0); // clone enabled branches only (photon* excluded)
  tout->SetDirectory(&fout);

  // Add corrected branches
  double lightyieldx[NFIBX];
  double lightyieldy[NFIBY];
  tout->Branch("lightyieldx", lightyieldx, Form("lightyieldx[%d]/D", NFIBX));
  tout->Branch("lightyieldy", lightyieldy, Form("lightyieldy[%d]/D", NFIBY));

  // Re-enable photon branches for event reading
  for (int i = 0; i < NFIBX; ++i) tin->SetBranchStatus(Form("photonx%d", i + 1), 1);
  for (int i = 0; i < NFIBY; ++i) tin->SetBranchStatus(Form("photony%d", i + 1), 1);

  // ---- Random generator ----
  // Seed with a mix of time and file name hash to avoid identical results across files.
  UInt_t seed = (UInt_t)time(nullptr);
  for (char c : inBase) seed = seed * 131u + (UInt_t)(unsigned char)c;
  TRandom3 r(seed);

  // ---- Event loop ----
  const Long64_t n = tin->GetEntries();
  std::cout << "[INFO] Processing: " << inBase << " entries=" << n << "\n";

  for (Long64_t ie = 0; ie < n; ++ie) {
    tin->GetEntry(ie);

    // Compute corrected yields for X
    for (int ix = 0; ix < NFIBX; ++ix) {
      const int nIn = toNonNegativeInt(photonx_in[ix]);
      const int nAfter = applyNoiseAndCrosstalk(r, nIn);
      lightyieldx[ix] = (double)applyPixelSaturation(r, nAfter);
    }

    // Compute corrected yields for Y
    for (int iy = 0; iy < NFIBY; ++iy) {
      const int nIn = toNonNegativeInt(photony_in[iy]);
      const int nAfter = applyNoiseAndCrosstalk(r, nIn);
      lightyieldy[iy] = (double)applyPixelSaturation(r, nAfter);
    }

    // Fill output: cloned branches (except photon*) + new lightyieldx/y
    tout->Fill();

    if (ie % 2000 == 0) {
      std::cout << "  " << ie << " / " << n << "\n";
    }
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "[INFO] Wrote: " << outPath << "\n";
  return true;
}

// ---------------- Macro entry ----------------
void mppccorrection(std::string inDirArg = "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/mergedata/x", std::string outDirArg = "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehitafter/mppccorrection/x") {
  const std::string inDir  = inDirArg;
  const std::string outDir = outDirArg;

  const auto files = listRootFiles(inDir);
  if (files.empty()) {
    std::cerr << "[WARN] No .root files found in: " << inDir << "\n";
    return;
  }

  int nOK = 0, nNG = 0;
  for (const auto& f : files) {
    const bool ok = processOneFile(f, outDir);
    if (ok) ++nOK;
    else    ++nNG;
  }

  std::cout << "[INFO] Done. success=" << nOK << " failed=" << nNG << "\n";
}
