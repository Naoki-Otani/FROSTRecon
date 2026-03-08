// mppc_correction.cpp
// Standalone tool to apply MPPC effects correction to the "wls" TTree in a ROOT file,
// and then apply ADC saturation correction to the resulting light yield.
//
// What it does:
//   - Reads an input ROOT file.
//   - Copies all TTrees except the target tree (default: "wls") as-is into the output ROOT file.
//   - For the target tree (default: "wls"):
//       * Reads photonx1..132 and photony1..140 (Double_t assumed).
//       * Applies noise + crosstalk + pixel saturation to get "light yield" in p.e.
//       * Applies ADC saturation correction ONLY when light yield >= 80 p.e.:
//           y = p0 + p1 * (1 - exp( - (x/p2)^p3 ))
//         Otherwise, y = x.
//       * Writes a new tree that:
//           - Keeps all original branches except photonx*/photony* (dropped).
//           - Adds new branches: lightyieldx[132], lightyieldy[140] (Double_t) as FINAL yields.
//
// Build:
//   g++ -O2 -std=c++17 mppc_correction.cpp $(root-config --cflags --libs) -o mppc_correction
//
// Usage:
//   ./mppc_correction -i input.root -o output.root [options]
//
// Options:
//   --tree <name>          Target tree name (default: wls)
//   --noise-rate <Hz>      Noise rate in Hz (default: 5.0e5)
//   --crosstalk <prob>     Crosstalk probability (default: 0.07)
//   --npixel <N>           Number of MPPC pixels (default: 559)
//   --int-range <sec>      Integration range in seconds (default: 1.0e-6)
//
//   --sat-p0 <val>         ADC saturation parameter p0 (default: -3967)
//   --sat-p1 <val>         ADC saturation parameter p1 (default: 6531)
//   --sat-p2 <val>         ADC saturation parameter p2 (default: 325.3)
//   --sat-p3 <val>         ADC saturation parameter p3 (default: 0.0241)
//   --sat-threshold <val>  Apply saturation only for x >= threshold (default: 80)
//
//   -h, --help             Show help

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>

// ---------------- Geometry ----------------
static constexpr int NFIBX = 132;
static constexpr int NFIBY = 140;

// ---------------- Parameters (defaults) ----------------
struct Params {
  std::string inPath;
  std::string outPath;
  std::string treeName = "wls";

  double integralRange = 1.0e-6; // [s]
  int    nPixel        = 559;

  double noiseRate     = 5.0e5;  // [Hz]
  double crosstalkProb = 0.07;

  // ADC saturation correction parameters
  double sat_p0 = -3967.0;
  double sat_p1 =  6531.0;
  double sat_p2 =   325.3;
  double sat_p3 =     0.0241;

  double sat_threshold = 80.0; // Apply saturation only if x >= threshold
};

static void PrintUsage(const char* prog) {
  std::cerr
    << "Usage:\n"
    << "  " << prog << " -i <input.root> -o <output.root> [options]\n\n"
    << "Options:\n"
    << "  --tree <name>          Target tree name (default: wls)\n"
    << "  --noise-rate <Hz>      Noise rate in Hz (default: 5.0e5)\n"
    << "  --crosstalk <prob>     Crosstalk probability (default: 0.07)\n"
    << "  --npixel <N>           Number of MPPC pixels (default: 559)\n"
    << "  --int-range <sec>      Integration range in seconds (default: 1.0e-6)\n"
    << "\n"
    << "  --sat-p0 <val>         ADC saturation parameter p0 (default: -3967)\n"
    << "  --sat-p1 <val>         ADC saturation parameter p1 (default: 6531)\n"
    << "  --sat-p2 <val>         ADC saturation parameter p2 (default: 325.3)\n"
    << "  --sat-p3 <val>         ADC saturation parameter p3 (default: 0.0241)\n"
    << "  --sat-threshold <val>  Apply saturation only for x >= threshold (default: 80)\n"
    << "\n"
    << "  -h, --help             Show this help\n";
}

static bool ParseArgs(int argc, char** argv, Params& p) {
  if (argc <= 1) return false;

  auto needValue = [&](int i) {
    if (i + 1 >= argc) {
      std::cerr << "[ERROR] Missing value after: " << argv[i] << "\n";
      return false;
    }
    return true;
  };

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];

    if (a == "-h" || a == "--help") {
      return false;
    } else if (a == "-i") {
      if (!needValue(i)) return false;
      p.inPath = argv[++i];
    } else if (a == "-o") {
      if (!needValue(i)) return false;
      p.outPath = argv[++i];
    } else if (a == "--tree") {
      if (!needValue(i)) return false;
      p.treeName = argv[++i];
    } else if (a == "--noise-rate") {
      if (!needValue(i)) return false;
      p.noiseRate = std::stod(argv[++i]);
    } else if (a == "--crosstalk") {
      if (!needValue(i)) return false;
      p.crosstalkProb = std::stod(argv[++i]);
    } else if (a == "--npixel") {
      if (!needValue(i)) return false;
      p.nPixel = std::stoi(argv[++i]);
    } else if (a == "--int-range") {
      if (!needValue(i)) return false;
      p.integralRange = std::stod(argv[++i]);
    } else if (a == "--sat-p0") {
      if (!needValue(i)) return false;
      p.sat_p0 = std::stod(argv[++i]);
    } else if (a == "--sat-p1") {
      if (!needValue(i)) return false;
      p.sat_p1 = std::stod(argv[++i]);
    } else if (a == "--sat-p2") {
      if (!needValue(i)) return false;
      p.sat_p2 = std::stod(argv[++i]);
    } else if (a == "--sat-p3") {
      if (!needValue(i)) return false;
      p.sat_p3 = std::stod(argv[++i]);
    } else if (a == "--sat-threshold") {
      if (!needValue(i)) return false;
      p.sat_threshold = std::stod(argv[++i]);
    } else {
      std::cerr << "[ERROR] Unknown argument: " << a << "\n";
      return false;
    }
  }

  if (p.inPath.empty() || p.outPath.empty()) {
    std::cerr << "[ERROR] -i and -o are required.\n";
    return false;
  }
  return true;
}

// Convert a (non-negative) double count into an integer hit count
static inline int ToNonNegativeInt(double x) {
  if (!std::isfinite(x) || x <= 0.0) return 0;
  return (int)std::llround(x);
}

// Apply noise (Poisson) + crosstalk (Binomial)
static inline int ApplyNoiseAndCrosstalk(TRandom3& r, int nIn, double noiseRate, double integralRange, double xtProb) {
  const int nNoise = r.Poisson(noiseRate * integralRange);
  int n = nIn + nNoise;
  if (n < 0) n = 0;

  const int nXT = r.Binomial(n, xtProb);
  n += nXT;
  return n;
}

// Apply MPPC pixel saturation effect by counting unique fired pixels
static int ApplyPixelSaturation(TRandom3& r, int nPE, int nPixel) {
  if (nPE <= 0) return 0;

  std::vector<int> pix;
  pix.reserve((size_t)nPE);

  for (int i = 0; i < nPE; ++i) {
    pix.push_back((int)r.Integer(nPixel)); // 0..nPixel-1
  }

  std::sort(pix.begin(), pix.end());
  const auto it = std::unique(pix.begin(), pix.end());
  return (int)std::distance(pix.begin(), it);
}

// ADC saturation correction:
// y = p0 + p1*(1 - exp(-(x/p2)^p3)), applied only when x >= threshold.
// For x < threshold, y = x.
static inline double ApplyAdcSaturation(double x, const Params& p) {
  if (x < p.sat_threshold) return x;

  const double t = std::pow(x / p.sat_p2, p.sat_p3);
  const double y = p.sat_p0 + p.sat_p1 * (1.0 - std::exp(-t));
  return y;
}

static uint32_t MakeSeedFromTimeAndName(const std::string& inPath) {
  uint32_t s = (uint32_t)time(nullptr);
  for (unsigned char c : inPath) s = s * 131u + (uint32_t)c;
  return s;
}

// Process the target tree: drop photon branches, add corrected lightyield branches.
static bool ProcessTargetTree(TFile& fin, TFile& fout, const Params& p) {
  TTree* tin = (TTree*)fin.Get(p.treeName.c_str());
  if (!tin) {
    std::cerr << "[ERROR] Target TTree '" << p.treeName << "' not found in input.\n";
    return false;
  }

  tin->SetBranchStatus("*", 1);

  // Input photon branches
  double photonx_in[NFIBX] = {0.0};
  double photony_in[NFIBY] = {0.0};

  for (int i = 0; i < NFIBX; ++i) {
    const std::string bn = "photonx" + std::to_string(i + 1);
    if (!tin->GetBranch(bn.c_str())) {
      std::cerr << "[ERROR] Missing branch: " << bn << " in tree '" << p.treeName << "'\n";
      return false;
    }
    tin->SetBranchAddress(bn.c_str(), &photonx_in[i]);
  }
  for (int i = 0; i < NFIBY; ++i) {
    const std::string bn = "photony" + std::to_string(i + 1);
    if (!tin->GetBranch(bn.c_str())) {
      std::cerr << "[ERROR] Missing branch: " << bn << " in tree '" << p.treeName << "'\n";
      return false;
    }
    tin->SetBranchAddress(bn.c_str(), &photony_in[i]);
  }

  // Disable photon branches before cloning so they are not in the output tree.
  for (int i = 0; i < NFIBX; ++i) tin->SetBranchStatus(Form("photonx%d", i + 1), 0);
  for (int i = 0; i < NFIBY; ++i) tin->SetBranchStatus(Form("photony%d", i + 1), 0);

  TTree* tout = tin->CloneTree(0);
  tout->SetName(p.treeName.c_str());
  tout->SetDirectory(&fout);

  // Output branches (FINAL yields after ADC saturation correction)
  double lightyieldx[NFIBX] = {0.0};
  double lightyieldy[NFIBY] = {0.0};
  tout->Branch("lightyieldx", lightyieldx, Form("lightyieldx[%d]/D", NFIBX));
  tout->Branch("lightyieldy", lightyieldy, Form("lightyieldy[%d]/D", NFIBY));

  // Re-enable photon branches for reading entries
  for (int i = 0; i < NFIBX; ++i) tin->SetBranchStatus(Form("photonx%d", i + 1), 1);
  for (int i = 0; i < NFIBY; ++i) tin->SetBranchStatus(Form("photony%d", i + 1), 1);

  // Random generator seeded from time + input name
  const uint32_t seed = MakeSeedFromTimeAndName(p.inPath);
  TRandom3 r(seed);

  const Long64_t n = tin->GetEntries();
  std::cout << "[INFO] Processing tree '" << p.treeName << "' entries=" << n << "\n";

  for (Long64_t ie = 0; ie < n; ++ie) {
    tin->GetEntry(ie);

    // X fibers
    for (int ix = 0; ix < NFIBX; ++ix) {
      const int nIn = ToNonNegativeInt(photonx_in[ix]);
      const int nAfter = ApplyNoiseAndCrosstalk(r, nIn, p.noiseRate, p.integralRange, p.crosstalkProb);
      const double ly = (double)ApplyPixelSaturation(r, nAfter, p.nPixel); // MPPC-corrected yield
      lightyieldx[ix] = ApplyAdcSaturation(ly, p); // final yield
    }

    // Y fibers
    for (int iy = 0; iy < NFIBY; ++iy) {
      const int nIn = ToNonNegativeInt(photony_in[iy]);
      const int nAfter = ApplyNoiseAndCrosstalk(r, nIn, p.noiseRate, p.integralRange, p.crosstalkProb);
      const double ly = (double)ApplyPixelSaturation(r, nAfter, p.nPixel);
      lightyieldy[iy] = ApplyAdcSaturation(ly, p);
    }

    tout->Fill();

    if (ie % 2000 == 0) {
      std::cout << "  " << ie << " / " << n << "\n";
    }
  }

  fout.cd();
  tout->Write("", TObject::kOverwrite);
  std::cout << "[INFO] Finished '" << p.treeName << "'\n";
  return true;
}

// Copy all TTrees except the target one as-is, and copy other objects as-is.
static bool CopyOtherObjects(TFile& fin, TFile& fout, const Params& p) {
  TIter nextKey(fin.GetListOfKeys());
  TKey* key = nullptr;

  while ((key = (TKey*)nextKey())) {
    const std::string name = key->GetName();

    // Skip target tree here; it will be handled separately.
    if (name == p.treeName) continue;

    TObject* obj = key->ReadObj();
    if (!obj) continue;

    if (obj->InheritsFrom(TTree::Class())) {
      TTree* t = (TTree*)obj;
      fout.cd();
      TTree* tclone = t->CloneTree(-1, "fast");
      tclone->SetName(t->GetName());
      tclone->Write("", TObject::kOverwrite);
      delete tclone;
      delete obj;
      continue;
    }

    // Other objects (histograms, etc.)
    fout.cd();
    obj->Write(name.c_str(), TObject::kOverwrite);
    delete obj;
  }
  return true;
}

int main(int argc, char** argv) {
  Params p;
  if (!ParseArgs(argc, argv, p)) {
    PrintUsage(argv[0]);
    return 2;
  }

  // Parameter sanity checks
  if (p.integralRange <= 0.0 || p.nPixel <= 0 ||
      p.noiseRate < 0.0 || p.crosstalkProb < 0.0 || p.crosstalkProb > 1.0 ||
      p.sat_threshold < 0.0) {
    std::cerr << "[ERROR] Invalid parameter values.\n";
    return 2;
  }

  TFile fin(p.inPath.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "[ERROR] Failed to open input: " << p.inPath << "\n";
    return 1;
  }

  // Create output directory if needed
  {
    const std::string outDir = gSystem->DirName(p.outPath.c_str());
    if (!outDir.empty()) gSystem->mkdir(outDir.c_str(), true);
  }

  TFile fout(p.outPath.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Failed to create output: " << p.outPath << "\n";
    return 1;
  }

  std::cout << "[INFO] Input : " << p.inPath << "\n";
  std::cout << "[INFO] Output: " << p.outPath << "\n";
  std::cout << "[INFO] Tree  : " << p.treeName << "\n";
  std::cout << "[INFO] Params: intRange=" << p.integralRange
            << " noiseRate=" << p.noiseRate
            << " xtProb=" << p.crosstalkProb
            << " nPixel=" << p.nPixel << "\n";
  std::cout << "[INFO] ADC sat: p0=" << p.sat_p0
            << " p1=" << p.sat_p1
            << " p2=" << p.sat_p2
            << " p3=" << p.sat_p3
            << " threshold=" << p.sat_threshold << "\n";

  if (!CopyOtherObjects(fin, fout, p)) {
    std::cerr << "[ERROR] Failed while copying non-target objects.\n";
    return 1;
  }

  if (!ProcessTargetTree(fin, fout, p)) {
    std::cerr << "[ERROR] Failed while processing target tree.\n";
    return 1;
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "[INFO] Done.\n";
  return 0;
}
