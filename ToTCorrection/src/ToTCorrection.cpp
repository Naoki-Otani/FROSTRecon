#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMath.h>
#include <TROOT.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

static const Int_t    NCH            = 272;
static const Int_t    NBUNCH         = 8;
static const Int_t    INT_RANGE      = 75;
static const Double_t BUNCH_INTERVAL = 43.5;

// -------------------------
// Helper: trim
// -------------------------
static inline std::string Trim(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) b++;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e-1]))) e--;
  return s.substr(b, e - b);
}

// -------------------------
// Read sampling_first_bunch rules from text file
// Format:
//   run_max  sampling_first_bunch
// with comments starting '#'
// -------------------------
static bool LoadSamplingRules(const std::string& rulePath,
                             std::vector<std::pair<int,int>>& rules)
{
  std::ifstream ifs(rulePath);
  if (!ifs) {
    std::cerr << "[ERROR] Cannot open rules file: " << rulePath << "\n";
    return false;
  }

  rules.clear();
  std::string line;
  while (std::getline(ifs, line)) {
    line = Trim(line);
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    int run_max = -1, sfb = -1;
    if (!(iss >> run_max >> sfb)) continue;
    rules.emplace_back(run_max, sfb);
  }

  // Ensure ascending by run_max
  std::sort(rules.begin(), rules.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });

  return !rules.empty();
}

// -------------------------
// Determine SAMPLING_FIRST_BUNCH for a given run using loaded rules
// For the first rule with run <= run_max, return its sampling index.
// -------------------------
static int GetSamplingFirstBunch(int run,
                                 const std::vector<std::pair<int,int>>& rules)
{
  for (const auto& r : rules) {
    if (run <= r.first) return r.second;
  }
  // If not found, use the last one (should not happen if last run_max is large)
  return rules.back().second;
}

// -------------------------
// Extract run number from filename
// Expect: run#####_#_####_lightyield.root (example: run00016_0_9999_lightyield.root)
// -------------------------
static bool ExtractRunNumber(const std::string& fname, int& run)
{
  // Capture "run" + 5 digits
  static const std::regex re(R"(run(\d{5})_.*_lightyield\.root$)");
  std::smatch m;
  if (!std::regex_search(fname, m, re)) return false;
  run = std::stoi(m[1].str());
  return true;
}

// -------------------------
// Find maximum value in vec that lies within [low, high].
// Return true if found.
// -------------------------
static bool MaxInWindow(const std::vector<double>& vec,
                        double low, double high,
                        double& outMax)
{
  bool found = false;
  double m = -std::numeric_limits<double>::infinity();
  for (double v : vec) {
    if (v >= low && v <= high) {
      found = true;
      if (v > m) m = v;
    }
  }
  if(vec.size()!=1) found=false;
  if (found) outMax = m;
  return found;
}

static bool MinInWindow(const std::vector<double>& vec,
                        double low, double high,
                        double& outMin)
{
  bool found = false;
  double m = std::numeric_limits<double>::infinity();
  for (double v : vec) {
    if (v >= low && v <= high) {
      found = true;
      if (v < m) m = v;
    }
  }
  if (found) outMin = m;
  return found;
}

// -------------------------
// Find t_max (sample index) where waveform is maximum within [Ts, Te].
// Ts/Te are in sampling-index units. We search waveform[ch][i] in that window.
// Return false if waveform is empty or the window is out of range.
// -------------------------
static bool FindTmaxInWaveform(const std::vector<double>& wf,
                              double Ts, double Te,
                              double& tmax_out)
{
  if (wf.empty()) return false;

  // Convert to inclusive integer sample window
  int iBeg = (int)std::ceil(std::min(Ts, Te));
  int iEnd = (int)std::floor(std::max(Ts, Te));
  if (iEnd < 0) return false;
  if (iBeg >= (int)wf.size()) return false;
  if (iBeg < 0) iBeg = 0;
  if (iEnd >= (int)wf.size()) iEnd = (int)wf.size() - 1;
  if (iBeg > iEnd) return false;

  double vmax = -std::numeric_limits<double>::infinity();
  int imax = iBeg;
  for (int i = iBeg; i <= iEnd; ++i) {
    if (wf[i] > vmax) {
      vmax = wf[i];
      imax = i;
    }
  }
  tmax_out = (double)imax;
  return true;
}



// Model: lightyield = A * exp(B*ToT + C*ToT^2)
static TF1* MakeFitFunc(const std::string& name, double xmin, double xmax)
{
  // TF1* f = new TF1(name.c_str(), "[0]*TMath::Exp([1]*x + [2]*x*x)", xmin, xmax);
  TF1* f = new TF1(name.c_str(), "[0]*(x/100.) + [1]*(x/100.)*(x/100.) + [2]*(x/100.)*(x/100.)*(x/100.)", xmin, xmax);
  f->SetParNames("A","B","C");
  return f;
}

// -------------------------
// Main
// -------------------------
int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <data directory>(eg:/group/nu/ninja/work/otani/FROST_beamdata/e71c)\n";
    return 1;
  }

  const std::string baseDir = argv[1];
  const std::string inputDir  = baseDir + "/rootfile_aftercalib";
  const std::string outDir    = baseDir + "/ToTcorrection";
  const std::string ruleFile  = baseDir + "/calibration/sampling_first_bunch/sampling_first_bunch_rules.txt";
  const std::string outPdf    = outDir + "/ToT_vs_LightYield_2D_allch.pdf";

  // Store cablenum label for each channel (filled when we first see it)
  std::vector<int> cableLabel(NCH, -1);

  // Create output directory if it doesn't exist
  gSystem->mkdir(outDir.c_str(), /*recursive=*/true);

  // Load sampling rules
  std::vector<std::pair<int,int>> samplingRules;
  if (!LoadSamplingRules(ruleFile, samplingRules)) {
    std::cerr << "[ERROR] Failed to load sampling rules.\n";
    return 2;
  }

  // Prepare histograms (one per channel)
  // Binning/range: adjust if needed.
  const int    nbinX = 200;   // ToT
  const double xMin  = 0.0;
  const double xMax  = 1500.0;
  const int    nbinY = 150;   // lightyield
  const double yMin  = 0.0;
  const double yMax  = 150.0;

  std::vector<std::unique_ptr<TH2D>> h2(NCH);
  for (int ch = 0; ch < NCH; ++ch) {
    std::ostringstream name, title;
    name  << "h2_ch" << ch;
    title << "ch " << ch << ";ToT (leading - trailing);light yield from ADC [p.e.]";
    h2[ch] = std::make_unique<TH2D>(name.str().c_str(), title.str().c_str(),
                                    nbinX, xMin, xMax,
                                    nbinY, yMin, yMax);
  }

  // List ROOT files in inputDir
  TSystemDirectory dir("inputDir", inputDir.c_str());
  TList* files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "[ERROR] Cannot list directory: " << inputDir << "\n";
    return 3;
  }

  // Match run#####_#_####_lightyield.root
  const std::regex fileRe(R"(run\d{5}_\d+_\d+_lightyield\.root$)");

  int nFilesUsed = 0;

  TIter next(files);
  while (TSystemFile* f = (TSystemFile*)next()) {
    std::string fname = f->GetName();
    if (f->IsDirectory()) continue;
    if (!std::regex_match(fname, fileRe)) continue;

    int run = -1;
    if (!ExtractRunNumber(fname, run)) continue;

    const int SAMPLING_FIRST_BUNCH = GetSamplingFirstBunch(run, samplingRules);

    const std::string fullPath = inputDir + "/" + fname;
    std::unique_ptr<TFile> tf(TFile::Open(fullPath.c_str(), "READ"));
    if (!tf || tf->IsZombie()) {
      std::cerr << "[WARN] Cannot open: " << fullPath << "\n";
      continue;
    }

    TTree* tree = dynamic_cast<TTree*>(tf->Get("tree"));
    if (!tree) {
      std::cerr << "[WARN] No TTree named 'tree' in: " << fullPath << "\n";
      continue;
    }

    // Disable all branches first to reduce I/O (waveform is huge).
    tree->SetBranchStatus("*", 0);

    // Branches
    Double_t lightyield[NCH][NBUNCH];
    Int_t    cablenum[NCH];
    std::vector<std::vector<double>>* leading  = nullptr;
    std::vector<std::vector<double>>* trailing = nullptr;
    std::vector<std::vector<double>>* waveform = nullptr;

    // Enable only lightweight branches for the first pass.
    tree->SetBranchStatus("lightyield", 1);
    tree->SetBranchStatus("cablenum",  1);
    tree->SetBranchStatus("leading",   1);
    tree->SetBranchStatus("trailing",  1);
    // Keep waveform disabled in normal GetEntry; read it only when needed.
    tree->SetBranchStatus("waveform",  0);

    tree->SetBranchAddress("lightyield", lightyield);
    tree->SetBranchAddress("cablenum",  cablenum);
    tree->SetBranchAddress("leading", &leading);
    tree->SetBranchAddress("trailing", &trailing);
    tree->SetBranchAddress("waveform", &waveform);

    // Pointer to waveform branch for selective reading.
    TBranch* brWaveform = tree->GetBranch("waveform");
    if (!brWaveform) {
      std::cerr << "[WARN] No branch named 'waveform' in: " << fullPath << "\n";
      continue;
    }

    const Long64_t nEnt = tree->GetEntries();
    for (Long64_t ie = 0; ie < nEnt; ++ie) {
      // Read only lightweight branches here (waveform is OFF).
      tree->GetEntry(ie);

      // At this point, waveform is not read yet. Do NOT check waveform here.
      if (!leading || !trailing) continue;
      if ((int)leading->size() < NCH)  continue;
      if ((int)trailing->size() < NCH) continue;

      // Quick pre-scan: if no lightyield>10 in this event, skip waveform read.
      bool needWaveform = false;
      for (int ch = 0; ch < NCH && !needWaveform; ++ch) {
        for (int bunch = 0; bunch < NBUNCH; ++bunch) {
          if (lightyield[ch][bunch] > 10.0) { needWaveform = true; break; }
        }
      }

      if (!needWaveform) continue;

      // std::cout<<"Event "<<ie<<" has lightyield > 10, reading waveform...\n";

      // Read waveform branch ONLY for this event.
      // Temporarily enable waveform; otherwise the pointer may remain null.
      tree->SetBranchStatus("waveform", 1);
      brWaveform->GetEntry(ie);
      tree->SetBranchStatus("waveform", 0);

      // Now waveform should be available.
      if (!waveform) continue;
      if ((int)waveform->size() < NCH) continue;
      std::cout<<"Event "<<ie<<" has valid leading/trailing/waveform vectors.\n";

      for (int ch = 0; ch < NCH; ++ch) {

        // Cache cablenum for this channel (use the first observed value)
        if (cableLabel[ch] < 0) {
          cableLabel[ch] = cablenum[ch];
        }
        const auto& vLead = leading->at(ch);
        const auto& vTrail = trailing->at(ch);
        const auto& wf = waveform->at(ch);

        // If either empty, no ToT can be formed for any bunch
        // (We still check per bunch window; empty means invalid.)
        for (int bunch = 0; bunch < NBUNCH; ++bunch) {

          const double ly = lightyield[ch][bunch];
          std::cout<<ly<<std::endl;
          if (ly <= 10.0) continue; // Requirement: lightyield > 10 only

          const double Ts = SAMPLING_FIRST_BUNCH + BUNCH_INTERVAL * bunch;
          const double Te = Ts + INT_RANGE;

          std::cout<<Ts<<" "<<Te<<std::endl;

          // Determine t_max from waveform peak within [Ts, Te]
          double t_max = 0.0;
          if (!FindTmaxInWaveform(wf, Ts, Te, t_max)) continue;

          // Define asymmetric search windows in the TDC-time domain.
          // NOTE: TDC time axis is reversed (as you pointed out), so we do NOT change ToT sign here.
          //
          // leading:  [-16*t_max+7334-50 , -16*Ts+7334] -> pick MIN
          // trailing: [-16*Te+7334       , -16*t_max+7334+50] -> pick MAX
          const double lead_a = -16.0 * t_max + 7334.0 - 50.0;
          const double lead_b = -16.0 * Ts    + 7334.0;
          const double lead_low  = std::min(lead_a, lead_b);
          const double lead_high = std::max(lead_a, lead_b);

          const double trail_a = -16.0 * Te    + 7334.0;
          const double trail_b = -16.0 * t_max + 7334.0 + 50.0;
          const double trail_low  = std::min(trail_a, trail_b);
          const double trail_high = std::max(trail_a, trail_b);

          double leadMin = 0.0, trailMax = 0.0;
          const bool hasLead  = MinInWindow(vLead,  lead_low,  lead_high,  leadMin);
          const bool hasTrail = MaxInWindow(vTrail, trail_low, trail_high, trailMax);

          // If leading/trailing not filled in that time window, ToT invalid
          if (!hasLead || !hasTrail) continue;

          // Keep the original sign convention: ToT = leading - trailing
          const double ToT = leadMin - trailMax;
          h2[ch]->Fill(ToT, ly);
        }
      }
    }

    nFilesUsed++;
    if (nFilesUsed % 10 == 0) {
      std::cout << "[INFO] Processed " << nFilesUsed << " files...\n";
    }
  }

  std::cout << "[INFO] Total processed files: " << nFilesUsed << "\n";
  if (nFilesUsed == 0) {
    std::cerr << "[ERROR] No matching ROOT files were processed.\n";
    return 4;
  }

  // Build profiles and fit functions once per channel, keep them alive
  std::vector<TProfile*> profs(NCH, nullptr);
  std::vector<TF1*>      fits (NCH, nullptr);
  // Graphs made from profile points after applying y-range cut
  std::vector<TGraphErrors*> grs(NCH, nullptr);

  // Use mode (most probable value) in each ToT bin to define representative black points
  const double yModeMin = 10.0;   // optional: accept mode only within this LY range
  const double yModeMax = 100.0;
  const int    minEntriesPerXbin = 30; // require enough stats in each ToT bin
  // Fit x-range (you can keep as you like)
  const double fitMinX_user = xMin;
  const double fitMaxX_user = xMax;

  // Update histogram titles using cablenum
  for (int ch = 0; ch < NCH; ++ch) {
    // Fallback if cablenum was never observed (should not happen)
    const int cnum = (cableLabel[ch] >= 0) ? cableLabel[ch] : ch;

    // Title format: "cable <cablenum>;x;y"
    // (Title in ROOT histogram includes axis titles separated by ';')
    TString newTitle;
    // x-axis: ToT, y-axis: lightyield (matches Fill(ToT, ly))
    newTitle.Form("cablenum %d;ToT (leading - trailing);light yield from ADC [p.e.]", cnum);
    h2[ch]->SetTitle(newTitle.Data());
  }

  // Create profiles (for visualization) and perform fits using MODE points
  for (int ch = 0; ch < NCH; ++ch) {
    std::string pname = "prof_ch" + std::to_string(ch);
    profs[ch] = h2[ch]->ProfileX(pname.c_str());   // owned by ROOT; keep pointer
    profs[ch]->SetDirectory(nullptr);              // detach from any file
    profs[ch]->SetTitle(h2[ch]->GetTitle());
    profs[ch]->GetXaxis()->SetTitle(h2[ch]->GetXaxis()->GetTitle());
    profs[ch]->GetYaxis()->SetTitle(h2[ch]->GetYaxis()->GetTitle());
    profs[ch]->SetStats(1);

    // --- Build TGraphErrors from MODE of Y in each X(ToT) bin ---
    std::vector<double> xs, ys, exs, eys;
    const int nxb = h2[ch]->GetNbinsX();
    xs.reserve(nxb); ys.reserve(nxb); exs.reserve(nxb); eys.reserve(nxb);

    for (int ibx = 1; ibx <= nxb; ++ibx) {
      // Optional: restrict x range used for mode extraction
      const double xCenter = h2[ch]->GetXaxis()->GetBinCenter(ibx);
      if (xCenter < fitMinX_user || xCenter > fitMaxX_user) continue;

      // Project Y for this single X bin
      std::string pyname = "py_ch" + std::to_string(ch) + "_x" + std::to_string(ibx);
      std::unique_ptr<TH1D> py(h2[ch]->ProjectionY(pyname.c_str(), ibx, ibx));
      const double nent = py->GetEntries();
      if (nent < minEntriesPerXbin) continue;

      const int maxBin = py->GetMaximumBin(); // mode bin
      const double yMode = py->GetXaxis()->GetBinCenter(maxBin);

      // Apply mode y-range cut (optional but usually helpful)
      if (yMode < yModeMin || yMode > yModeMax) continue;

      // Estimate y uncertainty: use RMS of the projected distribution
      // (If distribution is very non-Gaussian, consider using bin width as a floor.)
      double yerr = py->GetRMS();
      const double yBinW = py->GetXaxis()->GetBinWidth(maxBin);
      if (yerr < 0.5 * yBinW) yerr = 0.5 * yBinW;

      xs.push_back(xCenter);
      ys.push_back(yMode);
      exs.push_back(0.5 * h2[ch]->GetXaxis()->GetBinWidth(ibx));
      eys.push_back(yerr);
    }

    // Create graph (keep alive)
    std::string gname = "gr_ch" + std::to_string(ch);
    grs[ch] = new TGraphErrors((int)xs.size());
    grs[ch]->SetName(gname.c_str());
    for (int i = 0; i < (int)xs.size(); ++i) {
      grs[ch]->SetPoint(i, xs[i], ys[i]);
      grs[ch]->SetPointError(i, exs[i], eys[i]);
    }

    // Fit range in x: use user range, but if there are no points, keep default
    double fitMinX = fitMinX_user;
    double fitMaxX = fitMaxX_user;
    if (xs.empty()) {
      fitMinX = h2[ch]->GetXaxis()->GetXmin();
      fitMaxX = h2[ch]->GetXaxis()->GetXmax();
    }
    fitMinX=400;
    fitMaxX=1000;


    std::string fname = "f_ch" + std::to_string(ch);
    fits[ch] = MakeFitFunc(fname, fitMinX, fitMaxX);
    fits[ch]->SetParameters(1, 1, 1);

    // Do NOT fit here. Fit will be done in the drawing loop to ensure stats/fit box is created.

  }


  // Draw and save multi-page PDF: 8x4 pads, 272ch => 9 pages
  gStyle->SetOptStat(1110);
  gStyle->SetNumberContours(50);
  // Show fit results in stats box: entries/mean/rms + params + chi2/ndf
  gStyle->SetOptFit(1111);

  TCanvas c("c", "ToT vs LightYield 2D", 2400, 1400);
  c.Divide(8, 4, 0.001, 0.001);

  // Open multi-page PDF
  c.Print((outPdf + "[").c_str());

  const int perPage = 32;
  const int nPages = (NCH + perPage - 1) / perPage;

  for (int page = 0; page < nPages; ++page) {
    c.Clear();
    c.Divide(8, 4, 0.001, 0.001);

    for (int i = 0; i < perPage; ++i) {
      int ch = page * perPage + i;
      c.cd(i + 1);
      if (ch < NCH) {
        gPad->SetRightMargin(0.12);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.10);

        // Draw 2D first (no stats from TH2)
        h2[ch]->SetStats(0);
        h2[ch]->Draw("COLZ");

        // Draw profile first (stats box for profile if you want)
        profs[ch]->SetMarkerStyle(20);
        profs[ch]->SetMarkerSize(0.6);
        profs[ch]->SetLineWidth(2);
        profs[ch]->Draw("SAME");

        // Draw the mode points used for fit (black points) - TGraph is the fitted object
        grs[ch]->SetMarkerStyle(24);
        grs[ch]->SetMarkerSize(0.5);
        grs[ch]->Draw("P SAME");  // draw AFTER profile

        // IMPORTANT: Fit inside the drawing loop so that ROOT actually creates the fit box.
        // "R" uses TF1 range; "S" stores result; "0" is NOT used because we want drawing/box.
        grs[ch]->Fit(fits[ch], "QRSR");

        // Draw fit curve on top
        fits[ch]->SetLineWidth(3);
        fits[ch]->Draw("SAME");
        // Force stats box to appear for the profile and include fit info
        gPad->Modified();
        gPad->Update();
        // --- 1) Profile stats box (entries/mean/rms) ---
        if (auto* stp = dynamic_cast<TPaveStats*>(profs[ch]->FindObject("stats"))) {
          stp->SetTextSize(0.035);
          stp->SetX1NDC(0.12);
          stp->SetX2NDC(0.55);
          stp->SetY1NDC(0.72);
          stp->SetY2NDC(0.92);
        }

        // --- 2) Fit-result box belongs to TGraph, NOT to TProfile ---
        // For TGraph, stats box is typically stored in GetListOfFunctions().
        gPad->Modified();
        gPad->Update();
        if (auto* stg = dynamic_cast<TPaveStats*>(
              grs[ch]->GetListOfFunctions()->FindObject("stats"))) {
          stg->SetTextSize(0.035);
          // Place box (adjust if it overlaps colz palette)
          stg->SetX1NDC(0.58);
          stg->SetX2NDC(0.95);
          stg->SetY1NDC(0.72);
          stg->SetY2NDC(0.92);
        }
        std::cout << "ch=" << ch
          << " N=" << grs[ch]->GetN()
          << " funcs=" << grs[ch]->GetListOfFunctions()->GetSize()
          << " has_stats=" << (grs[ch]->GetListOfFunctions()->FindObject("stats")!=nullptr)
          << "\n";

      } else {
        // Empty pad
        gPad->Clear();
      }
    }

    c.Print(outPdf.c_str());
  }

  // Close multi-page PDF
  c.Print((outPdf + "]").c_str());

  std::cout << "[INFO] Saved: " << outPdf << "\n";
  return 0;
}
