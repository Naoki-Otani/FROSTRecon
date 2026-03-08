// File: DrawLightyieldvsADC.C
// Usage:
//   root -l -q DrawLightyieldvsADC.C
// or
//   root -l
//   .x DrawLightyieldvsADC.C

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TString.h>
#include <TError.h>
#include <TH1F.h>
#include <TLatex.h>

#include <fstream>
#include <sstream>
#include <string>

void DrawSaturatedLightyield() {
  // Input file path
  const char* inFile = "/home/nu/notani/FROSTRecon/ProcessMC/ADCsaturation/data/LightyieldvsADC.txt";

  // Disable the statistics box
  gStyle->SetOptStat(0);

  // Create a graph with y-errors
  auto gr = new TGraphErrors();
  gr->SetName("grSaturatedLightyield");
  gr->SetTitle("True light yield vs Saturated light yield;True light yield [p.e.];Saturated light yield [p.e.]");

  // Open the text file
  std::ifstream fin(inFile);
  if (!fin.is_open()) {
    Error("DrawSaturatedLightyield", "Cannot open input file: %s", inFile);
    return;
  }

  // Read line-by-line and skip comments/empty lines
  std::string line;
  int ip = 0;

  while (std::getline(fin, line)) {
    // Skip empty lines
    if (line.empty()) continue;

    // Skip comment lines (starting with '#')
    // Also skip lines that begin with spaces then '#'
    {
      std::string trimmed = line;
      trimmed.erase(0, trimmed.find_first_not_of(" \t"));
      if (trimmed.empty()) continue;
      if (trimmed[0] == '#') continue;
    }

    // Parse: amp adc adc_err
    double amp = 0.0, adc = 0.0, adc_err = 0.0;
    std::istringstream iss(line);
    if (!(iss >> amp >> adc >> adc_err)) {
      // If a line cannot be parsed, just ignore it
      continue;
    }

    gr->SetPoint(ip, amp*65./250., adc*0.016055);
    gr->SetPointError(ip, 0., adc_err*0.016055); // xerr = 0, yerr = adc_err
    ++ip;
  }

  fin.close();

  if (gr->GetN() <= 0) {
    Error("DrawSaturatedLightyield", "No valid data points were read from: %s", inFile);
    return;
  }

  // Marker/line style
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  gr->SetLineWidth(2);

  // Draw
  auto c1 = new TCanvas("c1", "True light yield vs Saturated light yield", 900, 700);
  c1->SetGrid();
  c1->SetLeftMargin(0.15);

  // Create an empty frame to define axis ranges
  TH1F* frame = c1->DrawFrame(0.0, 0.0, 2000.0, 400.0); // xmin,ymin,xmax,ymax
  frame->SetTitle("True light yield vs Saturated light yield;True light yield [p.e.];Saturated light yield [p.e.]");

  gr->Draw("P same"); // Axes + Points

  // Fit with y = a*x in the range x=[0,80]
  TF1* f = new TF1("f", "[0]+[1]*(1-exp(-pow((x/[2]),[3])))",80,1200);
  f->SetParameter(0, -3981); // Initial guess (adjust if needed)
  f->SetParameter(1, 6500);
  f->SetParameter(2, 1015);
  f->SetParameter(3, 0.02);

  // "R" uses the function range, "Q" is quiet, remove "Q" if you want output
  gr->Fit(f, "RQ");

  f->SetRange(80,2000);
  // (Optional) Draw the fit line explicitly (Fit already draws by default unless you use "0")
  f->Draw("same");

  gStyle->SetOptFit(1111); // Show fit parameters on the plot

  // Move the stats/fit box (NDC coordinates: 0..1)
  gPad->Update(); // Needed to create the stats box

  auto st = (TPaveStats*)gr->GetListOfFunctions()->FindObject("stats");
  if (st) {
    st->SetX1NDC(0.6); // left
    st->SetX2NDC(0.9); // right
    st->SetY1NDC(0.7); // bottom
    st->SetY2NDC(0.9); // top
  }
  gPad->Modified();
  gPad->Update();

  // Draw the fit function formula in the top-left using LaTeX
  TLatex latex;
  latex.SetTextColor(kRed);
  latex.SetNDC(true);        // Use normalized coordinates (0..1)
  latex.SetTextSize(0.035);  // Adjust text size as you like
  latex.SetTextAlign(13);    // Left-top alignment

  // LaTeX string for: A(x) = p0 + p1 (1 - exp(-(x/p2)^p3))
  latex.DrawLatex(0.18, 0.85,
    "Fit:  y=p_{0}+p_{1}\\left(1#minus e^{#minus \\left(x/p_{2}\\right)^{p_{3}}}\\right)"
  );

  double xref = 80.0;
  double yref = 80.0;

  double f80 = f->Eval(xref);
  cout<<f80<<endl;
  // Save outputs
  c1->SaveAs("/home/nu/notani/FROSTRecon/ProcessMC/ADCsaturation/plot/SaturatedLightyield.pdf");
}
