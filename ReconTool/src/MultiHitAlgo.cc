#include "MultiHitAlgo.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "TH1D.h"

#include "FROSTConstants.h"

namespace FROST {

void DivideGroup1(const TH1D* h,
                  std::vector<std::vector<int>>& group,
                  const MultiHitAlgoConfig& cfg) {
  // Split bins above threshold1 into contiguous groups.
  // Single-bin groups may be discarded if they look like statistical fluctuation,

  group.clear();

  const int nb = h->GetNbinsX();
  const int nfib = (nb == kNfibX) ? kNfibX : ((nb == kNfibY) ? kNfibY : nb);

  std::vector<int> overthr;
  overthr.reserve(nfib);

  for (int l = 1; l <= nfib; ++l) {
    if (h->GetBinContent(l) > cfg.threshold1) {
      overthr.push_back(l);
    }
  }

  std::vector<std::vector<int>> groupall;
  int groupcounter = 0;

  if (!overthr.empty()) {
    groupall.push_back(std::vector<int>());
    groupall[0].push_back(overthr.at(0));

    for (std::size_t l = 1; l < overthr.size(); ++l) {
      if (overthr.at(l) != overthr.at(l - 1) + 1) {
        ++groupcounter;
        groupall.push_back(std::vector<int>());
      }
      groupall[groupcounter].push_back(overthr.at(l));
    }

    int counter = 0;
    for (int l = 0; l <= groupcounter; ++l) {
      if (groupall[l].size() == 1) {
        const int b = groupall[l][0];
        const double c = h->GetBinContent(b);
        const double s = std::sqrt(c);

        if (b == 1) {
          if (std::abs(h->GetBinContent(b + 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
          }
        } else if (h->GetNbinsX() == kNfibX && b == kNfibX) {
          if (std::abs(h->GetBinContent(b - 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
          }
        } else if (h->GetNbinsX() == kNfibY && b == kNfibY) {
          if (std::abs(h->GetBinContent(b - 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
          }
        } else {
          if (std::abs(h->GetBinContent(b - 1) - c) > s &&
              std::abs(h->GetBinContent(b + 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
          }
        }
      } else {
        group.push_back(std::vector<int>());
        for (std::size_t k = 0; k < groupall[l].size(); ++k) {
          group[counter].push_back(groupall[l][k]);
        }
        ++counter;
      }
    }
  }
}

void DivideGroup2(const TH1D* h,
                  std::vector<std::vector<int>>& group,
                  const MultiHitAlgoConfig& cfg) {
  // For each existing group, take up to 6 strongest bins and re-split them
  // using thr = min(threshold2, 0.25*max) (original logic).
  //

  if (group.empty()) return;

  std::vector<std::vector<int>> groupref = group;
  group.clear();

  const int nb = h->GetNbinsX();
  const int nfib = (nb == kNfibX) ? kNfibX : ((nb == kNfibY) ? kNfibY : nb);

  int counter = 0;

  for (std::size_t i = 0; i < groupref.size(); ++i) {
    // Build a temporary histogram that contains only the bins of this group.
    TH1D tmp("tmp", "tmp", nfib, -nfib * 5.0, nfib * 5.0);
    tmp.SetDirectory(nullptr);
    for (int j = 1; j <= nfib; ++j) tmp.SetBinContent(j, 0.0);
    for (std::size_t j = 0; j < groupref[i].size(); ++j) {
      const int b = groupref[i][j];
      if (1 <= b && b <= nfib) tmp.SetBinContent(b, h->GetBinContent(b));
    }

    int rank[6];
    double content[6];

    rank[0] = tmp.GetMaximumBin();
    content[0] = tmp.GetBinContent(rank[0]);
    tmp.SetBinContent(rank[0], 0);

    for (int r = 1; r < 6; ++r) {
      rank[r] = tmp.GetMaximumBin();
      content[r] = tmp.GetBinContent(rank[r]);
      if (content[r] <= 0) {
        rank[r] = 1000;
      } else {
        tmp.SetBinContent(rank[r], 0);
      }
    }

    std::sort(rank, rank + 6);

    double thr = 0.0;
    if (content[0] > cfg.threshold2 * 4.0) {
      thr = cfg.threshold2;
    } else {
      thr = content[0] * 0.25;
    }

    std::vector<int> overthr;
    overthr.reserve(6);
    for (int l = 0; l < 6; ++l) {
      if (rank[l] <= nfib) {
        if (h->GetBinContent(rank[l]) > thr) {
          overthr.push_back(rank[l]);
        }
      }
    }

    if (overthr.empty()) continue;

    std::vector<std::vector<int>> groupall;
    int groupcounter = 0;

    groupall.push_back(std::vector<int>());
    groupall[0].push_back(overthr.at(0));

    for (std::size_t l = 1; l < overthr.size(); ++l) {
      if (overthr.at(l) != overthr.at(l - 1) + 1) {
        ++groupcounter;
        groupall.push_back(std::vector<int>());
      }
      groupall[groupcounter].push_back(overthr.at(l));
    }

    int localcounter = 0;
    for (int l = 0; l <= groupcounter; ++l) {
      if (groupall[l].size() == 1) {
        const int b = groupall[l][0];
        const double c = h->GetBinContent(b);
        const double s = std::sqrt(c);

        if (b == 1) {
          if (std::abs(h->GetBinContent(b + 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
            ++localcounter;
          }
        } else if (h->GetNbinsX() == kNfibX && b == kNfibX) {
          if (std::abs(h->GetBinContent(b - 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
            ++localcounter;
          }
        } else if (h->GetNbinsX() == kNfibY && b == kNfibY) {
          if (std::abs(h->GetBinContent(b - 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
            ++localcounter;
          }
        } else {
          if (std::abs(h->GetBinContent(b - 1) - c) > s &&
              std::abs(h->GetBinContent(b + 1) - c) > s) {
            group.push_back(std::vector<int>());
            group[counter].push_back(b);
            ++counter;
            ++localcounter;
          }
        }
      } else {
        group.push_back(std::vector<int>());
        for (std::size_t k = 0; k < groupall[l].size(); ++k) {
          group[counter].push_back(groupall[l][k]);
        }
        ++counter;
        ++localcounter;
      }
    }

    // If everything was rejected as "fluctuation", keep the first bin anyway.
    if (localcounter == 0) {
      group.push_back(std::vector<int>());
      group[counter].push_back(groupall[0][0]);
      ++counter;
    }
  }
}

void MultiReconstructionX(const TH1D* h,
                          std::vector<double>& xgmulti,
                          std::vector<int>& rectype,
                          const std::vector<std::vector<int>>& groupx,
                          const MultiHitAlgoConfig& cfg) {
  // rectype semantics follow the original pasted version:
  //   - When groupx.size()==1 -> push two candidates with rectype {0,1}
  //   - When groupx.size()>1  -> push multiple candidates; candidates from split case use rectype {0,1},
  //                             candidates from merged case use rectype {2}
  //   - When groupx is empty  -> fall back to global xg, then push {0,1}

  xgmulti.clear();
  rectype.clear();

  const int nfibx = kNfibX;
  const double alpha = cfg.alpha;

  auto xpos = [&](int bin) -> double {
    return (-5.0 * nfibx - 5.0 + 10.0 * static_cast<double>(bin));
  };

  double light = 0.0;
  double xgmulti0 = 0.0;
  double xgmulti1 = 0.0;

  if (!groupx.empty()) {
    if (groupx.size() == 1) {
      const auto& g = groupx[0];

      if (g.size() == 1) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;
        for (int i = 1; i <= nfibx; ++i) {
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 = xgmulti0 / light;
        xgmulti1 = xgmulti0;

      } else if (g.size() == 2) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;

        for (int i = g.at(0); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        xgmulti0 = xgmulti0 / light;

        light = 0.0;
        xgmulti1 += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha) * xpos(g.at(0));
        light += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha);
        for (int i = g.at(1); i <= nfibx; ++i) {
          if (i == g.at(1) + 2) break;
          xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti1 = xgmulti1 / light;

      } else if (g.size() == 3) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;

        for (int i = g.at(0); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        xgmulti0 = xgmulti0 / light;

        light = 0.0;
        xgmulti1 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        for (int i = g.at(2); i <= nfibx; ++i) {
          if (i == g.at(2) + 2) break;
          xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti1 = xgmulti1 / light;

      } else if (g.size() == 4) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;

        for (int i = g.at(1); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 = xgmulti0 / light;

        light = 0.0;
        for (int i = g.at(2); i <= nfibx; ++i) {
          if (i == g.at(3) + 2) break;
          xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti1 = xgmulti1 / light;

      } else if (g.size() == 5) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;

        for (int i = g.at(1); i >= 1; --i) {
          if (i == g.at(0) - 3) break;
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * xpos(g.at(2));
        light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
        xgmulti0 = xgmulti0 / light;

        light = 0.0;
        for (int i = g.at(3); i <= nfibx; ++i) {
          if (i == g.at(4) + 3) break;
          xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti1 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * xpos(g.at(2));
        light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
        xgmulti1 = xgmulti1 / light;

      } else if (g.size() == 6) {
        light = 0.0;
        xgmulti0 = 0.0;
        xgmulti1 = 0.0;

        for (int i = g.at(2); i >= 1; --i) {
          if (i == g.at(0) - 3) break;
          xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti0 = xgmulti0 / light;

        light = 0.0;
        for (int i = g.at(3); i <= nfibx; ++i) {
          if (i == g.at(5) + 3) break;
          xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        xgmulti1 = xgmulti1 / light;
      }

      xgmulti.push_back(xgmulti0);
      rectype.push_back(0);
      xgmulti.push_back(xgmulti1);
      rectype.push_back(1);

    } else {
      std::vector<int> dgroup;
      std::vector<int> middle;

      for (std::size_t i = 0; i + 1 < groupx.size(); ++i) {
        dgroup.push_back(groupx[i + 1].at(0) - groupx[i].at(groupx[i].size() - 1) - 1);
        if (dgroup.back() % 2 == 0) middle.push_back(dgroup.back() / 2);
        else middle.push_back((dgroup.back() + 1) / 2);
      }

      for (std::size_t k = 0; k < groupx.size(); ++k) {
        int middleminus;
        int middleplus;
        int dgroupminus;
        int dgroupplus;

        if (k == 0) {
          middleminus = 1;
          middleplus = middle[k];
          dgroupminus = -1;
          dgroupplus = dgroup[k];
        } else if (k + 1 == groupx.size()) {
          middleminus = middle[k - 1];
          middleplus = nfibx;
          dgroupminus = dgroup[k - 1];
          dgroupplus = -1;
        } else {
          middleminus = middle[k - 1];
          middleplus = middle[k];
          dgroupminus = dgroup[k - 1];
          dgroupplus = dgroup[k];
        }

        TH1D tmp("tmpx", "tmpx", nfibx, -nfibx * 5.0, nfibx * 5.0);
        tmp.SetDirectory(nullptr);
        for (int j = 1; j <= nfibx; ++j) tmp.SetBinContent(j, 0.0);
        for (std::size_t j = 0; j < groupx[k].size(); ++j) {
          tmp.SetBinContent(groupx[k][j], h->GetBinContent(groupx[k][j]));
        }

        int maxbin = tmp.GetMaximumBin();
        tmp.SetBinContent(maxbin, 0);
        int secondbin;
        int thirdbin;

        if (tmp.GetBinContent(tmp.GetMaximumBin()) <= 0) {
          secondbin = -1;
        } else {
          secondbin = tmp.GetMaximumBin();
          tmp.SetBinContent(secondbin, 0);
        }

        if (tmp.GetBinContent(tmp.GetMaximumBin()) <= 0) {
          thirdbin = -1;
        } else {
          thirdbin = tmp.GetMaximumBin();
        }

        if (secondbin > 0 && thirdbin > 0 &&
            (std::abs(secondbin - maxbin) > 1 ||
             (std::abs(thirdbin - maxbin) > 1 && std::abs(thirdbin - secondbin) > 1))) {

          const auto& g = groupx[k];

          if (g.size() == 2) {
            light = 0.0; xgmulti0 = 0.0; xgmulti1 = 0.0;
            for (int i = g.at(0); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            xgmulti0 = xgmulti0 / light;

            light = 0.0;
            xgmulti1 += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha) * xpos(g.at(0));
            light += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha);
            for (int i = g.at(1); i <= nfibx; ++i) {
              if (i == g.at(1) + 2) break;
              xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti1 = xgmulti1 / light;

          } else if (g.size() == 3) {
            light = 0.0; xgmulti0 = 0.0; xgmulti1 = 0.0;
            for (int i = g.at(0); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            xgmulti0 = xgmulti0 / light;

            light = 0.0;
            xgmulti1 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * xpos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            for (int i = g.at(2); i <= nfibx; ++i) {
              if (i == g.at(2) + 2) break;
              xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti1 = xgmulti1 / light;

          } else if (g.size() == 4) {
            light = 0.0; xgmulti0 = 0.0; xgmulti1 = 0.0;
            for (int i = g.at(1); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti0 = xgmulti0 / light;

            light = 0.0;
            for (int i = g.at(2); i <= nfibx; ++i) {
              if (i == g.at(3) + 2) break;
              xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti1 = xgmulti1 / light;

          } else if (g.size() == 5) {
            light = 0.0; xgmulti0 = 0.0; xgmulti1 = 0.0;
            for (int i = g.at(1); i >= 1; --i) {
              if (i == g.at(0) - 3 || (dgroupminus > 0 && i == g.at(0) - dgroupminus - 1)) break;
              xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti0 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * xpos(g.at(2));
            light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
            xgmulti0 = xgmulti0 / light;

            light = 0.0;
            for (int i = g.at(3); i <= nfibx; ++i) {
              if (i == g.at(4) + 3 || (dgroupplus > 0 && i == g.at(4) + dgroupplus + 1)) break;
              xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti1 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * xpos(g.at(2));
            light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
            xgmulti1 = xgmulti1 / light;

          } else if (g.size() == 6) {
            light = 0.0; xgmulti0 = 0.0; xgmulti1 = 0.0;
            for (int i = g.at(2); i >= 1; --i) {
              if (i == g.at(0) - 3 || (dgroupminus > 0 && i == g.at(0) - dgroupminus - 1)) break;
              xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti0 = xgmulti0 / light;

            light = 0.0;
            for (int i = g.at(3); i <= nfibx; ++i) {
              if (i == g.at(5) + 3 || (dgroupplus > 0 && i == g.at(5) + dgroupplus + 1)) break;
              xgmulti1 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            xgmulti1 = xgmulti1 / light;
          }

          xgmulti.push_back(xgmulti0);
          rectype.push_back(0);
          xgmulti.push_back(xgmulti1);
          rectype.push_back(1);

        } else {
          light = 0.0;
          double xg = 0.0;

          const auto& g = groupx[k];
          for (int i = g.at(0); i <= g.at(static_cast<int>(g.size()) - 1); ++i) {
            xg += std::pow(h->GetBinContent(i), alpha) * xpos(i);
            light += std::pow(h->GetBinContent(i), alpha);
          }

          int j = 1;
          for (int i = g.at(static_cast<int>(g.size()) - 1) + 1;
               i <= g.at(static_cast<int>(g.size()) - 1) + middleplus; ++i) {
            if ((g.at(0) - j) < g.at(0) - middleminus) break;

            if (dgroupplus > 0 && dgroupplus % 2 != 0 &&
                i == g.at(static_cast<int>(g.size()) - 1) + middleplus) {
              xg += std::pow(h->GetBinContent(i) / 2.0, alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i) / 2.0, alpha);
            } else {
              xg += std::pow(h->GetBinContent(i), alpha) * xpos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }

            if (dgroupminus > 0 && dgroupminus % 2 != 0 &&
                g.at(0) - j == g.at(0) - middleminus) {
              xg += std::pow(h->GetBinContent(g.at(0) - j) / 2.0, alpha) * xpos(g.at(0) - j);
              light += std::pow(h->GetBinContent(g.at(0) - j) / 2.0, alpha);
              ++j;
            } else {
              xg += std::pow(h->GetBinContent(g.at(0) - j), alpha) * xpos(g.at(0) - j);
              light += std::pow(h->GetBinContent(g.at(0) - j), alpha);
              ++j;
            }
          }

          xg = xg / light;
          xgmulti.push_back(xg);
          rectype.push_back(2);
        }
      }
    }
  } else {
    light = 0.0;
    xgmulti0 = 0.0;
    xgmulti1 = 0.0;
    for (int i = 1; i <= nfibx; ++i) {
      xgmulti0 += std::pow(h->GetBinContent(i), alpha) * xpos(i);
      light += std::pow(h->GetBinContent(i), alpha);
    }
    xgmulti0 = xgmulti0 / light;
    xgmulti1 = xgmulti0;

    xgmulti.push_back(xgmulti0);
    rectype.push_back(0);
    xgmulti.push_back(xgmulti1);
    rectype.push_back(1);
  }
}

void MultiReconstructionY(const TH1D* h,
                          std::vector<double>& ygmulti,
                          std::vector<int>& rectype,
                          const std::vector<std::vector<int>>& groupy,
                          const MultiHitAlgoConfig& cfg) {

  ygmulti.clear();
  rectype.clear();

  const int nfiby = kNfibY;
  const double alpha = cfg.alpha;

  auto ypos = [&](int bin) -> double {
    return (-5.0 * nfiby - 5.0 + 10.0 * static_cast<double>(bin));
  };

  double light = 0.0;
  double ygmulti0 = 0.0;
  double ygmulti1 = 0.0;

  if (!groupy.empty()) {
    if (groupy.size() == 1) {
      const auto& g = groupy[0];

      if (g.size() == 1) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;
        for (int i = 1; i <= nfiby; ++i) {
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 = ygmulti0 / light;
        ygmulti1 = ygmulti0;

      } else if (g.size() == 2) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;

        for (int i = g.at(0); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        ygmulti0 = ygmulti0 / light;

        light = 0.0;
        ygmulti1 += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha) * ypos(g.at(0));
        light += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha);
        for (int i = g.at(1); i <= nfiby; ++i) {
          if (i == g.at(1) + 2) break;
          ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti1 = ygmulti1 / light;

      } else if (g.size() == 3) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;

        for (int i = g.at(0); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        ygmulti0 = ygmulti0 / light;

        light = 0.0;
        ygmulti1 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
        light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
        for (int i = g.at(2); i <= nfiby; ++i) {
          if (i == g.at(2) + 2) break;
          ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti1 = ygmulti1 / light;

      } else if (g.size() == 4) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;

        for (int i = g.at(1); i >= 1; --i) {
          if (i == g.at(0) - 2) break;
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 = ygmulti0 / light;

        light = 0.0;
        for (int i = g.at(2); i <= nfiby; ++i) {
          if (i == g.at(3) + 2) break;
          ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti1 = ygmulti1 / light;

      } else if (g.size() == 5) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;

        for (int i = g.at(1); i >= 1; --i) {
          if (i == g.at(0) - 3) break;
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * ypos(g.at(2));
        light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
        ygmulti0 = ygmulti0 / light;

        light = 0.0;
        for (int i = g.at(3); i <= nfiby; ++i) {
          if (i == g.at(4) + 3) break;
          ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti1 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * ypos(g.at(2));
        light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
        ygmulti1 = ygmulti1 / light;

      } else if (g.size() == 6) {
        light = 0.0;
        ygmulti0 = 0.0;
        ygmulti1 = 0.0;

        for (int i = g.at(2); i >= 1; --i) {
          if (i == g.at(0) - 3) break;
          ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti0 = ygmulti0 / light;

        light = 0.0;
        for (int i = g.at(3); i <= nfiby; ++i) {
          if (i == g.at(5) + 3) break;
          ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
          light += std::pow(h->GetBinContent(i), alpha);
        }
        ygmulti1 = ygmulti1 / light;
      }

      ygmulti.push_back(ygmulti0);
      rectype.push_back(0);
      ygmulti.push_back(ygmulti1);
      rectype.push_back(1);

    } else {
      std::vector<int> dgroup;
      std::vector<int> middle;

      for (std::size_t i = 0; i + 1 < groupy.size(); ++i) {
        dgroup.push_back(groupy[i + 1].at(0) - groupy[i].at(groupy[i].size() - 1) - 1);
        if (dgroup.back() % 2 == 0) middle.push_back(dgroup.back() / 2);
        else middle.push_back((dgroup.back() + 1) / 2);
      }

      for (std::size_t k = 0; k < groupy.size(); ++k) {
        int middleminus;
        int middleplus;
        int dgroupminus;
        int dgroupplus;

        if (k == 0) {
          middleminus = 1;
          middleplus = middle[k];
          dgroupminus = -1;
          dgroupplus = dgroup[k];
        } else if (k + 1 == groupy.size()) {
          middleminus = middle[k - 1];
          middleplus = nfiby;
          dgroupminus = dgroup[k - 1];
          dgroupplus = -1;
        } else {
          middleminus = middle[k - 1];
          middleplus = middle[k];
          dgroupminus = dgroup[k - 1];
          dgroupplus = dgroup[k];
        }

        TH1D tmp("tmpy", "tmpy", nfiby, -nfiby * 5.0, nfiby * 5.0);
        tmp.SetDirectory(nullptr);
        for (int j = 1; j <= nfiby; ++j) tmp.SetBinContent(j, 0.0);
        for (std::size_t j = 0; j < groupy[k].size(); ++j) {
          tmp.SetBinContent(groupy[k][j], h->GetBinContent(groupy[k][j]));
        }

        int maxbin = tmp.GetMaximumBin();
        tmp.SetBinContent(maxbin, 0);
        int secondbin;
        int thirdbin;

        if (tmp.GetBinContent(tmp.GetMaximumBin()) <= 0) {
          secondbin = -1;
        } else {
          secondbin = tmp.GetMaximumBin();
          tmp.SetBinContent(secondbin, 0);
        }

        if (tmp.GetBinContent(tmp.GetMaximumBin()) <= 0) {
          thirdbin = -1;
        } else {
          thirdbin = tmp.GetMaximumBin();
        }

        if (secondbin > 0 && thirdbin > 0 &&
            (std::abs(secondbin - maxbin) > 1 ||
             (std::abs(thirdbin - maxbin) > 1 && std::abs(thirdbin - secondbin) > 1))) {

          const auto& g = groupy[k];

          if (g.size() == 2) {
            light = 0.0; ygmulti0 = 0.0; ygmulti1 = 0.0;
            for (int i = g.at(0); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            ygmulti0 = ygmulti0 / light;

            light = 0.0;
            ygmulti1 += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha) * ypos(g.at(0));
            light += std::pow(h->GetBinContent(g.at(0)) / 2.0, alpha);
            for (int i = g.at(1); i <= nfiby; ++i) {
              if (i == g.at(1) + 2) break;
              ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti1 = ygmulti1 / light;

          } else if (g.size() == 3) {
            light = 0.0; ygmulti0 = 0.0; ygmulti1 = 0.0;
            for (int i = g.at(0); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti0 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            ygmulti0 = ygmulti0 / light;

            light = 0.0;
            ygmulti1 += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha) * ypos(g.at(1));
            light += std::pow(h->GetBinContent(g.at(1)) / 2.0, alpha);
            for (int i = g.at(2); i <= nfiby; ++i) {
              if (i == g.at(2) + 2) break;
              ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti1 = ygmulti1 / light;

          } else if (g.size() == 4) {
            light = 0.0; ygmulti0 = 0.0; ygmulti1 = 0.0;
            for (int i = g.at(1); i >= 1; --i) {
              if (i == g.at(0) - 2) break;
              ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti0 = ygmulti0 / light;

            light = 0.0;
            for (int i = g.at(2); i <= nfiby; ++i) {
              if (i == g.at(3) + 2) break;
              ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti1 = ygmulti1 / light;

          } else if (g.size() == 5) {
            light = 0.0; ygmulti0 = 0.0; ygmulti1 = 0.0;
            for (int i = g.at(1); i >= 1; --i) {
              if (i == g.at(0) - 3 || (dgroupminus > 0 && i == g.at(0) - dgroupminus - 1)) break;
              ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti0 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * ypos(g.at(2));
            light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
            ygmulti0 = ygmulti0 / light;

            light = 0.0;
            for (int i = g.at(3); i <= nfiby; ++i) {
              if (i == g.at(4) + 3 || (dgroupplus > 0 && i == g.at(4) + dgroupplus + 1)) break;
              ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti1 += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha) * ypos(g.at(2));
            light += std::pow(h->GetBinContent(g.at(2)) / 2.0, alpha);
            ygmulti1 = ygmulti1 / light;

          } else if (g.size() == 6) {
            light = 0.0; ygmulti0 = 0.0; ygmulti1 = 0.0;
            for (int i = g.at(2); i >= 1; --i) {
              if (i == g.at(0) - 3 || (dgroupminus > 0 && i == g.at(0) - dgroupminus - 1)) break;
              ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti0 = ygmulti0 / light;

            light = 0.0;
            for (int i = g.at(3); i <= nfiby; ++i) {
              if (i == g.at(5) + 3 || (dgroupplus > 0 && i == g.at(5) + dgroupplus + 1)) break;
              ygmulti1 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }
            ygmulti1 = ygmulti1 / light;
          }

          ygmulti.push_back(ygmulti0);
          rectype.push_back(0);
          ygmulti.push_back(ygmulti1);
          rectype.push_back(1);

        } else {
          light = 0.0;
          double yg = 0.0;

          const auto& g = groupy[k];
          for (int i = g.at(0); i <= g.at(static_cast<int>(g.size()) - 1); ++i) {
            yg += std::pow(h->GetBinContent(i), alpha) * ypos(i);
            light += std::pow(h->GetBinContent(i), alpha);
          }

          int j = 1;
          for (int i = g.at(static_cast<int>(g.size()) - 1) + 1;
               i <= g.at(static_cast<int>(g.size()) - 1) + middleplus; ++i) {
            if ((g.at(0) - j) < g.at(0) - middleminus) break;

            if (dgroupplus > 0 && dgroupplus % 2 != 0 &&
                i == g.at(static_cast<int>(g.size()) - 1) + middleplus) {
              yg += std::pow(h->GetBinContent(i) / 2.0, alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i) / 2.0, alpha);
            } else {
              yg += std::pow(h->GetBinContent(i), alpha) * ypos(i);
              light += std::pow(h->GetBinContent(i), alpha);
            }

            if (dgroupminus > 0 && dgroupminus % 2 != 0 &&
                g.at(0) - j == g.at(0) - middleminus) {
              yg += std::pow(h->GetBinContent(g.at(0) - j) / 2.0, alpha) * ypos(g.at(0) - j);
              light += std::pow(h->GetBinContent(g.at(0) - j) / 2.0, alpha);
              ++j;
            } else {
              yg += std::pow(h->GetBinContent(g.at(0) - j), alpha) * ypos(g.at(0) - j);
              light += std::pow(h->GetBinContent(g.at(0) - j), alpha);
              ++j;
            }
          }

          yg = yg / light;
          ygmulti.push_back(yg);
          rectype.push_back(2);
        }
      }
    }
  } else {
    light = 0.0;
    ygmulti0 = 0.0;
    ygmulti1 = 0.0;
    for (int i = 1; i <= nfiby; ++i) {
      ygmulti0 += std::pow(h->GetBinContent(i), alpha) * ypos(i);
      light += std::pow(h->GetBinContent(i), alpha);
    }
    ygmulti0 = ygmulti0 / light;
    ygmulti1 = ygmulti0;

    ygmulti.push_back(ygmulti0);
    rectype.push_back(0);
    ygmulti.push_back(ygmulti1);
    rectype.push_back(1);
  }
}
}  // namespace FROST
