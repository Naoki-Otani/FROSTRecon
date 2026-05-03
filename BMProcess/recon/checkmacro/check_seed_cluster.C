// check_seed_cluster.C

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

void print_set(const std::set<int>& s)
{
  std::cout << "{ ";
  for (int x : s) std::cout << x << " ";
  std::cout << "}";
}

int max_plane_gap(const std::set<int>& planes)
{
  if (planes.size() < 2) return 0;

  int max_gap = 0;
  auto prev = planes.begin();
  auto it = std::next(prev);

  for (; it != planes.end(); ++it, ++prev) {
    int gap = *it - *prev;
    if (gap > max_gap) max_gap = gap;
  }

  return max_gap;
}

void print_missing_planes(const std::set<int>& planes)
{
  if (planes.size() < 2) return;

  auto prev = planes.begin();
  auto it = std::next(prev);

  for (; it != planes.end(); ++it, ++prev) {
    int gap = *it - *prev;

    if (gap >= 3) {
      std::cout << "      gap from plane " << *prev
                << " to " << *it
                << " : missing planes = ";

      for (int p = *prev + 1; p < *it; p++) {
        std::cout << p << " ";
      }

      std::cout << "\n";
    }
  }
}

void check_seed_cluster(
    const char* filename="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-BMPMrecon/track_merged/BMPM_track_2026-02-15_20-40-12_Run0.root",
    // const char* filename="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-BMPMrecon/vertex/BMPM_vertex_2026-02-15_20-40-12_Run0_0.root",
    Long64_t max_entries = 10,
    Long64_t start_entry = 0,
    int target_detector = 5,
    int min_gap_to_print = 3)
{
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  TTree* tree = (TTree*)f->Get("tree");
  if (!tree) {
    std::cerr << "Cannot find tree named 'tree'" << std::endl;
    return;
  }

  TTreeReader reader(tree);

  // Hit-level branches
  TTreeReaderArray<UInt_t> hit_id(reader, "hits_.hit_id_");
  TTreeReaderArray<UInt_t> hit_parent_cluster_id(reader, "hits_.parent_cluster_id_");
  TTreeReaderArray<Int_t> hit_detector(reader, "hits_.detector_");
  TTreeReaderArray<Int_t> hit_view(reader, "hits_.view_");
  TTreeReaderArray<Int_t> hit_plane(reader, "hits_.plane_");

  // Recon cluster-level branches
  TTreeReaderValue<Int_t> num_recon_clusters(reader, "num_recon_clusters_");
  TTreeReaderArray<UInt_t> cluster_id(reader, "recon_clusters_.cluster_id_");
  TTreeReaderArray<Int_t> cluster_view(reader, "recon_clusters_.view_");
  TTreeReaderArray<Int_t> cluster_num_hits(reader, "recon_clusters_.num_hits_");

  const Long64_t nentries = tree->GetEntries();

  Long64_t end_entry = nentries;
  if (max_entries >= 0) {
    end_entry = std::min(nentries, start_entry + max_entries);
  }

  reader.SetEntriesRange(start_entry, end_entry);

  Long64_t entry = start_entry;
  Long64_t n_suspicious_entries = 0;

  while (reader.Next()) {
    // ------------------------------------------------------------
    // 1. Make BabyMIND input hit plane list by view
    // ------------------------------------------------------------
    std::map<int, std::set<int>> input_planes_by_view;
    std::map<int, int> detector_counts;
    std::map<UInt_t, int> bm_parent_cluster_id_counts;

    for (int ihit = 0; ihit < hit_id.GetSize(); ihit++) {
      detector_counts[hit_detector[ihit]]++;

      if (hit_detector[ihit] != target_detector) continue;

      input_planes_by_view[hit_view[ihit]].insert(hit_plane[ihit]);
      bm_parent_cluster_id_counts[hit_parent_cluster_id[ihit]]++;
    }

    // ------------------------------------------------------------
    // 2. Check whether BabyMIND input hits have 2-plane gap
    //    gap >= 3 means:
    //      e.g. plane 12 and 15 exist, but 13 and 14 are missing.
    // ------------------------------------------------------------
    bool suspicious = false;

    for (const auto& item : input_planes_by_view) {
      const std::set<int>& planes = item.second;
      int gap = max_plane_gap(planes);

      if (gap >= min_gap_to_print) {
        suspicious = true;
      }
    }

    if (!suspicious) {
      entry++;
      continue;
    }

    n_suspicious_entries++;

    // ------------------------------------------------------------
    // 3. Print suspicious entry
    // ------------------------------------------------------------
    std::cout << "\n========================================\n";
    std::cout << "Entry = " << entry << "\n";
    std::cout << "Number of all hits = " << hit_id.GetSize() << "\n";
    std::cout << "num_recon_clusters_ = " << *num_recon_clusters << "\n";
    std::cout << "recon_clusters_.cluster_id_.GetSize() = "
              << cluster_id.GetSize() << "\n";

    std::cout << "\n[Detector counts]\n";
    for (const auto& item : detector_counts) {
      std::cout << "  detector = " << item.first
                << " nhits = " << item.second << "\n";
    }

    std::cout << "\n[Input BabyMIND planes by view]\n";
    for (const auto& item : input_planes_by_view) {
      int view = item.first;
      const std::set<int>& planes = item.second;

      int gap = max_plane_gap(planes);

      std::cout << "  view = " << view << " planes = ";
      print_set(planes);
      std::cout << "  max_gap = " << gap << "\n";

      if (gap >= min_gap_to_print) {
        print_missing_planes(planes);
      }
    }

    // ------------------------------------------------------------
    // 4. Print parent_cluster_id values of BabyMIND hits
    //    This tells whether hits_.parent_cluster_id_ can be used.
    // ------------------------------------------------------------
    std::cout << "\n[BabyMIND hit parent_cluster_id counts]\n";
    for (const auto& item : bm_parent_cluster_id_counts) {
      std::cout << "  parent_cluster_id = " << item.first
                << " nhits = " << item.second << "\n";
    }

    // ------------------------------------------------------------
    // 5. Print recon cluster summary
    //    This does not tell cluster hit planes unless parent_cluster_id works.
    // ------------------------------------------------------------
    std::cout << "\n[ReconCluster summary]\n";

    for (int icl = 0; icl < cluster_id.GetSize(); icl++) {
      std::cout << "  cluster index = " << icl
                << "  cluster_id = " << cluster_id[icl]
                << "  view = " << cluster_view[icl]
                << "  num_hits = " << cluster_num_hits[icl]
                << "\n";
    }

    // ------------------------------------------------------------
    // 6. Try to match hit parent_cluster_id to cluster_id
    //    If this prints empty plane lists, parent_cluster_id is not usable.
    // ------------------------------------------------------------
    std::cout << "\n[Try matching BabyMIND hits by parent_cluster_id]\n";

    for (int icl = 0; icl < cluster_id.GetSize(); icl++) {
      UInt_t cid = cluster_id[icl];

      std::map<int, std::set<int>> matched_planes_by_view;
      std::vector<UInt_t> matched_hit_ids;

      for (int ihit = 0; ihit < hit_id.GetSize(); ihit++) {
        if (hit_detector[ihit] != target_detector) continue;
        if (hit_parent_cluster_id[ihit] != cid) continue;

        matched_planes_by_view[hit_view[ihit]].insert(hit_plane[ihit]);
        matched_hit_ids.push_back(hit_id[ihit]);
      }

      std::cout << "  cluster_id = " << cid
                << "  cluster_view = " << cluster_view[icl]
                << "  cluster_num_hits = " << cluster_num_hits[icl]
                << "\n";

      if (matched_planes_by_view.empty()) {
        std::cout << "    No BabyMIND hits matched by parent_cluster_id\n";
      } else {
        for (const auto& item : matched_planes_by_view) {
          std::cout << "    matched view = " << item.first
                    << " planes = ";
          print_set(item.second);
          std::cout << "\n";
        }

        std::cout << "    matched hit ids = ";
        for (UInt_t hid : matched_hit_ids) {
          std::cout << hid << " ";
        }
        std::cout << "\n";
      }
    }

    entry++;
  }

  std::cout << "\n========================================\n";
  std::cout << "Finished scan\n";
  std::cout << "Scanned entries: " << start_entry << " - " << end_entry - 1 << "\n";
  std::cout << "Suspicious entries with max plane gap >= "
            << min_gap_to_print << " : "
            << n_suspicious_entries << "\n";

  f->Close();
}
