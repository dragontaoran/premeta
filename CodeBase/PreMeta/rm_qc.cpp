// Date: April 2016
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "rm_qc.h"

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_raremetal_utils.h"
#include "PreMeta/write_raremetal_utils.h"
#include "StringUtils/string_utils.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using file_reader_utils::Chromosome;
using file_reader_utils::Nucleotide;
using namespace string_utils;
using namespace std;

namespace premeta {

bool RareMetalQc(
    const bool is_new_version_from, const double& rescale,
    const bool remove_monomorphic, const bool remove_multi_allelic,
    const bool remove_multi_allelic_two, const bool remove_non_snps,
    const bool strict_matching,
    const string& in_score_file, const string& in_cov_file,
    const string& out_score_file, const string& out_cov_file,
    set<Position>* monomorphic_snps, set<Position>* multi_allelic_snps,
    set<Position>* multi_allelic_snps_two, set<Position>* non_snps,
    set<Position>* non_matching, set<Position>* non_std_chr_snps) {
  // Read Score File, looking for monomorphic and multi-allelic SNPs.
  map<Position, int> snps_in_score_file;
  if (!GetMonoAndMultiSnpsFromRareMetalScoreFile(
          in_score_file,
          strict_matching ? &snps_in_score_file : nullptr,
          monomorphic_snps, multi_allelic_snps, multi_allelic_snps_two,
          non_snps, non_std_chr_snps)) {
    cout << "ERROR: Unable to read score file '" << in_score_file
         << "'. Aborting." << endl;
    return false;
  }

  // Read Covariance File, looking for multi-allelic SNPs.
  map<Position, int> snps_in_cov_file;
  if (!GetMultiAllelicSnpsFromRareMetalCovarianceFile(
          in_cov_file,
          strict_matching ? &snps_in_cov_file : nullptr,
          multi_allelic_snps, non_std_chr_snps)) {
    cout << "ERROR: Unable to read covariance file '" << in_score_file
         << "'. Aborting." << endl;
    return false;
  }

  // Merge all SNPs to exclude in a single set.
  set<Position> snps_to_exclude;
  if (remove_monomorphic) {
    snps_to_exclude.insert(monomorphic_snps->begin(), monomorphic_snps->end());
  }
  if (remove_multi_allelic) {
    snps_to_exclude.insert(multi_allelic_snps->begin(), multi_allelic_snps->end());
  }
  if (remove_multi_allelic_two) {
    snps_to_exclude.insert(multi_allelic_snps_two->begin(), multi_allelic_snps_two->end());
  }
  if (remove_non_snps) {
    snps_to_exclude.insert(non_snps->begin(), non_snps->end());
  }
  if (non_std_chr_snps != nullptr) {
    snps_to_exclude.insert(non_std_chr_snps->begin(), non_std_chr_snps->end());
  }

  // Return, if no Monomorphic nor Multi-allelic nor Wrong Chromosome SNPs
  // were encountered.
  // UPDATE: Ran wanted to print out copies of .score and .cov files, even
  // if no SNPs were removed; so we skip this early abort.
  // if (snps_to_exclude.empty()) return true;

  // Write score and covariance file, by copying the originals but removing
  // monomorphic and multi-allelic SNPs.
  return PrintRareMetaliWithMonoAndMultiSnpsRemoved(
      in_score_file, in_cov_file, out_score_file, out_cov_file,
      *non_std_chr_snps, snps_to_exclude, snps_in_score_file, snps_in_cov_file,
      strict_matching ? non_matching : nullptr);
}

}  // namespace premeta
