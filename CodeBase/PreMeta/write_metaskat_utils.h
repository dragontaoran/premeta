// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Using provided inputs, Writes .MSSD and .MINFO files,
//              which mimic format of 1st-phase (summary statistics) output
//              of MetaSKAT (and can be subsequently used by MetaSKAT
//              for meta-analysis.

#ifndef WRITE_METASKAT_UTILS_H
#define WRITE_METASKAT_UTILS_H


#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "FileReaderUtils/vcf_utils.h"

#include <map>
#include <string>
#include <vector>

using file_reader_utils::Nucleotide;

namespace premeta {

bool WriteMssdFile(
    const int study_num, const double& rescale,
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<pair<Position, Position>, double>& covariances,
    const string& outfile);
// Same as above, with default empty 'snps_to_skip'.
inline bool WriteMssdFile(
    const int study_num, const double& rescale,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<pair<Position, Position>, double>& covariances,
    const string& outfile) {
  return WriteMssdFile(
      study_num, rescale, set<Position>(), monomorphic_snps,
      snp_to_ref_alt_and_study, gene_to_snp_positions, covariances, outfile);
}

bool WriteMInfoFile(
    const int study_num, const int num_samples, const int num_snps,
    const double& rescale, const double& sigma_sq,
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const string& outfile);
// Same as above, with default empty 'snps_to_skip'.
inline bool WriteMInfoFile(
    const int study_num, const int num_samples, const int num_snps,
    const double& rescale, const double& sigma_sq,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const string& outfile) {
  return WriteMInfoFile(
      study_num, num_samples, num_snps, rescale, sigma_sq, set<Position>(),
      monomorphic_snps, snp_to_ref_alt_and_study, gene_to_snp_positions,
      snp_info, outfile);
}

}  // namespace premeta
#endif
