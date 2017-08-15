// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Using provided inputs, Writes Score and Covariance files,
//              which mimic format of 1st-phase (summary statistics) output
//              of RareMetalWorker (and can be subsequently used by RareMetal
//              for meta-analysis.

#ifndef WRITE_RAREMETAL_UTILS_H
#define WRITE_RAREMETAL_UTILS_H


#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "FileReaderUtils/vcf_utils.h"

#include <set>
#include <string>
#include <vector>

using file_reader_utils::Nucleotide;

namespace premeta {

bool PrintToRareMetal(
    const int study_num, const int num_samples, const double& rescale,
    const double& cov_rescale, const double& u_stat_rescale,
    const set<Position>& monomorphic_snps,
    const set<Position>& snps_to_skip,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, Position>& windows,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    const string& out_score_file, const string& out_cov_file);
// Same as above, but uses default empty set for snps_to_skip.
inline bool PrintToRareMetal(
    const int study_num, const int num_samples, const double& rescale,
    const double& cov_rescale, const double& u_stat_rescale,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, Position>& windows,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    const string& out_score_file, const string& out_cov_file) {
  return PrintToRareMetal(
      study_num, num_samples, rescale, cov_rescale, u_stat_rescale,
      monomorphic_snps, set<Position>(), snp_to_ref_alt_and_study, windows, snp_info,
      covariances, out_score_file, out_cov_file);
}
// Same as above, but for use-case: RareMetal -> RareMetal (so don't
// have to recompute windows).  
bool PrintToRareMetal(
    const int study_num, const int num_samples,
    const double& rescale, const double& u_stat_rescale,
    const set<Position>& monomorphic_snps,
    const set<Position>& snps_without_score_info,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    const string& out_score_file, const string& out_cov_file);

// Prints Score and Covariance files as copies as the input files, except
// Monomorphic and Multi-Allelic SNPs have been removed.
bool PrintRareMetaliWithMonoAndMultiSnpsRemoved(
    const string& in_score_file, const string& in_cov_file,
    const string& out_score_file, const string& out_cov_file,
    const set<Position>& non_std_chr_snps,
    const set<Position>& snps_to_exclude,
    const map<Position, int>& snps_in_score_file,
    const map<Position, int>& snps_in_cov_file,
    set<Position>* non_matching);

bool PrintCovariance(
    const int study_num, const int num_samples, const double& rescale,
    const Position& snp_pos,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    ofstream& cov_outfile, double* self_covariance);

// Appends score_file to score_files and cov_file to cov_files.
bool AddFilesToList(
    const string& score_files, const string& cov_files,
    const string& score_file, const string& cov_file);

}  // namespace premeta
#endif
