// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Reads 1st-phase (summary statistics) output of RAREMETAL's
// ScoreSeq function, populating generic data structures.

#ifndef READ_RAREMETAL_UTILS_H
#define READ_RAREMETAL_UTILS_H

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <set>
#include <string>
#include <vector>

namespace premeta {

// Reads RareMetal's score file.
bool GetInfoFromRareMetalScoreFile(
    const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& file_info,
    int* num_samples, int* num_snps,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<Position, SnpInfo>* scores);

// Reads RareMetal's score file, identifying Monomorphic and Multi-Allelic Snps.
bool GetMonoAndMultiSnpsFromRareMetalScoreFile(
    const string& score_file,
    map<Position, int>* snps_in_score_file,
    set<Position>* monomorphic_snps,
    set<Position>* multi_allelic_snps,
    set<Position>* multi_allelic_snps_two,
    set<Position>* non_snps,
    set<Position>* non_std_chr_snps);

// Reads RareMetal's Gene Grouping file.
bool GetInfoFromRareMetalGroupFile(
    const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const FileInfo& file_info,
    set<Position>* snps_without_score_info,
    map<string, vector<Position>>* gene_to_snp_positions);

// Reads RareMetal's Covariance file.
bool GetInfoFromRareMetalCovarianceFile(
    const bool is_new_version, const bool keep_snps_without_score_info,
    const int study_num, const int num_samples, const FileInfo& file_info,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    set<Position>* snps_without_score_info,
    map<Position, map<Position, double>>* snp_to_cov_w_neighbors);

// Reads RareMetal's covariance file, identifying Multi-Allelic Snps.
bool GetMultiAllelicSnpsFromRareMetalCovarianceFile(
    const string& cov_file,
    map<Position, int>* snps_in_cov_file,
    set<Position>* multi_allelic_snps, set<Position>* non_std_chr_snps);

// Combines data from RareMetal's Covariance file with data from score file.
bool GetSigmaFromRareMetalFiles(
    const int num_samples,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    const map<Position, SnpInfo>& scores,
    double* sigma);

// Some Input formats require rescaling of some of the score variables.
void RescaleSnpInfo(
    const double& u_stat_rescale, const double& sqrt_v_stat_rescale,
    map<Position, SnpInfo>* scores);
// Some Input formats require rescaling of the covariances.
void RescaleCovariances(const double& rescale_amount,
                        map<pair<Position, Position>, double>* covariances);

// Converts snp_to_cov_w_neighbors (Sliding Window info on covariances) to
// gene-centric info on covariances.
bool ComputeCovarianceByGene(
    const int num_samples, const double& sigma_sq,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances);

// Make sure we don't have partially info on any snps (i.e. a SNP appears
// in Score or Covariance file, but no self-covariance is available).
bool SanityCheckSnpsToSkip(
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, SnpInfo>& scores,
    const map<pair<Position, Position>, double>& covariances);

}  // namespace premeta
#endif
