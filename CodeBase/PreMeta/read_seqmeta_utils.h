// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Reads 1st-phase (summary statistics) output of SEQMETA's
// ScoreSeq function, populating generic data structures.

#ifndef READ_SEQMETA_UTILS_H
#define READ_SEQMETA_UTILS_H

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <set>
#include <string>
#include <vector>

namespace premeta {

bool ParseSeqMetaRDataGeneScores(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    int* num_snps,
    const RDataNode* node,
    const string& gene_name,
    set<Position>* snps_to_skip,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores);
bool ParseSeqMetaRDataGeneCovarianceOld(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode* node,
    const string& gene_name,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances);
bool ParseSeqMetaRDataGeneCovariance(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode* node,
    const string& gene_name,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances);
bool ParseSeqMetaRDataGeneN(const RDataNode* node, int* num_samples);
bool ParseSeqMetaRDataGeneMaf(
    const int study_num,
    const int num_samples,
    const RDataNode* node,
    const string& gene_name,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* monomorphic_snps, map<Position, SnpInfo>* scores);
bool ParseSeqMetaRDataGeneSey(const RDataNode* node, double* sey);
bool ParseSeqMetaRDataMainObjectAttributeOne(
    const RDataPairList* first_attr, int* num_genes);
bool ParseSeqMetaRDataMainObjectAttributeTwo(
    const int num_genes, const string& tag_name,
    const RDataPairList* second_attr, vector<string>* genes);
bool ParseSeqMetaRDataMainObjectAttributeThree(
    const RDataPairList* third_attr);
bool ParseSeqMetaRDataMainObjectAttributeFour(
    const RDataPairList* fourth_attr);
bool ParseSeqMetaRDataMainObjectAttributes(
    const RDataNode& main_attr, bool* is_skat_cohort_version, vector<string>* genes);
bool ParseSeqMetaRDataGene(
    const int study_num,
    const bool is_skat_cohort_version,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const string& gene_name,
    const RDataNode* item, int* num_samples, int* num_snps, double* sey,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores,
    map<pair<Position, Position>, double>* covariances);
bool ParseSeqMetaRDataTree(
    const bool require_snp_in_allele_file, const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode& root,
    int* num_samples, int* num_snps, double* sey,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores,
    map<pair<Position, Position>, double>* covariances);

}  // namespace premeta
#endif
