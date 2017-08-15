// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Using provided inputs, Writes .RData file, which mimics format
//              of 1st-phase (summary statistics) output of SeqMeta (and can
//              be subsequently used by SeqMeta for meta-analysis.

#ifndef WRITE_SEQMETA_UTILS_H
#define WRITE_SEQMETA_UTILS_H


#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "FileReaderUtils/vcf_utils.h"

#include <set>
#include <string>
#include <vector>

namespace premeta {

bool ComputeDsCMatrixFromCovarianceMatrix(
    const int study_num, const double& sigma_sq, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<pair<Position, Position>, double>& covariances,
    vector<int>* i_values, vector<int>* p_values, vector<double>* x_values);
bool ConstructSeqMetaTreeGeneCovI(
    vector<int>* i_values, vector<string*>* tags, RDataPairList* i);
bool ConstructSeqMetaTreeGeneCovP(
    vector<int>* p_values, vector<string*>* tags, RDataPairList* p);
bool ConstructSeqMetaTreeGeneCovDim(
    const int num_positions, vector<string*>* tags, RDataPairList* dim);
bool ConstructSeqMetaTreeGeneCovDimnames(
    const vector<Position>& gene_positions,
    vector<string*>* tags, RDataPairList* dimnames);
bool ConstructSeqMetaTreeGeneCovX(
    vector<double>* x_values, vector<string*>* tags, RDataPairList* x);
bool ConstructSeqMetaTreeGeneCovUplo(
    vector<string*>* tags, RDataPairList* uplo);
bool ConstructSeqMetaTreeGeneCovFactors(
    vector<string*>* tags, RDataPairList* factors);
bool ConstructSeqMetaTreeGeneCovClass(
    vector<string*>* tags, RDataPairList* class_name);
bool ConstructSeqMetaTreeGeneScores(
    const int study_num, const double& sey, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    vector<string*>* tags, RDataNode* gene_scores);
bool ConstructSeqMetaTreeGeneCov(
    const int study_num, const double& sigma_sq, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* gene_cov);
bool ConstructSeqMetaTreeGeneN(
    const int num_samples, RDataNode* gene_n);
bool ConstructSeqMetaTreeGeneMaf(
    const int study_num,
    const vector<Position>& gene_positions,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    vector<string*>* tags, RDataNode* gene_maf);
bool ConstructSeqMetaTreeGeneSey(
    const double& sey, RDataNode* gene_sey);
bool ConstructSeqMetaTreeGene(
    const int study_num, const int num_samples,
    const double& rescale, const double& sey,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* gene_info);
bool ConstructSeqMetaTreeRootObject(
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* root);
bool ConstructSeqMetaTreeRootAttribute(
    const map<string, vector<Position>>& gene_to_snp_positions,
    vector<string*>* tags, RDataNode* root_attr);
bool ConstructSeqMetaTree(
    const string& object_name,
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* root);

// Print SeqMeta outfile.
bool PrintToSeqMeta(
    const RDataNode& root, const string& out_file);

}  // namespace premeta
#endif
