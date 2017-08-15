// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Reads 1st-phase (summary statistics) output of MASS's
// ScoreSeq function, populating generic data structures.

#ifndef READ_MASS_UTILS_H
#define READ_MASS_UTILS_H

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "FileReaderUtils/vcf_utils.h"

#include <set>
#include <string>
#include <vector>

using file_reader_utils::Nucleotide;

namespace premeta {

// For the SNP positions in this gene, transfer information from snp_info
// into covariances (i.e. create covariance matrices for this gene).
bool UpdateGeneCovariancesFromMassFile(
    const vector<Position>& gene_snp_positions,
    const map<Position, vector<double>>& gene_covariances,
    map<pair<Position, Position>, double>* covariances);

// Read in Mass input file into covariances.
bool GetInfoFromMassFile(
    const bool require_snp_in_allele_file, const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    int* num_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<Position, SnpInfo>* snp_info,
    map<pair<Position, Position>, double>* covariances);

// Parse sigma_sq from first line comment in mass_file.
bool GetNumSamplesFromMassFile(
    const FileInfo& mass_file_info, int* num_samples);

}  // namespace premeta
#endif
