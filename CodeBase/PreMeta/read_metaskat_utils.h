// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Reads 1st-phase (summary statistics) output of METASKAT's
// ScoreSeq function, populating generic data structures.

#ifndef READ_METASKAT_UTILS_H
#define READ_METASKAT_UTILS_H

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <set>
#include <string>
#include <vector>

namespace premeta {

uint32_t xcrc32(const unsigned char* buf, int len);

bool ParseMInfoHeader(
    const FileInfo& file_info,
    int* num_samples, int* num_genes, int* num_snps,
    int* num_unique_snps, int* num_header_lines);

bool ParseMInfoContents(
    const int study_num,
    const FileInfo& file_info, const map<Position, SnpInfo>& allele_info,
    const int num_header_lines, const int num_samples,
    const int num_genes, const int num_snps, const int num_unique_snps,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<string, uint64_t>* start_lines,
    map<Position, SnpInfo>* snp_info);

bool ParseMInfoFile(
    const int study_num,
    const FileInfo& file_info, const map<Position, SnpInfo>& allele_info,
    int* num_samples,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<string, uint64_t>* start_lines,
    map<Position, SnpInfo>* snp_info);

bool ParseMssdFile(
    const int study_num, const string& input_file,
    const map<Position, set<int>>& snp_to_excluding_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<string, uint64_t>& start_lines,
    map<pair<Position, Position>, double>* covariances);

}  // namespace premeta
#endif
