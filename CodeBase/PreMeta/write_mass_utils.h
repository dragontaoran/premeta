// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Using provided inputs, Writes a file that matches the format
//              of 1st-phase (summary statistics) output of MASS (ScoreSeq)
//              (and can be subsequently used by MASS for meta-analysis).

#ifndef WRITE_MASS_UTILS_H
#define WRITE_MASS_UTILS_H


#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "FileReaderUtils/vcf_utils.h"

#include <set>
#include <string>
#include <vector>

using file_reader_utils::Nucleotide;

namespace premeta {

// Uses the input maps to construct the MASS (ScoreSeq) output file, and
// prints it to 'out_file'.
bool PrintToMass(
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const string& out_file,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& scores,
    const map<pair<Position, Position>, double>& covariances);

bool PrintMassScriptFile(const string& mass_score_file, const string& outfile);

}  // namespace premeta
#endif
