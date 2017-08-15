// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Converts 1st-phase output of MASS, RAREMETAL, MetaSKAT, or
// SeqMeta to any of the other formats, so the 2nd-phase (meta-analysis) can
// be performed on data aggregated across studies that had utilized different
// software for the first phase.

#ifndef RM_QC_H
#define RM_QC_H

#include "FileReaderUtils/vcf_utils.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <cfloat>   // For DBL_MIN
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;
using file_reader_utils::Chromosome;
using file_reader_utils::Nucleotide;

namespace premeta {

// Goes through RareMetal score and covariance files, looking for and removing
// monomorphic and multi-allelic SNPs. Keeps track of removed SNPs in the
// provided containers.
extern bool RareMetalQc(
    const bool is_new_version_from, const double& rescale,
    const bool monomorphic_only, const bool multi_allelic_only,
    const bool remove_multi_allelic_two, const bool remove_non_snps,
    const bool strict_matching,
    const string& in_score_file, const string& in_cov_file,
    const string& out_score_file, const string& out_cov_file,
    set<Position>* monomorphic_snps, set<Position>* multi_allelic_snps,
    set<Position>* multi_allelic_snps_two, set<Position>* non_snps,
    set<Position>* non_matching, set<Position>* non_std_chr_snps);

}  // namespace premeta

#endif
