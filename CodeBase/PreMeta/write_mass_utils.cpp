// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "write_mass_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_mass_utils.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using file_reader_utils::Chromosome;
using file_reader_utils::Nucleotide;
using file_reader_utils::CsvUtils;
using file_reader_utils::VcfUtils;
using file_reader_utils::GenericDataHolder;
using file_reader_utils::GenericDataType;
using namespace map_utils;
using namespace math_utils;
using namespace std;

namespace premeta {

bool PrintToMass(
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const string& out_file,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& scores,
    const map<pair<Position, Position>, double>& covariances) {
  if (FloatEq(rescale, 0.0) || FloatEq(sigma_sq, 0.0)) {
    cout << "ERROR: Invalid value for rescale (" << rescale
         << ") or sigma_sq (" << sigma_sq << "): they must be non-zero. "
         << "Aborting." << endl;
    return false;
  }

  ofstream outfile;
  outfile.open(out_file.c_str());
  if (!outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }
  outfile << "#Samples = " << num_samples << endl;
  for (const auto& gene_itr : gene_to_snp_positions) {
    const string& gene = gene_itr.first;

    // Iterate through SNP positions in this gene.
    for (const Position& pos : gene_itr.second) {
      // Check if SNP is monomorphic.
      const bool is_monomorphic =
          monomorphic_snps.find(pos) != monomorphic_snps.end();

      // Check if SNP has REF/ALT swapped from golden file.
      const bool is_ref_alt_swapped =
          IsSnpRefAltSwapped(study_num, pos, snp_to_ref_alt_and_study);

      // Print gene.
      outfile << gene << kMassDefaultDelimiter;

      // Print Position.
      outfile << PrintPosition(pos) << kMassDefaultDelimiter;
      
      // Print Score info.
      map<Position, SnpInfo>::const_iterator score_itr = scores.find(pos);
      if (score_itr == scores.end()) {
        cout << "ERROR in Printing Mass File: unable to find score information "
             << "for Position " << PrintPosition(pos) << " on gene '" << gene
             << "'. Aborting." << endl;
        return false;
      } else {
        const SnpInfo& info = score_itr->second;
        // For num_samples, use num_non_missing_ or missing_rate_ info, if available.
        const bool has_hwe_info =
            info.n_het_ >= 0 && info.n_alt_ >= 0 && info.n_ref_ >= 0;
        const int num_samples_for_snp =
            info.num_non_missing_ >= 0 ? info.num_non_missing_ :
            has_hwe_info ? info.n_het_ + info.n_alt_ + info.n_ref_ :
            info.missing_rate_ >= 0.0 ?
            round(num_samples * (1.0 - info.missing_rate_)) : num_samples;
        double maf = info.maf_ < 0 ? 
            info.mac_ / (2.0 * num_samples_for_snp) : info.maf_;
        if (is_ref_alt_swapped) {
          maf = 1.0 - maf;
        }
        int mac = info.mac_ < 0 ?
            round(info.maf_ * 2.0 * num_samples_for_snp) : info.mac_;
        int n_alt = info.n_alt_ < 0 ?
            num_samples_for_snp * pow(info.maf_, 2.0) : info.n_alt_;
        const int n_het = info.n_het_ < 0 ? mac - (2 * n_alt) : info.n_het_;
        int n_ref = info.n_ref_ < 0 ?
            num_samples_for_snp - (n_het + n_alt) : info.n_ref_;
        if (is_ref_alt_swapped) {
          int temp = n_ref;
          n_ref = n_alt;
          n_alt = temp;
        }
        if (is_ref_alt_swapped) {
          mac = n_het + 2 * n_alt;
        }
        outfile << maf << kMassDefaultDelimiter << mac
                << kMassDefaultDelimiter << num_samples_for_snp
                << kMassDefaultDelimiter;
        outfile << n_ref << kMassDefaultDelimiter
                << n_het << kMassDefaultDelimiter
                << n_alt << kMassDefaultDelimiter;
        if (is_monomorphic) {
          outfile << 0.0;
        } else if (is_ref_alt_swapped) {
          outfile << (-1.0 * info.u_stat_ / (sigma_sq * rescale));
        } else {
          outfile << (info.u_stat_ / (sigma_sq * rescale));
        }
      }

      // Go through other SNPs on this gene, printing covariate info.
      for (const Position& pos_two : gene_itr.second) {
        if (is_monomorphic ||
            monomorphic_snps.find(pos_two) != monomorphic_snps.end()) {
          outfile << kMassDefaultDelimiter << 0.0;
          continue;
        }
        const bool is_ref_alt_swapped_two =
          IsSnpRefAltSwapped(study_num, pos_two, snp_to_ref_alt_and_study);
        const double swapped_multiplier =
            (is_ref_alt_swapped && is_ref_alt_swapped_two) ? 1.0 :
            (is_ref_alt_swapped || is_ref_alt_swapped_two) ? -1.0 : 1.0;
        double covariance = 0.0;
        if (!LookupCovariance(pos, pos_two, covariances, &covariance)) {
          cout << "ERROR in Printing Mass file: unable to find covariance for "
               << "positions " << PrintPosition(pos) << " and "
               << PrintPosition(pos_two) << endl;
          return false;
        }
        outfile << kMassDefaultDelimiter
                << (swapped_multiplier * covariance / (sigma_sq * rescale * rescale));
      }
      outfile << endl;
    }
  }
  outfile.close();
  return true;
}

bool PrintMassScriptFile(const string& mass_score_file, const string& outfile) {
  ofstream out_file;
  out_file.open(outfile.c_str(), ios::out | ios::app);
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  out_file << endl << "## === STUDY INFORMATION === ##" << endl;
  out_file << "FILE = " << mass_score_file << endl;
  out_file << "SKIP = 1" << endl;
  out_file << "GENE_ID_COLUMN = 1" << endl;
  out_file << "GVAR_ID_COLUMN = 2" << endl;
  // NOTE(11/20/15): The following changes are to keep consistent with changes
  // made to MASS v7.1.
  //out_file << "MAC_COLUMN = 4" << endl;
  //out_file << "N_OBS_COLUMN = 5" << endl;
  out_file << "MAF_COLUMN = 3" << endl;
  out_file << "SCORE_COLUMN = 9" << endl;

  out_file.close();
  return true;
}

}  // namespace premeta
