// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "write_raremetal_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/snp_hwe.h"
#include "StringUtils/string_utils.h"

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
using namespace string_utils;
using namespace std;

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
    const string& out_score_file, const string& out_cov_file) {
  // Get a set containing all SNP positions (will be used to fetch all
  // SNP's in a given window).
  const set<Position> snp_positions = Keys(windows);

  // Open score file for printing.
  ofstream score_outfile;
  score_outfile.open(out_score_file.c_str());
  if (!score_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_score_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Open covariates file for printing.
  ofstream cov_outfile;
  cov_outfile.open(out_cov_file.c_str());
  if (!cov_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_cov_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Print metadata to score file.
  const string version_num = "4.13.5";
  score_outfile << "##ProgramName=RareMetalWorker" << endl
                << "##Version=" << version_num << endl
                << "##Samples=" << num_samples << endl
                << "##AnalyzedSamples=" << num_samples << endl
                << "##Families=" << num_samples << endl
                << "##AnalyzedFamilies=" << num_samples << endl
                << "##Founders=" << num_samples << endl
                << "##AnalyzedFounders=" << num_samples << endl;

  // Print Header line to score file.
  score_outfile << "#CHROM" << kRareMetalScoreFileDefaultDelimiter << "POS"
                << kRareMetalScoreFileDefaultDelimiter << "REF"
                << kRareMetalScoreFileDefaultDelimiter << "ALT"
                << kRareMetalScoreFileDefaultDelimiter << "N_INFORMATIVE"
                << kRareMetalScoreFileDefaultDelimiter << "FOUNDER_AF"
                << kRareMetalScoreFileDefaultDelimiter << "ALL_AF"
                << kRareMetalScoreFileDefaultDelimiter << "INFORMATIVE_ALT_AC"
                << kRareMetalScoreFileDefaultDelimiter << "CALL_RATE"
                << kRareMetalScoreFileDefaultDelimiter << "HWE_PVALUE"
                << kRareMetalScoreFileDefaultDelimiter << "N_REF"
                << kRareMetalScoreFileDefaultDelimiter << "N_HET"
                << kRareMetalScoreFileDefaultDelimiter << "N_ALT"
                << kRareMetalScoreFileDefaultDelimiter << "U_STAT"
                << kRareMetalScoreFileDefaultDelimiter << "SQRT_V_STAT"
                << kRareMetalScoreFileDefaultDelimiter << "ALT_EFFSIZE"
                << kRareMetalScoreFileDefaultDelimiter << "PVALUE"
                << endl;

  // Print header line to covariates file.
  cov_outfile << "#CHROM" << kRareMetalCovFileDefaultDelimiter << "CURRENT_POS"
              << kRareMetalCovFileDefaultDelimiter << "MARKERS_IN_WINDOW"
              << kRareMetalCovFileDefaultDelimiter << "COV_MATRICES\n";

  // Print contents of both files.
  for (const pair<Position, SnpInfo>& info : snp_info) {
    const Position& snp_pos = info.first;
    const SnpInfo& snp_info = info.second;
    
    const bool is_monomorphic =
        monomorphic_snps.find(snp_pos) != monomorphic_snps.end();

    // Check if SNP has REF/ALT swapped from golden file.
    const bool is_ref_alt_swapped =
        IsSnpRefAltSwapped(study_num, snp_pos, snp_to_ref_alt_and_study);

    // Skip flagged SNPs.
    if (snps_to_skip.find(snp_pos) != snps_to_skip.end()) continue;

    // Print CHROM and POS to both files.
    score_outfile << snp_pos.chr_ << kRareMetalScoreFileDefaultDelimiter
                  << snp_pos.pos_ << kRareMetalScoreFileDefaultDelimiter;
    if (!is_monomorphic) {
      cov_outfile << snp_pos.chr_ << kRareMetalCovFileDefaultDelimiter
                  << snp_pos.pos_ << kRareMetalCovFileDefaultDelimiter;
    }

    // Print REF and ALT (Reference and Minor Allele types).
    if (is_ref_alt_swapped) {
      if (snp_info.minor_allele_.empty()) {
        score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
      } else {
        score_outfile << snp_info.minor_allele_
                      << kRareMetalScoreFileDefaultDelimiter;
      }
    } else if (snp_info.major_allele_.empty()) {
      score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.major_allele_
                    << kRareMetalScoreFileDefaultDelimiter;
    }
    if (is_ref_alt_swapped) {
      if (snp_info.major_allele_.empty()) {
        score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
      } else {
        score_outfile << snp_info.major_allele_
                      << kRareMetalScoreFileDefaultDelimiter;
      }
    } else if (snp_info.minor_allele_.empty()) {
      score_outfile << kDummyAltAllele << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.minor_allele_
                    << kRareMetalScoreFileDefaultDelimiter;
    }
    
    // For num_samples_for_gene, use num_non_missing_ or missing_rate_ info, if available.
    const bool has_hwe_info =
        snp_info.n_het_ >= 0 && snp_info.n_alt_ >= 0 && snp_info.n_ref_ >= 0;
    const int num_samples_for_snp =
        snp_info.num_non_missing_ >= 0 ? snp_info.num_non_missing_ :
        has_hwe_info ? snp_info.n_het_ + snp_info.n_alt_ + snp_info.n_ref_ :
        snp_info.missing_rate_ >= 0.0 ?
        round(num_samples * (1.0 - snp_info.missing_rate_)) : num_samples;
    const int mac = snp_info.mac_ < 0 ?
        round(snp_info.maf_ * 2.0 * num_samples_for_snp) : snp_info.mac_;
    // We intentionally round down (by casting as an int instead of using
    // 'round' function), as rounding up may force n_het to be negative.
    const int n_alt = snp_info.n_alt_ < 0 ?
         num_samples_for_snp * pow(snp_info.maf_, 2.0) : snp_info.n_alt_;
    const int n_het = snp_info.n_het_ < 0 ?
        mac - (2 * n_alt) : snp_info.n_het_;
    const int n_ref = snp_info.n_ref_ < 0 ?
        num_samples_for_snp - (n_het + n_alt) : snp_info.n_ref_;
    const int info_alt_ac = (snp_info.n_het_ >= 0 && snp_info.n_alt_ >= 0) ?
        snp_info.n_het_ + 2 * snp_info.n_alt_ : mac;
    const double founder_af = has_hwe_info ?
        static_cast<double>(info_alt_ac) / (2.0 * num_samples_for_snp) :
        snp_info.maf_;
    const double call_rate =
        static_cast<double>(num_samples_for_snp) / num_samples;

    // Print N_INFORMATIVE, FOUNDER_AF, ALL_AF, and INFORMATIVE_ALT_AC.
    score_outfile << num_samples << kRareMetalScoreFileDefaultDelimiter;
    if (is_ref_alt_swapped) {
      score_outfile << (1.0 - founder_af) << kRareMetalScoreFileDefaultDelimiter
                    << (1.0 - founder_af) << kRareMetalScoreFileDefaultDelimiter
                    << (n_het + 2 * n_ref)
                    << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.maf_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.maf_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.mac_ << kRareMetalScoreFileDefaultDelimiter;
    }

    // Print CALL_RATE (if not available, just print samples_for_gene).
    score_outfile << call_rate << kRareMetalScoreFileDefaultDelimiter;
    // Print HWE_PVALUE (if not available, just print 1).
    const double hwe = has_hwe_info ?
        GetSnpHWE(snp_info.n_het_, snp_info.n_ref_, snp_info.n_alt_) : 1.0;
    score_outfile << (hwe < 0.0 ? 1.0 : hwe)
                  << kRareMetalScoreFileDefaultDelimiter;

    // Print N_REF (estimate as: Num_samples - MAC if not available),
    // N_HET (estimate with MAC if not available), and
    // N_ALT (set to 0 if not available).
    if (is_ref_alt_swapped) {
      score_outfile << n_alt << kRareMetalScoreFileDefaultDelimiter
                    << n_het << kRareMetalScoreFileDefaultDelimiter
                    << n_ref << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << n_ref << kRareMetalScoreFileDefaultDelimiter
                    << n_het << kRareMetalScoreFileDefaultDelimiter
                    << n_alt << kRareMetalScoreFileDefaultDelimiter;
    }

    // Don't print out this SNP to covariance file if it is a Monomorphic SNP
    // (also print 'NA' in the remaining 4 score file columns: U_STAT,
    // SQRT_V_STAT, ALT_EFFSIZE, PVALUE).
    if (is_monomorphic) {
      score_outfile << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA\n";
      continue;
    }

    // Print the remainder of info for this SNP to the cov_outfile.
    // Find window for this SNP.
    double self_covariance = 0.0;
    if (!LookupCovariance(snp_pos, snp_pos, covariances, &self_covariance)) {
      cout << "ERROR: Unable to find self-covariance for SNP Position: "
           << PrintPosition(snp_pos) << endl;
      return false;
    }
    if (!FloatEq(rescale, 1.0) || !FloatEq(cov_rescale, 1.0)) {
      self_covariance *= (cov_rescale / (rescale * rescale));
    }

    map<Position, Position>::const_iterator window_itr = windows.find(snp_pos);
    set<Position>::const_iterator snp_itr = snp_positions.find(snp_pos);
    if (window_itr == windows.end() || snp_itr == snp_positions.end()) {
      cout << "ERROR Unable to find covariance info for SNP Position "
           << PrintPosition(snp_pos) << ". Aborting." << endl;
      return false;
    }
    const Position& window_end = window_itr->second;
    string cov_matrices = "";
    bool is_first = true;
    while (snp_itr != snp_positions.end()) {
      // Skip flagged SNPs.
      if (snps_to_skip.find(*snp_itr) != snps_to_skip.end()) continue;

      // Add separating comma.
      if (is_first) {
        is_first = false;
      } else {
        cov_outfile << ",";
        cov_matrices += ",";
      }

      // Print Position.
      cov_outfile << (*snp_itr).pos_;

      // Print Covariance (we ignore the return value of LookupCovariance,
      // as there are valid times it will fail to find a covariance: no
      // covariate info will be available for a given SNP and all of the
      // SNPs in its window that aren't in its gene.
      double covariance = 0.0;
      if (monomorphic_snps.find(*snp_itr) == monomorphic_snps.end()) {
        LookupCovariance(snp_pos, *snp_itr, covariances, &covariance);
        if (!FloatEq(rescale, 1.0) || !FloatEq(cov_rescale, 1.0)) {
          covariance *= (cov_rescale / (rescale * rescale));
        }
        covariance *= 1.0 / static_cast<double>(num_samples);
      }
      const bool is_ref_alt_swapped_two =
        IsSnpRefAltSwapped(study_num, *snp_itr, snp_to_ref_alt_and_study);
      const double swapped_multiplier =
          (is_ref_alt_swapped && is_ref_alt_swapped_two) ? 1.0 :
          (is_ref_alt_swapped || is_ref_alt_swapped_two) ? -1.0 : 1.0;

      cov_matrices += Itoa(swapped_multiplier * covariance);

      // Check if we've reached the window's end Position.
      if ((*snp_itr).chr_ == window_end.chr_ &&
          (*snp_itr).pos_ == window_end.pos_) {
        break;
      }
      snp_itr++;
    }
    cov_outfile << kRareMetalCovFileDefaultDelimiter << cov_matrices << "\n";

    // Print U-Stat and Sqrt_V-Stat.
    const double sqrt_v_stat = sqrt(self_covariance) / (rescale * rescale);
    const double u_stat_multiplier = is_ref_alt_swapped ? -1.0 : 1.0;
    const double u_stat =
        u_stat_multiplier * u_stat_rescale * snp_info.u_stat_ / rescale;
    score_outfile << u_stat << kRareMetalScoreFileDefaultDelimiter
                  << sqrt_v_stat << kRareMetalScoreFileDefaultDelimiter;

    // Print ALT_EFFSIZE.
    const double alt_effsize = u_stat / self_covariance;
    score_outfile << alt_effsize << kRareMetalScoreFileDefaultDelimiter;

    // Print p-value.
    const double u_stat_v_stat_ratio = abs(u_stat / sqrt_v_stat);

    // Want to take \Phi(u_stat_v_stat_ratio), where \Phi represents the CDF
    // for the standard normal distribution. C++ doesn't have an out-of-the-
    // box \Phi function, but it does have the error function erf and its
    // complement erfc, which are related to \Phi via:
    //    \Phi(x) = 0.5 * erfc(-x / sqrt(2))
    //cout << "\nPHB u_stat_v_stat_ratio: " << u_stat_v_stat_ratio << endl;
    const double u_stat_norm =
        0.5 * erfc(-1.0 * u_stat_v_stat_ratio / sqrt(2.0));
    score_outfile << (2.0 - 2.0 * u_stat_norm);
    score_outfile << "\n";
  }

  score_outfile.close();
  cov_outfile.close();
  return true;
}

bool PrintToRareMetal(
    const int study_num, const int num_samples,
    const double& rescale, const double& u_stat_rescale,
    const set<Position>& monomorphic_snps,
    const set<Position>& snps_without_score_info,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    const string& out_score_file, const string& out_cov_file) {
  // Open score file for printing.
  ofstream score_outfile;
  score_outfile.open(out_score_file.c_str());
  if (!score_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_score_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Open covariates file for printing.
  ofstream cov_outfile;
  cov_outfile.open(out_cov_file.c_str());
  if (!cov_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_cov_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Print metadata to score file.
  const string version_num = "4.13.5";
  score_outfile << "##ProgramName=RareMetalWorker" << endl
                << "##Version=" << version_num << endl
                << "##Samples=" << num_samples << endl
                << "##AnalyzedSamples=" << num_samples << endl
                << "##Families=" << num_samples << endl
                << "##AnalyzedFamilies=" << num_samples << endl
                << "##Founders=" << num_samples << endl
                << "##AnalyzedFounders=" << num_samples << endl;

  // Print Header line to score file.
  score_outfile << "#CHROM" << kRareMetalScoreFileDefaultDelimiter << "POS"
                << kRareMetalScoreFileDefaultDelimiter << "REF"
                << kRareMetalScoreFileDefaultDelimiter << "ALT"
                << kRareMetalScoreFileDefaultDelimiter << "N_INFORMATIVE"
                << kRareMetalScoreFileDefaultDelimiter << "FOUNDER_AF"
                << kRareMetalScoreFileDefaultDelimiter << "ALL_AF"
                << kRareMetalScoreFileDefaultDelimiter << "INFORMATIVE_ALT_AC"
                << kRareMetalScoreFileDefaultDelimiter << "CALL_RATE"
                << kRareMetalScoreFileDefaultDelimiter << "HWE_PVALUE"
                << kRareMetalScoreFileDefaultDelimiter << "N_REF"
                << kRareMetalScoreFileDefaultDelimiter << "N_HET"
                << kRareMetalScoreFileDefaultDelimiter << "N_ALT"
                << kRareMetalScoreFileDefaultDelimiter << "U_STAT"
                << kRareMetalScoreFileDefaultDelimiter << "SQRT_V_STAT"
                << kRareMetalScoreFileDefaultDelimiter << "ALT_EFFSIZE"
                << kRareMetalScoreFileDefaultDelimiter << "PVALUE"
                << endl;

  // Print header line to covariates file.
  cov_outfile << "#CHROM" << kRareMetalCovFileDefaultDelimiter << "CURRENT_POS"
              << kRareMetalCovFileDefaultDelimiter << "MARKERS_IN_WINDOW"
              << kRareMetalCovFileDefaultDelimiter << "COV_MATRICES\n";

  // Print contents of both files.
  for (const pair<Position, SnpInfo>& info : snp_info) {
    const Position& snp_pos = info.first;
    const SnpInfo& snp_info = info.second;

    const bool is_monomorphic =
        monomorphic_snps.find(snp_pos) != monomorphic_snps.end();

    // Check if SNP has REF/ALT swapped from golden file.
    const bool is_ref_alt_swapped =
        IsSnpRefAltSwapped(study_num, snp_pos, snp_to_ref_alt_and_study);

    // Print CHROM and POS to both files.
    score_outfile << snp_pos.chr_ << kRareMetalScoreFileDefaultDelimiter
                  << snp_pos.pos_ << kRareMetalScoreFileDefaultDelimiter;
    if (!is_monomorphic) {
      cov_outfile << snp_pos.chr_ << kRareMetalCovFileDefaultDelimiter
                  << snp_pos.pos_ << kRareMetalCovFileDefaultDelimiter;
    }

    // Print REF and ALT (Reference and Minor Allele types).
    if (is_ref_alt_swapped) {
      if (snp_info.minor_allele_.empty()) {
        score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
      } else {
        score_outfile << snp_info.minor_allele_
                      << kRareMetalScoreFileDefaultDelimiter;
      }
    } else if (snp_info.major_allele_.empty()) {
      score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.major_allele_
                    << kRareMetalScoreFileDefaultDelimiter;
    }
    if (is_ref_alt_swapped) {
      if (snp_info.major_allele_.empty()) {
        score_outfile << kDummyRefAllele << kRareMetalScoreFileDefaultDelimiter;
      } else {
        score_outfile << snp_info.major_allele_
                      << kRareMetalScoreFileDefaultDelimiter;
      }
    } else if (snp_info.minor_allele_.empty()) {
      score_outfile << kDummyAltAllele << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.minor_allele_
                    << kRareMetalScoreFileDefaultDelimiter;
    }
    
    // For num_samples_for_gene, use num_non_missing_ or missing_rate_ info,
    // if available.
    if (snp_info.n_het_ < 0 || snp_info.n_alt_ < 0 || snp_info.n_ref_ < 0 ||
        snp_info.mac_ < 0 || snp_info.num_non_missing_ < 0) {
      cout << "ERROR: N_HET, N_ALT, N_REF, MAC, or CALL_RATE information "
           << "missing for SNP " << PrintPosition(snp_pos) << endl;
      return false;
    }
    const double call_rate =
        static_cast<double>(snp_info.num_non_missing_) / num_samples;

    // Print N_INFORMATIVE, FOUNDER_AF, ALL_AF, INFORMATIVE_ALT_AC,
    // and CALL_RATE.
    score_outfile << num_samples << kRareMetalScoreFileDefaultDelimiter;
    if (is_ref_alt_swapped) {
      score_outfile << (1.0 - snp_info.maf_) << kRareMetalScoreFileDefaultDelimiter
                    << (1.0 - snp_info.maf_) << kRareMetalScoreFileDefaultDelimiter
                    << (snp_info.n_het_ + 2 * snp_info.n_ref_)
                    << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.maf_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.maf_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.mac_ << kRareMetalScoreFileDefaultDelimiter;
    }
    score_outfile << call_rate << kRareMetalScoreFileDefaultDelimiter;

    // Print HWE_PVALUE (if not available, just print 1).
    const double hwe =
        GetSnpHWE(snp_info.n_het_, snp_info.n_ref_, snp_info.n_alt_);
    score_outfile << (hwe < 0.0 ? 1.0 : hwe)
                  << kRareMetalScoreFileDefaultDelimiter;

    // Print N_REF (estimate as: Num_samples - MAC if not available),
    // N_HET (estimate with MAC if not available), and
    // N_ALT (set to 0 if not available).
    if (is_ref_alt_swapped) {
      score_outfile << snp_info.n_alt_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.n_het_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.n_ref_ << kRareMetalScoreFileDefaultDelimiter;
    } else {
      score_outfile << snp_info.n_ref_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.n_het_ << kRareMetalScoreFileDefaultDelimiter
                    << snp_info.n_alt_ << kRareMetalScoreFileDefaultDelimiter;
    }

    // Don't print out this SNP to covariance file if it is a Monomorphic SNP
    // (also print 'NA' in the remaining 4 score file columns: U_STAT,
    // SQRT_V_STAT, ALT_EFFSIZE, PVALUE).
    if (is_monomorphic) {
      score_outfile << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA" << kRareMetalScoreFileDefaultDelimiter
                    << "NA\n";
      continue;
    }

    // Print Covariance.
    double self_covariance;
    if (!PrintCovariance(
            study_num, num_samples, rescale, snp_pos, snp_to_ref_alt_and_study,
            snp_to_cov_w_neighbors, cov_outfile, &self_covariance)) {
      return false;
    }

    // Print U-Stat and Sqrt_V-Stat.
    const double sqrt_v_stat =
        sqrt(self_covariance) / (rescale * rescale);
    const double u_stat_multiplier = is_ref_alt_swapped ? -1.0 : 1.0;
    const double u_stat =
        u_stat_multiplier * u_stat_rescale * snp_info.u_stat_ / rescale;
    score_outfile << u_stat << kRareMetalScoreFileDefaultDelimiter
                  << sqrt_v_stat << kRareMetalScoreFileDefaultDelimiter;

    // Print ALT_EFFSIZE.
    const double alt_effsize = u_stat / (self_covariance);
    score_outfile << alt_effsize << kRareMetalScoreFileDefaultDelimiter;

    // Print p-value.
    const double u_stat_v_stat_ratio = abs(u_stat / sqrt_v_stat);


    // Want to take \Phi(u_stat_v_stat_ratio), where \Phi represents the CDF
    // for the standard normal distribution. C++ doesn't have an out-of-the-
    // box \Phi function, but it does have the error function erf and its
    // complement erfc, which are related to \Phi via:
    //    \Phi(x) = 0.5 * erfc(-x / sqrt(2))
    const double u_stat_norm =
        0.5 * erfc(-1.0 * u_stat_v_stat_ratio / sqrt(2.0));
    score_outfile << (2.0 - 2.0 * u_stat_norm);
    score_outfile << "\n";
  }

  // Now loop through all SNPs that had covariate information but no score
  // information, and print out the covariate information.
  for (const Position& snp_pos : snps_without_score_info) {
    if (monomorphic_snps.find(snp_pos) != monomorphic_snps.end()) {
      continue;
    }
    cov_outfile << snp_pos.chr_ << kRareMetalCovFileDefaultDelimiter
                << snp_pos.pos_ << kRareMetalCovFileDefaultDelimiter;
    double self_covariance;
    if (!PrintCovariance(
            study_num, num_samples, rescale, snp_pos, snp_to_ref_alt_and_study,
            snp_to_cov_w_neighbors, cov_outfile, &self_covariance)) {
      return false;
    }
  }

  score_outfile.close();
  cov_outfile.close();
  return true;
}

bool PrintRareMetaliWithMonoAndMultiSnpsRemoved(
    const string& in_score_file, const string& in_cov_file,
    const string& out_score_file, const string& out_cov_file,
    const set<Position>& non_std_chr_snps,
    const set<Position>& snps_to_exclude,
    const map<Position, int>& snps_in_score_file,
    const map<Position, int>& snps_in_cov_file,
    set<Position>* non_matching) {
  // UPDATE: Ran wanted to print out copies of .score and .cov file, even
  // if they are unchanged. So skip this check.
  /*
  if (snps_to_exclude.empty()) {
    cout << "ERROR: No Monomorphic nor Multi-Allelic SNPs to remove." << endl;
    return false;
  }
  */

  // Process Score File first. Need to simultaneously read from the old one and
  // write to the new one. This is straightforward copying, just omit lines
  // whose SNPs are in snps_to_exclude.
  //   - Open Score File for reading.
  ifstream score_file(in_score_file.c_str());
  if (!score_file.is_open()) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to open file."
         << endl;
    return false;
  }
  //   - Open score file for printing.
  ofstream score_outfile;
  score_outfile.open(out_score_file.c_str());
  if (!score_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_score_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }
  //   - Read line of score file
  set<string> sep;
  sep.insert(kRareMetalScoreFileDefaultDelimiter);
  int line_num = 0;
  string orig_line;
  while(getline(score_file, orig_line)) {
    line_num++;
    string line = orig_line;
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (HasPrefixString(line, kRareMetalScoreFileDefaultCommentChar)) {
      score_outfile << orig_line << endl;
      continue;
    } else if (RemoveLeadingWhitespace(line).empty()) {
      score_outfile << orig_line << endl;
      continue;
    } else if (HasPrefixString(line, "CHROM")) {
      score_outfile << orig_line << endl;
      continue;
    } else {
      // Non-comment line. Process Position, check if it is Monomorphic or
      // Multi-Allelic; if so, skip it. Otherwise, print line.
      vector<string> columns;
      Split(line, sep, false, &columns);
      if (columns.size() < 2) {
        cout << "ERROR: Unable to parse line " << line_num
             << " of score file '" << in_score_file
             << "': Line doesn't have expected number of columns:\n\t'"
             << orig_line << "'" << endl;
        return false;
      }
      // Position.
      uint64_t pos;
      if (!Stoi(columns[1], &pos)) {
        cout << "ERROR: Unable to parse line " << line_num
             << " of score file '" << in_score_file
             << "': Position (" << columns[1] << ") is not parsable as an "
             << "integer. Line:\n\t'" << orig_line << "'" << endl;
        return false;
      }
      // Chromosome.
      Chromosome chr;
      if (!VcfUtils::ParseChromosome(columns[0], &chr)) {
        Position chr_pos;
        chr_pos.chr_ = chr;
        chr_pos.pos_ = pos;
        if (non_std_chr_snps.find(chr_pos) == non_std_chr_snps.end()) {
          cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
               << "row " << line_num << " in " << in_score_file
               << ": Unrecognized chromosome '" << columns[0]
               << "'. Aborting.\n";
          return false;
        }
      }
      Position chr_pos;
      chr_pos.chr_ = chr;
      chr_pos.pos_ = pos;

      // Skip printing if this SNP is in snps_to_exclude.
      if (snps_to_exclude.find(chr_pos) != snps_to_exclude.end()) {
        continue;
      }

      // Skip printing this SNP if the number of times it appears in
      // snps_in_score_file doesn't match the number of times it appears
      // in snps_in_cov_file.
      if (non_matching != nullptr) {
        map<Position, int>::const_iterator score_snp_itr =
            snps_in_score_file.find(chr_pos);
        if (score_snp_itr == snps_in_score_file.end()) {
          non_matching->insert(chr_pos);
          continue;
        }
        map<Position, int>::const_iterator cov_snp_itr =
            snps_in_cov_file.find(chr_pos);
        if (cov_snp_itr == snps_in_cov_file.end()) {
          non_matching->insert(chr_pos);
          continue;
        }
        if (score_snp_itr->second != cov_snp_itr->second) {
          non_matching->insert(chr_pos);
          continue;
        }
      }

      // Print SNP.
      score_outfile << orig_line << endl;
    }
  }
  score_file.close();
  score_outfile.close();

  // Process Covariance File. Need to simultaneously read from the old one and
  // write to the new one. Two changes must be made:
  //   1) Lines corresponding to Monomorphic or Multi-Allelic SNPs should not
  //      be printed
  //   2) Within a window, whichever SNPs are removed due to (1) above, we need
  //      to update the covariance "matrices" by removing the covariances
  //      corresponding to the removed SNPs.
  // Step (2) is done using the 'MARKERS_IN_WINDOW' (for RAREMETALWORKER) or
  // 'MARKER_POS' (for RVTESTS) column.
  //   - Open Covariance File for reading.
  ifstream cov_file(in_cov_file.c_str());
  if (!cov_file.is_open()) {
    cout << "ERROR in getting info from RareMetal Covariance file: Unable to open file."
         << endl;
    return false;
  }
  //   - Open covariates file for printing.
  ofstream cov_outfile;
  cov_outfile.open(out_cov_file.c_str());
  if (!cov_outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_cov_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }
  //   - Read line of cov file
  line_num = 0;
  bool is_rvtests = false;
  while(getline(cov_file, orig_line)) {
    line_num++;
    string line = orig_line;
    string line_suffix = "";
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
      line_suffix = line.substr(line.length() - 1);
    }

    if (HasPrefixString(line, kRareMetalScoreFileDefaultCommentChar)) {
      cov_outfile << orig_line << endl;
      continue;
    } else if (HasPrefixString(line, "CHROM")) {
      is_rvtests = true;
      cov_outfile << orig_line << endl;
      continue;
    } else {
      // Non-comment line. Process Position, check if it is Monomorphic or
      // Multi-Allelic; if so, skip it. Otherwise, print line.
      vector<string> columns;
      Split(line, sep, false, &columns);
      if ((is_rvtests && columns.size() != 6) ||
           !is_rvtests && columns.size() != 4) {
        cout << "ERROR: Unable to parse line " << line_num
             << " of cov file '" << in_cov_file
             << "': Line doesn't have expected number of columns:\n\t'"
             << orig_line << "'" << endl;
        return false;
      }
      // Position.
      uint64_t pos;
      if (!Stoi(columns[1], &pos)) {
        cout << "ERROR: Unable to parse line " << line_num
             << " of score file '" << in_score_file
             << "': Position (" << columns[1] << ") is not parsable as an "
             << "integer. Line:\n\t'" << orig_line << "'" << endl;
        return false;
      }
      // Chromosome.
      Chromosome chr;
      if (!VcfUtils::ParseChromosome(columns[0], &chr)) {
        Position chr_pos;
        chr_pos.chr_ = chr;
        chr_pos.pos_ = pos;
        if (non_std_chr_snps.find(chr_pos) == non_std_chr_snps.end()) {
          cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
               << "row " << line_num << " in " << in_score_file
               << ": Unrecognized chromosome '" << columns[0]
               << "'. Aborting.\n";
          return false;
        }
      }
      Position chr_pos;
      chr_pos.chr_ = chr;
      chr_pos.pos_ = pos;
      if (snps_to_exclude.find(chr_pos) != snps_to_exclude.end()) {
        // Monomorphic, Multi-Allelic, or Wrong-Chr SNP; skip line.
        continue;
      }

      // Skip printing this SNP if the number of times it appears in
      // snps_in_score_file doesn't match the number of times it appears
      // in snps_in_cov_file.
      if (non_matching != nullptr) {
        map<Position, int>::const_iterator cov_snp_itr =
            snps_in_cov_file.find(chr_pos);
        if (cov_snp_itr == snps_in_cov_file.end()) {
          non_matching->insert(chr_pos);
          continue;
        }
        map<Position, int>::const_iterator score_snp_itr =
            snps_in_score_file.find(chr_pos);
        if (score_snp_itr == snps_in_score_file.end()) {
          non_matching->insert(chr_pos);
          continue;
        }
        if (score_snp_itr->second != cov_snp_itr->second) {
          non_matching->insert(chr_pos);
          continue;
        }
      }

      // SNP is not Monomorphic or Multi-Allelic; but perhaps its neighbors are.
      vector<string> orig_positions_in_window, positions_in_window;
      set<int> indices_to_skip;
      const int markers_in_window_col_index = is_rvtests ? 4 : 2;
      Split(columns[markers_in_window_col_index], ",", &orig_positions_in_window);
      for (int i = 0; i < orig_positions_in_window.size(); ++i) {
        // Parse position.
        uint64_t marker_pos;
        if (!Stoi(orig_positions_in_window[i], &marker_pos)) {
          cout << "ERROR: Unable to parse line " << line_num
               << " of score file '" << in_score_file
               << "': In MARKERS_IN_WINDOW, item " << i + 1
               << " Position (" << orig_positions_in_window[i]
               << ") is not parsable as an integer. Line:\n\t'"
               << orig_line << "'" << endl;
          return false;
        }
        Position neighbor_pos;
        neighbor_pos.chr_ = chr;
        neighbor_pos.pos_ = marker_pos;
        if (snps_to_exclude.find(neighbor_pos) != snps_to_exclude.end()) {
          // Monomorphic or Multi-Allelic SNP.
          indices_to_skip.insert(i);
        } else if (non_matching != nullptr) {
          map<Position, int>::const_iterator cov_snp_itr =
              snps_in_cov_file.find(neighbor_pos);
          if (cov_snp_itr == snps_in_cov_file.end()) {
            indices_to_skip.insert(i);
            continue;
          }
          map<Position, int>::const_iterator score_snp_itr =
              snps_in_score_file.find(neighbor_pos);
          if (score_snp_itr == snps_in_score_file.end()) {
            indices_to_skip.insert(i);
            continue;
          }
          if (score_snp_itr->second != cov_snp_itr->second) {
            indices_to_skip.insert(i);
            continue;
          }
          positions_in_window.push_back(orig_positions_in_window[i]);
        } else {
          positions_in_window.push_back(orig_positions_in_window[i]);
        }
      }

      // If no neighbors are Monomorphic or Multi-Allelic, just print original line.
      if (indices_to_skip.empty()) {
        cov_outfile << orig_line << endl;
        continue;
      }
      // There was at least one monomorphic/multi-allelic SNP. Reconstruct the
      // line, picking out the relevant parts.
      cov_outfile << columns[0] << "\t" << columns[1] << "\t";
      if (is_rvtests) {
        cov_outfile << columns[2] << "\t"
                    << Itoa(positions_in_window.size()) << "\t";
      }
      cov_outfile << Join(positions_in_window, ",") << "\t";
      vector<string> neighbor_covariances;
      const int covariances_col_index = is_rvtests ? 5 : 3;
      Split(columns[covariances_col_index], ",", &neighbor_covariances);
      if (neighbor_covariances.size() != orig_positions_in_window.size()) {
        cout << "ERROR: Unable to parse line " << line_num
             << " of score file '" << in_score_file
             << "': Number of positions in MARKERS_IN_WINDOW column ("
             << orig_positions_in_window.size() << ") doesn't match "
             << "number of covariances in COV_MATRICES column ("
             << neighbor_covariances.size() << "). Line:\n\t"
             << orig_line << "'" << endl;
        return false;
      }
      bool first_cov_printed = false;
      for (int i = 0; i < neighbor_covariances.size(); ++i) {
        if (indices_to_skip.find(i) != indices_to_skip.end()) continue;
        if (first_cov_printed) {
          cov_outfile << ",";
        }
        first_cov_printed = true;
        cov_outfile << neighbor_covariances[i];
      }
      cov_outfile << line_suffix << endl;
    }
  }

  cov_file.close();
  cov_outfile.close();

  return true;
}

bool PrintCovariance(
    const int study_num, const int num_samples, const double& rescale,
    const Position& snp_pos,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    ofstream& cov_outfile, double* self_covariance) {
  // Print the remainder of info for this SNP to the cov_outfile.
  map<Position, map<Position, double>>::const_iterator cov_itr =
      snp_to_cov_w_neighbors.find(snp_pos);
  if (cov_itr == snp_to_cov_w_neighbors.end()) {
    cout << "ERROR: Unable to find covariance info for SNP Position: "
         << PrintPosition(snp_pos) << endl;
    return false;
  }
  const map<Position, double>& neighbor_covariances = cov_itr->second;
  map<Position, double>::const_iterator self_cov_itr =
      neighbor_covariances.find(snp_pos);
  if (self_cov_itr == neighbor_covariances.end()) {
    cout << "ERROR: Unable to find self-covariance for SNP Position: "
         << PrintPosition(snp_pos) << endl;
    return false;
  }
  *self_covariance = self_cov_itr->second * num_samples;
  if (!FloatEq(rescale, 1.0)) {
    *self_covariance *= (1.0 / (rescale * rescale));
  }

  // Check if SNP has REF/ALT swapped from golden file.
  const bool is_ref_alt_swapped =
      IsSnpRefAltSwapped(study_num, snp_pos, snp_to_ref_alt_and_study);

  string cov_matrices = "";
  bool is_first = true;
  for (map<Position, double>::const_iterator neighbor_cov_itr =
           neighbor_covariances.begin();
       neighbor_cov_itr != neighbor_covariances.end(); ++neighbor_cov_itr) {
    // Add separating comma.
    if (is_first) {
      is_first = false;
    } else {
      cov_outfile << ",";
      cov_matrices += ",";
    }

    // Print Position.
    cov_outfile << (neighbor_cov_itr->first).pos_;

    // Print Covariance (we ignore the return value of LookupCovariance,
    // as there are valid times it will fail to find a covariance: no
    // covariate info will be available for a given SNP and all of the
    // SNPs in its window that aren't in its gene.
    double covariance = neighbor_cov_itr->second * num_samples;
    if (!FloatEq(rescale, 1.0)) {
      covariance *= (1.0 / (rescale * rescale));
    }
    covariance *= 1.0 / static_cast<double>(num_samples);

    const bool is_ref_alt_swapped_two =
      IsSnpRefAltSwapped(study_num, neighbor_cov_itr->first, snp_to_ref_alt_and_study);
    const double swapped_multiplier =
        (is_ref_alt_swapped && is_ref_alt_swapped_two) ? 1.0 :
        (is_ref_alt_swapped || is_ref_alt_swapped_two) ? -1.0 : 1.0;
    cov_matrices += Itoa(swapped_multiplier * covariance);
  }

  cov_outfile << kRareMetalCovFileDefaultDelimiter << cov_matrices << "\n";
  return true;
}

bool AddFilesToList(
    const string& score_files, const string& cov_files,
    const string& score_file, const string& cov_file) {
  ofstream score_out_file;
  score_out_file.open(score_files.c_str(), ios::out | ios::app);
  if (!score_out_file.is_open()) {
    cout << "ERROR in trying to write file '" << score_files
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  score_out_file << "## === STUDY INFORMATION === ##" << endl;
  score_out_file << score_file << endl;
  score_out_file.close();

  ofstream cov_out_file;
  cov_out_file.open(cov_files.c_str(), ios::out | ios::app);
  if (!cov_out_file.is_open()) {
    cout << "ERROR in trying to write file '" << cov_files
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  cov_out_file << "## === STUDY INFORMATION === ##" << endl;
  cov_out_file << cov_file << endl;
  cov_out_file.close();

  return true;
}

}  // namespace premeta
