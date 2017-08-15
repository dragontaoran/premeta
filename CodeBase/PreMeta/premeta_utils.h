// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Converts 1st-phase output of MASS, RAREMETAL, MetaSKAT, or
// SeqMeta to any of the other formats, so the 2nd-phase (meta-analysis) can
// be performed on data aggregated across studies that had utilized different
// software for the first phase.

#ifndef PREMETA_UTILS_H
#define PREMETA_UTILS_H

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

class PreMeta {
 public:
  static bool WriteSeqMetaBinary(const string& outfile);

  static bool WriteGoldenAlleleFile(const string& outfile);
  static bool ParseGoldenSnpFile(const string& golden_file);
  static void SetSnpRefAltCheckField(const bool value) {
    check_snp_ref_alt_order_ = value;
  }

  // Prints snp_to_excluding_study_ to file.
  static void PrintInconsistentSnps(const string& outfile);

  // Reads alleles for each SNP position from file of form:
  //   CHR:POSITION  REF  ALT
  static bool GetAllelesFromAlleleFile(
      const bool is_rm_or_metaskat, const int study_num,
      const FileInfo& file_info,
      map<Position, SnpInfo>* allele_info);

  // DEPRECATED.
  // Reads alleles for each SNP position.
  //static bool GetAllelesFromAlleleInfo(
  //    const AlleleInfo& info,
  //    map<Position, pair<string, string>>* allele_info);

  // =========================================================================== 
  // MASS -> X
  // =========================================================================== 
  // =========================================================================== 
  //     MASS -> MASS
  // =========================================================================== 
  static bool MassToMass(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& mass_file_info,
      const string& mass_script_file, const string& out_file);
  // =========================================================================== 
  //     MASS -> MetaSKAT
  // =========================================================================== 
  static bool MassToMetaSkat(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& mass_file_info,
      const string& minfo_outfile, const string& mssd_outfile);
  // =========================================================================== 
  //     MASS -> RareMetal
  // =========================================================================== 
  static bool MassToRareMetal(
      const int study_num, const int window_size, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& mass_file_info,
      const string& score_files, const string& cov_files,
      const string& out_score_file, const string& out_cov_file);
  // =========================================================================== 
  //     MASS -> SeqMeta
  // =========================================================================== 
  static bool MassToSeqMeta(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& mass_file_info,
      const string& out_file);
 
  // =========================================================================== 
  // MetaSKAT -> X
  // =========================================================================== 
  // =========================================================================== 
  //     MetaSKAT -> MASS
  // =========================================================================== 
  static bool MetaSkatToMass(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& meta_skat_m_info,
      const string& meta_skat_mssd_file,
      const string& mass_script_file, const string& out_file);
  // =========================================================================== 
  //     MetaSKAT -> MetaSKAT
  // =========================================================================== 
  static bool MetaSkatToMetaSkat(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& meta_skat_m_info,
      const string& meta_skat_mssd_file,
      const string& minfo_outfile, const string& mssd_outfile);
  // =========================================================================== 
  //     MetaSKAT -> RareMetal
  // =========================================================================== 
  static bool MetaSkatToRareMetal(
      const int study_num, const int window_size, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& meta_skat_m_info,
      const string& meta_skat_mssd_file,
      const string& score_files, const string& cov_files,
      const string& out_score_file, const string& out_cov_file);
  // =========================================================================== 
  //     MetaSKAT -> SeqMeta
  // =========================================================================== 
  static bool MetaSkatToSeqMeta(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& meta_skat_m_info,
      const string& meta_skat_mssd_file,
      const string& out_file);

  // =========================================================================== 
  // RareMetal -> X
  // =========================================================================== 
  // =========================================================================== 
  //     RareMetal -> MASS
  // =========================================================================== 
  static bool RareMetalToMass(
      const int study_num, const bool is_new_version,
      const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& group_file_info,
      const FileInfo& score_file_info,
      const FileInfo& covariance_file_info,
      const string& mass_script_file, const string& out_file);
  // =========================================================================== 
  //     RareMetal -> MetaSKAT
  // =========================================================================== 
  static bool RareMetalToMetaSkat(
      const int study_num, const bool is_new_version, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& group_file_info,
      const FileInfo& score_file_info,
      const FileInfo& covariance_file_info,
      const string& minfo_outfile, const string& mssd_outfile);
  // =========================================================================== 
  //     RareMetal -> RareMetal
  // =========================================================================== 
  static bool RareMetalToRareMetal(
      const int study_num, const bool is_new_version_from,
      const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& score_file_info,
      const FileInfo& covariance_file_info,
      const string& score_files, const string& cov_files,
      const string& out_score_file, const string& out_cov_file);
  // =========================================================================== 
  //     RareMetal -> SeqMeta
  // =========================================================================== 
  static bool RareMetalToSeqMeta(
      const int study_num, const bool is_new_version, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const FileInfo& group_file_info,
      const FileInfo& score_file_info,
      const FileInfo& covariance_file_info,
      const string& out_file);

  // =========================================================================== 
  // SeqMeta -> X
  // =========================================================================== 
  // =========================================================================== 
  //     SeqMeta -> MASS
  // =========================================================================== 
  static bool SeqMetaToMass(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const string& seq_meta_rdata_file,
      const string& mass_script_file, const string& out_file);
  // =========================================================================== 
  //     SeqMeta -> MetaSKAT
  // =========================================================================== 
  static bool SeqMetaToMetaSkat(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const string& seq_meta_rdata_file,
      const string& minfo_outfile, const string& mssd_outfile);
  // =========================================================================== 
  //     SeqMeta -> RareMetal
  // =========================================================================== 
  static bool SeqMetaToRareMetal(
      const int study_num, const int window_size, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const string& seq_meta_rdata_file,
      const string& score_files, const string& cov_files,
      const string& out_score_file, const string& out_cov_file);
  // =========================================================================== 
  //     SeqMeta -> SeqMeta
  // =========================================================================== 
  static bool SeqMetaToSeqMeta(
      const int study_num, const double& rescale,
      const map<Position, SnpInfo>& allele_info,
      const string& seq_meta_rdata_file,
      const string& out_file);

 private:
  // Member fields.
  //   - Golden Set of REF/ALT alleles for each SNP, and the set of studies
  //     that have these flipped.
  static map<Position, tuple<string, string, set<int>>> snp_to_ref_alt_and_study_;
  //   - Set of (SNP, List of Studies) for which the SNP was excluded from analysis
  //     due to not matching (up to permutation of REF <-> ALT) the golden file.
  static map<Position, set<int>> snp_to_excluding_study_;
  //   - Whether to check each SNP's Ref/Alt alleles with the "golden file".
  static bool check_snp_ref_alt_order_;

  // Extracts just the filename (i.e. strips path prefix and ".RData" suffix).
  static string GetSeqMetaObjectNameFromFile(const string& out_file);

  // Prints all snps in 'snps_to_skip' to file.
  static void PrintSkippedSnps(
      const string& orig_out_file, const set<Position>& snps_to_skip);

  // For each SNP, finds the position of the LAST SNP in its window.
  // Only used for translating MASS->RareMetal.
  static bool GetWindows(
      const int window_size, const set<Position>& snp_positions,
      map<Position, Position>* windows);

  // DEPRECATED.
  // Reads alleles for each SNP position from RareMetalWorker Score files.
  //static bool GetAllelesFromRareMetalFiles(
  //    const vector<FileInfo>& files,
  //    map<Position, pair<string, string>>* allele_info);

  // DEPRECATED.
  // Reads alleles for each SNP position from MetalSKAT MInfo files.
  //static bool GetAllelesFromMetaSkatFiles(
  //    const vector<FileInfo>& files,
  //    map<Position, pair<string, string>>* allele_info);

};

}  // namespace premeta
#endif
