// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "read_raremetal_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
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
using file_reader_utils::ReadCsvInput;
using file_reader_utils::ReadCsvOutput;
using namespace map_utils;
using namespace math_utils;
using namespace string_utils;
using namespace std;

namespace premeta {

bool GetInfoFromRareMetalScoreFile(
    const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& file_info,
    int* num_samples, int* num_snps,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<Position, SnpInfo>* scores) {
  // Score File has a bunch of header info. Parse it.
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to open file."
         << endl;
    return false;
  }
  *num_samples = 0;
  string line;
  bool is_rv_test_format = false;
  while(getline(file, line)) {
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (HasPrefixString(line, kRareMetalScoreFileSampleNum)) {
      vector<string> num_samples_parts;
      Split(line, "=", &num_samples_parts);
      if (num_samples_parts.size() != 2) {
        cout << "ERROR in getting info from RareMetal Score file: Unable to parse "
             << "number of samples: '" << line << "'. Aborting.\n";
        file.close();
        return false;
      }
      string num_samples_str;
      RemoveAllWhitespace(num_samples_parts[1], &num_samples_str);
      if (!Stoi(num_samples_str, num_samples)) {
        cout << "ERROR in getting info from RareMetal Score file: Unable to parse "
             << "number of samples: '" << num_samples_str << "'. Aborting.\n";
        file.close();
        return false;
      }
    } else if (HasPrefixString(line, file_info.comment_char_)) {
    } else if (HasPrefixString(line, "CHROM")) {
      is_rv_test_format = true;
      break;
    } else {
      // Reached a non-comment line.
      break;
    }
  }
  if (*num_samples == 0) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to find number of "
         << "samples in file (there should be a line prefixed by "
         << kRareMetalScoreFileSampleNum << ". Aborting.\n";
    file.close();
    return false;
  }
  file.close();

  // Now read Score File, skipping metadata at the top.
  const int toggle = is_rv_test_format ? -1 : 0;
  ReadCsvInput input;
  input.filename_ = file_info.name_;
  input.delimiters_.insert(file_info.delimiter_);
  input.comment_char_ = file_info.comment_char_;
  input.has_header_ = is_rv_test_format;
  vector<pair<int, GenericDataType>>& column_types = input.columns_to_read_;
  column_types.push_back(make_pair(1, GenericDataType::STRING));  // Chromosome.
  column_types.push_back(make_pair(2, GenericDataType::UINT_64));  // Position.
  column_types.push_back(make_pair(3, GenericDataType::STRING));  // REF.
  column_types.push_back(make_pair(4, GenericDataType::STRING));  // ALT.
  column_types.push_back(make_pair(5, GenericDataType::INT));  // N_INFORMATIVE.
  column_types.push_back(make_pair(7 + toggle, GenericDataType::DOUBLE));  // ALL_AF (MAF).
  column_types.push_back(make_pair(11 + toggle, GenericDataType::INT));  // N_REF.
  column_types.push_back(make_pair(12 + toggle, GenericDataType::INT));  // N_HET.
  column_types.push_back(make_pair(13 + toggle, GenericDataType::INT));  // N_ALT.
  column_types.push_back(make_pair(14 + toggle, GenericDataType::DOUBLE));  // U_STAT.
  column_types.push_back(make_pair(15 + toggle, GenericDataType::DOUBLE));  // SQRT_V_STAT.
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to read score file "
         << "'" << file_info.name_ << "'. Aborted with error:\n" << output.error_msg_
         << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_score_file = output.output_;

  // Iterate through parsed_score_file, populating scores.
  *num_snps = parsed_score_file.size();
  for (int i = 0; i < parsed_score_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_score_file[i];

    if (current_row.size() != column_types.size()) {
      cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
           << "row " << i << " in " << file_info.name_
           << ": Unable to find all requisite columns (using tab delimiter). "
           << "Found " << current_row.size() << ", expected 9. Check that "
           << "score file has proper format." << endl;
      return false;
    }

    // Chromosome.
    Chromosome chr;
    if (!VcfUtils::ParseChromosome(current_row[0].str_, &chr)) {
      cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
           << "row " << i << " in " << file_info.name_
           << ": Unrecognized chromosome '" << current_row[0].str_
           << "'. Aborting.\n";
      return false;
    }

    // Position.
    uint64_t pos = current_row[1].uint64_;
    Position chr_pos;
    chr_pos.chr_ = chr;
    chr_pos.pos_ = pos;
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(chr_pos);
      if (allele_info_itr != allele_info.end()) {
        chr_pos = allele_info_itr->second.pos_;
      }
    }

    if (scores->find(chr_pos) != scores->end()) {
      cout << "ERROR in getting info from RareMetal Score file on row "
           << i << " in " << file_info.name_ << ": Already have seen "
           << "Chromosome " << current_row[0].str_ << " and position "
           << pos << " on a previous line. Aborting.\n";
      return false;
    }

    // Score Info.
    SnpInfo info;
    info.major_allele_ = current_row[2].str_;
    info.minor_allele_ = current_row[3].str_;
    info.n_info_ = current_row[4].int_;
    info.maf_ = current_row[5].dbl_;
    info.n_ref_ = current_row[6].int_;
    info.n_het_ = current_row[7].int_;
    info.n_alt_ = current_row[8].int_;
    info.mac_ = info.n_het_ + (2 * info.n_alt_);
    info.num_non_missing_ = info.n_ref_ + info.n_het_ + info.n_alt_;
    info.u_stat_ = current_row[9].dbl_;
    info.sqrt_v_stat_ = current_row[10].dbl_;
    if (info.sqrt_v_stat_ <= 0.0) {
      cout << "ERROR in getting info from RareMetal Score file: Invalid SQRT_V_STAT "
           << "value found on row " << i << " in " << file_info.name_
           << ": " << info.sqrt_v_stat_ << ". Aborting.\n";
      return false;
    }
    if (info.num_non_missing_ == 0 ||
        info.maf_ < 1.0 / (4.0 * info.num_non_missing_) ||
        info.maf_ > 1.0 - (1.0 / (4.0 * info.num_non_missing_))) {
      monomorphic_snps->insert(chr_pos);
    }

    // Check that REF/ALT alleles match existing entry for this SNP. Otherwise,
    // either add this SNP to snp_to_ref_alt_and_study_ or snp_to_excluding_study_,
    // as appropriate.
    tuple<string, string, set<int>>* snp_info =
        FindOrInsert(
            chr_pos, *snp_to_ref_alt_and_study,
            make_tuple(info.major_allele_, info.minor_allele_, set<int>()));
    if (get<0>(*snp_info) != info.major_allele_ ||
        get<1>(*snp_info) != info.minor_allele_) {
      if (get<1>(*snp_info) != info.major_allele_ ||
          get<0>(*snp_info) != info.minor_allele_) {
        // This is a more aggregious difference than just a swapping of
        // REF <-> ALT. Add snp to snp_to_excluding_study_.
        set<int>* excluding_studies =
            FindOrInsert(chr_pos, *snp_to_excluding_study, set<int>());
        excluding_studies->insert(study_num);
        continue;
      } else {
        get<2>(*snp_info).insert(study_num);
      }
    }

    // The other fields of SnpInfo require covariance data to complete;
    // this will be done after reading in covariance data.
    scores->insert(make_pair(chr_pos, info));
  }

  return true;
}

bool GetMonoAndMultiSnpsFromRareMetalScoreFile(
    const string& score_file,
    map<Position, int>* snps_in_score_file,
    set<Position>* monomorphic_snps,
    set<Position>* multi_allelic_snps,
    set<Position>* multi_allelic_snps_two,
    set<Position>* non_snps,
    set<Position>* non_std_chr_snps) {
  // Score File has a bunch of header info. Parse it.
  ifstream file(score_file.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to open file."
         << endl;
    return false;
  }
  string line;
  bool is_rv_test_format = false;
  while(getline(file, line)) {
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (HasPrefixString(line, kRareMetalScoreFileDefaultCommentChar)) {
    } else if (HasPrefixString(line, "CHROM")) {
      is_rv_test_format = true;
      break;
    } else {
      // Reached a non-comment line.
      break;
    }
  }
  file.close();

  // Now read Score File, skipping metadata at the top.
  const int toggle = is_rv_test_format ? -1 : 0;
  ReadCsvInput input;
  input.filename_ = score_file;
  input.delimiters_.insert(kRareMetalScoreFileDefaultDelimiter);
  input.comment_char_ = kRareMetalScoreFileDefaultCommentChar;
  input.has_header_ = is_rv_test_format;
  vector<pair<int, GenericDataType>>& column_types = input.columns_to_read_;
  column_types.push_back(make_pair(1, GenericDataType::STRING));  // Chromosome.
  column_types.push_back(make_pair(2, GenericDataType::UINT_64));  // Position.
  column_types.push_back(make_pair(3, GenericDataType::STRING));  // REF.
  column_types.push_back(make_pair(4, GenericDataType::STRING));  // ALT.
  column_types.push_back(make_pair(7 + toggle, GenericDataType::DOUBLE));  // ALL_AF (MAF).
  column_types.push_back(make_pair(11 + toggle, GenericDataType::INT));  // N_REF.
  column_types.push_back(make_pair(12 + toggle, GenericDataType::INT));  // N_HET.
  column_types.push_back(make_pair(13 + toggle, GenericDataType::INT));  // N_ALT.
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from RareMetal Score file: Unable to read score file "
         << "'" << score_file << "'. Aborted with error:\n" << output.error_msg_
         << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_score_file = output.output_;

  // Iterate through parsed_score_file, populating scores.
  set<Position> snps_encountered;
  for (int i = 0; i < parsed_score_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_score_file[i];

    if (current_row.size() != column_types.size()) {
      cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
           << "row " << i << " in " << score_file
           << ": Unable to find all requisite columns (using tab delimiter). "
           << "Found " << current_row.size() << ", expected 9. Check that "
           << "score file has proper format." << endl;
      return false;
    }

    // Chromosome.
    Chromosome chr;
    if (!VcfUtils::ParseChromosome(current_row[0].str_, &chr)) {
      // We allow unparsable chromosomes, if these are going to be removed anyway,
      // which happens iff non_std_chr_snps is non-null. So don't abort in this case.
      if (non_std_chr_snps == nullptr) {
        cout << "ERROR in getting info from RareMetal Score file: Unable to parse data "
             << "row " << i << " in " << score_file
             << ": Unrecognized chromosome '" << current_row[0].str_
             << "'. Aborting.\n";
        return false;
      }
    }

    // Position.
    uint64_t pos = current_row[1].uint64_;
    Position chr_pos;
    chr_pos.chr_ = chr;
    chr_pos.pos_ = pos;

    // Update counter for how many times this SNP has been encountered.
    if (snps_in_score_file != nullptr) {
      int* num_occurences = FindOrInsert(chr_pos, *snps_in_score_file, 0);
      (*num_occurences)++;
    }

    // If appropriate (non_std_chr_snps is non-null), track SNPs that aren't on
    // Chromosome 1-22 nor X.
    const int chr_num = static_cast<int>(chr);
    if (non_std_chr_snps != nullptr && (chr_num < 1 || chr_num > 23)) {
      non_std_chr_snps->insert(chr_pos);
    }

    // Set SNP as one that has been seen in Score file; if already present, add
    // SNP to multi-allelic SNPs.
    if (!snps_encountered.insert(chr_pos).second) {
      multi_allelic_snps->insert(chr_pos);
    }

    // Check for Non-SNPs and Multi-Allelic (type 2) SNPs.
    const string& ref = current_row[2].str_;
    if (ref.find(",") != string::npos) multi_allelic_snps_two->insert(chr_pos);
    else if (ref.length() != 1) non_snps->insert(chr_pos);
    const string& alt = current_row[3].str_;
    if (alt.find(",") != string::npos) multi_allelic_snps_two->insert(chr_pos);
    else if (alt.length() != 1) non_snps->insert(chr_pos);

    // Check for Monomorphic SNP:
    //   - NM := Number Non-Missing SNPs = NHET + NALT + NREF
    //   - Monomorphic if: NM = 0, or MAF < 1/4NM or MAF > 1 - 1/4NM
    SnpInfo info;
    info.maf_ = current_row[4].dbl_;
    info.n_ref_ = current_row[5].int_;
    info.n_het_ = current_row[6].int_;
    info.n_alt_ = current_row[7].int_;
    info.num_non_missing_ = info.n_ref_ + info.n_het_ + info.n_alt_;
    if (info.num_non_missing_ == 0 ||
        info.maf_ < 1.0 / (4.0 * info.num_non_missing_) ||
        info.maf_ > 1.0 - (1.0 / (4.0 * info.num_non_missing_))) {
      monomorphic_snps->insert(chr_pos);
    }
  }

  return true;
}

bool GetInfoFromRareMetalGroupFile(
    const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const FileInfo& file_info,
    set<Position>* snps_without_score_info,
    map<string, vector<Position>>* gene_to_snp_positions) {
  ReadCsvInput input;
  input.filename_ = file_info.name_;
  input.delimiters_.insert(file_info.delimiter_);
  input.comment_char_ = file_info.comment_char_;
  input.has_header_ = false;
  // Read all columns (format [0, 1, 0] encodes range [1, \end]; see csv_utils.h)
  vector<int> group_columns;
  group_columns.push_back(0);
  group_columns.push_back(1);
  group_columns.push_back(0);
  vector<pair<vector<int>, GenericDataType>>& column_types = input.range_columns_to_read_;
  column_types.push_back(make_pair(group_columns, GenericDataType::STRING));
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from RareMetal Group file: Unable to read group "
         << "file '" << file_info.name_ << "'. Aborted with error:\n"
         << output.error_msg_ << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_snps_by_gene_file = output.output_;

  // Iterate through parsed_snps_by_gene_file, populating gene_to_snp_positions.
  for (int i = 0; i < parsed_snps_by_gene_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_snps_by_gene_file[i];

    // Fetch gene name (first entry on line).
    const string& gene_name = current_row[0].str_;
    if (gene_to_snp_positions->find(gene_name) != gene_to_snp_positions->end()) {
      cout << "ERROR in getting info from RareMetal Group file: Unable to parse line "
           << (i + 1) << ". Gene name '" << gene_name << "' has already appeared "
           << "earlier in the file. Aborting\n";
      return false;
    }

    if (current_row.size() < 2) {
      cout << "ERROR in getting info from RareMetal Group file: Unable to parse line "
           << (i + 1) << ". Only found one column '" << gene_name
           << "' on this line (using delimiter '" << file_info.delimiter_
           << "'). Aborting.\n";
      return false;
    }

    map<string, vector<Position>>::iterator itr =
        gene_to_snp_positions->insert(
            make_pair(gene_name, vector<Position>())).first;
    vector<Position>& snps = itr->second;

    // Go through rest of line, extracting SNP position info.
    for (int j = 1; j < current_row.size(); ++j) {
      const string& pos_info = current_row[j].str_;
      vector<string> pos_parts;
      Split(pos_info, ":", &pos_parts);
      if (pos_parts.size() != 4) {
        cout << "ERROR in getting info from RareMetal Group file: Unable to parse SNP "
             << j << " on line " << (i + 1) << ". Expected three colon "
             << "dividers, found " << (pos_parts.size() - 1) << endl;
        return false;
      }

      // Parse Chromosome.
      const string& chr_str = pos_parts[0];
      Chromosome chr;
      if (!VcfUtils::ParseChromosome(chr_str, &chr)) {
        cout << "ERROR in getting info from RareMetal Group file: Unable to parse SNP "
             << j << " on line " << (i + 1) << ". Unrecognized chromosome '"
             << chr_str << "'. Aborting.\n";
        return false;
      }

      // Parse Position.
      uint64_t pos;
      if (!Stoi(pos_parts[1], &pos)) {
        cout << "ERROR in getting info from RareMetal Group file: Unable to parse SNP "
             << j << " on line " << (i + 1) << ". Unable to parse position '"
             << pos_parts[1] << "' as an integer. Aborting.\n";
        return false;
      }

      // Add Position and Chromosome.
      Position chr_pos;
      chr_pos.chr_ = chr;
      chr_pos.pos_ = pos;

      // Use allele_info to get format of this SNP in the target software.
      if (!allele_info.empty()) {
        map<Position, SnpInfo>::const_iterator allele_info_itr =
            allele_info.find(chr_pos);
        if (allele_info_itr != allele_info.end()) {
          chr_pos = allele_info_itr->second.pos_;
        }
      }

      // Check that this SNP shouldn't be skipped (due to mismatched REF/ALT
      // alleles).
      if (ShouldExcludeSnp(study_num, chr_pos, snp_to_excluding_study)) {
          continue;
      }

      // Make sure this SNP appeared in the .score file. If not, keep track of
      // it separately (via snps_without_score_info), but don't add it to
      // gene_to_snp_positions.
      if (snp_info.find(chr_pos) == snp_info.end()) {
        snps_without_score_info->insert(chr_pos);
        continue;
      }

      snps.push_back(chr_pos);
    }

    // Remove this gene if it has no SNPs.
    if (snps.empty()) {
      gene_to_snp_positions->erase(gene_name);
    }
  }

  return true;
}

bool GetInfoFromRareMetalCovarianceFile(
    const bool is_new_version, const bool keep_snps_without_score_info,
    const int study_num, const int num_samples, const FileInfo& file_info,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, SnpInfo>& snp_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    set<Position>* snps_without_score_info,
    map<Position, map<Position, double>>* snp_to_cov_w_neighbors) {
  // Look at the Metadata at the top of the .cov.txt file to see if the format
  // is RAREMETALWORKER or RVTEST (the latter doesn't have a comment char at
  // the start of the Header line).
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting info from RareMetal Covariance file: Unable to open file."
         << endl;
    return false;
  }
  string line;
  bool is_rv_test_format = false;
  while(getline(file, line)) {
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (HasPrefixString(line, file_info.comment_char_)) {
    } else if (HasPrefixString(line, "CHROM")) {
      is_rv_test_format = true;
      break;
    } else {
      // Reached a non-comment line.
      break;
    }
  }

  const int toggle = is_rv_test_format ? 2 : 0;
  ReadCsvInput input;
  input.filename_ = file_info.name_;
  input.delimiters_.insert(file_info.delimiter_);
  input.comment_char_ = file_info.comment_char_;
  input.has_header_ = is_rv_test_format;
  vector<pair<int, GenericDataType>>& column_types = input.columns_to_read_;
  column_types.push_back(make_pair(1, GenericDataType::STRING));  // Chromosome.
  column_types.push_back(make_pair(2, GenericDataType::UINT_64));  // Position.
  column_types.push_back(make_pair(3 + toggle, GenericDataType::STRING));  // MARKERS_IN_WINDOW.
  column_types.push_back(make_pair(4 + toggle, GenericDataType::STRING));  // COV_MATRICES.
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from RareMetal Covariance file: Unable to read group "
         << "file '" << file_info.name_ << "'. Aborted with error:\n"
         << output.error_msg_ << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_cov_file = output.output_;

  // Iterate through parsed_cov_file, populating snp_to_cov_w_neighbors.
  for (int i = 0; i < parsed_cov_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_cov_file[i];

    // Chromosome.
    Chromosome chr;
    if (!VcfUtils::ParseChromosome(current_row[0].str_, &chr)) {
      cout << "ERROR in getting info from RareMetal Covariance file: Unable to parse "
           << "data row " << (i + 1) << " in " << file_info.name_
           << ": Unrecognized chromosome '" << current_row[0].str_
           << "'. Aborting.\n";
      return false;
    }

    // Position.
    uint64_t pos = current_row[1].uint64_;
    Position chr_pos;
    chr_pos.chr_ = chr;
    chr_pos.pos_ = pos;
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(chr_pos);
      if (allele_info_itr != allele_info.end()) {
        chr_pos = allele_info_itr->second.pos_;
      }
    }

    // Check that this SNP shouldn't be skipped (due to mismatched REF/ALT
    // alleles).
    if (ShouldExcludeSnp(study_num, chr_pos, snp_to_excluding_study)) {
      continue;
    }

    // Make sure this SNP appeared in the .score file. If not, keep track of
    // it separately (via snps_without_score_info), but don't add it to
    // snp_to_cov_w_neighbors.
    if (snp_info.find(chr_pos) == snp_info.end()) {
      snps_without_score_info->insert(chr_pos);
      if (!keep_snps_without_score_info) {
        continue;
      }
    }

    if (snp_to_cov_w_neighbors->find(chr_pos) != snp_to_cov_w_neighbors->end()) {
      cout << "ERROR in getting info from RareMetal Covariance file on row "
           << (i + 1) << " in " << file_info.name_ << ": Already have seen "
           << "Chromosome " << current_row[0].str_ << " and position "
           << pos << " on a previous line. Aborting.\n";
      return false;
    }

    map<Position, map<Position, double>>::iterator itr =
        snp_to_cov_w_neighbors->insert(make_pair(
            chr_pos, map<Position, double>())).first;
    map<Position, double>& neighbor_snps_to_cov = itr->second;

    // Markers in Window.
    vector<string> markers;
    Split(current_row[2].str_, ",", &markers);

    // Covariances.
    vector<string> covariances;
    Split(current_row[3].str_, ",", &covariances);

    if (markers.size() != covariances.size()) {
      cout << "ERROR in getting info from RareMetal Covariance file on row "
           << (i + 1) << " in " << file_info.name_ << ": Number of markers ("
           << markers.size() << ") does not match the number of Covariance "
           << "Matrices (" << covariances.size() << "). Aborting.\n";
      return false;
    }

    for (int j = 0; j < markers.size(); ++j) {
      uint64_t position;
      if (!Stoi(markers[j], &position)) {
        cout << "ERROR in getting info from RareMetal Covariance file: Unable to parse "
             << "Marker " << (j + 1) << " on line " << (i + 1) << " as an integer position"
             << ": " << markers[j] << ". Aborting.\n";
        return false;
      }
      double cov;
      if (!Stod(covariances[j], &cov)) {
        cout << "ERROR in getting info from RareMetal Covariance file: Unable to parse "
             << "Covariance " << (j + 1) << " on line " << (i + 1) << " as a double value"
             << ": " << covariances[j] << ". Aborting.\n";
        return false;
      }
      //if (is_new_version) {
      //  cov *= static_cast<double>(num_samples);
      //}
      Position neighbor_pos;
      neighbor_pos.chr_ = chr;
      neighbor_pos.pos_ = position;
      // Use allele_info to get format of this SNP in the target software.
      if (!allele_info.empty()) {
        map<Position, SnpInfo>::const_iterator allele_info_itr =
            allele_info.find(neighbor_pos);
        if (allele_info_itr != allele_info.end()) {
          neighbor_pos = allele_info_itr->second.pos_;
        }
      }

      // Check that this SNP shouldn't be skipped (due to mismatched REF/ALT
      // alleles).
      if (ShouldExcludeSnp(study_num, neighbor_pos, snp_to_excluding_study)) {
        continue;
      }

      // Make sure this SNP appeared in the .score file. If not, keep track of
      // it separately (via snps_without_score_info), but don't add it to
      // snp_to_cov_w_neighbors.
      if (snp_info.find(neighbor_pos) == snp_info.end()) {
        snps_without_score_info->insert(neighbor_pos);
        if (!keep_snps_without_score_info) {
          continue;
        }
      }

      const bool new_neighbor =
          neighbor_snps_to_cov.insert(make_pair(neighbor_pos, cov)).second;
      if (!new_neighbor) {
        cout << "ERROR in getting info from RareMetal Covariance file: Neighbor on chr "
             << chr << " with position " << position << " was found as the "
             << (j + 1) << "th item on line " << (i + 1) << ", but it also appeared "
             << "earlier. Aborting.\n";
        return false;
      }
    }
  }

  return true;
}

bool GetMultiAllelicSnpsFromRareMetalCovarianceFile(
    const string& cov_file,
    map<Position, int>* snps_in_cov_file,
    set<Position>* multi_allelic_snps, set<Position>* non_std_chr_snps) {
  // Look at the Metadata at the top of the .cov.txt file to see if the format
  // is RAREMETALWORKER or RVTEST (the latter doesn't have a comment char at
  // the start of the Header line).
  ifstream file(cov_file.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting info from RareMetal Covariance file: Unable to open file."
         << endl;
    return false;
  }
  string line;
  bool is_rv_test_format = false;
  while(getline(file, line)) {
    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (HasPrefixString(line, "#")) {
    } else if (HasPrefixString(line, "CHROM")) {
      is_rv_test_format = true;
      break;
    } else {
      // Reached a non-comment line.
      break;
    }
  }

  const int toggle = is_rv_test_format ? 2 : 0;
  ReadCsvInput input;
  input.filename_ = cov_file;
  input.delimiters_.insert(kRareMetalCovFileDefaultDelimiter);
  input.comment_char_ = kRareMetalCovFileDefaultCommentChar;
  input.has_header_ = is_rv_test_format;
  vector<pair<int, GenericDataType>>& column_types = input.columns_to_read_;
  column_types.push_back(make_pair(1, GenericDataType::STRING));  // Chromosome.
  column_types.push_back(make_pair(2, GenericDataType::UINT_64));  // Position.
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from RareMetal Covariance file: Unable to read group "
         << "file '" << cov_file << "'. Aborted with error:\n"
         << output.error_msg_ << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_cov_file = output.output_;

  // Iterate through parsed_cov_file, looking for multi-allelic SNPs.
  set<Position> encountered_snps;
  for (int i = 0; i < parsed_cov_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_cov_file[i];

    // Chromosome.
    Chromosome chr;
    if (!VcfUtils::ParseChromosome(current_row[0].str_, &chr)) {
      // We allow unparsable chromosomes, if these are going to be removed anyway,
      // which happens iff non_std_chr_snps is non-null. So don't abort in this case.
      if (non_std_chr_snps == nullptr) {
        cout << "ERROR in getting info from RareMetal Covariance file: Unable to parse "
             << "data row " << (i + 1) << " in " << cov_file
             << ": Unrecognized chromosome '" << current_row[0].str_
             << "'. Aborting.\n";
        return false;
      }
    }

    // Position.
    uint64_t pos = current_row[1].uint64_;
    Position chr_pos;
    chr_pos.chr_ = chr;
    chr_pos.pos_ = pos;

    // Update counter for how many times this SNP has been encountered.
    if (snps_in_cov_file != nullptr) {
      int* num_occurences = FindOrInsert(chr_pos, *snps_in_cov_file, 0);
      (*num_occurences)++;
    }

    // If appropriate (non_std_chr_snps is non-null), track SNPs that aren't on
    // Chromosome 1-22 nor X.
    const int chr_num = static_cast<int>(chr);
    if (non_std_chr_snps != nullptr && (chr_num < 1 || chr_num > 23)) {
      non_std_chr_snps->insert(chr_pos);
    }

    // Add this SNP to the set of already encountered SNPs, and if it was already
    // there, add it to multi_allelic_snps.
    if (!encountered_snps.insert(chr_pos).second) {
      multi_allelic_snps->insert(chr_pos);
    }
  }

  return true;
}

bool GetSigmaFromRareMetalFiles(
    const int num_samples,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    const map<Position, SnpInfo>& scores,
    double* sigma) {
  *sigma = -1.0;
  for (const auto& score : scores) {
    const Position& pos = score.first;
    const SnpInfo& pos_info = score.second;
    map<Position, map<Position, double>>::const_iterator cov_itr =
        snp_to_cov_w_neighbors.find(pos);
    if (cov_itr == snp_to_cov_w_neighbors.end()) continue;

    // snp_to_cov_w_neighbors should contain self-covariance.
    map<Position, double>::const_iterator self_cov_itr =
        cov_itr->second.find(pos);
    if (self_cov_itr ==  cov_itr->second.end()) continue;

    // Only need to compute sigma once, as it is the same for all positions.
    *sigma = pos_info.sqrt_v_stat_ / sqrt(num_samples * self_cov_itr->second);
    return true;
  }

  return false;
}

void RescaleCovariances(const double& rescale_amount,
                        map<pair<Position, Position>, double>* covariances) {
  for (map<pair<Position, Position>, double>::iterator itr = covariances->begin();
       itr != covariances->end(); ++itr) {
    itr->second *= rescale_amount;
  }
}

void RescaleSnpInfo(
    const double& u_stat_rescale, const double& sqrt_v_stat_rescale,
    map<Position, SnpInfo>* scores) {
  // Nothing to do if all rescale constants are 1.
  if (FloatEq(u_stat_rescale, 1.0) && FloatEq(sqrt_v_stat_rescale, 1.0)) {
    return;
  }

  for (map<Position, SnpInfo>::iterator itr = scores->begin();
       itr != scores->end(); ++itr) {
    SnpInfo& info = itr->second;
    info.u_stat_ *= u_stat_rescale;
    info.sqrt_v_stat_ *= sqrt_v_stat_rescale;
  }
}

bool ComputeCovarianceByGene(
    const int num_samples, const double& sigma_sq,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, map<Position, double>>& snp_to_cov_w_neighbors,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances) {
  for (const auto& gene_itr : gene_to_snp_positions) {
    const string& gene = gene_itr.first;
    const vector<Position>& positions = gene_itr.second;
    for (int i = 0; i < positions.size(); ++i) {
      const Position& row_pos = positions[i];
      map<Position, map<Position, double>>::const_iterator row_cov_itr =
          snp_to_cov_w_neighbors.find(row_pos);
      for (int j = i; j < positions.size(); ++j) {
        const Position& col_pos = positions[j];
        // First check to see if we already have a covariance for these
        // two SNP Positions.
        double covariance = DBL_MIN;
        if (LookupCovariance(row_pos, col_pos, *covariances, &covariance)) {
          continue;
        }
        map<Position, map<Position, double>>::const_iterator col_cov_itr =
            snp_to_cov_w_neighbors.find(col_pos);
        if (row_cov_itr == snp_to_cov_w_neighbors.end() &&
            col_cov_itr == snp_to_cov_w_neighbors.end()) {
          if (row_pos.chr_ == col_pos.chr_ && row_pos.pos_ == col_pos.pos_) {
            if (i != j) {
              cout << "ERROR in computing Covariance by gene: For gene '"
                   << gene << "', RareMetal's Gene Grouping file lists Position "
                   << PrintPosition(row_pos) << " multiple times. Aborting."
                   << endl;
              return false;
            }
            snps_to_skip->insert(row_pos);
            continue;
          }
          // No Covariance for this pair of Positions; enter 0.
          covariances->insert(make_pair(make_pair(row_pos, col_pos), 0.0));
          continue;
        }
        if (row_cov_itr != snp_to_cov_w_neighbors.end()) {
          const map<Position, double>& row_pos_covariances = row_cov_itr->second;
          map<Position, double>::const_iterator value_itr =
              row_pos_covariances.find(col_pos);
          if (value_itr != row_pos_covariances.end()) {
            covariances->insert(make_pair(
                  make_pair(row_pos, col_pos),
                  (value_itr->second * num_samples * sigma_sq)));
            continue;
          }
        }
        if (col_cov_itr != snp_to_cov_w_neighbors.end()) {
          const map<Position, double>& col_pos_covariances = col_cov_itr->second;
          map<Position, double>::const_iterator value_itr =
              col_pos_covariances.find(row_pos);
          if (value_itr != col_pos_covariances.end()) {
            covariances->insert(make_pair(
                  make_pair(row_pos, col_pos),
                  (value_itr->second * num_samples * sigma_sq)));
            continue;
          }
        }
        // If we reached this point, we failed to find the covariance; enter 0.
        if (row_pos.chr_ == col_pos.chr_ && row_pos.pos_ == col_pos.pos_) {
          if (i != j) {
            cout << "ERROR in computing Covariance by gene: For gene '"
                 << gene << "', RareMetal's Gene Grouping file lists Position "
                 << PrintPosition(row_pos) << " multiple times. Aborting."
                 << endl;
            return false;
          }
          snps_to_skip->insert(row_pos);
          continue;
        }
        covariances->insert(make_pair(make_pair(row_pos, col_pos), 0.0));
      }
    }
  }
  return true;
}

bool SanityCheckSnpsToSkip(
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, SnpInfo>& scores,
    const map<pair<Position, Position>, double>& covariances) {
  if (snps_to_skip.empty()) return true;
  for (const auto& scores_itr : scores) {
    if (monomorphic_snps.find(scores_itr.first) != monomorphic_snps.end()) {
      // No need to keep sanity-checking this SNP, since it is monomorphic.
      continue;
    }
    if (snps_to_skip.find(scores_itr.first) != snps_to_skip.end()) {
      cout << "ERROR: Found SNP "
           << PrintPosition(scores_itr.first)
           << ", but it appears in RareMetal Score file." << endl;
      return false;
    }
  }

  for (const auto& cov_itr : covariances) {
    if (snps_to_skip.find(cov_itr.first.first) != snps_to_skip.end() &&
        monomorphic_snps.find(cov_itr.first.first) == monomorphic_snps.end()) {
      cout << "ERROR: Found SNP "
           << PrintPosition(cov_itr.first.first)
           << " among the SNPs to skip, but it appears "
           << "in RareMetal Covariance file." << endl;
      return false;
    }
    if (snps_to_skip.find(cov_itr.first.second) != snps_to_skip.end() &&
        monomorphic_snps.find(cov_itr.first.second) == monomorphic_snps.end()) {
      cout << "ERROR: Found snp "
           << PrintPosition(cov_itr.first.second)
           << " among the SNPs to skip, but it appears "
           << "in RareMetal Covariance file." << endl;
      return false;
    }
  }

  return true;
}

}  // namespace premeta
