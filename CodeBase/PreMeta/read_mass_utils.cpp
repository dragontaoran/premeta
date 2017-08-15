// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "read_mass_utils.h"

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

bool UpdateGeneCovariancesFromMassFile(
    const vector<Position>& gene_snp_positions,
    const map<Position, vector<double>>& gene_covariances,
    map<pair<Position, Position>, double>* covariances) {
  for (int i = 0; i < gene_snp_positions.size(); ++i) {
    const Position& current_snp = gene_snp_positions[i];
    map<Position, vector<double>>::const_iterator cov_itr =
        gene_covariances.find(current_snp);
    if (cov_itr == gene_covariances.end()) {
      cout << "ERROR in reading covariance info from MASS file: "
           << "Unable to find covariance information for SNP at position "
           << PrintPosition(current_snp) << ". Aborting." << endl;
      return false;
    }
    const vector<double>& row_covariances = cov_itr->second;
    if (row_covariances.size() != gene_snp_positions.size()) {
      cout << "ERROR in reading covariance info from MASS file: Found "
           << gene_snp_positions.size() << " snp(s) for this gene, but row for "
           << "position " << PrintPosition(current_snp) 
           << " has " << row_covariances.size() << " covariances listed. "
           << "Aborting." << endl;
      return false;
    }
    // NOTE: The below line should be used if MASS's first step (ScoreSeq)
    // outputs its covariance matrices in UPPER-TRIANGULAR format. But as of
    // July 2015, ScoreSeq outputs covariance matrices in LOWER-TRIANGULAR
    // format, so we extract those positions from 'gene_covariances'.
    //for (int j = i; j < row_covariances.size(); ++j) {
    for (int j = 0; j <= i; ++j) {
      pair<map<pair<Position, Position>, double>::const_iterator, bool> in_itr =
          covariances->insert(make_pair(
              make_pair(current_snp, gene_snp_positions[j]),
              row_covariances[j]));
      if (!in_itr.second &&
          !FloatEq(in_itr.first->second, row_covariances[j])) {
        cout << "ERROR in reading covariance info from MASS file: For Positions "
             << PrintPosition(current_snp) << " and "
             << PrintPosition(gene_snp_positions[j])
             << ", multiple (and non-matching) covariances were found: "
             << in_itr.first->second << " and " << row_covariances[j]
             << ". Aborting." << endl;
        return false;
      }
    }
  }

  return true;
}

bool GetNumSamplesFromMassFile(
    const FileInfo& mass_file_info, int* num_samples) {
  ifstream file(mass_file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR: Unable to open MASS input file '"
         << mass_file_info.name_ << "'." << endl;
    file.close();
    return false;
  }
  string line = "";
  if (!getline(file, line)) {
    cout << "ERROR in reading first line of MASS file '"
         << mass_file_info.name_ << "'. Aborting." << endl;
    file.close();
    return false;
  }

  // Remove trailing character '13' (Carriage Return) at end of line, if present.
  if (line[line.length() - 1] == 13) {
    line = line.substr(0, line.length() - 1);
  }

  if (!HasPrefixString(line, mass_file_info.comment_char_)) {
    cout << "ERROR in reading first line ('#Samples') of MASS file: "
         << "Expected num Samples on first line of mass file '"
         << mass_file_info.name_
         << "', but first line (" << line << ") is not a comment. Aborting."
         << endl;
    file.close();
    return false;
  }
  if (line.find("Samples") == string::npos) {
    cout << "ERROR in reading first line ('#Samples') of MASS file: "
         << "Expected num samples on first line of mass file '"
         << mass_file_info.name_
         << "', but found on top line '" << line 
         << "', which does not specify " << "Samples" << ". Aborting."
         << endl;
    file.close();
    return false;
  }
  vector<string> num_samples_parts;
  Split(line, "=", &num_samples_parts);
  if (num_samples_parts.size() != 2) {
    cout << "ERROR in reading first line ('#Samples') of MASS file: "
         << "Expected num samples on first line of mass file '"
         << mass_file_info.name_
         << "', but found on first line '" << line 
         << "', which does not have the expected equals sign. Aborting."
         << endl;
    file.close();
    return false;
  }
  string num_samples_str;
  RemoveAllWhitespace(num_samples_parts[1], &num_samples_str);
  if (!Stoi(num_samples_str, num_samples)) {
    cout << "ERROR in reading first line ('#Samples') of MASS file: "
         << "Expected num samples on first line of mass file '"
         << mass_file_info.name_
         << "', but found on first line '" << line 
         << "', which cannot be parsed as an integer. Aborting."
         << endl;
    file.close();
    return false;
  }
  file.close();
  return true;
}

bool GetInfoFromMassFile(
    const bool require_snp_in_allele_file, const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    int* num_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<Position, SnpInfo>* snp_info,
    map<pair<Position, Position>, double>* covariances) {
  ReadCsvInput input;
  input.filename_ = mass_file_info.name_;
  input.delimiters_.insert(mass_file_info.delimiter_);
  input.comment_char_ = mass_file_info.comment_char_;
  input.has_header_ = false;
  vector<int> string_cols;
  string_cols.push_back(1);  // Gene Name.
  string_cols.push_back(2);  // Chr:Pos.
  vector<int> int_cols;
  int_cols.push_back(4);     // MAC
  int_cols.push_back(5);     // Num Non-Missing Samples.
  int_cols.push_back(6);     // N_REF.
  int_cols.push_back(7);     // N_HET.
  int_cols.push_back(8);     // N_ALT.
  // Vector [3, 0, 9, 0] indicates columns {3, 9 - END} are double; see csv_utils.h.
  vector<int> double_cols;
  double_cols.push_back(3);  // MAF
  double_cols.push_back(0);
  double_cols.push_back(9);  // Column 9 is U-Stat; columns 10-END are covariances.
  double_cols.push_back(0);
  vector<pair<vector<int>, GenericDataType>>& column_types = input.range_columns_to_read_;
  column_types.push_back(make_pair(string_cols, GenericDataType::STRING));
  column_types.push_back(make_pair(int_cols, GenericDataType::INT));
  column_types.push_back(make_pair(double_cols, GenericDataType::DOUBLE));
  ReadCsvOutput output;
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in getting info from Mass file: Unable to read MASS "
         << "file '" << mass_file_info.name_ << "'. Aborted with error:\n"
         << output.error_msg_ << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_cov_file = output.output_;

  *num_snps = parsed_cov_file.size();

  // Iterate through the lines of the Mass file, populating score and
  // covariance information.
  string current_gene = "";
  map<Position, vector<double>> gene_covariances;
  vector<Position>* gene_snp_positions;
  for (int i = 0; i < parsed_cov_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_cov_file[i];

    // Gene Name.
    const string& gene_name = current_row[0].str_;

    // Determine if we're processing a new gene, or the same as before.
    if (gene_name != current_gene) {
      // If it's a new gene. Merge covariances from previous gene.
      if (i != 0 &&
          !UpdateGeneCovariancesFromMassFile(
              *gene_snp_positions, gene_covariances, covariances)) {
        cout << "ERROR in getting info from Mass file: "
             << "Unable to combine covariances for gene '" << current_gene
             << "'. See row " << (i + 2) << endl;
        return false;
      }
      pair<map<string, vector<Position>>::iterator, bool> gene_insertion_itr =
          gene_to_snp_positions->insert(make_pair(gene_name, vector<Position>()));
      if (!gene_insertion_itr.second) {
        cout << "WARNING while getting info from Mass file: Gene' " << gene_name
             << "' appeared in a non-contiguous block in the input file. "
             << "This may lead to an error." << endl;
      }
      gene_snp_positions = &(gene_insertion_itr.first->second);
      gene_covariances.clear();
      current_gene = gene_name;
    }
    
    // Chromosome and Position.
    Position snp_pos;
    if (!ParsePosition(current_row[1].str_, &snp_pos)) {
      // Unable to parse Position into CHR:POS format; store the Position string
      // into the snp_id_ field.
      snp_pos.snp_id_ = current_row[1].str_;
    }

    // Make sure this SNP is not among those that should be skipped.
    if (ShouldExcludeSnp(study_num, snp_pos, snp_to_excluding_study)) {
      (*num_snps)--;
      continue;
    }

    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(snp_pos);
      if (allele_info_itr == allele_info.end()) {
        if (require_snp_in_allele_file) {
          // Depending on target software, we may or may not demand that the
          // SNP is present in the allele_file: If target software is MASS or
          // SeqMeta, then we don't need to know Major/Minor alleles; but if
          // target software is MetaSkat or RareMetal, then we need this info.
          snps_to_skip->insert(snp_pos);
          (*num_snps)--;
          continue;
        }
      } else {
        snp_pos = allele_info_itr->second.pos_;
      }
    } else if (require_snp_in_allele_file) {
      cout << "ERROR in getting info from Mass file: SNP_INFO file "
           << "not specified (or empty), so Major/Minor alleles are "
           << "not available (target software requires this information. "
           << "Aborting." << endl;
      return false;
    }

    gene_snp_positions->push_back(snp_pos);
    pair<map<Position, SnpInfo>::iterator, bool> insertion_itr =
        snp_info->insert(make_pair(snp_pos, SnpInfo()));
    SnpInfo& info = insertion_itr.first->second;
    // If this is first time seeing this SNP, add info for it.
    if (insertion_itr.second) {
      info.mac_ = current_row[2].int_;
      info.num_non_missing_ = current_row[3].int_;
      info.n_ref_ = current_row[4].int_;
      info.n_het_ = current_row[5].int_;
      info.n_alt_ = current_row[6].int_;
      info.maf_ = current_row[7].dbl_;
      info.u_stat_ = current_row[8].dbl_;
      if (info.num_non_missing_ == 0 ||
          info.maf_ < 1.0 / (4.0 * info.num_non_missing_) ||
          info.maf_ > 1.0 - (1.0 / (4.0 * info.num_non_missing_))) {
        monomorphic_snps->insert(snp_pos);
      }
    } else {
      // This SNP has already been seen before (in another gene). Sanity check
      // the current info matches the old info.
      const SnpInfo& info = insertion_itr.first->second;
      if (info.mac_ != current_row[2].int_ ||
          info.num_non_missing_ != current_row[3].int_ ||
          info.maf_ != current_row[7].dbl_ ||
          info.u_stat_ != current_row[8].dbl_ ||
          info.n_ref_ != current_row[4].int_ ||
          info.n_het_ != current_row[5].int_ ||
          info.n_alt_ != current_row[6].int_) {
        cout << "ERROR in getting info from Mass file on row "
             << (i + 2) << " in " << mass_file_info.name_ << ": Already have seen "
             << "SNP " << PrintPosition(snp_pos) 
             << " on a previous line, and the information about this "
             << "SNP is not consistent on these lines. Aborting.\n";
        return false;
      }
      /* Note: Having a Position (SNP) appear in multiple genes is ok; no
       * need to print out warning.
      cout << "WARNING in GetInfoFromMassFile on row "
           << (i + 2) << " in " << mass_file_info.name_ << ": Already have seen "
           << "Chromosome " << pos_parts[0] << " and position "
           << pos << " on a previous line.\n";
      */
    }

    // Covariances.
    map<Position, vector<double>>::iterator cov_itr =
        gene_covariances.insert(make_pair(snp_pos, vector<double>())).first;
    vector<double>& row_covariances = cov_itr->second;
    for (int j = 9; j < current_row.size(); ++j) {
      row_covariances.push_back(current_row[j].dbl_);
    }
  }

  // Finish processing the final gene.
  if (!UpdateGeneCovariancesFromMassFile(
          *gene_snp_positions, gene_covariances, covariances)) {
    cout << "ERROR in getting info from Mass file: Unable to combine "
         << "covariances for gene '" << current_gene << "'. See row "
         << parsed_cov_file.size() << endl;
    return false;
  }

  return true;
}

}  // namespace premeta
