// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "premeta_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_mass_utils.h"
#include "PreMeta/read_metaskat_utils.h"
#include "PreMeta/read_r_binary_utils.h"
#include "PreMeta/read_raremetal_utils.h"
#include "PreMeta/read_seqmeta_utils.h"
#include "PreMeta/write_mass_utils.h"
#include "PreMeta/write_metaskat_utils.h"
#include "PreMeta/write_raremetal_utils.h"
#include "PreMeta/write_r_binary_utils.h"
#include "PreMeta/write_seqmeta_utils.h"
#include "StringUtils/string_utils.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
//#include <rpc/xdr.h>
#include <set>
#include <string>
#include <vector>

// NOTE: The below commented out code works for calling methods in the
// R Internals library, but only allows calling methods in RApiSerialize.h,
// R.h, Rmath.h, and Rinternals.h. In particular, access to serialize.c is
// still blocked, and I'm not sure how to access it (it may be intentionally
// impossible to do so, see e.g.:
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-October/006739.html
//#include "RApiSerialize.h"
//#include <R.h>
//#include <Rmath.h>
//#include <Rinternals.h>
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

map<Position, tuple<string, string, set<int>>> PreMeta::snp_to_ref_alt_and_study_;
map<Position, set<int>> PreMeta::snp_to_excluding_study_;
bool PreMeta::check_snp_ref_alt_order_;

namespace {
// There are unfortunately a few valid keywords for the columns of the
// Allele (SNP_INFO) File:
//  (i)   Columns for Major and Minor Allele are optional (they are
//        mandatory if going from (MASS, SeqMeta) -> (MetaSkat, RM),
//        but since this file can also be used simply to translate
//        from one SNP_ID format to another, we don't demand
//        the Major/Minor Allele columns at this point)
//  (ii)  Column for the SNP_ID format for the target software is
//        mandatory; however, we allow header for this column to be
//        either kAlleleFilePositionHeader or kAlleleFileMetaSnpIdHeader
//  (iii) Column for the SNP_ID format for the study software is
//        optional (used if we want to translate SNP_ID format from
//        that used in Study to that required by target software). If
//        present, this column should have header
//        kAlleleFileStudySnpIdHeader
// To account for all of these possibilities, we scan the header in
// the following order:
//   1)  Assume SNP_INFO is used for *both* translating SNP_IDs
//       as well as identifying Major/Minor alleles; so all four
//       columns should be present
//         1a) Assume that (ii) uses the kAlleleFilePositionHeader keyword
//         1b) Assume that (ii) uses the kAlleleFileMetaSnpIdHeader keyword
//   2) Assume SNP_INFO is used only for translating SNP_IDs; so two
//      columns should be present
//         2a) Assume that (ii) uses the kAlleleFilePositionHeader keyword
//         2b) Assume that (ii) uses the kAlleleFileMetaSnpIdHeader keyword
//   3) Assume SNP_INFO is used only for identifying Major/Minor alleles; so
//      three columns should be present
//         3a) Assume that (ii) uses the kAlleleFilePositionHeader keyword
//         3b) Assume that (ii) uses the kAlleleFileMetaSnpIdHeader keyword
bool ReadAlleleFileHeader(
    const bool is_study_allele_file, const bool allow_format_two,
    const FileInfo& file_info, const string& col_one, const string& col_two,
    int* col_one_index, int* col_two_index, int* major_col, int* minor_col,
    int* max_column_index) {
  map<string, int> col_indices;
  bool found_columns = false;

  // Identify columns, assuming (1a).        
  col_indices.insert(make_pair(col_one, 1));
  col_indices.insert(make_pair(kAlleleFileMajorHeader, 2));
  col_indices.insert(make_pair(kAlleleFileMinorHeader, 3));
  col_indices.insert(make_pair(col_two, 4));
  if (GetColumnIndices(file_info, &col_indices)) {
    found_columns = true;
    // Subtract '1', because col_indices starts with index '1', while 'line_parts'
    // below starts with index '0'.
    *col_one_index = col_indices[col_one] - 1;
    *col_two_index = col_indices[col_two] - 1;
    *major_col = col_indices[kAlleleFileMajorHeader] - 1;
    *minor_col = col_indices[kAlleleFileMinorHeader] - 1;
  }

  if (is_study_allele_file) {
    // Identify columns, assuming (1b).        
    if (!found_columns) {
      // Reset col_indices, in case it was modified in GetColumnIndices
      // (which can happen if partial column name matches were found).
      col_indices.clear();
      col_indices.insert(make_pair(kAlleleFilePositionHeader, 1));
      col_indices.insert(make_pair(kAlleleFileMajorHeader, 2));
      col_indices.insert(make_pair(kAlleleFileMinorHeader, 3));
      col_indices.insert(make_pair(col_two, 4));
      if (GetColumnIndices(file_info, &col_indices)) {
        found_columns = true;
        // Subtract '1', because col_indices starts with index '1', while 'line_parts'
        // below starts with index '0'.
        *col_one_index = col_indices[kAlleleFilePositionHeader] - 1;
        *col_two_index = col_indices[col_two] - 1;
        *major_col = col_indices[kAlleleFileMajorHeader] - 1;
        *minor_col = col_indices[kAlleleFileMinorHeader] - 1;
      }
    }

    if (allow_format_two) {
      // Identify columns, assuming (2a).        
      if (!found_columns) {
        // Reset col_indices, in case it was modified in GetColumnIndices
        // (which can happen if partial column name matches were found).
        col_indices.clear();
        col_indices.insert(make_pair(col_one, 1));
        col_indices.insert(make_pair(col_two, 2));
        if (GetColumnIndices(file_info, &col_indices)) {
          found_columns = true;
          // Subtract '1', because col_indices starts with index '1', while 'line_parts'
          // below starts with index '0'.
          *col_one_index = col_indices[col_one] - 1;
          *col_two_index = col_indices[col_two] - 1;
        }
      }

      // Identify columns, assuming (2b).        
      if (!found_columns) {
        // Reset col_indices, in case it was modified in GetColumnIndices
        // (which can happen if partial column name matches were found).
        col_indices.clear();
        col_indices.insert(make_pair(kAlleleFilePositionHeader, 1));
        col_indices.insert(make_pair(col_two, 2));
        if (GetColumnIndices(file_info, &col_indices)) {
          found_columns = true;
          // Subtract '1', because col_indices starts with index '1', while 'line_parts'
          // below starts with index '0'.
          *col_one_index = col_indices[kAlleleFilePositionHeader] - 1;
          *col_two_index = col_indices[col_two] - 1;
        }
      }
    }
  }

  // Identify columns, assuming (3a).        
  if (!found_columns) {
    // Reset col_indices, in case it was modified in GetColumnIndices
    // (which can happen if partial column name matches were found).
    col_indices.clear();
    col_indices.insert(make_pair(col_one, 1));
    col_indices.insert(make_pair(kAlleleFileMajorHeader, 2));
    col_indices.insert(make_pair(kAlleleFileMinorHeader, 3));
    if (GetColumnIndices(file_info, &col_indices)) {
      found_columns = true;
      // Subtract '1', because col_indices starts with index '1', while 'line_parts'
      // below starts with index '0'.
      *col_one_index = col_indices[col_one] - 1;
      *major_col = col_indices[kAlleleFileMajorHeader] - 1;
      *minor_col = col_indices[kAlleleFileMinorHeader] - 1;
    }
  }

  // Identify columns, assuming (3b).        
  if (!found_columns) {
    // Reset col_indices, in case it was modified in GetColumnIndices
    // (which can happen if partial column name matches were found).
    col_indices.clear();
    col_indices.insert(make_pair(kAlleleFilePositionHeader, 1));
    col_indices.insert(make_pair(kAlleleFileMajorHeader, 2));
    col_indices.insert(make_pair(kAlleleFileMinorHeader, 3));
    if (GetColumnIndices(file_info, &col_indices)) {
      found_columns = true;
      // Subtract '1', because col_indices starts with index '1', while 'line_parts'
      // below starts with index '0'.
      if (is_study_allele_file) {
        *col_one_index = col_indices[kAlleleFilePositionHeader] - 1;
      } else {
        *col_two_index = col_indices[kAlleleFilePositionHeader] - 1;
      }
      *major_col = col_indices[kAlleleFileMajorHeader] - 1;
      *minor_col = col_indices[kAlleleFileMinorHeader] - 1;
    }
  }

  if (!found_columns) {
    // All attempts to parse header of SNP_INFO file failed.
    const string message_one = allow_format_two ? " (optionally)" : "";
    const string message_two = !is_study_allele_file ?
      "" : (" (optionally) \"" + string(kAlleleFileStudySnpIdHeader) + "\" and");
    cout << "ERROR: Unable to parse header of SNP_INFO File '"
         << file_info.name_ << "'; check that you specified the correct column "
         << "delimiter and comment character. Note that the Header should have "
         << "a columns titled \"" << col_one << "\" (or \""
         << kAlleleFilePositionHeader << "\"), as well as" << message_two
         << message_one << " the major/minor allele information in columns \""
         << kAlleleFileMajorHeader << "\", and \""
         << kAlleleFileMinorHeader << "\" (it may have other "
         << "columns as well)." << endl;
    return false;
  }

  for (const auto& col_itr : col_indices) {
    if (col_itr.second > *max_column_index) *max_column_index = col_itr.second;
  }

  return true;
}

}  // namespace

string PreMeta::GetSeqMetaObjectNameFromFile(const string& out_file) {
  string to_return = StripSuffixString(out_file, ".RData");
  to_return = StripSuffixString(to_return, ".Rdata");
  size_t dir_index = to_return.find("/");
  while (dir_index != string::npos) {
    to_return = to_return.substr(dir_index + 1);
    dir_index = to_return.find("/");
  }
  return to_return;
}

void PreMeta::PrintInconsistentSnps(const string& outfile) {
  if (snp_to_excluding_study_.empty()) return;

  cout << "WARNING: At least one SNP was skipped in at least one "
       << "study, due to SNP REF and/or ALT alleles not matching "
       << "another (already processed) study.\n"
       << "Check '" << outfile << "' for a list of SNPs that were "
       << "skipped, and the Studies that excluded them." << endl;

  ofstream out_file;
  out_file.open(outfile.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return;
  }

  // Print Header.
  out_file << "SNP\tStudies_Excluding_Snp" << endl;

  for (const pair<Position, set<int>>& snp_and_studies : snp_to_excluding_study_) {
    out_file << PrintPosition(snp_and_studies.first) << "\t"
             << Join(snp_and_studies.second, ",") << "\n";
  }

  out_file.close();
}

void PreMeta::PrintSkippedSnps(
    const string& orig_out_file, const set<Position>& snps_to_skip) {
  if (snps_to_skip.empty()) return;
  size_t index_to_replace = orig_out_file.find("raremetal_score_");
  if (index_to_replace == string::npos) {
    index_to_replace = orig_out_file.find("meta_skat_");
  }

  string outfile = "./snps_with_no_allele_info.txt";
  // Should always happen.
  if (index_to_replace != string::npos) {
    outfile =
        orig_out_file.substr(0, index_to_replace) +
        "snps_with_no_allele_info.txt";
  }

  cout << "WARNING: Not all SNP Positions have allele information. SNPs "
       << "without allele information will be omitted from the output; a "
       << "list of all such SNPs can be found in: '" << outfile << "'" << endl;

  ofstream out_file;
  out_file.open(outfile.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return;
  }

  for (const Position& pos : snps_to_skip) {
    out_file << PrintPosition(pos) << endl;
  }

  out_file.close();
}

bool PreMeta::GetWindows(
    const int window_size, const set<Position>& snp_positions,
    map<Position, Position>* windows) {
  if (snp_positions.empty()) {
    cout << "ERROR in Geting Windows: No SNPs specified for which to get "
         << "windows. Aborting." << endl;
    return false;
  }
  // Keep two iterators, one marking the SNP for which we're looking
  // for the end position, and one walking down positions, looking for
  // the end position.
  set<Position>::const_iterator current_snp_itr = snp_positions.begin();
  set<Position>::const_iterator end_window_itr = snp_positions.begin();
  Chromosome current_chr;
  uint64_t current_pos;
  while (current_snp_itr != snp_positions.end()) {
    current_chr = current_snp_itr->chr_;
    current_pos = current_snp_itr->pos_;
    if (current_chr == Chromosome::CHROMOSOME_UNKNOWN && current_pos == 0) {
      // This Position is specified by SNP_ID instead of CHR:POS, so window
      // cannot be computed for it. Return false.
      cout << "ERROR in GetWindows: Unable to compute a sliding window for "
           << "SNP '" << PrintPosition(*current_snp_itr)
           << "', which is not in CHR:POS format." << endl;
      return false;
    }
    while (end_window_itr != snp_positions.end() &&
           end_window_itr->chr_ == current_chr &&
           end_window_itr->pos_ < current_pos + window_size) {
      if (end_window_itr->chr_ == Chromosome::CHROMOSOME_UNKNOWN &&
          end_window_itr->pos_ == 0) {
        // This Position is specified by SNP_ID instead of CHR:POS, so window
        // cannot be computed for it. Return false.
        cout << "ERROR in GetWindows: Unable to compute a sliding window for "
             << "SNP '" << PrintPosition(*end_window_itr)
             << "', which is not in CHR:POS format." << endl;
        return false;
      }
      ++end_window_itr;
    }
    windows->insert(make_pair(*current_snp_itr, *(--end_window_itr)));
    if (end_window_itr->chr_ == Chromosome::CHROMOSOME_UNKNOWN &&
        end_window_itr->pos_ == 0) {
      // This Position is specified by SNP_ID instead of CHR:POS, so window
      // cannot be computed for it. Return false.
      cout << "ERROR in GetWindows: Unable to compute a sliding window for "
           << "SNP '" << PrintPosition(*end_window_itr)
           << "', which is not in CHR:POS format." << endl;
      return false;
    }
    ++end_window_itr;
    ++current_snp_itr;
  }
  return true;
}

bool PreMeta::WriteGoldenAlleleFile(const string& outfile) {
  ofstream out_file;
  out_file.open(outfile.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Print Header Line.
  out_file << "SnpId\tCHR:POS\tREF\tALT\tFlipped_Studies\n";
  for (const pair<Position, tuple<string, string, set<int>>>& snp :
       snp_to_ref_alt_and_study_) {
    const Position& snp_pos = snp.first;
    out_file << snp_pos.snp_id_ << "\t"
             << (snp_pos.chr_ == Chromosome::CHROMOSOME_UNKNOWN ?
                 "-" : PrintPosition(snp_pos))
             << "\t"
             << get<0>(snp.second) << "\t"
             << get<1>(snp.second) << "\t"
             << Join(get<2>(snp.second), ",") << "\n";
  }

  out_file.close();

  return true;
}

bool PreMeta::ParseGoldenSnpFile(const string& golden_file) {
  FileInfo file_info;
  file_info.name_ = golden_file;
  file_info.delimiter_ = "\t";
  file_info.comment_char_ = "#";

  // Read header, to determine the columns.
  int snp_id_col = -1;
  int pos_col = -1;
  int ref_col = -1;
  int alt_col = -1;
  int max_column_index = -1;
  if (!ReadAlleleFileHeader(
          false, false, file_info,
          kAlleleFileSnpIdHeader, kAlleleFilePositionHeader,
          &snp_id_col, &pos_col, &ref_col, &alt_col, &max_column_index)) {
    return false;
  }

  // Open Allele File.
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting Alleles from allele annotation file: "
         << "Unable to open file '" << file_info.name_ << "'." << endl;
    return false;
  }
  
  // Read Allele File.
  string line;
  int line_num = 0;
  while(getline(file, line)) {
    line_num++;

    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    // Skip Comment lines.
    if (HasPrefixString(line, file_info.comment_char_)) {
      continue;
    }

    // Parse line, splitting it into the columns.
    vector<string> line_parts;
    Split(line, file_info.delimiter_, &line_parts);
    if (line_parts.size() < max_column_index) {
      cout << "ERROR in getting Alleles from allele annotation file: "
           << "Comment char: '" << file_info.comment_char_ << "', delimiter: '"
           << file_info.delimiter_ << "'. "
           << "On line " << line_num << " of Allele file '"
           << file_info.name_ << "', unable to parse line: '"
           << line << "'. Expected at least 3 columns, found "
           << line_parts.size()
           << ":\n" << Join(line_parts, "\n") << endl;
      return false;
    }

    // Parse SNP_ID (the ID that will be used for Meta-Analysis).
    Position snp_pos;
    if (pos_col >= 0 && !ParsePosition(line_parts[pos_col], &snp_pos)) {
      cout << "ERROR: Unable to parse SNP's position on line "
           << line_num << " of the global snp info file '"
           << file_info.name_ << "'." << endl;
      return false;
    }

    // Parse SNP_ID (the ID that was used for Study-Level Analysis).
    if (snp_id_col >= 0) {
      snp_pos.snp_id_ = line_parts[pos_col];
    }

    const string major_str = line_parts[ref_col];
    const string minor_str = line_parts[alt_col];

    if (!snp_to_ref_alt_and_study_.insert(make_pair(
            snp_pos, make_tuple(major_str, minor_str, set<int>()))).second) {
      cout << "ERROR: When reading golden snp info file '"
           << file_info.name_ << "', encountered SNP '"
           << PrintPosition(snp_pos) << "' twice (second occurence on line "
           << line_num << ")." << endl;
      return false;
    }
  }

  return true;
}

bool PreMeta::GetAllelesFromAlleleFile(
    const bool is_rm_or_metaskat, const int study_num,
    const FileInfo& file_info,
    map<Position, SnpInfo>* allele_info) {
  int meta_pos_col = -1;
  int study_pos_col = -1;
  int major_col = -1;
  int minor_col = -1;
  int max_column_index = -1;
  if (!ReadAlleleFileHeader(
          true, is_rm_or_metaskat, file_info,
          kAlleleFileMetaSnpIdHeader, kAlleleFileStudySnpIdHeader,
          &meta_pos_col, &study_pos_col, &major_col, &minor_col, &max_column_index)) {
    return false;
  }

  // Open Allele File.
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR in getting Alleles from allele annotation file: "
         << "Unable to open file '" << file_info.name_ << "'." << endl;
    return false;
  }
  
  // Read Allele File.
  string line;
  int line_num = 0;
  while(getline(file, line)) {
    line_num++;

    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    // Skip Comment lines.
    if (HasPrefixString(line, file_info.comment_char_)) {
      continue;
    }

    // Parse line, splitting it into the columns.
    vector<string> line_parts;
    Split(line, file_info.delimiter_, &line_parts);
    if (line_parts.size() < max_column_index) {
      cout << "ERROR in getting Alleles from allele annotation file: "
           << "Comment char: '" << file_info.comment_char_ << "', delimiter: '"
           << file_info.delimiter_ << "'. "
           << "On line " << line_num << " of Allele file '"
           << file_info.name_ << "', unable to parse line: '"
           << line << "'. Expected at least 3 columns, found "
           << line_parts.size()
           << ":\n" << Join(line_parts, "\n") << endl;
      return false;
    }

    // Parse SNP_ID (the ID that will be used for Meta-Analysis).
    Position meta_pos;
    if (!ParsePosition(line_parts[meta_pos_col], &meta_pos)) {
      meta_pos.snp_id_ = line_parts[meta_pos_col];
    }

    // Parse SNP_ID (the ID that was used for Study-Level Analysis).
    Position study_pos;
    if (study_pos_col >= 0) {
      if (!ParsePosition(line_parts[study_pos_col], &study_pos)) {
        study_pos.snp_id_ = line_parts[study_pos_col];
      }
    }

    string major_str = "";
    string minor_str = "";
    if (major_col >= 0 && minor_col >= 0) {
      major_str = line_parts[major_col];
      minor_str = line_parts[minor_col];
    }

    // Insert an entry into allele_info using Key study_pos (if present),
    // otherwise meta_pos.
    const Position pos_key = study_pos_col >= 0 ? study_pos : meta_pos;
    pair<map<Position, SnpInfo>::iterator, bool> info_insert_itr =
        allele_info->insert(make_pair(pos_key, SnpInfo()));
    // If we already have Allele information for this SNP, check that the
    // current line is consistent with it.
    if (!info_insert_itr.second) {
      const SnpInfo& old_info = info_insert_itr.first->second;
      if (old_info.pos_.chr_ != meta_pos.chr_ ||
          old_info.pos_.pos_ != meta_pos.pos_ ||
          old_info.pos_.snp_id_ != meta_pos.snp_id_ ||
          (old_info.major_allele_ != major_str &&
           !old_info.major_allele_.empty() && !major_str.empty()) ||
          (old_info.minor_allele_ != minor_str &&
           !old_info.minor_allele_.empty() && !minor_str.empty())) {
        cout << "ERROR in getting Alleles from allele annotation file: "
             << "On line " << line_num << " of Allele file '"
             << file_info.name_ << "', already have SNP allele information "
             << "for Position " << PrintPosition(pos_key)
             << " which is inconsistent with present values. Aborting." << endl;
        return false;
      }
    }
    SnpInfo& current_info = info_insert_itr.first->second;
    current_info.pos_ = meta_pos;
    if (major_col >= 0 && minor_col >= 0) {
      current_info.major_allele_ = major_str;
      current_info.minor_allele_ = minor_str;
    }

    // Insert this SNP's info to snp_to_ref_alt_and_study_, if REF/ALT info
    // is available.
    if (!is_rm_or_metaskat) {
      if (major_col < 0 || minor_col < 0) {
        cout << "ERROR: Line " << line_num << " of SNP INFO file '"
             << file_info.name_ << " is missing requisite REF and/or ALT "
             << "column." << endl;
        return false;
      }
      tuple<string, string, set<int>>* snp_info =
          FindOrInsert(pos_key, snp_to_ref_alt_and_study_,
                       make_tuple(major_str, minor_str, set<int>()));
      // Check that REF/ALT alleles match existing entry for this SNP. Otherwise,
      // either add this SNP to snp_to_ref_alt_and_study_ or snp_to_excluding_study_,
      // as appropriate.
      if (get<0>(*snp_info) != major_str || get<1>(*snp_info) != minor_str) {
        if (get<1>(*snp_info) != major_str || get<0>(*snp_info) != minor_str) {
          // This is a more aggregious difference than just a swapping of
          // REF <-> ALT. Add snp to snp_to_excluding_study_.
          set<int>* excluding_studies =
              FindOrInsert(pos_key, snp_to_excluding_study_, set<int>());
          excluding_studies->insert(study_num);
        } else {
          get<2>(*snp_info).insert(study_num);
        }
      }
    }
  }

  return true;
}

/* DEPRECATED */
/*
bool PreMeta::GetAllelesFromRareMetalFiles(
    const vector<FileInfo>& files,
    map<Position, pair<string, string>>* allele_info) {
  for (const FileInfo& file_info : files) {
    // Read score_file, which has info about all SNPS, including Allele Info.
    int num_samples, num_snps;
    map<Position, SnpInfo> snp_info;
    if (!GetInfoFromRareMetalScoreFile(
            file_info, &num_samples, &num_snps, &snp_info)) {
      cout << "ERROR in getting Alleles from RareMetal files: "
           << "Unable to parse RareMetal score file '" << file_info.name_
           << "'. Aborting.\n";
      return false;
    }
    for (const auto& snp_itr : snp_info) {
      const Position& pos = snp_itr.first;
      const SnpInfo& info = snp_itr.second;
      pair<map<Position, pair<string, string>>::iterator, bool>
          info_insert_itr = allele_info->insert(make_pair(
              pos, make_pair(info.major_allele_, info.minor_allele_)));
      // If we already have Allele information for this SNP, check that the
      // current line is consistent with it.
      if (!info_insert_itr.second) {
        const pair<string, string>& old_pair = info_insert_itr.first->second;
        if (old_pair.first != info.major_allele_ ||
            old_pair.second != info.minor_allele_) {
        cout << "ERROR in getting Alleles from RareMetal files: "
               << "In RareMetalWorker's score file '"
               << file_info.name_ << "': Already had SNP allele information "
               << "for Position " << PrintPosition(pos)
               << " which is inconsistent with present values. Previously: "
               << "(Major = " << old_pair.first << ", Minor = "
               << old_pair.second << "); vs. Currently: (Major = "
               << info.major_allele_ << ", Minor = " << info.minor_allele_
               << "). Aborting." << endl;
          return false;
        }
      }
    }
  }
  return true;
}
*/ 
/* END DEPRECATED */

/* DEPRECATED */
/*
bool PreMeta::GetAllelesFromMetaSkatFiles(
    const vector<FileInfo>& files,
    map<Position, pair<string, string>>* allele_info) {
  for (const FileInfo& file_info : files) {
    // Parse MInfo File.
    int num_samples;
    map<string, vector<Position>> gene_to_snp_positions;
    map<pair<Position, string>, uint64_t> start_lines;
    map<Position, SnpInfo> snp_info;
    if (!ParseMInfoFile(
            file_info,
            &num_samples, &gene_to_snp_positions, &start_lines, &snp_info)) {
      cout << "ERROR in getting Alleles from MetaSkat files: "
           << "Unable to parse MetaSKAT MInfo file '" << file_info.name_
           << "'. Aborting.\n";
      return false;
    }
    for (const auto& snp_itr : snp_info) {
      const Position& pos = snp_itr.first;
      const SnpInfo& info = snp_itr.second;
      pair<map<Position, pair<string, string>>::iterator, bool>
          info_insert_itr = allele_info->insert(make_pair(
              pos, make_pair(info.major_allele_, info.minor_allele_)));
      // If we already have Allele information for this SNP, check that the
      // current line is consistent with it.
      if (!info_insert_itr.second) {
        const pair<string, string>& old_pair = info_insert_itr.first->second;
        if (old_pair.first != info.major_allele_ ||
            old_pair.second != info.minor_allele_) {
          cout << "ERROR in getting Alleles from MetaSkat files: "
               << "In MetaSKAT's MInfo file '"
               << file_info.name_ << "': Already had SNP allele information "
               << "for Position " << PrintPosition(pos)
               << " which is inconsistent with present values. Previously: "
               << "(Major = " << old_pair.first << ", Minor = "
               << old_pair.second << "); vs. Currently: (Major = "
               << info.major_allele_ << ", Minor = " << info.minor_allele_
               << "). Aborting." << endl;
          return false;
        }
      }
    }
  }
  return true;
}
*/
/* END DEPRECATED */

/* DEPRECATED */
/* 
bool PreMeta::GetAllelesFromAlleleInfo(
    const AlleleInfo& info,
    map<Position, pair<string, string>>* allele_info) {
  for (const FileInfo& file_info : info.supp_files_) {
    if (!GetAllelesFromAlleleFile(file_info, allele_info)) {
      return false;
    }
  }

  if (!GetAllelesFromRareMetalFiles(
          info.raremetal_files_, allele_info)) {
    return false;
  }

  if (!GetAllelesFromMetaSkatFiles(
          info.metaskat_files_, allele_info)) {
    return false;
  }
  return true;
}
*/
/* END DEPRECATED */

bool PreMeta::RareMetalToRareMetal(
    const int study_num, const bool is_new_version_from,
    const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& score_file_info,
    const FileInfo& covariance_file_info,
    const string& score_files, const string& cov_files,
    const string& out_score_file, const string& out_cov_file) {
  // Read score_file, which has info about all SNPS, including:
  // postion, MAF, u-stat, and v-stat.
  int num_samples, num_snps;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!GetInfoFromRareMetalScoreFile(
          study_num, allele_info, score_file_info,
          &num_samples, &num_snps, &monomorphic_snps,
          &snp_to_ref_alt_and_study_, &snp_to_excluding_study_, &snp_info)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to parse "
         << "score file '" << score_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Read covariance_file, which has for every SNP (in a sliding window)
  // the covariance of that SNP with it's neighbors.
  set<Position> snps_without_score_info;
  map<Position, map<Position, double>> snp_to_cov_w_neighbors;
  if (!GetInfoFromRareMetalCovarianceFile(
          is_new_version_from, true, study_num, num_samples, covariance_file_info,
          allele_info, snp_info, snp_to_excluding_study_,
          &snps_without_score_info, &snp_to_cov_w_neighbors)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to parse covariance "
         << "file '" << covariance_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Get sigma, if using old version of RareMetalWorker.
  double sigma = 1.0;
  if (!is_new_version_from && !GetSigmaFromRareMetalFiles(
          num_samples, snp_to_cov_w_neighbors, snp_info, &sigma)) {
    cout << "ERROR in converting RareMetal to RareMetal: "
         << "Unable to compute ResidualVariance from covariances. "
         << "Aborting." << endl;
    return false;
  }

  // Modify snp_info by adjusting by an appropriate constant.
  const double u_stat_rescale = FloatEq(sigma, 0.0) ? 1.0 : 1.0 / sigma;

  // Print RareMetal Score and Covariance files.
  if (!PrintToRareMetal(
          study_num, num_samples, rescale, u_stat_rescale, monomorphic_snps,
          snps_without_score_info, snp_to_ref_alt_and_study_,
          snp_info, snp_to_cov_w_neighbors, out_score_file, out_cov_file)) {
    cout << "ERROR in printing RareMetal to RareMetal. Aborting.\n";
    return false;
  }

  // Add above files to Score and Covariance files list.
  return AddFilesToList(score_files, cov_files, out_score_file, out_cov_file);
}

bool PreMeta::RareMetalToMass(
    const int study_num, const bool is_new_version,
    const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& group_file_info,
    const FileInfo& score_file_info,
    const FileInfo& covariance_file_info,
    const string& mass_script_file, const string& out_file) {
  // Read score_file, which has info about all SNPS, including:
  // postion, MAF, u-stat, and v-stat.
  int num_samples, num_snps;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!GetInfoFromRareMetalScoreFile(
          study_num, allele_info, score_file_info,
          &num_samples, &num_snps, &monomorphic_snps,
          &snp_to_ref_alt_and_study_, &snp_to_excluding_study_, &snp_info)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to parse "
         << "score file '" << score_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Read group file, which has for each gene, the chr:pos of every SNP.
  // In paritcular, we populate the following map, which maps gene name
  // to an internal map, keyed by chromosome and value (a set of) the
  // corresponding positions within that (gene, chromosome) that have a SNP.
  set<Position> snps_without_score_info;
  map<string, vector<Position>> gene_to_snp_positions;
  if (!GetInfoFromRareMetalGroupFile(
          study_num, allele_info, snp_info, snp_to_excluding_study_,
          group_file_info, &snps_without_score_info, &gene_to_snp_positions)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to parse group file "
         << "'" << group_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Read covariance_file, which has for every SNP (in a sliding window)
  // the covariance of that SNP with it's neighbors.
  map<Position, map<Position, double>> snp_to_cov_w_neighbors;
  if (!GetInfoFromRareMetalCovarianceFile(
          is_new_version, false, study_num, num_samples, covariance_file_info,
          allele_info, snp_info, snp_to_excluding_study_,
          &snps_without_score_info, &snp_to_cov_w_neighbors)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to parse covariance "
         << "file '" << covariance_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Get sigma, if using old version of RareMetalWorker.
  double sigma = 1.0;
  if (!is_new_version && !GetSigmaFromRareMetalFiles(
          num_samples, snp_to_cov_w_neighbors, snp_info, &sigma)) {
    cout << "ERROR in converting RareMetal to Mass: "
         << "Unable to compute ResidualVariance from covariances. "
         << "Aborting." << endl;
    return false;
  }

  // Modify snp_info by adjusting by an appropriate constant.
  const double v_stat_rescale = FloatEq(sigma, 0.0) ? 1.0 : 1.0 / pow(sigma, 2.0);
  if (!FloatEq(v_stat_rescale, 1.0)) {
    RescaleSnpInfo(1.0, v_stat_rescale, &snp_info);
  }

  // Compute Covariances (by Gene).
  set<Position> snps_to_skip;
  map<pair<Position, Position>, double> covariances;
  if (!ComputeCovarianceByGene(
          num_samples, sigma,
          gene_to_snp_positions, snp_to_cov_w_neighbors,
          &snps_to_skip, &covariances)) {
    cout << "ERROR in converting RareMetal to Mass: Unable to get "
         << "gene covariances. Aborting.\n";
    return false;
  }
  
  // Make sure we don't have partially info on any snps (i.e. a SNP appears
  // in Score or Covariance file, but no self-covariance is available).
  if (!SanityCheckSnpsToSkip(
          snps_to_skip, monomorphic_snps, snp_info, covariances)) {
    return false;
  }

  // Print Mass script file.
  if (!PrintMassScriptFile(out_file, mass_script_file)) return false;

  // Print to file.
  if (!PrintToMass(study_num, num_samples, rescale, sigma, out_file,
                   monomorphic_snps, snp_to_ref_alt_and_study_,
                   gene_to_snp_positions, snp_info, covariances)) {
    cout << "ERROR in converting RareMetal to Mass: Failed to Print. Aborting."
         << endl;
    return false;
  }

  return true;
}

bool PreMeta::RareMetalToMetaSkat(
    const int study_num, const bool is_new_version, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& group_file_info,
    const FileInfo& score_file_info,
    const FileInfo& covariance_file_info,
    const string& minfo_outfile, const string& mssd_outfile) {
  // Read score_file, which has info about all SNPS, including:
  // postion, MAF, u-stat, and v-stat.
  int num_samples, num_snps;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!GetInfoFromRareMetalScoreFile(
          study_num, allele_info, score_file_info,
          &num_samples, &num_snps, &monomorphic_snps,
          &snp_to_ref_alt_and_study_, &snp_to_excluding_study_, &snp_info)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to parse "
         << "score file '" << score_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Read group file, which has for each gene, the chr:pos of every SNP.
  // In paritcular, we populate the following map, which maps gene name
  // to an internal map, keyed by chromosome and value (a set of) the
  // corresponding positions within that (gene, chromosome) that have a SNP.
  set<Position> snps_without_score_info;
  map<string, vector<Position>> gene_to_snp_positions;
  if (!GetInfoFromRareMetalGroupFile(
          study_num, allele_info, snp_info, snp_to_excluding_study_, group_file_info,
          &snps_without_score_info, &gene_to_snp_positions)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to parse group file '"
         << group_file_info.name_ << "'. Aborting.\n";
    return false;
  }

  // Read covariance_file, which has for every SNP (in a sliding window)
  // the covariance of that SNP with it's neighbors.
  map<Position, map<Position, double>> snp_to_cov_w_neighbors;
  if (!GetInfoFromRareMetalCovarianceFile(
          is_new_version, false, study_num, num_samples, covariance_file_info,
          allele_info, snp_info, snp_to_excluding_study_,
          &snps_without_score_info, &snp_to_cov_w_neighbors)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to parse covariance file '"
         << covariance_file_info.name_ << "'. Aborting.\n";
    return false;
  }

  // Get sigma, if using old version of RareMetalWorker.
  double sigma = 1.0;
  if (!is_new_version && !GetSigmaFromRareMetalFiles(
          num_samples, snp_to_cov_w_neighbors, snp_info, &sigma)) {
    cout << "ERROR in converting RareMetal to MetaSkat::GetSigmaFromRareMetalFiles. "
         << "Aborting.\n";
    return false;
  }

  // Modify snp_info by adjusting by an appropriate constant.
  const double v_stat_rescale = FloatEq(sigma, 0.0) ? 1.0 : 1.0 / pow(sigma, 2.0);
  if (!FloatEq(v_stat_rescale, 1.0)) {
    RescaleSnpInfo(1.0, v_stat_rescale, &snp_info);
  }

  // Compute Covariances (by Gene).
  set<Position> snps_to_skip;
  map<pair<Position, Position>, double> covariances;
  if (!ComputeCovarianceByGene(
          num_samples, sigma,
          gene_to_snp_positions, snp_to_cov_w_neighbors,
          &snps_to_skip, &covariances)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to get "
         << "gene covariances. Aborting.\n";
    return false;
  }
  
  // Make sure we don't have partially info on any snps (i.e. a SNP appears
  // in Score or Covariance file, but no self-covariance is available).
  if (!SanityCheckSnpsToSkip(
          snps_to_skip, monomorphic_snps, snp_info, covariances)) {
    return false;
  }

  if (!WriteMInfoFile(
          study_num, num_samples, num_snps, rescale, sigma, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions, snp_info, minfo_outfile)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to write .Minfo "
         << "file. Aborting.\n";
    return false;
  }
  
  if (!WriteMssdFile(
          study_num, rescale, monomorphic_snps, snp_to_ref_alt_and_study_,
          gene_to_snp_positions, covariances, mssd_outfile)) {
    cout << "ERROR in converting RareMetal to MetaSkat: Unable to write .MSSD "
         << "file. Aborting.\n";
    return false;
  }

  return true;
}

bool PreMeta::RareMetalToSeqMeta(
    const int study_num, const bool is_new_version, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& group_file_info,
    const FileInfo& score_file_info,
    const FileInfo& covariance_file_info,
    const string& out_file) {
  // Read score_file, which has info about all SNPS, including:
  // postion, MAF, u-stat, and v-stat.
  int num_samples, num_snps;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!GetInfoFromRareMetalScoreFile(
          study_num, allele_info, score_file_info,
          &num_samples, &num_snps, &monomorphic_snps,
          &snp_to_ref_alt_and_study_, &snp_to_excluding_study_, &snp_info)) {
    cout << "ERROR in converting RareMetal to SeqMeta: Unable to parse "
         << "score file '" << score_file_info.name_ << "'. Aborting." << endl;
    return false;
  }

  // Read group file, which has for each gene, the chr:pos of every SNP.
  // In paritcular, we populate the following map, which maps gene name
  // to an internal map, keyed by chromosome and value (a set of) the
  // corresponding positions within that (gene, chromosome) that have a SNP.
  set<Position> snps_without_score_info;
  map<string, vector<Position>> gene_to_snp_positions;
  if (!GetInfoFromRareMetalGroupFile(
          study_num, allele_info, snp_info, snp_to_excluding_study_, group_file_info,
          &snps_without_score_info, &gene_to_snp_positions)) {
    cout << "ERROR in converting RareMetal to SeqMeta: Unable to parse group file '"
         << group_file_info.name_ << "'. Aborting.\n";
    return false;
  }

  // Read covariance_file, which has for every SNP (in a sliding window)
  // the covariance of that SNP with it's neighbors.
  map<Position, map<Position, double>> snp_to_cov_w_neighbors;
  if (!GetInfoFromRareMetalCovarianceFile(
          is_new_version, false, study_num, num_samples, covariance_file_info,
          allele_info, snp_info, snp_to_excluding_study_,
          &snps_without_score_info, &snp_to_cov_w_neighbors)) {
    cout << "ERROR in converting RareMetal to SeqMeta: Unable to parse covariance file '"
         << covariance_file_info.name_ << "'. Aborting.\n";
    return false;
  }

  // Get sigma, if using old version of RareMetalWorker.
  double sigma = 1.0;
  if (!is_new_version && !GetSigmaFromRareMetalFiles(
          num_samples, snp_to_cov_w_neighbors, snp_info, &sigma)) {
    cout << "ERROR in converting RareMetal to SeqMeta: Unable to get "
         << "Residual Variance from covariances. Aborting.\n";
    return false;
  }

  // Modify snp_info by adjusting by an appropriate constant.
  const double v_stat_rescale = FloatEq(sigma, 0.0) ? 1.0 : 1.0 / pow(sigma, 2.0);
  if (!FloatEq(v_stat_rescale, 1.0)) {
    RescaleSnpInfo(1.0, v_stat_rescale, &snp_info);
  }

  // Compute Covariances (by Gene).
  set<Position> snps_to_skip;
  map<pair<Position, Position>, double> covariances;
  if (!ComputeCovarianceByGene(
          num_samples, sigma,
          gene_to_snp_positions, snp_to_cov_w_neighbors,
          &snps_to_skip, &covariances)) {
    cout << "ERROR in converting RareMetal to SeqMeta: Unable to transform "
         << "sliding-window covariances to gene-based covariances. Aborting."
         << endl;
    return false;
  }
  
  // Make sure we don't have partially info on any snps (i.e. a SNP appears
  // in Score or Covariance file, but no self-covariance is available).
  if (!SanityCheckSnpsToSkip(
          snps_to_skip, monomorphic_snps, snp_info, covariances)) {
    return false;
  }

  // Convert info into a SeqMeta RDataTree.
  vector<string*>* tags = new vector<string*>();
  RDataNode* root = new RDataNode();
  const string obj_name = GetSeqMetaObjectNameFromFile(out_file);
  if (!ConstructSeqMetaTree(
          obj_name, study_num, num_samples, rescale, sigma, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions, snp_info,
          covariances, tags, root)) {
    cout << "ERROR in printing RareMetal to SeqMeta. Aborting.\n";
    return false;
  }

  // Print SeqMeta .Rdata file.
  if (!PrintToSeqMeta(*root, out_file)) {
    cout << "ERROR in Printing RareMetal to SeqMeta. Aborting.\n";
    return false;
  }

  // Delete the RDataTree.
  set<string*> sexptype_tags;
  for (string* tag : *tags) sexptype_tags.insert(tag);
  DeleteRDataTree(root, sexptype_tags);
  DeleteRDataTags(*tags);
  delete tags;

  return true;
}

bool PreMeta::MassToRareMetal(
    const int study_num, const int window_size, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const string& score_files, const string& cov_files,
    const string& out_score_file, const string& out_cov_file) {
  // Read Residual Variance from Mass input file.
  int num_samples = 0;
  if (!GetNumSamplesFromMassFile(mass_file_info, &num_samples)) {
    return false;
  }

  // Read Mass input file.
  int num_snps = 0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  // TODO(PHB): Determine if this is the proper time to do 'rescale'ing.
  if (!GetInfoFromMassFile(
          true, study_num, allele_info, mass_file_info, snp_to_excluding_study_,
          &num_snps, &gene_to_snp_positions, &snps_to_skip, &monomorphic_snps,
          &snp_info, &covariances)) {
    cout << "ERROR in converting Mass to RareMetal: Unable to parse MASS score "
         << "file. Aborting.\n";
    return false;
  }

  // Read Allele file.
  if (!CopyAlleleInfoToSnpInfo(allele_info, &snp_info)) {
    cout << "ERROR in converting Mass to RareMetal: Unable to get Allele "
         << "information. Aborting.\n";
    return false;
  }
  
  // Print to file Snps that will be skipped due to no Allele Info.
  PrintSkippedSnps(out_score_file, snps_to_skip);

  // Compute Windows.
  const int window = window_size == -1 ? kRareMetalWindowSize : window_size;
  map<Position, Position> windows;
  if (!GetWindows(window, Keys(snp_info), &windows)) {
    cout << "ERROR in converting Mass to RareMetal: Unable to compute windows. "
         << "Aborting.\n";
    return false;
  }

  // Print RareMetal Score and Covariance files.
  if (!PrintToRareMetal(
          study_num, num_samples, rescale, 1.0, 1.0,
          monomorphic_snps, snps_to_skip, snp_to_ref_alt_and_study_,
          windows, snp_info, covariances,
          out_score_file, out_cov_file)) {
    cout << "ERROR in printing Mass to RareMetal. Aborting.\n";
    return false;
  }

  // Add above files to Score and Covariance files list.
  return AddFilesToList(score_files, cov_files, out_score_file, out_cov_file);
}

bool PreMeta::SeqMetaToMass(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const string& seq_meta_rdata_file,
    const string& mass_script_file, const string& out_file) {
  // Parse RData File.
  vector<string*> sexptype_tags;
  RDataNode* root = new RDataNode();
  if (!ReadRDataFile(seq_meta_rdata_file, &sexptype_tags, root)) {
    /*
    // PHB Temp: Print SeqMeta tree.
    ofstream PHB_outfile;
    const string& PHB_outname = "seq_meta_bad_input_to_printed_tree_full.txt";
    PHB_outfile.open(PHB_outname.c_str());
    cout << "\nPHB Printing SeqMeta input from '" << seq_meta_rdata_file
         << "' to: " << PHB_outname << endl;
    PrintRDataTree(*root, "", PHB_outfile);
    PHB_outfile.close();
    //END PHB Temp
    */

    set<string*> tags;
    for (string* tag : sexptype_tags) tags.insert(tag);
    DeleteRDataTree(root, tags);
    DeleteRDataTags(sexptype_tags);
    return false;
  }

  /*
  // PHB Temp: Print SeqMeta tree.
  ofstream PHB_outfile;
  const string& PHB_outname = "seq_meta_input_to_printed_tree_full.txt";
  PHB_outfile.open(PHB_outname.c_str());
  cout << "\nPHB Printing SeqMeta input from '" << seq_meta_rdata_file
       << "' to: " << PHB_outname << endl;
  PrintRDataTree(*root, "", PHB_outfile);
  PHB_outfile.close();
  //END PHB Temp
  */

  // Convert the parsed RData File, now stored by an RData Tree,
  // into the expected structure of SeqMeta output.
  int num_samples = -1;
  int num_snps = 0;
  double sey = -1.0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  bool parsed_tree =
      ParseSeqMetaRDataTree(
          false, study_num, allele_info, snp_to_excluding_study_,
          *root, &num_samples, &num_snps, &sey,
          &snps_to_skip, &monomorphic_snps, &gene_to_snp_positions,
          &snp_info, &covariances);

  // Delete the RDataTree.
  set<string*> tags;
  for (string* tag : sexptype_tags) tags.insert(tag);
  DeleteRDataTree(root, tags);
  DeleteRDataTags(sexptype_tags);

  if (!parsed_tree) return false;

  // Print Mass script file.
  if (!PrintMassScriptFile(out_file, mass_script_file)) return false;

  return PrintToMass(study_num, num_samples, rescale, pow(sey, 2.0),
                     out_file, monomorphic_snps, snp_to_ref_alt_and_study_,
                     gene_to_snp_positions, snp_info, covariances);
}

bool PreMeta::SeqMetaToRareMetal(
    const int study_num, const int window_size, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const string& seq_meta_rdata_file,
    const string& score_files, const string& cov_files,
    const string& out_score_file, const string& out_cov_file) {
  // Parse RData File.
  vector<string*> sexptype_tags;
  RDataNode* root = new RDataNode();
  if (!ReadRDataFile(seq_meta_rdata_file, &sexptype_tags, root)) {
    set<string*> tags;
    for (string* tag : sexptype_tags) tags.insert(tag);
    DeleteRDataTree(root, tags);
    DeleteRDataTags(sexptype_tags);
    return false;
  }

  // Convert the parsed RData File, now stored by an RData Tree,
  // into the expected structure of SeqMeta output.
  int num_samples = -1;
  int num_snps = 0;
  double sey = -1.0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  bool parsed_tree =
      ParseSeqMetaRDataTree(
          true, study_num, allele_info, snp_to_excluding_study_,
          *root, &num_samples, &num_snps, &sey,
          &snps_to_skip, &monomorphic_snps, &gene_to_snp_positions,
          &snp_info, &covariances);

  // Delete the RDataTree.
  set<string*> tags;
  for (string* tag : sexptype_tags) tags.insert(tag);
  DeleteRDataTree(root, tags);
  DeleteRDataTags(sexptype_tags);

  if (!parsed_tree) return false;

  // Read Allele file.
  if (!CopyAlleleInfoToSnpInfo(allele_info, &snp_info)) {
    cout << "ERROR in converting SeqMeta to RareMetal: Failure in reading "
         << "allele information in annotation file. Aborting.\n";
    return false;
  }
  
  // Print to file Snps that will be skipped due to no Allele Info.
  PrintSkippedSnps(out_score_file, snps_to_skip);

  // Compute Windows.
  const int window = window_size == -1 ? kRareMetalWindowSize : window_size;
  map<Position, Position> windows;
  if (!GetWindows(window, Keys(snp_info), &windows)) {
    cout << "ERROR in converting SeqMeta to RareMetal: Unable to get windows. "
         << "Aborting.\n";
    return false;
  }

  // Print RareMetal Score and Covariance files.
  const double sigma_sq = pow(sey, 2.0);
  const double cov_rescale = FloatEq(sigma_sq, 0.0) ? 1.0 : 1.0 / sigma_sq;
  if (!PrintToRareMetal(
          study_num, num_samples, rescale, cov_rescale, cov_rescale,
          monomorphic_snps, snps_to_skip, snp_to_ref_alt_and_study_,
          windows, snp_info, covariances, out_score_file, out_cov_file)) {
    cout << "ERROR in printing SeqMeta to RareMetal. Aborting.\n";
    return false;
  }

  // Add above files to Score and Covariance files list.
  return AddFilesToList(score_files, cov_files, out_score_file, out_cov_file);
}

bool PreMeta::SeqMetaToMetaSkat(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const string& seq_meta_rdata_file,
    const string& minfo_outfile, const string& mssd_outfile) {
  // Parse RData File.
  vector<string*> sexptype_tags;
  RDataNode* root = new RDataNode();
  if (!ReadRDataFile(seq_meta_rdata_file, &sexptype_tags, root)) {
    set<string*> tags;
    for (string* tag : sexptype_tags) tags.insert(tag);
    DeleteRDataTree(root, tags);
    DeleteRDataTags(sexptype_tags);
    return false;
  }

  // Convert the parsed RData File, now stored by an RData Tree,
  // into the expected structure of SeqMeta output.
  int num_samples = -1;
  int num_snps = 0;
  double sey = -1.0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  bool parsed_tree =
      ParseSeqMetaRDataTree(
          true, study_num, allele_info, snp_to_excluding_study_,
          *root, &num_samples, &num_snps, &sey,
          &snps_to_skip, &monomorphic_snps, &gene_to_snp_positions,
          &snp_info, &covariances);

  // Delete the RDataTree.
  set<string*> tags;
  for (string* tag : sexptype_tags) tags.insert(tag);
  DeleteRDataTree(root, tags);
  DeleteRDataTags(sexptype_tags);

  if (!parsed_tree) return false;

  // Rescale scores and covariances.
  const double sigma_sq = pow(sey, 2.0);
  const double cov_rescale = FloatEq(sigma_sq, 0.0) ? 1.0 : 1.0 / sigma_sq;
  if (!FloatEq(cov_rescale, 1.0)) {
    RescaleSnpInfo(cov_rescale, cov_rescale, &snp_info);
    RescaleCovariances(cov_rescale, &covariances);
  }

  // Read Allele file.
  if (!CopyAlleleInfoToSnpInfo(allele_info, &snp_info)) {
    cout << "ERROR in converting SeqMeta to MetaSkat: failed parsing "
         << "allele info in the annotation file. Aborting.\n";
    return false;
  }
  
  // Print to file Snps that will be skipped due to no Allele Info.
  PrintSkippedSnps(minfo_outfile, snps_to_skip);

  if (!WriteMInfoFile(
          study_num, num_samples, num_snps, rescale, pow(sey, 2.0), snps_to_skip,
          monomorphic_snps, snp_to_ref_alt_and_study_, gene_to_snp_positions,
          snp_info, minfo_outfile)) {
    cout << "ERROR in converting SeqMeta to MetaSkat: Unable to write "
         << ".MInfo file. Aborting.\n";
    return false;
  }
  
  if (!WriteMssdFile(
          study_num, rescale, snps_to_skip, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions,
          covariances, mssd_outfile)) {
    cout << "ERROR in converting SeqMeta to MetaSkat: Unable to write "
         << ".MSSD file. Aborting.\n";
    return false;
  }

  return true;
}

bool PreMeta::SeqMetaToSeqMeta(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const string& seq_meta_rdata_file,
    const string& out_file) {
  // Parse RData File.
  vector<string*> sexptype_tags;
  RDataNode* root = new RDataNode();
  if (!ReadRDataFile(seq_meta_rdata_file, &sexptype_tags, root)) {
    set<string*> tags;
    for (string* tag : sexptype_tags) tags.insert(tag);
    DeleteRDataTree(root, tags);
    DeleteRDataTags(sexptype_tags);
    return false;
  }

  // Convert the parsed RData File, now stored by an RData Tree,
  // into the expected structure of SeqMeta output.
  int num_samples = -1;
  int num_snps = 0;
  double sey = -1.0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  bool parsed_tree =
      ParseSeqMetaRDataTree(
          false, study_num, allele_info, snp_to_excluding_study_,
          *root, &num_samples, &num_snps, &sey,
          &snps_to_skip, &monomorphic_snps, &gene_to_snp_positions,
          &snp_info, &covariances);

  // Delete the RDataTree.
  set<string*> tags;
  for (string* tag : sexptype_tags) tags.insert(tag);
  DeleteRDataTree(root, tags);
  DeleteRDataTags(sexptype_tags);

  if (!parsed_tree) return false;

  // Convert info into a SeqMeta RDataTree.
  vector<string*>* new_tags = new vector<string*>();
  RDataNode* new_root = new RDataNode();
  const string obj_name = GetSeqMetaObjectNameFromFile(out_file);
  if (!ConstructSeqMetaTree(
          obj_name, study_num, num_samples, rescale, pow(sey, 2.0),
          monomorphic_snps, snp_to_ref_alt_and_study_,
          gene_to_snp_positions, snp_info, covariances, new_tags, new_root)) {
    cout << "ERROR in printing output for Mass to SeqMeta. Aborting.\n";
    return false;
  }

  // Print SeqMeta .Rdata file.
  if (!PrintToSeqMeta(*new_root, out_file)) {
    cout << "ERROR in Printing output for Mass to SeqMeta. Aborting.\n";
    return false;
  }

  // Delete the RDataTree.
  set<string*> new_sexptype_tags;
  for (string* tag : *new_tags) new_sexptype_tags.insert(tag);
  DeleteRDataTree(new_root, new_sexptype_tags);
  DeleteRDataTags(*new_tags);
  delete new_tags;

  return true;
}

bool PreMeta::MassToMass(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const string& mass_script_file, const string& out_file) {
  int num_samples = 0;
  if (!GetNumSamplesFromMassFile(mass_file_info, &num_samples)) {
    return false;
  }

  // Read Mass input file.
  int num_snps = 0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  // TODO(PHB): Determine if this is the proper time to do 'rescale'ing.
  if (!GetInfoFromMassFile(
          false, study_num, allele_info, mass_file_info, snp_to_excluding_study_,
          &num_snps, &gene_to_snp_positions, &snps_to_skip, &monomorphic_snps,
          &snp_info, &covariances)) {
    cout << "ERROR in converting Mass to Mass: Unable to parse MASS score "
         << "file. Aborting.\n";
    return false;
  }

  // Print Mass script file.
  if (!PrintMassScriptFile(out_file, mass_script_file)) return false;

  return PrintToMass(study_num, num_samples, rescale, 1.0,
                     out_file, monomorphic_snps, snp_to_ref_alt_and_study_,
                     gene_to_snp_positions, snp_info, covariances);
}

bool PreMeta::MassToSeqMeta(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const string& out_file) {
  double sigma_sq = 1.0;
  int num_samples = 0;
  if (!GetNumSamplesFromMassFile(mass_file_info, &num_samples)) {
    return false;
  }

  // Read Mass input file.
  int num_snps = 0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  // TODO(PHB): Determine if this is the proper time to do 'rescale'ing.
  if (!GetInfoFromMassFile(
          false, study_num, allele_info, mass_file_info, snp_to_excluding_study_,
          &num_snps, &gene_to_snp_positions, &snps_to_skip, &monomorphic_snps,
          &snp_info, &covariances)) {
    cout << "ERROR in converting Mass to SeqMeta: Unable to parse MASS score "
         << "file. Aborting.\n";
    return false;
  }

  // Convert info into a SeqMeta RDataTree.
  vector<string*>* tags = new vector<string*>();
  RDataNode* root = new RDataNode();
  const string obj_name = GetSeqMetaObjectNameFromFile(out_file);
  if (!ConstructSeqMetaTree(
          obj_name, study_num, num_samples, rescale, sigma_sq, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions, snp_info,
          covariances, tags, root)) {
    cout << "ERROR in printing output for Mass to SeqMeta. Aborting.\n";
    return false;
  }

  // Print SeqMeta .Rdata file.
  if (!PrintToSeqMeta(*root, out_file)) {
    cout << "ERROR in Printing output for Mass to SeqMeta. Aborting.\n";
    return false;
  }

  // Delete the RDataTree.
  set<string*> sexptype_tags;
  for (string* tag : *tags) sexptype_tags.insert(tag);
  DeleteRDataTree(root, sexptype_tags);
  DeleteRDataTags(*tags);
  delete tags;

  return true;
}

bool PreMeta::MetaSkatToMass(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& meta_skat_m_info,
    const string& meta_skat_mssd_file,
    const string& mass_script_file, const string& out_file) {
  // Parse MInfo File.
  int num_samples;
  map<string, vector<Position>> gene_to_snp_positions;
  map<string, uint64_t> start_lines;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!ParseMInfoFile(
          study_num, meta_skat_m_info, allele_info,
          &num_samples, &monomorphic_snps, &snp_to_ref_alt_and_study_,
          &snp_to_excluding_study_, &gene_to_snp_positions,
          &start_lines, &snp_info)) {
    return false;
  }

  // Parse MSSD File.
  map<pair<Position, Position>, double> covariances;
  if (!ParseMssdFile(
          study_num, meta_skat_mssd_file, snp_to_excluding_study_,
          gene_to_snp_positions, start_lines, &covariances)) {
    return false;
  }

  // Print Mass script file.
  if (!PrintMassScriptFile(out_file, mass_script_file)) return false;

  const double sigma_sq = 1.0;
  return PrintToMass(study_num, num_samples, rescale, sigma_sq,
                     out_file, monomorphic_snps, snp_to_ref_alt_and_study_,
                     gene_to_snp_positions, snp_info, covariances);
}

bool PreMeta::MetaSkatToMetaSkat(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& meta_skat_m_info,
    const string& meta_skat_mssd_file,
    const string& minfo_outfile, const string& mssd_outfile) {
  // Parse MInfo File.
  int num_samples;
  map<string, vector<Position>> gene_to_snp_positions;
  map<string, uint64_t> start_lines;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!ParseMInfoFile(
          study_num, meta_skat_m_info, allele_info,
          &num_samples, &monomorphic_snps, &snp_to_ref_alt_and_study_,
          &snp_to_excluding_study_, &gene_to_snp_positions,
          &start_lines, &snp_info)) {
    return false;
  }

  // Parse MSSD File.
  map<pair<Position, Position>, double> covariances;
  if (!ParseMssdFile(
          study_num, meta_skat_mssd_file, snp_to_excluding_study_,
          gene_to_snp_positions, start_lines, &covariances)) {
    return false;
  }

  // Write MInfo and MSSD files.
  if (!WriteMInfoFile(
          study_num, num_samples, snp_info.size(), rescale, 1.0, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions, snp_info, minfo_outfile)) {
    cout << "ERROR in converting Mass to MetaSkat: Failure in writing "
         << ".MInfo file. Aborting.\n";
    return false;
  }
  
  if (!WriteMssdFile(
          study_num, rescale, monomorphic_snps, snp_to_ref_alt_and_study_,
          gene_to_snp_positions, covariances, mssd_outfile)) {
    cout << "ERROR in converting Mass to MetaSkat: Failure in writing "
         << ".MSSD file. Aborting.\n";
    return false;
  }

  return true;
}

bool PreMeta::MetaSkatToRareMetal(
    const int study_num, const int window_size, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& meta_skat_m_info,
    const string& meta_skat_mssd_file,
    const string& score_files, const string& cov_files,
    const string& out_score_file, const string& out_cov_file) {
  // Parse MInfo File.
  int num_samples;
  map<string, vector<Position>> gene_to_snp_positions;
  map<string, uint64_t> start_lines;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!ParseMInfoFile(
          study_num, meta_skat_m_info, allele_info,
          &num_samples, &monomorphic_snps, &snp_to_ref_alt_and_study_,
          &snp_to_excluding_study_, &gene_to_snp_positions,
          &start_lines, &snp_info)) {
    return false;
  }

  // Parse MSSD File.
  map<pair<Position, Position>, double> covariances;
  if (!ParseMssdFile(
          study_num, meta_skat_mssd_file, snp_to_excluding_study_,
          gene_to_snp_positions, start_lines, &covariances)) {
    return false;
  }

  // Compute Windows.
  const int window = window_size == -1 ? kRareMetalWindowSize : window_size;
  map<Position, Position> windows;
  if (!GetWindows(window, Keys(snp_info), &windows)) {
    cout << "ERROR in converting MetaSkat to RareMetal: Unable to get windows. "
         << "Aborting.\n";
    return false;
  }

  // Print RareMetal Score and Covariance files.
  // MetaSkat doesn't provide sigma_sq, so we use a value of 1.0, which
  // is the best we can do. Note that if printing to the old version of
  // R, then U_STAT and V_STAT are going to be off by a factor of
  // sigma_sq (and SE will be off by sqrt(sigma_sq)).
  if (!PrintToRareMetal(
          study_num, num_samples, rescale, 1.0, 1.0,
          monomorphic_snps, snp_to_ref_alt_and_study_, windows, snp_info,
          covariances, out_score_file, out_cov_file)) {
    cout << "ERROR in printing MetaSkat to RareMetal. Aborting.\n";
    return false;
  }

  // Add above files to Score and Covariance files list.
  return AddFilesToList(score_files, cov_files, out_score_file, out_cov_file);
}

bool PreMeta::MetaSkatToSeqMeta(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& meta_skat_m_info,
    const string& meta_skat_mssd_file,
    const string& out_file) {
  // Parse MInfo File.
  int num_samples;
  map<string, vector<Position>> gene_to_snp_positions;
  map<string, uint64_t> start_lines;
  map<Position, SnpInfo> snp_info;
  set<Position> monomorphic_snps;
  if (!ParseMInfoFile(
          study_num, meta_skat_m_info, allele_info,
          &num_samples, &monomorphic_snps, &snp_to_ref_alt_and_study_,
          &snp_to_excluding_study_, &gene_to_snp_positions,
          &start_lines, &snp_info)) {
    return false;
  }

  // Parse MSSD File.
  map<pair<Position, Position>, double> covariances;
  if (!ParseMssdFile(
          study_num, meta_skat_mssd_file, snp_to_excluding_study_,
          gene_to_snp_positions, start_lines, &covariances)) {
    return false;
  }


  // Convert info into a SeqMeta RDataTree.
  const double sigma_sq = 1.0;
  vector<string*>* tags = new vector<string*>();
  RDataNode* root = new RDataNode();
  const string obj_name = GetSeqMetaObjectNameFromFile(out_file);
  if (!ConstructSeqMetaTree(
          obj_name, study_num, num_samples, rescale, sigma_sq,
          monomorphic_snps, snp_to_ref_alt_and_study_,
          gene_to_snp_positions, snp_info, covariances, tags, root)) {
    cout << "ERROR in printing MetaSKAT to SeqMeta. Aborting.\n";
    return false;
  }

  // Print SeqMeta .Rdata file.
  if (!PrintToSeqMeta(*root, out_file)) {
    cout << "ERROR in Printing MetaSKAT to SeqMeta. Aborting.\n";
    return false;
  }

  // Delete the RDataTree.
  set<string*> sexptype_tags;
  for (string* tag : *tags) sexptype_tags.insert(tag);
  DeleteRDataTree(root, sexptype_tags);
  DeleteRDataTags(*tags);
  delete tags;

  return true;
}


bool PreMeta::MassToMetaSkat(
    const int study_num, const double& rescale,
    const map<Position, SnpInfo>& allele_info,
    const FileInfo& mass_file_info,
    const string& minfo_outfile, const string& mssd_outfile) {
  // Read Residual Variance from Mass input file.
  double sigma_sq = 1.0;
  int num_samples = 0;
  if (!GetNumSamplesFromMassFile(mass_file_info, &num_samples)) {
    return false;
  }

  // Read Mass input file.
  int num_snps = 0;
  map<string, vector<Position>> gene_to_snp_positions;
  map<Position, SnpInfo> snp_info;
  map<pair<Position, Position>, double> covariances;
  set<Position> snps_to_skip;
  set<Position> monomorphic_snps;
  // TODO(PHB): Determine if this is the proper time to do 'rescale'ing.
  if (!GetInfoFromMassFile(
          true, study_num, allele_info, mass_file_info, snp_to_excluding_study_,
          &num_snps, &gene_to_snp_positions, &snps_to_skip, &monomorphic_snps,
          &snp_info, &covariances)) {
    cout << "ERROR in converting Mass to MetaSkat: Unable to parse MASS score "
         << "file. Aborting.\n";
    return false;
  }

  // Read Allele file.
  if (!CopyAlleleInfoToSnpInfo(allele_info, &snp_info)) {
    cout << "ERROR in converting Mass to MetaSkat: Failed to parse allele "
         << "info from annotation file. Aborting.\n";
    return false;
  }
  
  // Print to file Snps that will be skipped due to no Allele Info.
  PrintSkippedSnps(minfo_outfile, snps_to_skip);

  if (!WriteMInfoFile(
          study_num, num_samples, num_snps, rescale, sigma_sq, snps_to_skip,
          monomorphic_snps, snp_to_ref_alt_and_study_, gene_to_snp_positions,
          snp_info, minfo_outfile)) {
    cout << "ERROR in converting Mass to MetaSkat: Failure in writing "
         << ".MInfo file. Aborting.\n";
    return false;
  }
  
  if (!WriteMssdFile(
          study_num, rescale, snps_to_skip, monomorphic_snps,
          snp_to_ref_alt_and_study_, gene_to_snp_positions,
          covariances, mssd_outfile)) {
    cout << "ERROR in converting Mass to MetaSkat: Failure in writing "
         << ".MSSD file. Aborting.\n";
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
// Test: Create binary .rdata file to be read by R.
// -----------------------------------------------------------------------------
/*
// Input:
//   type: R object type, or 'SEXPTYPE'; See:
//         http://cran.r-project.org/doc/manuals/r-release/R-ints.html#SEXPs
int pack_flags(int type, int levels, int is_object, int has_attr, int has_tag) {
  int flags;

  if (type == 9) {
    // scalar string type, used for symbol names
    levels &= (~((1 << 5) | 1));
  }

  flags = type | (levels << 12);
  if (is_object)
    flags |= (1 << 8);
  if (has_attr)
    flags |= (1 << 9);
  if (has_tag)
    flags |= (1 << 10);

  return flags;
}

void encode_integer(int i, char* buf) {
  XDR xdrs;
  int success;

  xdrmem_create(&xdrs, buf, kRXDRIntegerBytes, XDR_ENCODE);
  success = xdr_int(&xdrs, &i);
  xdr_destroy(&xdrs);
  if (!success) {
    printf("encode_integer failed\n");
    exit(1);
  }
}

void encode_double(double d, char* buf) {
  XDR xdrs;
  int success;

  xdrmem_create(&xdrs, buf, kRXDRDoubleBytes, XDR_ENCODE);
  success = xdr_double(&xdrs, &d);
  xdr_destroy(&xdrs);
  if (!success) {
    printf("encode_double failed\n");
    exit(1);
  }
}

void write_data(const char* buf, int len, FILE* fp) {
  int res;
  res = fwrite(buf, sizeof(char), len, fp);
  if (res != len) {
    printf("Write failed\n");
    exit(1);
  }
}

void write_integer(int i, char* buf, FILE* fp) {
  encode_integer(i, buf);
  write_data(buf, kRXDRIntegerBytes, fp);
}

void write_double(double d, char* buf, FILE* fp) {
  encode_double(d, buf);
  write_data(buf, kRXDRDoubleBytes, fp);
}

bool PreMeta::WriteSeqMetaBinary(const string& outfile) {
  FILE* fp;
  int res;
  char buf[128];

  fp = fopen(outfile.c_str(), "w");
  if (fp == NULL) {
    printf("Couldn't open file for writing\n");
    return false;
  }

  // Write magic: XDR_V2
  char xdr_version[6] = "RDX2\n";
  write_data(xdr_version, 5, fp);

  // Write format
  char xdr_format[3] = "X\n";
  write_data(xdr_format, 2, fp);

  // Write R version information
  write_integer(2, buf, fp);       // serialization version: 2
  write_integer(133633, buf, fp);  // Current R version (2.10.1 in this case)
  write_integer(131840, buf, fp);  // Version number for R 2.3.0

  // The saved R objects are wrapped in a list of dotted pairs before saving.
  // Next we write out flags needed for this list.
  write_integer(pack_flags(2, 0, 0, 0, 1), buf, fp);

  // Write the name of the variable we're storing
  write_integer(1, buf, fp);                          // symbol type
  // When type (first parameter) is 9, we will AND level (second parameter)
  // with ~33 = 11011110. So making 2nd parameter 33 makes all of this a no-op.
  // So, given that type is 9, setting level (2nd parameter) to any of the
  // following are equivalent: 0, 1, 32, 33. Similarly, if type is 9, then
  // for any level X, flipping the 1st and/or 6th bits of X has no effect.
  write_integer(pack_flags(9, 0, 0, 0, 0), buf, fp); // symbol flags
  //write_integer(pack_flags(9, 33, 0, 0, 0), buf, fp); // symbol flags
  int var_name_length;
  const char* var_name =
      StringToCharPointer("foo_foo", &var_name_length);
      //StringToCharPointer("seqMeta_" + outfile, &var_name_length);
  write_integer(var_name_length, buf, fp);
  write_data(var_name, var_name_length, fp);

  // Now write the actual variable data
  write_integer(pack_flags(14, 0, 0, 0, 0), buf, fp); // vector of reals
  write_integer(3, buf, fp);   // length of vector
  write_double(1.1, buf, fp);  // first value
  write_double(2.2, buf, fp);  // second value
  write_double(-3.1, buf, fp);  // third value

  // Tell R we're done
  write_integer(254, buf, fp);
  fclose(fp);

  return true;
}

// Input:
//   type: R object type, or 'SEXPTYPE'; See:
//         http://cran.r-project.org/doc/manuals/r-release/R-ints.html#SEXPs
bool unpack_flags(
    int flags,
    int* type, int* levels, bool* is_object, bool* has_attr, bool* has_tag) {
//int unpack_flags(int type, int levels, int is_object, int has_attr, int has_tag) {
  cout << "\nPHB flags: " << flags << endl;
  *type = flags & 31;
  *is_object = (flags >> 8 & 1);
  *has_attr = (flags >> 9 & 1);
  *has_tag = (flags >> 10 & 1);
  if (*type == 9) {
    *levels = (flags >> 12);
  } else {
    *levels = (flags >> 12);
  }

  if (type == 9) {
    // scalar string type, used for symbol names
    levels &= (~((1 << 5) | 1));
  }

  flags = type | (levels << 12);
  if (is_object)
    flags |= (1 << 8);
  if (has_attr)
    flags |= (1 << 9);
  if (has_tag)
    flags |= (1 << 10);

  return flags;
}

void decode_integer(int i, char* buf) {
  XDR xdrs;
  int success;

  xdrmem_create(&xdrs, buf, kRXDRIntegerBytes, XDR_ENCODE);
  success = xdr_int(&xdrs, &i);
  xdr_destroy(&xdrs);
  if (!success) {
    printf("encode_integer failed\n");
    exit(1);
  }
}

void decode_double(double d, char* buf) {
  XDR xdrs;
  int success;

  xdrmem_create(&xdrs, buf, kRXDRDoubleBytes, XDR_ENCODE);
  success = xdr_double(&xdrs, &d);
  xdr_destroy(&xdrs);
  if (!success) {
    printf("encode_double failed\n");
    exit(1);
  }
}

void read_data(FILE* fp, char* buf) {
  int res;
  res = fwrite(buf, sizeof(char), len, fp);
  if (res != len) {
    printf("Write failed\n");
    exit(1);
  }
}
*/

}  // namespace premeta
