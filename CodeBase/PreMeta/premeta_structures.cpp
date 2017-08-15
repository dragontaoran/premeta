// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "premeta_structures.h"

#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "StringUtils/string_utils.h"

#include <fstream>
#include <iostream>
#include <map>

using file_reader_utils::VcfUtils;
using namespace string_utils;
using namespace map_utils;

namespace premeta {

bool ParsePosition(const string& input, Chromosome* chr, uint64_t* pos) {
  vector<string> pos_parts;
  Split(input, ":", &pos_parts);
  if (pos_parts.size() != 2) return false;

  // Parse Chromosome.
  const string& chr_str = pos_parts[0];
  if (!VcfUtils::ParseChromosome(chr_str, chr)) return false;

  // Parse Position.
  if (!Stoi(pos_parts[1], pos)) return false;

  return true;
}

bool ParsePosition(const string& input, Position* pos) {
  return ParsePosition(input, &(pos->chr_), &(pos->pos_));
}

string PrintPosition(const Position& pos) {
  // See if we have a valid CHR:POS to print; otherwise, print snp_id_.
  if (pos.chr_ == Chromosome::CHROMOSOME_UNKNOWN && pos.pos_ == 0) {
    return pos.snp_id_;
  }
  return
      VcfUtils::PrintChromosome(pos.chr_) + ":" + Itoa(pos.pos_);
}

bool ParseSoftwareVersion(const string& input, SoftwareVersion* version) {
  if (version == nullptr) return false;
  
  version->base_ = -1;
  version->sub_version_ = -1;
  version->sub_sub_version_ = -1;

  string no_whitespace;
  RemoveAllWhitespace(input, &no_whitespace);
  if (no_whitespace.empty()) return false;
  vector<string> parts;
  Split(no_whitespace, ".", &parts);
  if (parts.size() > 3) return false;

  for (int i = 0; i < parts.size(); ++i) {
    int version_num;
    if (!Stoi(parts[i], &version_num)) return false;
    if (i == 0) version->base_ = version_num;
    if (i == 1) version->sub_version_ = version_num;
    if (i == 2) version->sub_sub_version_ = version_num;
  }

  return true;
}

string PrintSoftwareVersion(const DefaultSoftwareVersion version) {
  if (version == SOFTWARE_VERSION_MASS) {
    return "7.0";
  }
  if (version == SOFTWARE_VERSION_METASKAT) {
    return "0.40";
  }
  if (version == SOFTWARE_VERSION_RAREMETAL_NEW) {
    return "4.13.5";
  }
  if (version == SOFTWARE_VERSION_RAREMETAL_OLD) {
    return "0.4.0";
  }
  if (version == SOFTWARE_VERSION_SEQMETA) {
    return "1.5";
  }
  return "Unknwn";
}

string PrintSoftwareVersion(const SoftwareVersion& version) {
  if (version.base_ < 0) return "Unknown";
  string to_return = "";
  to_return += Itoa(version.base_);
  if (version.sub_version_ < 0) return to_return;
  to_return += "." + Itoa(version.sub_version_);
  if (version.sub_sub_version_ < 0) return to_return;
  to_return += "." + Itoa(version.sub_sub_version_);
  return to_return;
}

bool GetColumnIndices(
    const FileInfo& file_info, map<string, int>* column_indices) {
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR. Unable to find file '" << file_info.name_
         << "'. Aborting." << endl;
    return false;
  }
  string line = "";
  int line_num = 0;
  bool should_break = false;
  while (getline(file, line) && !should_break) {
    line_num++;

    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    if (!HasPrefixString(line, file_info.comment_char_)) {
      // For most use-cases, we could just 'break' here. But MetaSKAT's .MInfo
      // file actually has the header on the first non-comment row. So we set
      // 'should_break' to true, so that this will be the last line processed.
      should_break = true;
    }

    line = StripPrefixString(line, file_info.comment_char_);

    vector<string> cols;
    Split(line, file_info.delimiter_, &cols);
    if (cols.size() < column_indices->size()) continue;
    int num_matches = 0;
    for (int i = 0; i < cols.size(); ++i) {
      const string& col = cols[i];
      map<string, int>::iterator col_itr = column_indices->find(col);
      if (col_itr == column_indices->end()) continue;
      col_itr->second = (i + 1);
      num_matches++;
    }
    if (num_matches == column_indices->size()) {
      return true;
    }
    if (num_matches > column_indices->size()) {
      continue;
    }
    if (num_matches > 0 && num_matches < column_indices->size()) {
      continue;
    }
  }

  return false;
}

bool LookupCovariance(
    const Position& pos_one, const Position& pos_two,
    const map<pair<Position, Position>, double>& covariances,
    double* covariance) {
  // First check if (pos_one, pos_two) is a Key in covariances.
  map<pair<Position, Position>, double>::const_iterator orientation_one_itr =
      covariances.find(make_pair(pos_one, pos_two));
  const bool found_pair = orientation_one_itr != covariances.end();

  if (found_pair) {
    *covariance = orientation_one_itr->second;
  }

  // If pos_one == pos_two, no need to check opposite orientation.
  if (pos_one.chr_ == pos_two.chr_ && pos_one.pos_ == pos_two.pos_ &&
      pos_one.snp_id_ == pos_two.snp_id_) {
    if (!found_pair) {
      // In some places (e.g. ComputeCovarianceByGene, where this will
      // always happen), we are expected to hit this. So don't print out a
      // WARNING here; instead, all LookupCovariance failure messages
      // are handled by the calling function.
    }
    return found_pair;
  }

  // Check the other orientation.
  map<pair<Position, Position>, double>::const_iterator orientation_two_itr =
      covariances.find(make_pair(pos_two, pos_one));
  const bool found_pair_two = orientation_two_itr != covariances.end();

  if (found_pair_two) {
    if (found_pair) {
      if (*covariance != orientation_two_itr->second) {
        cout << "ERROR in Looking up Covariance: two different covariance "
             << "values were found for Positions " << PrintPosition(pos_one)
             << " and " << PrintPosition(pos_two) << ". Aborting." << endl;
        return false;
      }
    } else {
      *covariance = orientation_two_itr->second;
    }
  }

  return (found_pair || found_pair_two);
}

bool ShouldExcludeSnp(
    const int study_num, const Position& pos,
    const map<Position, set<int>>& snp_to_excluding_study) {
  map<Position, set<int>>::const_iterator exclude_snp_itr =
      snp_to_excluding_study.find(pos);
  if (exclude_snp_itr == snp_to_excluding_study.end()) return false;
  const set<int>& excluding_studies = exclude_snp_itr->second;
  return excluding_studies.find(study_num) != excluding_studies.end();
}

bool IsSnpRefAltSwapped(
    const int study_num, const Position& pos,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study) {
  map<Position, tuple<string, string, set<int>>>::const_iterator swapped_itr =
      snp_to_ref_alt_and_study.find(pos);
  if (swapped_itr == snp_to_ref_alt_and_study.end()) return false;

  const set<int>& swapped_studies = get<2>(swapped_itr->second);
  return swapped_studies.find(study_num) != swapped_studies.end();
}

bool CopyAlleleInfoToSnpInfo(
    const map<Position, SnpInfo>& allele_info,
    map<Position, SnpInfo>* snp_info) {
  for (const pair<Position, SnpInfo>& input_snp : allele_info) {
    const Position& meta_pos = input_snp.second.pos_;
    map<Position, SnpInfo>::iterator snp_itr = snp_info->find(meta_pos);
    if (snp_itr == snp_info->end()) {
      // NOTE: Updated this 11/18 so that new SNPs are not added to SNP_INFO.
      // This was necessary because sometimes the SNP_INFO file may have
      // information on SNPs that are not in a given study, and then
      // write_raremetal_utils.cpp loops through the SNPs in snp_info, and
      // complains if it can't find score information for that SNP (which
      // would happen in the present example, since SNPs in the SNP_INFO
      // file that aren't in the study used to get inserted into snp_info
      // here, and obviously they won't have score information.
      // Thus, the below code was commented out 11/18.
      /*
      SnpInfo current_snp_info;
      current_snp_info.major_allele_ = input_snp.second.major_allele_;
      current_snp_info.minor_allele_ = input_snp.second.minor_allele_;
      snp_info->insert(make_pair(meta_pos, current_snp_info));
      */
    } else {
      SnpInfo& current_snp_info = snp_itr->second;
      current_snp_info.major_allele_ = input_snp.second.major_allele_;
      current_snp_info.minor_allele_ = input_snp.second.minor_allele_;
    }
  }

  return true;
}

}  // namespace premeta
