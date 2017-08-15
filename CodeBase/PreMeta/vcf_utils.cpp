// Author: paulbunn@email.unc.edu (Paul Bunn)
// Last Updated: March 2015

#include "vcf_utils.h"

#include "StringUtils/string_utils.h"

#include <ctime>     // For printing out current time
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using string_utils::StringUtils;

namespace file_reader_utils {

static bool kInitializeGtSeparatorsOnce = true;
const int64_t MAX_BUFFER_SIZE = 2147483648; // 2 GiB

bool VcfUtils::IsSegmentDone(const ProcessLineResponse& response_type,
                             const int num_variants_to_keep,
                             const int num_variants_taken,
                             const OutputOptions& output_options) {
  if (output_options.segment_type_ != SegmentType::GENE) {
    if (num_variants_taken >= num_variants_to_keep) {
      return true;
    }
    return false;
  } else if (output_options.segment_type_ == SegmentType::GENE) {
    if (num_variants_taken > 0 &&
        response_type == SNP_NOT_IN_SEGMENT) {
      // Since SNPs that are part of a gene will be a contiguous block
      // of lines in the vcf file, and since we processed at least one
      // SNP in the gene already, we must have reached the end of this
      // contiguous block.
      return true;
    }
    return false;
  } else {
    // Unexpected SegmentType. Print warning, and return false.
    cout << "ERROR: Unexpected SegmentType: "
         << output_options.segment_type_
         << ". Ignoring OutputOptions SegmentType (process whole file). "
         << endl;
  }
  return false;
}

ProcessLineResponse VcfUtils::ProcessLine(
    const string& line,
    const int expected_num_cols,
    const vector<string>& header_row,
    const GroupingFilter& grouping_filter,
    const FilteringOptions& filtering_options,
    const SnpOutputType& output_type,
    FilteringCounts* filtering_counts,
    vector<SnpType>* snp_values_per_sample,
    vector<SnpInfo>* snp_info) {
  // Split line into its columns.
  vector<string> columns;
  StringUtils::Split(
      line, "\t", false /* Do Not Skip/Flatten Empty Columns */, &columns);
  if (expected_num_cols != columns.size()) {
    cout << "ERROR Reading line: number of columns (" << columns.size()
         << ") does not match number of columns in header ("
         << expected_num_cols << "). Line:\n" << line << endl;
    return ProcessLineResponse::FAILURE;
  }

  // Create the item that will hold all info for this Variant; populate
  // its fields as we process this line.
  SnpInfo info;

  // Process Chromesome Number.
  if (!ParseChromosome(columns[0], &info.chromosome_)) {
    cout << "ERROR: Unable to get chromosome number from line:\n"
         << line << endl;
    return ProcessLineResponse::FAILURE;
  }

  // Process Position.
  if (!GetPosition(columns[1], &info.position_)) {
    cout << "ERROR: Unable to get position from line:\n"
         << line << endl;
    return ProcessLineResponse::FAILURE;
  }

  if (IsPositionFiltered(info.chromosome_, info.position_,
                         filtering_options.position_filter_)) {
    filtering_counts->num_position_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Process SNP ID.
  info.id_ = columns[2];

  // To save space, don't store default SNP ID ".".
  if (info.id_ == ".") info.id_ = "";

  // Filter by SNP ID, if appropriate.
  if (IsSnpFiltered(info.id_, filtering_options.snps_to_discard_,
                    filtering_options.snps_to_keep_)) {
    filtering_counts->num_snp_id_column_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Apply GroupingFilter to see if this Variant belongs int the output.
  if (!grouping_filter.variants_by_position_.empty()) {
    const map<Chromosome, set<int64_t>>::const_iterator grouping_itr =
        grouping_filter.variants_by_position_.find(info.chromosome_);
    if (grouping_itr == grouping_filter.variants_by_position_.end()) {
      return ProcessLineResponse::SNP_NOT_IN_SEGMENT;
    }
    const set<int64_t>& variants_to_keep = grouping_itr->second;
    if (variants_to_keep.find(info.position_) == variants_to_keep.end()) {
      return ProcessLineResponse::SNP_NOT_IN_SEGMENT;
    }
  } else if (get<0>(grouping_filter.gene_position_) !=
             Chromosome::CHROMOSOME_UNKNOWN) {
    if (get<0>(grouping_filter.gene_position_) != info.chromosome_ ||
        get<1>(grouping_filter.gene_position_) > info.position_ ||
        get<2>(grouping_filter.gene_position_) < info.position_) {
      return ProcessLineResponse::SNP_NOT_IN_SEGMENT;
    }
  } else if (!grouping_filter.snp_ids_.empty()) {
    if (grouping_filter.snp_ids_.find(info.id_) ==
        grouping_filter.snp_ids_.end()) {
      return ProcessLineResponse::SNP_NOT_IN_SEGMENT;
    }
  } else {
    cout << "ERROR: Exactly one GroupingFilter field should be set. "
         << "Aborting." << endl;
    return ProcessLineResponse::FAILURE;
  }

  // Process Major and Minor Alleles.
  if(!GetAlleleTypes(
         columns[3], columns[4],
         &info.major_allele_, &info.minor_allele_, &info.triallele_)) {
    return ProcessLineResponse::FAILURE;
  }

  // Filter by Quality, if appropriate.
  if (IsQualityFiltered(columns[5], filtering_options.min_quality_threshold_)) {
    filtering_counts->num_low_quality_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Filter by Filter Column, if appropriate.
  if (IsFilterFiltered(columns[6], filtering_options.filter_column_tags_)) {
    filtering_counts->num_filter_column_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Apply INFO filters, if appropriate.
  if (IsInfoFiltered(columns[7], filtering_options.info_column_filters_,
                     &info.minor_allele_freq_, &info.triallele_freq_)) {
    filtering_counts->num_info_column_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Apply Allele filters, if appropriate.
  if (IsAlleleFiltered(info.minor_allele_, info.triallele_,
                       info.minor_allele_freq_, info.triallele_freq_,
                       filtering_options.allele_filter_)) { 
    filtering_counts->num_allele_type_exclusions_++;
    return ProcessLineResponse::FILTERED;
  }

  // Determine which entry in data columns corresponds to the genotype.
  // TODO(PHB): Determine if any of the other data (besides GT) is could
  // be of interest, and if so, generalize below to search of other
  // possibly interesting data types.
  const string& format = columns[8];
  vector<string> format_cols;
  StringUtils::Split(format, ":", &format_cols);
  int gt_index = -1;
  for (int i = 0; i < format_cols.size(); ++i) {
    if (format_cols[i] == "GT") {
      gt_index = i;
      break;
    }
  }
  if (gt_index < 0) {
    cout << "ERROR: Row does not have GT info: " << format << endl;
    return ProcessLineResponse::FAILURE;
  }

  // Update num_sample_id_exclusions_ only if this is the first row
  // for this block (i.e. we don't need to count every row, since
  // they will all be the same; just count for the first row).
  bool should_update_num_samples_filtered =
      filtering_counts->num_sample_id_exclusions_ == 0;
  // Process the data.
  for (int i = 9; i < expected_num_cols; ++i) {
    if (!filtering_options.samples_to_keep_.empty() ||
        !filtering_options.samples_to_discard_.empty()) {
      if (!filtering_options.samples_to_keep_.empty()) {
        const set<string>& to_keep = filtering_options.samples_to_keep_;
        if (to_keep.find(header_row[i]) == to_keep.end()) {
          if (should_update_num_samples_filtered) {
            filtering_counts->num_sample_id_exclusions_++;
          }
          continue;
        }
      } else {
        const set<string>& to_discard = filtering_options.samples_to_discard_;
        if (to_discard.find(header_row[i]) != to_discard.end()) {
          if (should_update_num_samples_filtered) {
            filtering_counts->num_sample_id_exclusions_++;
          }
          continue;
        }
      }
    }
    snp_values_per_sample->push_back(SnpType());
    if (!ProcessDataColumn(
            columns[i], gt_index, filtering_options.actual_data_filters_,
            output_type, filtering_counts, &(snp_values_per_sample->back()))) {
      cout << "ERROR: Unable to process data column " << i
           << ": " << columns[i] << endl;
      return ProcessLineResponse::FAILURE;
    }
  }

  // Add SnpInfo.
  snp_info->push_back(info);

  return ProcessLineResponse::SNP_IN_SEGMENT;
}

bool VcfUtils::ParseChromosome(const string& chr_str, Chromosome* chr) {
  if (chr_str == "X") *chr = Chromosome::CHROMOSOME_X;
  else if (chr_str == "Y") *chr = Chromosome::CHROMOSOME_Y;
  else if (chr_str == "1") *chr = Chromosome::CHROMOSOME_ONE;
  else if (chr_str == "2") *chr = Chromosome::CHROMOSOME_TWO;
  else if (chr_str == "3") *chr = Chromosome::CHROMOSOME_THREE;
  else if (chr_str == "4") *chr = Chromosome::CHROMOSOME_FOUR;
  else if (chr_str == "5") *chr = Chromosome::CHROMOSOME_FIVE;
  else if (chr_str == "6") *chr = Chromosome::CHROMOSOME_SIX;
  else if (chr_str == "7") *chr = Chromosome::CHROMOSOME_SEVEN;
  else if (chr_str == "8") *chr = Chromosome::CHROMOSOME_EIGHT;
  else if (chr_str == "9") *chr = Chromosome::CHROMOSOME_NINE;
  else if (chr_str == "10") *chr = Chromosome::CHROMOSOME_TEN;
  else if (chr_str == "11") *chr = Chromosome::CHROMOSOME_ELEVEN;
  else if (chr_str == "12") *chr = Chromosome::CHROMOSOME_TWELVE;
  else if (chr_str == "13") *chr = Chromosome::CHROMOSOME_THIRTEEN;
  else if (chr_str == "14") *chr = Chromosome::CHROMOSOME_FOURTEEN;
  else if (chr_str == "15") *chr = Chromosome::CHROMOSOME_FIFTEEN;
  else if (chr_str == "16") *chr = Chromosome::CHROMOSOME_SIXTEEN;
  else if (chr_str == "17") *chr = Chromosome::CHROMOSOME_SEVENTEEN;
  else if (chr_str == "18") *chr = Chromosome::CHROMOSOME_EIGHTEEN;
  else if (chr_str == "19") *chr = Chromosome::CHROMOSOME_NINETEEN;
  else if (chr_str == "20") *chr = Chromosome::CHROMOSOME_TWENTY;
  else if (chr_str == "21") *chr = Chromosome::CHROMOSOME_TWENTYONE;
  else if (chr_str == "22") *chr = Chromosome::CHROMOSOME_TWENTYTWO;
  else {
    *chr = Chromosome::CHROMOSOME_UNKNOWN;
    return false;
  }
  return true;
}

string VcfUtils::PrintChromosome(const Chromosome chr) {
  switch (chr) {
    case Chromosome::CHROMOSOME_UNKNOWN:
    {
      return "Unknown";
    }
    case Chromosome::CHROMOSOME_ONE:
    {
      return "1";
    }
    case Chromosome::CHROMOSOME_TWO:
    {
      return "2";
    }
    case Chromosome::CHROMOSOME_THREE:
    {
      return "3";
    }
    case Chromosome::CHROMOSOME_FOUR:
    {
      return "4";
    }
    case Chromosome::CHROMOSOME_FIVE:
    {
      return "5";
    }
    case Chromosome::CHROMOSOME_SIX:
    {
      return "6";
    }
    case Chromosome::CHROMOSOME_SEVEN:
    {
      return "7";
    }
    case Chromosome::CHROMOSOME_EIGHT:
    {
      return "8";
    }
    case Chromosome::CHROMOSOME_NINE:
    {
      return "9";
    }
    case Chromosome::CHROMOSOME_TEN:
    {
      return "10";
    }
    case Chromosome::CHROMOSOME_ELEVEN:
    {
      return "11";
    }
    case Chromosome::CHROMOSOME_TWELVE:
    {
      return "12";
    }
    case Chromosome::CHROMOSOME_THIRTEEN:
    {
      return "13";
    }
    case Chromosome::CHROMOSOME_FOURTEEN:
    {
      return "14";
    }
    case Chromosome::CHROMOSOME_FIFTEEN:
    {
      return "15";
    }
    case Chromosome::CHROMOSOME_SIXTEEN:
    {
      return "16";
    }
    case Chromosome::CHROMOSOME_SEVENTEEN:
    {
      return "17";
    }
    case Chromosome::CHROMOSOME_EIGHTEEN:
    {
      return "18";
    }
    case Chromosome::CHROMOSOME_NINETEEN:
    {
      return "19";
    }
    case Chromosome::CHROMOSOME_TWENTY:
    {
      return "20";
    }
    case Chromosome::CHROMOSOME_TWENTYONE:
    {
      return "21";
    }
    case Chromosome::CHROMOSOME_TWENTYTWO:
    {
      return "22";
    }
    case Chromosome::CHROMOSOME_X:
    {
      return "X";
    }
    case Chromosome::CHROMOSOME_Y:
    {
      return "Y";
    }
    default:
      // Unrecognized/Unsupported Chromosome type.
      return "";
  }
  // Shouldn't reach here.
  return "";
}

string VcfUtils::GetNucleotideString(const Nucleotide n) {
  switch (n) {
    case Nucleotide::NUCLEOTIDE_UNKNOWN:
    {
      return "Unknown";
    }
    case Nucleotide::NUCLEOTIDE_NA:
    {
      return "NA";
    }
    case Nucleotide::A:
    {
      return "A";
    }
    case Nucleotide::C:
    {
      return "C";
    }
    case Nucleotide::G:
    {
      return "G";
    }
    case Nucleotide::T:
    {
      return "T";
    }
    case Nucleotide::U:
    {
      return "U";
    }
    case Nucleotide::R:
    {
      return "R";
    }
    case Nucleotide::Y:
    {
      return "Y";
    }
    case Nucleotide::W:
    {
      return "W";
    }
    case Nucleotide::S:
    {
      return "S";
    }
    case Nucleotide::M:
    {
      return "M";
    }
    case Nucleotide::K:
    {
      return "K";
    }
    case Nucleotide::N:
    {
      return "N";
    }
    case Nucleotide::NUCLEOTIDE_INSERTION:
    {
      return "Insertion";
    }
    case Nucleotide::NUCLEOTIDE_DELETION:
    {
      return "Deletion";
    }
    case Nucleotide::NUCLEOTIDE_SUBSTITUTION:
    {
      return "Substitution";
    }
    case Nucleotide::NUCLEOTIDE_MULTI:
    {
      return "Sequence";
    }
    case Nucleotide::NUCLEOTIDE_NONE:
    {
      return "None";
    }
    case Nucleotide::NUCLEOTIDE_OTHER:
    {
      return "Other";
    }
    default:
      // Unrecognized/Unsupported Nucleotide type.
      return "";
  }
  // Shouldn't reach here.
  return "";
}

bool VcfUtils::AddVariantToFilter(const string& line, GroupingFilter* filter) {
  // Lines with spaces or tabs or not valid.
  if (line.find(" ") != string::npos || line.find("\t") != string::npos) {
    cout << "ERROR: Invalid variant: " << line << " has wrong format. " << endl;
    return false;
  }
  if (line.find(":") == string::npos) {
    // Format Two: Snp_ID.
    filter->snp_ids_.insert(line);
  } else {
    // Format One: CHROMOSOME:POSITION.
    int delimeter_pos = line.find(":");
    Chromosome chromosome;
    if (!ParseChromosome(line.substr(0, delimeter_pos), &chromosome)) {
      cout << "ERROR: Unable to parse Chromosome from line: " << line << endl;
      return false;
    }
    int64_t position;
    if (!StringUtils::Stoi(line.substr(delimeter_pos + 1), &position)) {
      cout << "ERROR: Unable to parse Position from line: " << line << endl;
      return false;
    }
    map<Chromosome, set<int64_t>>::iterator chr_itr =
        filter->variants_by_position_.find(chromosome);
    if (chr_itr == filter->variants_by_position_.end()) {
      // Don't yet have a Value yet for this chromosome. Create set and add it.
      set<int64_t> variants;
      variants.insert(position);
      filter->variants_by_position_.insert(make_pair(chromosome, variants));
    } else {
      // Already have variants for this chromosome. Add this one to the set.
      set<int64_t>& variants = chr_itr->second;
      variants.insert(position);
    }
  }
  return true;
}

bool VcfUtils::GetVariantsToKeep(
    const OutputOptions& output_options,
    GroupingFilter* filter, int* num_variants_to_keep) {
  if (output_options.segment_type_ == SegmentType::LIST) {
    const ListGroupingInfo& list_info = output_options.by_list_;
    // Sanity check input.
    if (list_info.current_file_index_ < 0 ||
        list_info.current_file_index_ >= list_info.list_files_.size()) {
      cout << "ERROR: Current file index (" << list_info.current_file_index_
           << ") is not compatible with number of Variant files provided ("
           << list_info.list_files_.size() << "). No output will be generated."
           << endl;
      return false;
    }

    // Read variants from file into filter->snp_ids_.
    const string& variants_filename =
        list_info.list_files_[list_info.current_file_index_];
    ifstream variants_file(variants_filename, ios::in);
    if (!variants_file.is_open()) {
      cout << "ERROR: Unable to open current variants file: "
           << variants_filename << ". No output will be generated." << endl;
      return false;
    }
    string line;
    while (getline(variants_file, line)) {
      (*num_variants_to_keep)++;
      if (!AddVariantToFilter(line, filter)) {
        cout << "ERROR: Unable to process line of Variant File '"
             << variants_filename << ": " << line << endl
             << "No output will be generated." << endl;
        return false;
      }
    }
    variants_file.close();
  } else if (output_options.segment_type_ == SegmentType::GENE) {
    const GeneGroupingInfo grouping_info = output_options.by_gene_;
    map<string, tuple<Chromosome, int64_t, int64_t>>::const_iterator gene_itr =
        grouping_info.name_to_position_.find(grouping_info.current_gene_);
    if (gene_itr == grouping_info.name_to_position_.end()) {
      cout << "ERROR: Unable to find gene: " << grouping_info.current_gene_
           << " in gene position map. Aborting." << endl;
      return false;
    }
    filter->gene_position_ = gene_itr->second;
    return true;
  } else if (output_options.segment_type_ == SegmentType::WINDOW) {
    const WindowGroupingInfo info = output_options.by_window_;
    *num_variants_to_keep = info.window_size_;
  } else {
    cout << "ERROR: Unrecognized segment_type: "
         << output_options.segment_type_ << ". Aborting." << endl;
    return false;
  }
  return true;
}

bool VcfUtils::GetPosition(const string& column, int64_t* position) {
  return StringUtils::Stoi(column, position);
}

Nucleotide VcfUtils::GetNucleotide(const char c) {
  if (c == 'A') return Nucleotide::A;
  if (c == 'C') return Nucleotide::C;
  if (c == 'G') return Nucleotide::G;
  if (c == 'T') return Nucleotide::T;
  if (c == 'U') return Nucleotide::U;
  if (c == 'R') return Nucleotide::R;
  if (c == 'Y') return Nucleotide::Y;
  if (c == 'W') return Nucleotide::W;
  if (c == 'S') return Nucleotide::S;
  if (c == 'M') return Nucleotide::M;
  if (c == 'K') return Nucleotide::K;
  if (c == 'N') return Nucleotide::N;
  if (c == '.') return Nucleotide::NUCLEOTIDE_INSERTION;
  // Unexpected nucelotide.
  cout << "ERROR: Unexpected Nucleotide: " << c << endl;
  return Nucleotide::NUCLEOTIDE_NONE;
}

bool VcfUtils::GetAlleleTypes(
    const string& ref, const string& alt,
    Nucleotide* major, Nucleotide* minor, Nucleotide* triallele) {
  if (ref.empty() || alt.empty()) {
    cout << "ERROR Empty alt and/or ref allele: " << ref << ", " << alt << endl;
    return false;
  }

  // Parse Major Allele.
  if (ref.length() != 1) {
    *major = NUCLEOTIDE_MULTI;
  } else {
    *major = GetNucleotide(ref.at(0));
  }

  // Parse Minor Allele(s).
  vector<string> minors;
  StringUtils::Split(alt, ",", &minors);
  // TODO(PHB): Should we handle more than 3?
  if (minors.size() >= 3) {
    cout << "ERROR: Too many minor alleles (can handle at most 2) for alt: "
         << alt << ". Only the first 2 are taken." << endl;
  } else if (minors.size() == 1) {
    *triallele = NUCLEOTIDE_NONE;
  }
  for (int i = 0; i < minors.size() && i < 2; ++i) {
    const string& allele_str = minors[i];
    Nucleotide* allele = i == 0 ? minor : triallele;
    if (allele_str.length() != 1) {
      if (*major != NUCLEOTIDE_MULTI) {
        *allele = NUCLEOTIDE_INSERTION;
      } else if (ref.length() == allele_str.length()) {
        *allele = NUCLEOTIDE_SUBSTITUTION;
      } else if (ref.length() > allele_str.length()) {
        *allele = NUCLEOTIDE_DELETION;
      } else {
        *allele = NUCLEOTIDE_INSERTION;
      }
    } else if (*major == NUCLEOTIDE_MULTI) {
      *allele = NUCLEOTIDE_DELETION;
    } else {
      *allele = GetNucleotide(allele_str.at(0));
    }
  }
  return true;
}

bool VcfUtils::IsPositionInRange(
    const int64_t position, const map<int64_t, int64_t>& ranges) {
  if (ranges.empty()) return false;
  map<int64_t, int64_t>::const_iterator lower_bound =
      ranges.lower_bound(position);
  if (lower_bound == ranges.end()) {
    // All Range Start positions are smaller than position; check the last one.
    lower_bound--;
    return position <= lower_bound->second;
  } else if (lower_bound != ranges.begin()) {
    // Current range's start position is bigger than position; check previous.
    lower_bound--;
    return position <= lower_bound->second;
  } else {
    // All range start positions are bigger than position.
    return false;
  }
}

bool VcfUtils::IsPositionFiltered(
    const Chromosome chromosome, const int64_t& position,
    const PositionFilter& position_filter) {
  if (!position_filter.chromosomes_to_filter_.empty()) {
    const bool in_set =
        position_filter.chromosomes_to_filter_.find(chromosome) !=
        position_filter.chromosomes_to_filter_.end();
    return (position_filter.filter_type_ == FilterType::EXCLUDE && in_set) ||
           (position_filter.filter_type_ == FilterType::INCLUDE && !in_set);
  }
  if (!position_filter.chromosome_position_pairs_to_filter_.empty()) {
    map<Chromosome, set<int64_t>>::const_iterator chr_itr =
        position_filter.chromosome_position_pairs_to_filter_.find(chromosome);
    if (chr_itr == position_filter.chromosome_position_pairs_to_filter_.end()) {
      return position_filter.filter_type_ == FilterType::INCLUDE;
    }
    const set<int64_t>& positions = chr_itr->second;
    const bool in_set = positions.find(position) != positions.end();
    return (position_filter.filter_type_ == FilterType::EXCLUDE && in_set) ||
           (position_filter.filter_type_ == FilterType::INCLUDE && !in_set);
  }
  if (position_filter.chromosome_range_pairs_to_filter_.empty()) {
    // No position filters set; no filtering should be done.
    return false;
  }
  map<Chromosome, map<int64_t, int64_t>>::const_iterator range_itr =
      position_filter.chromosome_range_pairs_to_filter_.find(chromosome);
  if (range_itr == position_filter.chromosome_range_pairs_to_filter_.end()) {
    // No filtering on this chromosome, return.
    return false;
  }
  const map<int64_t, int64_t>& ranges = range_itr->second;
  if (position_filter.filter_type_ == FilterType::INCLUDE) {
    return !IsPositionInRange(position, ranges);
  } else {
    return IsPositionInRange(position, ranges);
  }
}

bool VcfUtils::IsSnpFiltered(const string& snp_id,
                             const set<string>& snps_to_discard,
                             const set<string>& snps_to_keep) {
  if (!snps_to_discard.empty()) {
    return snps_to_discard.find(snp_id) != snps_to_discard.end();
  } else if (!snps_to_keep.empty()) {
    return snps_to_keep.find(snp_id) == snps_to_keep.end();
  }
  return false;
}

// TODO(PHB): Currently, only allow filtering on minor and trialleles (i.e.
// not allowed for ref allele). Consider expanding to allow filtering based
// on ref allele, as well as based on more general nucleotide settings (e.g.
// Purine or Weak or Keto or Indele, etc.). For some of these, the input
// has enough info, just need to modify the function below, e.g. if allele_filter
// specifies to filter Purine, should check if (alt/tri) is A, G, or R.
bool VcfUtils::IsAlleleFiltered(
    const Nucleotide& alt, const Nucleotide& triallele,
    const double& minor_freq, const double& tri_freq,
    const AlleleTypeFilter& allele_filter) {
  // Filter by AlleleType (INSERTION, DELETION, SUBSTITUION).
  if (!allele_filter.type_.empty()) {
    if (IsAlleleTypeFiltered(
            alt, allele_filter.filter_type_, allele_filter.type_) ||
        (triallele != NUCLEOTIDE_NONE &&
         IsAlleleTypeFiltered(
             triallele, allele_filter.filter_type_, allele_filter.type_))) {
      return true;
    }
  }

  // Filter by NumAlleleVariant (BIALLELIC, TRIALLELIC, QUADALLELIC).
  if (!allele_filter.number_.empty()) {
    NumAlleleVariant num = triallele == NUCLEOTIDE_NONE ?
        NumAlleleVariant::BIALLELIC : NumAlleleVariant::TRIALLELIC;
    const bool in_set =
        allele_filter.number_.find(num) != allele_filter.number_.end();
    if ((allele_filter.filter_type_ == FilterType::INCLUDE && !in_set) ||
        (allele_filter.filter_type_ == FilterType::EXCLUDE && in_set)) {
      return true;
    }
  }

  // Filter by Minor/Tri-Allele Frequency.
  if (allele_filter.filter_type_ == FilterType::INCLUDE) {
    if (minor_freq < allele_filter.minor_min_ || 
        (allele_filter.minor_max_ > 0.0 &&
         minor_freq > allele_filter.minor_max_)) {
      return true;
    }
    if (triallele != NUCLEOTIDE_NONE &&
        (tri_freq < allele_filter.tri_min_ || 
         (allele_filter.tri_max_ > 0.0 &&
          tri_freq > allele_filter.tri_max_))) {
      return true;
    }
  } else if (allele_filter.filter_type_ == FilterType::EXCLUDE) {
    if (minor_freq < allele_filter.minor_max_ || 
        (allele_filter.minor_min_ > 0.0 &&
         minor_freq > allele_filter.minor_min_)) {
      return true;
    }
    if (triallele != NUCLEOTIDE_NONE &&
        (tri_freq < allele_filter.tri_max_ || 
         (allele_filter.tri_min_ > 0.0 &&
          tri_freq > allele_filter.tri_min_))) {
      return true;
    }
  }
  return false;
}

AlleleType VcfUtils::GetAlleleType(const Nucleotide n) {
  if (n == Nucleotide::NUCLEOTIDE_INSERTION) return AlleleType::INSERTION;
  if (n == Nucleotide::NUCLEOTIDE_DELETION) return AlleleType::DELETION;
  return AlleleType::SUBSTITUION;
}

bool VcfUtils::IsAlleleTypeFiltered(
    const Nucleotide& alt, const FilterType& filter_type,
    const set<AlleleType>& type_filter) {
  if (type_filter.empty()) return false;
  AlleleType type = GetAlleleType(alt);
  const bool in_set = type_filter.find(type) != type_filter.end();
  return (filter_type == FilterType::INCLUDE && !in_set) ||
         (filter_type == FilterType::EXCLUDE && in_set);

}

bool VcfUtils::IsQualityFiltered(
    const string& quality_str, const double& min_quality_threshold) {
  double quality;
  if (!StringUtils::Stod(quality_str, &quality)) {
    cout << "ERROR: Unable to parse quality: " << quality_str
         << " as a double." << endl;
    if (min_quality_threshold > 0.0) return true;
  }
  return quality < min_quality_threshold;
}

bool VcfUtils::IsFilterFiltered(
    const string& column, const set<string>& filter_tags) {
  return (filter_tags.find(column) != filter_tags.end()) ||
         (column != "PASS" && filter_tags.find("ALL") != filter_tags.end());
}

bool VcfUtils::IsKeyValueFiltered(
    const string& key, const string& value,
    const vector<KeyValueFilter>& filters) {
  vector<string> parts;
  StringUtils::Split(value, ",", &parts);
  for (const KeyValueFilter& filter : filters) {
    if (key != filter.id_) continue;
    bool in_set = false;
    bool has_match = false;
    for (const string& part : parts) {
      if (filter.type_ == KeyType::KT_STRING) {
        in_set = filter.strings_to_filter_.find(part) == filter.strings_to_filter_.end();
      } else if (filter.type_ == KeyType::KT_INTEGER) {
        int v;
        if (!StringUtils::Stoi(part, &v)) {
          cout << "ERROR: Unable to parse Info part " << part
               << " from info column: " << value
               << " as an integer." << endl;
          return true;
        }
        in_set =
            (!filter.filter_zero_ || v != 0) &&
            (filter.int_value_ == 0 || v == filter.int_value_) &&
            ((filter.min_threshold_ == 0.0 || v >= filter.min_threshold_) &&
             (filter.max_threshold_ == 0.0 || v <= filter.max_threshold_)); 
      } else if (filter.type_ == KeyType::KT_FLOAT) {
        double v;
        if (!StringUtils::Stod(part, &v)) {
          cout << "ERROR: Unable to parse Info part " << part
               << " from info column: " << value
               << " as an double." << endl;
          return true;
        }
        in_set =
            (!filter.filter_zero_ || v != 0.0) &&
            (filter.float_value_ == 0.0 || v == filter.float_value_) &&
            ((filter.min_threshold_ == 0.0 || v >= filter.min_threshold_) &&
             (filter.max_threshold_ == 0.0 || v <= filter.max_threshold_)); 
      }
      if (filter.filter_type_ == FilterType::INCLUDE && in_set) {
        // Passed filter. Break.
        has_match = true;
        break;
      }
      if (filter.filter_type_ == FilterType::EXCLUDE && in_set) {
        return true;
      }
    }
    if (filter.filter_type_ == FilterType::INCLUDE && !has_match) {
      // None of the parts passed the filter (otherwise we would've broken
      // out of for loop above). Thus, filter imposed.
      return true;
    }
  }
  return false;
}

bool VcfUtils::IsInfoFiltered(
    const string& info, const vector<KeyValueFilter>& info_column_filters,
    double* minor_freq, double* tri_freq) {
  vector<string> fields;
  StringUtils::Split(info, ";", &fields);
  for (const string& field : fields) {
    vector<string> parts;
    StringUtils::Split(field, "=", &parts);
    if (parts.size() != 2) {
      cout << "ERROR: Unexpected INFO column field part: " << parts.size()
           << " in column field: " << field
           << " in column: " << info << endl;
      return true;
    }
    const string& key = parts[0];
    const string& value = parts[1];
    if (IsKeyValueFiltered(key, value, info_column_filters)) return true;
    // Record allele frequencies.
    if (key == "AF") {
      vector<string> frequencies;
      StringUtils::Split(value, ",", &frequencies);
      StringUtils::Stod(frequencies[0], minor_freq);
      if (frequencies.size() > 1) {
        StringUtils::Stod(frequencies[1], tri_freq);
      }
    }
  }
  return false;
}

SnpType VcfUtils::GetSnpType(const char c) {
  if (c == '0') return SnpType::SNP_TYPE_MAJOR;
  if (c == '1') return SnpType::SNP_TYPE_MINOR;
  if (c == '2') return SnpType::SNP_TYPE_TRI;
  return SnpType::SNP_TYPE_UNKNOWN;
}

// TODO(PHB): Implement usage of data_filters and populating filtering_counts.
// TODO(PHB): Are there ever NA values, or missing values, for genotype?
// If so, how to handle these gracefully (right now, I abort if there are
// any invalid/missing genotype values).
bool VcfUtils::ProcessDataColumn(
    const string& column, const int gt_index,
    const vector<KeyValueFilter>& data_filters,
    const SnpOutputType& output_type,
    FilteringCounts* filtering_counts,
    SnpType* snp_type) {
  // Get Genotype portion of the Column.
  vector<string> column_parts;
  StringUtils::Split(column, ":", &column_parts);
  if (column_parts.size() <= gt_index) {
    cout << "ERROR: Found " << column_parts.size() << " items in column "
         << column << ", but gt_index is: " << gt_index << endl;
    return false;
  }

  // Break Genotype into info about each allele.
  const string& gt = column_parts[gt_index];
  const int gt_length = gt.length();
  // Only support Haploids and Diploids; abort if we don't have this.
  if (gt_length != 1 && gt_length != 3) {
    cout << "ERROR: Unexpected Genotype: " << gt << endl;
    return false;
  }

  // Handle Haploid
  if (gt_length == 1) {
    SnpType type = GetSnpType(gt.at(0));
    if (type == SnpType::SNP_TYPE_UNKNOWN) {
      cout << "ERROR: Unexpected Haploid Genotype: " << gt << endl;
      return false;
    }
    if (output_type == SnpOutputType::COUNT_ONLY) {
      *snp_type = (type == SnpType::SNP_TYPE_MAJOR) ?
          SnpType::SNP_TYPE_NONE : SnpType::SNP_TYPE_ONE;
      return true;
    } else if (output_type == SnpOutputType::PHASED ||
               output_type == SnpOutputType::UNPHASED) {
      *snp_type = type;
      return true;
    } else {
      cout << "ERROR: Unrecognized output_type: " << output_type << endl;
      return false;
    }
  }

  // Handle Diploid.
  if (gt.at(1) != '/' && gt.at(1) != '|') {
    cout << "ERROR: Unrecognized phase separator in genotype: " << gt << endl;
    return false;
  }
  SnpType first_allele = GetSnpType(gt.at(0));
  SnpType second_allele = GetSnpType(gt.at(2));
  if (first_allele == SnpType::SNP_TYPE_UNKNOWN ||
      second_allele == SnpType::SNP_TYPE_UNKNOWN) {
    cout << "ERROR: Unexpected Haploid Genotype: " << gt << endl;
    return false;
  }
  if (output_type == SnpOutputType::COUNT_ONLY) {
    int num_minor =
        (first_allele == SnpType::SNP_TYPE_MAJOR ? 0 : 1) +
        (second_allele == SnpType::SNP_TYPE_MAJOR ? 0 : 1);
    *snp_type = num_minor == 0 ?
      SnpType::SNP_TYPE_NONE :
      (num_minor == 1 ? SnpType::SNP_TYPE_ONE : SnpType::SNP_TYPE_TWO);
    return true;
  } else if (output_type == SnpOutputType::UNPHASED) {
    if (first_allele == SnpType::SNP_TYPE_MAJOR) {
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        *snp_type = SnpType::SNP_TYPE_MAJOR_MAJOR;
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        *snp_type = SnpType::SNP_TYPE_MAJOR_U_MINOR;
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        *snp_type = SnpType::SNP_TYPE_MAJOR_U_TRI;
      }
    } else if (first_allele == SnpType::SNP_TYPE_MINOR) {
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        *snp_type = SnpType::SNP_TYPE_MAJOR_U_MINOR;
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        *snp_type = SnpType::SNP_TYPE_MINOR_MINOR;
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        *snp_type = SnpType::SNP_TYPE_MINOR_U_TRI;
      }
    } else {  // first_allele == SnpType::SNP_TYPE_TRI
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        *snp_type = SnpType::SNP_TYPE_MAJOR_U_TRI;
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        *snp_type = SnpType::SNP_TYPE_MINOR_U_TRI;
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        *snp_type = SnpType::SNP_TYPE_TRI_TRI;
      }
    }
  } else if (output_type == SnpOutputType::PHASED) {
    if (first_allele == SnpType::SNP_TYPE_MAJOR) {
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        *snp_type = SnpType::SNP_TYPE_MAJOR_MAJOR;
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_MAJOR_P_MINOR;
        } else {
          *snp_type = SnpType::SNP_TYPE_MAJOR_U_MINOR;
        }
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_MAJOR_P_TRI;
        } else {
          *snp_type = SnpType::SNP_TYPE_MAJOR_U_TRI;
        }
      }
    } else if (first_allele == SnpType::SNP_TYPE_MINOR) {
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_MINOR_P_MAJOR;
        } else {
          *snp_type = SnpType::SNP_TYPE_MAJOR_U_MINOR;
        }
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        *snp_type = SnpType::SNP_TYPE_MINOR_MINOR;
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_MINOR_P_TRI;
        } else {
          *snp_type = SnpType::SNP_TYPE_MINOR_U_TRI;
        }
      }
    } else {  // first_allele == SnpType::SNP_TYPE_TRI
      if (second_allele == SnpType::SNP_TYPE_MAJOR) {
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_TRI_P_MAJOR;
        } else {
          *snp_type = SnpType::SNP_TYPE_MAJOR_U_TRI;
        }
      } else if (second_allele == SnpType::SNP_TYPE_MINOR) {
        if (gt.at(1) == '/') {
          *snp_type = SnpType::SNP_TYPE_TRI_P_MINOR;
        } else {
          *snp_type = SnpType::SNP_TYPE_MINOR_U_TRI;
        }
      } else {  // second_allele == SnpType::SNP_TYPE_TRI
        *snp_type = SnpType::SNP_TYPE_TRI_TRI;
      }
    }
  } else {
    cout << "ERROR: Unrecognized output_type: " << output_type << endl;
    return false;
  }
  return true;
}

void VcfUtils::SetSamplesToFilterFromFile(
    const string& filename, set<string>* samples_to_filter) {
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << file << ". Aborting." << endl;
    return;
  }
  string line;
  while (getline(file, line)) {
    samples_to_filter->insert(line);
  }
  file.close();
}

void VcfUtils::SetVariantsToFilterFromFile(
    const string& filename, set<string>* snps_to_filter) {
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << file << ". Aborting." << endl;
    return;
  }
  string line;
  while (getline(file, line)) {
    snps_to_filter->insert(line);
  }
  file.close();
}

void VcfUtils::SetVariantsToFilterFromFile(const string& filename,
                                           PositionFilter* filter_positions) {
  // Open File.
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << file << ". Aborting." << endl;
    return;
  }

  // Read Lines of file.
  string line;
  int line_index = 1;
  while (getline(file, line)) {
    vector<string> line_parts;
    StringUtils::Split(line, "\t", &line_parts);
    // Check Line has appropriate format.
    if (line_parts.size() != 2 && line_parts.size() != 3) {
      cout << "ERROR in parsing file " << filename
           << ": Unexpected line " << line_index << ":\n" << line
           << "\nDoes not have expected format:\nCHROMOSOME POSITION" << endl;
      return;
    }
    // Parse Chromosome.
    Chromosome chr;
    if (!ParseChromosome(line_parts[0], &chr)) {
      cout << "ERROR in parsing file " << filename
           << ": Unable to parse line " << line_index
           << " which has unparsable Chromosome:\n"
           << line << endl;
      return;
    }
    // Parse Position(s).
    if ((line_parts.size() == 2)) {
      // Line has format: CHROMOSOME POSITION. Parse Position.
      int64_t position;
      if (!StringUtils::Stoi(line_parts[1], &position)) {
        cout << "ERROR in parsing file " << filename
             << ": Unable to parse line " << line_index
             << " which has unparsable position:\n"
             << line << endl;
        return;
      }
      // Add Line to filter_positions->chromosome_position_pairs_to_filter_.
      map<Chromosome, set<int64_t>>::iterator itr =
          filter_positions->chromosome_position_pairs_to_filter_.find(chr);
      if (itr == filter_positions->chromosome_position_pairs_to_filter_.end()) {
        // Haven't seen this chromosome yet. Create map and insert.
        set<int64_t> to_add;
        to_add.insert(position);
        filter_positions->chromosome_position_pairs_to_filter_.insert(make_pair(
            chr, to_add));
      } else {
        // Already created a map for this chromosome. Add to it.
        itr->second.insert(position);
      }
    } else {
      // Line has format: CHROMOSOME RANGE_START RANGE_END.
      // Parse range start and end.
      int64_t range_start, range_end;
      if (!StringUtils::Stoi(line_parts[1], &range_start)) {
        cout << "ERROR in parsing file " << filename
             << ": Unable to parse line " << line_index
             << " which has unparsable range start:\n"
             << line << endl;
        return;
      }
      if (!StringUtils::Stoi(line_parts[2], &range_end)) {
        cout << "ERROR in parsing file " << filename
             << ": Unable to parse line " << line_index
             << " which has unparsable range end:\n"
             << line << endl;
        return;
      }
      // Add Line to filter_positions->chromosome_range_pairs_to_filter_.
      map<Chromosome, map<int64_t, int64_t>>::iterator itr =
          filter_positions->chromosome_range_pairs_to_filter_.find(chr);
      if (itr == filter_positions->chromosome_range_pairs_to_filter_.end()) {
        // Haven't seen this chromosome yet. Create map and insert.
        map<int64_t, int64_t> to_add;
        to_add.insert(make_pair(range_start, range_end));
        filter_positions->chromosome_range_pairs_to_filter_.insert(make_pair(
            chr, to_add));
      } else {
        // Already created a map for this chromosome. Add to it.
        itr->second.insert(make_pair(range_start, range_end));
      }
    }
    line_index++;
  }
  file.close();
}

void VcfUtils::GetGenePositionMapFromFile(
    const string& filename,
    map<string, tuple<Chromosome, int64_t, int64_t>>* gene_name_to_position) {
  // Open File.
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << file << ". Aborting." << endl;
    return;
  }

  // Read Lines of file.
  string line;
  int line_index = 1;
  while (getline(file, line)) {
    vector<string> line_parts;
    StringUtils::Split(line, "\t", &line_parts);
    // Check Line has appropriate format.
    if (line_parts.size() != 4 && line_parts.size() != 5) {
      cout << "ERROR in parsing file " << filename
           << ": Unexpected line " << line_index << ":\n" << line
           << "\nDoes not have expected format:\nCHROMOSOME POSITION" << endl;
      return;
    }
    // Parse Gene_ID.
    const string& gene_id = line_parts[0];
    // Parse Chromosome.
    Chromosome chr;
    if (!ParseChromosome(line_parts[1], &chr)) {
      cout << "ERROR in parsing file " << filename
           << ": Unable to parse line " << line_index
           << " which has unparsable Chromosome:\n"
           << line << endl;
      return;
    }
    // Parse range start and end.
    int64_t range_start, range_end;
    if (!StringUtils::Stoi(line_parts[2], &range_start)) {
      cout << "ERROR in parsing file " << filename
           << ": Unable to parse line " << line_index
           << " which has unparsable range start:\n"
           << line << endl;
      return;
    }
    if (!StringUtils::Stoi(line_parts[3], &range_end)) {
      cout << "ERROR in parsing file " << filename
           << ": Unable to parse line " << line_index
           << " which has unparsable range end:\n"
           << line << endl;
      return;
    }
    // Parse Strand, if present.
    if (line_parts.size() == 5) {
      if (line_parts[4] == "-") {
        cout << "ERROR in parsing file " << filename
             << ": Line " << line_index
             << " indicates Negative Strand, which is currently unsupported:\n"
             << line << endl;
        return;
      } else if (line_parts[4] != "+") {
        cout << "ERROR in parsing file " << filename
             << ": Line " << line_index
             << " has unrecognized value in Strand Column:\n"
             << line << endl;
        return;
      }
    }
    // Add Line to gene_name_to_position.
    tuple<Chromosome, int64_t, int64_t> to_add =
        make_tuple(chr, range_start, range_end);
    gene_name_to_position->insert(make_pair(gene_id, to_add));

    line_index++;
  }
  file.close();
}

bool VcfUtils::IsMetadataLine(const string& line) {
  // Metadata lines are prefixed with a repeated Comment Char '#'.
  return StringUtils::HasPrefixString(line, "##");
}

bool VcfUtils::ProcessHeaderLine(
    const string& line, const FilteringOptions& filtering_options,
    vector<string>* header_row) {
  // Check that this line is the Header Line.
  if (!StringUtils::HasPrefixString(line, "#CHROM")) return false;
  vector<string> temp_header;
  if (!StringUtils::Split(line, "\t", &temp_header)) return false;
  // Check header_row has some data columns.
  if (temp_header.size() < 10) return false;
  string& first_col = temp_header[0];
  // Get rid of Comment character '#' at start of Header Line.
  first_col = first_col.substr(1);
  // Get rid of Samples in the filter list.
  const set<string>& samples = filtering_options.samples_to_discard_.empty() ?
      filtering_options.samples_to_keep_ :
      filtering_options.samples_to_discard_;
  const bool to_keep = filtering_options.samples_to_discard_.empty();
  for (int i = 0; i < temp_header.size(); ++i) {
    if (i < 9) {
      header_row->push_back(temp_header[i]);
    } else if (samples.empty()) {
      header_row->push_back(temp_header[i]);
    } else if (samples.find(temp_header[i]) == samples.end()) {
      if (to_keep) {
        // Sample is not in the set of Samples to keep, skip it.
        continue;
      } else {
        header_row->push_back(temp_header[i]);
      }
    } else if (to_keep) {
      header_row->push_back(temp_header[i]);
    }
    // Remaining case is that Sample is in samples and to_keep is false, so this
    // Sample is in to_discard, so it should be skipped; i.e. do nothing here.
  }
  return true;
}

void VcfUtils::GetHeaderRow(
    const string& filename,
    const FilteringOptions& filtering_options,
    vector<string>* header_row) {
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << file << ". Aborting." << endl;
    return;
  }
  string line;
  while (getline(file, line)) {
    if (IsMetadataLine(line)) {
      // TODO(PHB): Determine if we need to store any info from metadata, and if
      // so, do it here.
      continue;
    } else if (!ProcessHeaderLine(line, filtering_options, header_row)) {
      cout << "ERROR: Unrecognized line in vcf file is neither Metadata nor "
           << "Header Row. Aborting." << endl;
    }
    break;
  }
  file.close();
}

void VcfUtils::ProcessVcfFile(
    const string& filename,
    const FilteringOptions& filtering_options,
    const OutputOptions& output_options,
    vector<string>* header_row,
    int64_t* next_char_to_read_pos,
    int* num_lines_processed,
    FilteringCounts* filtering_counts,
    vector<vector<SnpType>>* snp_values_per_sample,
    vector<SnpInfo>* snp_info) {
  // Sanity check input.
  if (header_row == nullptr || next_char_to_read_pos == nullptr ||
      filtering_counts == nullptr || snp_values_per_sample == nullptr ||
      snp_info == nullptr) {
    cout << "ERROR: Null parameters passed to ProcessVcfFile. Aborting." << endl;
    return;
  }
  
  // Get header_row, if not already populated.
  if (header_row->empty()) {
    GetHeaderRow(filename, filtering_options, header_row);
  }
  if (header_row->empty()) return;
  const int num_expected_cols = header_row->size();

  // Open file at end (ios::ate) so we can first grab its size.
  ifstream file(filename, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << filename << ". Aborting." << endl;
    return;
  }

  // If output_options indicate Variant Segments are specified by file
  // (SegmentType::LIST), store the list of variants in a local variable.
  GroupingFilter grouping_filter;
  int num_variants_to_keep;
  if (!GetVariantsToKeep(
          output_options, &grouping_filter, &num_variants_to_keep)) {
    return;
  }

  // Set total file_size, and from it, block_size.
  int64_t file_size = file.tellg();
  if (next_char_to_read_pos != nullptr && *next_char_to_read_pos == file_size) {
    return;
  }
  const int64_t block_size = min(MAX_BUFFER_SIZE, file_size);
  char* memblock = new char[block_size];

  // Return to beginning of file, or *next_char_to_read_pos if set.
  int64_t start_pos =
      (next_char_to_read_pos == nullptr || *next_char_to_read_pos < 0) ?
      0 : *next_char_to_read_pos;
  file.seekg(start_pos, ios::beg);
  int64_t iteration = 1;
  int64_t next_block_start_pos = file.tellg();
  int64_t prev_start_pos = -1;
  int num_variants_taken = 0;
  string final_line;
  // Garbage local variable to simplify code below, used if caller of
  // this function has set next_char_to_read_pos to null.
  int64_t next_start_pos;
  if (next_char_to_read_pos == nullptr) next_char_to_read_pos = &next_start_pos;
  while (true) {
    // Sanity check we've made progress since the last block was read.
    int64_t current_file_pos = file.tellg();
    if (prev_start_pos == next_block_start_pos && prev_start_pos >= 0) {
      // We enter here when either:
      //  a) The last line from the previous block was so long that there
      //     wasn't a single line-break in the entire block; OR
      //  b) We reached the last line of the file.
      // In either case, we store the last line, and abort.
      if (!final_line.empty()) {
        snp_values_per_sample->push_back(vector<SnpType>());
        ProcessLineResponse response_type = ProcessLine(
            final_line, num_expected_cols, *header_row, grouping_filter,
            filtering_options,
            output_options.output_type_, filtering_counts,
            &(snp_values_per_sample->back()), snp_info);
        if (response_type == ProcessLineResponse::FAILURE) {
          cout << "ERROR processing line number " << *num_lines_processed
               << ":\n" << final_line << endl;
          return;
        }
        (*num_lines_processed)++;
      }
      if (current_file_pos + final_line.length() != file_size) {
        // Case (a): This line is too long (fills up the entire buffer before
        // reaching the NewLine character.
        cout << "ERROR: Line starting at position " << current_file_pos
             << " is too long. Aborting." << endl;
      }
      *next_char_to_read_pos = file_size;
      break;
    } else if (current_file_pos != next_block_start_pos) {
      file.seekg(next_block_start_pos);
    }

    prev_start_pos = file.tellg();
    const int64_t remaining_bytes = file_size - prev_start_pos;
    // Check if we're done reading file.
    if (remaining_bytes <= 0) {
      *next_char_to_read_pos = file_size;
      break;
    }

    // Resize memblock, if this is the final block and there are fewer
    // than block_size bytes left.
    int64_t current_block_size = min(remaining_bytes, block_size);
    if (current_block_size < block_size) {
      delete[] memblock;
      memblock = new char[current_block_size];
    }

    // Read block.
    file.read(memblock, current_block_size);

    // Process block.
    stringstream ss;
    ss.rdbuf()->pubsetbuf(memblock, current_block_size);

    // Determine if this block ends right at the end of line, as
    // determined by the New Line character '10'.
    ss.seekg(-1, ios::end);
    int end_char = ss.peek();
    bool last_char_is_newline = end_char == 10;
    ss.seekg(0, ios::beg);

    // Process Lines.
    string str;
    int64_t current_stream_pos;
    bool at_least_one_line = false;
    bool is_first_line_of_block = true;
    bool should_break = false;
    while(getline(ss, str)) {
      if (is_first_line_of_block) {
        is_first_line_of_block = false;
        str = final_line + str;
      }
      // Check if we're processing the last line of the block.
      current_stream_pos = ss.tellg();
      if (current_stream_pos == -1) {
        // Final linee of the block. Don't process it now, but store it
        // in 'final_line', so it can processed with the next block.
        final_line = str;
      } else {
        // Not the last line of this block; process it.
        at_least_one_line = true;
        snp_values_per_sample->push_back(vector<SnpType>());
        ProcessLineResponse response_type = ProcessLine(
            str, num_expected_cols, *header_row, grouping_filter,
            filtering_options,
            output_options.output_type_, filtering_counts,
            &(snp_values_per_sample->back()), snp_info);
        if (response_type == ProcessLineResponse::FAILURE) {
          cout << "ERROR processing line number " << *num_lines_processed
               << ":\n" << str << endl;
          return;
        } else if (response_type == ProcessLineResponse::SNP_IN_SEGMENT) {
          ++num_variants_taken;
        }
        (*num_lines_processed)++;
        should_break =
            IsSegmentDone(response_type, num_variants_to_keep,
                          num_variants_taken, output_options);
        if (should_break) {
          // We've read the last relevant row for this segment; mark
          // the position of the next character in the next segment by
          // setting next_char_to_read_pos appropriately.
          if (response_type == ProcessLineResponse::SNP_NOT_IN_SEGMENT) {
            // This line was not a part of this Segment. Set
            // next_char_to_read_pos to be the index representing the
            // first character on this line.
            *next_char_to_read_pos =
                prev_start_pos + current_stream_pos - str.length();
          } else {
            // This line was part of this Segment. Set next_char_to_read_pos to
            // be the index representing the first character on the next line.
            *next_char_to_read_pos = prev_start_pos + current_stream_pos;
          }
          break;
        }
      }
    }

    // Segment ended already, no need to keep reading file.
    if (should_break) break;

    // If last character in the block was NL character, then final_line
    // will not get reset from its value from the previous block. Reset
    // it now.
    if (last_char_is_newline) {
      final_line = "";
    }

    // Set next block start position.
    if (at_least_one_line) {
      next_block_start_pos = file.tellg();
    }
    iteration++;
  }

  file.close();
  delete[] memblock;
}

}  // namespace file_reader_utils
