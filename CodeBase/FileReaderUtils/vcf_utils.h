// Author: paulbunn@email.unc.edu (Paul Bunn)
// Last Updated: March 2015
//
// Description: Class with helper functions for reading/filtering vcf files.

#include "StringUtils/string_utils.h"

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifndef VCF_UTILS_H
#define VCF_UTILS_H

using namespace string_utils;
using namespace std;

namespace file_reader_utils {

// Possible types (String, Integer, Float, etc.) seen in vcf file for
// the keys in the INFO, FILTER, and FORMAT columns.
enum KeyType {
  KT_STRING,
  KT_INTEGER,
  KT_FLOAT,
};

enum SegmentType {
  WINDOW,
  GENE,
  LIST,
};

// Contents of the final output table will contain data of this type, which
// represents the SNP information for that (SNP, Sample) entry.
// Quad-allelic variants, and triploid and higher, are not supported.
enum SnpType {
  SNP_TYPE_UNKNOWN,      // Missing
  SNP_TYPE_NA,            // Same as UNKNOWN
  // Counts only.
  SNP_TYPE_NONE,          // 0 minor alleles (both alleles dominant)
  SNP_TYPE_ONE,           // 1 minor allele
  SNP_TYPE_TWO,           // 2 minor alleles (i.e. both alleles minor)
  // Haploid.
  SNP_TYPE_MAJOR,         // Genotype: '0'
  SNP_TYPE_MINOR,         // Genotype: '1'
  SNP_TYPE_TRI,           // Genotype: '2'
  // Diploid; phased.
  SNP_TYPE_MAJOR_MAJOR,   // Genotype: '0/0'
  SNP_TYPE_MAJOR_P_MINOR, // Genotype: '0/1'
  SNP_TYPE_MINOR_P_MAJOR, // Genotype: '1/0'
  SNP_TYPE_MINOR_MINOR,   // Genotype: '1/1'
  SNP_TYPE_MAJOR_P_TRI,   // Genotype: '0/2'
  SNP_TYPE_TRI_P_MAJOR,   // Genotype: '2/0'
  SNP_TYPE_TRI_TRI,       // Genotype: '2/2'
  SNP_TYPE_MINOR_P_TRI,   // Genotype: '1/2'
  SNP_TYPE_TRI_P_MINOR,   // Genotype: '2/1'
  // Diploid; unphased (also includes SNP_TYPE_MAJOR_MAJOR, SNP_TYPE_MINOR_MINOR,
  // and SNP_TYPE_TRI_TRI above).
  SNP_TYPE_MAJOR_U_MINOR, // Genotype: '0|1'
  SNP_TYPE_MAJOR_U_TRI,   // Genotype: '0|2'
  SNP_TYPE_MINOR_U_TRI,   // Genotype: '1|2'
};

enum SnpOutputType {
  COUNT_ONLY,  // Use only SnpTypes: NA, NONE, ONE, and TWO
  UNPHASED,    // Use only SnpTypes: NA, MAJOR_MAJOR, MINOR_MINOR, TRI_TRI,
               //                    MAJOR_U_MINOR, MAJOR_U_TRI, and MINOR_U_TRI
  PHASED,      // Use all SnpTypes, EXCEPT: NONE, ONE, and TWO
};

// Number of Possible Variants.
enum NumAlleleVariant {
  BIALLELIC,
  TRIALLELIC,
  QUADALLELIC,
};

// Type of Variant.
enum AlleleType {
  SUBSTITUION,
  INSERTION,
  DELETION,
};

enum FilterType {
  INCLUDE,
  EXCLUDE,
  // Only available for filters applied to a single data cell
  // (e.g. acutal_data_filters_), not for filters applied to an
  // entire row (e.g. position_filters, info_filters, etc.)
  MARK_NA,
};

enum Chromosome {
  CHROMOSOME_UNKNOWN,
  CHROMOSOME_ONE,
  CHROMOSOME_TWO,
  CHROMOSOME_THREE,
  CHROMOSOME_FOUR,
  CHROMOSOME_FIVE,
  CHROMOSOME_SIX,
  CHROMOSOME_SEVEN,
  CHROMOSOME_EIGHT,
  CHROMOSOME_NINE,
  CHROMOSOME_TEN,
  CHROMOSOME_ELEVEN,
  CHROMOSOME_TWELVE,
  CHROMOSOME_THIRTEEN,
  CHROMOSOME_FOURTEEN,
  CHROMOSOME_FIFTEEN,
  CHROMOSOME_SIXTEEN,
  CHROMOSOME_SEVENTEEN,
  CHROMOSOME_EIGHTEEN,
  CHROMOSOME_NINETEEN,
  CHROMOSOME_TWENTY,
  CHROMOSOME_TWENTYONE,
  CHROMOSOME_TWENTYTWO,
  CHROMOSOME_X,
  CHROMOSOME_Y,
};

enum ReferenceGenome {
  HG_UNKNONWN,
  HG18,
  HG19,
};

enum Nucleotide {
  NUCLEOTIDE_UNKNOWN,
  NUCLEOTIDE_NA,    // Same as unknown
  A,     // Adenine
  C,     // Cytosine
  G,     // Guanine
  T,     // Thymine
  U,     // Uracil
  R,     // Purine (A or G)
  Y,     // Pyrimidine (C or T)
  W,     // Weak (A or T)
  S,     // Strong (G or C)
  M,     // Amino (A or C)
  K,     // Keto (G or T)
  N,     // Any nucleotide
  NUCLEOTIDE_INSERTION,       // Not a nucleotide (indele of type insertion)
  NUCLEOTIDE_DELETION,        // Not a nucleotide (indele of type deletion)
  NUCLEOTIDE_SUBSTITUTION,    // Used for when REF and ALT have same length > 1
  NUCLEOTIDE_MULTI,           // A sequence (length > 1) of nucelotides
  NUCLEOTIDE_NONE,            // No Nucleotide present (used for e.g. triallele)
  NUCLEOTIDE_OTHER,           // Not a nucleotide, for some other reason
};

enum ProcessLineResponse {
  FAILURE,
  FILTERED,
  SNP_IN_SEGMENT,
  SNP_NOT_IN_SEGMENT,
};

// NOTE: SnpInfo below uses separate fields for Chromosome and Position instead
// of this structure, because this structure wasn't available when those fields
// of SnpInfo were introduced, and updating those fields to instead use Position
// requires some work...
struct SnpPosition {
  Chromosome chr_;
  int64_t pos_;

  SnpPosition() {
    chr_ = Chromosome::CHROMOSOME_UNKNOWN;
    pos_ = -1;
  }

  SnpPosition(Chromosome chr, const int64_t& pos) {
    chr_ = chr;
    pos_ = pos;
  }

  bool operator <(const SnpPosition& x) const {
    return std::tie(chr_, pos_) < std::tie(x.chr_, x.pos_);
  }
};

// A structure to hold information about a SNP.
struct SnpInfo {
  // SNP Name/Id (if it has one).
  string id_;

  // SNP Position.
  Chromosome chromosome_;
  int64_t position_;

  // In case this SNP corresponds to a row in some data file, this field
  // stores the row index.
  int row_index_;

  // Major/Minor Allele for this SNP.
  Nucleotide major_allele_;
  Nucleotide minor_allele_;
  Nucleotide triallele_;
  double minor_allele_freq_;
  double triallele_freq_;

  // In case we don't know what the major/minor allele is, so a data file
  // refers to the alleles via generic numbers. This fields give a mapping
  // from the generic numbers to the Nucleotide. Note that the field type is
  // 'string' instead of 'Nucleotide' to allow more flexibility (e.g. Plink
  // allows alleles to be any non-zero character).
  string allele_one_;
  string allele_two_;
  string allele_three_;
  string allele_four_;

  SnpInfo() {
    id_ = "";
    chromosome_ = Chromosome::CHROMOSOME_UNKNOWN;
    position_ = -1;
    row_index_ = -1;
    major_allele_ = Nucleotide::NUCLEOTIDE_UNKNOWN;
    minor_allele_ = Nucleotide::NUCLEOTIDE_UNKNOWN;
    triallele_ = Nucleotide::NUCLEOTIDE_UNKNOWN;
    minor_allele_freq_ = 0.0;
    triallele_freq_ = 0.0;
    allele_one_ = "";
    allele_two_ = "";
    allele_three_ = "";
    allele_four_ = "";
  }

  // So that we can put SnpInfo objects in a set or map, we need a comparator.
  bool operator <(const SnpInfo& x) const {
    return std::tie(chromosome_, position_, id_, major_allele_, minor_allele_,
                    triallele_, minor_allele_freq_, triallele_freq_,
                    allele_one_, allele_two_, allele_three_, allele_four_) <
           std::tie(x.chromosome_, x.position_, x.id_, x.major_allele_, x.minor_allele_,
                    x.triallele_, x.minor_allele_freq_, x.triallele_freq_,
                    x.allele_one_, x.allele_two_, x.allele_three_, x.allele_four_);
  }

  // Print out all fields.
  void Print(string* output) const {
    if (output == nullptr) return;
    *output += "id_: " + id_ + "\n";
    *output += "chromosome_: " + Itoa(chromosome_) + "\n";
    *output += "position_: " + Itoa(position_) + "\n";
    *output += "row_index_: " + Itoa(row_index_) + "\n";
    *output += "major_allele_: " + Itoa(major_allele_) + "\n";
    *output += "minor_allele_: " + Itoa(minor_allele_) + "\n";
    *output += "triallele_: " + Itoa(triallele_) + "\n";
    *output += "minor_allele_freq_: " + Itoa(minor_allele_freq_) + "\n";
    *output += "triallele_freq_: " + Itoa(triallele_freq_) + "\n";
    *output += "allele_one_: " + allele_one_ + "\n";
    *output += "allele_two_: " + allele_two_ + "\n";
    *output += "allele_three_: " + allele_three_ + "\n";
    *output += "allele_four_: " + allele_four_ + "\n";
  }
};

// A structure to hold filtering options based on allele-type.
struct AlleleTypeFilter {
  // Whether to include or exclude these alleles.
  FilterType filter_type_;
  // The categories of alleles to filter on. Determination of whether a given
  // Variant lies in the filter will be based on:
  //   - 'OR' within the type_ and number_ sets.
  //   - 'AND' between the two sets (unless empty, in which case it is ignored).
  set<AlleleType> type_;
  set<NumAlleleVariant> number_;

  // Filter based on minor allele frequency.
  double minor_min_;
  double minor_max_;
  // Filter based on tri-allele frequency.
  double tri_min_;
  double tri_max_;
};

// A structure to hold filtering options based on position.
struct PositionFilter {
  // Whether to include or exclude variants based on the below position criteria.
  FilterType filter_type_;

  // Filter by chromosome.
  set<Chromosome> chromosomes_to_filter_;

  // Filter by (chromosome, position).
  map<Chromosome, set<int64_t>> chromosome_position_pairs_to_filter_;
  map<Chromosome, map<int64_t, int64_t>> chromosome_range_pairs_to_filter_;
};

// A structure to hold filtering options based on vcf Key, Values
// (e.g. for INFO column).
struct KeyValueFilter {
  // The (typically 2 character) key/id for this info item.
  string id_;
  // The KeyType associated to this info item.
  KeyType type_;
  // Filter type.
  FilterType filter_type_;

  // Possible filters for KT_STRING types:
  set<string> strings_to_filter_;

  // Possible filters for KT_INTEGER and KT_FLOAT types. Can either filter by
  // exact value or via a (1- or 2-sided) range.
  // Discard any values that lie below this value.
  float min_threshold_;
  // Discard any values that lie above this value.
  float max_threshold_;
  // Filter based on exact value.
  // NOTE: Don't use these if the filter value is zero! Must use filter_zero_.
  int int_value_;
  float float_value_;
  bool filter_zero_;
};

// A structure to hold overall statistics (counts) for each filter.
struct FilteringCounts {
  int num_snp_id_column_exclusions_;
  int num_sample_id_exclusions_;
  int num_filter_column_exclusions_;
  int num_info_column_exclusions_;
  int num_na_values_found_;
  int num_na_values_created_;
  int num_low_quality_exclusions_;
  int num_allele_type_exclusions_;
  int num_position_exclusions_;
};

// A structure to hold all possible filtering options.
struct FilteringOptions {
  // Filter based on SNP Id (ID column of vcf file). Only one of the two sets
  // below should be populated (i.e. either keep only snps in the first set,
  // or keep all snps EXCEPT those in the second set).
  set<string> snps_to_keep_;
  set<string> snps_to_discard_;

  // Filter samples based on sample_id. Only one of the two sets
  // below should be populated (i.e. either keep only samples in the first set,
  // or keep all samples EXCEPT those in the second set).
  set<string> samples_to_keep_;
  set<string> samples_to_discard_;

  // The tags that can appear in vcf's filter column that should be removed.
  // Can specify the string "ALL" to discard variants that fail any filter.
  set<string> filter_column_tags_;

  // Filter based on vcf's INFO column.
  vector<KeyValueFilter> info_column_filters_;

  // Specific Cells to skip (e.g. put as N/A) based on rules governing
  // actual data in the cell.
  vector<KeyValueFilter> actual_data_filters_;

  // Filter based on vcf's QUALITY column.
  double min_quality_threshold_;

  // Filter based on Allele type (from vcf's ALT column).
  AlleleTypeFilter allele_filter_;

  // Filter based on position.
  PositionFilter position_filter_;
};

// A structure to hold the various criteria that will identify if the current
// Variant/SNP being processed belongs in the current group.
struct GroupingFilter {
  // The next two fields will be populated if SegmentType is LIST.
  // Key is Chromosome, Value is the set of positions that should be taken.
  map<Chromosome, set<int64_t>> variants_by_position_;
  set<string> snp_ids_;
  // First Coordinate is Chromosome, Second is gene start, Third is gene end.
  // TODO(PHB): Determine if any genes lie on multiple chromosomes, and/or have
  // multiple segments on the same chromosome. If so, I'll have to modify
  // data structure below.
  // This field will be populated if SegmentType is GENE.
  tuple<Chromosome, int64_t, int64_t> gene_position_; 
  // If SegmentType == WINDOW, we use window_size_ directly; so no need
  // to have a field here for this grouping type.
};

// A structure to hold all information required to perform grouping
// by gene.
struct GeneGroupingInfo {
  // Holds the positions of all genes. This data structure can be created
  // from a genepos file; see GetGenePositionMapFromFile() below.
  const map<string, tuple<Chromosome, int64_t, int64_t>>& name_to_position_;
  // The name of the current gene to be processed.
  string current_gene_;
  // TODO(PHB): Determine if this is needed/used by any API.
  // The index (line/row) of the gene to be processed, based on its location
  // in genepos_file_.
  //int current_gene_index_;
};

// A structure to hold all information required to perform grouping
// by list. The list can either specify a SNP (variant) by its
// (Chromosome, Position) or by it's SNP ID.
struct ListGroupingInfo {
  // Each element should be a filename; each file should have lines with format:
  //   CHROMOSOME:POSITION
  // or
  //   SNP_ID
  // Each file corresponds to one group.
  vector<string> list_files_;
  // Which file (index within list_files_) should be processed in this iteration.
  int current_file_index_;
};

// A structure to hold all information required to perform grouping
// by window.
struct WindowGroupingInfo {
  int window_size_;
  // Which window should be processed in current iteration.
  // TODO(PHB): Currently, this field is not used, as 'next_char_to_read_pos' is
  // used instead to skip to the appropriate line in the vcf file.
  //int current_window_index_;
};

// A structure for specifying how variants should be grouped.
struct OutputOptions {
  // The kind of data to store in the final Output table.
  SnpOutputType output_type_;

  // The kind of grouping mechanism to use.
  SegmentType segment_type_;

  // At most one of the following should be set (in correspondence to
  // segment_type_ above).
  GeneGroupingInfo by_gene_;
  WindowGroupingInfo by_window_;
  ListGroupingInfo by_list_;
};

// A structure holding the relevant information for each vcf file (e.g.
// the filename, filtering options, grouping options, and header_row).
struct VcfGroup {
  string file_;
  FilteringOptions filtering_options_;
  OutputOptions output_options_;
  vector<string> header_row_;
};


class VcfUtils {
 public:
  // Input:
  //  - file: The location of the vcf file
  //  - filtering_options: Filtering that should be done on rows/columns
  //  - header_row: The first line of the vcf file (title names of the columns);
  //                can be pointer to empty vector if Header has not yet been
  //                read (e.g. on a previous iteration.
  //  - next_char_to_read_pos: (Character) Position to begin reading from. If
  //                           negative or NULL, will read from beginning of
  //                           file (i.e. same as this being equal to zero).
  //                           If non-null, will also be set as an output
  //                           parameter, as the character to start from on the
  //                           next iteration.
  //  - output_options: How rows should be grouped/organized for generating
  //                    output tables/summary statistics; also includes the
  //                    kind of data (SnpOutputType) to be stored in each cell
  //                    Passed as pointer, as it will be modified (see Output)
  // Output:
  //  - snp_values_per_sample: Primary output table, w/ Rows=Samples, Cols=SNP.
  //                           Cell entries are of type SnpType (specific format
  //                           determined by output_options.output_type_)
  //  - next_char_to_read_pos: If only a chunk of data from the vcf file was
  //                           read (based on grouping/filtering_options), this
  //                           will store the character index of the first
  //                           character on the next line to be read on the
  //                           next iteration.
  //  - num_lines_processed: The number of lines processed in current iteration.
  //  - header_row: If header_row is empty when passed in, it will be populated,
  //                and returned so that future iterations can re-use it.
  //  - snp_info: Contains extra info (e.g. Major/Minor alleles, MAF, average
  //              Read-Depth, etc.) about each SNP. The index of a SNP in this
  //              vector matches the order of the SNPs in snp_values_per_sample,
  //              i.e. this also serves as the Header Column.
  static void ProcessVcfFile(const string& filename,
                      const FilteringOptions& filtering_options,
                      const OutputOptions& output_options,
                      vector<string>* header_row,
                      int64_t* next_char_to_read_pos,
                      int* num_lines_processed,
                      FilteringCounts* filtering_counts,
                      vector<vector<SnpType>>* snp_values_per_sample,
                      vector<SnpInfo>* snp_info);
  // Same as above, with no start/stop file read positions (read whole file).
  static void ProcessVcfFile(const string& filename,
                      const FilteringOptions& filtering_options,
                      const OutputOptions& output_options,
                      vector<string>* header_row,
                      int* num_lines_processed,
                      FilteringCounts* filtering_counts,
                      vector<vector<SnpType>>* snp_values_per_sample,
                      vector<SnpInfo>* snp_info) {
    ProcessVcfFile(
        filename, filtering_options, output_options, header_row, nullptr,
        num_lines_processed, filtering_counts, snp_values_per_sample, snp_info);
  }
  // Same as above, different API (file, filtering_options, and output_options
  // are condensed into a single 'vcf_group').
  static void ProcessVcfFile(VcfGroup* vcf_group,
                      int64_t* next_char_to_read_pos,
                      int* num_lines_processed,
                      FilteringCounts* filtering_counts,
                      vector<vector<SnpType>>* snp_values_per_sample,
                      vector<SnpInfo>* snp_info) {
    ProcessVcfFile(
        vcf_group->file_, vcf_group->filtering_options_,
        vcf_group->output_options_, &(vcf_group->header_row_),
        next_char_to_read_pos, num_lines_processed,
        filtering_counts, snp_values_per_sample, snp_info);
  }
  // Same as above, with no start/stop file read positions (read whole file).
  static void ProcessVcfFile(VcfGroup* vcf_group,
                      int* num_lines_processed,
                      FilteringCounts* filtering_counts,
                      vector<vector<SnpType>>* snp_values_per_sample,
                      vector<SnpInfo>* snp_info) {
    ProcessVcfFile(
        vcf_group->file_, vcf_group->filtering_options_,
        vcf_group->output_options_, &(vcf_group->header_row_),
        num_lines_processed, filtering_counts,
        snp_values_per_sample, snp_info);
  }
  // Same as above, but for multiple input vcf files.
  // TODO(PHB): Think about how to merge info from multiple files into one
  // output table; i.e. gonna have to do more than just iterate over each
  // element of vcf_groups, as is currently done below.
  static void ProcessVcfFile(const vector<VcfGroup*>& vcf_groups,
                      int64_t* next_char_to_read_pos,
                      int* num_lines_processed,
                      FilteringCounts* filtering_counts,
                      vector<vector<SnpType>>* snp_values_per_sample,
                      vector<SnpInfo>* snp_info) {
    for (VcfGroup* group : vcf_groups) {
      ProcessVcfFile(
          group, next_char_to_read_pos, num_lines_processed,
          filtering_counts, snp_values_per_sample, snp_info);
    }
  }
  // Same as above, with no start/stop file read positions (read whole file).
  static void ProcessVcfFile(const vector<VcfGroup*>& vcf_groups,
                             int* num_lines_processed,
                             FilteringCounts* filtering_counts,
                             vector<vector<SnpType>>* snp_values_per_sample,
                             vector<SnpInfo>* snp_info) {
    for (VcfGroup* group : vcf_groups) {
      ProcessVcfFile(group, num_lines_processed, filtering_counts,
                     snp_values_per_sample, snp_info);
    }
  }

  // Populates name_to_position (Key=GeneName, Value=(Chrom, Start_Pos, End_Pos))
  // based on genepos file, which is a file that optionally has Header Line
  // (which is ignored), followed by lines of format:
  //   GENE_ID CHROMOSOME START_POSITION END_POSITION [STRAND]
  // where last column (STRAND) is optional (defaults to '+' if absent).
  // All positions should be with respect to same reference genome as vcf file.
  // TODO(PHB): Currently, I only handle '+' strand. Determine if support for
  // '-' strand is needed, and if so, update all data structures and code to
  // offer flexibility to specify strand.
  static void GetGenePositionMapFromFile(
      const string& filename,
      map<string, tuple<Chromosome, int64_t, int64_t>>* gene_name_to_position);

  // File 'filename' has lines of format:
  //   Sample_ID
  // Populates samples_to_filter.
  // NOTE: How samples_to_filter is used (positive or negative filter, i.e.
  // whether these are the samples that are saved or the samples that are ignored)
  // is determined by calling function.
  static void SetSamplesToFilterFromFile(const string& filename,
                                         set<string>* samples_to_filter);

  // File 'filename' has lines of format:
  //   SNP_ID
  // Populates snps_to_filter.
  // NOTE: How snps_to_filter is used (positive or negative filter, i.e.
  // whether these are the snps that are saved or the snps that are ignored)
  // is determined by calling function.
  static void SetVariantsToFilterFromFile(const string& filename,
                                          set<string>* snps_to_filter);

  // File 'filename' has lines of format:
  //   CHROME POSITION
  // Or
  //   CHROME START_RANGE END_RANGE
  // Populates filter_positions.chromosome_position_pairs_to_filter_ from lines
  // of the first form, and filter_positions.chromosome_range_pairs_to_filter_
  // from lines of the second form.
  // NOTE: File should have lines only of one of the two forms, as at most of
  // chromosome_position_pairs_to_filter_ and chromosome_range_pairs_to_filter_
  // should be non-empty.
  // NOTE: How filter_positions is used (positive or negative filter, i.e.
  // whether these are the snps that are saved or the snps that are ignored)
  // is determined by calling function.
  static void SetVariantsToFilterFromFile(const string& filename,
                                          PositionFilter* filter_positions);

  // Attempts to parse chr_str as one of the Chromosome types; returns true
  // upon success, false otherwise.
  static bool ParseChromosome(const string& chr_str, Chromosome* chr);

  // Prints a string representation of chr.
  static string PrintChromosome(const Chromosome chr);

  // Attempts to parse input as a Nucleotide type.
  static bool ParseNucleotide(const string& input, Nucleotide* n);

  // Converts a Nucleotide enum into its string representative.
  static string GetNucleotideString(const Nucleotide n);

 private:
  // Processes one line (variant/SNP) in the vcf file:
  //   - Checks filtering_options to see if the variant should be kept
  //     (if not, updates filtering_counts appropriately)
  //   - Makes sure line is processed successfully (otherwise, error is
  //     recorded by the returned ProcessLineResponse)
  //   - Checks that this line belongs in the output, as specified by
  //     output_options
  //   - If line should be kept, updates snp_values_per_sample and snp_info
  //     appropriately.
  static ProcessLineResponse ProcessLine(
      const string& line,
      const int expected_num_cols,
      const vector<string>& header_row,
      const GroupingFilter& grouping_filter,
      const FilteringOptions& filtering_options,
      const SnpOutputType& output_type,
      FilteringCounts* filtering_counts,
      vector<SnpType>* snp_values_per_sample,
      vector<SnpInfo>* snp_info);

  // Determines if the current Segment (Window, Gene, or Variant List) is
  // complete, and returns true if so.
  static bool IsSegmentDone(const ProcessLineResponse& response_type,
                     const int num_variants_to_keep,
                     const int num_variants_taken,
                     const OutputOptions& output_options);

  // Attempts to parse line in one of two possible formats:
  //   CHROMOSOME:POSITION
  // or
  //   SNP_ID
  // If format one, populates filter->variants_by_position_; if format two,
  // populates filter->snp_ids_. If neither format, returns false.
  static bool AddVariantToFilter(const string& line, GroupingFilter* filter);

  // Populates variants_to_keep with the contents of the i^th file of
  // list_info.list_files_, where i = list_info.current_file_index_.
  static bool GetVariantsToKeep(const OutputOptions& output_options,
                         GroupingFilter* filter, int* num_variants_to_keep);

  // Parses SNP position from POSITION column of the vcf file.
  static bool GetPosition(const string& column, int64_t* position);

  // Returns the Nucleotide type that c represents.
  static Nucleotide GetNucleotide(const char c);

  // Parses Allele types from ALT and REF columns. Note that quad-allelic
  // variants are not handled (will return false).
  static bool GetAlleleTypes(
      const string& ref, const string& alt,
      Nucleotide* major, Nucleotide* minor, Nucleotide* triallele);

  // Returns whether position lies in any of the ranges in 'ranges', where
  // the ranges in 'ranges' are determined by (Start_Range, End_Range) :=
  // (Key, Value).
  static bool IsPositionInRange(
      const int64_t position, const map<int64_t, int64_t>& ranges);

  // Determines if Variant should be skipped, based on its CHROM and
  // POSITION columns, and on the PositionFilter options.
  static bool IsPositionFiltered(
      const Chromosome chromosome, const int64_t& position,
      const PositionFilter& position_filter);

  // Determines if Variant should be skipped, based on its REF and ALT columns
  // and allele_filter.
  static bool IsAlleleFiltered(
    const Nucleotide& alt, const Nucleotide& triallele,
    const double& minor_freq, const double& tri_freq,
    const AlleleTypeFilter& allele_filter);

  // Determines if Variant should be skipped, based on its AlleleType.
  static bool IsAlleleTypeFiltered(
      const Nucleotide& alt, const FilterType& filter_type,
      const set<AlleleType>& type_filter);

  // Returns the appropriate AlleleType for n.
  static AlleleType GetAlleleType(const Nucleotide n);

  // Determines if Variant should be skipped, based on its SNP_ID column
  // and sets snps_to_discard and snps_to_keep.
  static bool IsSnpFiltered(const string& snp_id,
                     const set<string>& snps_to_discard,
                     const set<string>& snps_to_keep);

  // Determines if Variant should be skipped, based on its QUAL column
  // and min_quality_threshold.
  static bool IsQualityFiltered(
      const string& quality, const double& min_quality_threshold);

  // Determines if Variant should be skipped, based on its FILTER column
  // value and the filtering options in filter_column_tags.
  static bool IsFilterFiltered(
      const string& column, const set<string>& filter_column_tags);

  // Determines if Variant should be skipped, based on its Key and Value
  // (from e.g. INFO Column) and the KeyValueFilters in filters.
  static bool IsKeyValueFiltered(
      const string& key, const string& value,
      const vector<KeyValueFilter>& filters);

  // Determines if Variant should be skipped, based on its INFO column
  // and the info_column_filters.
  // Also populates minor_freq and tri_freq.
  static bool IsInfoFiltered(
      const string& info, const vector<KeyValueFilter>& info_column_filters,
      double* minor_freq, double* tri_freq);

  // Returns SNP_TYPE_MAJOR if c is '0', SNP_TYPE_MINOR if c is '1', SNP_TYPE_TRI
  // if c is '2', and SNP_TYPE_UNKNOWN otherwise.
  static SnpType GetSnpType(const char c);

  // Processes the given data column (one of the Sample columns from vcf file).
  static bool ProcessDataColumn(
      const string& column, const int gt_index,
      const vector<KeyValueFilter>& data_filters,
      const SnpOutputType& output_type,
      FilteringCounts* filtering_counts,
      SnpType* snp_type);


  // Reads vcf file to find and extract its Header Row.
  static void GetHeaderRow(
      const string& filename,
      const FilteringOptions& filtering_options,
      vector<string>* header_row);

  // Determines if the given line has the format of a METADATA row of a vcf file.
  static bool IsMetadataLine(const string& line);

  // Attempts to process the current line as the Header Row. Returns true
  // if successfully parsed, false otherwise.
  static bool ProcessHeaderLine(const string& line,
                         const FilteringOptions& filtering_options,
                         vector<string>* header_row);

};

}  // namespace file_reader_utils

#endif
