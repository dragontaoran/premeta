// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "read_metaskat_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "StringUtils/string_utils.h"
//PHB#include "TestUtils/test_utils.h"

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
//PHBusing namespace test_utils;
using namespace std;

namespace premeta {

// Performs a standard Cyclic Redundancy Check (CRC) on the first len characters
// of buf. This matches the CRC used by MetaSKAT (the first 4 characters of each
// gene's 'block' of the MSSD file).
uint32_t xcrc32(const unsigned char* buf, int len) {
  uint32_t init = 0xffffffff;
  uint32_t crc = init;
  while (len--) {
    crc = (crc << 8) ^ crc32_table[((crc >> 24) ^ *buf) & 255];
    buf++;
  }
  return crc;
}

bool ParseMInfoHeader(
    const FileInfo& file_info,
    int* num_samples, int* num_genes, int* num_snps,
    int* num_unique_snps, int* num_header_lines) {
  ifstream file(file_info.name_.c_str());
  if (!file.is_open()) {
    cout << "ERROR in parsing MInfo Header: Unable to open '" << file_info.name_
         << "'." << endl;
    return false;
  }
  *num_header_lines = 0;
  *num_samples = -1;
  *num_genes = -1;
  *num_snps = -1;
  *num_unique_snps = -1;
  int n_all = -1;
  string line;
  while(getline(file, line)) {
    if (!HasPrefixString(line, file_info.comment_char_)) break;
    (*num_header_lines)++;
    string suffix;
    if (StripPrefixString(line, "#N.ALL=", &suffix)) {
      if (!Stoi(suffix, &n_all)) {
        cout << "ERROR in parsing MInfo Header: Unable to parse "
             << "number of samples: '" << suffix << "'. Aborting.\n";
        return false;
      }
      if (*num_samples != -1 && *num_samples != n_all) {
        cout << "WARNING in ParseMInfoHeader: Found mismatching N.ALL ("
             << n_all << ") and N (" << *num_samples << ")." << endl;
      }
    } else if (StripPrefixString(line, "#N=", &suffix)) {
      if (!Stoi(suffix, num_samples)) {
        cout << "ERROR in parsing MInfo Header: Unable to parse "
             << "number of samples: '" << suffix << "'. Aborting.\n";
        return false;
      }
      if (n_all != -1 && *num_samples != n_all) {
        cout << "WARNING while parsing MInfo Header: Found mismatching N.ALL ("
             << n_all << ") and N (" << *num_samples << ")." << endl;
      }
    } else if (StripPrefixString(line, "#nSNPs=", &suffix)) {
      if (!Stoi(suffix, num_snps)) {
        cout << "ERROR in parsing MInfo Header: Unable to parse "
             << "number of SNPs: '" << suffix << "'. Aborting.\n";
        return false;
      }
    } else if (StripPrefixString(line, "#nSNPs.unique=", &suffix)) {
      if (!Stoi(suffix, num_unique_snps)) {
        cout << "ERROR in parsing MInfo Header: Unable to parse "
             << "number of Unique SNPs: '" << suffix << "'. Aborting.\n";
        return false;
      }
    } else if (StripPrefixString(line, "#nSets=", &suffix)) {
      if (!Stoi(suffix, num_genes)) {
        cout << "ERROR in parsing MInfo Header: Unable to parse "
             << "number of genes: '" << suffix << "'. Aborting.\n";
        return false;
      }
    } else {
      cout << "WARNING while parsing MInfo Header: Unrecognized comment line:\n\t"
           << line << "\nwill be ignored." << endl;
    }
  }

  // Sanity-check we parsed requisite info (num samples, num snps, and
  // num unique snps) from header.
  if (*num_samples == -1) {
    cout << "ERROR in parsing MInfo Header: Unable to find number of "
         << "samples in file (there should be a line prefixed by "
         << "#N=" << ". Aborting.\n";
    return false;
  }
  if (*num_genes == -1) {
    cout << "ERROR in parsing MInfo Header: Unable to find number of "
         << "genes in file (there should be a line prefixed by "
         << "#nSets=" << ". Aborting.\n";
    return false;
  }
  if (*num_snps == -1) {
    cout << "ERROR in parsing MInfo Header: Unable to find number of "
         << "SNPs in file (there should be a line prefixed by "
         << "#nSNPs=" << ". Aborting.\n";
    return false;
  }
  if (*num_unique_snps == -1) {
    cout << "ERROR in parsing MInfo Header: Unable to find number of "
         << "Unique SNPs in file (there should be a line prefixed by "
         << "#nSNPs.unique=" << ". Aborting.\n";
    return false;
  }
  file.close();
  return true;
}

bool ParseMInfoContents(
    const int study_num,
    const FileInfo& file_info, const map<Position, SnpInfo>& allele_info,
    const int num_header_lines, const int num_samples,
    const int num_genes, const int num_snps, const int num_unique_snps,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<string, uint64_t>* start_lines,
    map<Position, SnpInfo>* snp_info) {
  string allele_one = kMInfoHeaderMajorAllele;
  string allele_two = kMInfoHeaderMinorAllele;
  map<string, int> col_indices;
  col_indices.insert(make_pair(kMInfoHeaderSetId, -1));
  col_indices.insert(make_pair(kMInfoHeaderSetIdNumeric, -1));
  col_indices.insert(make_pair(kMInfoHeaderSnpId, -1));
  col_indices.insert(make_pair(kMInfoHeaderScore, -1));
  col_indices.insert(make_pair(kMInfoHeaderMaf, -1));
  col_indices.insert(make_pair(kMInfoHeaderMissingRate, -1));
  col_indices.insert(make_pair(kMInfoHeaderMajorAllele, -1));
  col_indices.insert(make_pair(kMInfoHeaderMinorAllele, -1));
  col_indices.insert(make_pair(kMInfoHeaderPass, -1));
  col_indices.insert(make_pair(kMInfoHeaderStartPos, -1));
  if (!GetColumnIndices(file_info, &col_indices)) {
    // Try old MetaSKAT format.
    col_indices.clear();
    allele_one = "Allele1";
    allele_two = "Allele2";
    col_indices.insert(make_pair(kMInfoHeaderSetId, -1));
    col_indices.insert(make_pair(kMInfoHeaderSetIdNumeric, -1));
    col_indices.insert(make_pair(kMInfoHeaderSnpId, -1));
    col_indices.insert(make_pair(kMInfoHeaderScore, -1));
    col_indices.insert(make_pair(kMInfoHeaderMaf, -1));
    col_indices.insert(make_pair(kMInfoHeaderMissingRate, -1));
    col_indices.insert(make_pair(allele_one, -1));
    col_indices.insert(make_pair(allele_two, -1));
    col_indices.insert(make_pair(kMInfoHeaderPass, -1));
    col_indices.insert(make_pair(kMInfoHeaderStartPos, -1));
    if (!GetColumnIndices(file_info, &col_indices)) {
      cout << "ERROR. Unable to parse MetaSKAT's .MInfo file '"
           << file_info.name_ << "': Could not find all the expected "
           << "columns: (" << kMInfoHeaderSetId << ", "
           << kMInfoHeaderSetIdNumeric << ", " << kMInfoHeaderSnpId
           << ", " << kMInfoHeaderScore << ", " << kMInfoHeaderMaf
           << ", " << kMInfoHeaderMissingRate << ", "
           << kMInfoHeaderMajorAllele << ", " << kMInfoHeaderMinorAllele
           << ", " << kMInfoHeaderPass << ", " << kMInfoHeaderStartPos
           << "). Check version used to generate .MInfo file "
           << "(PreMeta is compatible with version "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_METASKAT)
           << endl;
      return false;
    }
  }

  ReadCsvInput input;
  input.filename_ = file_info.name_;
  input.delimiters_.insert(file_info.delimiter_);
  input.comment_char_ = file_info.comment_char_;
  input.has_header_ = true;
  vector<pair<int, GenericDataType>>& column_types = input.columns_to_read_;
  column_types.push_back(make_pair(  // Gene Name.
      col_indices[kMInfoHeaderSetId], GenericDataType::STRING));
  column_types.push_back(make_pair(  // Gene Index.
      col_indices[kMInfoHeaderSetIdNumeric], GenericDataType::INT));
  column_types.push_back(make_pair(  // Position.
      col_indices[kMInfoHeaderSnpId], GenericDataType::STRING));
  column_types.push_back(make_pair(  // Score (U_STAT).
      col_indices[kMInfoHeaderScore], GenericDataType::DOUBLE));
  column_types.push_back(make_pair(  // MAF.
      col_indices[kMInfoHeaderMaf], GenericDataType::DOUBLE));
  column_types.push_back(make_pair(  // MISSING_RATE.
      col_indices[kMInfoHeaderMissingRate], GenericDataType::DOUBLE));
  column_types.push_back(make_pair(  // Major Allele.
      col_indices[allele_one], GenericDataType::STRING));
  column_types.push_back(make_pair(  // Minor Allele.
      col_indices[allele_two], GenericDataType::STRING));
  column_types.push_back(make_pair(  // PASS.
      col_indices[kMInfoHeaderPass], GenericDataType::STRING));
  column_types.push_back(make_pair(  // StartPos.
      col_indices[kMInfoHeaderStartPos], GenericDataType::UINT_64));
  ReadCsvOutput output; 
  if (!CsvUtils::ReadCsv(input, &output)) {
    cout << "ERROR in parsing MInfo file: Unable to read MInfo file "
         << "'" << file_info.name_ << "'. Aborted with error:\n"
         << output.error_msg_ << endl;
    return false;
  }
  const vector<vector<GenericDataHolder>>& parsed_minfo_file = output.output_;

  // Sanity-check we read expected number of SNP positions.
  if (parsed_minfo_file.size() != num_snps) {
    cout << "ERROR in parsing MInfo file: Expected " << num_snps
         << "rows, but found " << parsed_minfo_file.size() << endl;
    return false;
  }

  // Iterate through parsed_minfo_file, populating scores.
  string current_gene_name = "";
  int current_gene_id = -1;
  int num_duplicate_snps = 0;
  for (int i = 0; i < parsed_minfo_file.size(); ++i) {
    const vector<GenericDataHolder>& current_row = parsed_minfo_file[i];

    // Parse Position.
    Position pos;
    if (!ParsePosition(current_row[2].str_, &pos)) {
      // Unable to parse Position into CHR:POS format; store the Position string
      // into the snp_id_ field.
      pos.snp_id_ = current_row[2].str_;
    }
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(pos);
      if (allele_info_itr != allele_info.end()) {
        pos = allele_info_itr->second.pos_;
      }
    }

    // Check that REF/ALT alleles match existing entry for this SNP. Otherwise,
    // either add this SNP to snp_to_ref_alt_and_study_ or snp_to_excluding_study_,
    // as appropriate.
    tuple<string, string, set<int>>* snp_ref_alt_info =
        FindOrInsert(
            pos, *snp_to_ref_alt_and_study,
            make_tuple(current_row[6].str_, current_row[7].str_, set<int>()));
    if (get<0>(*snp_ref_alt_info) != current_row[6].str_ ||
        get<1>(*snp_ref_alt_info) != current_row[7].str_) {
      if (get<1>(*snp_ref_alt_info) != current_row[6].str_ ||
          get<0>(*snp_ref_alt_info) != current_row[7].str_) {
        // This is a more aggregious difference than just a swapping of
        // REF <-> ALT. Add snp to snp_to_excluding_study_.
        set<int>* excluding_studies =
            FindOrInsert(pos, *snp_to_excluding_study, set<int>());
        excluding_studies->insert(study_num);
        continue;
      } else {
        get<2>(*snp_ref_alt_info).insert(study_num);
      }
    }

    // Add/insert this position for this gene.
    const string& gene_name = current_row[0].str_;
    const int gene_id = current_row[1].int_;
    if ((gene_name == current_gene_name) != (gene_id == current_gene_id)) {
      cout << "ERROR in parsing MInfo file on line " << (i + num_header_lines + 2)
           << ": Previous gene (" << current_gene_name << ") and current gene ("
           << gene_name << (gene_name == current_gene_name ? ") " : ") don't ")
           << "match, but previous gene id (" << current_gene_id
           << (gene_id == current_gene_id ? ") matches " : ") doesn't match")
           << "current gene id (" << gene_id << "). Aborting." << endl;
      return false;
    }
    vector<Position>* snps_on_gene =
        FindOrInsert(gene_name, *gene_to_snp_positions, vector<Position>());
    snps_on_gene->push_back(pos);

    // Add SNP Info for this position.
    pair<map<Position, SnpInfo>::iterator, bool> snp_info_itr =
        snp_info->insert(make_pair(pos, SnpInfo()));
    SnpInfo& info = snp_info_itr.first->second;
    const int mac =
        round(current_row[4].dbl_ * 2.0 * num_samples * (1.0 - current_row[5].dbl_));
    if (!snp_info_itr.second) {
      // We've already seen this position. Sanity check values match.
      if (!FloatEq(info.u_stat_, current_row[3].dbl_) ||
          !FloatEq(info.maf_, current_row[4].dbl_) ||
          info.mac_ != mac ||
          info.missing_rate_ != current_row[5].dbl_ ||
          info.major_allele_ != current_row[6].str_ ||
          info.minor_allele_ != current_row[7].str_) {
        cout << "ERROR in parsing MInfo file on line " << (i + num_header_lines + 2)
             << ": Already seen Position " << PrintPosition(pos) << " before, but "
             << "previous values don't match current values: Old score ("
             << info.u_stat_ << ") vs. New score (" << current_row[3].dbl_
             << "); Old MAF (" << info.maf_ << ") vs. New MAF ("
             << current_row[4].dbl_ << "); Old MAC (" << info.mac_
             << ") vs. New MAC (" << mac << "); New Missing Rate ("
             << current_row[5].dbl_ << ") vs Old Missing Rate ("
             << info.missing_rate_ << "); New Major Allele ("
             << current_row[6].str_ << ") vs Old Major Allele ("
             << info.major_allele_ << "); New Minor Allele ("
             << current_row[7].str_ << ") vs Old Minor Allele ("
             << info.minor_allele_ << ")." << endl;
        return false;
      }
      num_duplicate_snps++;
    } else {
      info.u_stat_ = current_row[3].dbl_;
      info.maf_ = current_row[4].dbl_;
      info.mac_ = mac;
      info.missing_rate_ = current_row[5].dbl_;
      info.major_allele_ = current_row[6].str_;
      info.minor_allele_ = current_row[7].str_;
      if (num_snps == 0 ||
          info.maf_ < 1.0 / (4.0 * num_snps) ||
          info.maf_ > 1.0 - (1.0 / (4.0 * num_snps))) {
        monomorphic_snps->insert(pos);
      }      
    }

    // Get start position (line number w.r.t. MSSD file) for this gene.
    // Note that MInfo file has 1 SNP per line, and each line contains the
    // gene's start position (line number), which is thus repeated on each
    // line (SNP) of the MInfo file that belongs to the same gene. We only
    // need to record the gene's start position once (e.g. on the first
    // SNP we read for the gene); and then sanity-check all subsequent
    // SNPs on that gene indicate the same start position.
    pair<map<string, uint64_t>::iterator, bool> start_lines_inserter =
        start_lines->insert(make_pair(gene_name, current_row[9].uint64_));
    if (!start_lines_inserter.second &&
        start_lines_inserter.first->second != current_row[9].uint64_) {
      cout << "ERROR in parsing MInfo file on line " << (i + num_header_lines + 2)
           << ": Position " << PrintPosition(pos) << " on gene '"
           << gene_name << "' indicates start line "
           << current_row[9].uint64_ << ", but a previous SNP Position "
           << "indicated that gene starts on line "
           << start_lines_inserter.first->second << ". Aborting." << endl;
      return false;
    }

    // TODO(PHB): Columns 7-10, 12 aren't needed for generating MASS output. Determine
    // if they are needed for any of the output types, and if not, remove reading
    // these columns above.

    current_gene_name = gene_name;
    current_gene_id = gene_id;
  }

  // Sanity-check the right number of unique SNP positions were found.
  if (num_duplicate_snps != num_snps - num_unique_snps) {
    cout << "ERROR in parsing MInfo file: Expected " << num_unique_snps
         << " Unique SNP positions, but found " << num_snps - num_duplicate_snps
         << endl;
    return false;
  }

  // Sanity-check the right number of Genes were found.
  if (num_genes != gene_to_snp_positions->size()) {
    cout << "ERROR in parsing MInfo file: Expected " << num_genes
         << " genes, but found " << gene_to_snp_positions->size() << endl;
    return false;
  }

  return true;
}

bool ParseMInfoFile(
    const int study_num,
    const FileInfo& file_info, const map<Position, SnpInfo>& allele_info,
    int* num_samples,
    set<Position>* monomorphic_snps,
    map<Position, tuple<string, string, set<int>>>* snp_to_ref_alt_and_study,
    map<Position, set<int>>* snp_to_excluding_study,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<string, uint64_t>* start_lines,
    map<Position, SnpInfo>* snp_info) {
  int num_header_lines = -1;
  int num_snps = -1;
  int num_unique_snps = -1;
  int num_genes = -1;
  if (!ParseMInfoHeader(
          file_info, num_samples,
          &num_genes, &num_snps, &num_unique_snps, &num_header_lines)) {
    return false;
  }

  if (!ParseMInfoContents(
        study_num, file_info, allele_info, num_header_lines,
        *num_samples, num_genes, num_snps, num_unique_snps,
        monomorphic_snps, snp_to_ref_alt_and_study, snp_to_excluding_study,
        gene_to_snp_positions, start_lines, snp_info)) {
    return false;
  }

  return true;
}

bool ParseMssdFile(
    const int study_num, const string& input_file,
    const map<Position, set<int>>& snp_to_excluding_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<string, uint64_t>& start_lines,
    map<pair<Position, Position>, double>* covariances) {

  ifstream file(input_file, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR in parsing Mssd file: Unable to open file: " << input_file
         << ". Aborting." << endl;
    return false;
  }

  /*
  // PHB Temp: writing MSSD file, as a sanity check.
  const string PHB_outfile = "foo_mssd_input_as_chars.txt";
  cout << "\nPrinting MSSD file as characters (int) to '"
       << PHB_outfile << "'" << endl;
  int64_t file_size = file.tellg();
  char* in_buffer = new char[file_size];
  file.seekg(0, ios::beg);
  // Read block.
  file.read(in_buffer, file_size);
  ofstream outfile;
  outfile.open(PHB_outfile.c_str());
  for (int64_t i = 0; i < file_size; ++i) {
    unsigned char c = in_buffer[i];
    outfile << "Char " << i << ": " << (int)c << endl;
  }
  outfile.close();
  cout << "\nDone Printing MSSD file." << endl;
  // END PHB Temp.
  */

  // Parse MSSD file by parsing contiguous gene blocks.
  int num_genes_processed = 0;
  const int percentage = start_lines.size() / 100;
  for (map<string, uint64_t>::const_iterator itr = start_lines.begin();
       itr != start_lines.end(); ++itr) {
    ++num_genes_processed;
    const uint64_t& start_pos = itr->second;
    const string& current_gene = itr->first;

    // Get SNP positions in this gene.
    map<string, vector<Position>>::const_iterator gene_pos_itr =
        gene_to_snp_positions.find(current_gene);
    if (gene_pos_itr == gene_to_snp_positions.end()) {
      cout << "ERROR in parsing MSSD file: Unable to find SNP positions for gene '"
           << current_gene << ". Aborting." << endl;
      return false;
    }
    const vector<Position>& snp_positions = gene_pos_itr->second;
    const int num_snps = snp_positions.size();

    // For upper-triangular (num_snps x num_snps) matrix, we expect
    // (num_snps * (num_snps + 1)) / 2 entries, with each entry a float.
    const int num_upper_tri_entries = (num_snps * (num_snps + 1)) / 2; 
    const int num_upper_tri_bytes = num_upper_tri_entries * sizeof(float);

    // Based on comment in MetaSKAT's MatFile.h, it seems like if you want
    // to restrict to 5,000 x 5,000 matrices, then the number of upper-triangle
    // entries is 5k * (5k + 1) / 2 = 12,502,500, and then the number of bytes
    // needed to represent this (4 bytes per entry) is 50,010,000. So I should
    // probably set kMaxMatrixSize to 50010000 and compare num_upper_tri_bytes
    // to this, but I'll mimic the original code in MatFile.cpp.
    if (num_upper_tri_bytes > kMaxMatrixSize){
      cout << "ERROR in parsing Mssd file: Unable to handle gene '" << current_gene
           << "', which has more than 5,000 SNPs (found " << num_snps
           << " SNPs). Aborting." << endl;
      return false;
    } else if (num_upper_tri_bytes > kMaxBufferSize) {
      cout << "ERROR in parsing Mssd file: Unable to handle gene '" << current_gene
           << "', which has too many SNPs (" << num_upper_tri_bytes
           << ") and thus the gene block in MSSD exceeds kMaxBufferSize ("
           << kMaxBufferSize << ") bytes. Aborting." << endl;
      return false;
    }
    
    file.seekg(start_pos, ios::beg);

    // Read first 4 bytes, which is the CRC of the covariance matrix.
    char crc_buffer[10];
    file.read(crc_buffer, sizeof(unsigned int)); // first 4 bytes
    uint32_t crc1;
    memcpy(&crc1, crc_buffer, sizeof(unsigned int));

    // Read gene block.
    char* gene_block = new char[num_upper_tri_bytes];
    file.read(gene_block, num_upper_tri_bytes);
    const uint32_t crc2 =
        xcrc32((unsigned char*) gene_block, num_upper_tri_bytes);

    // Check CRC.
    if (crc1 != crc2) {
      cout << "ERROR in parsing Mssd file: Failed CRC for gene '"
           << current_gene << "' at gene group number " << num_genes_processed
           << ", which has start position "
           << start_pos << " and num_snps " << num_snps << ". First 4 bytes: "
           << crc1 << ", CRC of Covariance Matrix: " << crc2
           << ". Aborting." << endl;
      return false;
    }
   
    // Parse bytes as Float values.
    // NOTE: this is not strictly necessary, as I can just parse them as I use
    // them (see commented-out code below) when populating 'covariances', but
    // doing so here doesn't cost anything, and makes the code cleaner.
    float* covariance_values = new float[num_upper_tri_bytes];
    memcpy(covariance_values, gene_block, num_upper_tri_bytes);

    // Copy covariances to 'covariances'.
    int k = 0;
    for (int i = 0; i < num_snps; i++) {
      const Position& row_pos = snp_positions[i];
      // Check that this SNP shouldn't be skipped (due to mismatched REF/ALT
      // alleles).
      if (ShouldExcludeSnp(study_num, row_pos, snp_to_excluding_study)) {
          continue;
      }
      for (int j = i; j < num_snps; j++) {
        const Position& col_pos = snp_positions[j];
        // Check that this SNP shouldn't be skipped (due to mismatched REF/ALT
        // alleles).
        if (ShouldExcludeSnp(study_num, col_pos, snp_to_excluding_study)) {
            continue;
        }
        // NOTE: can uncomment below code if you don't want to use the
        // temporary 'covariance_values' structure.
        // const double cov = static_cast<double>(ParseFloat(gene_block, k));
        // k += 4;
        const double cov = static_cast<double>(covariance_values[k]);
        k++;

        pair<map<pair<Position, Position>, double>::iterator, bool> cov_itr =
            covariances->insert(make_pair(make_pair(row_pos, col_pos), cov));
        if (!cov_itr.second && !FloatEq(cov_itr.first->second, cov)) {
          cout << "ERROR in parsing Mssd file: Covariance for Positions "
               << PrintPosition(row_pos) << ", " << PrintPosition(col_pos)
               << " appears twice, with different values. Aborting." << endl;
          return false;
        }
      }
    }
    delete gene_block;
    delete covariance_values;
  }

  return true;
}

}  // namespace premeta
