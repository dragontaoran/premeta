// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "write_metaskat_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_metaskat_utils.h"

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
 
void FilterFlaggedSnpPositions(
    const vector<Position>& orig_gene_snp_positions,
    const set<Position>& snps_to_skip,
    vector<Position>* gene_snp_positions) {
  for (const Position& pos : orig_gene_snp_positions) {
    if (snps_to_skip.find(pos) == snps_to_skip.end()) {
      gene_snp_positions->push_back(pos);
    }
  }
}

bool WriteMssdFile(
    const int study_num, const double& rescale,
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<pair<Position, Position>, double>& covariances,
    const string& outfile) {
  ofstream out_file;
  out_file.open(outfile.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // TODO(PHB): Figure out why character '17' is written first.
  char PHB_foo[1];
  int PHB_bar = 17;
  memcpy(PHB_foo, &PHB_bar, sizeof(int));
  out_file.write(PHB_foo, sizeof(char));

  // Loop through all genes, writing covariance matrices for each.
  for (const auto& gene_itr : gene_to_snp_positions) {
    const vector<Position>& orig_gene_snp_positions = gene_itr.second;
    vector<Position> gene_snp_positions;
    FilterFlaggedSnpPositions(
        orig_gene_snp_positions, snps_to_skip, &gene_snp_positions);
    const int num_snps_on_gene = gene_snp_positions.size();
   
    // Loop through all pairs of SNPs on this gene, creating the
    // (lower-triangular) covariance matrix via 'covariances'.
    const int num_upper_tri_entries =
        (num_snps_on_gene * (1 + num_snps_on_gene)) / 2;
    float m_buffer_float[num_upper_tri_entries];
    int k = 0;
    for (int i = 0; i < num_snps_on_gene; ++i) {
      const Position& col_pos = gene_snp_positions[i];

      // Check if SNP has REF/ALT swapped from golden file.
      const bool is_ref_alt_swapped =
          IsSnpRefAltSwapped(study_num, col_pos, snp_to_ref_alt_and_study);

      for (int j = i; j < num_snps_on_gene; ++j) {
        const Position& row_pos = gene_snp_positions[j];

        const bool is_ref_alt_swapped_two =
          IsSnpRefAltSwapped(study_num, row_pos, snp_to_ref_alt_and_study);
        const double swapped_multiplier =
            (is_ref_alt_swapped && is_ref_alt_swapped_two) ? 1.0 :
            (is_ref_alt_swapped || is_ref_alt_swapped_two) ? -1.0 : 1.0;

        double covariance;
 
        const bool is_monomorphic =
            monomorphic_snps.find(col_pos) != monomorphic_snps.end() ||
            monomorphic_snps.find(row_pos) != monomorphic_snps.end();
        if (is_monomorphic ||
            !LookupCovariance(col_pos, row_pos, covariances, &covariance)) {
          // There may be times when we don't have covariance, e.g. if
          // converting from RareMetalWorker to MetaSKAT, and if the distance
          // between these 2 SNPs is greater than the Window size. For such
          // cases, indicate zero covariance.
          covariance = 0.0;
        }
        // Convert double to float (MetaSKAT does this to save space), and
        // store in m_buffer_float.
        m_buffer_float[k] =
            (float)(swapped_multiplier * covariance / (rescale * rescale));
        k++;
      }
    }
    
    int length = sizeof(float) * num_upper_tri_entries;
    char m_buffer_crc[length];
    memcpy(m_buffer_crc, m_buffer_float, length);

    uint32_t crc = xcrc32((unsigned char*) m_buffer_crc, length);
    char crc_c[4];
    memcpy(crc_c, &crc, 4);

    out_file.write(crc_c, sizeof(uint32_t));
    out_file.write((char *) m_buffer_crc, length);
  }

  out_file.close();
  return true;
}

bool WriteMInfoFile(
    const int study_num, const int num_samples, const int num_snps,
    const double& rescale, const double& sigma_sq,
    const set<Position>& snps_to_skip,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const string& outfile) {
  ofstream out_file;
  out_file.open(outfile.c_str());
  if (!out_file.is_open()) {
    cout << "ERROR in trying to write file '" << outfile
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  // Print Meta-data.
  out_file << "#N.ALL=" << num_samples << endl;
  out_file << "#N=" << num_samples << endl;
  out_file << "#nSets=" << gene_to_snp_positions.size() << endl;
  out_file << "#nSNPs=" << num_snps << endl;
  out_file << "#nSNPs.unique=" << snp_info.size() << endl;

  // Print Header row.
  out_file << "SetID" << kMInfoDefaultDelimiter << "SetID_numeric"
           << kMInfoDefaultDelimiter << "SNPID Score"
           << kMInfoDefaultDelimiter << "MAF" << kMInfoDefaultDelimiter
           << "MissingRate" << kMInfoDefaultDelimiter
           << "MajorAllele" << kMInfoDefaultDelimiter << "MinorAllele"
           << kMInfoDefaultDelimiter << "PASS" << kMInfoDefaultDelimiter
           << "StartPOS" << kMInfoDefaultDelimiter << "StartPOSPermu" << endl;

  // Print Rest of File.
  int gene_index = 1;
  int start_pos = 1;
  for (const auto& gene_itr : gene_to_snp_positions) {
    const vector<Position>& orig_gene_snp_positions = gene_itr.second;
    vector<Position> gene_snp_positions;
    FilterFlaggedSnpPositions(
        orig_gene_snp_positions, snps_to_skip, &gene_snp_positions);
    const int num_snps_on_gene = gene_snp_positions.size();
    for (int i = 0; i < num_snps_on_gene; ++i) {
      const Position& snp_pos = gene_snp_positions[i];

      // Check if SNP has REF/ALT swapped from golden file.
      const bool is_ref_alt_swapped =
          IsSnpRefAltSwapped(study_num, snp_pos, snp_to_ref_alt_and_study);

      map<Position, SnpInfo>::const_iterator snp_info_itr =
          snp_info.find(snp_pos);
      if (snp_info_itr == snp_info.end()) {
        cout << "ERROR in writing MInfo file: For SNP number " << (i + 1)
             << " on gene '" << gene_itr.first << "', unable to find any "
             << "information for this SNP. Aborting." << endl;
        return false;
      }

      // Gene name.
      out_file << gene_itr.first << kMInfoDefaultDelimiter;

      // Gene index.
      out_file << gene_index << kMInfoDefaultDelimiter;

      // SNP Position.
      out_file << PrintPosition(snp_pos) << kMInfoDefaultDelimiter;

      // Score.
      if (monomorphic_snps.find(snp_pos) != monomorphic_snps.end()) {
        out_file << 0.0;
      } else if (is_ref_alt_swapped) {
        out_file << (-1.0 * snp_info_itr->second.u_stat_ / rescale);
      } else {
        out_file << (snp_info_itr->second.u_stat_ / rescale);
      }
      out_file << kMInfoDefaultDelimiter;

      // MAF.
      if (is_ref_alt_swapped) {
        out_file << (1.0 - snp_info_itr->second.maf_) << kMInfoDefaultDelimiter;
      } else {
        out_file << snp_info_itr->second.maf_ << kMInfoDefaultDelimiter;
      }

      // MissingRate.
      if (snp_info_itr->second.missing_rate_ < 0.0) {
        if (snp_info_itr->second.num_non_missing_ >= 0) {
          // Set Missing Rate from N_INFO, if available.
          out_file << (1.0 - static_cast<double>(snp_info_itr->second.num_non_missing_) /
                       static_cast<double>(num_samples))
                   << kMInfoDefaultDelimiter;
        } else {
          // Set default value of 0 when we don't have missing_rate info.
          out_file << "0" << kMInfoDefaultDelimiter;
        }
      } else {
        out_file << snp_info_itr->second.missing_rate_
                 << kMInfoDefaultDelimiter;
      }

      // Major and Minor Alleles.
      if (is_ref_alt_swapped) {
        out_file << snp_info_itr->second.minor_allele_
                 << kMInfoDefaultDelimiter
                 << snp_info_itr->second.major_allele_
                 << kMInfoDefaultDelimiter;
      } else {
        out_file << snp_info_itr->second.major_allele_
                 << kMInfoDefaultDelimiter
                 << snp_info_itr->second.minor_allele_
                 << kMInfoDefaultDelimiter;
      }

      // PASS
      out_file << "PASS" << kMInfoDefaultDelimiter;

      // StartPOS
      out_file << start_pos << kMInfoDefaultDelimiter;

      // StartPOSPermu
      out_file << "0" << endl;
    }

    start_pos += 4 + (2 * (num_snps_on_gene) * (num_snps_on_gene + 1));
    ++gene_index;
  }
  out_file.close();
  return true;
}

}  // namespace premeta
