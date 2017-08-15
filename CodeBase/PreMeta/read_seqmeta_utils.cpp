// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "read_seqmeta_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_r_binary_utils.h"
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

bool ParseSeqMetaRDataMainObjectAttributeOne(
    const RDataPairList* first_attr, int* num_genes) {
  if (first_attr == nullptr || first_attr->tag_ == nullptr ||
      first_attr->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "first Attribute has unexpected format. Aborting." << endl;
    return false;
  }
  if (*(first_attr->tag_) != "dim") {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "first Attribute has tag '" << *(first_attr->tag_)
         << "' (expected 'dim'). Aborting." << endl;
    return false;
  }
  const RDataNode& first_attr_obj_node = *(first_attr->obj_);
  if (first_attr_obj_node.tag_ != nullptr ||
      first_attr_obj_node.attr_ != nullptr ||
      first_attr_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "first Attribute's object has unexpected format. Aborting." << endl;
    return false;
  }
  const RDataObject& first_attr_obj = *(first_attr_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, first_attr_obj) ||
      first_attr_obj.int_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "first Attribute's Object has unexpected format. Aborting." << endl;
    return false;
  }
  *num_genes = (*(first_attr_obj.int_vec_))[0];
  return true;
}

bool ParseSeqMetaRDataMainObjectAttributeTwo(
    const bool old_version, const int num_genes, const string& tag_name,
    const RDataPairList* second_attr, vector<string>* genes) {
  if (second_attr == nullptr || second_attr->tag_ == nullptr ||
      second_attr->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Object representing "
         << "the Main Attribute's gene names has unexpected format. Aborting."
         << endl;
    return false;
  }
  if (*(second_attr->tag_) != tag_name) {
    cout << "\nERROR in parsing SeqMeta RData object: Object representing "
         << "the Main Attribute's gene names has wrong tag: Expected '"
         << tag_name << "', found '" << *(second_attr->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& second_attr_obj_node = *(second_attr->obj_);
  if (second_attr_obj_node.tag_ != nullptr ||
      second_attr_obj_node.attr_ != nullptr ||
      second_attr_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Object representing "
         << "the Main Attribute's gene names has unexpected format. Aborting."
         << endl;
    return false;
  }

  // Older versions of SeqMeta store the gene names directly as the Pair-List's
  // object, while newer versions (i.e. those that PreMeta supports) have an
  // extra layer (Pair-List's object is a list (of one element), and that one
  // element is the object that stores the gene names).
  RDataObject second_attr_list_obj;
  if (!old_version) {
    const RDataObject& second_attr_obj = *(second_attr_obj_node.obj_);
    if (!ExpectObjectType(SEXPTYPE_LIST, second_attr_obj) ||
        second_attr_obj.list_vec_->size() != 1) {
      cout << "\nERROR in parsing SeqMeta RData object: Object representing "
           << "the Main Attribute's gene names has unexpected format. Aborting."
           << endl;
      return false;
    }
    const RDataNode* second_attr_list = (*(second_attr_obj.list_vec_))[0];

    if (second_attr_list == nullptr || second_attr_list->tag_ != nullptr ||
        second_attr_list->attr_ != nullptr || second_attr_list->obj_ == nullptr) {
      cout << "\nERROR in parsing SeqMeta RData object: Object representing "
           << "the Main Attribute's gene names has unexpected format. Aborting."
           << endl;
      return false;
    }
    second_attr_list_obj = *(second_attr_list->obj_);
  } else {
    second_attr_list_obj = *(second_attr_obj_node.obj_);
  }

  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, second_attr_list_obj) ||
      (num_genes >= 0 && second_attr_list_obj.str_vec_->size() != num_genes)) {
    cout << "\nERROR in parsing SeqMeta RData object: Object representing "
         << "the Main Attribute's gene names is a list with unexpected "
         << "number of genes (expected " << num_genes << " genes based on "
         << "Main Attribute's 2nd field, found "
         << second_attr_list_obj.str_vec_->size()
         << "). Aborting." << endl;
    return false;
  }
  for (const string& gene : *(second_attr_list_obj.str_vec_)) {
    genes->push_back(gene);
  }
  return true;
}

bool ParseSeqMetaRDataMainObjectAttributeThree(
    const RDataPairList* third_attr) {
  if (third_attr == nullptr || third_attr->tag_ == nullptr ||
      third_attr->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "third Attribute has unexpected format. Aborting." << endl;
    return false;
  }
  if (*(third_attr->tag_) != "family") {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "third Attribute has tag '" << *(third_attr->tag_)
         << "' (expected 'family'). Aborting." << endl;
    return false;
  }
  const RDataNode& third_attr_obj_node = *(third_attr->obj_);
  if (third_attr_obj_node.tag_ != nullptr ||
      third_attr_obj_node.attr_ != nullptr ||
      third_attr_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "third Attribute's object has unexpected format. Aborting." << endl;
    return false;
  }
  const RDataObject& third_attr_obj = *(third_attr_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, third_attr_obj) ||
      third_attr_obj.str_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "third Attribute's Object has unexpected format. Aborting." << endl;
    return false;
  }
  if ((*(third_attr_obj.str_vec_))[0] != "gaussian") {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "third Attribute's Object should be 'gaussian'; but found '"
         << (*(third_attr_obj.str_vec_))[0] << "'. Aborting." << endl;
    return false;
  }
  return true;
}

bool ParseSeqMetaRDataMainObjectAttributeFour(
    const RDataPairList* fourth_attr, bool* is_skat_cohort_version) {
  if (fourth_attr == nullptr || fourth_attr->tag_ == nullptr ||
      fourth_attr->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "fourth Attribute has unexpected format. Aborting." << endl;
    return false;
  }
  if (*(fourth_attr->tag_) != "class") {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "fourth Attribute has tag '" << *(fourth_attr->tag_)
         << "' (expected 'class'). Aborting." << endl;
    return false;
  }
  const RDataNode& fourth_attr_obj_node = *(fourth_attr->obj_);
  if (fourth_attr_obj_node.tag_ != nullptr ||
      fourth_attr_obj_node.attr_ != nullptr ||
      fourth_attr_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "fourth Attribute's object has unexpected format. Aborting." << endl;
    return false;
  }
  const RDataObject& fourth_attr_obj = *(fourth_attr_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, fourth_attr_obj) ||
      fourth_attr_obj.str_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "4th Attribute's Object has unexpected format. Aborting." << endl;
    return false;
  }
  if ((*(fourth_attr_obj.str_vec_))[0] == "skatCohort") {
    cout << "\nWARNING: While parsing SeqaMeta RData object, main attribute's "
         << "tag should be 'seqMeta', but found 'skatCohort'. This indicates "
         << "the RData object was constructed using an older (unsupported by "
         << "PreMeta) version of seqMeta. PreMeta will attempt to parse "
         << "this RData object as if it came from a supported version of "
         << "seqMeta." << endl;
    *is_skat_cohort_version = true;
  } else if ((*(fourth_attr_obj.str_vec_))[0] != "seqMeta") {
    cout << "\nERROR in parsing SeqMeta RData object: Main Attribute's "
         << "fourth Attribute's Object should be 'seqMeta'; but found '"
         << (*(fourth_attr_obj.str_vec_))[0] << "'. Aborting." << endl;
    return false;
  }
  return true;
}

bool ParseSeqMetaRDataMainObjectAttributes(
    const RDataNode& main_attr, bool* is_skat_cohort_version, vector<string>* genes) {
  // Sanity check Main Object's Attributes does not have a Tag or Attributes.
  if (main_attr.attr_ != nullptr || main_attr.tag_ != nullptr ||
      main_attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format "
         << "for main_attr. Aborting." << endl;
    return false;
  }

  // Sanity check Main Object's Attribute type is a List.
  const RDataObject& main_attr_obj = *(main_attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, main_attr_obj) ||
      main_attr_obj.pair_list_vec_->empty()) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format "
         << "for Main Attribute's Object. Aborting." << endl;
    return false;
  }

  const bool is_old_version = main_attr_obj.pair_list_vec_->size() != 4;
  // Prompt warning about unsupported version if List length is not 4.
  if (is_old_version) {
    cout << "\nWARNING in parsing SeqMeta RData object: Expected 4 attributes "
         << "for the seqMeta object, found "
         << main_attr_obj.pair_list_vec_->size()
         << ". This likely means the version of seqMeta that was used to "
         << "create this RData Object is not a version supported by PreMeta. "
         << "PreMeta will attempt to parse this RData object as if it came "
         << "from a supported version of seqMeta." << endl;
  }
  const vector<RDataPairList*>& main_attr_items = *(main_attr_obj.pair_list_vec_);

  // Parse first Attribute: sanity check it has Tag "dim", and record num_genes.
  int num_genes = -1;
  if (!is_old_version &&
      !ParseSeqMetaRDataMainObjectAttributeOne(main_attr_items[0], &num_genes)) {
    return false;
  }

  // Parse second Attribute: sanity check it has Tag "dimnames" and store genes.
  const int gene_names_index = is_old_version ? 0 : 1;
  const string tag_name = is_old_version ? "names" : "dimnames";
  if (!ParseSeqMetaRDataMainObjectAttributeTwo(
          is_old_version, num_genes, tag_name,
          main_attr_items[gene_names_index], genes)) {
    return false;
  }

  // Parse third Attribute: sanity check it has Tag "family" and value "gaussian"
  if (!is_old_version &&
      !ParseSeqMetaRDataMainObjectAttributeThree(main_attr_items[2])) {
    return false;
  }

  // Parse fourth Attribute: sanity check it has Tag "class" and value "seqMeta"
  const int class_name_index = is_old_version ? 1 : 3;
  if (!ParseSeqMetaRDataMainObjectAttributeFour(
          main_attr_items[class_name_index], is_skat_cohort_version)) {
    return false;
  }

  return true;
}

bool ParseSeqMetaRDataGeneScores(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    int* num_snps,
    const RDataNode* node,
    const string& gene_name,
    set<Position>* snps_to_skip,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores) {
  // Sanity check we haven't seen this geen before.
  if (gene_to_snp_positions->find(gene_name) != gene_to_snp_positions->end()) {
    cout << "\nERROR in parsing SeqMeta RData object: Already have processed "
         << "gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }

  // Sanity check node has expected format.
  if (node == nullptr || node->tag_ != nullptr ||
      node->attr_ == nullptr || node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object's Gene scores: Node has bad format. "
         << "Aborting." << endl;
    return false;
  }

  // Parse Attribute.
  const RDataNode& attr = *(node->attr_);
  if (attr.tag_ != nullptr || attr.attr_ != nullptr || attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Attr has bad format. "
         << "Aborting." << endl;
    return false;
  }
  const RDataObject& attr_obj = *(attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, attr_obj) ||
      attr_obj.pair_list_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Attr is not a list "
         << "with 1 element. Aborting." << endl;
    return false;
  }
  const RDataPairList* attr_list = (*(attr_obj.pair_list_vec_))[0];
  if (attr_list == nullptr || attr_list->tag_ == nullptr ||
      attr_list->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Attribute List "
         << "has bad format. Aborting." << endl;
    return false;
  }
  if (*(attr_list->tag_) != "names") {
    cout << "\nERROR in parsing SeqMeta RData object: Attribute Tag is '"
         << *(attr_list->tag_) << "' (expected 'names'). Aborting." << endl;
    return false;
  }
  const RDataNode& attr_list_obj_node = *(attr_list->obj_);
  if (attr_list_obj_node.tag_ != nullptr ||
      attr_list_obj_node.attr_ != nullptr ||
      attr_list_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Attribute List's "
         << "object has bad format. Aborting." << endl;
    return false;
  }
  const RDataObject& attr_names = *(attr_list_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, attr_names)) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected a String "
         << "vector of attribute names. Aborting." << endl;
    return false;
  }
  const vector<string>& names = *(attr_names.str_vec_);
  *num_snps += names.size();

  // Parse Scores.
  const RDataObject& scores_obj = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_REAL, scores_obj) ||
      scores_obj.real_vec_->size() != names.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Scores Object "
         << "has bad format. Aborting." << endl;
    return false;
  }
  const vector<double>& score_values = *(scores_obj.real_vec_);

  // Store scores and names.
  vector<Position>& positions = gene_to_snp_positions->insert(
      make_pair(gene_name, vector<Position>())).first->second;
  for (int i = 0; i < names.size(); ++i) {
    const double& value = score_values[i];
    const string& name = names[i];
    Position pos;
    if (!ParsePosition(name, &pos)) {
      pos.snp_id_ = name;
      pos.chr_ = Chromosome::CHROMOSOME_UNKNOWN;
      pos.pos_ = 0;
    }
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(pos);
      if (allele_info_itr == allele_info.end()) {
        if (require_snp_in_allele_file) {
          (*num_snps)--;
          snps_to_skip->insert(pos);
          continue;
        }
      } else {
        pos = allele_info_itr->second.pos_;
      }
    } else if (require_snp_in_allele_file) {
      cout << "\nERROR in getting info from Mass file: SNP_INFO file "
           << "not specified (or empty), so Major/Minor alleles are "
           << "not available (target software requires this information. "
           << "Aborting." << endl;
      return false;
    }

    // Make sure this SNP shouldn't be excluded from this study (based
    // on mismatching REF/ALT alleles).
    if (ShouldExcludeSnp(study_num, pos, snp_to_excluding_study)) {
      (*num_snps)--;
      continue;
    }

    positions.push_back(pos);

    pair<map<Position, SnpInfo>::iterator, bool> insertion_itr =
        scores->insert(make_pair(pos, SnpInfo()));
    if (!insertion_itr.second &&
        insertion_itr.first->second.u_stat_ != value) {
      cout << "\nERROR in parsing SeqMeta RData object: Position '"
           << name << "' has different score for gene '" << gene_name
           << "than when it was previously seen (" << value
           << " now vs. " << insertion_itr.first->second.u_stat_
           << "before). Aborting." << endl;
      return false;
    }
    // If we already have SnpInfo for this position, proceed to next.
    if (!insertion_itr.second) continue;
    insertion_itr.first->second.u_stat_ = value;
  }

  return true;
}

bool ParseSeqMetaRDataGeneCovarianceOld(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode* node,
    const string& gene_name,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances) {
  // Sanity check node has expected format: Old versions store the covariance
  // matrix as an (Object, Attribute), where the Object is simply a Real Vector
  // of all of the values (Row 1, Row 2, ...), and the Attribute is a Pair-List
  // with 2 elements: The first gives the dimensions of the matrix, the second
  // gives the col/row names (i.e. the gene ids).
  if (node == nullptr || node->tag_ != nullptr ||
      node->attr_ == nullptr || node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object's Covariance: Node has bad format. "
         << "Aborting." << endl;
    return false;
  }

  // First, parse the Attribute for Matrix dimensions and Row/Column Names.
  const RDataNode& object_attr = *(node->attr_);
  if (object_attr.tag_ != nullptr || object_attr.attr_ != nullptr ||
      object_attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute for "
         << "Covariance matrix. Aborting." << endl;
    return false;
  }
  const RDataObject& attr_obj = *(object_attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, attr_obj) ||
      attr_obj.pair_list_vec_->size() != 2) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute for "
         << "Covariance matrix: Expected two attribute descriptors (for "
         << "covariance matrix dimensions and row/col names), but found "
         << attr_obj.pair_list_vec_->size() << ". Aborting." << endl;
    return false;
  }
  // Parse Covariance Matrix size (number of dimensions).
  const RDataPairList* dim_list = (*(attr_obj.pair_list_vec_))[0];
  if (dim_list == nullptr || dim_list->tag_ == nullptr ||
      dim_list->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected "
         << "Covariance matrix's 'dim' attribute to have a Tag and 2-D "
         << "Integer vector (representing the dimensions). Aborting." << endl;
    return false;
  }
  if (*(dim_list->tag_) != "dim") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected "
         << "Covariance matrix's first attribute to be the Dimensions "
         << "of the matrix; but found tag '" << *(dim_list->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& dim_list_obj_node = *(dim_list->obj_);
  if (dim_list_obj_node.tag_ != nullptr || dim_list_obj_node.attr_ != nullptr ||
      dim_list_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute for "
         << "Covariance matrix. Aborting." << endl;
    return false;
  }
  const RDataObject& dim_vector = *(dim_list_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, dim_vector) ||
      dim_vector.int_vec_->size() != 2) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix "
         << "Dimensions are not an integer vector of size 2. Aborting."
         << endl;
    return false;
  }
  const int num_rows = (*(dim_vector.int_vec_))[0];
  if (num_rows != (*(dim_vector.int_vec_))[1]) {
    cout << "\nERROR in parsing SeqMeta RData object: Non-square Matrix "
         << "has dimensions: (" << num_rows << ", "
         << (*(dim_vector.int_vec_))[1] << "). Aborting." << endl;
    return false;
  }
  // Parse Covariance Matrix row/col names.
  const RDataPairList* dimnames_list_node = (*(attr_obj.pair_list_vec_))[1];
  if (dimnames_list_node == nullptr || dimnames_list_node->tag_ == nullptr ||
      dimnames_list_node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected "
         << "Covariance matrix's 'dimnames' attribute to have a Tag and two "
         << "lists (representing the row/col names). Aborting." << endl;
    return false;
  }
  if (*(dimnames_list_node->tag_) != "dimnames") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected "
         << "Covariance matrix's 2nd attribute to be the Row/Col names "
         << "but found tag '" << *(dimnames_list_node->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& dimnames_list_obj_node = *(dimnames_list_node->obj_);
  if (dimnames_list_obj_node.tag_ != nullptr ||
      dimnames_list_obj_node.attr_ != nullptr ||
      dimnames_list_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute for "
         << "Covariance matrix. Aborting." << endl;
    return false;
  }
  const RDataObject& dimnames_list = *(dimnames_list_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_LIST, dimnames_list) ||
      dimnames_list.list_vec_->size() != 2) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix "
         << "Dimnames is not a list of size 2. Aborting."
         << endl;
    return false;
  }
  // Parse row names.
  const RDataNode* rownames_node = (*(dimnames_list.list_vec_))[0];
  if (rownames_node == nullptr || rownames_node->tag_ != nullptr ||
      rownames_node->attr_ != nullptr || rownames_node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix's "
         << "rownames_node has wrong format. Aborting." << endl;
    return false;
  }
  const RDataObject& rownames_obj = *(rownames_node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, rownames_obj) ||
      rownames_obj.str_vec_->size() != num_rows) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix's "
         << "rownames_obj is not a string vector of length " << num_rows
         << ". Aborting." << endl;
    return false;
  }
  const vector<string>& rownames_str = *(rownames_obj.str_vec_);
  // Parse column names.
  const RDataNode* colnames_node = (*(dimnames_list.list_vec_))[1];
  if (colnames_node == nullptr || colnames_node->tag_ != nullptr ||
      colnames_node->attr_ != nullptr || colnames_node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix's "
         << "colnames_node has wrong format. Aborting." << endl;
    return false;
  }
  const RDataObject& colnames_obj = *(colnames_node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, colnames_obj) ||
      colnames_obj.str_vec_->size() != num_rows) {
    cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix's "
         << "colnames_obj is not a string vector of length " << num_rows
         << ". Aborting." << endl;
    return false;
  }
  const vector<string>& colnames_str = *(colnames_obj.str_vec_);
  if (!equal(rownames_str.begin(), rownames_str.begin() + num_rows,
             colnames_str.begin())) {
    cout << "\nERROR in parsing SeqMeta RData object: Rownames don't match "
         << "colnames:\n\t[" << Join(rownames_str, ", ") << "]\n\t["
         << Join(colnames_str, ", ") << "]\nAborting." << endl;
    return false;
  }
  // Convert rownames from String to Position.
  vector<Position> rownames;
  for (const string& name : rownames_str) {
    rownames.push_back(Position());
    Position& pos = rownames.back();
    if (!ParsePosition(name, &pos)) {
      pos.snp_id_ = name;
      pos.chr_ = Chromosome::CHROMOSOME_UNKNOWN;
      pos.pos_ = 0;
    }
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(pos);
      if (allele_info_itr == allele_info.end()) {
        if (require_snp_in_allele_file) {
          snps_to_skip->insert(pos);
          continue;
        }
      } else {
        pos = allele_info_itr->second.pos_;
      }
    } else if (require_snp_in_allele_file) {
      cout << "\nERROR in getting info from Mass file: SNP_INFO file "
           << "not specified (or empty), so Major/Minor alleles are "
           << "not available (target software requires this information. "
           << "Aborting." << endl;
      return false;
    }
  }
  // Sanity check Positions for this gene are consistent to what was read
  // earlier.
  map<string, vector<Position>>::const_iterator gene_to_pos_itr =
      gene_to_snp_positions.find(gene_name);
  if (gene_to_pos_itr == gene_to_snp_positions.end()) {
    cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
         << "new gene '" << gene_name << "'. Aborting." << endl; return false;
  }
  const vector<Position>& positions = gene_to_pos_itr->second;
  if (rownames.size() != positions.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Found "
         << rownames.size() << " positions for gene '" << gene_name
         << "', but earlier only found " << positions.size()
         << " positions. Aborting." << endl;
    return false;
  }
  for (int i = 0; i < num_rows; ++i) {
    if (rownames[i].chr_ != positions[i].chr_ ||
        rownames[i].pos_ != positions[i].pos_ ||
        rownames[i].snp_id_ != positions[i].snp_id_) {
      cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
           << "new position " << PrintPosition(rownames[i])
           << " in gene '" << gene_name << "'. Aborting." << endl;
      return false;
    }
  }
  
  // Now Parse Covariance Matrix.
  const RDataObject& values_vector = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_REAL, values_vector) ||
      values_vector.real_vec_->size() != num_rows * num_rows) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected Covariance "
         << "Matrix to have size (" << num_rows << ", " << num_rows
         << "), but found " << values_vector.real_vec_->size()
         << "entries. Aborting." << endl;
    return false;
  }
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_rows; ++col) {
      const double& value = (*values_vector.real_vec_)[row * num_rows + col];
      // Covariance matrix should be symmetric. Only need to populate values
      // once (e.g. for lower-triangular part).
      if (col <= row) {
        const Position& pos_one = rownames[row];
        const Position& pos_two = rownames[col];
        if (ShouldExcludeSnp(study_num, pos_one, snp_to_excluding_study) ||
            ShouldExcludeSnp(study_num, pos_two, snp_to_excluding_study)) {
          continue;
        }
        pair<map<pair<Position, Position>, double>::iterator, bool> insertion_itr =
            covariances->insert(make_pair(make_pair(pos_one, pos_two), value));
        if (!insertion_itr.second && insertion_itr.first->second != value) {
          cout << "\nERROR in parsing SeqMeta RData object: covariance for pair ("
               << pos_one.chr_ << ":" << pos_one.pos_ << ", " << pos_two.chr_
               << ":" << pos_two.pos_ << ") was previously found to be "
               << insertion_itr.first->second << ", but in gene '" << gene_name
               << "' it was found to be " << value << ". Aborting." << endl;
          return false;
        }
      } else {
        // Sanity-check Covariance Matrix was symmetric.
        if (!FloatEq(value, (*values_vector.real_vec_)[col * num_rows + row])) {
          cout << "\nERROR in parsing SeqMeta RData object: Covariance Matrix "
               << "is not symmetric at position (" << row + 1 << ", " << col + 1
               << "): " << value << " is not equal to position ("
               << col + 1 << ", " << row + 1 << "): "
               << (*values_vector.real_vec_)[col * num_rows + row]
               << ". Aborting." << endl;
          return false;
        }
      }
    }
  }

  return true;
}

bool ParseSeqMetaRDataGeneCovariance(
    const int study_num,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode* node,
    const string& gene_name,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* snps_to_skip,
    map<pair<Position, Position>, double>* covariances) {
  // Sanity check node has expected format.
  if (node == nullptr || node->tag_ != nullptr ||
      node->attr_ != nullptr || node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object's Covariance: Node has bad format. "
         << "Aborting." << endl;
    return false;
  }
  const RDataObject& class_object = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_CLASS, class_object) ||
      class_object.class_->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected a Class Object "
         << "for the first Covariance block. Aborting." << endl;
    return false;
  }
  
  // There should be 8 elements in a list:
  //   {i, p, Dim, Dimnames, x, uplo, factors, class}.
  const RDataObject& obj = *(class_object.class_->obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, obj) ||
      obj.pair_list_vec_->size() != 8) {
    cout << "\nERROR in parsing SeqMeta RData object: Object's list does not "
         << "have the expected 8 elements. Aborting." << endl;
    return false;
  }

  // Parse 'i'.
  const RDataPairList* first_item = (*(obj.pair_list_vec_))[0];
  if (first_item == nullptr || first_item->tag_ == nullptr ||
      first_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's first item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(first_item->tag_) != "i") {
    cout << "\nERROR in parsing SeqMeta RData object: List's first item "
         << "tag is not 'i': '" << *(first_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& first_item_obj_node = *(first_item->obj_);
  if (first_item_obj_node.tag_ != nullptr ||
      first_item_obj_node.attr_ != nullptr ||
      first_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's first item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& first_obj = *(first_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, first_obj)) {
    cout << "\nERROR in parsing SeqMeta RData object: List's first item's "
         << "object is not an integer vector. Aborting."
         << endl;
    return false;
  }
  const vector<int>& i_info = *(first_obj.int_vec_);

  // Parse 'p'.
  const RDataPairList* second_item = (*(obj.pair_list_vec_))[1];
  if (second_item == nullptr || second_item->tag_ == nullptr ||
      second_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's second item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(second_item->tag_) != "p") {
    cout << "\nERROR in parsing SeqMeta RData object: List's second item "
         << "tag is not 'p': '" << *(second_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& second_item_obj_node = *(second_item->obj_);
  if (second_item_obj_node.tag_ != nullptr ||
      second_item_obj_node.attr_ != nullptr ||
      second_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's second item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& second_obj = *(second_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, second_obj)) {
    cout << "\nERROR in parsing SeqMeta RData object: List's second item's "
         << "object is not an integer vector. Aborting."
         << endl;
    return false;
  }
  const vector<int>& p_info = *(second_obj.int_vec_);
  const int num_nonzero_elements = p_info.back();

  // Parse 'Dim'.
  const RDataPairList* third_item = (*(obj.pair_list_vec_))[2];
  if (third_item == nullptr || third_item->tag_ == nullptr ||
      third_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's third item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(third_item->tag_) != "Dim") {
    cout << "\nERROR in parsing SeqMeta RData object: List's third item "
         << "tag is not 'Dim': '" << *(third_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& third_item_obj_node = *(third_item->obj_);
  if (third_item_obj_node.tag_ != nullptr ||
      third_item_obj_node.attr_ != nullptr ||
      third_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's third item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& third_obj = *(third_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, third_obj) ||
      third_obj.int_vec_->size() != 2) {
    cout << "\nERROR in parsing SeqMeta RData object: List's third item's "
         << "object is not an integer vector of size 2. Aborting."
         << endl;
    return false;
  }
  const int num_rows = (*(third_obj.int_vec_))[0];
  if (num_rows != (*(third_obj.int_vec_))[1]) {
    cout << "\nERROR in parsing SeqMeta RData object: Non-square Matrix "
         << "has dimensions: (" << num_rows << ", "
         << (*(third_obj.int_vec_))[1] << "). Aborting." << endl;
    return false;
  }
  if (p_info.size() != num_rows + 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected p to have "
         << "(num_cols + 1) elements, but num_cols is " << num_rows
         << " and size p is " << p_info.size() << ". Aborting." << endl;
    return false;
  }

  // Parse 'Dimnames'.
  const RDataPairList* fourth_item = (*(obj.pair_list_vec_))[3];
  if (fourth_item == nullptr || fourth_item->tag_ == nullptr ||
      fourth_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(fourth_item->tag_) != "Dimnames") {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item "
         << "tag is not 'Dimnames': '" << *(fourth_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& fourth_item_obj_node = *(fourth_item->obj_);
  if (fourth_item_obj_node.tag_ != nullptr ||
      fourth_item_obj_node.attr_ != nullptr ||
      fourth_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& fourth_obj = *(fourth_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_LIST, fourth_obj) ||
      fourth_obj.list_vec_->size() != 2) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "object is not a list of size 2. Aborting."
         << endl;
    return false;
  }
  // Parse row names.
  const RDataNode* rownames_node = (*(fourth_obj.list_vec_))[0];
  if (rownames_node == nullptr || rownames_node->tag_ != nullptr ||
      rownames_node->attr_ != nullptr || rownames_node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "rownames_node has wrong format. Aborting." << endl;
    return false;
  }
  const RDataObject& rownames_obj = *(rownames_node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, rownames_obj) ||
      rownames_obj.str_vec_->size() != num_rows) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "rownames_obj is not a string vector of length " << num_rows
         << ". Aborting." << endl;
    return false;
  }
  const vector<string>& rownames_str = *(rownames_obj.str_vec_);
  // Parse column names.
  const RDataNode* colnames_node = (*(fourth_obj.list_vec_))[1];
  if (colnames_node == nullptr || colnames_node->tag_ != nullptr ||
      colnames_node->attr_ != nullptr || colnames_node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "colnames_node has wrong format. Aborting." << endl;
    return false;
  }
  const RDataObject& colnames_obj = *(colnames_node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, colnames_obj) ||
      colnames_obj.str_vec_->size() != num_rows) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fourth item's "
         << "colnames_obj is not a string vector of length " << num_rows
         << ". Aborting." << endl;
    return false;
  }
  const vector<string>& colnames_str = *(colnames_obj.str_vec_);
  if (!equal(rownames_str.begin(), rownames_str.begin() + num_rows,
             colnames_str.begin())) {
    cout << "\nERROR in parsing SeqMeta RData object: Rownames don't match "
         << "colnames:\n\t[" << Join(rownames_str, ", ") << "]\n\t["
         << Join(colnames_str, ", ") << "]\nAborting." << endl;
    return false;
  }
  // Convert rownames from String to Position.
  vector<Position> rownames;
  for (const string& name : rownames_str) {
    rownames.push_back(Position());
    Position& pos = rownames.back();
    if (!ParsePosition(name, &pos)) {
      pos.snp_id_ = name;
      pos.chr_ = Chromosome::CHROMOSOME_UNKNOWN;
      pos.pos_ = 0;
    }
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(pos);
      if (allele_info_itr == allele_info.end()) {
        if (require_snp_in_allele_file) {
          snps_to_skip->insert(pos);
          continue;
        }
      } else {
        pos = allele_info_itr->second.pos_;
      }
    } else if (require_snp_in_allele_file) {
      cout << "\nERROR in getting info from Mass file: SNP_INFO file "
           << "not specified (or empty), so Major/Minor alleles are "
           << "not available (target software requires this information. "
           << "Aborting." << endl;
      return false;
    }
  }
  // Sanity check Positions for this gene are consistent to what was read
  // earlier.
  map<string, vector<Position>>::const_iterator gene_to_pos_itr =
      gene_to_snp_positions.find(gene_name);
  if (gene_to_pos_itr == gene_to_snp_positions.end()) {
    cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
         << "new gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }
  const vector<Position>& positions = gene_to_pos_itr->second;
  if (rownames.size() != positions.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Found "
         << rownames.size() << " positions for gene '" << gene_name
         << "', but earlier only found " << positions.size()
         << " positions. Aborting." << endl;
    return false;
  }
  for (int i = 0; i < num_rows; ++i) {
    if (rownames[i].chr_ != positions[i].chr_ ||
        rownames[i].pos_ != positions[i].pos_ ||
        rownames[i].snp_id_ != positions[i].snp_id_) {
      cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
           << "new position " << PrintPosition(rownames[i])
           << " in gene '" << gene_name << "'. Aborting." << endl;
      return false;
    }
  }

  // Parse 'x'.
  const RDataPairList* fifth_item = (*(obj.pair_list_vec_))[4];
  if (fifth_item == nullptr || fifth_item->tag_ == nullptr ||
      fifth_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fifth item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(fifth_item->tag_) != "x") {
    cout << "\nERROR in parsing SeqMeta RData object: List's fifth item "
         << "tag is not 'x': '" << *(fifth_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& fifth_item_obj_node = *(fifth_item->obj_);
  if (fifth_item_obj_node.tag_ != nullptr ||
      fifth_item_obj_node.attr_ != nullptr ||
      fifth_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fifth item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& fifth_obj = *(fifth_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_REAL, fifth_obj) ||
      fifth_obj.real_vec_->size() != num_nonzero_elements) {
    cout << "\nERROR in parsing SeqMeta RData object: List's fifth item's "
         << "object is not a real vector of size " << num_nonzero_elements
         << ". Aborting." << endl;
    return false;
  }
  const vector<double>& values = *(fifth_obj.real_vec_);

  // Parse 'uplo'.
  const RDataPairList* sixth_item = (*(obj.pair_list_vec_))[5];
  if (sixth_item == nullptr || sixth_item->tag_ == nullptr ||
      sixth_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's sixth item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(sixth_item->tag_) != "uplo") {
    cout << "\nERROR in parsing SeqMeta RData object: List's sixth item "
         << "tag is not 'uplo': '" << *(sixth_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& sixth_item_obj_node = *(sixth_item->obj_);
  if (sixth_item_obj_node.tag_ != nullptr ||
      sixth_item_obj_node.attr_ != nullptr ||
      sixth_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's sixth item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& sixth_obj = *(sixth_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, sixth_obj) ||
      sixth_obj.str_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: List's sixth item's "
         << "object is not a String vector of size 1. Aborting."
         << endl;
    return false;
  }
  const string& uplo_type = (*(sixth_obj.str_vec_))[0];
  if (uplo_type != "U") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected uplo type "
         << "'U', but found '" << uplo_type << "'. Aborting." << endl;
    return false;
  }

  // Parse 'factors'.
  const RDataPairList* seventh_item = (*(obj.pair_list_vec_))[6];
  if (seventh_item == nullptr || seventh_item->tag_ == nullptr ||
      seventh_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's seventh item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(seventh_item->tag_) != "factors") {
    cout << "\nERROR in parsing SeqMeta RData object: List's seventh item "
         << "tag is not 'uplo': '" << *(seventh_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  // We don't care about factors, so we won't parse them.

  // Parse 'class'.
  const RDataPairList* eighth_item = (*(obj.pair_list_vec_))[7];
  if (eighth_item == nullptr || eighth_item->tag_ == nullptr ||
      eighth_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's eighth item "
         << "does not have expected format. Aborting." << endl;
    return false;
  }
  if (*(eighth_item->tag_) != "class") {
    cout << "\nERROR in parsing SeqMeta RData object: List's eighth item "
         << "tag is not 'class': '" << *(eighth_item->tag_) << "'. Aborting."
         << endl;
    return false;
  }
  const RDataNode& eighth_item_obj_node = *(eighth_item->obj_);
  if (eighth_item_obj_node.tag_ != nullptr ||
      eighth_item_obj_node.attr_ == nullptr ||
      eighth_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's eighth item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& eighth_obj = *(eighth_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, eighth_obj) ||
      eighth_obj.str_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: List's eighth item's "
         << "object is not a String vector of size 1. Aborting."
         << endl;
    return false;
  }
  const string& class_name = (*(eighth_obj.str_vec_))[0];
  if (class_name != "dsCMatrix") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class name "
         << "'dsCMatrix', but found '" << class_name << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& class_attr = *(eighth_item_obj_node.attr_);
  if (class_attr.tag_ != nullptr || class_attr.attr_ != nullptr ||
      class_attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad format for class "
         << "attribute. Aborting." << endl;
    return false;
  }
  const RDataObject& class_obj = *(class_attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, class_obj) ||
      class_obj.pair_list_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class' object "
         << "to be a list of size 1. Aborting." << endl;
    return false;
  }
  const RDataPairList* class_list = (*(class_obj.pair_list_vec_))[0];
  if (class_list == nullptr || class_list->tag_ == nullptr ||
      class_list->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class' list "
         << "to have a Tag and string vector. Aborting." << endl;
    return false;
  }
  if (*(class_list->tag_) != "package") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class' list "
         << "to have tag 'package', but found '" << *(class_list->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& class_list_obj_node = *(class_list->obj_);
  if (class_list_obj_node.tag_ != nullptr ||
      class_list_obj_node.attr_ != nullptr ||
      class_list_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: List's first item's "
         << "object node does not have expected format. Aborting." << endl;
    return false;
  }
  const RDataObject& class_package_obj = *(class_list_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, class_package_obj) ||
      class_package_obj.str_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class' "
         << "package to be a string vector with one element. Aborting." << endl;
    return false;
  }
  if ((*(class_package_obj.str_vec_))[0] != "Matrix") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected class' "
         << "package to be 'Matrix', but found '"
         << (*(class_package_obj.str_vec_))[0] << "'. Aborting." << endl;
    return false;
  }

  // Now actually parse the covariance matrix.
  int i_itr = 0;
  int p_itr = 0;
  int current_col_num_nonzeros_total = 0;
  int current_col_num_nonzeros_used = 0;
  int prev_row = 0;
  int prev_col = 0;
  for (const double& value : values) {
    while ((current_col_num_nonzeros_total == current_col_num_nonzeros_used) &&
           (num_rows >= p_itr + 1)) {
      current_col_num_nonzeros_total = p_info[p_itr + 1] - p_info[p_itr];
      current_col_num_nonzeros_used = 0;
      p_itr++;
    }
    if (current_col_num_nonzeros_total == current_col_num_nonzeros_used) {
      cout << "\nERROR in parsing SeqMeta RData object: Did not find all "
           << "non-zero values before reaching the end of the cov matrix. "
           << "num_rows: " << num_rows << ", i_itr: " << i_itr << ", p_itr: "
           << p_itr << ", current_col_num_nonzeros_total: "
           << current_col_num_nonzeros_total
           << ", current_col_num_nonzeros_used: "
           << current_col_num_nonzeros_used
           << "Aborting." << endl;
      return false;
    }
    const int row = i_info[i_itr];
    const int col = p_itr - 1;
    // Fill in all zero-variances between last non-zero value and this one.
    for (int c = prev_col; c <= col; ++c) {
      if (ShouldExcludeSnp(study_num, rownames[c], snp_to_excluding_study)) {
        continue;
      }
      for (int r = (c == prev_col ? (prev_row + 1) : 0);
           (r <= c && (c < col || r < row)); ++r) {
        if (ShouldExcludeSnp(study_num, rownames[r], snp_to_excluding_study)) {
          continue;
        }
        pair<map<pair<Position, Position>, double>::iterator, bool> insertion_itr =
            covariances->insert(make_pair(make_pair(rownames[r], rownames[c]), 0.0));
        if (!insertion_itr.second && insertion_itr.first->second != 0.0) {
          cout << "\nERROR in parsing SeqMeta RData object: covariance for pair ("
               << rownames[r].chr_ << ":" << rownames[r].pos_ << ", "
               << rownames[c].chr_ << ":" << rownames[c].pos_
               << ") was previously found to be "
               << insertion_itr.first->second << ", but in gene '" << gene_name
               << "' it was found to be 0. Aborting." << endl;
          return false;
        }
      }
    }
    const Position& pos_one = (row < col) ? rownames[row] : rownames[col];
    const Position& pos_two = (row < col) ? rownames[col] : rownames[row];
    if (ShouldExcludeSnp(study_num, pos_one, snp_to_excluding_study) ||
        ShouldExcludeSnp(study_num, pos_two, snp_to_excluding_study)) {
      continue;
    }
    pair<map<pair<Position, Position>, double>::iterator, bool> insertion_itr =
        covariances->insert(make_pair(make_pair(pos_one, pos_two), value));
    if (!insertion_itr.second && insertion_itr.first->second != value) {
      cout << "\nERROR in parsing SeqMeta RData object: covariance for pair ("
           << pos_one.chr_ << ":" << pos_one.pos_ << ", " << pos_two.chr_
           << ":" << pos_two.pos_ << ") was previously found to be "
           << insertion_itr.first->second << ", but in gene '" << gene_name
           << "' it was found to be " << value << ". Aborting." << endl;
      return false;
    }
    prev_row = row;
    prev_col = col;
    current_col_num_nonzeros_used++;
    i_itr++;
  }

  return true;
}

bool ParseSeqMetaRDataGeneN(const RDataNode* node, int* num_samples) {
  if (node == nullptr || node->tag_ != nullptr || node->attr_ != nullptr ||
      node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: node has bad format. Aborting.\n";
    return false;
  }
  const RDataObject& obj = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_INT, obj) ||
      obj.int_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: bad object. Aborting." << endl;
    return false;
  }
  const int n = (*(obj.int_vec_))[0];
  if (*num_samples == -1) {
    *num_samples = n;
  }
  if (*num_samples != n) {
    cout << "\nERROR in parsing SeqMeta RData object: Previous value for n ("
         << *num_samples << ") does not match current value (" << n
         << "). Aborting." << endl;
    return false;
  }

  return true;
}

bool ParseSeqMetaRDataGeneMaf(
    const int study_num,
    const int num_samples,
    const RDataNode* node,
    const string& gene_name,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    set<Position>* monomorphic_snps, map<Position, SnpInfo>* scores) {
  // Sanity check Positions for this gene are consistent to what was read
  // earlier.
  map<string, vector<Position>>::const_iterator gene_to_pos_itr =
      gene_to_snp_positions.find(gene_name);
  if (gene_to_pos_itr == gene_to_snp_positions.end()) {
    cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
         << "new gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }
  const vector<Position>& positions = gene_to_pos_itr->second;

  // Sanity Check node's format.
  if (node == nullptr || node->tag_ != nullptr ||
      node->attr_ == nullptr || node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: bad node format."
         << "Aborting." << endl;
    return false;
  }

  // Parse Attribute 'names': all the Positions associated to this gene.
  const RDataNode& attr = *(node->attr_);
  if (attr.tag_ != nullptr || attr.attr_ != nullptr || attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute format. "
         << "Aborting." << endl;
    return false;
  }
  const RDataObject& attr_obj = *(attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, attr_obj) ||
      attr_obj.pair_list_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: Attribute should be a list "
         << "with one object. Aborting." << endl;
    return false;
  }
  const RDataPairList* attr_list = (*(attr_obj.pair_list_vec_))[0];
  if (attr_list == nullptr || attr_list->tag_ == nullptr ||
      attr_list->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute list format. "
         << "Aborting." << endl;
    return false;
  }
  if (*(attr_list->tag_) != "names") {
    cout << "\nERROR in parsing SeqMeta RData object: Expected Attribute's Tag "
         << "to be 'names', but found '" << *(attr_list->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& attr_list_obj_node = *(attr_list->obj_);
  if (attr_list_obj_node.tag_ != nullptr ||
      attr_list_obj_node.attr_ != nullptr ||
      attr_list_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Bad Attribute format in "
         << "list's object. Aborting." << endl;
    return false;
  }
  const RDataObject& names_obj = *(attr_list_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, names_obj) ||
      names_obj.str_vec_->size() != positions.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected Attribute "
         << "'names' to be a vector with " << positions.size()
         << "elements. Aborting." << endl;
    return false;
  }
  const vector<string>& names = *(names_obj.str_vec_);

  // Parse Object (Maf values).
  const RDataObject& maf_obj = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_REAL, maf_obj) ||
      maf_obj.real_vec_->size() != names.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Expected MAF object to be a "
         << "real vector of length " << names.size() << ". Aborting." << endl;
    return false;
  }
  const vector<double> maf = *(maf_obj.real_vec_);
  for (int i = 0; i < names.size(); ++i) {
    Position pos;
    if (!ParsePosition(names[i], &pos)) {
      // Unable to parse Position into CHR:POS format.
      // NOTE(11/23/15): Removed the following 3 lines and added the 4th-6th
      // lines, as we now support Positions that are not in CHR:POS format.
      //cout << "\nERROR in parsing SeqMeta RData object: Unable to parse name '"
      //     << names[i] << "' as a Position. Aborting." << endl;
      //return false;
      pos.snp_id_ = names[i];
      pos.chr_ = Chromosome::CHROMOSOME_UNKNOWN;
      pos.pos_ = 0;
    }
    if (ShouldExcludeSnp(study_num, pos, snp_to_excluding_study)) {
      continue;
    }
    // Use allele_info to get format of this SNP in the target software.
    if (!allele_info.empty()) {
      map<Position, SnpInfo>::const_iterator allele_info_itr =
          allele_info.find(pos);
      if (allele_info_itr != allele_info.end()) {
        pos = allele_info_itr->second.pos_;
      }
    }
    if (positions[i].chr_ != pos.chr_ || positions[i].pos_ != pos.pos_ ||
        positions[i].snp_id_ != pos.snp_id_) {
      cout << "\nERROR in parsing SeqMeta RData object: Found unexpected "
           << "new position " << pos.chr_ << ":" << pos.pos_
           << " in gene '" << gene_name << "'. Aborting." << endl;
      return false;
    }
    map<Position, SnpInfo>::iterator snp_itr = scores->find(pos);
    if (snp_itr == scores->end()) {
      cout << "\nERROR in parsing SeqMeta RData object: No data (score) found "
           << "for Position " << pos.chr_ << ":" << pos.pos_
           << ". Aborting." << endl;
      return false;
    }
    SnpInfo& info = snp_itr->second;
    if (FloatEq(info.maf_, -1.0)) {
      info.maf_ = maf[i];
      if (num_samples == 0 ||
          info.maf_ < 1.0 / (4.0 * num_samples) ||
          info.maf_ > 1.0 - (1.0 / (4.0 * num_samples))) {
        monomorphic_snps->insert(pos);
      }
    }
    if (!FloatEq(info.maf_, maf[i])) {
      cout << "\nERROR in parsing SeqMeta RData object: MAF for Position: "
           << names[i] << " does not match previous value ("
           << info.maf_ << " vs. " << maf[i] << "). Aborting." << endl;
      return false;
    }
  }

  return true;
}

bool ParseSeqMetaRDataGeneSey(const RDataNode* node, double* sey) {
  if (node == nullptr || node->tag_ != nullptr || node->attr_ != nullptr ||
      node->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: node has bad format. Aborting.\n";
    return false;
  }
  const RDataObject& obj = *(node->obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_REAL, obj) ||
      obj.real_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: bad object. Aborting." << endl;
    return false;
  }
  const double& new_sey = (*(obj.real_vec_))[0];

  if (FloatEq(*sey, -1.0)) {
    *sey = new_sey;
  }
  if (!FloatEq(*sey, new_sey)) {
    cout << "\nERROR in parsing SeqMeta RData object: Previous value for sey ("
         << *sey << ") does not match current value (" << new_sey
         << "). Aborting." << endl;
    return false;
  }
  return true;
}

bool ParseSeqMetaRDataGene(
    const int study_num,
    const bool is_skat_cohort_version,
    const bool require_snp_in_allele_file,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const string& gene_name,
    const RDataNode* item, int* num_samples, int* num_snps, double* sey,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores,
    map<pair<Position, Position>, double>* covariances) {
  if (item == nullptr || item->tag_ != nullptr ||
      item->attr_ == nullptr || item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: item has unexpected format. "
         << "Aborting." << endl;
    return false;
  }

  // Sanity Check Attribute: Tag 'names' w/ Obj: [scores, cov, n, maf, sey]
  const RDataNode& item_attr = *(item->attr_);
  if (item_attr.tag_ != nullptr || item_attr.attr_ != nullptr ||
      item_attr.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: item has unexpected Attribute. "
         << "Aborting." << endl;
    return false;
  }

  // Parse Attributes.
  const RDataObject& item_attr_obj = *(item_attr.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, item_attr_obj) ||
      item_attr_obj.pair_list_vec_->size() != 1) {
    cout << "\nERROR in parsing SeqMeta RData object: item has unexpected Attribute "
         << "Object. Aborting." << endl;
    return false;
  }
  const RDataPairList* item_attr_list_item = (*(item_attr_obj.pair_list_vec_))[0];
  if (item_attr_list_item == nullptr || item_attr_list_item->tag_ == nullptr ||
      item_attr_list_item->obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: item has unexpected Attribute "
         << "Object format. Aborting." << endl;
    return false;
  }
  if (*(item_attr_list_item->tag_) != "names") {
    cout << "\nERROR in parsing SeqMeta RData object: item's Attribute tag is unexpected. "
         << "Expected 'names', found '" << *(item_attr_list_item->tag_)
         << "'. Aborting." << endl;
    return false;
  }
  const RDataNode& item_attr_list_item_obj_node = *(item_attr_list_item->obj_);
  if (item_attr_list_item_obj_node.tag_ != nullptr ||
      item_attr_list_item_obj_node.attr_ != nullptr ||
      item_attr_list_item_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: item has unexpected Attribute "
         << "Object Node format. Aborting." << endl;
    return false;
  }
  const RDataObject& attr_vec = *(item_attr_list_item_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_VECTOR_STRING, attr_vec) ||
      attr_vec.str_vec_->size() != 5) {
    cout << "\nERROR in parsing SeqMeta RData object: item's Attribute vector has wrong "
         << "size. Aborting." << endl;
    return false;
  }
  if ((*(attr_vec.str_vec_))[0] != "scores" ||
      (*(attr_vec.str_vec_))[1] != "cov" ||
      (*(attr_vec.str_vec_))[2] != "n" ||
      (*(attr_vec.str_vec_))[3] != "maf" ||
      (*(attr_vec.str_vec_))[4] != "sey") {
    cout << "\nERROR in parsing SeqMeta RData object: item's Attribute vector has "
         << "unexpected names: ["
         << Join(*(attr_vec.str_vec_), ", ") << "]. Aborting.\n";
    return false;
  }

  // Parse Item: It is a list of 5 items, corresponding to
  // [scores, cov, n, maf, sey].
  const RDataObject& item_obj = *(item->obj_);
  if (!ExpectObjectType(SEXPTYPE_LIST, item_obj) ||
      item_obj.list_vec_->size() != 5) {
    cout << "\nERROR in parsing SeqMeta RData object: item's Object has wrong size. "
         << "Expected 5 dimension ([scores, cov, n, maf, sey]), but found "
         << item_obj.list_vec_->size() << ". Aborting." << endl;
    return false;
  }

  // Parse scores.
  if (!ParseSeqMetaRDataGeneScores(
          study_num, require_snp_in_allele_file, allele_info, snp_to_excluding_study,
          num_snps, (*(item_obj.list_vec_))[0], gene_name,
          snps_to_skip, gene_to_snp_positions, scores)) {
    cout << "\nERROR in parsing SeqMeta RData object::ParseSeqMetaRDataGeneScores "
         << "for gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }

  // Parse cov.
  if (is_skat_cohort_version && !ParseSeqMetaRDataGeneCovarianceOld(
          study_num, require_snp_in_allele_file, allele_info,
          snp_to_excluding_study, (*(item_obj.list_vec_))[1],
          gene_name, *gene_to_snp_positions, snps_to_skip, covariances)) {
    cout << "\nERROR in parsing (Old version of) SeqMeta RData object: "
         << "Could not get covariance for gene '" << gene_name
         << "'. Aborting." << endl;
    return false;
  } else if (!is_skat_cohort_version && !ParseSeqMetaRDataGeneCovariance(
          study_num, require_snp_in_allele_file, allele_info,
          snp_to_excluding_study, (*(item_obj.list_vec_))[1],
          gene_name, *gene_to_snp_positions, snps_to_skip, covariances)) {
    cout << "\nERROR in parsing SeqMeta RData object: "
         << "Could not get covariance for gene '" << gene_name
         << "'. Aborting." << endl;
    return false;
  }

  // Parse n.
  if (!ParseSeqMetaRDataGeneN((*(item_obj.list_vec_))[2], num_samples)) {
    cout << "\nERROR in parsing SeqMeta RData object::ParseSeqMetaRDataGeneN "
         << "for gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }

  // Parse maf.
  if (!ParseSeqMetaRDataGeneMaf(
          study_num, *num_samples, (*(item_obj.list_vec_))[3], gene_name,
          allele_info, snp_to_excluding_study, *gene_to_snp_positions,
          monomorphic_snps, scores)) {
    cout << "\nERROR in parsing SeqMeta RData object::ParseSeqMetaRDataGeneMaf "
         << "for gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }

  // Parse sey.
  if (!ParseSeqMetaRDataGeneSey((*(item_obj.list_vec_))[4], sey)) {
    cout << "\nERROR in parsing SeqMeta RData object::ParseSeqMetaRDataGeneSey "
         << "for gene '" << gene_name << "'. Aborting." << endl;
    return false;
  }

  return true;
}

bool ParseSeqMetaRDataTree(
    const bool require_snp_in_allele_file, const int study_num,
    const map<Position, SnpInfo>& allele_info,
    const map<Position, set<int>>& snp_to_excluding_study,
    const RDataNode& root,
    int* num_samples, int* num_snps, double* sey,
    set<Position>* snps_to_skip, set<Position>* monomorphic_snps,
    map<string, vector<Position>>* gene_to_snp_positions,
    map<Position, SnpInfo>* scores,
    map<pair<Position, Position>, double>* covariances) {
  // Sanity check root.
  if (root.attr_ != nullptr || root.tag_ != nullptr || root.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format for root. "
         << "Aborting." << endl;
    return false;
  }

  // Sanity check root's object type is a pair list, with one element.
  const RDataObject& root_obj = *(root.obj_);
  if (!ExpectObjectType(SEXPTYPE_PAIR_LIST, root_obj) ||
      root_obj.pair_list_vec_->size() != 1 ||
      (*(root_obj.pair_list_vec_))[0] == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format for root's "
         << "Object. Aborting." << endl;
    return false;
  }

  const RDataPairList& main_node = *((*(root_obj.pair_list_vec_))[0]);

  // Sanity check main node.
  if (main_node.tag_ == nullptr || main_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format for main node. "
         << "Aborting." << endl;
    return false;
  }

  // Parse Main Object Name.
  const string main_obj_name = *(main_node.tag_);

  // Parse Main Object's Object.
  const RDataNode& main_node_obj_node = *(main_node.obj_);
  if (main_node_obj_node.tag_ != nullptr ||
      main_node_obj_node.attr_ == nullptr ||
      main_node_obj_node.obj_ == nullptr) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format for main node's "
         << "object. Aborting." << endl;
    return false;
  }

  // Parse Main Object's Attributes.
  const RDataNode& main_attr = *(main_node_obj_node.attr_);
  vector<string> genes;
  bool is_skat_cohort_version = false;
  if (!ParseSeqMetaRDataMainObjectAttributes(
          main_attr, &is_skat_cohort_version, &genes)) {
    return false;
  }

  // Parse Main Object's Items.
  const RDataObject& main_obj = *(main_node_obj_node.obj_);
  if (!ExpectObjectType(SEXPTYPE_LIST, main_obj) ||
      main_obj.list_vec_->empty()) {
    cout << "\nERROR in parsing SeqMeta RData object: Unexpected format for Main "
         << "Object. Aborting." << endl;
    return false;
  }
  const vector<RDataNode*>& items = *(main_obj.list_vec_);
  if (items.size() != genes.size()) {
    cout << "\nERROR in parsing SeqMeta RData object: Mismatching number of genes. "
         << "Expected " << genes.size() << ", but found " << items.size()
         << ". Aborting." << endl;
    return false;
  }
  for (int i = 0; i < items.size(); ++i) {
    if (!ParseSeqMetaRDataGene(
            study_num, is_skat_cohort_version, require_snp_in_allele_file,
            allele_info, snp_to_excluding_study,
            genes[i], items[i], num_samples, num_snps, sey,
            snps_to_skip, monomorphic_snps, gene_to_snp_positions,
            scores, covariances)) {
      cout << "\nERROR in parsing SeqMeta RData object for child " << i << endl;
      return false;
    }
  }

  return true;
}

}  // namespace premeta
