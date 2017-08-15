// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Structs and Enums used for PreMeta.

#ifndef PREMETA_STRUCTURES_H
#define PREMETA_STRUCTURES_H

#include "FileReaderUtils/vcf_utils.h"
#include "StringUtils/string_utils.h"

#include <cfloat>   // For DBL_MIN

using file_reader_utils::Chromosome;

namespace premeta {

// R has various String types: Byte, Latin-1, UTF-8, and ASCII.
// These are coded by the "levels" part of the pack flags, which has 16 bits, but
// only bits 1, 2, 3, 5, and 6 are relevant for determining String type (where bit
// indexing starts at '0', not '1'; i.e. relevant bits start at SECOND bit).
enum RStringType {
  R_STRING_NONE = 0,  // Level: 00000000 00000000 (shouldn't ever see this)
  R_STRING_NA,        // Level: 00000000 00000010
  R_STRING_BYTE,      // Level: 00000000 00000100
  R_STRING_LATIN,     // Level: 00000000 00001000
  R_STRING_UTF8,      // Level: 00000000 00010000
                      // Level: 00000000 00100000 (bit 5 doesn't determine type)
  R_STRING_ASCII,     // Level: 00000000 01000000
};

// R SEXPTYPES.
// See http://cran.r-project.org/doc/manuals/r-release/R-ints.html#SEXPTYPEs
// It's crucial that the order (and hence the actual enum values) are in
// correspondence with the names in the above link (in order to not have
// this dependence, I could write a function to translate
// (R SEXTYPE) int -> (premeta_utils::Sexptype) enum,
// but this seems not worth it, since (1) R SEXPTYPES have been stable for
// a number of years (20 or so), so the assumption they'll match seems ok;
// and (2) I'd have to update the function anyway; so might as well just
// update the order of the enum below and use that as the "function".
enum Sexptype {
  SEXPTYPE_UNKNOWN = 0,
  SEXPTYPE_SYMBOL,
  SEXPTYPE_PAIR_LIST,
  SEXPTYPE_CLOSURE,
  SEXPTYPE_ENVIRONMENT,
  SEXPTYPE_PROMISE,
  SEXPTYPE_LANG_OBJ,
  SEXPTYPE_SPECIAL,
  SEXPTYPE_BUILTIN,
  SEXPTYPE_STRING,
  SEXPTYPE_VECTOR_BOOL,
  SEXPTYPE_UNUSED_ELEVEN,
  SEXPTYPE_UNUSED_TWELVE,
  SEXPTYPE_VECTOR_INT,
  SEXPTYPE_VECTOR_REAL,
  SEXPTYPE_VECTOR_COMPLEX,
  SEXPTYPE_VECTOR_STRING,
  SEXPTYPE_DOTS,
  SEXPTYPE_ANY,
  SEXPTYPE_LIST,
  SEXPTYPE_VECTOR_EXPRESSION,
  SEXPTYPE_BYTE_CODE,
  SEXPTYPE_POINTER,
  SEXPTYPE_WEAK,
  SEXPTYPE_VECTOR_RAW,
  SEXPTYPE_CLASS,
  SEXPTYPE_PSUEDO,  // For Sexptype values 241-255, see:
                    // http://cran.r-project.org/doc/manuals/r-release/R-ints.html#FOOT7
};

// The versions of each software that PreMeta is equipped to handle. This enum
// if used primarily (exclusively?) in PrintVersion() below.
enum DefaultSoftwareVersion {
  SOFTWARE_VERSION_UNKNOWN,
  SOFTWARE_VERSION_MASS,
  SOFTWARE_VERSION_METASKAT,
  SOFTWARE_VERSION_RAREMETAL_NEW,
  SOFTWARE_VERSION_RAREMETAL_OLD,
  SOFTWARE_VERSION_SEQMETA,
};

// The Chromosome and Position index of a SNP; and/or the name of the SNP.
struct Position {
  Chromosome chr_;
  uint64_t pos_;
  string snp_id_;

  Position() {
    chr_ = Chromosome::CHROMOSOME_UNKNOWN;
    pos_ = 0;
    snp_id_ = "";
  }

  bool operator <(const Position& x) const {
    // Return false (indicating equality) if snp_id_ matches, irregardless
    // of Chromosome/Position.
    if (!x.snp_id_.empty() && x.snp_id_ == snp_id_) return false;

    // Handle case x's Chromosome/Position fields are not set.
    if (x.chr_ == Chromosome::CHROMOSOME_UNKNOWN && x.pos_ == 0) {
      // If this's Chromosome/Position fields are also not set, fallback
      // to snp_id_.
      if (chr_ == Chromosome::CHROMOSOME_UNKNOWN && pos_ == 0) {
        return snp_id_ < x.snp_id_;
      }
      // 'This' has Chromosome/Position set, but x doesn't. Go ahead and
      // consider x < this.
      return false;
    // 'This' doesn't have Chromosome/Position set, but x does. Go ahead and
    // consider this < x.
    } else if (chr_ == Chromosome::CHROMOSOME_UNKNOWN && pos_ == 0) {
      return true;
    }

    // Both x and 'this' have Chromosome/Position set; take whichever is
    // smaller according to these.
    return std::tie(chr_, pos_) < std::tie(x.chr_, x.pos_);
  }
};

// Structure for storing information about a given SNP.
// Relationships between some of the following:
//   - num_non_missing_ = Num_Samples * (1 - MissingRate)
//   - mac_ = num_non_missing_ * 2 * maf_
//   - mac_ = n_het_ + (2 * n_alt_)
//   - num_non_missing_ = n_ref_ + n_het_ + n_alt_
// Additional equalities for RAREMETAL terminology:
//   - CALL_RATE = num_non_missing_ / Num_Samples
//   - FOUNDER_AF = ALL_AF = maf_
//   - N_INFO = Num_Samples
//   - INFO_ALT_AC = mac_
//
struct SnpInfo {
  Position pos_;
  int n_info_;
  int num_non_missing_;
  int mac_;
  int n_ref_;
  int n_het_;
  int n_alt_;
  double missing_rate_;
  double maf_;
  double u_stat_;
  double sqrt_v_stat_;
  double var_;
  string major_allele_;
  string minor_allele_;

  SnpInfo() {
    n_info_ = -1;
    num_non_missing_ = -1;
    mac_ = -1;
    n_ref_ = -1;
    n_het_ = -1;
    n_alt_ = -1;
    maf_ = -1.0;
    missing_rate_ = -1.0;
    u_stat_ = -1.0;
    sqrt_v_stat_ = -1.0;
    var_ = -1.0;
    major_allele_ = "A";
    minor_allele_ = "T";
  }
};

// Describes the type of R object being described; See:
// http://cran.r-project.org/doc/manuals/r-release/R-ints.html
struct SexptypeInfo {
  Sexptype type_;
  bool is_obj_;
  bool has_attr_;
  bool has_tag_;
  int level_;

  SexptypeInfo() {
    type_ = SEXPTYPE_UNKNOWN;
    is_obj_ = false;
    has_attr_ = false;
    has_tag_ = false;
    level_ = 0;
  }
};

// Basic properties of a file: File name (location), comment character,
// and delimiter.
struct FileInfo {
  string name_;
  string delimiter_;
  string comment_char_;

  FileInfo() {
    name_ = "";
    delimiter_ = "\t";
    comment_char_ = "#";
  }
};

// Indicates which files are available to get allele information from.
// Since two of the four software's output summary statistics contain
// allele information, we can use those outputs to get allele info. We
// also allow the user to supply supplementary file containing allele
// information (in this case, the supplementary file MUST have form:
//   CHR:POS  REF  ALT
struct AlleleInfo {
  vector<FileInfo> supp_files_;
  vector<FileInfo> metaskat_files_;
  vector<FileInfo> raremetal_files_;
};

// Holds Version information.
struct SoftwareVersion {
  int base_;
  int sub_version_;
  int sub_sub_version_;

  SoftwareVersion() {
    base_ = -1;
    sub_version_ = -1;
    sub_sub_version_ = -1;
  }
};

// Forward declaration of the Object field of a RDataNode (definition of this
// struct is below; forward-declared here due to the cyclic dependency.
struct RDataObject;

// Holds the three features of an RData node: Object, Tag, and Attribute.
struct RDataNode {
  RDataObject* obj_;
  RDataNode* attr_;
  string* tag_;
  bool is_object_;

  RDataNode() {
    obj_ = nullptr;
    attr_ = nullptr;
    tag_ = nullptr;
    is_object_ = false;
  }
};

// Holds features (Tag and Object) for an RData Pair-List.
struct RDataPairList {
  string* tag_;
  RDataNode* obj_;

  RDataPairList() {
    tag_ = nullptr;
    obj_ = nullptr;
  }
};

// The Object part of an RData node; exactly one of the fields should be
// non-null.
struct RDataObject {
  string* str_;
  vector<int>* int_vec_;
  vector<bool>* bool_vec_;
  vector<string>* str_vec_;
  vector<double>* real_vec_;
  vector<RDataNode*>* list_vec_;
  vector<RDataPairList*>* pair_list_vec_;
  RDataNode* class_;

  RDataObject() {
    str_ = nullptr;
    int_vec_ = nullptr;
    bool_vec_ = nullptr;
    str_vec_ = nullptr;
    real_vec_ = nullptr;
    list_vec_ = nullptr;
    pair_list_vec_ = nullptr;
    class_ = nullptr;
  }
};

bool ParsePosition(const string& input, Chromosome* chr, uint64_t* pos);

bool ParsePosition(const string& input, Position* pos);

string PrintPosition(const Position& pos);

bool ParseSoftwareVersion(const string& input, SoftwareVersion* version);

string PrintSoftwareVersion(const DefaultSoftwareVersion version);
string PrintSoftwareVersion(const SoftwareVersion& version);

// Finds the column index of each string in column_indices, by searching
// through comment lines in the input file. Returns true if all column
// indices were found.
bool GetColumnIndices(
    const FileInfo& file_info, map<string, int>* column_indices);
 
// Searches covariances for (pos_one, pos_two) and (pos_two, pos_one).
bool LookupCovariance(
    const Position& pos_one, const Position& pos_two,
    const map<pair<Position, Position>, double>& covariances,
    double* covariance);

bool ShouldExcludeSnp(
    const int study_num, const Position& pos,
    const map<Position, set<int>>& snp_to_excluding_study);

bool IsSnpRefAltSwapped(
    const int study_num, const Position& pos,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study);

// The Key of allele_info is the Position of a SNP, and the corresponding
// Value (of type SnpInfo) only has 3 fields present: pos_, major_allele_,
// and minor_allele_. The pos_ field is in-case the input file uses a
// different format for SNP ID/Position than the target software.
// Then the snp_info of Key pos_ will copy major_allele_ and minor_allele_
// to the corresponding Value. Also, Positions in snp_info that are not
// present in allele_info are added to snps_to_skip.
bool CopyAlleleInfoToSnpInfo(
    const map<Position, SnpInfo>& allele_info,
    map<Position, SnpInfo>* snp_info);

}  // namespace premeta
#endif
