// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "read_r_binary_utils.h"

#include "FileReaderUtils/csv_utils.h"
#include "FileReaderUtils/vcf_utils.h"
#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <zlib.h>

using namespace std;
using file_reader_utils::Chromosome;
using file_reader_utils::Nucleotide;
using namespace map_utils;
using namespace math_utils;

namespace premeta {

// Overload cout's << operator so that I can print enum Sexptype's
// String representation instead of the enum (integer) value.
ostream& operator<<(std::ostream& out, const Sexptype value) {
  static map<Sexptype, string> sexptype_strings;
  if (sexptype_strings.size() == 0) {
/* NOTE: The following would work (after adding all enums, not just first three)
*  but I'd actually rather print a custom string, not the official enum name.
#define INSERT_ELEMENT(p) sexptype_strings[p] = #p
      INSERT_ELEMENT(SEXPTYPE_UNKNOWN);
      INSERT_ELEMENT(SEXPTYPE_SYMBOL);
      INSERT_ELEMENT(SEXPTYPE_PAIR_LIST);
#undef INSERT_ELEMENT
*/
    sexptype_strings.insert(make_pair(
        SEXPTYPE_UNKNOWN, "Unknown Sexptype"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_SYMBOL, "Symbol"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_PAIR_LIST, "Pair-List"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_CLOSURE, "Closure"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_ENVIRONMENT, "Environment"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_PROMISE, "Promise"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_LANG_OBJ, "Language Object"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_SPECIAL, "Special"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_BUILTIN, "Built-in"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_STRING, "String"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_BOOL, "Bool Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_UNUSED_ELEVEN, "Unused Sexptype 11"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_UNUSED_TWELVE, "Unused Sexptype 12"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_INT, "Integer Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_REAL, "Real Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_COMPLEX, "Complex Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_STRING, "String Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_DOTS, "Dots"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_ANY, "Any"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_LIST, "List"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_EXPRESSION, "Expression Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_BYTE_CODE, "Byte Code"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_POINTER, "Pointer"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_WEAK, "Weak"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_VECTOR_RAW, "Raw Vector"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_CLASS, "Class"));
    sexptype_strings.insert(make_pair(
        SEXPTYPE_PSUEDO, "Pseudo"));
  }
  if (value > 26) return out << "ERROR: Invalid Sexptype value: " << value;
  return out << sexptype_strings[value];
}

void AppendTagIfNotPresent(string* name, vector<string*>* tags) {
  for (const string* tag : *tags) {
    if (*tag == *name) return;
  }
  tags->push_back(name);
}

bool ExpectObjectType(
    const Sexptype type, const RDataObject& obj) {
  if (type == SEXPTYPE_STRING) {
    return (obj.str_ != nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_VECTOR_INT) {
    return (obj.str_ == nullptr && obj.int_vec_ != nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_VECTOR_BOOL) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ != nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_VECTOR_STRING) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ != nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_VECTOR_REAL) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ != nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_LIST) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ != nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ == nullptr);
  } else if (type == SEXPTYPE_PAIR_LIST) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ != nullptr && obj.class_ == nullptr);
  }  else if (type == SEXPTYPE_CLASS) {
    return (obj.str_ == nullptr && obj.int_vec_ == nullptr &&
            obj.bool_vec_ == nullptr && obj.str_vec_ == nullptr &&
            obj.real_vec_ == nullptr && obj.list_vec_ == nullptr &&
            obj.pair_list_vec_ == nullptr && obj.class_ != nullptr);
  }
  return false;
}

bool ReadRDataFile(
    const string& infile, vector<string*>* sexptype_tags, RDataNode* root) {
  // NOTE: See Note at top of file regarding using R Interals library.
  //SEXP foo, zed;
  //R_outpstream_t bar;
  //SEXP hungry = serializeToRaw(foo);
  // END NOTE.

  ifstream file(infile, ios::in|ios::binary|ios::ate);
  if (!file.is_open()) {
    cout << "ERROR: Unable to open file: " << infile << ". Aborting." << endl;
    return false;
  }
  // TODO(PHB): Handle case file_size > kMaxBufferSize.
  int64_t file_size = file.tellg();
  const unsigned long block_size =
      static_cast<unsigned long>(min(kMaxBufferSize, file_size));
  if (block_size != file_size) {
    cout << "ERROR: Unable to handle files with more than "
         << kMaxBufferSize << " characters." << endl;
    return false;
  }

  char* in_buffer = new char[block_size];
  file.seekg(0, ios::beg);
  // Read block.
  file.read(in_buffer, block_size);

  // We'll determine if the file is compressed (zipped) or not based on
  // reading the first 5 bytes (if file is uncompressed, first 5 bytes of
  // a .RData file should be 'RDX2\n').
  const char first_char = in_buffer[0];
  const char second_char = in_buffer[1];
  const char third_char = in_buffer[2];
  const char fourth_char = in_buffer[3];
  const char fifth_char = in_buffer[4];
  char* uncompressed_bytes;
  unsigned long num_bytes;
  // 82 = 'R', 68 = 'D', 88 = 'X', 50 = 2, 10 = '\n'; where "RDX2\n"
  // is the expected header at the top of any .Rdata binary file.
  if (first_char != 82 || second_char != 68 || third_char != 88 ||
       fourth_char != 50 || fifth_char != 10) {
    // Unzip File.
    if (!DecompressZippedBinaryFile(
            infile, in_buffer, block_size,
            &uncompressed_bytes, &num_bytes)) {
      free(in_buffer);
      free(uncompressed_bytes);
      file.close();
      return false;
    }
  } else {
    uncompressed_bytes = in_buffer;
    num_bytes = (unsigned long) block_size;
  }

  // Parse RData file.
  if (!ParseRDataFile(
          infile, uncompressed_bytes, num_bytes,
          sexptype_tags, root)) {
    if (in_buffer == uncompressed_bytes) {
      free(in_buffer);
    } else {
      free(in_buffer);
      free(uncompressed_bytes);
    }
    file.close();
    return false;
  }
  
  if (in_buffer == uncompressed_bytes) {
    free(in_buffer);
  } else {
    free(in_buffer);
    free(uncompressed_bytes);
  }
  file.close();

  return true;
}

bool DecompressZippedBinaryFile(
    const string& input_filename,
    char* in_buffer, const uint64_t block_size,
    char** uncompressed_bytes, unsigned long* num_bytes) {
  if (block_size == 0) {
    cout << "ERROR in decompressing zipped file: "
         << "Empty file '" << input_filename << "'. Aborting." << endl;
    return false;
  }
  
  unsigned full_length = block_size ;
  unsigned half_length = block_size / 2;
  
  unsigned uncompressed_length = full_length;
  *uncompressed_bytes = new char[uncompressed_length];
  
  z_stream strm;
  strm.next_in = (Bytef*) in_buffer;
  strm.avail_in = block_size;
  strm.total_out = 0;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  
  bool done = false;
  
  if (inflateInit2(&strm, (16 + MAX_WBITS)) != Z_OK) {
    free(*uncompressed_bytes);
    cout << "ERROR in decompressing zipped file: "
         << "Unable to decompress file: '" << input_filename
         << "'. Aborting." << endl;
    return false;
  }
  
  while (!done) {
    // Check if output buffer is too small.
    if (strm.total_out >= uncompressed_length) {
      // Increase size of output buffer.
      char* temp = new char[uncompressed_length + half_length];
      memcpy(temp, *uncompressed_bytes, uncompressed_length);
      uncompressed_length += half_length;
      free(*uncompressed_bytes);
      *uncompressed_bytes = temp;
    }

    strm.next_out = (Bytef *) (*uncompressed_bytes + strm.total_out);
    strm.avail_out = uncompressed_length - strm.total_out;
  
    // Inflate another chunk.
    int err = inflate (&strm, Z_SYNC_FLUSH);
    if (err == Z_STREAM_END) done = true;
    else if (err != Z_OK) {
      break;
    }
  }
  
  if (inflateEnd (&strm) != Z_OK) {
    free(*uncompressed_bytes);
    cout << "ERROR in decompressing zipped file: "
         << "Unable to decompress file: '" << input_filename
         << "'. Aborting." << endl;
    return false;
  }

  *num_bytes = strm.total_out;

  return true;
}

bool ValidCharPosition(
    const unsigned long& next_char, const unsigned long& max_char_pos) {
  // Sanity check next_char < max_char_pos and op is supported.
  if (next_char >= max_char_pos) {
    cout << "ERROR in parsing .RData file: end of file reached "
         << "unexpectedly (not closing 'END' block. Aborting." << endl;
    return false;
  }
  return true;
}

bool ParseRDataHeader(
    const char* uncompressed_bytes, string* obj_name) {
  // Parse Header.
  const string header(uncompressed_bytes, 7);
  if (header != "RDX2\nX\n") {
    cout << "ERROR in Parsing RData Header: expected file to start with "
         << "'RDX2\\nX\\n', but instead found '" << header << "'. Aborting.\n";
    return false;
  }

  // Parse Serialization.
  int serialization1 = uncompressed_bytes[7];
  int serialization2 = uncompressed_bytes[8];
  int serialization3 = uncompressed_bytes[9];
  int serialization4 = uncompressed_bytes[10];

  // Parse R Version Number used to write the RData file.
  int r_release1 = uncompressed_bytes[11];
  int r_release2 = uncompressed_bytes[12];
  int r_release3 = uncompressed_bytes[13];
  int r_release4 = uncompressed_bytes[14];

  // Parse R Version Number (the minimal version needed to read the RData file).
  int r_version1 = uncompressed_bytes[15];
  int r_version2 = uncompressed_bytes[16];
  int r_version3 = uncompressed_bytes[17];
  int r_version4 = uncompressed_bytes[18];

  cout << "Parsing SeqMeta file constructed under serialization: "
       << serialization1 << serialization2 << serialization3 << serialization4
       << ", under R Release: "
       << r_release1 << r_release2 << r_release3 << r_release4
       << ", and R Version: "
       << r_version1 << r_version2 << r_version3 << r_version4 << endl;

  // Ensure next bytes specify the global pair-list pack_flags.
  if ((int) uncompressed_bytes[19] != 0 || (int) uncompressed_bytes[20] != 0 ||
      (int) uncompressed_bytes[21] != 4 || (int) uncompressed_bytes[22] != 2) {
    cout << "ERROR in Parsing Header: Expected highest-level structure "
         << "to be a Pair-List with a Tag (encoding: '0042'), but found: '"
         << (int) uncompressed_bytes[19] << (int) uncompressed_bytes[20]
         << (int) uncompressed_bytes[21] << (int) uncompressed_bytes[22]
         << "'. Aborting." << endl;
    return false;
  }

  // Ensure next bytes specify '0001' for SEXPTYPE.
  if ((int) uncompressed_bytes[23] != 0 || (int) uncompressed_bytes[24] != 0 ||
      (int) uncompressed_bytes[25] != 0 || (int) uncompressed_bytes[26] != 1) {
    cout << "ERROR in Parsing Header: Expected type of highest level "
         << "object to be SEXPTYPE ('0001'), but instead found '"
         << (int) uncompressed_bytes[23] << (int) uncompressed_bytes[24]
         << (int) uncompressed_bytes[25] << (int) uncompressed_bytes[26]
         << "'. Aborting." << endl;
    return false;
  }

  // Next 4-bytes should specify a String. Ordinarily, it should be '0409',
  // which indicates a Latin string; but we demand here only a string '***9'.
  if ((int) uncompressed_bytes[30] != 9) {
    cout << "ERROR in Parsing Header: Expected String ('***9') descriptor for "
         << "tag of highest level object, but found: '"
         << (int) uncompressed_bytes[27] << (int) uncompressed_bytes[28]
         << (int) uncompressed_bytes[29] << (int) uncompressed_bytes[30]
         << "'. Aborting." << endl;
    return false;
  }

  // Get length of Object Name.
  int obj_length = ParseRDataInteger(uncompressed_bytes, 31);

  // Read Object Name.
  *obj_name = ParseRDataString(uncompressed_bytes, 35, obj_length);

  return true;
}

bool ParseRDataBoolean(
    const char* input, const unsigned long& start_pos) {
  return input[start_pos] != 0;
}

int ParseRDataInteger(
    const char* input, const unsigned long& start_pos) {
  unsigned char int_char[4];
  memcpy(int_char, &input[start_pos], 4);
  int lowest_bits = (int) int_char[3];
  int low_bits = (int) int_char[2];
  int high_bits = (int) int_char[1];
  int highest_bits = (int) int_char[0];
  return lowest_bits + (low_bits << 8) + (high_bits << 16) + (highest_bits << 24);
}

double ParseRDataFloat(
    const char* input, const unsigned long& start_pos) {
  double d;

  unsigned char b0, b1, b2, b3, b4, b5, b6, b7;
  memcpy(&b0, &input[start_pos], sizeof(b0));
  memcpy(&b1, &input[start_pos + 1], sizeof(b1));
  memcpy(&b2, &input[start_pos + 2], sizeof(b2));
  memcpy(&b3, &input[start_pos + 3], sizeof(b3));
  memcpy(&b4, &input[start_pos + 4], sizeof(b4));
  memcpy(&b5, &input[start_pos + 5], sizeof(b5));
  memcpy(&b6, &input[start_pos + 6], sizeof(b6));
  memcpy(&b7, &input[start_pos + 7], sizeof(b7));

  *((unsigned char*)(&d) + 7) = b0;
  *((unsigned char*)(&d) + 6) = b1;
  *((unsigned char*)(&d) + 5) = b2;
  *((unsigned char*)(&d) + 4) = b3;
  *((unsigned char*)(&d) + 3) = b4;
  *((unsigned char*)(&d) + 2) = b5;
  *((unsigned char*)(&d) + 1) = b6;
  *((unsigned char*)(&d) + 0) = b7;

  return d;
}

RStringType ParseStringType(const int level) {
  // 65696 = 11111111 10100001
  // Level should not have any bits NOT in positions 1, 2, 3, 4, or 6
  if (level & 65696) return R_STRING_NONE;
  const int is_na       = (level >> 1) & 1;
  const int is_byte     = (level >> 2) & 1;
  const int is_latin    = (level >> 3) & 1;
  const int is_utf      = (level >> 4) & 1;
  const int is_ascii    = (level >> 6) & 1;

  // Exactly one of the above bits should be set.
  if (is_na + is_byte + is_latin + is_utf + is_ascii != 1) {
    return R_STRING_NONE;
  }

  if (is_na) return R_STRING_NA;
  if (is_byte) return R_STRING_BYTE;
  if (is_latin) return R_STRING_LATIN;
  if (is_utf) return R_STRING_UTF8;
  if (is_ascii) return R_STRING_ASCII;

  return R_STRING_NONE;
}

string ParseRDataString(
    const char* input, const unsigned long& start_pos, const int length) {
  // Reserve one extra slot for the null-terminating character '\0'.
  char name[length + 1];
  memcpy(name, &input[start_pos], length);
  // Add null-terminating character, so conversion to string knows where
  // the char array terminates.
  name[length] = '\0';
  string to_return(name);
  return to_return;
}

bool ParseRDataSxptype(
    const char* uncompressed_bytes, const unsigned long& next_char,
    int* sxptype, int* second_byte) {
  const int flags = ParseRDataInteger(uncompressed_bytes, next_char);
  // A Pair-List's sxpinfo (the first block after a Pair-List ('0042') block)
  // should either be 1 or 255. We read it from the lowest byte.
  *sxptype = (flags & 255);
  // We may also need the next lowest byte.
  *second_byte = (flags >> 8 & 255);
  return true;
}

bool ParseRDataPackFlag(
    const char* uncompressed_bytes, const unsigned long& next_char,
    SexptypeInfo* info) {
  const int flags = ParseRDataInteger(uncompressed_bytes, next_char);
  info->type_ = static_cast<Sexptype>(flags & 31);
  info->is_obj_ = (flags >> 8 & 1);
  info->has_attr_ = (flags >> 9 & 1);
  info->has_tag_ = (flags >> 10 & 1);
  if (info->type_ == 9) {
    info->level_ = (flags >> 12) & 94;
  } else {
    info->level_ = (flags >> 12);
  }
  return true;
}

bool ParseRDataTag(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char, string* tag) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  SexptypeInfo info;
  // Parse Tag sxptype (should be STRING).
  if (!ParseRDataPackFlag(uncompressed_bytes, next_char, &info)) {
    return false;
  }
  // Check Tag sxptype is indeed STRING.
  if (info.type_ != SEXPTYPE_STRING) {
    cout << "ERROR in parsing RDatavTag at character: " << next_char
         << ". Expected sxptype STRING for Tag, instead found: "
         << info.type_ << ". Aborting." << endl;
    return false;
  }
  // Parse String (Tag).
  next_char += 4;
  if (!ParseRDataStringBlock(
          info, uncompressed_bytes, max_char_pos, next_char, tag)) {
    return false;
  }
  return true;
}

bool ParseRDataAttribute(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* attr) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  SexptypeInfo info;
  if (!ParseRDataPackFlag(uncompressed_bytes, next_char, &info)) {
    return false;
  }
  if (info.type_ != SEXPTYPE_PAIR_LIST) {
    cout << "ERROR in Parsing RData Attribute at character " << next_char
         << ". Expected Attribute to be a Pair-List, instead found sxptype: "
         << info.type_ << ". Aborting." << endl;
    return false;
  }

  next_char += 4;
  attr->obj_ = new RDataObject();
  return ParseRDataPairListBlock(
      info, uncompressed_bytes, max_char_pos, next_char, sexptype_tags, attr);
}

bool ParseRDataSxpInfoBlock(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataPairList* node) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  // Parse Type.
  int sxptype, second_byte;
  if (!ParseRDataSxptype(uncompressed_bytes, next_char,
                         &sxptype, &second_byte)) {
    return false;
  }

  if (sxptype == 1) {
    // Sxptype 0001 is the normal, expected type for a Pair-List.
    // It has Tag first, then the Pair-List Item. Parse Tag.
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
    sexptype_tags->push_back(node->tag_);
  } else if (sxptype == 255) {
    // Sxptype '0 0 X 255' is like '0001', except it indicates we should
    // copy the tag from the X^th sextype object.
    if (second_byte < 1 || second_byte > sexptype_tags->size()) {
      cout << "ERROR in parsing RData SxpInfoBlock at character: " << next_char
           << ". Cannot copy tag of the sxptype " << second_byte
           << " because this sxptype index is not available. Aborting." << endl;
      return false;
    } 
    node->tag_ = (*sexptype_tags)[second_byte - 1];
  }
  if (sxptype == 1 || sxptype == 255) {
    // Next for Sxptype = '0001', and first for Sxptype = '0 0 0 255', is
    // the Pair-List Item. Parse that.
    next_char += 4;
    node->obj_ = new RDataNode();
    if (!ParseRDataItemBlock(
            uncompressed_bytes, max_char_pos, next_char, sexptype_tags, node->obj_)) {
      return false;
    }
  } else {
    cout << "ERROR in parsing RData SxpInfoBlock: Uknown Pair-List "
         << "sxptype: " << sxptype << " at character: "
         << next_char << ". Aborting." << endl;
    return false;
  }
  return true;
}

bool ParseRDataItemBlock(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  SexptypeInfo info;
  if (!ParseRDataPackFlag(uncompressed_bytes, next_char, &info)) {
    return false;
  }

  node->obj_ = new RDataObject();
  next_char += 4;
  if (info.type_ == SEXPTYPE_PAIR_LIST) {
    return ParseRDataPairListBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_LIST) {
    return ParseRDataListBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_STRING) {
    node->obj_->str_ = new string();
    return ParseRDataStringBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        node->obj_->str_);
  } else if (info.type_ == SEXPTYPE_VECTOR_BOOL) {
    return ParseRDataVectorBoolBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_VECTOR_INT) {
    return ParseRDataVectorIntBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_VECTOR_REAL) {
    return ParseRDataVectorRealBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_VECTOR_STRING) {
    return ParseRDataVectorStringBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ == SEXPTYPE_CLASS) {
    return ParseRDataClassBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else if (info.type_ > 240 && info.type_ < 256) {
    return ParseRDataPseudoBlock(
        info, uncompressed_bytes, max_char_pos, next_char,
        sexptype_tags, node);
  } else {
    cout << "ERROR in parsing RData ItemBlock: Unexpected input parameters "
         << "(should never see this error) at character position: "
         << next_char << "; info.type_: " << info.type_ << ". Aborting" << endl;
    return false;
  }
}

bool ParseRDataStringBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    string* str) {
  RStringType string_type = ParseStringType(info.level_);
  // Currently, only R_STRING_ASCII String type is supported.
  if (string_type != R_STRING_ASCII) {
    cout << "ERROR in parsing RData StringBlock: Unexpected String type: "
         << string_type << " ( from level " << info.level_
         << ") at character position: " << next_char
         << ". Aborting." << endl;
    return false;
  }

  // Parse String Length.
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  const int str_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse String.
  next_char += 4;
  if (!ValidCharPosition(next_char + str_length, max_char_pos)) return false;
  *str = ParseRDataString(uncompressed_bytes, next_char, str_length);


  // Advance current position to four characters before the end of the string
  // (we subtract 4, because when all functions return, it is assumed that
  // next_char points to the start of the last block of 4 that they read).
  next_char += str_length - 4;
  return true;
}

bool ParseRDataVectorBoolBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  // Get Vector Length.
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  const int vec_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse Elements of the Vector.
  next_char += 4;
  // A double is expressed using 8 bits.
  if (!ValidCharPosition(next_char + 4 * vec_length, max_char_pos)) return false;
  node->obj_->bool_vec_ = new vector<bool>();
  for (int i = 0; i < vec_length; ++i) {
    node->obj_->bool_vec_->push_back(
        ParseRDataBoolean(uncompressed_bytes, next_char));
    next_char++;
  }

  // Advance current position to four characters before the end of the last
  // double in the vector (we subtract 4, because when all functions return,
  // it is assumed that next_char points to the start of the last block of
  // four that they read).
  next_char -= 4;

  // Parse Attributes, if relevant.
  if (info.has_attr_) {
    next_char += 4;
    node->attr_ = new RDataNode();
    if (!ParseRDataAttribute(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->attr_)) {
      return false;
    }
  }

  // Parse List Tag, if relevant.
  if (info.has_tag_) {
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
  }

  return true;
}

bool ParseRDataVectorIntBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  // Get Vector Length.
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  const int vec_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse Elements of the Vector.
  next_char += 4;
  // A double is expressed using 8 bits.
  if (!ValidCharPosition(next_char + 4 * vec_length, max_char_pos)) return false;
  node->obj_->int_vec_ = new vector<int>();
  for (int i = 0; i < vec_length; ++i) {
    node->obj_->int_vec_->push_back(
        ParseRDataInteger(uncompressed_bytes, next_char));
    next_char += 4;
  }

  // Advance current position to four characters before the end of the last
  // double in the vector (we subtract 4, because when all functions return,
  // it is assumed that next_char points to the start of the last block of
  // four that they read).
  next_char -= 4;

  // Parse Attributes, if relevant.
  if (info.has_attr_) {
    next_char += 4;
    node->attr_ = new RDataNode();
    if (!ParseRDataAttribute(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->attr_)) {
      return false;
    }
  }

  // Parse List Tag, if relevant.
  if (info.has_tag_) {
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
  }

  return true;
}

bool ParseRDataVectorRealBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  // Get Vector Length.
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  const int vec_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse Elements of the Vector.
  next_char += 4;
  // A double is expressed using 8 bits.
  if (!ValidCharPosition(next_char + 8 * vec_length, max_char_pos)) return false;
  node->obj_->real_vec_ = new vector<double>();
  for (int i = 0; i < vec_length; ++i) {
    node->obj_->real_vec_->push_back(
        ParseRDataFloat(uncompressed_bytes, next_char));
    next_char += 8;
  }

  // Advance current position to four characters before the end of the last
  // double in the vector (we subtract 4, because when all functions return,
  // it is assumed that next_char points to the start of the last block of
  // four that they read).
  next_char -= 4;

  // Parse Attributes, if relevant.
  if (info.has_attr_) {
    next_char += 4;
    node->attr_ = new RDataNode();
    if (!ParseRDataAttribute(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->attr_)) {
      return false;
    }
  }

  // Parse List Tag, if relevant.
  if (info.has_tag_) {
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
  }

  return true;
}

bool ParseRDataVectorStringBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  // Get Vector Length.
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  const int vec_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse Elements of the Vector.
  node->obj_->str_vec_ = new vector<string>();
  for (int i = 0; i < vec_length; ++i) {
    // Parse String Type.
    next_char += 4;
    if (!ValidCharPosition(next_char, max_char_pos)) return false;
    SexptypeInfo info;
    if (!ParseRDataPackFlag(uncompressed_bytes, next_char, &info)) {
      return false;
    }

    // Ensure type is String.
    if (info.type_ != SEXPTYPE_STRING) {
      cout << "ERROR in parsing RData VectorStringBlock at character position "
           << next_char << ": Expected String Type (" << SEXPTYPE_STRING
           << ") Pack_Flag block, but found type: " << info.type_
           << ". Aborting." << endl;
      return false;
    }

    next_char += 4;
    string value = "";
    if (!ParseRDataStringBlock(
            info, uncompressed_bytes, max_char_pos, next_char, &value)) {
      return false;
    }
    node->obj_->str_vec_->push_back(value);
  }

  // Parse Attributes, if relevant.
  if (info.has_attr_) {
    next_char += 4;
    node->attr_ = new RDataNode();
    if (!ParseRDataAttribute(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->attr_)) {
      return false;
    }
  }

  // Parse List Tag, if relevant.
  if (info.has_tag_) {
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
  }

  return true;
}

bool ParseRDataPairListBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  // Currently, for Pair-List SEXPTYPE, all non-low-level bits must be zero,
  // except for Tag, which should be true.
  if (info.is_obj_ || info.has_attr_ || !info.has_tag_ || info.level_ != 0) {
    cout << "ERROR in parsing RData PairListBlock: Unexpected SexptypeInfo "
         << "at character position " << next_char << ". Aborting." << endl;
    return false;
  }

  // Keep reading Pair-List siblings until reach an END block.
  int is_end = ParseRDataInteger(uncompressed_bytes, next_char);
  int i = 1;
  node->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  while (is_end != kItemEndMarker) {
    node->obj_->pair_list_vec_->push_back(new RDataPairList());

    // Next block: Sxptype.
    if (!ValidCharPosition(next_char, max_char_pos)) return false;
    if (!ParseRDataSxpInfoBlock(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->obj_->pair_list_vec_->back())) {
      return false;
    }
    next_char += 4;

    if (next_char >= max_char_pos) {
      cout << "ERROR in parsing RData PairListBlock: Reached end of file before "
           << "finding the End of Pair List marker '000254'. Current position: "
           << next_char << ". Aborting." << endl;
      return false;
    }
    is_end = ParseRDataInteger(uncompressed_bytes, next_char);
    if (is_end == kItemEndMarker) break;

    // To reach here, we must be in a set of sibling Pair-Lists. Make
    // sure this is indeed the case.
    SexptypeInfo info;
    if (!ParseRDataPackFlag(uncompressed_bytes, next_char, &info)) {
      return false;
    }
    if (info.type_ != SEXPTYPE_PAIR_LIST) {
      cout << "ERROR in parsing RData PairListBlock: On iteration " << i
           << " of Pair-List loop, expected next block to be END block "
           << "or another (sibling) Pair-List, but found instead block type: "
           << info.type_ << "; at position " << next_char << endl;
      return false;
    }
    next_char += 4;
    ++i;
  }

  return true;
}

bool ParseRDataListBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  // Parse List Length.
  int list_length = ParseRDataInteger(uncompressed_bytes, next_char);

  // Parse List Items.
  node->obj_->list_vec_ = new vector<RDataNode*>();
  for (int i = 0; i < list_length; ++i) {
    node->obj_->list_vec_->push_back(new RDataNode());
    if (info.is_obj_) node->obj_->list_vec_->back()->is_object_ = true;
    next_char += 4;
    if (!ParseRDataItemBlock(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->obj_->list_vec_->back())) {
      return false;
    }
  }

  // Parse Attributes, if relevant.
  if (info.has_attr_) {
    next_char += 4;
    node->attr_ = new RDataNode();
    if (!ParseRDataAttribute(
            uncompressed_bytes, max_char_pos, next_char,
            sexptype_tags, node->attr_)) {
      return false;
    }
  }

  // Parse List Tag, if relevant.
  if (info.has_tag_) {
    next_char += 4;
    node->tag_ = new string();
    if (!ParseRDataTag(
            uncompressed_bytes, max_char_pos, next_char, node->tag_)) {
      return false;
    }
  }

  return true;
}

bool ParseRDataClassBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  node->obj_ = new RDataObject();
  node->obj_->class_ = new RDataNode();
  if (info.is_obj_) node->obj_->class_->is_object_ = true;
  return ParseRDataItemBlock(
      uncompressed_bytes, max_char_pos, next_char, sexptype_tags,
      node->obj_->class_);
}

bool ParseRDataPseudoBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node) {
  if (!ValidCharPosition(next_char, max_char_pos)) return false;
  // I haven't seen blocks of this type used yet. Print error message.
  cout << "ERROR in parsing RData PseudoBlock at character " << next_char
       << ": Handling of blocks of sxptype 'Pseudo' not yet implemented. "
       << "Aborting." << endl;
  return false;
}

bool ParseRDataFile(
    const string& input_file,
    const char* uncompressed_bytes, unsigned long& num_bytes,
    vector<string*>* sexptype_tags, RDataNode* root) {

  /*
  // PHB Temp: writing decompressed file, as a sanity check.
  cout << "\nPrinting unzipped file as characters (int) to '"
       << "foo_seq_meta_input_as_chars.txt'" << endl;
  ofstream outfile;
  outfile.open("foo_seq_meta_input_as_chars.txt");
  for (int64_t i = 0; i < num_bytes; ++i) {
    unsigned char c = uncompressed_bytes[i];
    outfile << "Char " << i << ": " << (int)c << endl;
  }
  outfile.close();
  cout << "\nDone Printing unzipped file." << endl;
  // END PHB Temp.
  */

  if (num_bytes < 23) {
    cout << "ERROR in Parsing RData input file '" << input_file
         << "': (uncompressed) file has only "
         << num_bytes << " bytes. Aborting." << endl;
    return false;
  }

  // Make sure RData file has expected header.
  string obj_name;
  if (!ParseRDataHeader(uncompressed_bytes, &obj_name)) {
    return false;
  }

  // Parse rest of file. Start at the '0042' (Global Pair-List), 
  // which is the 20th character (char index 19).
  unsigned long next_char = 19;
  return ParseRDataItemBlock(
      uncompressed_bytes, num_bytes, next_char, sexptype_tags, root);
}

}  // namespace premeta
