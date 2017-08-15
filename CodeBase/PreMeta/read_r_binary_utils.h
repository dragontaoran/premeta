// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Reads RData files, parsing them as an RDataTree.

#ifndef READ_R_BINARY_UTILS_H
#define READ_R_BINARY_UTILS_H

#include "PreMeta/premeta_structures.h"

#include <set>
#include <string>
#include <vector>

namespace premeta {

ostream& operator<<(std::ostream& out, const Sexptype value);

void AppendTagIfNotPresent(string* name, vector<string*>* tags);

bool ExpectObjectType(const Sexptype type, const RDataObject& obj);

bool ReadRDataFile(
    const string& infile, vector<string*>* sexptype_tags, RDataNode* root);

bool ParseRDataFile(
    const string& input_file,
    const char* uncompressed_bytes, unsigned long& num_bytes,
    vector<string*>* sexptype_tags, RDataNode* root);

bool DecompressZippedBinaryFile(
    const string& input_filename,
    char* in_buffer, const uint64_t block_size,
    char** uncompressed_bytes, unsigned long* num_bytes);

bool ValidCharPosition(
    const unsigned long& next_char, const unsigned long& max_char_pos);
RStringType ParseStringType(const int level);

bool ParseRDataHeader(
    const char* uncompressed_bytes, string* obj_name);

bool ParseRDataBoolean(
    const char* input, const unsigned long& start_pos);
int ParseRDataInteger(
    const char* input, const unsigned long& start_pos);
double ParseRDataFloat(
    const char* input, const unsigned long& start_pos);
string ParseRDataString(
    const char* input, const unsigned long& start_pos, const int length);

bool ParseRDataSxptype(
    const char* uncompressed_bytes, const unsigned long& next_char,
    int* sxptype, int* second_byte);

bool ParseRDataPackFlag(
    const char* uncompressed_bytes, const unsigned long& next_char,
    SexptypeInfo* info);

bool ParseRDataTag(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char, string* tag);

bool ParseRDataAttribute(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* attr);

bool ParseRDataSxpInfoBlock(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataPairList* node);

bool ParseRDataItemBlock(
    const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataPairListBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataListBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataStringBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char, string* str);

bool ParseRDataVectorBoolBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataVectorIntBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataVectorRealBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataVectorStringBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataClassBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

bool ParseRDataPseudoBlock(
    const SexptypeInfo& info, const char* uncompressed_bytes,
    const unsigned long& max_char_pos, unsigned long& next_char,
    vector<string*>* sexptype_tags, RDataNode* node);

}  // namespace premeta
#endif
