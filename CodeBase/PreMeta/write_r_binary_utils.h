// Date: July 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)
//
// Description: Writes the RDataTree to a .RData object that can be loaded by R.

#ifndef WRITE_R_BINARY_UTILS_H
#define WRITE_R_BINARY_UTILS_H

#include "PreMeta/premeta_structures.h"

#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace premeta {

char* StringToCharPointer(const string& input, int* size);

void WriteBoolean(const bool to_write, ofstream& outfile);

void WriteInteger(const int to_write, ofstream& outfile);

void WriteFloat(const float& to_write, ofstream& outfile);

void WriteStringType(const RStringType type, ofstream& outfile);

void WriteStringBlocks(const string& to_write, ofstream& outfile);

void WriteString(const string& to_write, ofstream& outfile);

void WriteSxptype(const bool copy_tag, const int tag_index, ofstream& outfile);

void WritePackFlag(
    const Sexptype type, const bool is_object, const bool has_attribute,
    const bool has_tag, const int level, ofstream& outfile);

void WriteBoolVector(const vector<bool>& to_write, ofstream& outfile);

void WriteIntegerVector(const vector<int>& to_write, ofstream& outfile);

void WriteStringVector(const vector<string>& to_write, ofstream& outfile);

void WriteRealVector(const vector<double>& to_write, ofstream& outfile);

void WriteRDataHeader(ofstream& outfile);

void WriteRDataTree(
    vector<string>* tags, const RDataNode& root, ofstream& outfile);

void PrintRDataObject(
    const RDataObject& node, const string& tabs, ofstream& outfile);

void PrintRDataTree(
    const RDataNode& root, const string& tabs, ofstream& outfile);

void DeleteRDataObject(RDataObject* obj, const set<string*>& sexptype_tags);

void DeleteRDataTree(RDataNode* root, const set<string*>& sexptype_tags);

void DeleteRDataTags(vector<string*> tags);

}  // namespace premeta
#endif
