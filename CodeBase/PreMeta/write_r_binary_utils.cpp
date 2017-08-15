// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "write_r_binary_utils.h"

#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;

namespace premeta {
  
char* StringToCharPointer(const string& input, int* size) {
  vector<char> holder;
  for (int i = 0; i < input.length(); ++i) {
    const char c = input[i];
    if (c == '\0') continue;
    holder.push_back(c);
  }
  *size = holder.size();
  return holder.data();
}

int FindTag(const string& tag, const vector<string>& tags) {
  if (tag.empty()) return -1;
  for (int i = 0; i < tags.size(); ++i) {
    if (tag == tags[i]) return i;
  }
  return -1;
}

void WriteBoolean(const bool to_write, ofstream& outfile) {
  char bool_char[1];
  bool_char[0] = to_write;
  outfile.write(bool_char, 1);
}

void WriteInteger(const int to_write, ofstream& outfile) {
  char int_char[4];
  int_char[0] = (to_write >> 24 & 255);
  int_char[1] = (to_write >> 16 & 255);
  int_char[2] = (to_write >> 8 & 255);
  int_char[3] = (to_write & 255);
  outfile.write(int_char, 4);
}

void WriteDouble(const double& to_write, ofstream& outfile) {
  char* temp_float_char = (char*) &to_write;

  // The above puts the bytes in the wrong order (Endian). Swap them.
  char float_char[8];
  for (int i = 0; i < 8; ++i) {
    float_char[i] = temp_float_char[7 - i];
  }

  outfile.write(float_char, 8);
}

void WriteStringType(const RStringType type, ofstream& outfile) {
  // '9' is the SXPTYPE for String, and the other bit describes the
  // encoding of the string (ASCII, LATIN, etc.) See bottom section of:
  // http://cran.r-project.org/doc/manuals/r-release/R-ints.html#Rest-of-header
  if (type == R_STRING_NA) {
    const int level = 9 + (1 << 16);
    WriteInteger(level, outfile);
  } else if (type == R_STRING_BYTE) {
    const int level = 9 + (1 << 17);
    WriteInteger(level, outfile);
  } else if (type == R_STRING_LATIN) {
    const int level = 9 + (1 << 18);
    WriteInteger(level, outfile);
  } else if (type == R_STRING_UTF8) {
    const int level = 9 + (1 << 19);
    WriteInteger(level, outfile);
  } else if (type == R_STRING_ASCII) {
    const int level = 9 + (1 << 22);
    WriteInteger(level, outfile);
  } else {
    cout << "ERROR Unable to parse RStringType: "
         << type << endl;
  }
}

void WriteString(const string& to_write, ofstream& outfile) {
  for (int i = 0; i < to_write.length(); ++i) {
    const unsigned char current = to_write.at(i);
    outfile << current;
  }
}

void WriteStringBlocks(const string& to_write, ofstream& outfile) {
  // Write String Type (Latin).
  WriteStringType(R_STRING_LATIN, outfile);

  // Write String Length.
  WriteInteger(to_write.length(), outfile);

  // Write String.
  WriteString(to_write, outfile);
}

void WriteSxptype(const bool copy_tag, const int tag_index, ofstream& outfile) {
  const int level = copy_tag ? 255 + (tag_index << 8) : 1;
  WriteInteger(level, outfile);
}

void WritePackFlag(
    const Sexptype type, const bool is_object, const bool has_attribute,
    const bool has_tag, const int level, ofstream& outfile) {
  const int to_write =
      (static_cast<int>(type)) +
      (static_cast<int>(is_object) << 8) +
      (static_cast<int>(has_attribute) << 9) +
      (static_cast<int>(has_tag) << 10) +
      (level << 12);
  WriteInteger(to_write, outfile);
}

void WriteEndBlock(ofstream& outfile) {
  WriteInteger(254, outfile);
}

void WriteBoolVector(const vector<bool>& to_write, ofstream& outfile) {
  // Write Vector length.
  WriteInteger(to_write.size(), outfile);

  // Write Vector contents.
  for (const bool value : to_write) {
    WriteBoolean(value, outfile);
  }
}

void WriteIntegerVector(const vector<int>& to_write, ofstream& outfile) {
  // Write Vector length.
  WriteInteger(to_write.size(), outfile);

  // Write Vector contents.
  for (const int value : to_write) {
    WriteInteger(value, outfile);
  }
}

void WriteStringVector(const vector<string>& to_write, ofstream& outfile) {
  // Write Vector length.
  WriteInteger(to_write.size(), outfile);

  // Write Vector contents.
  for (const string& value : to_write) {
    WriteStringBlocks(value, outfile);
  }
}

void WriteRealVector(const vector<double>& to_write, ofstream& outfile) {
  // Write Vector length.
  WriteInteger(to_write.size(), outfile);

  // Write Vector contents.
  for (const double& value : to_write) {
    WriteDouble(value, outfile);
  }
}

void WriteRDataHeader(ofstream& outfile) {
  // Write header: RDX2.
  outfile.write("RDX2\n", 5);

  // Write format: X.
  outfile.write("X\n", 2);

  // Write R version information
  WriteInteger(2, outfile);       // Serialization version: 2
  WriteInteger(196865, outfile);  // Current R version (3.1.1 in this case)
  WriteInteger(131840, outfile);  // Version number for R 2.3.0
}

void WriteRDataTree(
    vector<string>* tags, const RDataNode& root, ofstream& outfile) {
  // Write Node's Object.
  if (root.obj_ != nullptr) {
    const RDataObject& node = *(root.obj_);
    if (node.str_ != nullptr) {
      WriteStringBlocks(*(node.str_), outfile);
    } else if (node.class_ != nullptr) {
      // Level 10,000 marks S4 Class objects; and all that I have seen
      // set 'is_object' and 'has_attribute' bits to true.
      WritePackFlag(SEXPTYPE_CLASS, true, true, false, 16, outfile);
      //WritePackFlag(SEXPTYPE_LIST, true, true, false, 10000, outfile);
      WriteRDataTree(tags, *node.class_, outfile);
    } else if (node.int_vec_ != nullptr) {
      WritePackFlag(
          SEXPTYPE_VECTOR_INT,
          root.is_object_, root.attr_ != nullptr,
          (root.tag_ != nullptr && !root.tag_->empty()),
          0, outfile);
      WriteIntegerVector(*(node.int_vec_), outfile);
    } else if (node.bool_vec_ != nullptr) {
      WritePackFlag(
          SEXPTYPE_VECTOR_BOOL,
          root.is_object_, root.attr_ != nullptr,
          (root.tag_ != nullptr && !root.tag_->empty()),
          0, outfile);
      WriteBoolVector(*(node.bool_vec_), outfile);
    } else if (node.str_vec_ != nullptr) {
      WritePackFlag(
        SEXPTYPE_VECTOR_STRING,
          root.is_object_, root.attr_ != nullptr,
          (root.tag_ != nullptr && !root.tag_->empty()),
          0, outfile);
      WriteStringVector(*(node.str_vec_), outfile);
    } else if (node.real_vec_ != nullptr) {
      WritePackFlag(
          SEXPTYPE_VECTOR_REAL,
          root.is_object_, root.attr_ != nullptr,
          (root.tag_ != nullptr && !root.tag_->empty()),
          0, outfile);
      WriteRealVector(*(node.real_vec_), outfile);
    } else if (node.list_vec_ != nullptr) {
      WritePackFlag(
          SEXPTYPE_LIST,
          root.is_object_, root.attr_ != nullptr,
          (root.tag_ != nullptr && !root.tag_->empty()),
          0, outfile);
      WriteInteger(node.list_vec_->size(), outfile);
      for (int i = 0; i < node.list_vec_->size(); ++i) {
        const RDataNode& list_item = *((*(node.list_vec_))[i]);
        WriteRDataTree(tags, list_item, outfile);
      }
    } else if (node.pair_list_vec_ != nullptr) {
      for (int i = 0; i < node.pair_list_vec_->size(); ++i) {
        const RDataPairList& list_item = *((*(node.pair_list_vec_))[i]);
        
        // Write Pair-List SEXPTYPE.
        WritePackFlag(SEXPTYPE_PAIR_LIST, false, false, true, 0, outfile);
        
        // Must know if the tag is original before knowing which Sxptype to write.
        string tag = "";
        if (list_item.tag_ != nullptr) tag = *(list_item.tag_);
        const int tag_index = FindTag(tag, *tags);

        // Write Sxptype ('0 0 0 1' or '0 0 0 255').
        WriteSxptype(tag_index >= 0, tag_index + 1, outfile);

        // Write Tag. (Pair-Lists write tag second, and only if tag has not
        // already appeared as a tag for an earlier item).
        if (tag_index == -1) WriteStringBlocks(tag, outfile);

        // Store tag, if appropriate.
        if (!tag.empty() && tag_index == -1) tags->push_back(tag);
        
        // Write the actual object.
        const RDataNode* list_item_obj = list_item.obj_;
        WriteRDataTree(tags, *list_item_obj, outfile);
      }
      WriteEndBlock(outfile);
    }
  }

  // Write Node's Attribute.
  if (root.attr_ != nullptr) {
    WriteRDataTree(tags, *(root.attr_), outfile);
  }

  // Write Node's Tag.
  if (root.tag_ != nullptr) {
    string tag = *(root.tag_);
    const int tag_index = FindTag(tag, *tags);
    if (tag_index == -1) WriteStringBlocks(tag, outfile);
    // Store tag, if appropriate.
    if (!tag.empty() && tag_index == -1) tags->push_back(tag);
  }
}

void PrintRDataObject(
    const RDataObject& node, const string& tabs, ofstream& outfile) {
  if (node.str_ != nullptr && !node.str_->empty()) {
    outfile << "String: " << *(node.str_) << endl;
  } else if (node.class_ != nullptr) {
     outfile << "S4 Class:" << endl;
     PrintRDataTree(*node.class_, (tabs + "  "), outfile);
  } else if (node.int_vec_ != nullptr && !node.int_vec_->empty()) {
    outfile << "Integer Vector: [";
    for (int i = 0; i < node.int_vec_->size(); ++i) {
      if (i > 0) outfile << ", ";
      outfile << (*node.int_vec_)[i];
    }
    outfile << "]" << endl;
  } else if (node.bool_vec_ != nullptr && !node.bool_vec_->empty()) {
    outfile << "Boolean Vector: [";
    for (int i = 0; i < node.bool_vec_->size(); ++i) {
      if (i > 0) outfile << ", ";
      outfile << (*node.bool_vec_)[i];
    }
    outfile << "]" << endl;
  } else if (node.str_vec_ != nullptr && !node.str_vec_->empty()) {
    outfile << "String Vector: [";
    for (int i = 0; i < node.str_vec_->size(); ++i) {
      if (i > 0) outfile << ", ";
      outfile << (*node.str_vec_)[i];
    }
    outfile << "]" << endl;
  } else if (node.real_vec_ != nullptr && !node.real_vec_->empty()) {
    outfile << "Real Vector: [";
    for (int i = 0; i < node.real_vec_->size(); ++i) {
      if (i > 0) outfile << ", ";
      outfile << (*node.real_vec_)[i];
    }
    outfile << "]" << endl;
  } else if (node.list_vec_ != nullptr) {
    if (node.list_vec_->empty()) {
      outfile << "List: Empty" << endl;
    } else {
      outfile << "List of " << node.list_vec_->size() << " Elements:" << endl;
    }
    for (int i = 0; i < node.list_vec_->size(); ++i) {
      outfile << tabs << "List Item " << (i + 1) << " of "
              << node.list_vec_->size() << ":" << endl;
      PrintRDataTree(*((*(node.list_vec_))[i]), (tabs + "  "), outfile);
    }
  } else if (node.pair_list_vec_ != nullptr) {
    if (node.pair_list_vec_->empty()) {
      outfile << "Pair-List: Empty" << endl;
    } else {
      outfile << "Pair-List of " << node.pair_list_vec_->size() << " Elements:"
              << endl;
    }
    for (int i = 0; i < node.pair_list_vec_->size(); ++i) {
      outfile << tabs << "Pair-List Item " << (i + 1) << " of "
              << node.pair_list_vec_->size() << " (Tag: "
              << *(((*(node.pair_list_vec_))[i])->tag_) << "):" << endl;
      PrintRDataTree(*(((*(node.pair_list_vec_))[i])->obj_), (tabs + "  "), outfile);
    }
  } else {
    // This Object is empty. Simply end line, and return.
    outfile << "<Empty>" << endl;
  }
}

void PrintRDataTree(
    const RDataNode& root, const string& tabs, ofstream& outfile) {
  if (root.is_object_) outfile << tabs << "Is Object: True" << endl;
  if (root.tag_ != nullptr) outfile << tabs << "Tag: " << *(root.tag_) << endl;
  if (root.obj_ != nullptr) {
    outfile << tabs << "Object: ";
    PrintRDataObject(*(root.obj_), (tabs + "  "), outfile);
  }
  if (root.attr_ != nullptr) {
    outfile << tabs << "Attribute: " << endl;
    PrintRDataTree(*(root.attr_), (tabs + "  "), outfile);
  }
}

void DeleteRDataObject(
    RDataObject* obj, const set<string*>& sexptype_tags) {
  if (obj->str_ != nullptr) delete obj->str_;
  if (obj->int_vec_ != nullptr) delete obj->int_vec_;
  if (obj->bool_vec_ != nullptr) delete obj->bool_vec_;
  if (obj->str_vec_ != nullptr) delete obj->str_vec_;
  if (obj->real_vec_ != nullptr) delete obj->real_vec_;
  if (obj->list_vec_ != nullptr) {
    for (int i = 0; i < obj->list_vec_->size(); ++i) {
      DeleteRDataTree((*(obj->list_vec_))[i], sexptype_tags);
    }
    delete obj->list_vec_;
  }
  if (obj->pair_list_vec_ != nullptr) {
    for (int i = 0; i < obj->pair_list_vec_->size(); ++i) {
      DeleteRDataTree(((*(obj->pair_list_vec_))[i])->obj_, sexptype_tags);
    }
    delete obj->pair_list_vec_;
  }
}

void DeleteRDataTree(
    RDataNode* root, const set<string*>& sexptype_tags) {
  if (root->tag_ != nullptr &&
      sexptype_tags.find(root->tag_) == sexptype_tags.end()) {
    delete root->tag_;
  }
  if (root->attr_ != nullptr) DeleteRDataTree(root->attr_, sexptype_tags);
  if (root->obj_ != nullptr) DeleteRDataObject(root->obj_, sexptype_tags);
  delete root;
}

void DeleteRDataTags(vector<string*> tags) {
  for (string* tag : tags) delete tag;
}

}  // namespace premeta
