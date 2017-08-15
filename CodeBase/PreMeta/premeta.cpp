// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "premeta_constants.h"
#include "premeta_structures.h"
#include "premeta_utils.h"

#include "StringUtils/string_utils.h"

#include <cstdlib>
#include <iostream>
#include <string>

using premeta::AlleleInfo;
using premeta::DefaultSoftwareVersion;
using premeta::FileInfo;
using premeta::Position;
using premeta::PreMeta;
using premeta::SnpInfo;
using premeta::SoftwareVersion;
using namespace string_utils;
using namespace std;

enum Software {
  SOFTWARE_UNKNOWN,
  MASS,
  METASKAT,
  RAREMETAL,
  SEQMETA,
};

struct BlockInfo {
  FileInfo input_one_;
  FileInfo input_two_;
  FileInfo input_three_;
  FileInfo snp_info_;

  Software from_;
  SoftwareVersion v_from_;

  double rescale_;

  int window_size_;

  string output_one_;
  string output_two_;

  BlockInfo() {
    from_ = Software::SOFTWARE_UNKNOWN;
    output_one_ = "";
    output_two_ = "";

    rescale_ = 1.0;
    window_size_ = 1000000;  // 1 Million.
  }
};

void PrintUsage() {
  cout << "Expected usage:\n\t" << "premeta.exe --script SCRIPT_FILE "
       << "--format [MASS | RAREMETAL | seqMeta | MetaSKAT] "
       << "[--version NUM] [--output_dir PATH_TO_DIR]" << endl;
}

string GetDirectory(const string& path) {
  size_t last_dir = path.find_last_of("/");
  if (last_dir == string::npos) return "";
  return path.substr(0, last_dir + 1);
}

Software StringToSoftware(const string& software) {
  if (software == "MASS" || software == "Mass" ||
      software == "mass") {
    return Software::MASS;
  } else if (software == "MetaSkat" || software == "METASKAT" ||
             software == "MetaSKAT" || software == "metaskat" ||
             software == "Metaskat") {
    return Software::METASKAT;
  } else if (software == "SeqMeta" || software == "Seqmeta" ||
             software == "SEQMETA" || software == "seqMeta" ||
             software == "SeqMETA" || software == "seqmeta") {
    return Software::SEQMETA;
  } else if (software == "RAREMETAL" || software == "Raremetal" ||
             software == "RareMetal" || software == "RareMETAL" ||
             software == "raremetal") {
    return Software::RAREMETAL;
  }
  cout << "\nERROR Reading Command: After '--format' argument, expected "
       << "one of the valid softwares:\nMASS, RareMetal, seqMeta, or MetaSKAT.\n";
  return Software::SOFTWARE_UNKNOWN;
}

// Updates the input string by appending the given integer. If the string
// contains a '.' (representing a file suffix), then the integer is appeneded
// just before the first '.'.
string AppendIndex(const string& input, const int index) {
  size_t last_slash = input.rfind("/");
  string filename = input;
  string directory = "";
  if (last_slash != string::npos) {
    filename = input.substr(last_slash + 1);
    directory = input.substr(0, last_slash + 1);
  }
  size_t first_dot = filename.find(".");
  if (first_dot == string::npos) {
    return input + Itoa(index);
  } else {
    return (directory + filename.substr(0, first_dot) +
            Itoa(index) + filename.substr(first_dot));
  }
}

string SetGoldenSnpFilename(
    const string& output_dir, const string& global_snp_info_file) {
  if (global_snp_info_file.empty()) {
    return output_dir + "global_snp_info.txt";
  }

  size_t first_dot = global_snp_info_file.find(".");
  if (first_dot == string::npos) {
    return global_snp_info_file + "_new.txt";
  }

  size_t last_slash = global_snp_info_file.rfind("/");
  string filename = global_snp_info_file;
  string directory = "";
  if (last_slash != string::npos) {
    filename = global_snp_info_file.substr(last_slash + 1);
    directory = global_snp_info_file.substr(0, last_slash + 1);
  }
  return (directory + filename.substr(0, first_dot) +
          "_new" + filename.substr(first_dot));
}

bool ShouldCheckSnpRefAltOrdering(
    const Software to_s, const string& global_snp_info_file,
    const vector<BlockInfo> conversions,
    bool* check_snp_ref_alt_order) {
  // If target Software is RareMetal or MetaSkat, demand allele checking.
  if (to_s == Software::RAREMETAL || to_s == Software::METASKAT) {
    *check_snp_ref_alt_order = true;
  }

  // If Global Snp File was provided, demand allele checking.
  if (!global_snp_info_file.empty()) {
    *check_snp_ref_alt_order = true;
  }

  // If any of the (MASS or SEQMETA) studies included an allele file, demand
  // allele file. Further, keep track if any such studies *don't* include
  // an allele file, so we can prompt an error in case some do and some don't.
  bool exist_studies_w_o_allele_file = false;
  for (const BlockInfo& info : conversions) {
    // Nothing to check if source study is from RareMetal or MetaSkat.
    if (info.from_ == Software::RAREMETAL || info.from_ == Software::METASKAT) {
      continue;
    }

    if (info.snp_info_.name_.empty()) {
      exist_studies_w_o_allele_file = true;
    } else {
      *check_snp_ref_alt_order = true;
    }
  }

  // If check_snp_ref_alt_order is true, make sure all studies have the
  // requisite SNP Ref/Alt labelling.
  if (*check_snp_ref_alt_order && exist_studies_w_o_allele_file) {
    cout << "ERROR: All studies must have SNP info that labels "
         << "REF and ALT alleles (studies from RareMetal and MetaSkat "
         << "automatically include this information in the study output "
         << "files, but this information must be explicitly included "
         << "by the user for MASS and SeqMeta studies, via the "
         << "SNP_INFO keyword in the script file), make sure all "
         << "MASS and SeqMeta studies in the script file use the "
         << "SNP_INFO keyword that references a file that contains "
         << "REF/ALT information for all the SNPs.\n"
         << "NOTE: the requirement that all studies have SNP info that "
         << "labels REF and ALT alleles arises because one of the following "
         << "conditions is met:\n\t1) The target software is RareMetal or "
         << "MetaSkat\n\t2) A global snp info file was specified via the "
         << "--global_snp_info command-line argument\n\t3) At least one "
         << "of the studies in the script file specified a Snp Info file "
         << "via the SNP_INFO keyword\nIf you do not wish to do allele "
         << "checking, you can make sure all three of the above conditions "
         << "are not met." << endl;
    return false;
  }

  return true;
}

// Check to see if a file already exists of this name. If so, modify name by
// appending an integer. Returns false if name had to be modified.
bool EnsureNameIsUnique(string* input) {
  int index_to_append = 0;
  const string orig_input = *input;
  while (FILE* file = fopen(input->c_str(), "r")) {
    index_to_append++;
    *input = AppendIndex(orig_input, index_to_append);
    fclose(file);
  }
  return index_to_append == 0;
}

bool ParseCommandLineArgs(
    const int argc, char* argv[],
    string* script_file, string* output_dir, string* global_snp_info_file,
    Software* to, SoftwareVersion* v_to) {
  // Loop through command line arguments, starting at 1 (first arg is the command
  // to run the program).
  for (int i = 1; i < argc; ++i) {
    string arg = string(argv[i]);
    if (arg == "--script") {
      if (i == argc - 1) {
        cout << "\nERROR Reading Command: Expected argument after '--script'.\n";
        return false;
      }
      ++i;
      *script_file = string(argv[i]);
    } else if (arg == "--software" || arg == "--format") {
      if (i == argc - 1) {
        cout << "\nERROR Reading Command: Expected argument after '--format'.\n";
        return false;
      }
      ++i;
      const string software = string(argv[i]);
      *to = StringToSoftware(software);
      if (*to == Software::SOFTWARE_UNKNOWN) {
        return false;
      }
    } else if (arg == "--version") {
      if (i == argc - 1) {
        cout << "\nERROR Reading Command: Expected argument after '--version'.\n";
        return false;
      }
      ++i;
      const string version = string(argv[i]);
      if (!premeta::ParseSoftwareVersion(version, v_to)) {
        cout << "\nERROR Reading Command: Unable to parse version '"
             << version << "'" << endl;
        return false;
      }
    } else if (arg == "--output_dir") {
      if (i == argc - 1) {
        cout << "\nERROR Reading Command: Expected argument after '--output_dir'.\n";
        return false;
      }
      ++i;
      *output_dir = string(argv[i]);
      if (!HasSuffixString(*output_dir, "/")) {
        *output_dir = *output_dir + "/";
      }
    } else if (arg == "--global_snp_info") {
      if (i == argc - 1) {
        cout << "\nERROR Reading Command: Expected argument after '--output_dir'.\n";
        return false;
      }
      ++i;
      *global_snp_info_file = string(argv[i]);
    } else {
      cout << "\nERROR Reading Command: unrecognized argument: " << arg << endl;
      return false;
    }
  }

  return true;
}

bool IsBlockComplete(const Software to_s, const BlockInfo* info) {
  // On the first pass of the while loop in ParseScriptFile, info will be null.
  if (info == nullptr) return true;

  if (info->from_ == Software::SOFTWARE_UNKNOWN) {
    cout << "ERROR: Software type not specified for block." << endl;
    return false;
  } else if (info->from_ == Software::MASS) {
    if (info->input_one_.name_.empty()) {
      cout << "ERROR: Missing 'FILE' in Mass block." << endl;
      return false;
    }
    if (to_s == Software::RAREMETAL || to_s == Software::METASKAT) {
      if (info->snp_info_.name_.empty()) {
        cout << "ERROR: Missing 'SNP_INFO' in Mass block, which "
             << "is mandatory when converting MASS to RareMetal or "
             << "MetaSKAT." << endl;
        return false;
      }
    }
    if (!info->input_three_.name_.empty()) {
      cout << "ERROR: At most two input files expected for MASS." << endl;
      return false;
    }
    return true;
  } else if (info->from_ == Software::RAREMETAL) {
    if (to_s != Software::RAREMETAL && info->input_one_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_GROUP' in RareMetal block." << endl;
      return false;
    }
    if (info->input_two_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_SCORE' in RareMetal block." << endl;
      return false;
    }
    if (info->input_three_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_COV' in RareMetal block." << endl;
      return false;
    }
    return true;
  } else if (info->from_ == Software::SEQMETA) {
    if (info->input_one_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_RDATA' in seqMeta block." << endl;
      return false;
    }
    if (to_s == Software::RAREMETAL || to_s == Software::METASKAT) {
      if (info->snp_info_.name_.empty()) {
        cout << "ERROR: Missing 'SNP_INFO' in SeqMeta block, which "
             << "is mandatory when converting SeqMeta to RareMetal or "
             << "MetaSKAT." << endl;
        return false;
      }
    }
    if (!info->input_three_.name_.empty()) {
      cout << "ERROR: At most two input files expected for seqMeta." << endl;
      return false;
    }
    return true;
  } else if (info->from_ == Software::METASKAT) {
    if (info->input_one_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_MInfo' in MetaSKAT block." << endl;
      return false;
    }
    if (info->input_two_.name_.empty()) {
      cout << "ERROR: Missing 'FILE_MSSD' in MetaSKAT block." << endl;
      return false;
    }
    if (!info->input_three_.name_.empty()) {
      cout << "ERROR: Only two input files expected for MetaSKAT." << endl;
      return false;
    }
    return true;
  }
  cout << "ERROR: Unrecognized block in script file." << endl;
  return false;
}

bool IsAlleleBlockComplete(const FileInfo* info) {
  // On the first pass of the while loop in ParseScriptFile, info will be null.
  if (info == nullptr) return true;

  if (info->name_.empty()) {
    cout << "ERROR: Missing 'FILE' in Allele File block." << endl;
    return false;
  }
  return true;
}

void FillDefaultFileSettings(BlockInfo* info) {
  if (info->from_ == Software::MASS) {
    info->input_one_.delimiter_ = premeta::kMassDefaultDelimiter;
    info->input_one_.comment_char_ = premeta::kMassDefaultCommentChar;
    info->input_two_.delimiter_ = premeta::kAlleleFileDefaultDelimiter;
    info->input_two_.comment_char_ = premeta::kAlleleFileDefaultCommentChar;
  } else if (info->from_ == Software::METASKAT) {
    info->input_one_.delimiter_ = premeta::kMInfoDefaultDelimiter;
    info->input_one_.comment_char_ = premeta::kMInfoDefaultCommentChar;
  } else if (info->from_ == Software::RAREMETAL) {
    info->input_one_.delimiter_ =
        premeta::kRareMetalGroupFileDefaultDelimiter;
    info->input_one_.comment_char_ =
        premeta::kRareMetalGroupFileDefaultCommentChar;
    info->input_two_.delimiter_ =
        premeta::kRareMetalScoreFileDefaultDelimiter;
    info->input_two_.comment_char_ =
        premeta::kRareMetalScoreFileDefaultCommentChar;
    info->input_three_.delimiter_ =
        premeta::kRareMetalCovFileDefaultDelimiter;
    info->input_three_.comment_char_ =
        premeta::kRareMetalCovFileDefaultCommentChar;
  } else if (info->from_ == Software::SEQMETA) {
    info->input_two_.delimiter_ = premeta::kAlleleFileDefaultDelimiter;
    info->input_two_.comment_char_ = premeta::kAlleleFileDefaultCommentChar;
  }
}

string PreMetaStripQuotes(const string& input) {
  string to_return = input;
  if (input.length() >= 2) {
    if (HasPrefixString(input, "\"") &&
        HasSuffixString(input, "\"")) {
      to_return = input.substr(1, input.length() - 2);
    }
    if (HasPrefixString(input, "'") &&
        HasSuffixString(input, "'")) {
      to_return = input.substr(1, input.length() - 2);
    }
  }

  // Because we Stripped all whitespace from script file lines,
  // we may have inadvertently deleted the 'space' seperator.
  // We'll go ahead and assume this is the case, and return " "
  // here; if we're wrong, the user likely misformatted the
  // script file anyway.
  if (to_return.empty()) return " ";

  // Handle special characters.
  if (to_return == "\\t") return "\t";
  if (to_return == "\\n") return "\n";
  if (to_return == "\\s") return " ";

  // No other special characters handled; print warning if we
  // see escape character.
  if (HasPrefixString(to_return, "\\")) {
    cout << endl << "WARNING: Found string '" << to_return
         << "' in script file; the escape character is unsupported "
         << "by PreMeta, and is unlikely to work correctly. Please "
         << "update script file (and make appropriate changes to "
         << "other input files) to either not use escaped characters, "
         << "or use one of the supported escape characters: '"
         << "'\\s' (space), '\\t' (tab), or '\\n' (new-line)." << endl;
  }

  return to_return;
}

bool IsNewRareMetalVersion(const SoftwareVersion& version) {
  return
      version.base_ > 4 ||
      (version.base_ == 4 && version.sub_version_ > 13) ||
      (version.base_ == 4 && version.sub_version_ == 13 &&
       version.sub_sub_version_ >= 2);
}

bool CheckVersion(const SoftwareVersion& version,
                  const Software from_s) {
  if (from_s == Software::MASS) {
    if (version.base_ != 7 ||
        version.sub_version_ != 0 ||
        version.sub_sub_version_ != -1) {
      if (version.base_ < 7) {
        cout << endl << "ERROR: MASS Versions earlier than "
             << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_MASS)
             << " are not supported." << endl;
        return false;
      }
      cout << endl << "WARNING: PreMeta assumes Mass version "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_MASS)
           << ". Converting version " << PrintSoftwareVersion(version)
           << " is not guarenteed to work..." << endl;
    }
  } else if (from_s == Software::METASKAT) {
    if (version.base_ != 0 || version.sub_version_ != 40 ||
        version.sub_sub_version_ != -1) {
      cout << endl << "WARNING: PreMeta assumes MetaSKAT version "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_METASKAT)
           << ". Converting version " << PrintSoftwareVersion(version)
           << " is not guarenteed to work..." << endl;
    }
  } else if (from_s == Software::RAREMETAL) {
    if ((version.base_ != 4 ||
         version.sub_version_ != 13 ||
         version.sub_sub_version_ != 5) &&
        (version.base_ != 0 ||
         version.sub_version_ != 4 ||
         version.sub_sub_version_ != 0)) {
      cout << endl << "WARNING: PreMeta assumes RareMetalWorker version "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_RAREMETAL_OLD)
           << " or "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_RAREMETAL_NEW)
           << ". Converting version " << PrintSoftwareVersion(version)
           << " is not guarenteed to work..." << endl;
    }
  } else if (from_s == Software::SEQMETA) {
    if (version.base_ != 1 || version.sub_version_ != 5 ||
        version.sub_sub_version_ != -1) {
      cout << endl << "WARNING: PreMeta assumes seqMeta version "
           << PrintSoftwareVersion(
                 DefaultSoftwareVersion::SOFTWARE_VERSION_SEQMETA)
           << ". Converting version " << PrintSoftwareVersion(version)
           << " is not guarenteed to work..." << endl;
    }
  }
  return true;
}

bool ParseScriptFile(const string& script_file,
                     const Software to_s,
                     vector<BlockInfo>* conversions) {
  ifstream file(script_file.c_str());
  if (!file.is_open()) {
    cout << "ERROR: Unable to open script file '" << script_file
         << "'. Aborting." << endl;
    return false;
  }

  BlockInfo* current_info = nullptr;
  string input_line, line;
  int line_index = 0;
  int block_index = 0;
  while(getline(file, input_line)) {
    line_index++;

    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (input_line[input_line.length() - 1] == 13) {
      input_line = input_line.substr(0, input_line.length() - 1);
    }

    line.clear();
    RemoveAllWhitespace(input_line, &line);
    // Skip empty lines and comment lines.
    if (line.empty() || HasPrefixString(line, "#")) continue;
    if (HasPrefixString(line, "SOFTWARE=") ||
        HasPrefixString(line, "Software=") ||
        HasPrefixString(line, "software=") ||
        HasPrefixString(line, "FORMAT=") ||
        HasPrefixString(line, "Format=") ||
        HasPrefixString(line, "format=")) {
      block_index++;
      if (!IsBlockComplete(to_s, current_info)) {
        cout << endl << "ERROR in reading line " << line_index
             << " of Software block " << block_index << " of script file: "
             << "Encountered new 'SOFTWARE' tag, but not enough information "
             << "was provided in the previous block. Aborting." << endl;
        return false;
      }
      const size_t start_index = 1 + line.find("=");
      conversions->push_back(BlockInfo());
      current_info = &conversions->back();
      current_info->from_ = StringToSoftware(line.substr(start_index));
      if (current_info->from_ == Software::SOFTWARE_UNKNOWN) {
        cout << "ERROR in reading line " << line_index << " of Software block "
             << block_index << " of script file: "
             << "Unable to parse Software '" << line.substr(start_index)
             << "'. Options are: MASS, METASKAT, RAREMETAL, or SEQMETA."
             << "Aborting." << endl;
        return false;
      }
      FillDefaultFileSettings(current_info);
      continue;
    }
    if (HasPrefixString(line, "VERSION=")) {
      const size_t start_index = 1 + line.find("=");
      if (!premeta::ParseSoftwareVersion(
              line.substr(start_index), &current_info->v_from_)) {
        cout << "ERROR in reading line " << line_index << " of Software block "
             << block_index << " of script file: "
             << "Unable to parse Version '" << line.substr(start_index)
             << "'. Aborting." << endl;
        return false;
      }
      if (!CheckVersion(current_info->v_from_, current_info->from_)) {
        cout << "ERROR on line " << line_index << " of Software block "
             << block_index << " of script file: "
             << "Version " << line.substr(start_index)
             << " is not supported. Aborting." << endl;
        return false;
      }
    } else if (HasPrefixString(line, "SNP_INFO=") ||
               HasPrefixString(line, "Snp_Info=") ||
               HasPrefixString(line, "Snp_info=") ||
               HasPrefixString(line, "snp_info=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->snp_info_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "SNP_INFO_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->snp_info_.delimiter_ =
        PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "SNP_INFO_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->snp_info_.comment_char_ =
        PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "RESCALE=")) {
      const size_t start_index = 1 + line.find("=");
      double rescale = -1.0;
      if (!Stod(PreMetaStripQuotes(line.substr(start_index)), &rescale)) {
        cout << "ERROR in reading line " << line_index << " of Software block "
             << block_index << " of script file: "
             << "Unable to parse RESCALE '" << line.substr(start_index)
             << "' as a numeric value." << endl;
        return false;
      }
      current_info->rescale_ = rescale;
    } else if (HasPrefixString(line, "OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_one_ =
          PreMetaStripQuotes(line.substr(start_index));
      current_info->output_two_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "MASS_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_one_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "SCORE_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_one_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "COV_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_two_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "MINFO_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_one_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "MSSD_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_two_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "RDATA_OUT_FILE_PREFIX=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->output_one_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "WINDOW_SIZE=")) {
      const size_t start_index = 1 + line.find("=");
      int window_size;
      if (!Stoi(line.substr(start_index), &window_size)) {
        cout << "ERROR in reading line " << line_index << " of Software block "
             << block_index << " of script file: "
             << "Unable to parse WINDOW_SIZE: '" << line.substr(start_index)
             << "' as an integer." << endl;
        return false;
      }
      current_info->window_size_ = window_size;
    } else if (HasPrefixString(line, "FILE_GROUP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_GROUP_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_GROUP_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_SCORE=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_SCORE_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_SCORE_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_COV=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_three_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_COV_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_three_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_COV_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_three_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MINFO=") ||
               HasPrefixString(line, "FILE_MInfo=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MINFO_SEP=") ||
               HasPrefixString(line, "FILE_MInfo_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MINFO_COMMENT=") ||
               HasPrefixString(line, "FILE_MInfo_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MSSD=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MSSD_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_MSSD_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_two_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_RDATA=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.name_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_RDATA_SEP=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.delimiter_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else if (HasPrefixString(line, "FILE_RDATA_COMMENT=")) {
      const size_t start_index = 1 + line.find("=");
      current_info->input_one_.comment_char_ =
          PreMetaStripQuotes(line.substr(start_index));
    } else {
      cout << "ERROR in reading line " << line_index << " of Software block "
           << block_index << " of script file: "
           << "Unable to parse line: '" << line << "'"
           << "Aborting." << endl;
      return false;
    }
  }

  return true;
}

void PrintScriptFileMessage() {
  cout << endl
       << "See documentation at:\n\thttp://dlin.web.unc.edu/software/premeta/\n"
       << "for proper format of script file." << endl;
}

int main(int argc, char* argv[]) {
  cout << endl;
  
  // Parse Command-Line.
  string script_file = "";
  string output_dir = "";
  string global_snp_info_file = "";
  Software to_s = Software::SOFTWARE_UNKNOWN;
  SoftwareVersion to_v;
  if (!ParseCommandLineArgs(
          argc, argv, &script_file, &output_dir, &global_snp_info_file,
          &to_s, &to_v)) {
    PrintUsage();
    return -1;
  }

  // Handle global allele file.
  string golden_allele_file =
      SetGoldenSnpFilename(output_dir, global_snp_info_file);
  if (!global_snp_info_file.empty() &&
      !PreMeta::ParseGoldenSnpFile(global_snp_info_file)) {
    cout << "ERROR: Unable to parse global snp info file '"
         << global_snp_info_file << "'." << endl;
    return -1;
  }

  // Sanity-check command-line had required components.
  if (to_s == Software::SOFTWARE_UNKNOWN) {
    cout << "ERROR: You must specify the target software." << endl;
    PrintUsage();
    return -1;
  }
  if (script_file.empty()) {
    cout << "ERROR: You must specify a script file." << endl;
    PrintUsage();
    return -1;
  }

  // Set and check target version; either use default values, or sanity-check
  // user specified appropriate values.
  if (to_v.base_ == -1 && to_v.sub_version_ == -1 &&
      to_v.sub_sub_version_ == -1) {
    string software_str = "";
    string version_str = "";
    if (to_s == Software::MASS) {
      software_str = "MASS";
      version_str = PrintSoftwareVersion(
          DefaultSoftwareVersion::SOFTWARE_VERSION_MASS);
      to_v.base_ = 7;
      to_v.sub_version_ = 0;
    } else if (to_s == Software::METASKAT) {
      software_str = "MetaSKAT";
      version_str = PrintSoftwareVersion(
          DefaultSoftwareVersion::SOFTWARE_VERSION_METASKAT);
    } else if (to_s == Software::RAREMETAL) {
      software_str = "RareMetal";
      version_str = PrintSoftwareVersion(
          DefaultSoftwareVersion::SOFTWARE_VERSION_RAREMETAL_NEW);
      // Update to_v to indicate the New version of RareMetal.
      to_v.base_ = 4;
      to_v.sub_version_ = 13;
      to_v.sub_sub_version_ = 5;
    } else if (to_s == Software::SEQMETA) {
      software_str = "seqMeta";
      version_str = PrintSoftwareVersion(
          DefaultSoftwareVersion::SOFTWARE_VERSION_SEQMETA);
    }
    cout << endl << "WARNING: No output version specified. Printing output "
         << "format consistent with " << software_str << " version "
         << version_str << endl;
  }

  // Sanity-check target version is compatible with the specified software.
  // Abort if version is not compatible; output warning if version is newer
  // than the default ones.
  if (to_s == Software::RAREMETAL &&
      (to_v.base_ < 4 ||
       (to_v.base_ ==  4 && to_v.sub_version_ < 13) ||
       (to_v.base_ ==  4 && to_v.sub_version_ == 13 &&
        to_v.sub_sub_version_ < 5))) {
    cout << endl << "ERROR: PreMeta does not support conversion to versions "
         << "of RareMetal older than "
         << PrintSoftwareVersion(
                DefaultSoftwareVersion::SOFTWARE_VERSION_RAREMETAL_NEW)
         << ". Aborting." << endl;
    return -1;
  } else if (to_s == Software::RAREMETAL &&
             (to_v.base_ != 4 || to_v.sub_version_ != 13 ||
              to_v.sub_sub_version_ != 5)) {
    cout << endl << "WARNING: PreMeta currently supports RareMetal version "
         << PrintSoftwareVersion(
               DefaultSoftwareVersion::SOFTWARE_VERSION_RAREMETAL_NEW)
         << ", and may not be compatible with higher versions." << endl;
  } else if (to_s == Software::MASS && to_v.base_ < 7) {
    cout << endl << "ERROR: PreMeta does not support conversion to versions "
         << "of MASS older than "
         << PrintSoftwareVersion(
                DefaultSoftwareVersion::SOFTWARE_VERSION_MASS)
         << ". Aborting." << endl;
    return -1;
  } else if (to_s == Software::MASS &&
             (to_v.base_ > 7 || to_v.sub_version_ > 0)) {
    cout << endl << "WARNING: PreMeta currently supports MASS version "
         << PrintSoftwareVersion(
               DefaultSoftwareVersion::SOFTWARE_VERSION_MASS)
         << ", and may not be compatible with higher versions." << endl;
  } else if (to_s == Software::METASKAT &&
             (to_v.base_ != 0 || to_v.sub_version_ != 40)) {
    cout << endl << "WARNING: PreMeta currently supports METASKAT version "
         << PrintSoftwareVersion(
               DefaultSoftwareVersion::SOFTWARE_VERSION_METASKAT)
         << ", and may not be compatible with other versions." << endl;
  } else if (to_s == Software::SEQMETA &&
             (to_v.base_ != 1 || to_v.sub_version_ != 5)) {
    cout << endl << "WARNING: PreMeta currently supports SeqMETA version "
         << PrintSoftwareVersion(
               DefaultSoftwareVersion::SOFTWARE_VERSION_SEQMETA)
         << ", and may not be compatible with other versions." << endl;
  }

  // Parse Script File.
  vector<BlockInfo> conversions;
  if (!ParseScriptFile(script_file, to_s, &conversions)) {
    PrintScriptFileMessage();
    return -1;
  }

  // Only need to populate allele information if target software is RareMetal or
  // MetaSKAT.
  /* UPDATE: As of July 25, 2015, we decided to demand the user provide
   * an allele file explicitly, and NOT use other (RareMetalWorker,
   * MetaSKAT) files. This decision was made because:
   *   1) Whoever is running PreMeta to convert summary-statistic
   *      output data should have access to Allele info
   *   2) It is possible that Major/Minor Alleles are different in
   *      different studies, or perhaps just for a single patient.
   *      In either case, meta-analysis will give the wrong results
   *      if we ignore the difference and merge allele information.
   * Thus, we skip populating a global allele list, and instead populate
   * allele list for each block of output files in the script file.     
  AlleleInfo allele_info;
  map<Position, pair<string, string>> allele_map;
  if (to_s == Software::RAREMETAL || to_s == Software::METASKAT) {
    // Add all allele files explicitly given in script file.
    for (const FileInfo& info : allele_files) {
      allele_info.supp_files_.push_back(FileInfo());
      FileInfo& new_info = allele_info.supp_files_.back();
      new_info.name_ = info.name_;
      new_info.delimiter_ = info.delimiter_;
      new_info.comment_char_ = info.comment_char_;
    }

    // Add all allele info that can be parsed from RareMetal's score file and
    // MetaSKAT's MInfo file.
    for (const BlockInfo& info : conversions) {
      if (info.from_ == Software::RAREMETAL) {
        allele_info.raremetal_files_.push_back(FileInfo());
        FileInfo& new_info = allele_info.raremetal_files_.back();
        new_info.name_ = info.input_two_.name_;
        new_info.delimiter_ = info.input_two_.delimiter_;
        new_info.comment_char_ = info.input_two_.comment_char_;
      } else if (info.from_ == Software::METASKAT) {
        allele_info.metaskat_files_.push_back(FileInfo());
        FileInfo& new_info = allele_info.metaskat_files_.back();
        new_info.name_ = info.input_one_.name_;
        new_info.delimiter_ = info.input_one_.delimiter_;
        new_info.comment_char_ = info.input_one_.comment_char_;
      }
    }

    if (!PreMeta::GetAllelesFromAlleleInfo(allele_info, &allele_map)) {
      cout << "ERROR: Unable to Get Allele Information from the provided "
           << "files. Aborting." << endl;
      return -1;
    }
  }
  */

  // Make sure that all studies provide a Snp Info file, if required.
  bool check_snp_ref_alt_order = false;
  if (!ShouldCheckSnpRefAltOrdering(
          to_s, global_snp_info_file, conversions, &check_snp_ref_alt_order)) {
    return -1;
  }
  PreMeta::SetSnpRefAltCheckField(check_snp_ref_alt_order);

  // Process each block in the script file, performing the necessary
  // changes to the input files to produce the appropriate output files.
  int block_itr = 1;
  int from_and_to_match = 0;
  bool deleted_mass_script_file = false;
  vector<pair<set<string>, set<string>>> file_conversion_names;
  for (const BlockInfo& info : conversions) {
    if (to_s == Software::MASS) {
      const string block_output_dir = info.output_one_.empty() ?
          output_dir : GetDirectory(info.output_one_).empty() ?
          output_dir : "";
      string outfile_one = info.output_one_.empty() ?
          block_output_dir + "mass_score_file" +
          Itoa(block_itr - from_and_to_match) + ".txt" :
          block_output_dir + info.output_one_ + ".txt";
      string warning_message = "";
      if (!EnsureNameIsUnique(&outfile_one)) {
        warning_message = "\nWARNING: filename updated to '" + outfile_one +
                          "' to avoid collision with an existing file.\n";
      }
      const string mass_script_file =
          info.output_one_.empty() ? output_dir + "mass_script.txt" :
          GetDirectory(info.output_one_).empty() ?
          output_dir + "mass_script.txt" :
          GetDirectory(info.output_one_) + "mass_script.txt";
      // Delete existing script file, if it exists (PreMeta outputs a mass
      // script file, which is used in the Meta-Analysis phase of PreMeta.
      // In order to account for multiple studies, PreMeta appends to this
      // script file; so we make sure that it is empty before appending
      // for the first time).
      if (!deleted_mass_script_file) {
        ofstream out_file;
        out_file.open(mass_script_file.c_str(), ios::out);
        out_file.close();
        deleted_mass_script_file = true;
      }
      // Parse allele_file
      map<Position, SnpInfo> allele_map;
      if (!info.snp_info_.name_.empty()) {
        if (!PreMeta::GetAllelesFromAlleleFile(
                (info.from_ == Software::RAREMETAL ||
                 info.from_ == Software::METASKAT),
                block_itr, info.snp_info_, &allele_map)) {
          cout << "ERROR: Unable to Get Allele Information from the provided "
               << "file '" << info.snp_info_.name_ << "'. Aborting." << endl;
          return -1;
        }
      }
      if (info.from_ == Software::MASS && info.rescale_ == 1.0) {
        cout << endl << "Skipping conversion from MASS to MASS for file block "
             << block_itr << " of script file (since RESCALE not specified)."
             << endl;
        block_itr++;
        from_and_to_match++;
        continue;
      } else if (info.from_ == Software::MASS) {
        cout << endl << "Converting MASS files in block " << block_itr
             << " of script file to MASS format..."
             << warning_message << endl;
        if (!PreMeta::MassToMass(
                block_itr, info.rescale_, allele_map, info.input_one_,
                mass_script_file, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting MASS files in block " << block_itr
             << " of script file to MASS format." << endl;
      } else if (info.from_ == Software::METASKAT) {
        cout << endl << "Converting MetaSKAT files in block " << block_itr
             << " of script file to MASS format..."
             << warning_message << endl;
        if (!PreMeta::MetaSkatToMass(
                block_itr, info.rescale_, allele_map, info.input_one_,
                info.input_two_.name_, mass_script_file, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting MetaSKAT files in block " << block_itr
             << " of script file to MASS format." << endl;
      } else if (info.from_ == Software::RAREMETAL) {
        cout << endl << "Converting RAREMETAL files in block " << block_itr
             << " of script file to MASS format..."
             << warning_message << endl;
        if (!PreMeta::RareMetalToMass(
                block_itr, IsNewRareMetalVersion(info.v_from_),
                info.rescale_, allele_map, info.input_one_, info.input_two_,
                info.input_three_, mass_script_file,
                outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        from_files.insert(info.input_three_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting RAREMETAL files in block " << block_itr
             << " of script file to MASS format." << endl;
      } else if (info.from_ == Software::SEQMETA) {
        cout << endl << "Converting SeqMeta files in block " << block_itr
             << " of script file to MASS format..."
             << warning_message << endl;
        if (!PreMeta::SeqMetaToMass(
                block_itr, info.rescale_, allele_map, info.input_one_.name_,
                mass_script_file, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting SeqMeta files in block " << block_itr
             << " of script file to MASS format." << endl;
      }
    } else if (to_s == Software::METASKAT) {
      const string block_output_dir = info.output_one_.empty() ?
          output_dir : GetDirectory(info.output_one_).empty() ?
          output_dir : "";
      string outfile_one = info.output_one_.empty() ?
          block_output_dir + "meta_skat_" +
          Itoa(block_itr - from_and_to_match) + ".MInfo" :
          block_output_dir + info.output_one_ + ".MInfo";
      string outfile_two = info.output_two_.empty() ?
          block_output_dir + "meta_skat_" +
          Itoa(block_itr - from_and_to_match) + ".MSSD" :
          block_output_dir + info.output_two_ + ".MSSD";
      string warning_message = "";
      if (!EnsureNameIsUnique(&outfile_one)) {
        warning_message = "\nWARNING: filename updated to '" + outfile_one +
                          "' to avoid collision with an existing file.\n";
      }
      if (!EnsureNameIsUnique(&outfile_two)) {
        warning_message += "\nWARNING: filename updated to '" + outfile_two +
                           "' to avoid collision with an existing file.\n";
      }
      // Parse allele_file
      map<Position, SnpInfo> allele_map;
      if (!info.snp_info_.name_.empty()) {
        if (!PreMeta::GetAllelesFromAlleleFile(
                (info.from_ == Software::RAREMETAL ||
                 info.from_ == Software::METASKAT),
                block_itr, info.snp_info_, &allele_map)) {
          cout << "ERROR: Unable to Get Allele Information from the provided "
               << "file '" << info.snp_info_.name_ << "'. Aborting." << endl;
          return -1;
        }
      }
      if (info.from_ == Software::MASS) {
        cout << endl << "Converting MASS files in block " << block_itr
             << " of script file to METASKAT format..."
             << warning_message << endl;
        if (info.snp_info_.name_.empty()) {
          cout << "ERROR: Must specify SNP_INFO file (for Major/Minor Allele "
               << "information) when using PreMeta to convert from MASS "
               << "to METASKAT. Aborting." << endl;
          return -1;
        }
        if (!PreMeta::MassToMetaSkat(
                block_itr, info.rescale_, allele_map, info.input_one_,
                outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting MASS files in block " << block_itr
             << " of script file to METASKAT format." << endl;
      } else if (info.from_ == Software::METASKAT && info.rescale_ == 1.0) {
        cout << endl << "Skipping conversion from METASKAT to METASKAT for "
             << "file block " << block_itr << " of script file (since "
             << "RESCALE not specified)." << endl;
        block_itr++;
        from_and_to_match++;
        continue;
      } else if (info.from_ == Software::METASKAT) {
        cout << endl << "Converting METASKAT files in block " << block_itr
             << " of script file to METASKAT format..."
             << warning_message << endl;
        if (!PreMeta::MetaSkatToMetaSkat(
                block_itr, info.rescale_, allele_map, info.input_one_,
                info.input_two_.name_, outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        from_files.insert(info.input_three_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting METASKAT files in block " << block_itr
             << " of script file to METASKAT format." << endl;
      } else if (info.from_ == Software::RAREMETAL) {
        cout << endl << "Converting RAREMETAL files in block " << block_itr
             << " of script file to METASKAT format..."
             << warning_message << endl;
        if (!PreMeta::RareMetalToMetaSkat(
                block_itr, IsNewRareMetalVersion(info.v_from_),
                info.rescale_, allele_map, info.input_one_, info.input_two_,
                info.input_three_, outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        from_files.insert(info.input_three_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting RAREMETAL files in block " << block_itr
             << " of script file to METASKAT format." << endl;
      } else if (info.from_ == Software::SEQMETA) {
        cout << endl << "Converting SeqMeta files in block " << block_itr
             << " of script file to METASKAT format..."
             << warning_message << endl;
        if (info.snp_info_.name_.empty()) {
          cout << "ERROR: Must specify SNP_INFO file (for Major/Minor Allele "
               << "information) when using PreMeta to convert from SEQMETA "
               << "to METASKAT. Aborting." << endl;
          return -1;
        }
        if (!PreMeta::SeqMetaToMetaSkat(
                block_itr, info.rescale_, allele_map, info.input_one_.name_,
                outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting SeqMeta files in block " << block_itr
             << " of script file to METASKAT format." << endl;
      }
    } else if (to_s == Software::RAREMETAL) {
      const string block_output_dir = info.output_one_.empty() ?
          output_dir : GetDirectory(info.output_one_).empty() ?
          output_dir : "";
      string outfile_one = info.output_one_.empty() ?
          block_output_dir + "raremetal_" +
          Itoa(block_itr - from_and_to_match) + ".score.txt" :
          block_output_dir + info.output_one_ + ".score.txt";
      string outfile_two = info.output_two_.empty() ?
          block_output_dir + "raremetal_" +
          Itoa(block_itr- from_and_to_match) + ".cov.txt" :
          block_output_dir + info.output_two_ + ".cov.txt";
      const string score_files =
          info.output_one_.empty() ? output_dir + "score_files.txt" :
          GetDirectory(info.output_one_).empty() ?
          output_dir + "score_files.txt" :
          GetDirectory(info.output_one_) + "score_files.txt";
      const string cov_files =
          info.output_one_.empty() ? output_dir + "cov_files.txt" :
          GetDirectory(info.output_one_).empty() ?
          output_dir + "cov_files.txt" :
          GetDirectory(info.output_one_) + "cov_files.txt";
      string warning_message = "";
      if (!EnsureNameIsUnique(&outfile_one)) {
        warning_message = "\nWARNING: filename updated to '" + outfile_one +
                          "' to avoid collision with an existing file.\n";
      }
      if (!EnsureNameIsUnique(&outfile_two)) {
        warning_message += "\nWARNING: filename updated to '" + outfile_two +
                           "' to avoid collision with an existing file.\n";
      }
      // Parse allele_file
      map<Position, SnpInfo> allele_map;
      if (!info.snp_info_.name_.empty()) {
        if (!PreMeta::GetAllelesFromAlleleFile(
                (info.from_ == Software::RAREMETAL ||
                 info.from_ == Software::METASKAT),
                block_itr, info.snp_info_, &allele_map)) {
          cout << "ERROR: Unable to Get Allele Information from the provided "
               << "file '" << info.snp_info_.name_ << "'. Aborting." << endl;
          return -1;
        }
      }
      if (info.from_ == Software::MASS) {
        cout << endl << "Converting MASS files in block " << block_itr
             << " of script file to RAREMETAL format..."
             << warning_message << endl;
        if (info.snp_info_.name_.empty()) {
          cout << "ERROR: Must specify SNP_INFO file (for Major/Minor Allele "
               << "information) when using PreMeta to convert from MASS "
               << "to RAREMETAL. Aborting." << endl;
          return -1;
        }
        if (!PreMeta::MassToRareMetal(
                block_itr, info.window_size_, info.rescale_,
                allele_map, info.input_one_, score_files, cov_files,
                outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting MASS files in block " << block_itr
             << " of script file to RAREMETAL format." << endl;
      } else if (info.from_ == Software::METASKAT) {
        cout << endl << "Converting METASKAT files in block " << block_itr
             << " of script file to RAREMETAL format..."
             << warning_message << endl;
        if (!PreMeta::MetaSkatToRareMetal(
                block_itr, info.window_size_, info.rescale_, allele_map,
                info.input_one_, info.input_two_.name_,
                score_files, cov_files, outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting METASKAT files in block " << block_itr
             << " of script file to RAREMETAL format." << endl;
      } else if (info.from_ == Software::RAREMETAL &&
                 info.rescale_ == 1.0 &&
                 (IsNewRareMetalVersion(info.v_from_) ||
                  !IsNewRareMetalVersion(to_v))) {
        cout << endl << "Skipping conversion from RAREMETAL to RAREMETAL for "
             << "file block " << block_itr << " of script file (since "
             << "RESCALE not specified, and not converting from old version "
             << "0.4.0 to newer version 4.13.5)." << endl;
        block_itr++;
        from_and_to_match++;
        continue;
      } else if (info.from_ == Software::RAREMETAL) {
        cout << endl << "Converting RAREMETAL files in block " << block_itr
             << " of script file to RAREMETAL format..."
             << warning_message << endl;
        if (!PreMeta::RareMetalToRareMetal(
                block_itr, IsNewRareMetalVersion(info.v_from_), info.rescale_, allele_map,
                info.input_two_, info.input_three_,
                score_files, cov_files, outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting RAREMETAL files in block " << block_itr
             << " of script file to RAREMETAL format." << endl;
      } else if (info.from_ == Software::SEQMETA) {
        cout << endl << "Converting SeqMeta files in block " << block_itr
             << " of script file to RAREMETAL format..."
             << warning_message << endl;
        if (info.snp_info_.name_.empty()) {
          cout << "ERROR: Must specify SNP_INFO file (for Major/Minor Allele "
               << "information) when using PreMeta to convert from SEQMETA "
               << "to RAREMETAL. Aborting." << endl;
          return -1;
        }
        if (!PreMeta::SeqMetaToRareMetal(
                block_itr, info.window_size_, info.rescale_,
                allele_map, info.input_one_.name_, score_files, cov_files,
                outfile_one, outfile_two)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        to_files.insert(outfile_two);
        cout << "Done converting SeqMeta files in block " << block_itr
             << " of script file to RAREMETAL format." << endl;
      }
    } else if (to_s == Software::SEQMETA) {
      const string block_output_dir = info.output_one_.empty() ?
          output_dir : GetDirectory(info.output_one_).empty() ?
          output_dir : "";
      string outfile_one = info.output_one_.empty() ?
          block_output_dir + "seqMeta_" +
          Itoa(block_itr - from_and_to_match) + ".RData" :
          block_output_dir + info.output_one_ + ".RData";
      string warning_message = "";
      if (!EnsureNameIsUnique(&outfile_one)) {
        warning_message = "\nWARNING: filename updated to '" + outfile_one +
                          "' to avoid collision with an existing file.\n";
      }
      // Parse allele_file
      map<Position, SnpInfo> allele_map;
      if (!info.snp_info_.name_.empty()) {
        if (!PreMeta::GetAllelesFromAlleleFile(
                (info.from_ == Software::RAREMETAL ||
                 info.from_ == Software::METASKAT),
                block_itr, info.snp_info_, &allele_map)) {
          cout << "ERROR: Unable to Get Allele Information from the provided "
               << "file '" << info.snp_info_.name_ << "'. Aborting." << endl;
          return -1;
        }
      }
      if (info.from_ == Software::MASS) {
        cout << endl << "Converting MASS files in block " << block_itr
             << " of script file to SEQMETA format..."
             << warning_message << endl;
        if (!PreMeta::MassToSeqMeta(
                block_itr, info.rescale_, allele_map, info.input_one_, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting MASS files in block " << block_itr
             << " of script file to SEQMETA format." << endl;
      } else if (info.from_ == Software::METASKAT) {
        cout << endl << "Converting METASKAT files in block " << block_itr
             << " of script file to SEQMETA format..."
             << warning_message << endl;
        if (!PreMeta::MetaSkatToSeqMeta(
                block_itr, info.rescale_, allele_map, info.input_one_,
                info.input_two_.name_, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting METASKAT files in block " << block_itr
             << " of script file to SEQMETA format." << endl;
      } else if (info.from_ == Software::RAREMETAL) {
        cout << endl << "Converting RAREMETAL files in block " << block_itr
             << " of script file to SEQMETA format..."
             << warning_message << endl;
        if (!PreMeta::RareMetalToSeqMeta(
                block_itr, IsNewRareMetalVersion(info.v_from_),
                info.rescale_, allele_map, info.input_one_, info.input_two_,
                info.input_three_, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        from_files.insert(info.input_two_.name_);
        from_files.insert(info.input_three_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting RAREMETAL files in block " << block_itr
             << " of script file to SEQMETA format." << endl;
      } else if (info.from_ == Software::SEQMETA && info.rescale_ == 1.0) {
        cout << endl << "Skipping conversion from SEQMETA to SEQMETA for "
             << "file block " << block_itr << " of script file (since "
             << "RESCALE not specified)." << endl;
        block_itr++;
        from_and_to_match++;
        continue;
      } else if (info.from_ == Software::SEQMETA) {
        cout << endl << "Converting SEQMETA files in block " << block_itr
             << " of script file to SEQMETA format..."
             << warning_message << endl;
        if (!PreMeta::SeqMetaToSeqMeta(
                block_itr, info.rescale_, allele_map, info.input_one_.name_, outfile_one)) {
          cout << "ERROR in converting file block " << block_itr
               << " of script file. Aborting." << endl;
          return -1;
        }
        file_conversion_names.push_back(pair<set<string>, set<string>>());
        pair<set<string>, set<string>>& files = file_conversion_names.back();
        set<string>& from_files = files.first;
        set<string>& to_files = files.second;
        from_files.insert(info.input_one_.name_);
        to_files.insert(outfile_one);
        cout << "Done converting SEQMETA files in block " << block_itr
             << " of script file to SEQMETA format." << endl;
      }
    }
    block_itr++;
  }

  cout << "\nPrinting Global Allele File...\n";
  PreMeta::WriteGoldenAlleleFile(golden_allele_file);

  // Print any SNPs that were skipped by a study, due to mismatching
  // REF/ALT alleles.
  PreMeta::PrintInconsistentSnps(output_dir + "inconsistent_snp_file.txt");

  // Done! Print summary.
  cout << "\n\nDONE! Converted:\n\n";
  for (const pair<set<string>, set<string>>& pairs : file_conversion_names) {
    cout << "\t" << Join(pairs.first, ", ") << "  ->" << endl
         << "\t" << Join(pairs.second, ", ") << endl << endl;
  }
  return 0;
}
