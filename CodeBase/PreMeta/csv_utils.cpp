// Author: paulbunn@email.unc.edu (Paul Bunn)
// Last Updated: March 2015

#include "csv_utils.h"

#include "StringUtils/string_utils.h"

#include <climits>  // For INT_MIN
#include <cfloat>   // For DBL_MIN
#include <iostream> // For cout (for debugging; can remove if not needed).
#include <fstream>  // For sprintf
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;
using string_utils::StringUtils;

namespace file_reader_utils {

namespace {

bool IsMissingValue(const string& input) {
  if (input.empty() || input == "NA" || input == "N/A" ||
      input == "na" || input == "n/a") {
    return true;
  }
  return false;
}

// Determines if line_num should be skipped, based on the input parameters.
// Note that there are some ambiguities, e.g. if line_num is in both
// lines_to_keep and lines_to_skip. The caller is responsible for making
// sure this doesn't happen; but these are handled by preferentially
// taking inclusion over exclusion.
bool ShouldSkipLine(
    const string& line, const string& comment_id,
    const int line_num, const bool has_header,
    const set<int>* lines_to_skip, const set<int>* lines_to_keep,
    const int range_to_keep_start, const int range_to_keep_end,
    const int range_to_skip_start, const int range_to_skip_end,
    bool* read_header_line) {
  // Skip empty lines, or lines that begin with comment_id.
  string trimmed_line;
  StringUtils::RemoveLeadingWhitespace(line, &trimmed_line);
  if (trimmed_line.empty()) {
    return true;
  }
  if (!comment_id.empty() &&
      StringUtils::HasPrefixString(trimmed_line, comment_id)) {
    return true;
  }

  // Skip header line.
  if (!(*read_header_line) && has_header) {
    *read_header_line = true;
    return true;
  }

  // Keep line if it is in lines_to_keep.
  if (lines_to_keep != nullptr &&
      (lines_to_keep->find(line_num) != lines_to_keep->end())) {
    return false;
  } else if (lines_to_keep != nullptr) {
    return true;
  }

  // Keep line if it is within range specified by
  //   [range_to_keep_start, range_to_keep_end].
  if (range_to_keep_start >= 0) {
    if (range_to_keep_end >= 0) {
      return !(line_num >= range_to_keep_start &&
               line_num <= range_to_keep_end);
    } else {
      return !(line_num >= range_to_keep_start);
    }
  } else if (range_to_keep_end >= 0) {
    return !(line_num <= range_to_keep_end);
  }

  // Skip line if it is in lines_to_skip.
  if (lines_to_skip != nullptr &&
      (lines_to_skip->find(line_num) != lines_to_skip->end())) {
    return true;
  } else if (lines_to_skip != nullptr) {
    return false;
  }

  // Skip line if it is within range specified by
  //   [range_to_skip_start, range_to_skip_end].
  if (range_to_skip_start >= 0) {
    if (range_to_skip_end >= 0) {
      return (line_num >= range_to_skip_start ||
              line_num <= range_to_skip_end);
    } else {
      return line_num >= range_to_skip_start;
    }
  } else if (range_to_skip_end >= 0) {
    return line_num <= range_to_skip_end;
  }
  
  // No filtering applies to this line_num.
  return false;
}

}  // namespace

bool CsvUtils::ReadCsv(
    const string& filename, const string& delimiter, const string& comment_id,
    const bool has_header,
    const double& na_double, const int na_int, const string& na_string,
    const set<int>* lines_to_skip, const set<int>* lines_to_keep,
    const int range_to_keep_start, const int range_to_keep_end,
    const int range_to_skip_start, const int range_to_skip_end,
    const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
    vector<vector<GenericDataHolder>>* output,
    vector<pair<int, string>>* error_lines, char* error_msg) {
  if (output == nullptr || columns_to_read.empty()) {
    sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Bad input.");
    return false;
  }
  output->clear();

  // Open File.
  ifstream file(filename.c_str());
  if (!file.is_open()) {
    sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Unable to open file '%s'.",
            filename.c_str());
    return false;
  }

  // Loop through lines of file.
  bool read_header_line = false;
  string line;
  int line_index = 1;
  int last_line_to_keep = -1;
  if (range_to_keep_end >= 0) {
    last_line_to_keep = range_to_keep_end;
  } else if (lines_to_keep != nullptr && !lines_to_keep->empty()) {
    last_line_to_keep = *((lines_to_keep->end())--);
  }
  while(getline(file, line)) {
    // Optimization to abort early if we've reached the end of the range
    // of line numbers to keep.
    if (last_line_to_keep >= 0 && line_index > last_line_to_keep) {
      break;
    }

    // Remove trailing character '13' (Carriage Return) at end of line, if present.
    if (line[line.length() - 1] == 13) {
      line = line.substr(0, line.length() - 1);
    }

    // Skip lines that shouldn't be processed (e.g. header line).
    if (ShouldSkipLine(line, comment_id, line_index, has_header,
                       lines_to_skip, lines_to_keep,
                       range_to_keep_start, range_to_keep_end,
                       range_to_skip_start, range_to_skip_end,
                       &read_header_line)) {
      line_index++;
      continue;
    }

    // Split line by column; initially reading it into a vector of strings.
    vector<string> columns;
    StringUtils::Split(
        line, delimiter, false /* Do Not Skip/Flatten Empty Columns */, &columns);

    // Pick out the relevant columns, and parse into the appropriate data type.
    vector<GenericDataHolder> row_values;
    for (const pair<vector<int>, GenericDataType>& itr : columns_to_read) {
      const GenericDataType& type = itr.second;
      bool within_range = false;
      int start_range = -1;
      int end_range = -1;
      for (const int columns_itr : itr.first) {
        // Special parsing needed for range-notation. See comments in csv_utils.h
        // for details.
        if (columns_itr == 0) {
          if (within_range) {
            // Nothing to do, ready to proceed to loop below.
          } else {
            within_range = true;
            continue;
          }
        } else if (within_range) {
          if (start_range == -1) {
            start_range = columns_itr;
            continue;
          } else if (end_range == -1) {
            end_range = columns_itr;
            continue;
          } else {
            sprintf(error_msg, "Error reading line %d: columns specifications "
                               "cannot be read. Specifications: [%s]. See "
                               "csv_utils for proper usage. Aborting.",
                    line_index, StringUtils::Join(itr.first, ", ").c_str());
            return false;
          }
        } else {
          start_range = columns_itr;
          end_range = columns_itr;
        }
        if (start_range >= 0 && end_range == -1) {
          end_range = columns.size();
        } else if (start_range == -1 || start_range > end_range) {
          sprintf(error_msg, "Error reading line %d: columns specifications "
                             "cannot be read. Specifications: [%s]. See "
                             "csv_utils for proper usage. Aborting.",
                  line_index, StringUtils::Join(itr.first, ", ").c_str());
          return false;
        }
        for (int column = start_range; column <= end_range; column++) {
          if (columns.size() < column) {
            sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Attempting to read "
                               "column %d of line %d: '%s', but only parsed %d "
                               "columns.",
                    column, line_index, line.c_str(), columns.size());
            return false;
          }

          // Parse the column according to it's type.
          GenericDataHolder holder;
          const string& col_value = columns[column - 1];
          if (IsMissingValue(col_value)) {
            if (type == GenericDataType::STRING) holder.str_ = na_string;
            else if (type == GenericDataType::INT) holder.int_ = na_int;
            else if (type == GenericDataType::INT_64) holder.int_ = na_int;
            else if (type == GenericDataType::UINT_64) holder.int_ = na_int;
            else if (type == GenericDataType::DOUBLE) holder.dbl_ = na_double;
          } else if (type == GenericDataType::STRING) {
            holder.str_ = col_value;
          } else if (type == GenericDataType::INT) {
            if (!StringUtils::Stoi(col_value, &holder.int_)) {
              if (error_lines == nullptr) {
                sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Attempting to read "
                                   "column %d of line %d ('%s'): expected INT value, "
                                   "but observed '%s'.",
                      column, line_index, line.c_str(), col_value.c_str());
                return false;
              } else {
                error_lines->push_back(make_pair(line_index, col_value));
                holder.int_ = na_int;
              }
            }
          } else if (type == GenericDataType::INT_64) {
            if (!StringUtils::Stoi(col_value, &holder.int64_)) {
              if (error_lines == nullptr) {
                sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Attempting to read "
                                   "column %d of line %d ('%s'): expected INT value, "
                                   "but observed '%s'.",
                      column, line_index, line.c_str(), col_value.c_str());
                return false;
              } else {
                error_lines->push_back(make_pair(line_index, col_value));
                holder.int_ = na_int;
              }
            }
          } else if (type == GenericDataType::UINT_64) {
            if (!StringUtils::Stoi(col_value, &holder.uint64_)) {
              if (error_lines == nullptr) {
                sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Attempting to read "
                                   "column %d of line %d ('%s'): expected INT value, "
                                   "but observed '%s'.",
                      column, line_index, line.c_str(), col_value.c_str());
                return false;
              } else {
                error_lines->push_back(make_pair(line_index, col_value));
                holder.int_ = na_int;
              }
            }
          } else if (type == GenericDataType::DOUBLE) {
            if (!StringUtils::Stod(col_value, &holder.dbl_)) {
              if (error_lines == nullptr) {
                sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Attempting to read "
                                   "column %d of line %d ('%s'): expected DOUBLE "
                                   "value, but observed '%s'.",
                      column, line_index, line.c_str(), col_value.c_str());
                return false;
              } else {
                error_lines->push_back(make_pair(line_index, col_value));
                holder.dbl_ = na_double;
              }
            }
          } else {
            sprintf(error_msg, "ERROR in CsvUtils::ReadCsv. Unsupported "
                               "GenericDataType: %d.", type);
            return false;
          }
          row_values.push_back(holder);
        }
        within_range = false;
        start_range = -1;
        end_range = -1;
      }
    }
    output->push_back(row_values);
    line_index++;
  }
  return true;
}

bool CsvUtils::ReadCsv(
    const string& filename, const string& delimiter, const string& comment_id,
    const bool has_header,
    const double& na_double, const int na_int, const string& na_string,
    const set<int>* lines_to_skip, const set<int>* lines_to_keep,
    const int range_to_keep_start, const int range_to_keep_end,
    const int range_to_skip_start, const int range_to_skip_end,
    const vector<pair<int, GenericDataType>>& columns_to_read,
    vector<vector<GenericDataHolder>>* output,
    vector<pair<int, string>>* error_lines, char* error_msg) {
  // Convert the intermediate data structure to the one expected by the
  // other ReadCsv function.
  vector<pair<vector<int>, GenericDataType>> new_format_cols_to_read;
  for (pair<int, GenericDataType> item : columns_to_read) {
    vector<int> vector_from_int;
    vector_from_int.push_back(item.first);
    new_format_cols_to_read.push_back(make_pair(
          vector_from_int, item.second));
  }
  return ReadCsv(
      filename, delimiter, comment_id, has_header, na_double, na_int, na_string,
      lines_to_skip, lines_to_keep, range_to_keep_start, range_to_keep_end,
      range_to_skip_start, range_to_skip_end, new_format_cols_to_read,
      output, error_lines, error_msg);
}
}  // namespace file_reader_utils
