// Author: paulbunn@email.unc.edu (Paul Bunn)
// Last Updated: March 2015
//
// Description: Helper functions for reading in .csv files (or more
// generally, and file in table format).

#ifndef CSV_UTILS_H
#define CSV_UTILS_H

#include "StringUtils/string_utils.h"

#include <climits>  // For INT_MIN
#include <cfloat>   // For DBL_MIN
#include <set>
#include <string.h>
#include <vector>

using namespace std;
using string_utils::StringUtils;

namespace file_reader_utils {

// An enum representing the various data types that a column may take.
enum GenericDataType {
  STRING,
  DOUBLE,
  INT,
  INT_64,
  UINT_64,
};

// Holds a data value of one of the three basic types (string, int, double).
// This generic struct allows us to construct a data structure (e.g. a vector)
// of values of different underlying types (e.g. create a vector of GenericDataHolder).
struct GenericDataHolder {
 public:
  string str_;
  double dbl_;
  int int_;
  int64_t int64_;
  uint64_t uint64_;
};

// A class defining basic helper functions for getting string/double values
// out of the various data holder objects above.
class CsvUtils {
 public:
  // Reads filename (an input file in csv (or table) format by reading in
  // the columns specified by columns_to_read, and populating 'output' with
  // the results. For example, on input file:
  //   AGE  CITY  IQ  HEIGHT
  //   34   L.A.  125 6.1
  //   22   NYC   65  6.5
  //   9    D.C.  87  3.9
  // and columns_to_read = [([1, 3], INT), ([4], DOUBLE)]
  // output would be: [ [int_ = 34, int_ = 125, dbl_ = 6.1],
  //                    [int_ = 22, int_ = 65,  dbl_ = 6.5],
  //                    [int_ = 9,  int_ = 87,  dbl_ = 3.9] ]
  // columns_to_read allows specifying ranges. To do this, pad each desired
  // range with '0's (as a hack to indicate the values in between are a range);
  // where you can indicate no end to the range by having only one value in
  // between the '0'-padding. For example, if you were interested in columns:
  //   1-10, 14, 16, 18, 20-22, 24-28, and 30 - END_OF_COLUMS,
  // (say all these columns have same type; columns of other types are handled
  // similarly) then you would specify for this GenericDataType the vector:
  //   [0, 1, 10, 0, 14, 16, 18, 0, 20, 22, 0, 0, 24, 28, 0, 0, 30, 0].
  //
  // Additional options include column delimiter in input file, whether
  // the input file has a header row, symbol marking the start of comment lines
  // (these lines are ignored) and values to use for each data type when a row
  // is missing a value.
  // Also, if error_lines is NULL, then returns false whenever a column cannot
  // be parsed into the desired data type; but if non-null, then error lines,
  // together with the (first) column string that caused the problem, are put
  // into error_lines and ReadCsv still returns true (assuming no other errors).
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const double& na_double, const int na_int, const string& na_string,
      const set<int>* lines_to_skip, const set<int>* lines_to_keep,
      const int range_to_keep_start, const int range_to_keep_end,
      const int range_to_skip_start, const int range_to_skip_end,
      const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output,
      vector<pair<int, string>>* error_lines, char* error_msg);
  // Same as above, but uses default values for all the ranges.
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const double& na_double, const int na_int, const string& na_string,
      const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output,
      vector<pair<int, string>>* error_lines, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header,
        na_double, na_int, na_string, nullptr, nullptr, -1, -1, -1, -1,
        columns_to_read, output, error_lines, error_msg);
  }
  // Same as above using default values for na_double (DBL_MIN),
  // na_string ("NA"), and na_int (INT_MIN).
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output,
      vector<pair<int, string>>* error_lines, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header, DBL_MIN, INT_MIN, "NA",
        columns_to_read, output, error_lines, error_msg);
  }
  // Same as above with nullptr for error_lines.
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header, columns_to_read,
        output, nullptr, error_msg);
  }
  // Same as above with empty string for comment_id.
  static bool ReadCsv(
      const string& filename, const string& delimiter, const bool has_header,
      const vector<pair<vector<int>, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output, char* error_msg) {
    return ReadCsv(
        filename, delimiter, "", has_header, columns_to_read, output, error_msg);
  }
  // Same as above, but with less flexibility: must explicitly specify every
  // column that you desire to read (the above API allows a more general
  // mechanism for specifying columns of interest, e.g. by allowing ranges).
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const double& na_double, const int na_int, const string& na_string,
      const set<int>* lines_to_skip, const set<int>* lines_to_keep,
      const int range_to_keep_start, const int range_to_keep_end,
      const int range_to_skip_start, const int range_to_skip_end,
      const vector<pair<int, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output,
      vector<pair<int, string>>* error_lines, char* error_msg);
  // Same as above, but uses default values for all the ranges.
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const double& na_double, const int na_int, const string& na_string,
      const vector<pair<int, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output,
      vector<pair<int, string>>* error_lines, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header,
        na_double, na_int, na_string, nullptr, nullptr, -1, -1, -1, -1,
        columns_to_read, output, error_lines, error_msg);
  }
  // Same as above, but passes in nullptr for error_lines (so any unparsable
  // column in any row will cause ReadCsv to return false).
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const double& na_double, const int na_int, const string& na_string,
      const vector<pair<int, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header, na_double, na_int,
        na_string, columns_to_read, output, nullptr, error_msg);
  }
  // Same as above using default values for na_double (DBL_MIN),
  // na_string ("NA"), and na_int (INT_MIN).
  static bool ReadCsv(
      const string& filename, const string& delimiter, const string& comment_id,
      const bool has_header,
      const vector<pair<int, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output, char* error_msg) {
    return ReadCsv(
        filename, delimiter, comment_id, has_header, DBL_MIN, INT_MIN, "NA",
        columns_to_read, output, error_msg);
  }
  // Same as above using default values for delimiter (\t), has_header (true),
  // comment_id (""), na_double (DBL_MIN), na_string ("NA"), and na_int (INT_MIN).
  static bool ReadCsv(
      const string& filename,
      const vector<pair<int, GenericDataType>>& columns_to_read,
      vector<vector<GenericDataHolder>>* output, char* error_msg) {
    return ReadCsv(
        filename, "\t", "", true, DBL_MIN, INT_MIN, "NA",
        columns_to_read, output, error_msg);
  }
};

}  // namespace file_reader_utils

#endif
