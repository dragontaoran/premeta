#include "read_table_with_header.h"

#include "FileReaderUtils/read_input.h"
#include "MathUtils/data_structures.h"
#include "StringUtils/string_utils.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <fstream>   // For sprintf and fstream.
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using file_reader_utils::ReadInput;
using math_utils::CensoringData;
using math_utils::DataHolder;
using math_utils::LinearTerm;
using math_utils::VariableTerm;
using string_utils::StringUtils;

namespace file_reader_utils {

namespace {

bool AddToken(
    const string& token, const int current_col_index, const int line_num,
    vector<DataHolder>* values, map<int, set<string> >* nominal_columns,
    char* error_msg) {
  
  DataHolder item;
  
  // Check if this column represents a nominal variable (by seeing if column was
  // already identified as such, or if the value cannot be parsed as a double).
  map<int, set<string> >::iterator itr =
      nominal_columns->find(current_col_index);
  if (itr != nominal_columns->end() ||
      !StringUtils::Stod(token, &(item.value_))) {
    // Set item.name_, to indicate this value is nominal (not numeric).
    item.name_ = token;

    // Check if this is the first time we've identified this column as nominal.
    if (token == "NA" || token == "N/A") {
      // Don't need to do anything; having this condition here makes it cleaner
      // than wrapping the following code (which should be skipped if token
      // indicates N/A) in an 'if' condition.
    } else if (itr == nominal_columns->end()) {
      // If this is not the first row (case) of data values, log an error.
      // (NOTE: This could happen if the user didn't mark the column as nominal
      // via '$' and the first value(s) read from the column could be parsed
      // as a double. For simplicity, we reject such input formats, and demand
      // that the user either remove ambiguous data values, and/or append '$'
      // to the title/header row of the data to indicate all values in that column
      // are nominal).
      if (line_num != 2) {
        sprintf(error_msg, "on line %d, column %d:\n\tExpected numerical "
                           "value, but found '%s'. Must update data file, "
                           "either by indicating column %d represents a "
                           "nominal value (by appending '$' to the column name "
                           "or removing all numerical-appearing entries from "
                           "this column) or update %s to be a numerical value.",
                line_num, (current_col_index + 1), token.c_str(),
                (current_col_index + 1), token.c_str());
        return false;
      }

      // Add this column to nominal_columns.
      set<string> value_to_add;
      value_to_add.insert(token);
      nominal_columns->insert(make_pair(current_col_index, value_to_add));
    } else {
      // This column was already known to be nominal. Add the new value to the
      // set of values for this column.
      set<string>& existing_values = itr->second;
      existing_values.insert(token);
    }
  }

  // Add item to values.
  if (values->size() > current_col_index) {
    (*values)[current_col_index] = item;
  } else {
    values->push_back(item);
  }
  return true;
}

}  // namespace

// Returns info on the appropriate format for the input file.
string ReadTableWithHeader::PrintInputFormat() {
  return  "\nNote that input.txt should have format:\n  "
          "\tX_1\tX_2$\t...\tX_p\n"
          "\tX_1,1\tX_2,1\t...\tX_p,1\n"
          "\tX_1,2\tX_2,2\t...\tX_p,2\n"
          "\t...\t...\t...\t...\n"
          "\tX_1,n\tX_2,n\t...\tX_p,n\n"
          "Where the first row is the 'TITLE' row, with titles for "
          "each of the (in)dependent variables. Titles that have a "
          "'$' suffix are considered NOMINAL variables. The distinction "
          "between independent vs. dependent variables, as well as the "
          "model (linear formula) specifying the relationship between the "
          "variables, can be specified once input data file has been read. "
          "\n  You can use simply X_1, ..., X_p for the title names in "
          "the first row, if you don't want to explicitly name them; "
          "but the first row MUST be the 'TITLE' row (i.e. it must not "
          "contain the actual data). The next n rows should be filled "
          "with the appropriate data."
          "\n  The input/output file names cannot contain spaces.\n";
}

bool ReadTableWithHeader::GetNominalColumns(
    const vector<string>& variable_names,
    vector<vector<DataHolder> >& data_values,
    map<int, set<string> >* nominal_columns, char* error_msg) {
  // Sanity Check input.
  if (variable_names.empty() || nominal_columns == nullptr) {
    sprintf(error_msg, "Invalid API call to GetNominalColumns. Aborting.");
    return false;
  }

  // Get the columns that have been inferred to be NOMINAL (based on
  // presence of '$' or a non-double value).
  string nominal_columns_detected;
  if (!nominal_columns->empty()) {
    string nominal_column_names;
    for (map<int, set<string> >::const_iterator itr = nominal_columns->begin();
         itr != nominal_columns->end(); ++itr) {
      if (itr->first >= variable_names.size()) {
        sprintf(error_msg, "Error in processing header line: column index is "
                           "out of bounds. Aborting.");
        return false;
      }
      if (!nominal_column_names.empty()) {
        nominal_column_names += ", ";
      }
      nominal_column_names += variable_names[itr->first];
    }
    nominal_columns_detected =
        "\n############################################################"
        "\nNominal Columns detected:\n" + nominal_column_names + "\n";
  } else {
    nominal_columns_detected =
        "\n############################################################"
        "\nNo nominal columns detected.\n";
  }

  // Display inferred NOMINAL columns to user.
  const string nominal_col_input_request = 
       "\nTo accept the above listed columns as the nominal ones, press "
       "<ENTER>.\nOtherwise, enter (comma-seperated) nominal columns, "
       "or enter 'NONE' if all variables are numeric\n(column "
       "name must exactly match input file, and is case-sensitive):\n> ";
  cout << nominal_columns_detected << nominal_col_input_request;

  // Prompt user for actual NOMINAL columns.
  set<int> new_nominal_columns;
  bool valid_input = false;
  string user_response;
  char temp;
  while (!valid_input) {
    user_response = "";
    cin.get(temp);
    while (temp != '\n') {
      user_response += temp;
      cin.get(temp);
    }
    char input_error[512] = "";
    if (user_response.empty()) {
      // User has accepted the detected nominal columns.
      valid_input = true;
    } else if (ReadInput::ProcessUserEnteredNominalColumns(
          user_response, variable_names, &new_nominal_columns, input_error)) {
      // Update set of nominal columns based on user input.
      if (!UpdateNominalColumns(
              data_values, new_nominal_columns, nominal_columns, input_error)) {
        cout << "\nUnable to honor your selection for NOMINAL variables: "
             << input_error << "\n"
             << nominal_columns_detected << nominal_col_input_request; 
      } else {
        valid_input = true;
      }
    } else {
      cout << "\nUnable to process your response: " << input_error << "\n"
           << nominal_columns_detected << nominal_col_input_request; 
    }
  }

  return true;
}

bool ReadTableWithHeader::ReadLine(
    const string& line, const int expected_num_columns, const int line_num,
    vector<DataHolder>* values, map<int, set<string> >* nominal_columns,
    set<int>* na_cols_for_line, char* error_msg) {
  if (values == nullptr || nominal_columns == nullptr) {
    sprintf(error_msg, "NULL 'values' or 'nominal_columns'; Check API call "
                       "to ReadLine(). Aborting.");
    return false;
  }
  if (line.empty()) {
    return true;
  }

  if (expected_num_columns >= 0) {
    values->resize(expected_num_columns);
  }

  // Iterate through line, stopping at white space.
  string token = "";
  int current_col_index = 0;
  for (string::const_iterator itr = line.begin(); itr != line.end(); ++itr) {
    // Append current non-whitespace character, and proceed to next character.
    if (*itr != ' ' && *itr != '\t') {
      token += *itr;
      continue;
    }

    // Sanity check line doesn't have too many entries or end in extra
    // whitespace.
    if (expected_num_columns >= 0 &&
        current_col_index > expected_num_columns - 1) {
      sprintf(error_msg,
              "Unable to parse line %d. Line either line ends in extra "
              "whitespace (space or tab), or too many data values).\n",
              line_num);
      return false;
    }

    // Nothing to do if token is empty (can happen e.g. with consecutive
    // whitespace characters in input file).
    if (token.empty()) {
      continue;
    }

    // Current char is whitespace (separator). Add token to 'values'.
    if (!AddToken(token, current_col_index, line_num,
                  values, nominal_columns, error_msg)) {
      return false;
    }
    
    // Add column to na_cols_for_line if necessary.
    if (token == "NA" || token == "N/A") {
      na_cols_for_line->insert(current_col_index);
    }

    current_col_index++;
    token = "";
  }

  // Reached end of line. Make sure line had proper number of entries.
  if (token.empty()) {
    sprintf(error_msg, "Unable to parse line %d:  Line ends in "
                       "extraneous whitespace (space or tab).\n", line_num);
    return false;
  }
  if (expected_num_columns >=0 &&
      current_col_index != expected_num_columns - 1) {
    sprintf(error_msg, "Unable to parse line %d: Wrong number of entries"
            ": expected %d entries (based on TITLE line), but found %d.\n",
            line_num, expected_num_columns, current_col_index + 1);
    return false;
  }
  
  // Store final token on the line.
  if (!AddToken(token, current_col_index, line_num,
                values, nominal_columns, error_msg)) {
    return false;
  }

  // Add column to na_cols_for_line if necessary.
  if (token == "NA" || token == "N/A") {
    na_cols_for_line->insert(current_col_index);
  }

  return true;
}

// Reads all data values (from file) for input and output variables;
// stores in indep_vars and dep_vars, respectively. Returns true if
// file is successfully parsed; returns false otherwise.
bool ReadTableWithHeader::ReadInputData(
    ifstream& file, const int expected_num_columns,
    vector<vector<DataHolder> >* values,
    map<int, set<string> >* nominal_columns,
    map<int, set<int> >* na_columns,
    char* error_msg) {
  // Sanity check input.
  if (values == nullptr || nominal_columns == nullptr) {
    sprintf(error_msg, "Null 'values' and/or 'nominal_columns'; "
                       "Check API call to ReadInputData(). Aborting.");
    return false;
  }

  values->clear();
  string line;
  int line_num = 2;
  int actual_num_columns = expected_num_columns;
  while (getline(file, line)) {
    vector<DataHolder> sample_values;
    set<int> na_columns_for_line;
    if (!ReadLine(line, actual_num_columns, line_num,
                  &sample_values, nominal_columns, &na_columns_for_line,
                  error_msg)) {
      return false;
    }
    if (!na_columns_for_line.empty()) {
      na_columns->insert(make_pair(values->size(), na_columns_for_line));
    }
    actual_num_columns = sample_values.size();
    ++line_num;
    values->push_back(sample_values);
  }
  return true;
}

bool ReadTableWithHeader::GetTitles(
    const string& title_line,
    vector<string>* titles, map<int, set<string> >* nominal_columns) {
  // Iterate through characters on line, using whitespace as a delimeter.
  // Treat consecutive non-whitespace characters as the 'token' representing
  // the title for that column.
  titles->clear();
  string temp_title = "";
  int column_num = 0;
  for (string::const_iterator itr = title_line.begin(); itr != title_line.end(); ++itr) {
    if (*itr != ' ' && *itr != '\t') {
      // Non-whitespace character, append it to current string of consecutive
      // characters (representing the current token).
      temp_title += *itr;
    } else if (!temp_title.empty()) {
      // Whitespace character; add current token to 'titles'.
      titles->push_back(temp_title);
      // Check whether the current title has a '$' suffix
      // (indicates a NOMINAL variable).
      if (nominal_columns != nullptr &&
          temp_title.substr(temp_title.length() - 1) == "$") {
        nominal_columns->insert(make_pair(column_num, set<string>()));
      }
      temp_title = "";
      column_num++;
    }
  }

  // Store final word on the first line.
  if (!temp_title.empty()) {
    titles->push_back(temp_title);
    // Check whether the current title has a '$' suffix
    // (indicates a NOMINAL variable).
    if (nominal_columns != nullptr &&
        temp_title.substr(temp_title.length() - 1) == "$") {
      nominal_columns->insert(make_pair(column_num, set<string>()));
    }
  }

  // Sanity check that titles are distinct.
  map<string, int> names;
  for (int i = 0; i < titles->size(); ++i) {
    if (names.find((*titles)[i]) != names.end()) return false;
    names.insert(make_pair((*titles)[i], i));
  }

  return true;
}

bool ReadTableWithHeader::ParseStrataAndSubgroup(
      const string& strata_str, const string& subgroup_str,
      const vector<string>& header,
      string* subgroup_rhs,
      set<int>* strata_cols, vector<int>* subgroup_cols) {
  // Get indices of the columns indicated by strata_str.
  if (!ReadInput::GetColumnIndicesFromString(
          strata_str, header, strata_cols)) {
    return false;
  }

  // Get indices of the columns indicated by subgroup_str.
  if (!ReadInput::GetSubgroupColumns(
          subgroup_str, header, subgroup_rhs, subgroup_cols)) {
      return false;
  }
  return true;
}

bool ReadTableWithHeader::UpdateNominalColumns(
    vector<vector<DataHolder> >& data_values,
    const set<int>& new_columns,
    map<int, set<string> >* columns,
    char* error_msg) {
  // Go through 'columns', checking for indices that are not present in
  // 'new_columns'. If any are found, make sure that all the values in
  // the corresponding set can be treated as a double: if so, remove the
  // (index, set) from 'columns'; otherwise, log an error and return false.
  set<int> columns_to_delete;
  for (map<int, set<string> >::const_iterator old_itr = columns->begin();
       old_itr != columns->end(); ++old_itr) {
    if (new_columns.find(old_itr->first) == new_columns.end()) {
      // Old index does not appear in new_columns. Make sure all values can
      // be represented as a double, and update data_values accordingly.
      for (set<string>::const_iterator value_itr = (old_itr->second).begin();
           value_itr != (old_itr->second).end(); ++value_itr) {
        double temp;
        if (!StringUtils::Stod(*value_itr, &temp)) {
          sprintf(error_msg, "Column '%d' was not entered as a NOMINAL column, "
                             "but value '%s' (which is not parsable as a "
                             "numerical value) appears in one of the data rows "
                             "for this column.",
                  (1 + old_itr->first), value_itr->c_str());
          return false;
        }
      }
      // All values in the column are parsable as a double. Update data_values
      // to take the double versions of all values in this column.
      for (vector<vector<DataHolder> >::iterator data_itr = data_values.begin();
           data_itr != data_values.end(); ++data_itr) {
        DataHolder& holder = (*data_itr)[old_itr->first];
        if (holder.name_.empty() ||
            !StringUtils::Stod(holder.name_, &(holder.value_))) {
          sprintf(error_msg, "Unknown error: unable to interpret '%s' from "
                             "Column '%d' as a numerical value.",
                  holder.name_.c_str(), (1 + old_itr->first));
          return false;
        }
        holder.name_ = "";
      }
      columns_to_delete.insert(old_itr->first);
    }
  }

  // Remove extra columns from 'old_columns' (could've done this in above loop,
  // but didn't want to mess up the iteration through 'old_columns' be deleting
  // elements while iterating through it).
  for (set<int>::const_iterator itr = columns_to_delete.begin();
       itr != columns_to_delete.end(); ++itr) {
    columns->erase(*itr);
  }

  // As of now, 'old_columns' is a subset of 'columns'. If they are equal, we're done,
  // and can return.
  if (columns->size() == new_columns.size()) return true;

  // 'old_columns' is a proper subset of 'columns'. For all the new columns, we
  // have to populate 'old_columns' with this index, together with a set that
  // represents all of the distinct values in that column.
  // First, get a list of all the new columns to add.
  set<int> new_columns_to_add;
  for (set<int>::const_iterator new_itr = new_columns.begin();
       new_itr != new_columns.end(); ++new_itr) {
    if (columns->find(*new_itr) == columns->end()) {
      new_columns_to_add.insert(*new_itr);
    }
  }

  // For each column in 'new_columns_to_add', create a set of all the distinct
  // values in this column, and populate 'columns' with it and column index.
  for (vector<vector<DataHolder> >::iterator data_itr = data_values.begin();
       data_itr != data_values.end(); ++ data_itr) {
    for (set<int>::const_iterator itr = new_columns_to_add.begin();
         itr != new_columns_to_add.end(); ++itr) {
      // Sanity check, should never be true.
      if (*itr >= data_itr->size()) {
        sprintf(error_msg, "Unknown error: column %d is outside of the range "
                           "of input variables (%d).", *itr, data_itr->size());
        return false;
      }

      // Get data value for this (row, column).
      DataHolder& holder = (*data_itr)[*itr];
      if (!holder.name_.empty()) {
        sprintf(error_msg, "Unknown error: column %d has a stored data value "
                           "that is a string (%s), even though column was "
                           "hitherto considered a non-NOMINAL column.",
                           *itr, holder.name_.c_str());
        return false;
      }

      // Update data value to be a string (instead of a numerical value).
      stringstream ss;
      ss << holder.value_;
      holder.name_ = ss.str();
      holder.value_ = 0.0;

      // Update 'columns' by adding the new string value to the set.
      map<int, set<string> >::iterator old_itr = columns->find(*itr);
      if (old_itr == columns->end()) {
        // This is the first value being added for this column.
        set<string> value_to_add;
        value_to_add.insert(holder.name_);
        columns->insert(make_pair(*itr, value_to_add));
      } else {
        // This column was already known to be nominal. Add the new value to the
        // set of values for this column.
        set<string>& existing_values = old_itr->second;
        existing_values.insert(holder.name_);
      }
    }
  }

  return true;
}

/*
// Prints Variable titles together with the estimated regression coefficients
// and covariance; output is print to file indicated by 'outfile'.
bool ReadTableWithHeader::Print(
    const string& outfile,
    const vector<string>& titles,
    const double& variance,
    const MatrixXd& inverse_of_indep_vars,
    const MatrixXd& regression_coefficients,
    char* error_msg) {
  ofstream out_file;
  out_file.open(outfile.c_str());
  // Write title line.
  out_file << "Variable_Name\t\tEstimate\tVariance\tS.E.E.  \tT-Statistic\tp-value\n";
  char var_name[512] = "";
  for (int i = 0; i < regression_coefficients.size(); ++i) {
    if (i == 0) {
      sprintf(var_name, "Constant (B_0)");
    } else {
      sprintf(var_name, "%s (B_%d)", titles[i].c_str(), i);
    }
    double est_val = regression_coefficients(i, 0);
    double est_var = variance * inverse_of_indep_vars(i, i);
    double see = -1.0;
    if (est_var >= 0.0) {
      see = sqrt(est_var);
    } else {
      sprintf(error_msg, "Negative estimated variance (should never happen): %.04f\n", est_var);
    }
    double test_stat = -1.0;
    if (see > 0.0) {
      test_stat = est_val / see;
    } else {
      sprintf(error_msg, "Negative S.E.E. (should never happen). est_var: %.04f,"
              "S.E.E.: %.04f.\n", est_var, see);
    }
    double p_value = -1.0;
    p_value = ReverseIncompleteGammaFunction(0.5, (0.5 * test_stat * test_stat));
    char out_line[512] = "";
    sprintf(out_line, "%s\t\t%0.06f\t%0.06f\t%0.06f\t%0.06f\t%0.06f\n",
            var_name, est_val, est_var, see, test_stat, p_value);
    out_file << out_line;
  }
  out_file.close();
  return true;
}
*/
bool ReadTableWithHeader::ReadTitleLine(
    const string& filename, const string& working_dir,
    vector<string>* titles, map<int, set<string> >* nominal_columns,
    char* error_msg) {
  // Open input file.
  ifstream input_file(filename.c_str());
  if (!input_file.is_open()) {
    const string dir = working_dir.empty() ? "" :
                       (": (" + working_dir + ").");
    sprintf(error_msg, "Unable to find file. Make sure that it "
                       "is present in your current directory%s.",
                       dir.c_str());
    return false;
  }

  // Read Title Line.
  string title_line;
  if (!getline(input_file, title_line)) {
    sprintf(error_msg, "Error: Empty file.");
    return false;
  }
  if (!GetTitles(title_line, titles, nominal_columns)) {
    sprintf(error_msg, "Improper Title Line.\n%s",
                       PrintInputFormat().c_str());
    return false;
  }

  // Sanity check that titles are distinct.
  map<string, int> names;
  for (int i = 0; i < titles->size(); ++i) {
    if (names.find((*titles)[i]) != names.end()) {
      sprintf(error_msg, "'%s' appears as the variable name for two different"
                         "columns (%d and %d). Aborting.",
              (*titles)[i].c_str(), (1 + names[(*titles)[i]]), (1 + i));
      return false;
    }
    names.insert(make_pair((*titles)[i], i));
  }

  return true;
}

bool ReadTableWithHeader::ReadFile(
    const string& filename, const string& working_dir,
    vector<string>* header, vector<vector<DataHolder> >* data_values,
    map<int, set<string> >* nominal_columns, map<int, set<int> >* na_columns,
    char* error_msg) {
  if (header == nullptr || data_values == nullptr) {
    sprintf(error_msg, "Null header and/or data_values; check API for call "
                       "to ReadFile(). Aborting.");
    return false;
  }

  // Open input file.
  ifstream input_file(filename.c_str());
  if (!input_file.is_open()) {
    const string dir = working_dir.empty() ? "" :
                       (": (" + working_dir + ").");
    sprintf(error_msg, "Unable to find file. Make sure that it "
                       "is present in your current directory: %s.",
                       dir.c_str());
    return false;
  }
  
  // Read Title line of input file.
  string title_line;
  if (!getline(input_file, title_line)) {
    sprintf(error_msg, "Error: Empty input file.");
    return false;
  }
  if (!GetTitles(title_line, header, nominal_columns)) {
    sprintf(error_msg, "Improper file format.\n%s",
                       PrintInputFormat().c_str());
    return false;
  }

  // Read the rest (actual data) of input file.
  char new_error_msg[512] = "";
  if (!ReadInputData(
          input_file, header->size(),
          data_values, nominal_columns, na_columns, new_error_msg)) {
    string new_error_msg_str(new_error_msg);
    sprintf(error_msg, "Improper file format %s\n%s",
            new_error_msg_str.c_str(), PrintInputFormat().c_str());
    return false;
  }
  return true;
}

}  // namespace file_reader_utils
