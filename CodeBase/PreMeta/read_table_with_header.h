#include "FileReaderUtils/read_input.h"
#include "MathUtils/data_structures.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifndef READ_TABLE_WITH_HEADER_H
#define READ_TABLE_WITH_HEADER_H

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using math_utils::CensoringData;
using math_utils::DataHolder;
using math_utils::LinearTerm;
using math_utils::VariableTerm;

namespace file_reader_utils {

class ReadTableWithHeader {
 public:
  // Returns a string describing the expected format of the
  // input file to be read.
  static string PrintInputFormat();

  // Parses Strata and Subgroup.
  static bool ParseStrataAndSubgroup(
      const string& strata_str, const string& subgroup_str,
      const vector<string>& header,
      string* subgroup_rhs,
      set<int>* strata_cols, vector<int>* subgroup_cols);

  // Reads in filename, prompts user to enter model, and generates the
  // dependent/independent variables used in that model.
  // Populates dep_var and indep_vars with these values, and populates
  // legend with all of the linear terms (for the independent variables)
  // in the model. Returns true if successful, false otherwise (in which
  // case an error message containing the error is printed).
  template<typename T>
  // typename T (dep_var) should either be:
  //   1) vector<double> (for linear, logistic regression models)
  //   2) vector<CensoringData> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool ReadFileAndGetModel(
      const bool include_error_term, const string& filename, const string& working_dir,
      vector<string>* legend, T* dep_var, vector<VectorXd>* indep_vars,
      string* final_model) {
    // Sanity check input.
    if (legend == nullptr || dep_var == nullptr || indep_vars == nullptr) {
      cout << "\nERROR in ReadFileAndGetModel: null input parameters." << endl;
      return false;
    }

    // Open and read input file.
    vector<string> header;
    vector<vector<DataHolder> > data_values;
    map<int, set<string> > nominal_columns;
    map<int, set<int> > na_columns;
    char error_msg[1024] = "";
    if (!ReadFile(
            filename, working_dir,
            &header, &data_values, &nominal_columns, &na_columns, error_msg)) {
      string error_msg_str(error_msg);
      cout << "\nERROR: Unable to ReadFile (" << filename << "): "
           << error_msg_str << "\n";
      return false;
    }

    // Print the detected nominal columns, and get input from user for
    // the actual nominal columns.
    if (!ReadTableWithHeader::GetNominalColumns(
            header, data_values, &nominal_columns, error_msg)) {
      string error_msg_str(error_msg);
      cout << "\nERROR in getting nominal columns: " << error_msg_str << "\n";
      return false;
    }

    return CallGetAndParseModel(
        include_error_term, header, data_values, nominal_columns, na_columns,
        legend, dep_var, indep_vars, final_model);
  }

  // Same as above, passing in empty string for working_dir, and true for
  // 'include_error_term'.
  template<typename T>
  // typename T (dep_var) should either be:
  //   1) vector<double> (for linear, logistic regression models)
  //   2) vector<CensoringData> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool ReadFileAndGetModel(
      const string& filename,
      vector<string>* legend, T* dep_var, vector<VectorXd>* indep_vars,
      string* final_model) {
    return ReadFileAndGetModel(
        true, filename, "", legend, dep_var, indep_vars, final_model);
  }

  // Same as above, with slightly different API. Main differences:
  //   - Use NOMINAL variables detected (via '$' in the title header, or
  //     based on input data that is not parsable as a numeric value), and
  //     don't prompt user to clarify/update
  //   - Instead of prompting user for a Model, it is input as a string
  //   - There is an ofstream to log errors/comments to
  //   - Populate Eigen Matrices and Vector instead of (C++) vectors
  // Note that ofstream should be opened by caller, and closed by caller
  // upon function return.
  template<typename T>
  // typename T (dep_var) should either be:
  //   1) VectorXd (for linear, logistic regression models)
  //   2) vector<CensoringData> (for Cox Proportional Hazards Model)
  static bool ReadFileAndParseModel(
      const bool include_error_term, const bool subgroup_as_cov,
      const string& filename, const string& model,
      const string& strata_str, const string& subgroup_str,
      ofstream& log_file,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend,
      T* dep_var, MatrixXd* indep_vars, string* final_model) {
    // Sanity check input.
    if (legend == nullptr || dep_var == nullptr || indep_vars == nullptr) {
      log_file << "\nERROR in ReadFileAndParseModel: null input parameters.\n";
      return false;
    }

    // Open and read input file.
    vector<string> header;
    vector<vector<DataHolder> > data_values;
    map<int, set<string> > nominal_columns;
    map<int, set<int> > na_columns;
    char error_msg[1024] = "";
    if (!ReadFile(
            filename,
            &header, &data_values, &nominal_columns, &na_columns, error_msg)) {
      string error_msg_str(error_msg);
      log_file << "\nERROR in ReadFileAndParseModel: Unable to ReadFile ("
               << filename << "): " << error_msg_str << "\n";
      return false;
    }

    set<int> strata_cols;
    vector<int> subgroup_cols;
    string subgroup_rhs;
    if (!ParseStrataAndSubgroup(
          strata_str, subgroup_str, header,
          &subgroup_rhs, &strata_cols, &subgroup_cols)) {
      log_file << "ERROR in ParseStrataAndSubgroup." << endl;
      return false;
    }

    // Parse Model.
    return CallParseInputData(
        include_error_term, subgroup_as_cov,
        model, subgroup_rhs, log_file, header, data_values,
        strata_cols, subgroup_cols, nominal_columns, na_columns,
        strata, subgroup_sizes, legend, dep_var, indep_vars, final_model);
  }
  // Same as above, passing in 'true' for subgroup_as_cov.
  template<typename T>
  static bool ReadFileAndParseModel(
      const bool include_error_term,
      const string& filename, const string& model,
      const string& strata_str, const string& subgroup_str,
      ofstream& log_file,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend, T* dep_var, MatrixXd* indep_vars, string* final_model) {
    return ReadFileAndParseModel<T>(
        include_error_term, !subgroup_str.empty(),
        filename, model, strata_str, subgroup_str, log_file, strata,
        subgroup_sizes, legend, dep_var, indep_vars, final_model);
  }
  // Same as above, but passes in empty subgroup_str and null subgroup_sizes.
  template<typename T>
  static bool ReadFileAndParseModel(
      const bool include_error_term,
      const string& filename, const string& model,
      const string& strata_str, ofstream& log_file,
      vector<pair<int, int>>* strata, vector<string>* legend,
      T* dep_var, MatrixXd* indep_vars, string* final_model) {
    return ReadFileAndParseModel<T>(
        include_error_term, filename, model, strata_str, "" /* subgroup_str */,
        log_file, strata, nullptr, legend, dep_var, indep_vars, final_model);
  }
  // Same as above, but passes in empty strata_cols and null strata
  // (these will be ignored).
  template<typename T>
  static bool ReadFileAndParseModel(
      const bool include_error_term, const string& filename,
      const string& model, ofstream& log_file,
      vector<string>* legend, T* dep_var, MatrixXd* indep_vars,
      string* final_model) {
    return ReadFileAndParseModel<T>(
        include_error_term, filename, model, "" /* strata_str */, log_file,
        nullptr /* strata */, legend, dep_var, indep_vars, final_model);
  }

  // Reads in the form of the model and the raw data, and constructs the
  // values for dep_var and indep_vars.
  template<typename T, typename U>
  // typename T (dep_var) should either be:
  //   1) VectorXd              (for linear, logistic regression models)
  //   2) vector<CensoringData> (for Cox Proportional Hazards Model)
  // typename U should either be:
  //   1) VariableTerm          (for Linear/Logistic Regression)
  //   2) vector<VariableTerm>  (for Cox PH Model)
  static bool ParseInputData(
      const bool include_error_term,
      const bool subgroup_as_cov,
      const string& model,
      const string& subgroup_rhs,
      ofstream& log_file,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const set<int>& strata_cols,
      const vector<int>& subgroup_cols,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend,
      T* dep_var, MatrixXd* indep_vars, string* final_model) {
    // Read user-entered model, linking the names to those in the
    // header row of the input file, and interpretting the RHS of
    // the user entered model into concrete variables and operations.
    set<int> input_cols_used;
    U dependent_var;
    vector<LinearTerm> linear_vars;
    if (!ParseModel<U>(model, log_file, header, nominal_columns,
                       &input_cols_used, &dependent_var, &linear_vars)) {
      log_file << "\nERROR in ParseInputData: unable to ParseModel.\n";
      return false;
    }
    // Add strata_cols and subgroup_cols to input_cols_used.
    for (const int subgroup_col : subgroup_cols) {
      input_cols_used.insert(subgroup_col);
    }
    for (const int strata_col : strata_cols) {
      input_cols_used.insert(strata_col);
    }

    // Read in the input data, using the parsed model to compute each
    // (in)dependent term.
    char error_msg[512] = "";
    if (!ReadInput::ParseModel<T, MatrixXd, U>(
            include_error_term, subgroup_as_cov, subgroup_rhs, header,
            strata_cols, subgroup_cols, nominal_columns, data_values,
            dependent_var, linear_vars, na_columns, input_cols_used,
            dep_var, indep_vars, strata, subgroup_sizes,
            legend, final_model, error_msg)) {
      string error_msg_str(error_msg);
      log_file << "\nERROR in parsing model: " << error_msg_str << "\n";
      return false;
    }

    return true;
  }

 private:
  // Read the first line of 'file', placing each read item into 'titles'.
  // Returns true if file was successfully parsed, false otherwise.
  static bool GetTitles(const string& title_line,
                        vector<string>* titles,
                        map<int, set<string> >* nominal_columns);
  
  // Prompts user for desired nominal columns, displaying first the columns that
  // were already detected as nominal (based on Title line having a '$' next
  // to that variable name, and/or non-double values in the column). Updates
  // nominal_columns with the actual columns the user wants to be nominal.
  static bool GetNominalColumns(
      const vector<string>& variable_names,
      vector<vector<DataHolder> >& data_values,
      map<int, set<string> >* nominal_columns, char* error_msg);

  // If expected_num_columns >=0, then this function sets 'values' to be a
  // vector of size 'expected_num_columns', and then attempts to read this many
  // values from the given 'line'. Otherwise, reads however many items exist
  // in the line into 'values'. Also populates 'nominal_columns' with the
  // columns (and all distinct values that appeared in that column) that were
  // indicated as 'nominal' via the header file (via '$') or because there was
  // a value in that column that could not be parsed as a double.
  // Also populates 'na_cols_for_line' with the indices of columns for which
  // this line indicates 'NA'.
  // Returns true if successful, otherwise returns false (and 'error_msg'
  // contains details of the failure).
  static bool ReadLine(
      const string& line, const int expected_num_columns, const int line_num,
      vector<DataHolder>* values, map<int, set<string> >* nominal_columns,
      set<int>* na_cols_for_line, char* error_msg);

  // Parses the first line of working_dir/filename into 'titles'.
  // Also, if nominal_columns is non-null, adds (col_index, empty set) to
  // nominal_columns when that column represents a NOMINAL variable (as
  // determined by '$' suffix in the column title).
  // Returns true if successful, otherwise returns false (and 'error_msg'
  // contains details of the failure).
  static bool ReadTitleLine(
      const string& filename, const string& working_dir,
      vector<string>* titles, map<int, set<string> >* nominal_columns,
      char* error_msg);

  // Calls ReadTitleLine above with NULL for nominal_columns.
  static bool ReadTitleLine(
      const string& filename, const string& working_dir,
      vector<string>* titles, char* error_msg) {
    return ReadTitleLine(filename, working_dir, titles, NULL, error_msg);
  }

  // Calls ReadTitleLine above with an empty working dir.
  static bool ReadTitleLine(
      const string& filename,
      vector<string>* titles, map<int, set<string> >* nominal_columns,
      char* error_msg) {
    return ReadTitleLine(filename, "", titles, nominal_columns, error_msg);
  }
  
  // Calls ReadTitleLine above with an empty working dir and
  // NULL for nominal_columns.
  static bool ReadTitleLine(
      const string& filename, vector<string>* titles, char* error_msg) {
    return ReadTitleLine(filename, "", titles, NULL, error_msg);
  }

  // Parses the data in 'file' into 'data_values'. More specifically, 
  // 'file' is all the data from an input file (minus the title line); each
  // row is expected to have 'expected_num_columns' (to match title line);
  // unless 'expected_num_columns' is -1, in which case this function only
  // demands that all rows have the same number of columns.
  // Also populates 'nominal_columns' with the columns that (had at least
  // one row that) could not be parsed as a double.
  // Returns true if file was successfully parsed, false otherwise (in which
  // case 'error_msg' will contain info about why the parsing failed).
  static bool ReadInputData(
      ifstream& file, const int expected_num_columns,
      vector<vector<DataHolder> >* data_values,
      map<int, set<string> >* nominal_columns,
      map<int, set<int> >* na_columns,
      char* error_msg);
  
  // Same as above, but passes in 'expected_num_columns' = -1.
  static bool ReadInputData(
      ifstream& file,
      vector<vector<DataHolder> >* data_values,
      map<int, set<string> >* nominal_columns,
      map<int, set<int> >* na_columns,
      char* error_msg) {
    return ReadInputData(file, -1,
                         data_values, nominal_columns, na_columns, error_msg);
  }

  // Parses the first line of working_dir/filename into 'header', and parses the
  // rest of the file into data_values. Populates 'nominal_columns' with the
  // columns that (had at least one row that) could not be parsed as a double,
  // and populates 'na_columns' such that the map is keyed by the input data row,
  // and the corresponding set is all indices for which that row indicated 'NA'.
  // Returns true if successful, otherwise returns false (and 'error_msg'
  // contains details of the failure).
  static bool ReadFile(
      const string& filename, const string& working_dir,
      vector<string>* header, vector<vector<DataHolder> >* data_values,
      map<int, set<string> >* nominal_columns, map<int, set<int> >* na_columns,
      char* error_msg);

  // Calls ReadFile above with an empty working dir.
  static bool ReadFile(
      const string& filename,
      vector<string>* header, vector<vector<DataHolder> >* data_values,
      map<int, set<string> >* nominal_columns, map<int, set<int> >* na_columns,
      char* error_msg) {
    return
        ReadFile(filename, "",
                 header, data_values, nominal_columns, na_columns, error_msg);
  }

  // Updates 'columns' to reflect the indices in 'new_columns' by:
  //  (A) For indices in 'columns' that are NOT present in 'new_columns',
  //      check to make sure all data values in the corresponding set can
  //      be parsed as a double; otherwise return false
  //  (B) For indices in 'new_columns' that are not in 'columns', add the
  //      index to 'columns', with the corresponding set being all the
  //      (distinct) data values that appear among the data values for
  //      that column
  //  Returns true as long as (A) above does not fail.
  static bool UpdateNominalColumns(
      vector<vector<DataHolder> >& data_values,
      const set<int>& new_columns,
      map<int, set<string> >* columns,
      char* error_msg);

  // Populates dep_var and indep_vars with values from data_values and populates
  // legend with all of the linear terms (for the independent variables)
  // in the model. Returns true if successful, false otherwise (in which
  // case an error message containing the error is printed).
  template<typename U>
  // typename U should either be:
  //   1) VariableTerm          (for Linear/Logistic Regression)
  //   2) vector<VariableTerm>  (for Cox PH Model)
  static bool ParseModel(
      const string& model, ofstream& log_file,
      const vector<string>& header,
      const map<int, set<string> >& nominal_columns,
      set<int>* input_cols_used,
      U* dependent_var, vector<LinearTerm>* linear_vars) {
    // Parse passed-in model.
    char error_msg[1024] = "";
    if (!ReadInput::ProcessUserEnteredModel<U>(
            model, dependent_var, linear_vars, error_msg)) {
      string error_msg_str(error_msg);
      log_file << "\nERROR: Unable to parse model:\n" << model << endl
               << "Error message:\n" << error_msg_str << endl;
      return false;
    }

    // Sanity Check Model.
    map<string, int> titles;
    for (int i = 0; i < header.size(); ++i) {
      titles.insert(make_pair(header[i], i));
    }
    if (!ReadInput::SanityCheckModel<U>(
            titles, nominal_columns,
            dependent_var, linear_vars, input_cols_used, error_msg)) {
      string error_msg_str(error_msg);
      log_file << "\nERROR: " << error_msg_str << "\n";
      return false;
    }
    return true;
  }

  // Dummy function to call ParseInputData with the proper template type:
  // use typename = VariableTerm for Linear/Logistic Regression, and typname =
  // vector<VariableTerm> for Cox PH Model. The below is for Linear/Logistic
  // Regression (since it's dep_var vector stores a single dependent variable
  // per sample; Cox PH Model requires multiple dependent variables for
  // Survival/Censoring Time and Status); the one below that is for Cox.
  static bool CallParseInputData(
      const bool include_error_term,
      const bool subgroup_as_cov,
      const string& model,
      const string& subgroup_rhs,
      ofstream& log_file,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const set<int>& strata_cols,
      const vector<int>& subgroup_cols,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend,
      VectorXd* dep_var, MatrixXd* indep_vars, string* final_model) {
    return ParseInputData<VectorXd, VariableTerm>(
        include_error_term, subgroup_as_cov,
        model, subgroup_rhs, log_file, header, data_values,
        strata_cols, subgroup_cols, nominal_columns, na_columns,
        strata, subgroup_sizes, legend, dep_var, indep_vars, final_model);
  }
  // Same as above, for the Cox API.
  static bool CallParseInputData(
      const bool include_error_term,
      const bool subgroup_as_cov,
      const string& model,
      const string& subgroup_rhs,
      ofstream& log_file,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const set<int>& strata_cols,
      const vector<int>& subgroup_cols,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend,
      vector<CensoringData>* dep_var, MatrixXd* indep_vars, string* final_model) {
    return ParseInputData<vector<CensoringData>, vector<VariableTerm>>(
        include_error_term, subgroup_as_cov,
        model, subgroup_rhs, log_file, header, data_values,
        strata_cols, subgroup_cols, nominal_columns, na_columns,
        strata, subgroup_sizes, legend, dep_var, indep_vars, final_model);
  }

  // Prompts user for the model, and populates legend, dep_var, and indep_vars.
  template<typename T, typename U>
  // typename T (dep_var) should either be:
  //   1) VectorXd              (for linear, logistic regression models)
  //   2) vector<CensoringData> (for Cox Proportional Hazards Model)
  // typename U should either be:
  //   1) VariableTerm          (for Linear/Logistic Regression)
  //   2) vector<VariableTerm>  (for Cox PH Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool GetAndParseModel(
      const bool include_error_term,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<string>* legend,
      T* dep_var, vector<VectorXd>* indep_vars, string* final_model) {
    // Get Model Specifications from User.
    U dependent_var;
    vector<LinearTerm> linear_vars;
    set<int> input_cols_used;
    char error_msg[512] = "";
    if (!ReadInput::GetModel<U>(
            header, nominal_columns,
            &dependent_var, &linear_vars, &input_cols_used, error_msg)) {
      string error_msg_str(error_msg);
      cout << "\nERROR in getting model: " << error_msg_str << "\n";
      return false;
    }

    // Parse model to create vectors for dependent variable and each linear term
    // in the model.
    if (!ReadInput::ParseModel<T, vector<VectorXd>, U>(
            include_error_term, header, nominal_columns, data_values,
            dependent_var, linear_vars, na_columns, input_cols_used,
            dep_var, indep_vars, legend, final_model, error_msg)) {
      string error_msg_str(error_msg);
      cout << "\nERROR in parsing model: " << error_msg_str << "\n";
      return false;
    }

    return true; 
  }

  // Dummy function to call GetAndParseModel with the proper template type:
  // use typename = VariableTerm for Linear/Logistic Regression, and typname =
  // vector<VariableTerm> for Cox PH Model. The below is for Linear/Logistic
  // Regression (since it's dep_var vector stores a single dependent variable
  // per sample; Cox PH Model requires multiple dependent variables for
  // Survival/Censoring Time and Status); the one below that is for Cox.
  static bool CallGetAndParseModel(
      const bool include_error_term,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<string>* legend,
      vector<double>* dep_var, vector<VectorXd>* indep_vars, string* final_model) {
    return GetAndParseModel<vector<double>, VariableTerm>(
        include_error_term, header, data_values, nominal_columns, na_columns,
        legend, dep_var, indep_vars, final_model);
  }
  // Same as above, for the Cox API.
  static bool CallGetAndParseModel(
      const bool include_error_term,
      const vector<string>& header,
      const vector<vector<DataHolder> >& data_values,
      const map<int, set<string> >& nominal_columns,
      const map<int, set<int> >& na_columns,
      vector<string>* legend,
      vector<CensoringData>* dep_var, vector<VectorXd>* indep_vars,
      string* final_model) {
    return GetAndParseModel<vector<CensoringData>, vector<VariableTerm>>(
        include_error_term, header, data_values, nominal_columns, na_columns,
        legend, dep_var, indep_vars, final_model);
  }
};

}  // namespace file_reader_utils

#endif
