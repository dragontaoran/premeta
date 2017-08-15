#include "MathUtils/data_structures.h"
#include "StringUtils/string_utils.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifndef READ_INPUT_H
#define READ_INPUT_H

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using math_utils::CensoringData;
using math_utils::DataHolder;
using math_utils::LinearTerm;
using math_utils::Operation;
using math_utils::SamplingParameters;
using math_utils::VariableTerm;
using string_utils::StringUtils;

namespace {

static const char kSubgroupString[] = "Subgroup";

void CopyVariables(const vector<double>& in, vector<CensoringData>* out) {
  // This API needed only because Template doesn't know that this method
  // should never be called with this type for 'out'.
  // Do nothing, but record error.
  cout << "ERROR in CopyVariables: Unexpected call with CensoringData API."
       << endl;
}

void CopyVariables(const vector<double>& in, vector<double>* out) {
  if (in.empty() || out == nullptr) return;
  out->clear();
  for (const double& value : in) {
    out->push_back(value);
  }
}

void CopyVariables(const vector<double>& in, VectorXd* out) {
  if (in.empty() || out == nullptr) return;
  out->resize(in.size());
  for (int i = 0; i < in.size(); ++i) {
    (*out)[i] = in[i];
  }
}

void CopyVariables(const vector<VectorXd>& in, vector<VectorXd>* out) {
  if (in.empty() || out == nullptr) return;
  out->clear();
  for (const VectorXd& row : in) {
    out->push_back(row);
  }
}

void CopyVariables(const vector<VectorXd>& in, MatrixXd* out) {
  if (in.empty() || out == nullptr) return;
  const VectorXd& first_row = in[0];
  out->resize(in.size(), first_row.size());
  for (int i = 0; i < in.size(); ++i) {
    out->row(i) = in[i];
  }
}

}  // namespace

namespace file_reader_utils {

class ReadInput {
 public:
  // Returns true if value is any of:
  //   '1', '1.0', 'T', 'true', 'True', or 'TRUE'
  static bool IsTrueIndicator(const string& value);
  // Returns true if value is any of:
  //   '0', '0.0', 'F', 'false', 'False', or 'FALSE'
  static bool IsFalseIndicator(const string& value);

  static string GetErrorTermId();

  // Prompts user with input request 'query_text'. If user-entered response
  // is parsable as an int, return that int. Otherwise, display error to
  // user and re-prompt.
  static int GetIntegerFromUser(const string& query_text);
  // Ditto above, except tries to parse response as a double.
  static double GetDoubleFromUser(const string& query_text);
  // Asks user to enter the desired number of data rows to simulate ("n").
  static int GetNumDataRowsPerSimulationFromUser();
  // Asks user to enter the desired number of simulations to run ("k").
  static int GetNumSimulationsFromUser();

  // Attempts to parse strata by finding the referenced strings in 'header'
  // and populates 'strata_cols' with the corresponding indices.
  static bool GetColumnIndicesFromString(
      const string& strata, const vector<string>& header, set<int>* strata_cols);
  // Attempts to parse subgroup by finding the referenced strings in 'header'
  // and populates 'subgroup_cols' with the corresponding indices. This is same
  // as above, except we use a vector rather than a set, in case the original
  // order of terms in subgroup string are important (for subgroup they are,
  // since there is a RHS of the --subgroup parameter, which specifies values
  // to match, in the same order as the columns were listed on the LHS).
  static bool GetColumnIndicesFromString(
      const string& subgroup, const vector<string>& header,
      vector<int>* subgroup_cols);

  // Parses subgroup_str by splitting aroung the '=', parsing the LHS into
  // subgroup_cols (based on matching the terms that appear to their
  // corresponding index in 'header'), and putting the RHS in subgroup_rhs.
  static bool GetSubgroupColumns(
    const string& subgroup_str, const vector<string>& header,
    string* subgroup_rhs, vector<int>* subgroup_cols);

  // Parses the RHS of the user-entered --subgroup flag.
  // Example:
  //   RHS: {(WHITE, SQUARE), (WHITE, CIRCLE), (BROWN, SQUARE)}
  // Then subgroups would be a vector of length 3, and each inner-vector
  // would have length two:
  //  subgroups: [[White, Square], [White, Circle], [Brown, Square]]
  static bool ParseSubgroups(
      const bool subgroup_as_cov,
      const string& subgroup_rhs, const int num_expected_subgroups,
      vector<vector<string>>* subgroups);

  // Reads in a string, and populates SamplingParameters appropriately.
  // Returns true if string was successfully parsed.
  static bool ParseSamplingParameter(const string& value,
                                     SamplingParameters* params);

  // For each LinearTerm in linear_vars, asks user how each variable should
  // be sampled, and records user's responses in 'sampling_params'. Upon
  // return, sampling_params->size() should equal the number of variables
  // that appear across all linear terms plus '1' (an extra SamplingParameter
  // is needed to specify how to sample the error term; this special
  // SamplingParameter is stored under key 'kErrorTermId').
  static bool GetSamplingParametersFromUser(
      const vector<LinearTerm>& linear_vars,
      map<string, SamplingParameters>* sampling_params, char* error_msg);

  // Makes sure that the model specified by dependent_var is consistent:
  //   - The syntax is recognizable
  //   - The operations specified are supported (currently +, *, log, exp)
  //   - The NOMINAL variables never have log or exp applied to them
  // If so, populates input_cols_used with all the columns referred to by
  // dependent_var and returns true. Otherwise, returns false.
  static bool SanityCheckDependentVariable(
      const map<string, int>& titles,
      const map<int, set<string> >& nominal_columns,
      const VariableTerm* dependent_var,
      set<int>* input_cols_used,
      char* error_msg);
  // Same as above, with different API for dependent_var (represents Cox
  // Hazards model, so it has two LinearTerms).
  static bool SanityCheckDependentVariable(
      const map<string, int>& titles,
      const map<int, set<string> >& nominal_columns,
      const vector<VariableTerm>* dependent_var,
      set<int>* input_cols_used,
      char* error_msg);

  // Makes sure that the model specified by linear_vars is consistent:
  //   - The syntax is recognizable
  //   - The operations specified are supported (currently +, *, log, exp)
  //   - The NOMINAL variables never have log or exp applied to them
  // If so, populates input_cols_used with all the columns referred to by
  // linear_vars and returns true. Otherwise, returns false.
  static bool SanityCheckIndependentVariables(
      const map<string, int>& titles,
      const map<int, set<string> >& nominal_columns,
      const vector<LinearTerm>* linear_vars,
      set<int>* input_cols_used,
      char* error_msg);

  // Makes sure that the model specified by dependent_var and linear_vars
  // is consistent:
  //   - The syntax is recognizable
  //   - The operations specified are supported (currently +, *, log, exp)
  //   - The NOMINAL variables never have log or exp applied to them
  // If so, populates input_cols_used with all the columns referred to by
  // dependent_var and linear_vars, and returns true. Otherwise, returns false
  // and populates error_msg with the reason for failure.
  template<typename T>
  // typename T (dependent_var) should either be:
  //   1) VariableTerm (for linear, logistic regression models)
  //   2) vector<VariableTerm> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool SanityCheckModel(
      const map<string, int>& titles,
      const map<int, set<string> >& nominal_columns,
      const T* dependent_var, const vector<LinearTerm>* linear_vars,
      set<int>* input_cols_used, char* error_msg) {
    // Sanity check inputs.
    if (dependent_var == nullptr || linear_vars == nullptr ||
        input_cols_used == nullptr) {
      sprintf(error_msg, "Null pointer passed to SanityCheckModel. Check API.");
      return false;
    }

    if (!SanityCheckDependentVariable(
            titles, nominal_columns, dependent_var, input_cols_used, error_msg)) {
      return false;
    }    
    if (!SanityCheckIndependentVariables(
            titles, nominal_columns, linear_vars, input_cols_used, error_msg)) {
      return false;
    }
    return true;
  }

  // Returns appropriate text (strings) that will be used to prompt user for
  // appropriate input, depending on whether the Model is Linear/Logistic or
  // Cox Proportional Hazards (we distinguish among these cases based on
  // whether the dependent variable is a single VariableTerm (linear/logistic)
  // or a vector of VariableTerms (one of each of the time variables, and one
  // for the status ("Delta") variable; for Cox).
  static string PromptUserForModelText(const VariableTerm* dep_var) {
    return
        "\nEnter your linear Model."
         "\n\nYour expression can involve '+', '*', 'Log', 'exp', 'sqrt', and "
         "'pow' terms, and should have one equals sign with the dependent "
         "variable on the left and linear terms (involving the independent "
         "variables) on the right. Regression coefficients, explicit mention "
         "of a constant term (assumed implicitly), and the error term should "
         "not be specified. Nominal variables can optionally be specified by a "
         "'$' suffix (otherwise they will be automatically detected based on t"
         "he data in the input file). Variable names should not contain spaces";
  }
  static string PromptUserForModelText(vector<VariableTerm>* dep_var) {
    return
        "\nEnter your Cox Model."
         "\n\nYour expression can involve '+', '*', 'Log', 'exp', 'sqrt', and "
         "'pow' terms, and should have one equals sign with the dependent "
         "variables (e.g. time, status) on the left and linear terms "
         "(involving the independent variables) on the right. Regression "
         "coefficients and the error term should not be specified; and no "
         "constant term is assumed. Nominal variables can optionally be "
         "specified by a '$' suffix (otherwise they will be automatically "
         "detected based on the data in the input file). Variable names "
         "should not contain spaces";
  }
  static string ExampleModel(const VariableTerm* dep_var) {
    return
        "Example:\n\tLog(Y) = "
        "RACE$ * Log(AGE) + RACE$ + AGE * pow(HEIGHT, 0.5) + WEIGHT\n";
  }
  static string ExampleModel(const vector<VariableTerm>* dep_var) {
    return
        "Example:\n\t(Log(Survival_Time), Log(Censoring_Time), Status) = "
        "RACE$ * Log(AGE) + RACE$ + AGE * pow(HEIGHT, 0.5) + WEIGHT\n";
  }

  // Prompts user for the Linear Model to be solved. Makes sure that the model
  // is parsable, and if so, populates dependent_var and linear_vars with the
  // model and returns true. Returns false otherwise.
  template<typename T>
  // typename T (dependent_var) should either be:
  //   1) VariableTerm (for linear, logistic regression models)
  //   2) vector<VariableTerm> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool PromptForAndParseModel(
      const bool from_file,
      T* dependent_var, vector<LinearTerm>* linear_vars,
      char* error_msg) {
    // Sanity check input.
    if (dependent_var == nullptr || linear_vars == nullptr) {
      sprintf(error_msg, "Null input variables in PromptForAndParseModel(). "
                         "Check API in call stack and try again. Aborting.");
      return false;
    }

    // Request Linear Model from user.
    const string divider =
        "\n############################################################";
    const string data_from_file = !from_file ? "" :
        ", and should correspond to (a subset of) the names "
        "used in the input data file";
    string user_input_request = PromptUserForModelText(dependent_var);
    user_input_request += data_from_file + ".\n" + ExampleModel(dependent_var);
    cout << divider << user_input_request << " $> [Linear Model]: ";

    // Read in Linear Model.
    bool valid_input = false;
    char temp;
    while (!valid_input) {
      string user_response = "";
      cin.get(temp);
      while (temp != '\n') {
        user_response += temp;
        cin.get(temp);
      }
      char input_error[512] = "";
      if (!ProcessUserEnteredModel<T>(
              user_response, dependent_var, linear_vars, input_error)) {
        cout << divider << "\nUnable to process your model: " << input_error
             << "\n" << user_input_request << " $> [Linear Model]: ";
        continue;
      }
      valid_input = true;
    }

    return true;  
  }
  
  // Prompts user for the Linear Model to be solved. Makes sure that the model
  // specified is consistent:
  //   - The syntax is recognizable
  //   - The operations specified are supported (currently +, *, log, exp)
  //   - The NOMINAL variables never have log or exp applied to them
  // If any of the above fail, user is re-prompted to enter a valid model.
  // Populates dependent_var and linear_vars with the model; e.g. for model:
  //   Log(Y) = c_0 + c_1 * Log(X_1) * X_2 + c_2 * exp(X_2),
  // user should enter (notice contants and error term are not entered):
  //   Log(Y) = Log(X_1) * X_2 + exp(X_2)
  // and then dependent_var will represent Log(Y), and linear_vars will be a
  // vector of size two, with first term representing Log(X_1) * X_2, and
  // second term representing exp(X_2).
  // Also populates input_cols_used with the indices of all columns (from
  // the input data file) that are used in the model.
  // Returns true unless an unexpected error is encountered.
  template<typename T>
  // typename T (dependent_var) should either be:
  //   1) VariableTerm (for linear, logistic regression models)
  //   2) vector<VariableTerm> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool GetModel(
      const vector<string>& variable_names,
      const map<int, set<string> >& nominal_columns,
      T* dependent_var, vector<LinearTerm>* linear_vars,
      set<int>* input_cols_used, char* error_msg) {
    // Sanity check input.
    if (dependent_var == nullptr || linear_vars == nullptr) {
      sprintf(error_msg, "Null input variables in GetModel(). Check API in "
                         "call stack and try again. Aborting.");
      return false;
    }

    // Will need the following set of titles to make sure all the user-
    // entered variables used in the model match one of the variables
    // from the header line of the input file.
    map<string, int> titles;
    for (int i = 0; i < variable_names.size(); ++i) {
      titles.insert(make_pair(variable_names[i], i));
    }

    // Read User-Entered model.
    bool valid_model = false;
    while (!valid_model) {
      char err_msg[1024] = "";
      if (!PromptForAndParseModel<T>(true, dependent_var, linear_vars, err_msg)) {
        string err_msg_str(err_msg);
        cout << "ERROR: " << err_msg_str << "\n";
        continue;
      }

      // Sanity Check User-Entered Model.
      if (!SanityCheckModel<T>(
              titles, nominal_columns, dependent_var, linear_vars,
              input_cols_used, err_msg)) {
        string err_msg_str(err_msg);
        cout << "ERROR: " << err_msg_str << "\n";
        continue;
      }
      valid_model = true;
    }

    return true; 
  }
  // Same as above, but with no a-priori assumptions about valid variable
  // names nor which variables are nominal (only relevant when sanity checking
  // model; i.e. we won't be verifying that the variable names user specifies
  // are already present, e.g. in a data file).
  template<typename T>
  static bool GetModel(
      T* dependent_var, vector<LinearTerm>* linear_vars, char* error_msg) {
    // Sanity check input.
    if (dependent_var == nullptr || linear_vars == nullptr) {
      sprintf(error_msg, "Null input variables in GetModel(). Check API in "
                         "call stack and try again. Aborting.");
      return false;
    }

    // Read User-Entered model.
    bool valid_model = false;
    while (!valid_model) {
      char err_msg[1024] = "";
      if (!PromptForAndParseModel<T>(true, dependent_var, linear_vars, err_msg)) {
        string err_msg_str(err_msg);
        cout << "ERROR: " << err_msg_str << "\n";
        continue;
      }
      valid_model = true;
    }

    return true; 
  }

  // Given an input term (expected form 'Log(foo)', 'Exp(foo)', or simply 'foo')
  // sets term->term_title_ to either be all of input (in case 'Log' or 'Exp'
  // not detected) with term->op_ = IDENTITY, or breaks off the 'Log' or 'Exp'
  // and sets term->term_title_ to what is inside the parentheses, and sets
  // term->op_ accordingly. Returns false if 'Log' or 'Exp' was found but the
  // rest of the term was not enclosed in parentheses, returns true otherwise.
  static bool ParseLinearTerm(
      const string& input, VariableTerm* term, char* error_msg);

  // Calls ParseLinearTerm on the dep_term.
  static bool ParseDependentTerm(
      const string& input, VariableTerm* dep_term, char* error_msg) {
    return ParseLinearTerm(input, dep_term, error_msg);
  }
  // Same as above, but with different API (for Cox Proportional Hazards Model,
  // where there are two parts to the dependent variable: Time, and Status).
  static bool ParseDependentTerm(
      const string& input,
      vector<VariableTerm>* dep_term, char* error_msg);

  // Reads dependent_var into dep_var; where multiple values are parsed
  // in the case the variable is of Nominal type.
  static bool ParseDependentVariable(
      const vector<DataHolder>& sample_values,
      const map<string, int>& name_to_column,
      const map<int, set<string> >& nominal_columns,
      const VariableTerm& dependent_var,
      vector<double>* dep_var, char* error_msg);
  // Same as above, with different API for dependent_var (a vector, to hold
  // the multiple dependent variables in the Cox PH Model).
  static bool ParseDependentVariable(
      const vector<DataHolder>& sample_values,
      const map<string, int>& name_to_column,
      const map<int, set<string> >& nominal_columns,
      const vector<VariableTerm>& dependent_var,
      vector<double>* dep_var, char* error_msg);

  // 'dep_var_temp' should have two or three elements, with the last element
  // representing the Status (e.g. Delta_i, a.k.a. I(T <= C)), and the first
  // one/two elements representing the (survival, censoring) time(s). Construct
  // a CensoringData object from dep_var_temp, and push to back of dep_var.
  static bool ParseDependentVariables(
      const vector<double>& dep_var_temp,
      vector<CensoringData>* dep_var, char* error_msg);
  // Same as above, but for different API for dep_var. Note that no actual code
  // flow should come here, so we return false and log error (this is only here
  // at all because the place that calls this function is templated and doesn't
  // know that it will never actually be called; so not supporting it causes
  // compilation error.
  static bool ParseDependentVariables(
      const vector<double>& dep_var_temp,
      vector<double>* dep_var, char* error_msg) {
    sprintf(error_msg, 
            "ERROR in ParseDependentVariables: Unexpected call with wrong "
            "API (vector<double>) for dep_var. Aborting.");
    return false;
  }
  static bool ParseDependentVariables(
      const vector<double>& dep_var_temp,
      VectorXd* dep_var, char* error_msg) {
    sprintf(error_msg, 
            "ERROR in ParseDependentVariables: Unexpected call with wrong "
            "API (VectorXd) for dep_var. Aborting.");
    return false;
  }

  // Parses the string representing the RHS of a linear/logistic/cox regression
  // equation into the provided linear_vars.
  static bool ProcessIndependentVariables(
      const string& eq_rhs,
      double** constant_term, vector<LinearTerm>* linear_vars, char* error_msg);

  // Base on user-entered 'model', parses the model into the dependent_var
  // and linear_vars; returns true if model could be successfully parses,
  // and false otherwise (with error_msg containing the reason for failure).
  template<typename T>
  // typename T (dependent_var) should either be:
  //   1) VariableTerm (for linear, logistic regression models)
  //   2) vector<VariableTerm> (for Cox Proportional Hazards Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool ProcessUserEnteredModel(
      const string& model, double** constant_term,
      T* dependent_var, vector<LinearTerm>* linear_vars,
      char* error_msg) {
    // Sanity check input.
    if (dependent_var == nullptr || linear_vars == nullptr) {
      sprintf(error_msg, "NULL dependent_var or linear_vars. Check API in call "
                         "to ProcessUserEnteredModel.");
      return false;
    }
    linear_vars->clear();

    // Simplify model by removing whitespace.
    string parsed_model;
    StringUtils::RemoveAllWhitespace(model, &parsed_model);
    if (parsed_model.empty()) {
      sprintf(error_msg, "Empty model received from user input.");
      return false;
    }

    // Split model around the equality (separating (in)dependent variables).
    if (parsed_model.find("=") == string::npos) {
      sprintf(error_msg, "Input model contains no '=' sign.");
      return false;
    }
    vector<string> dep_indep_split;
    StringUtils::Split(parsed_model, "=", &dep_indep_split);
    if (dep_indep_split.size() > 2) {
      sprintf(error_msg, "Input model contains multiple '=' signs.");
      return false;
    }

    // Parse dependent variable.
    const string& dep_var = dep_indep_split[0];
    if (dep_var.empty()) {
      sprintf(error_msg, "Dependent variable (expression on left side of '=' "
                         "sign) is empty (or contains only whitespace).");
      return false;
    }
    if (!ParseDependentTerm(dep_var, dependent_var, error_msg)) return false; 

    // Parse Independent Variables.
    if (dep_indep_split.size() == 1) {
      // In some use cases (e.g. running stratified cox with subgroup, where
      // --subgroup_as_cov was true), it is valid to have the RHS of the model
      // empty (the covariates are obtained elsewhere). In this case, we just
      // want to leave linear_vars empty.
      return true;
    }
    return ProcessIndependentVariables(
        dep_indep_split[1], constant_term, linear_vars, error_msg);
  }
  // Same as above, but uses dummy variable for constant term (caller
  // will not have access to it).
  template<typename T>
  static bool ProcessUserEnteredModel(
      const string& model,
      T* dependent_var, vector<LinearTerm>* linear_vars, char* error_msg) {
    double dummy_value;
    double* dummy_ptr = &dummy_value;
    return ProcessUserEnteredModel<T>(
        model, &dummy_ptr, dependent_var, linear_vars, error_msg);
  }

  // Returns a string describing the expected format of the
  // input file to be read.
  static bool ProcessUserEnteredNominalColumns(
      const string& user_input, const vector<string>& variable_names,
      set<int>* nominal_columns, char* error_msg);
 
  // Given a model as described by dependent_var and linear_vars, populates
  // dep_var and indep_vars with the appropriate values, and prints the final
  // titles in 'legend'. For example, if model is:
  //   Log(Y) = X_1 * Log(X_2) + X_2 + X_1,
  // and X_1 is a NOMINAL variable with three possible values (so there will
  // be two corresponding indicator variables for it), then final formula is:
  //   Log(Y) = X_1,1 * Log(X_2) + X_1,2 * Log(X_2) + X_2 + X_1,1 + X_1,2
  // and legend is:
  //   [Log(Y), (X_1,1 * Log(X_2)), (X_1,2 * Log(X_2)), X_2, X_1,1, X_1,2]
  // Uses na_columns and input_cols_used to determine which rows from the input
  // data file to skip, based on that row having 'NA' in a relevant column.
  // Returns true unless as unexpected error is encountered.
  template<typename T, typename U, typename V>
  // typename T (dep_var) should either be:
  //   1) vector<CensoringData> (for Cox PH Model)
  //   2) vector<double>        (for Linear/Logistic Regression)
  //   3) VectorXd              (for Linear/Logistic Regression)
  // typename U (indep_vars) should either be:
  //   1) vector<VectorXd>
  //   2) MatrixXd
  // typename V (dependent_var) should be either:
  //   1) VariableTerm          (for Linear/Logistic Regression)
  //   2) vector<VariableTerm>  (for Cox PH Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool ParseModel(
      const bool include_error_term,
      const bool subgroup_as_cov,
      const string& subgroup_rhs,
      const vector<string>& variable_names,
      const set<int>& strata_cols,
      const vector<int>& subgroup_cols,
      const map<int, set<string> >& nominal_columns,
      const vector<vector<DataHolder> >& data_values,
      const V& dependent_var,
      const vector<LinearTerm>& linear_vars,
      const map<int, set<int> >& na_columns,
      const set<int>& input_cols_used,
      T* dep_var, U* indep_vars,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend, string* final_model, char* error_msg) {
    // Sanity Check input.
    if (dep_var == nullptr || indep_vars == nullptr || legend == nullptr) {
      sprintf(error_msg, "Null input value pased to ParseModel. Check API.");
      return false;
    }
    
    // Get Dependent and Independent variables.
    vector<VectorXd> indep_vars_temp;
    if (!ParseModel<T, V>(
          include_error_term, subgroup_as_cov, subgroup_rhs, variable_names,
          strata_cols, subgroup_cols, nominal_columns, data_values,
          dependent_var, linear_vars, na_columns, input_cols_used, dep_var, 
          &indep_vars_temp, strata, subgroup_sizes,
          legend, final_model, error_msg)) {
      return false;
    }

    // Now copy values into the provided containers.
    CopyVariables(indep_vars_temp, indep_vars);

    return true;
  }
  // Same as above, but doesn't have input parameters for strata or subgroup.
  // Instead, passes dummy values for these to the main ParseModel
  // function (which will ignore them).
  template<typename T, typename U, typename V>
  static bool ParseModel(
      const bool include_error_term,
      const vector<string>& variable_names,
      const map<int, set<string> >& nominal_columns,
      const vector<vector<DataHolder> >& data_values,
      const V& dependent_var,
      const vector<LinearTerm>& linear_vars,
      const map<int, set<int> >& na_columns,
      const set<int>& input_cols_used,
      T* dep_var, U* indep_vars,
      vector<string>* legend, string* final_model, char* error_msg) {
    set<int> dummy_strata_cols;
    vector<int> dummy_subgroup_cols;
    return ParseModel<T, U, V>(
        include_error_term, false, "", variable_names,
        dummy_strata_cols, dummy_subgroup_cols, nominal_columns, data_values,
        dependent_var, linear_vars, na_columns, input_cols_used,
        dep_var, indep_vars, nullptr, nullptr, legend, final_model, error_msg);
  }

 private:
  // Cox PH model has more than one dependent variable (Survival/Censoring time
  // and Status); for templated functions that support either, the format of the
  // input parameter 'dependent_var' determines the type of regression model.
  static bool IsCox(const VariableTerm& dependent_var) { return false; }
  static bool IsCox(const vector<VariableTerm>& dependent_var) { return true; }

  // Returns a string expression for the dependent variable(s).
  static string GetDependentVarString(const VariableTerm& dependent_var) {
    return math_utils::DataStructures::GetTermString(dependent_var);
  }
  static string GetDependentVarString(
      const vector<VariableTerm>& dependent_var) {
    string dep_var_str = "(";
    for (int i = 0; i < dependent_var.size(); ++i) {
      const VariableTerm& term = dependent_var[i];
      if (i > 0) dep_var_str += ", ";
      dep_var_str += math_utils::DataStructures::GetTermString(term);
    }
    dep_var_str += ")";
    return dep_var_str;
  }

  // Same as ParseModel() in 'public' section above, except removes one layer
  // of templatization by forcing indep_vars to be of type vector<VectorXd>.
  template<typename T, typename V>
  // typename T (dep_var) should either be:
  //   1) vector<CensoringData> (for Cox PH Model)
  //   2) vector<double>        (for Linear/Logistic Regression)
  //   3) VectorXd              (for Linear/Logistic Regression)
  // typename V (dependent_var) should be either:
  //   1) VariableTerm          (for Linear/Logistic Regression)
  //   2) vector<VariableTerm>  (for Cox PH Model)
  // This is implemented here, as Template Methods must appear in the header
  // file to avoid linking errors.
  static bool ParseModel(
      const bool include_error_term,
      const bool subgroup_as_cov,
      const string& subgroup_rhs,
      const vector<string>& variable_names,
      const set<int>& strata_cols,
      const vector<int>& subgroup_cols,
      const map<int, set<string> >& nominal_columns,
      const vector<vector<DataHolder> >& data_values,
      const V& dependent_var,
      const vector<LinearTerm>& linear_vars,
      const map<int, set<int> >& na_columns,
      const set<int>& input_cols_used,
      T* dep_var, vector<VectorXd>* indep_vars,
      vector<pair<int, int>>* strata, map<int, int>* subgroup_sizes,
      vector<string>* legend, string* final_model,
      char* error_msg) {
    // Get a mapping from the variable title to the column in the input data
    // that the title corresponds to.
    map<string, int> name_to_column;
    for (int i = 0; i < variable_names.size(); ++i) {
      name_to_column.insert(make_pair(variable_names[i], i));
    }

    // Linear and Logistic Regression uses a constant term on the RHS of the
    // Regression model, Cox Proportional Hazards Model does not. Here, we
    // distinguish them based on how many dependent variables there are
    // (Cox has at least two: one for Time, one for Status (Delta)).
    const bool is_cox = IsCox(dependent_var);

    vector<vector<string>> subgroups;
    if (!ParseSubgroups(subgroup_as_cov, subgroup_rhs, subgroup_cols.size(),
                        &subgroups)) {
      return false;
    }

    // Create a structure that will be needed to expand the 'subgroup' covariate
    // (for when subgroup_as_cov is true and there are m > 2 subgroups provided,
    // it becomes a nominal variable that is expanded into m - 1 covariates).
    set<string> subgroup_nominal_values;
    map<int, set<string>> subgroup_nominal_column;
    map<string, int> subgroup_name_to_dummy_column;
    if (subgroup_as_cov && !subgroup_cols.empty()) { 
      for (int i = 0; i < subgroups.size(); ++i) {
        subgroup_nominal_values.insert(StringUtils::Itoa(i));
      }
      subgroup_nominal_column.insert(make_pair(0, subgroup_nominal_values));
      subgroup_name_to_dummy_column.insert(make_pair(kSubgroupString, 0));
    }

    // Go through all the data values (by sample), and compute each linear term.
    legend->clear();
    vector<double> dep_var_temp;
    vector<vector<string>> strata_names;
    for (int i = 0; i < data_values.size(); ++i) {
      const vector<DataHolder>& sample_values = data_values[i];
      // Skip data rows that have 'NA' in at least one column that is needed
      // for the model.
      bool keep_row = true;
      map<int, set<int> >::const_iterator na_columns_itr = na_columns.find(i);
      if (na_columns_itr != na_columns.end()) {
        for (const int na_column : na_columns_itr->second) {
          if (input_cols_used.find(na_column) != input_cols_used.end()) {
            keep_row = false;
            break;
          }
        }
      }
      if (!keep_row) continue;

      // Fetch the subgroup for this row, if necessary. Also determines the
      // subgroup_index (which of the user-entered subgroups the data in this
      // row matches, if any), which will be used (if subgroup_as_cov is true)
      // determine the value of the Subgroup Indicator covariate.
      int subgroup_index = -1;
      if (!ParseSubgroup(
              i, subgroup_cols, subgroups, sample_values, &subgroup_index)) {
        return false;
      }
      // If subgroup_index was not set, then this row does not belong to
      // the subgroup.
      if (!subgroup_cols.empty() && subgroup_index == -1) {
        continue;
      } else if (!subgroup_cols.empty()) {
        map<int, int>::iterator size_itr = subgroup_sizes->find(subgroup_index);
        if (size_itr == subgroup_sizes->end()) {
          subgroup_sizes->insert(make_pair(subgroup_index, 1));
        } else {
          (size_itr->second)++;
        }
      }

      // Fetch the strata for this row, if necessary.
      int strata_index = -1;
      if (!ParseStrata(i, strata_cols, sample_values,
                       &strata_index, &strata_names, strata)) {
        return false;
      }

      // Generate the values for the dependent variable(s).
      if (!ParseDependentVariable(
              sample_values, name_to_column, nominal_columns, dependent_var,
              &dep_var_temp, error_msg)) {
        sprintf(error_msg, 
                "ERROR in ParseDependentVariable, row: %d", i);
        return false;
      }
      if (is_cox) {
        // Parse the dependent variables (representing Survival/Censoring Time and
        // Status/Delta) as a CensoringData object, and store in dep_var.
        if (!ParseDependentVariables(dep_var_temp, dep_var, error_msg)) {
          sprintf(error_msg, 
                  "ERROR in ParseDependentVariables, row: %d", i);
          return false;
        }
        dep_var_temp.clear();
      }

      // Generate the vector of values for the Independent Variables.
      // For Linear/Logistic Regression Models, there is a constant term
      // assumed for RHS of model. For Cox model, if subgroup_as_cov is true,
      // then there will be extra term(s) for the "Indicator" of the subgroup.
      vector<double> first_indep_values;
      vector<double>* first_indep_values_ptr = nullptr;
      vector<string> first_legend_titles;
      if (!is_cox) {
        first_indep_values.push_back(1.0);
        first_indep_values_ptr = &first_indep_values;
        first_legend_titles.push_back("Constant");
      }
      if (!subgroup_cols.empty() && subgroup_as_cov) {
        // Need to construct m - 1 (where m = number of subgroup categories)
        // covariates for the nominal category of which subgroup this row is in.
        first_indep_values_ptr = &first_indep_values;
        char tmp_error_msg[512];
        VariableTerm info_holder;
        info_holder.term_title_ = kSubgroupString;
        DataHolder subgroup_value;
        subgroup_value.name_ = StringUtils::Itoa(subgroup_index);
        vector<DataHolder> subgroup_value_holder;
        subgroup_value_holder.push_back(subgroup_value);
        VariableTerm dummy_subgroup_term;
        dummy_subgroup_term.op_ = math_utils::Operation::IDENTITY;
        dummy_subgroup_term.term_title_ = kSubgroupString;
        vector<VariableTerm> dummy_subgroup_term_holder;
        dummy_subgroup_term_holder.push_back(dummy_subgroup_term);
        if (!ComputeSampleValues(
                math_utils::Operation::MULT, "", 0, 1.0,
                subgroup_name_to_dummy_column, subgroup_nominal_column,
                subgroup_value_holder, dummy_subgroup_term_holder,
                first_indep_values_ptr, &first_legend_titles,
                tmp_error_msg)) {
          string tmp_error(tmp_error_msg);
          sprintf(error_msg, "%s\n", tmp_error.c_str());
          return false;
        }
      }
      if (!ComputeIndependentValues(
              first_legend_titles, nominal_columns, linear_vars, name_to_column,
              sample_values, first_indep_values_ptr, legend,
              indep_vars, error_msg)) {
          string tmp_error(error_msg);
          sprintf(error_msg,
                 "ERROR in ComputeSampleValues for subgroup, row: %d\n%s\n",
                 i, tmp_error.c_str());
        return false;
      }
    }

    // For Linear/Logistic Regression, we temporarily stored dependent variables
    // in dep_var_temp; copy them now into the provided container.
    // Also, start preparing strings used to represent the model.
    if (!is_cox) {
      CopyVariables(dep_var_temp, dep_var);
    }

    // Print out final model.
    string dep_var_str = GetDependentVarString(dependent_var);
    StoreFinalModel(
        is_cox, include_error_term, dep_var_str, legend, final_model);
    return true;
  }

  // Populates indep_vars with the appropriate values.
  static bool ComputeIndependentValues(
      const vector<string>& first_legend_titles,
      const map<int, set<string> >& nominal_columns,
      const vector<LinearTerm>& linear_vars,
      const map<string, int>& name_to_column,
      const vector<DataHolder>& sample_values,
      const vector<double>* first_indep_values,
      vector<string>* legend,
      vector<VectorXd>* indep_vars,
      char* error_msg);
  // Same as above, different API for first_[legend_title | indep_value]
  static bool ComputeIndependentValues(
      const string& first_legend_title,
      const map<int, set<string> >& nominal_columns,
      const vector<LinearTerm>& linear_vars,
      const map<string, int>& name_to_column,
      const vector<DataHolder>& sample_values,
      const double* first_indep_value,
      vector<string>* legend,
      vector<VectorXd>* indep_vars,
      char* error_msg) {
    vector<string> first_legend_titles;
    if (!first_legend_title.empty()) {
      first_legend_titles.push_back(first_legend_title);
    }
    vector<double>* first_indep_values = nullptr;
    vector<double> first_indep_values_holder;
    if (first_indep_value != nullptr) {
      first_indep_values_holder.push_back(*first_indep_value);
      first_indep_values = &first_indep_values_holder;
    }
    return ComputeIndependentValues(
        first_legend_titles, nominal_columns, linear_vars, name_to_column,
        sample_values, first_indep_values, legend, indep_vars, error_msg);
  }

  // Prints the final model to screen, as well as storing it in final_model.
  static void StoreFinalModel(
      const bool is_cox, const bool include_error_term,
      const string& dep_var_str, const vector<string>* legend,
      string* final_model);

  // Populate a vector with the values of the columns in sample_values, as
  // indexed by subgroup_cols. Try to match the resulting vector with one
  // of those in subgroups; if a match is found, set subgroup_index to
  // the corresponding index of subgroups, and add row_index to subgroup.
  static bool ParseSubgroup(
      const int row_index, const vector<int>& subgroup_cols,
      const vector<vector<string>>& subgroups,
      const vector<DataHolder>& sample_values,
      int* subgroup_index);

  // If strata_cols is non-empty, for indices in strata_cols, looks for the
  // corresponding values in sample_values, and puts them in a vector; then
  // sees if such a vector already exists in strata_names. If so, marks the
  // index 'j' of this vector in strata_names, and sets strata_index = j,
  // and also adds row_index to the j^th set in 'strata'. If not, pushes
  // the vector to the back of strata_names, and creates a new set (with
  // one element 'j') that it pushes to the back of strata.
  static bool ParseStrata(
      const int row_index, const set<int>& strata_cols,
      const vector<DataHolder>& sample_values,
      int* strata_index,
      vector<vector<string>>* strata_names, vector<pair<int, int>>* strata);

  // Asks user how variable 'var_name' should be sampled, and records user's
  // responses in 'sampling_params'.
  static void GetSamplingParameterFromUser(
      const string& var_name, SamplingParameters* params);

  // 'subterms' represents the individual terms of a linear term, e.g.
  // in the model:
  //   Y = Log(X_1) * X_2 + X_3
  // then for the first linear term, the subterms would be Log(X_1) and X_2,
  // and the 'op' would be multiplication (currently, this is the only op
  // that is recognized for combining subterms of a linear term).
  // This method computes each of the values needed (by looking up each
  // variable's value in data_values), and puts the value in indep_vars.
  // Note that for any NOMINAL variable, we may need to create multiple
  // indicator variables; so e.g. if X_2 is nominal with 3 distinct values,
  // the linear term Log(X_1) * X_2 becomes Log(X_1) * X_2,1 + Log(X_1) * X_2,2.
  // 'legend' will be the string representation of this linear term (or terms,
  // in the case of NOMINAL variable expansion to multiple indicator variables).
  // nominal_columns is necessary to determine the NOMINAL variables, and
  // name_to_column provides a mapping from the variable name (as recorded by
  // subterms) to the column that has the data value for it (in data_values).
  // This is a recursive method, iterating on both the number of subterms
  // in this linear term (i.e. size of 'subterms'), as well as looping over
  // all indicator variables for any NOMINAL variable types. 'depth' records
  // the current depth of the recursion (i.e. the current index in 'subterms').
  static bool ComputeSampleValues(
      const Operation op, const string& existing_title, const int depth,
      const double& existing_value, const map<string, int>& name_to_column,
      const map<int, set<string> >& nominal_columns,
      const vector<DataHolder>& data_values,
      const vector<VariableTerm>& subterms,
      vector<double>* indep_vars, vector<string>* legend, char* error_msg);
  // Same as above, but with a single subterm.
  static bool ComputeSampleValues(
      const Operation op, const string& existing_title, const int depth,
      const double& existing_value, const map<string, int>& name_to_column,
      const map<int, set<string> >& nominal_columns,
      const vector<DataHolder>& data_values,
      const VariableTerm& subterm,
      vector<double>* indep_vars, vector<string>* legend, char* error_msg) {
    vector<VariableTerm> subterms;
    subterms.push_back(subterm);
    return ComputeSampleValues(
        op, existing_title, depth, existing_value, name_to_column,
        nominal_columns, data_values, subterms, indep_vars, legend, error_msg);
  }

  // A helper function for ComputeSampleValues() above. See description there
  // for details.
  static bool UpdateValueAndIterate(
      const Operation op, const string& existing_title, const int depth,
      const double& existing_value, const double& term_value,
      const map<string, int>& name_to_column,
      const map<int, set<string> >& nominal_columns,
      const vector<DataHolder>& data_values,
      const vector<VariableTerm>& subterms,
      vector<double>* indep_vars, vector<string>* legend, char* error_msg);

};

}  // namespace file_reader_utils

#endif
