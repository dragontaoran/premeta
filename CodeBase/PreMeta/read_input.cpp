#include "read_input.h"
#include "MathUtils/data_structures.h"
#include "StringUtils/string_utils.h"

#include <Eigen/Dense>
#include <errno.h>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using math_utils::CensoringData;
using math_utils::DataHolder;
using math_utils::Distribution;
using math_utils::LinearTerm;
using math_utils::Operation;
using math_utils::SamplingParameters;
using math_utils::VariableTerm;
using string_utils::StringUtils;

static const char kErrorTermId[] = "PHB_error_term_id_PHB";

namespace file_reader_utils {

string ReadInput::GetErrorTermId() { return kErrorTermId; }

bool ReadInput::IsTrueIndicator(const string& value) {
  return (value == "1") || (value == "1.0") || (value == "T") ||
         (value == "True") || (value == "true") || (value == "TRUE");
}

bool ReadInput::IsFalseIndicator(const string& value) {
  return (value == "0") || (value == "0.0") || (value == "F") ||
         (value == "False") || (value == "false") || (value == "FALSE");
}

int ReadInput::GetIntegerFromUser(const string& query_text) {
  int n;
  bool valid_input = false;
  char temp;
  cout << query_text;
  while (!valid_input) {
    string user_response = "";
    cin.get(temp);
    while (temp != '\n') {
      user_response += temp;
      cin.get(temp);
    }
    if (!StringUtils::Stoi(user_response, &n)) {
      cout << "\nUnable to process your input '" << user_response
           << "' as an integer.\n" 
           << query_text;
      continue;
    }
    valid_input = true;
  }
  return n;
}

double ReadInput::GetDoubleFromUser(const string& query_text) {
  double d;
  bool valid_input = false;
  char temp;
  cout << query_text;
  while (!valid_input) {
    string user_response = "";
    cin.get(temp);
    while (temp != '\n') {
      user_response += temp;
      cin.get(temp);
    }
    if (!StringUtils::Stod(user_response, &d)) {
      cout << "\nUnable to process your input '" << user_response
           << "' as a numerical value.\n" 
           << query_text;
      continue;
    }
    valid_input = true;
  }
  return d;
}

int ReadInput::GetNumDataRowsPerSimulationFromUser() {
  bool valid_input = false;
  char temp;
  const string user_input_request =
      "\nPlease enter the number of sample rows (input data points) 'n'\n"
      " $> [Enter value for 'n']: ";
  return GetIntegerFromUser(user_input_request);
}

int ReadInput::GetNumSimulationsFromUser() {
  bool valid_input = false;
  char temp;
  const string user_input_request =
      "\nPlease enter the number of simulations 'k' you want to run\n"
      " $> [Enter value for 'k']: ";
  return GetIntegerFromUser(user_input_request);
}

bool ReadInput::GetColumnIndicesFromString(
    const string& strata, const vector<string>& header, set<int>* strata_cols) {
  if (strata.empty()) return true;
  vector<string> titles;
  string strata_stripped = StringUtils::StripPrefixString(strata, "(");
  strata_stripped = StringUtils::StripSuffixString(strata_stripped, ")");
  StringUtils::Split(strata_stripped, ",", &titles);
  for (const string& title : titles) {
    bool found_match = false;
    for (int i = 0; i < header.size(); ++i) {
      if (header[i] == title) {
        strata_cols->insert(i);
        found_match = true;
        break;
      }
    }
    if (!found_match) {
      cout << "ERROR in GetColumnIndicesFromString: "
           << "Unable to find strata column '"
           << title << "' among the titles in the Header row of the input file:"
           << endl << StringUtils::Join(header) << endl;
      return false;
    }
  }
  return true;
}

bool ReadInput::GetColumnIndicesFromString(
    const string& subgroup, const vector<string>& header,
    vector<int>* subgroup_cols) {
  vector<string> titles;
  string subgroup_stripped = StringUtils::StripPrefixString(subgroup, "(");
  subgroup_stripped = StringUtils::StripSuffixString(subgroup_stripped, ")");
  StringUtils::Split(subgroup_stripped, ",", &titles);
  for (const string& title : titles) {
    bool found_match = false;
    for (int i = 0; i < header.size(); ++i) {
      if (header[i] == title) {
        subgroup_cols->push_back(i);
        found_match = true;
        break;
      }
    }
    if (!found_match) {
      cout << "ERROR in GetColumnIndicesFromString: "
           << "Unable to find subgroup column '"
           << title << "' among the titles in the Header row of the input file:"
           << endl << StringUtils::Join(header) << endl;
      return false;
    }
  }
  return true;
}

bool ReadInput::GetSubgroupColumns(
    const string& subgroup_str, const vector<string>& header,
    string* subgroup_rhs, vector<int>* subgroup_cols) {
  // Nothing to do if no subgroup to parse.
  if (subgroup_str.empty()) return true;

  vector<string> subgroup_terms;
  StringUtils::Split(subgroup_str, "=", &subgroup_terms);
  if (subgroup_terms.size() != 2) {
    cout << "ERROR in ReadFileAndParseModel parsing subgroup '"
         << subgroup_str << "': Expected Column Title(s) on LHS of an "
         << "equals sign, and a list of value-tuples on the RHS."
         << endl;
    return false;
  }
  StringUtils::RemoveAllWhitespace(subgroup_terms[1], subgroup_rhs);
  if (!GetColumnIndicesFromString(subgroup_terms[0], header, subgroup_cols)) {
    cout << "ERROR in ReadFileAndParseModel getting subgroup columns: "
         << "unable to parse spec '" << subgroup_str << "'." << endl;
    return false;
  }
  return true;
}

bool ReadInput::ParseSubgroups(
    const bool subgroup_as_cov,
    const string& subgroup_rhs, const int num_expected_subgroups,
    vector<vector<string>>* subgroups) {
  if (num_expected_subgroups == 0) return true;
  string subgroup_rhs_stripped;
  StringUtils::RemoveAllWhitespace(subgroup_rhs, &subgroup_rhs_stripped);
  subgroup_rhs_stripped =
      StringUtils::StripPrefixString(subgroup_rhs_stripped, "{");
  subgroup_rhs_stripped =
      StringUtils::StripSuffixString(subgroup_rhs_stripped, "}");
  if (subgroup_rhs_stripped.empty()) {
    cout << "\nERROR in ParseSubgroups: No RHS found for subgroup: "
         << subgroup_rhs << endl;
    return false;
  }
  vector<string> rhs_terms;
  StringUtils::Split(subgroup_rhs_stripped, "),(", &rhs_terms);
  for (string& term : rhs_terms) {
    term = StringUtils::StripPrefixString(term, "(");
    term = StringUtils::StripSuffixString(term, ")");
    subgroups->push_back(vector<string>());
    vector<string>& subgroup = subgroups->back();
    StringUtils::Split(term, ",", &subgroup);
    if (subgroup.size() != num_expected_subgroups) {
      cout << "\nERROR in ParseSubgroups: subgroup '" << term
           << "' does not have the expected number of column values ("
           << num_expected_subgroups << ")." << endl;
      return false;
    }
  }

  // If there is only one Subgroup category, we cannot use subgroup as
  // a covariate (all rows will have the same value for this covariate).
  if (subgroup_as_cov && subgroups->size() < 2) {
    cout << "\nERROR in ParseSubgroups: Only one subgroup type was provided, "
         << "so subgroup cannot be used as a covariate (non-invertible "
         << "information matrix). You must either not use subgroup as a "
         << "covariate (e.g. option --nosubgroup_as_cov), or specify two "
         << "or more subgroup types." << endl;
    return false;
  }

  // Sanity check all subgroups are distinct.
  set<string> unique_subgroups;
  const string& kDummySeperator = "PHB_FOO_PHB";
  for (const vector<string>& subgroup_terms : *subgroups) {
    const string subgroup_concat =
        StringUtils::Join(subgroup_terms, kDummySeperator);
    if (unique_subgroups.find(subgroup_concat) == unique_subgroups.end()) {
      unique_subgroups.insert(subgroup_concat);
    } else {
      cout << "\nERROR in ParseSubgroups: Detected duplicate subgroup "
           << " categories ('" << StringUtils::Join(subgroup_terms, ",")
           << "') on RHS of --subgroup_model expression: " << subgroup_rhs
           << endl;
      return false;
    }
  }

  return true;
}

bool ReadInput::ParseSubgroup(
    const int row_index, const vector<int>& subgroup_cols,
    const vector<vector<string>>& subgroups,
    const vector<DataHolder>& sample_values,
    int* subgroup_index) {
  // It is valid that some calls to this shouldn't do anything; such
  // calls will have NULL subgroup_cols and NULL subgroup.
  if (subgroup_cols.empty()) return true;

  // Go through sample_values, picking out the values in the relevant
  // (subgroup) columns, and storing them in values_for_subgroup_cols.
  vector<string> values_for_subgroup_cols;
  for (int i = 0; i < subgroup_cols.size(); ++i) {
    const int subgroup_col = subgroup_cols[i];
    if (subgroup_col >= sample_values.size()) {
      cout << "\nERROR in ParseSubgroup: subgroup_col specifies a column index ("
           << subgroup_col << ") that is bigger than the number of columns in "
           << "input data (" << sample_values.size() << "). Aborting." << endl;
      return false;
    }
    const DataHolder& data = sample_values[subgroup_col];
    const string& col_value =
        data.name_.empty() ? StringUtils::Itoa(data.value_) : data.name_;
    values_for_subgroup_cols.push_back(col_value);
  }

  // Try to match the values in values_for_subgroup_cols to one of the
  // indices in subgroups. If a match is found, set subgroup index
  // accordingly.
  for (int i = 0; i < subgroups.size(); ++i) {
    const vector<string>& subgroup_type = subgroups[i];
    if (equal(values_for_subgroup_cols.begin(), values_for_subgroup_cols.end(),
              subgroup_type.begin())) {
      *subgroup_index = i;
      return true;
    }
  }
  return true;
}

bool ReadInput::ParseStrata(
    const int row_index, const set<int>& strata_cols,
    const vector<DataHolder>& sample_values,
    int* strata_index,
    vector<vector<string>>* strata_names, vector<pair<int, int>>* strata) {
  // It is valid that some calls to this shouldn't do anything; such
  // calls will have empty strata_cols and have strata = NULL.
  if (strata_cols.empty()) {
    strata->push_back(make_pair(row_index, 0));
    return true;
  }

  vector<string> strata_values;
  for (const int col_index : strata_cols) {
    // Sanity check it is a valid column index.
    if (col_index >= sample_values.size()) {
      cout << "\nERROR in ParseStrata: strata_cols specifies a column index ("
           << col_index << ") that is bigger than the number of columns in "
           << "input data (" << sample_values.size() << "). Aborting." << endl;
      return false;
    }
    const DataHolder& data = sample_values[col_index];
    strata_values.push_back(
        data.name_.empty() ? StringUtils::Itoa(data.value_) : data.name_);
  }
  

  // Look to see if this strata_name is already present in strata_names.
  bool found_match = false;
  for (int i = 0; i < strata_names->size(); ++i) {
    const vector<string>& strata_name = (*strata_names)[i];
    if (equal(strata_name.begin(), strata_name.end(), strata_values.begin())) {
      // A strata with this name already exists. Add 'row_index' to it.
      *strata_index = i;
      strata->push_back(make_pair(row_index, i));
      found_match = true;
      break;
    }
  }
  if (!found_match) {
    // This name not in strata_names. Create a new strata with this name,
    // and add 'row_index' to it.
    strata->push_back(make_pair(row_index, strata_names->size()));
    strata_names->push_back(strata_values);
  }
  
  return true;
}

bool ReadInput::ParseSamplingParameter(const string& value,
                                       SamplingParameters* params) {
  if (value.empty()) return false;

  // Distribution description should end in closing parentheses.
  string temp;
  if (!StringUtils::StripSuffixString(value, ")", &temp)) return false;
  string clipped_value;
  // Normal Distribution.
  if (StringUtils::StripPrefixString(temp, "N(", &clipped_value)) {
    size_t seperator = clipped_value.find(",");
    // Make sure there is a comma with characters before and after it.
    if (seperator == string::npos || seperator == 0 ||
        seperator == clipped_value.length() - 1) {
      return false;
    }
    double mean, std_dev;
    if (!StringUtils::Stod(clipped_value.substr(0, seperator), &mean) ||
        !StringUtils::Stod(clipped_value.substr(seperator + 1), &std_dev)) {
      return false;
    }
    params->type_ = Distribution::NORMAL;
    params->mean_ = mean;
    params->std_dev_ = std_dev;
    return true;
  } else if (StringUtils::StripPrefixString(temp, "U(", &clipped_value)) {
    size_t seperator = clipped_value.find(",");
    // Make sure there is a comma with characters before and after it.
    if (seperator == string::npos || seperator == 0 ||
        seperator == clipped_value.length() - 1) {
      return false;
    }
    double start, end;
    if (!StringUtils::Stod(clipped_value.substr(0, seperator), &start) ||
        !StringUtils::Stod(clipped_value.substr(seperator + 1), &end)) {
      return false;
    }
    params->type_ = Distribution::UNIFORM;
    params->range_start_ = start;
    params->range_end_ = end;
    return true;
  } else if (StringUtils::StripPrefixString(temp, "P(", &clipped_value)) {
    double mean;
    if (!StringUtils::Stod(clipped_value, &mean)) {
      return false;
    }
    params->type_ = Distribution::POISSON;
    params->mean_ = mean;
    return true;
  }
  return false;
}

void ReadInput::GetSamplingParameterFromUser(
    const string& var_name, SamplingParameters* params) {
  if (params == nullptr || var_name.empty()) return;

  // Get Distribution type.
  bool valid_type = false;
  char temp;
  string user_input_request =
      "\nPlease enter the Distribution Type from which you want to sample "
      "for variable '";
  user_input_request += var_name + "'.\nFormat should be: U(a,b); N(a,b); " +
                        "or P(a)\nWhere U(a,b) is for Uniform on [a,b]; " +
                        "N(a,b) is for Normal on [a,b]; and P(a) is Poisson." +
                        "Default (if nothing entered) is N(0,1)\n  " +
                        "$> [Enter Distribution Type]: ";
  cout << user_input_request;
  while (!valid_type) {
    string user_response = "";
    cin.get(temp);
    while (temp != '\n') {
      user_response += temp;
      cin.get(temp);
    }
    if (user_response.empty()) {
      params->type_ = Distribution::NORMAL;
      params->mean_ = 0.0;
      params->std_dev_ = 1.0;
    } else if (!ParseSamplingParameter(user_response, params)) {
      cout << "\nUnable to process your input '" << user_response
           << "' as one of the available distribution types.\n" 
           << user_input_request;
      continue;
    }
    valid_type = true;
  }
}

bool ReadInput::GetSamplingParametersFromUser(
    const vector<LinearTerm>& linear_vars,
    map<string, SamplingParameters>* sampling_params, char* error_msg) {
  // Sanity check input.
  if (sampling_params == nullptr) {
    sprintf(error_msg, "Null input to function GetSamplingParametersFromUser. "
                       "Check calling method.");
    return false;
  }

  // Get sampling parameters for error term.
  sampling_params->clear();
  SamplingParameters error_params;
  GetSamplingParameterFromUser("Error Term", &error_params);
  sampling_params->insert(make_pair(kErrorTermId, error_params));

  // Loop through all linear terms, asking user to input sampling parameters
  // for each variable encountered.
  for (const LinearTerm& term : linear_vars) {
    // Loop through all subterms of each linear term.
    for (const VariableTerm& subterm : term.terms_) {
      map<string, SamplingParameters>::const_iterator sampling_itr =
          sampling_params->find(subterm.term_title_);
      if (sampling_itr == sampling_params->end()) {
        SamplingParameters new_params;
        GetSamplingParameterFromUser(subterm.term_title_, &new_params);
        sampling_params->insert(make_pair(subterm.term_title_, new_params));
      }
    }
  }
  return true;
}

bool ReadInput::SanityCheckDependentVariable(
    const map<string, int>& titles,
    const map<int, set<string> >& nominal_columns,
    const VariableTerm* dependent_var,
    set<int>* input_cols_used, char* error_msg) {
  // Sanity-check dep variable name matches one of the input variables.
  map<string, int>::const_iterator title_itr =
      titles.find(dependent_var->term_title_);
  if (title_itr == titles.end()) {
    string valid_titles;
    for (const auto& title : titles) {
      if (!valid_titles.empty()) valid_titles += ", ";
      valid_titles += title.first;
    }
    sprintf(error_msg, "Unable to find '%s' among the variable names found "
                       "on the header line of the input file.\nVariable "
                       "names found on Header Line: %s\n",
            dependent_var->term_title_.c_str(), valid_titles.c_str());
    return false;
  }
  const int dep_var_column = title_itr->second;
  input_cols_used->insert(dep_var_column);

  // Sanity check operation specified for the dependent variable:
  // should not apply 'log' or 'exp' to a nominal variable.
  map<int, set<string> >::const_iterator itr =
      nominal_columns.find(dep_var_column);
  if (itr != nominal_columns.end() &&
      dependent_var->op_ != Operation::IDENTITY) {
    sprintf(error_msg, "Dependent variable '%s' is of NOMINAL type, so you "
                       "cannot apply operation %d to it.\n",
            dependent_var->term_title_.c_str(), dependent_var->op_);
    return false;
  }
  if (itr != nominal_columns.end() && itr->second.size() != 2) {
    sprintf(error_msg, "Dependent variable '%s' is of NOMINAL type, are you "
                       "sure you intended this? Note that this is only valid "
                       "if this column has exactly two distinct values across "
                       "all samples, but the input data has %d values.\n",
            dependent_var->term_title_.c_str(), itr->second.size());
    return false;
  }
  return true;
}

bool ReadInput::SanityCheckDependentVariable(
    const map<string, int>& titles,
    const map<int, set<string> >& nominal_columns,
    const vector<VariableTerm>* dependent_var,
    set<int>* input_cols_used, char* error_msg) {
  // Sanity-check dependent_var is either size two or three (one or two
  // dependent variables, plus the status variable).
  if (dependent_var->size() != 2 && dependent_var->size() != 3) {
    sprintf(error_msg, "ERROR in SanityCheckDependentVariable: improper "
                       "number of elements (%d) in dependent_var.",
            dependent_var->size());
    return false;
  }

  // Sanity-check dep variable name matches one of the input variables.
  vector<int> columns;
  for (const VariableTerm& term : *dependent_var) {
    const string& term_title = term.term_title_;
    map<string, int>::const_iterator title_itr = titles.find(term_title);
    if (title_itr == titles.end()) {
      string valid_titles;
      for (const auto& title : titles) {
        if (!valid_titles.empty()) valid_titles += ", ";
        valid_titles += title.first;
      }
      sprintf(error_msg, "Unable to find '%s' among the variable "
                         "names found on the header line of the input file.\n"
                         "Variable names found on Header Line: %s\n",
              term_title.c_str(), valid_titles.c_str());
      return false;
    }
    columns.push_back(title_itr->second);
    input_cols_used->insert(title_itr->second);
  }

  // Sanity check operation specified for the dependent variable:
  // should not apply 'log' or 'exp' to a nominal variable.
  // Also, dependent variable can only be a nominal variable if it
  // is the Status variable, and there are only two values (True, False).
  const int status_index = dependent_var->size() - 1;
  for (int i = 0; i < columns.size(); ++i) {
    int column = columns[i];
    const VariableTerm& term = (*dependent_var)[i];
    map<int, set<string> >::const_iterator itr = nominal_columns.find(column);
    if (itr != nominal_columns.end() && term.op_ != Operation::IDENTITY) {
      sprintf(error_msg, "Dependent variable '%s' is of NOMINAL type, so you "
                         "cannot apply operation %d to it.\n",
              term.term_title_.c_str(), term.op_);
      return false;
    }
    if (itr != nominal_columns.end() &&
        (i < status_index || itr->second.size() != 2)) {
      sprintf(error_msg, "Dependent variable '%s' is of NOMINAL type, are you "
                         "sure you intended this? Note that this is only valid "
                         "if this column has exactly two distinct values across "
                         "all samples, but the input data has %d values.\n",
              term.term_title_.c_str(), itr->second.size());
      return false;
    }
  }
  return true;
}

bool ReadInput::SanityCheckIndependentVariables(
    const map<string, int>& titles,
    const map<int, set<string> >& nominal_columns,
    const vector<LinearTerm>* linear_vars,
    set<int>* input_cols_used, char* error_msg) {
  for (const LinearTerm& itr : *linear_vars) {
    // Currently, the only supported operation for LinearTerm is MULT.
    if (itr.op_ != Operation::MULT) {
      sprintf(error_msg, "Only Multiplication is currently supported "
                         "for combining subterms of a linear term.\n");
      return false;
    }
    // Check all sub-terms of this linear term.
    for (const VariableTerm& term : itr.terms_) {
      // Sanity-check indep variable name matches one of the input variables.
      map<string, int>::const_iterator term_title_itr =
          titles.find(term.term_title_);
      if (term_title_itr == titles.end()) {
        string valid_titles;
        for (const auto& title : titles) {
          if (!valid_titles.empty()) valid_titles += ", ";
          valid_titles += title.first;
        }
        sprintf(error_msg, "Cannot find '%s' among the variable names found "
                           "on the header line of the input file.\nVariable "
                           "names found on Header Line: %s\n",
                term.term_title_.c_str(), valid_titles.c_str());
        return false;
      }
      const int indep_var_column = term_title_itr->second;
      input_cols_used->insert(indep_var_column);

      // Sanity check operation specified for the independent variable:
      // should not apply any operation to a nominal variable.
      if (term.op_ != Operation::IDENTITY &&
          nominal_columns.find(indep_var_column) != nominal_columns.end()) {
        const string op_name =
            term.op_ == Operation::MULT ? "Multiplication" :
            term.op_ == Operation::ADD ? "Addition" :
            term.op_ == Operation::EXP ? "Exponential" :
            term.op_ == Operation::SQRT ? "Square Root" :
            term.op_ == Operation::POW ? "Exponent" :
            term.op_ == Operation::LOG ? "Log" : "Unknown Operation";
        sprintf(error_msg, "Independent variable '%s' is of NOMINAL type, so "
                           "you cannot apply %s to it.\n",
                term.term_title_.c_str(), op_name.c_str());
        return false;
      }
    }
  }
  return true;
}

bool ReadInput::ParseDependentVariable(
    const vector<DataHolder>& sample_values,
    const map<string, int>& name_to_column,
    const map<int, set<string> >& nominal_columns,
    const VariableTerm& dependent_var,
    vector<double>* dep_var, char* error_msg) {
  int old_size = dep_var->size();
  // Parse Dependent variable.
  if (!ComputeSampleValues(
        Operation::MULT, "", 0, 1.0, name_to_column, nominal_columns,
        sample_values, dependent_var, dep_var, nullptr, error_msg)) {
    return false;
  }
  // Dependent Variable should not be NOMINAL, and so size of dep_var should
  // have grown by one (NOMINAL variables expand the number of variables).
  if (dep_var->size() - old_size != 1) {
    sprintf(error_msg, "ERROR in ParseDependentVariable: Too many variables "
                       "required for the dependent variable (found %d)",
            dep_var->size() - old_size);
    return false;
  }
  return true;
}

bool ReadInput::ParseDependentVariable(
    const vector<DataHolder>& sample_values,
    const map<string, int>& name_to_column,
    const map<int, set<string> >& nominal_columns,
    const vector<VariableTerm>& dependent_var,
    vector<double>* dep_var, char* error_msg) {
  // Sanity check input.
  if (dependent_var.size() != 2 && dependent_var.size() != 3) {
    sprintf(error_msg, "ERROR in ParseDependentVariable: Too many terms (%d) "
                       "in dependent_var.", dependent_var.size());
    return false;
  }

  // Parse each Dependent variable (one or two 'time' variables, plus the
  // 'status' variable).
  for (const VariableTerm& term : dependent_var) {
    int old_size = dep_var->size();
    if (!ComputeSampleValues(
          Operation::MULT, "", 0, 1.0, name_to_column, nominal_columns,
          sample_values, term, dep_var, nullptr, error_msg)) {
      return false;
    }
    // Dependent Variables should not be NOMINAL, and so number of dependent
    // variables computed should match dependent_var.size() (NOMINAL variables
    // will expand the number of variables).
    if (dep_var->size() - old_size != 1) {
      sprintf(error_msg, "ERROR in ParseDependentVariable: Too many variables "
                         "required for the dependent variable (found %d)",
              dep_var->size() - old_size);
      return false;
    }
  }
  return true;
}

bool ReadInput::ParseDependentVariables(
    const vector<double>& dep_var_temp,
    vector<CensoringData>* dep_var, char* error_msg) {
  // Sanity check input.
  if (dep_var_temp.size() != 2 && dep_var_temp.size() != 3) {
    sprintf(error_msg, "ERROR in ParseDependentVariables: Too many terms (%d) "
                       "in dependent_var.", dep_var_temp.size());
    return false;
  }

  // Sanity-check Status variable is True or False.
  const int status_index = dep_var_temp.size() == 2 ? 1 : 2;
  const double& status = dep_var_temp[status_index];
  if (status != 0.0 && status != 1.0) {
    sprintf(error_msg, "Status should be 0.0 (false) or 1.0 (true); found: "
                       "%0.04f", status);
    return false;
  }
  const bool is_alive = status == 0.0;

  // Parse the Dependent variables based on their values.
  CensoringData data;
  data.is_alive_ = is_alive;
  if (dep_var_temp.size() == 2) {
    if (is_alive) {
      data.censoring_time_ = dep_var_temp[0];
      // Artificially set survival time to be bigger than censoring_time_,
      // to ensure it is not used: All functions using CensoringData should
      // look at is_alive_ to determine which (among survival_time_ and
      // censoring_time_) is used; this is just a safeguard, in case there
      // are accidentally any places that use the minimum of the two, in
      // which case we don't want the default value of 0.0 for an
      // uninitialized survival_time_ to be the min.
      data.survival_time_ = dep_var_temp[0] + 1.0;
    } else {
      data.survival_time_ = dep_var_temp[0];
      // Artificially set censoring_time_ to be bigger than survival_time_,
      // for same reason as above (in reverse).
      data.censoring_time_ = dep_var_temp[0] + 1.0;
    }
  } else {
    data.survival_time_ = dep_var_temp[0];
    data.censoring_time_ = dep_var_temp[1];
  }
  dep_var->push_back(data);
  return true;
}

void ReadInput::StoreFinalModel(
    const bool is_cox, const bool include_error_term,
    const string& dep_var_str, const vector<string>* legend,
    string* final_model) {
  string constant_term_in_title;
  string constant_term_in_model;
  if (!is_cox) {
    constant_term_in_title = ", constant term (B_0)";
    constant_term_in_model = "B_0";
  }
  string indep_vars_str = " = ";
  for (int i = 0; i < legend->size(); ++i) {
    const string& legend_term = (*legend)[i];
    int j = is_cox ? i + 1 : i;
    if (j == 0) {
      // Constant term.
      indep_vars_str += constant_term_in_model;
      continue;
    }
    // If we already have at least one term (the constant term), add
    // the proper '+' sepearator (here, length 3 is because indep_vars_str
    // is initialized with having 3 characters for ' = ').
    if (indep_vars_str.length() > 3) {
      indep_vars_str += " + ";
    }
    indep_vars_str += "B_" + StringUtils::Itoa(j) + " * " + legend_term;
  }
  const string error_term_in_title = include_error_term ?
      ", and error term (e)):\n\t" : "):\n\t";
  const string error_term_in_model =
      include_error_term ? " + e\n" : "\n";
  const string title = "Model";
  const string model_msg =
      title + " (includes indicator variables for NOMINAL variables" +
      constant_term_in_title + error_term_in_title + dep_var_str +
      indep_vars_str + error_term_in_model;
  if (final_model != nullptr) *final_model = model_msg;
}

bool ReadInput::ComputeIndependentValues(
    const vector<string>& first_legend_titles,
    const map<int, set<string> >& nominal_columns,
    const vector<LinearTerm>& linear_vars,
    const map<string, int>& name_to_column,
    const vector<DataHolder>& sample_values,
    const vector<double>* first_indep_values,
    vector<string>* legend,
    vector<VectorXd>* indep_vars,
    char* error_msg) {
  // Only update legend once (on the first sample row of data). Note
  // this isn't strictly necessary, I could re-write it every time,
  // or write it once seperately outside of the loop over all sample
  // data rows; but the former seems extraneous and the latter would
  // involve duplicating work already being done here.
  vector<string>* do_legend_once = nullptr;
  if (legend != nullptr && legend->empty()) {
    for (const string& title : first_legend_titles) {
      legend->push_back(title);
    }
    do_legend_once = legend;
  }

  vector<double> indep_values;
  if (first_indep_values != nullptr) {
    for (const double& value : *first_indep_values) {
      indep_values.push_back(value);
    }
  }
  for (const LinearTerm& linear_term : linear_vars) {
    if (!ComputeSampleValues(
          linear_term.op_, "" /* Term Name */, 0 /* Initial Depth */,
          linear_term.constant_, name_to_column,
          nominal_columns, sample_values, linear_term.terms_,
          &indep_values, do_legend_once, error_msg)) {
      return false;
    }
  }

  VectorXd new_sample_row;
  new_sample_row.resize(indep_values.size());
  for (int j = 0; j < indep_values.size(); ++j) {
    new_sample_row[j] = indep_values[j];
  }
  indep_vars->push_back(new_sample_row);
  return true;
}


bool ReadInput::ComputeSampleValues(
    const Operation op, const string& existing_title, const int depth,
    const double& existing_value, const map<string, int>& name_to_column,
    const map<int, set<string> >& nominal_columns,
    const vector<DataHolder>& data_values,
    const vector<VariableTerm>& subterms,
    vector<double>* indep_vars, vector<string>* legend, char* error_msg) {
  // Sanity check input.
  if (indep_vars == nullptr) {
    sprintf(error_msg, "Null inputs to ComputeSampleValues. Check API.");
    return false;
  }
  if (subterms.size() <= depth) {
    sprintf(error_msg, "Depth (%d) exceeds size of 'subterms' (%d) inside "
                       "call to ComputeSampleValues().",
            depth, subterms.size());
    return false;
  }

  // Sanity check that the current term (variable) has a valid title.
  VariableTerm current_term = subterms[depth];
  map<string, int>::const_iterator name_to_col_itr =
      name_to_column.find(current_term.term_title_);
  if (name_to_col_itr == name_to_column.end()) {
    sprintf(error_msg, "Unrecognized variable '%s' does not appear in the "
                       "header line of the input file. Aborting.",
            current_term.term_title_.c_str());
    return false;
  }

  // Set local constants that will be needed below.
  const int index = name_to_col_itr->second;
  const string current_title =
      math_utils::DataStructures::GetTermString(current_term);
  string new_title = existing_title.empty() ? current_title :
      existing_title + math_utils::DataStructures::GetOpString(op) +
      current_title;

  // Check if the current term corresponds to a NOMINAL variable.
  map<int, set<string> >::const_iterator nominal_column_itr =
      nominal_columns.find(index);
  if (nominal_column_itr == nominal_columns.end()) {
    // Current term is NOT a nominal variable.
    if (!UpdateValueAndIterate(
          op, new_title, depth, existing_value, data_values[index].value_,
          name_to_column, nominal_columns, data_values, subterms,
          indep_vars, legend, error_msg)) {
      return false;
    }
    return true;
  } else {
    // Current term is a NOMINAL variable. Determine how many indicator variables
    // are necessary (based on how many distinct values appear for this variable),
    // and create this many indicator variables.
    const int indicators_needed = (nominal_column_itr->second.size() - 1);
    // Create k - 1 indicator values (k = number distinct values for this
    // NOMINAL variable).
    bool first_nominal_var = true;
    for (const string& indicator : nominal_column_itr->second) {
      double indicator_value;
      string new_title_with_subscript = new_title;
      if (indicators_needed == 0) {
        // There is only one distinct value that appears for this NOMINAL
        // variable. Handle this special case seperately.
        indicator_value = 1.0;
        new_title_with_subscript = existing_title;
        // Print out warning message. Just do this once, no need to
        // do it for every sample row of data; we use whether legend ==
        // nullptr as a proxy for doing this only once.
        if (legend != nullptr) {
          cout << "\nWARNING: Variable '" << current_term.term_title_
               << "' is of NOMINAL type, with only one distinct value ("
               << indicator << ").\nIgnoring it in the model.\n";
        }
      } else if (first_nominal_var) {
        // Should loop through one fewer time than number of values in
        // nominal_column_itr->second.
        first_nominal_var = false;
        continue;
      } else if (data_values[index].name_ == indicator) {
        indicator_value = 1.0;
        new_title_with_subscript += "_" + indicator;
      } else {
        indicator_value = 0.0;
        new_title_with_subscript += "_" + indicator;
      }
      if (!UpdateValueAndIterate(
            op, new_title_with_subscript, depth, existing_value,
            indicator_value, name_to_column, nominal_columns,
            data_values, subterms, indep_vars, legend, error_msg)) {
        return false;
      }
    }
    return true;
  }
}

bool ReadInput::UpdateValueAndIterate(
    const Operation op, const string& existing_title, const int depth,
    const double& existing_value, const double& term_value,
    const map<string, int>& name_to_column,
    const map<int, set<string> >& nominal_columns,
    const vector<DataHolder>& data_values,
    const vector<VariableTerm>& subterms,
    vector<double>* indep_vars, vector<string>* legend, char* error_msg) {
  // Update existing_value.
  double new_value;
  if (!math_utils::DataStructures::ComputeGroupOperation(
          op, existing_value,
          math_utils::DataStructures::ComputeSelfOperation(subterms[depth],
                                                           term_value),
          &new_value)) {
    sprintf(error_msg, "Unable to ComputeGroupOperation. Aborting.");
    return false;
  }
  // Check if this is the last subterm.
  if (depth + 1 == subterms.size()) {
    indep_vars->push_back(new_value);
    if (legend != nullptr) legend->push_back(existing_title);
    return true;
  } else {
    // There are more subterms. Call ComputeSampleValues iteratively on
    // the remaining ones.
    return ComputeSampleValues(
        op, existing_title, depth + 1, new_value, name_to_column,
        nominal_columns, data_values, subterms, indep_vars, legend, error_msg);
  }
}

bool ReadInput::ProcessIndependentVariables(
    const string& eq_rhs, double** constant_term,
    vector<LinearTerm>* linear_vars, char* error_msg) {
  // Parse independent variables.
  // Loop over all linear terms (linear terms are those seperated by '+').
  vector<string> linear_terms;
  StringUtils::Split(eq_rhs, "+", &linear_terms);
  bool constant_term_present = false;
  for (const string& term : linear_terms) {
    vector<VariableTerm> mult_terms;
    LinearTerm to_add;
    // Coefficients for each linear term is assumed to be 1.0, unless otherwise
    // specified (i.e. this can be updated below by user-entered value).
    to_add.constant_ = 1.0;
    vector<string> sub_terms;
    // Parse each subterm: loop over the 'subterms' of the current linear term
    // (subterms are seperated by '*' e.g. 2.0*X_1*X_2*Log(X_3) has 4 subterms).
    StringUtils::Split(term, "*", &sub_terms);
    // Check if this linear term is the constant term (a numerical value that
    // appears by itself).
    double temp_double;
    if (sub_terms.size() == 1 &&
        StringUtils::Stod(sub_terms[0], &temp_double)) {
      if (constant_term_present) {
         sprintf(error_msg,
                 "Found multiple constant terms in model. Aborting.");
         return false;
      }
      **constant_term = temp_double;
      constant_term_present = true;
      continue;
    }
    // Process linear term.
    for (const string& sub_term : sub_terms) {
      // Check if this sub_term can be parsed as a value, in which case it is
      // treated as the coefficient.
      double coefficient;
      if (StringUtils::Stod(sub_term, &coefficient)) {
        to_add.constant_ = coefficient;
        continue;
      }
      VariableTerm sub_term_to_add;
      if (!ParseLinearTerm(sub_term, &sub_term_to_add, error_msg)) return false;
      mult_terms.push_back(sub_term_to_add);
    }
    // Add linear term to linear_vars.
    to_add.terms_ = mult_terms;
    to_add.op_ = Operation::MULT;
    linear_vars->push_back(to_add);
  }

  // To indicate no constant term was found, set passed-in parameter to NULL.
  if (!constant_term_present) *constant_term = nullptr;

  return true;
}

bool ReadInput::ParseLinearTerm(
    const string& input, VariableTerm* term, char* error_msg) {
  // Sanity check input.
  if (term == nullptr) {
    sprintf(error_msg, "NULL term in call to ParseLinearTerm. Check API.");
    return false;
  }

  // Return false if input is empty.
  if (input.empty()) {
    sprintf(error_msg, "Empty input to ParseLinearTerm().");
    return false;
  }

  // Check for 'Log', 'LOG', or 'log' term.
  string stripped_input;
  if (StringUtils::StripPrefixString(input, "Log(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "LOG(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "log(", &stripped_input)) {
    // "Log" detected; ensure that anticipated opening parentheses '('
    // and closing parentheses ')' exist.
    const int input_length = stripped_input.length();
    if (input_length < 2 || stripped_input[input_length - 1] != ')') {
      sprintf(error_msg, "Detected 'Log' in linear term %s, but there is "
                         "no closing parenthese at the end of the term.",
              input.c_str());
      return false;
    }
    term->term_title_ = stripped_input.substr(0, input_length - 1);
    term->op_ = Operation::LOG;
    return true;
  }
    
  // Check for 'Exp', 'EXP', or 'exp' term.
  if (StringUtils::StripPrefixString(input, "Exp(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "EXP(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "exp(", &stripped_input)) {
    // "Exp" detected; ensure that anticipated opening parentheses '('
    // and closing parentheses ')' exist.
    const int input_length = stripped_input.length();
    if (input_length < 2 || stripped_input[input_length - 1] != ')') {
      sprintf(error_msg, "Detected 'exp(' in linear term %s, but there is "
                         "no closing parenthese at the end of the term.",
              input.c_str());
      return false;
    }
    term->term_title_ = stripped_input.substr(0, input_length - 1);
    term->op_ = Operation::EXP;
    return true;
  }

  // Check for 'Sqrt', 'SQRT', or 'sqrt' term.
  if (StringUtils::StripPrefixString(input, "Sqrt(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "SQRT(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "sqrt(", &stripped_input)) {
    // "Sqrt" detected; ensure that anticipated opening parentheses '('
    // and closing parentheses ')' exist.
    const int input_length = stripped_input.length();
    if (input_length < 2 || stripped_input[input_length - 1] != ')') {
      sprintf(error_msg, "Detected 'sqrt(' in linear term %s, but there is "
                         "no closing parenthese at the end of the term.",
              input.c_str());
      return false;
    }
    term->term_title_ = stripped_input.substr(0, input_length - 1);
    term->op_ = Operation::SQRT;
    return true;
  }

  // Check for 'Pow', 'POW', or 'pow' term.
  if (StringUtils::StripPrefixString(input, "Pow(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "POW(", &stripped_input) ||
      StringUtils::StripPrefixString(input, "pow(", &stripped_input)) {
    // "Exp" detected; ensure that anticipated opening parentheses '('
    // and closing parentheses ')' exist.
    const int input_length = stripped_input.length();
    if (input_length < 2 || stripped_input[input_length - 1] != ')') {
      sprintf(error_msg, "Detected 'pow(' in linear term %s, but there is "
                         "no closing parenthese at the end of the term.",
              input.c_str());
      return false;
    }
    size_t seperator = stripped_input.find(",");
    if (seperator == string::npos) {
      sprintf(error_msg, "Detected 'pow(' in linear term %s, but there is "
                         "no ',' marking the end of the variable name and "
                         "the desired exponent.\n Example usage: pow(X_1, "
                         "0.5), where X_1 is the variable name and 0.5 is "
                         "the desired exponent.", input.c_str());
      return false;
    }
    string exponent_str = stripped_input.substr(seperator + 1);
    // Remove trailing ')' from exponent.
    exponent_str = exponent_str.substr(0, exponent_str.length() - 1);
    double exponent;
    if (exponent_str.empty() || !StringUtils::Stod(exponent_str, &exponent)) {
      sprintf(error_msg, "Detected 'pow(' in linear term %s, but the text "
                         "following the ',' (which should be the exponent) "
                         "cannot be parsed as a numeric value.\nExample "
                         "usage: pow(X_1, 0.5), where X_1 is the variable "
                         "name and 0.5 is the desired exponent.",
              input.c_str());
      return false;
    }
    term->term_title_ = stripped_input.substr(0, seperator);
    term->op_ = Operation::POW;
    term->exponent_ = exponent;
    return true;
  }

  // No 'Exp', 'Log', 'Sqrt', or 'Pow' detected.
  // Interpret entire string as the term.
  term->term_title_ = input;
  term->op_ = Operation::IDENTITY;
  return true;
}

bool ReadInput::ParseDependentTerm(
    const string& input,
    vector<VariableTerm>* dep_term, char* error_msg) {
  // Sanity check input.
  if (dep_term == nullptr) {
    sprintf(error_msg, "NULL term in call to ParseLinearTerm. Check API.");
    return false;
  }

  // Return false if input is empty.
  if (input.empty()) {
    sprintf(error_msg, "Empty input to ParseLinearTerm().");
    return false;
  }

  vector<string> parts;
  StringUtils::Split(input, ",", &parts);

  if (parts.size() != 2 && parts.size() != 3) {
    sprintf(error_msg, "Unable to parse LHS of Regression equation: %s",
            input.c_str());
    return false;
  }

  // Process first dep variable.
  string dep_var, stripped_dep_var;
  StringUtils::StripPrefixString(parts[0], "(", &dep_var);
  StringUtils::RemoveExtraWhitespace(dep_var, &stripped_dep_var);
  VariableTerm term;
  if (!ParseLinearTerm(stripped_dep_var, &term, error_msg)) return false;
  dep_term->clear();
  dep_term->push_back(term);

  // Process second dep variable (if present).
  int status_var_index = 1;
  if (parts.size() == 3) {
    string dep_var_two;
    StringUtils::RemoveExtraWhitespace(parts[1], &dep_var_two);
    VariableTerm term_two;
    if (!ParseLinearTerm(dep_var_two, &term_two, error_msg)) return false;
    dep_term->push_back(term_two);
    status_var_index++;
  }

  // Process last dep variable (Status).
  string status_var, stripped_status_var;
  StringUtils::StripSuffixString(parts[status_var_index], ")", &status_var);
  StringUtils::RemoveExtraWhitespace(status_var, &stripped_status_var);
  VariableTerm status_term;
  status_term.term_title_ = stripped_status_var;
  status_term.op_ = Operation::IDENTITY;
  dep_term->push_back(status_term);
  return true;
}

bool ReadInput::ProcessUserEnteredNominalColumns(
    const string& user_input, const vector<string>& variable_names,
    set<int>* nominal_columns, char* error_msg) {
  // This should never happen based on current calls to this function.
  if (user_input.empty()) return false;

  // Check the special tokens that represent the user has indicated
  // that all input variables are numeric. In this case, nothing to do.
  if (user_input == "NONE" || user_input == "None" || user_input == "none") {
    return true;
  }
  
  // Sanity check inputs.
  if (variable_names.empty()) {
    sprintf(error_msg, "No variable names found. Check that the first line "
                       "in your input file is appropriately formatted.");
    return false;
  }
  if (nominal_columns == nullptr) {
    sprintf(error_msg, "NULL 'nominal_columns'. Check API in call to "
                       "ProcessUserEnteredNominalColumns.");
    return false;
  }

  // Clear nominal_columns.
  nominal_columns->clear();

  // Parse user_input.
  vector<string> new_nominal_columns;
  if (!StringUtils::Split(user_input, ", ", &new_nominal_columns)) {
    sprintf(error_msg, "Unreadable input.");
    return false;
  }

  // Iterate through the user-input column names, and match them
  // to the column indices of the header line of the input file.
  for (int i = 0; i < new_nominal_columns.size(); ++i) {
    bool found_match = false;
    for (int j = 0; j < variable_names.size(); ++j) {
      if (new_nominal_columns[i] == variable_names[j]) {
        nominal_columns->insert(j);
        found_match = true;
        break;
      }
    }
    if (!found_match) {
      sprintf(error_msg, "Unable to find '%s' as one of the variable "
                         "names in the header line of the input file.",
              new_nominal_columns[i].c_str());
      return false;
    }
  }

  return true;
}

}  // namespace file_reader_utils
