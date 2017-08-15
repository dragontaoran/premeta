// Date: June 2015
// Author: paulbunn@email.unc.edu (Paul Bunn)

#include "write_seqmeta_utils.h"

#include "MapUtils/map_utils.h"
#include "MathUtils/number_comparison.h"
#include "PreMeta/premeta_constants.h"
#include "PreMeta/premeta_structures.h"
#include "PreMeta/read_r_binary_utils.h"
#include "PreMeta/write_r_binary_utils.h"

#include <cfloat>   // For DBL_MIN
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace map_utils;
using namespace math_utils;
using namespace std;

namespace premeta {

bool ComputeDsCMatrixFromCovarianceMatrix(
    const int study_num, const double& sigma_sq, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<pair<Position, Position>, double>& covariances,
    vector<int>* i_values, vector<int>* p_values, vector<double>* x_values) {
  const int dim = gene_positions.size();
  int nonzero_values_counter = 0;
  // For some reason, the p-vector always starts off with a '0' in the
  // first coordinate.
  p_values->push_back(0);
  for (int i = 0; i < dim; ++i) {
    const Position& row = gene_positions[i];
    const bool row_is_swapped =
        IsSnpRefAltSwapped(study_num, row, snp_to_ref_alt_and_study);
    for (int j = 0; j <= i; ++j) {
      const Position& col = gene_positions[j];
      const bool col_is_swapped =
          IsSnpRefAltSwapped(study_num, col, snp_to_ref_alt_and_study);
      const double swap_multiplier =
          (row_is_swapped && col_is_swapped) ? 1.0 :
          (row_is_swapped || col_is_swapped) ? -1.0 : 1.0;
      double covariance = DBL_MIN;
      if (monomorphic_snps.find(row) != monomorphic_snps.end() ||
          monomorphic_snps.find(col) != monomorphic_snps.end()) {
        covariance = 0.0;
      } else {
        if (!LookupCovariance(row, col, covariances, &covariance) ||
            covariance == DBL_MIN) {
          cout << "ERROR: Unable to construct covariance matrix for SeqMeta. "
               << "Couldn't find covariance for Positions (" << PrintPosition(row)
               << ", " << PrintPosition(col) << ")." << endl;
          return false;
        }
      }
      if (covariance == 0.0) {
        continue;
      }
      i_values->push_back(j);
      nonzero_values_counter++;
      x_values->push_back(swap_multiplier * sigma_sq * covariance / (rescale * rescale));
    }
    p_values->push_back(nonzero_values_counter);
  }
  return true;
}

bool ConstructSeqMetaTreeGeneCovI(
    vector<int>* i_values, vector<string*>* tags, RDataPairList* i) {
  i->tag_ = new string("i");
  AppendTagIfNotPresent(i->tag_, tags);
  i->obj_ = new RDataNode();
  i->obj_->obj_ = new RDataObject();
  i->obj_->obj_->int_vec_ = i_values;
  return true;
}

bool ConstructSeqMetaTreeGeneCovP(
    vector<int>* p_values, vector<string*>* tags, RDataPairList* p) {
  p->tag_ = new string("p");
  AppendTagIfNotPresent(p->tag_, tags);
  p->obj_ = new RDataNode();
  p->obj_->obj_ = new RDataObject();
  p->obj_->obj_->int_vec_ = p_values;
  return true;
}

bool ConstructSeqMetaTreeGeneCovDim(
    const int num_positions, vector<string*>* tags, RDataPairList* dim) {
  dim->tag_ = new string("Dim");
  AppendTagIfNotPresent(dim->tag_, tags);
  dim->obj_ = new RDataNode();
  dim->obj_->obj_ = new RDataObject();
  dim->obj_->obj_->int_vec_ = new vector<int>();
  dim->obj_->obj_->int_vec_->push_back(num_positions);
  dim->obj_->obj_->int_vec_->push_back(num_positions);
  return true;
}

bool ConstructSeqMetaTreeGeneCovDimnames(
    const vector<Position>& gene_positions,
    vector<string*>* tags, RDataPairList* dimnames) {
  dimnames->tag_ = new string("Dimnames");
  AppendTagIfNotPresent(dimnames->tag_, tags);
  dimnames->obj_ = new RDataNode();
  dimnames->obj_->obj_ = new RDataObject();
  dimnames->obj_->obj_->list_vec_ = new vector<RDataNode*>();
  
  RDataNode* rownames = new RDataNode();
  rownames->obj_ = new RDataObject();
  rownames->obj_->str_vec_ = new vector<string>();
  for (const Position& pos : gene_positions) {
    rownames->obj_->str_vec_->push_back(PrintPosition(pos));
  }
  dimnames->obj_->obj_->list_vec_->push_back(rownames);

  RDataNode* colnames = new RDataNode();
  colnames->obj_ = new RDataObject();
  colnames->obj_->str_vec_ = new vector<string>();
  *(colnames->obj_->str_vec_) = *(rownames->obj_->str_vec_);
  dimnames->obj_->obj_->list_vec_->push_back(colnames);
  return true;
}

bool ConstructSeqMetaTreeGeneCovX(
    vector<double>* x_values, vector<string*>* tags, RDataPairList* x) {
  x->tag_ = new string("x");
  AppendTagIfNotPresent(x->tag_, tags);
  x->obj_ = new RDataNode();
  x->obj_->obj_ = new RDataObject();
  x->obj_->obj_->real_vec_ = x_values;
  return true;
}

bool ConstructSeqMetaTreeGeneCovUplo(
    vector<string*>* tags, RDataPairList* uplo) {
  uplo->tag_ = new string("uplo");
  AppendTagIfNotPresent(uplo->tag_, tags);
  uplo->obj_ = new RDataNode();
  uplo->obj_->obj_ = new RDataObject();
  uplo->obj_->obj_->str_vec_ = new vector<string>();
  uplo->obj_->obj_->str_vec_->push_back("U");
  return true;
}

bool ConstructSeqMetaTreeGeneCovFactors(
    vector<string*>* tags, RDataPairList* factors) {
  factors->tag_ = new string("factors");
  AppendTagIfNotPresent(factors->tag_, tags);
  // Create an empty list, which is the Object for 'factors'.
  factors->obj_ = new RDataNode();
  factors->obj_->obj_ = new RDataObject();
  factors->obj_->obj_->list_vec_ = new vector<RDataNode*>();
  return true;
}

bool ConstructSeqMetaTreeGeneCovClass(
    vector<string*>* tags, RDataPairList* class_name) {
  // Class Tag.
  class_name->tag_ = new string("class");
  AppendTagIfNotPresent(class_name->tag_, tags);

  // Class Object.
  class_name->obj_ = new RDataNode();
  class_name->obj_->obj_ = new RDataObject();
  class_name->obj_->obj_->str_vec_ = new vector<string>();
  class_name->obj_->obj_->str_vec_->push_back("dsCMatrix");

  // Class Attribute.
  class_name->obj_->attr_ = new RDataNode();
  class_name->obj_->attr_->obj_ = new RDataObject();
  class_name->obj_->attr_->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  RDataPairList* package = new RDataPairList();
  package->tag_ = new string("package");
  AppendTagIfNotPresent(package->tag_, tags);
  package->obj_ = new RDataNode();
  package->obj_->obj_ = new RDataObject();
  package->obj_->obj_->str_vec_ = new vector<string>();
  package->obj_->obj_->str_vec_->push_back("Matrix");
  class_name->obj_->attr_->obj_->pair_list_vec_->push_back(package);

  return true;
}

bool ConstructSeqMetaTreeGeneScores(
    const int study_num, const double& sey, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    vector<string*>* tags, RDataNode* gene_scores) {
  // Parse Scores Object.
  gene_scores->obj_ = new RDataObject();
  gene_scores->obj_->real_vec_ = new vector<double>();
  // Will also populate the Score Attribute names in this loop, so we don't
  // have to iterate through the loop a second time below when processing
  // the Attribute.
  vector<string>* score_attr_names = new vector<string>();
  const double sigma_sq = sey * sey;
  for (const Position& pos : gene_positions) {
    score_attr_names->push_back(PrintPosition(pos));
    map<Position, SnpInfo>::const_iterator snp_itr = snp_info.find(pos);
    if (snp_itr == snp_info.end()) {
      cout << "ERROR: Unable to find score for Position "
           << PrintPosition(pos) << endl;
      return false;
    }
    if (monomorphic_snps.find(pos) != monomorphic_snps.end()) {
      gene_scores->obj_->real_vec_->push_back(0.0);
    } else if (IsSnpRefAltSwapped(study_num, pos, snp_to_ref_alt_and_study)) {
      gene_scores->obj_->real_vec_->push_back(
          -1.0 * (snp_itr->second.u_stat_) * sigma_sq / rescale);
    } else {
      gene_scores->obj_->real_vec_->push_back(
          (snp_itr->second.u_stat_) * sigma_sq / rescale);
    }
  }

  // Parse Scores Attribute.
  gene_scores->attr_ = new RDataNode();
  gene_scores->attr_->obj_ = new RDataObject();
  gene_scores->attr_->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  RDataPairList* scores_attr = new RDataPairList();
  scores_attr->tag_ = new string("names");
  AppendTagIfNotPresent(scores_attr->tag_, tags);
  scores_attr->obj_ = new RDataNode();
  scores_attr->obj_->obj_ = new RDataObject();
  scores_attr->obj_->obj_->str_vec_ = score_attr_names;
  gene_scores->attr_->obj_->pair_list_vec_->push_back(scores_attr);

  return true;
}

bool ConstructSeqMetaTreeGeneCov(
    const int study_num, const double& sigma_sq, const double& rescale,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* gene_cov) {
  gene_cov->obj_ = new RDataObject();
  gene_cov->obj_->class_ = new RDataNode();
  gene_cov->obj_->class_->is_object_ = true;
  gene_cov->obj_->class_->obj_ = new RDataObject();
  gene_cov->obj_->class_->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  vector<RDataPairList*>* pair_list_vec =
      gene_cov->obj_->class_->obj_->pair_list_vec_;

  // Construct the (i, p, x) fields of dsCMatrix class from 'covariances'.
  vector<int>* i_values = new vector<int>();
  vector<int>* p_values = new vector<int>();
  vector<double>* x_values = new vector<double>();
  if (!ComputeDsCMatrixFromCovarianceMatrix(
        study_num, sigma_sq, rescale, gene_positions, monomorphic_snps,
        snp_to_ref_alt_and_study, covariances, i_values, p_values, x_values)) {
    return false;
  }

  // Parse gene$cov$i.
  RDataPairList* i = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovI(i_values, tags, i)) {
    return false;
  }
  pair_list_vec->push_back(i);

  // Parse gene$cov$p.
  RDataPairList* p = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovP(p_values, tags, p)) {
    return false;
  }
  pair_list_vec->push_back(p);

  // Parse gene$cov$Dim.
  RDataPairList* dim = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovDim(gene_positions.size(), tags, dim)) {
    return false;
  }
  pair_list_vec->push_back(dim);

  // Parse gene$cov$Dimnames.
  RDataPairList* dimnames = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovDimnames(gene_positions, tags, dimnames)) {
    return false;
  }
  pair_list_vec->push_back(dimnames);

  // Parse gene$cov$x.
  RDataPairList* x = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovX(x_values, tags, x)) {
    return false;
  }
  pair_list_vec->push_back(x);

  // Parse gene$cov$uplo.
  RDataPairList* uplo = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovUplo(tags, uplo)) {
    return false;
  }
  pair_list_vec->push_back(uplo);

  // Parse gene$cov$factors.
  RDataPairList* factors = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovFactors(tags, factors)) {
    return false;
  }
  pair_list_vec->push_back(factors);

  // Parse gene$cov$class.
  RDataPairList* class_name = new RDataPairList();
  if (!ConstructSeqMetaTreeGeneCovClass(tags, class_name)) {
    return false;
  }
  pair_list_vec->push_back(class_name);

  return true;
}

bool ConstructSeqMetaTreeGeneN(
    const int num_samples, RDataNode* gene_n) {
  gene_n->obj_ = new RDataObject();
  gene_n->obj_->int_vec_ = new vector<int>();
  gene_n->obj_->int_vec_->push_back(num_samples);
  return true;
}

bool ConstructSeqMetaTreeGeneMaf(
    const int study_num,
    const vector<Position>& gene_positions,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    vector<string*>* tags, RDataNode* gene_maf) {
  // Parse MAF Object.
  gene_maf->obj_ = new RDataObject();
  gene_maf->obj_->real_vec_ = new vector<double>();
  // Will also populate the MAF Attribute names in this loop, so we don't
  // have to iterate through the loop a second time below when processing
  // the Attribute.
  vector<string>* maf_attr_names = new vector<string>();
  for (const Position& pos : gene_positions) {
    maf_attr_names->push_back(PrintPosition(pos));
    map<Position, SnpInfo>::const_iterator snp_itr = snp_info.find(pos);
    if (snp_itr == snp_info.end()) {
      cout << "ERROR: Unable to find maf for Position "
           << PrintPosition(pos) << endl;
      return false;
    }
    if (IsSnpRefAltSwapped(study_num, pos, snp_to_ref_alt_and_study)) {
      gene_maf->obj_->real_vec_->push_back(1.0 - snp_itr->second.maf_);
    } else {
      gene_maf->obj_->real_vec_->push_back(snp_itr->second.maf_);
    }
  }

  // Parse MAF Attribute.
  gene_maf->attr_ = new RDataNode();
  gene_maf->attr_->obj_ = new RDataObject();
  gene_maf->attr_->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  RDataPairList* maf_attr = new RDataPairList();
  maf_attr->tag_ = new string("names");
  AppendTagIfNotPresent(maf_attr->tag_, tags);
  maf_attr->obj_ = new RDataNode();
  maf_attr->obj_->obj_ = new RDataObject();
  maf_attr->obj_->obj_->str_vec_ = maf_attr_names;
  gene_maf->attr_->obj_->pair_list_vec_->push_back(maf_attr);

  return true;
}

bool ConstructSeqMetaTreeGeneSey(
    const double& sey, RDataNode* gene_sey) {
  gene_sey->obj_ = new RDataObject();
  gene_sey->obj_->real_vec_ = new vector<double>();
  gene_sey->obj_->real_vec_->push_back(sey);
  return true;
}

bool ConstructSeqMetaTreeGene(
    const int study_num, const int num_samples,
    const double& rescale, const double& sey,
    const vector<Position>& gene_positions,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* gene_info) {
  // Parse Gene's Object: List of 5 items:
  //   {$scores, $cov, $n, $maf, $sey}
  gene_info->obj_ = new RDataObject();
  gene_info->obj_->list_vec_ = new vector<RDataNode*>();
  
  // Parse gene$scores.
  RDataNode* gene_scores = new RDataNode();
  if (!ConstructSeqMetaTreeGeneScores(
          study_num, sey, rescale, gene_positions, monomorphic_snps,
          snp_to_ref_alt_and_study, snp_info, tags, gene_scores)) {
    return false;
  }
  gene_info->obj_->list_vec_->push_back(gene_scores);

  // Parse gene$cov.
  RDataNode* gene_cov = new RDataNode();
  if (!ConstructSeqMetaTreeGeneCov(
          study_num, sey * sey, rescale, gene_positions, monomorphic_snps,
          snp_to_ref_alt_and_study, covariances, tags, gene_cov)) {
    return false;
  }
  gene_info->obj_->list_vec_->push_back(gene_cov);

  // Parse gene$n.
  RDataNode* gene_n = new RDataNode();
  if (!ConstructSeqMetaTreeGeneN(num_samples, gene_n)) {
    return false;
  }
  gene_info->obj_->list_vec_->push_back(gene_n);

  // Parse gene$maf.
  RDataNode* gene_maf = new RDataNode();
  if (!ConstructSeqMetaTreeGeneMaf(
          study_num, gene_positions, snp_to_ref_alt_and_study, snp_info, tags, gene_maf)) {
    return false;
  }
  gene_info->obj_->list_vec_->push_back(gene_maf);

  // Parse gene$sey
  RDataNode* gene_sey = new RDataNode();
  if (!ConstructSeqMetaTreeGeneSey(sey, gene_sey)) {
    return false;
  }
  gene_info->obj_->list_vec_->push_back(gene_sey);

  // Parse Gene's Attribute.
  gene_info->attr_ = new RDataNode();
  gene_info->attr_->obj_ = new RDataObject();
  gene_info->attr_->obj_->pair_list_vec_ = new vector<RDataPairList*>();
  RDataPairList* gene_info_attr = new RDataPairList();
  gene_info_attr->tag_ = new string("names");
  AppendTagIfNotPresent(gene_info_attr->tag_, tags);
  gene_info_attr->obj_ = new RDataNode();
  gene_info_attr->obj_->obj_ = new RDataObject();
  gene_info_attr->obj_->obj_->str_vec_ = new vector<string>();
  gene_info_attr->obj_->obj_->str_vec_->push_back("scores");
  gene_info_attr->obj_->obj_->str_vec_->push_back("cov");
  gene_info_attr->obj_->obj_->str_vec_->push_back("n");
  gene_info_attr->obj_->obj_->str_vec_->push_back("maf");
  gene_info_attr->obj_->obj_->str_vec_->push_back("sey");
  gene_info->attr_->obj_->pair_list_vec_->push_back(gene_info_attr);
  return true;
}

bool ConstructSeqMetaTreeRootObject(
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* root) {
  // Set root's is_object_ field.
  root->is_object_ = true;

  // Set root's attr_ field.
  root->attr_ = new RDataNode();
  if (!ConstructSeqMetaTreeRootAttribute(
        gene_to_snp_positions, tags, root->attr_)) {
    return false;
  }

  // Set root's obj_ field: A list, with one item per gene.
  root->obj_ = new RDataObject();
  RDataObject* obj = root->obj_;
  obj->list_vec_ = new vector<RDataNode*>();
  const double sey = pow(sigma_sq, 0.5);
  for (map<string, vector<Position>>::const_iterator gene_itr =
           gene_to_snp_positions.begin();
       gene_itr != gene_to_snp_positions.end(); ++gene_itr) {
    const vector<Position>& gene_positions = gene_itr->second;
    RDataNode* gene_info = new RDataNode();
    if (!ConstructSeqMetaTreeGene(
            study_num, num_samples, rescale, sey, gene_positions,
            monomorphic_snps, snp_to_ref_alt_and_study,
            snp_info, covariances, tags, gene_info)) {
      cout << "ERROR in generating seqMeta output for gene "
           << gene_itr->first << "'. Aborting." << endl;
      return false;
    }
    obj->list_vec_->push_back(gene_info);
  }
  return true;
}

bool ConstructSeqMetaTreeRootAttribute(
    const map<string, vector<Position>>& gene_to_snp_positions,
    vector<string*>* tags, RDataNode* root_attr) {
  root_attr->obj_ = new RDataObject();
  root_attr->obj_->pair_list_vec_ = new vector<RDataPairList*>();

  // First Attribute item is 'dim', which gives number of genes.
  RDataPairList* dim = new RDataPairList();
  dim->tag_ = new string("dim");
  AppendTagIfNotPresent(dim->tag_, tags);
  dim->obj_ = new RDataNode();
  dim->obj_->obj_ = new RDataObject();
  dim->obj_->obj_->int_vec_ = new vector<int>();
  dim->obj_->obj_->int_vec_->push_back(gene_to_snp_positions.size());
  root_attr->obj_->pair_list_vec_->push_back(dim);

  // Second Attribute item is 'dimnames'.
  RDataPairList* dimnames = new RDataPairList();
  dimnames->tag_ = new string("dimnames");
  AppendTagIfNotPresent(dimnames->tag_, tags);
  dimnames->obj_ = new RDataNode();
  dimnames->obj_->obj_ = new RDataObject();
  dimnames->obj_->obj_->list_vec_ = new vector<RDataNode*>();
  RDataNode* gene_names = new RDataNode();
  gene_names->obj_ = new RDataObject();
  gene_names->obj_->str_vec_ = new vector<string>();
  for (const string& gene_name : Keys(gene_to_snp_positions)) {
    gene_names->obj_->str_vec_->push_back(gene_name);
  }
  dimnames->obj_->obj_->list_vec_->push_back(gene_names);
  root_attr->obj_->pair_list_vec_->push_back(dimnames);

  // Third Attribute item is 'family'.
  RDataPairList* family = new RDataPairList();
  family->tag_ = new string("family");
  AppendTagIfNotPresent(family->tag_, tags);
  family->obj_ = new RDataNode();
  family->obj_->obj_ = new RDataObject();
  family->obj_->obj_->str_vec_ = new vector<string>();
  family->obj_->obj_->str_vec_->push_back("gaussian");
  root_attr->obj_->pair_list_vec_->push_back(family);

  // Fourth Attribute item is 'class'.
  RDataPairList* class_name = new RDataPairList();
  class_name->tag_ = new string("class");
  AppendTagIfNotPresent(class_name->tag_, tags);
  class_name->obj_ = new RDataNode();
  class_name->obj_->obj_ = new RDataObject();
  class_name->obj_->obj_->str_vec_ = new vector<string>();
  class_name->obj_->obj_->str_vec_->push_back("seqMeta");
  root_attr->obj_->pair_list_vec_->push_back(class_name);

  return true;
}

bool ConstructSeqMetaTree(
    const string& object_name,
    const int study_num, const int num_samples,
    const double& rescale, const double& sigma_sq,
    const set<Position>& monomorphic_snps,
    const map<Position, tuple<string, string, set<int>>>& snp_to_ref_alt_and_study,
    const map<string, vector<Position>>& gene_to_snp_positions,
    const map<Position, SnpInfo>& snp_info,
    const map<pair<Position, Position>, double>& covariances,
    vector<string*>* tags, RDataNode* root) {
  // Set root's obj_ field.
  root->obj_ = new RDataObject();
  root->obj_->pair_list_vec_ = new vector<RDataPairList*>();

  RDataPairList* root_pair_list = new RDataPairList();
  root_pair_list->tag_ = new string(object_name);
  AppendTagIfNotPresent(root_pair_list->tag_, tags);
  root_pair_list->obj_ = new RDataNode();
  root->obj_->pair_list_vec_->push_back(root_pair_list);

  // Construct Root Pair-List's Object.
  if (!ConstructSeqMetaTreeRootObject(
        study_num, num_samples, rescale, sigma_sq, monomorphic_snps,
        snp_to_ref_alt_and_study, gene_to_snp_positions,
        snp_info, covariances, tags, root_pair_list->obj_)) {
    return false;
  }

  return true;
}

bool PrintToSeqMeta(
    const RDataNode& root, const string& out_file) {
  
  ofstream outfile(out_file.c_str(), ofstream::binary);
  if (!outfile.is_open()) {
    cout << "ERROR in trying to write file '" << out_file
         << "': Unable to open file. Make sure directory exists." << endl;
    return false;
  }

  WriteRDataHeader(outfile);

  vector<string> tags;
  WriteRDataTree(&tags, root, outfile);

  outfile.close();
  return true;
}

}  // namespace premeta
