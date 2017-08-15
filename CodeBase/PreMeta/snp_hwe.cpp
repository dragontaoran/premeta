#include "snp_hwe.h"

#include <cstdlib>
#include <iostream>

using namespace std;

namespace premeta {

double GetSnpHWE(const int n_het, const int n_ref, const int n_alt) {
  if (n_ref < 0 || n_alt < 0 || n_het < 0) {
    cout << "ERROR in computing HWE: all counts should be positive. n_het: "
         << n_het << ", n_ref: " << n_ref << ", n_alt: " << n_alt << endl;
    return -1.0;
  }

  const int obs_homc = n_ref < n_alt ? n_alt : n_ref;
  const int obs_homr = n_ref < n_alt ? n_ref : n_alt;
  const int rare_copies = 2 * obs_homr + n_het;
  const int genotypes = n_het + obs_homc + obs_homr;

  double* het_probs =
      (double*) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == nullptr) {
    cout << "ERROR in computing HWE: Unable to allocate array for heterozygote "
         << " probabilities." << endl;
    return -1.0;
  }
  
  // Initialize het_probs to all be 0.0.
  for (int i = 0; i <= rare_copies; i++) {
     het_probs[i] = 0.0;
  }

  // Start at midpoint.
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

  // Check to ensure that midpoint and rare alleles have same parity.
  if ((rare_copies & 1) ^ (mid & 1)) {
     mid++;
  }

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] =
        het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) /
        (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    // 2 fewer heterozygotes for next iteration ->
    //   add one rare, one common homozygote.
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] =
        het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /
        ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    // Add 2 heterozygotes for next iteration ->
    //   subtract one rare, one common homozygote.
    curr_homr--;
    curr_homc--;
  }

  for (int i = 0; i <= rare_copies; i++) {
    het_probs[i] /= sum;
  }

  /*
  // Alternate p-value calculation for p_hi/p_lo
  double p_hi = het_probs[n_het];
  for (i = n_het + 1; i <= rare_copies; i++)
    p_hi += het_probs[i];
  
  double p_lo = het_probs[n_het];
  for (i = n_het - 1; i >= 0; i--)
     p_lo += het_probs[i];

  
  double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */

  double p_hwe = 0.0;
  // p-value calculation for p_hwe.
  for (int i = 0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[n_het]) continue;
    p_hwe += het_probs[i];
  }
  
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}

}  // namespace premeta
