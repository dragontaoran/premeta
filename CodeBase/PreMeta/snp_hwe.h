// Aurthor: Jan Wigginton
//
// Description: 
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  

#ifndef SNP_HWE_H
#define SNP_HWE_H

namespace premeta {

// Computes and returns HWE from N_HET, N_REF, and N_ALT values.
double GetSnpHWE(const int n_het, const int n_ref, const int n_alt);

}  // namespace premeta
#endif
