/* Copyright 2016-2018 Lingfei Wang
 * 
 * This file is part of Findr.
 * 
 * Findr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Findr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with Findr.  If not, see <http://www.gnu.org/licenses/>.
 */
/* This part contains the main interface function of genotype assisted pij inference.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_H_
#define _HEADER_LIB_PIJ_GASSIST_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Estimates the p-value of A  B against A->B from genotype and expression data with 5 tests.
 * E is always the best eQTL of A. Full data is required.
 * g:	(ng,ns) Genotype data, =0,1,...,nv-1. Each is the best eQTL of the corresponding gene in t.
 * t:	(ng,ns) Expression data of A.
 * t2:	(nt,ns) Expression data of B. Can be A or a superset of A.
 * p1:	(ng) P-values of step 1. Tests E->A v.s. E  A.
 * p2:	(ng,nt) P-values of step 2. Tests E->B v.s. E  B.
 * p3:	(ng,nt) P-values of step 3. Tests E->A->B v.s. E->A->B with E->B.
 * p4:	(ng,nt) P-values of step 4. Tests E->A->B with E->B v.s. E->A  B.
 * p5:	(ng,nt) P-values of step 5. Tests E->A->B with E->B v.s. A<-E->B.
 * nv:	Number of possible values each genotype entry may take, =number of alleles+1.
 * memlimit:	The function is able to split very large datasets (ng and nt) into smaller chunks for inference. This variable specifies the approximate memory usage limit. Note: For large datasets, a too small memory limit can fail the function. For unlimited memory, set memlimit=-1.
 * Return:	0 on sucess
 * Appendix:
 * 		ng:	Number of genes with best eQTL.
 * 		nt:	Number of genes with expression data for B
 * 		ns:	Number of samples.
 */
int pijs_gassist_pv(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,size_t memlimit);

/* Estimates the probability of A->B from genotype and expression data with 5 tests.
 * E is always the best eQTL of A. Full data is required.
 * g:	(ng,ns) Genotype data, =0,1,...,nv-1. Each is the best eQTL of the corresponding gene in t.
 * t:	(ng,ns) Expression data of A.
 * t2:	(nt,ns) Expression data of B. Can be A or a superset of A.
 * p1:	(ng) Probabilities of step 1. Tests E->A v.s. E  A. For nodiag=0, because the function expects significant eQTLs, p1 always return 1. For nodiag=1, uses diagonal elements of p2. Consider replacing p1 with your own (1-FDR) from eQTL discovery.
 * p2:	(ng,nt) Probabilities of step 2. Tests E->B v.s. E  B.
 * p3:	(ng,nt) Probabilities of step 3. Tests E->A->B v.s. E->A->B with E->B.
 * p4:	(ng,nt) Probabilities of step 4. Tests E->A->B with E->B v.s. E->A  B.
 * p5:	(ng,nt) Probabilities of step 5. Tests E->A->B with E->B v.s. A<-E->B.
 * nv:	Number of possible values each genotype entry may take, =number of alleles+1.
 * nodiag:	When the top ng rows of t2 is exactly t, diagonals of p2 and p3 are meaningless. In this case, set nodiag to 1 to avoid inclusion of NANs. For nodiag=0, t and t2 should not have any identical genes.
 * memlimit:	The function is able to split very large datasets (ng and nt) into smaller chunks for inference. This variable specifies the approximate memory usage limit. Note: For large datasets, a too small memory limit can fail the function. For unlimited memory, set memlimit=-1.
 * Return:	0 on sucess
 * Appendix:
 * 		ng:	Number of genes with best eQTL.
 * 		nt:	Number of genes with expression data for B
 * 		ns:	Number of samples.
 */
int pijs_gassist(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag,size_t memlimit);

/* Estimates the probability of A->B from genotype and expression data with defaults combination of tests. Uses results from pijs_gassist. Variables have the same definitions except:
 * ans:	(ng,nt) Predicted probability of A->B based on default combination of 5 tests. The default combination is (p2*p5+p4)/2. Note: this combination does not include p1.
 * Return:	0 on sucess
 */
int pij_gassist(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag,size_t memlimit);

/* Estimates the probability of A->B from genotype and expression data with traditional causal inference method.
 * NOTE:	This is not and is not intended as a loyal reimplementation of the Trigger R package. Instead, it aims at reusing methods and tests of Findr to produce inferences that mimicks the three tests performed by Trigger. Many implementational details are different between this function and Trigger, althrough a significant (but not full) overlap has been observed in existing studies. This method does not include p1.
 * Inputs and ouputs are the same as function pij_gassist.
 */
int pij_gassist_trad(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag,size_t memlimit);

#ifdef __cplusplus
}
#endif
#endif

























