/* Copyright 2016 Lingfei Wang
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
/* This part contains the log likelihood ratio calculations.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_LLR_H_
#define _HEADER_LIB_PIJ_GASSIST_LLR_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

//Struct of info pack for each llr thread
struct pij_gassist_llr_threadinfo{
	const MATRIXG*	g;				//genotype matrix
	const MATRIXF*	t;				//transcript matrix for A
	const MATRIXF*	t2;				//transcript matrix for B
	size_t	nv;				//Number of possible values for each genotype
	const VECTORF*	vbuff1;			//Vector (ns) buffer of all 1.
	//(ng) Output for Log likelihood calculation step 1
	VECTORF			*llr1;
	//(ng,nt) Output for Log likelihood calculation step 2,3
	MATRIXF			*llr2b,*llr2c,*llr3;	
};

/* Calculates the ratio and mean of all transcripts (t) for genes (g) with existing buffers.
 * g:		MATRIXG (ng,ns) genotype data, for multiple SNP and samples. Each element takes the value 0 to nv-1
 * t1:		MATRIXF (ng,ns) of transcript A data.
 * t2:		MATRIXF (ng,ns) of transcript B data.
 * mratio:	MATRIXF (ng,nv) of categorical ratio of samples for each SNP type. For return purpose.
 * mmean1:	MATRIXF (ng,nv) of categorical means of corresponding transcript A. For return purpose.
 * mmean2:	MATRIXF (ng,nv) of categorical means of corresponding transcript B. For return purpose.
 * nv:		number of possible values of g.
 */
void pij_gassist_llr_ratioandmean_1v1(const MATRIXG* g,const MATRIXF* t1,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF* mmean2,size_t nv);

/* Calculates the ratio and mean of all transcripts (t) for genes (g) with existing buffers.
 * g:		MATRIXG (ng,ns) genotype data, for multiple SNP and samples. Each element takes the value 0 to nv-1
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns) of transcript data for B
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For return purpose.
 * mmean1:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * nv:		number of possible values of g.
 * mb1:		MATRIXF[nv] (ng,ns) Buffer matrix
 * vb:		buffer. const VECTORF (ns). Must be set to 1 for all elements.
 */
void pij_gassist_llr_ratioandmean_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,MATRIXF** mb1,const VECTORF* vb);

/* Calculates the ratio and mean of all transcripts (t) for genes (g) with existing buffers.
 * g:		MATRIXG (ng,ns) genotype data, for multiple SNP and samples. Each element takes the value 0 to nv-1
 * t:		MATRIXF (nt,ns) of transcript data.
 * t2:		MATRIXF (nt,ns) of transcript data for B
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For return purpose.
 * mmean1:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * nv:		number of possible values of g.
 * vb:		buffer. const VECTORF (ns). Must be set to 1 for all elements.
 * Return:	0 on success
 */
int pij_gassist_llr_ratioandmean(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,const VECTORF* vb);

//Struct for parameters of pij_gassist_llr_block_buffed (below)
struct pij_gassist_llr_block_buffed_params{
	//Number of (E,A) pairs
	size_t			ng;
	//Number of possible genotypes per SNP
	size_t			nv;
	//(nv,ng) Genotype ratio matrix
	const MATRIXF*	mratio;
	//(nv,ng) Mean per genotype matrix
	const MATRIXF*	mmean1;
	//[nv](ng,nt) Mean per genotype matrix
	const MATRIXF**	mmean2;
	//(ng,nt) Correlation matrix
	MATRIXF*		mcov;		
	//(ng) Output vector for log likelihood ratio 1
	VECTORF* 		llr1;		
	//(ng,nt) Output matrix for log likelihood ratio 2 bold
	MATRIXF*		llr2b;		
	//(ng,nt) Output matrix for log likelihood ratio 2 conservative
	MATRIXF*		llr2c;		
	//(ng,nt) Output matrix for log likelihood ratio 3
	MATRIXF*		llr3;		
	//[nv](ng,nt) Buffer matrix
	MATRIXF**		mb1;		
};

/* Calculates the 3 log likelihood ratios of Trigger with nonpermuted data in block form:
 * 1. E->A v.s. E no relation with A
 * 2. A<-E->B with A--B v.s. E->A<-B
 * 3. Test E orthogonal with B | A | E->A and E->B. Two realizations:
 * 		1) Michoel:	E->A->B v.s. A<-E->B with A--B
 * 		2) Original: Define C=B-<B,A>, and test E no relation with C v.s. E->C (Not implemented)
 * Uses GSL BLAS.
 * Note: for each row, g must be the best eQTL of t of the same row.
 */
void pij_gassist_llr_block_buffed(const struct pij_gassist_llr_block_buffed_params* p);

/* Wrapper of pij_gassist_llr_block_buffed. Performs memory allocation and pre-calculations of
 * categorical mean and ratio, and the covariance matrix before invoking pij_gassist_llr_block_buffed
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix for A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix for B
 * nv:		Number of possible values for each genotype
 * llr1:	VECTORF (end-start). Log likelihood ratios for test 1.
 * llr2b:	MATRIXF (end-start,nt). Log likelihood ratios for test 2 bold.
 * llr2c:	MATRIXF (end-start,nt). Log likelihood ratios for test 2 conservative.
 * llr3:	MATRIXF (end-start,nt). Log likelihood ratios for test 3.
 * vb1:		VECTORF (ns). Constant buffer vector, must be set to 1 initially.
 * Return:	0 on success.
 * Notes:	1.	block range only applicable to A, all other transcripts as B are always considered.
 * 			2.	for each row, g must be the best eQTL of t of the same row.
 * 			3.	Definitions:	ng: number of SNPs=number of transcripts for A
 * 								nt: number of transcripts for B
 * 								ns: number of samples
 */
int pij_gassist_llr_block(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,VECTORF* llr1,MATRIXF* llr2b,MATRIXF* llr2c,MATRIXF* llr3,const VECTORF* vb1);

/* Calculates the log likelihood ratio of step 2 when the transcripts of those with eQTLs are different with those
 * to be tested against: A<-E->B with A--B v.s. E->A<-B. A and B can have different sets of transcripts here.
 * This is a more conservative version.
 * g:		MATRIXF (ng,ns) of genotype data
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns)	of transcript data for B
 * nv:		number of possible values of g.
 * llr2:	MATRIXF (ng,nt) of output
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For buffer purpose.
 * mmean:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For buffer purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For buffer purpose.
 * mmb1:	MATRIXF[nv] (ng,ns) Buffer matrix
 * vb1:		const VECTORF (ns). Must be set to 1 for all elements.
 */
void pij_gassist_llr2c_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,const VECTORF* vb1);

int pij_gassist_llr2c(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2);

/* Calculates the log likelihood ratio of step 2 when the transcripts of those with eQTLs are different with those
 * to be tested against: A<-E->B with A--B v.s. E->A. A and B can have different sets of transcripts here.
 * This is a less conservative version because it does not exclude B->A.
 * g:		MATRIXF (ng,ns) of genotype data
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns)	of transcript data for B
 * nv:		number of possible values of g.
 * llr2:	MATRIXF (ng,nt) of output
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For buffer purpose.
 * mmean:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For buffer purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For buffer purpose.
 * mmb1:	MATRIXF[nv] (ng,ns) Buffer matrix
 * mmb2:	MATRIXF[nv] (ng,nt) Buffer matrix
 * mb1:		MATRIXF (ng,nt). Buffer matrix
 * mb2:		MATRIXF (ng,nt). Buffer matrix
 * mb3:		MATRIXF (ng,nt). Buffer matrix
 * vb1:		const VECTORF (ns). Must be set to 1 for all elements.
 * vb2:		VECTORF (ng). Buffer vector
 */
void pij_gassist_llr2b_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,MATRIXF** mmb2,MATRIXF* mb1,MATRIXF* mb2,MATRIXF* mb3,const VECTORF* vb1,VECTORF* vb2);

/* Calculates the log likelihood ratio of step 2 when the transcripts of those with eQTLs are different with those
 * to be tested against: A<-E->B with A--B v.s. E->A<-B. A and B can have different sets of transcripts here.
 * g:		MATRIXF (ng,ns) of genotype data
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns)	of transcript data for B
 * nv:		number of possible values of g.
 * llr2:	MATRIXF (ng,nt) of output
 * Return:	0 on success
 */
int pij_gassist_llr2b(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2);

//#define pij_gassist_llr2 pij_gassist_llr2_bold
//#define pij_gassist_llr2_buffed pij_gassist_llr2_bold_buffed

/* Calculate log likelihood ratio of step 3 for single (E,A,B):
 * log likelihood(E->A->B)-log likelihood(A<-E->B with A--B)
 * fcov:	covariance (A,B) for each B
 * vratio:	(nv)	categorical ratio of each value of g
 * vmean1:	(nv)	categorical mean of A
 * vmean2:	(nv)	categorical mean of B
 * vb1:		(nv)	buffer vector f_alpha mu_{alpha 1}
 * Return:	log likelihood ratio
 */
FTYPE pij_gassist_llr3_one_from_stats_buffed(FTYPE fcov,const VECTORF* vratio,const VECTORF* vmean1,const VECTORF* vmean2,VECTORF* vb1);

/* Calculate log likelihood ratio of step 3 for single (E,A) but multiple B:
 * log likelihood(E->A->B)-log likelihood(A<-E->B with A--B)
 * g:		(ns) genotype data (E), each=0,1,...,nv-1
 * t1:		(ns) Supernormalizedtranscript data for A
 * t2:		(nt,ns) Supernormalized transcript data for Bs
 * nv:		Number of values each entry of g can take.
 * llr3:	(nt) output buffer for calculated llr
 * vb1:		(ns) constant buffer, must be set to all 1 initially
 * vratio:	(nv) categorical ratio of each value of g
 * vmean1:	(nv)	categorical mean of A
 * mmean:	(nt,nv)	categorical mean of each B for each value of g
 * vcov:	(nt) covariance (A,B) for each B
 * vb2:		(nv)	buffer vector f_alpha mu_{alpha 1}
 * vb3:		(nt)	buffer matrix for f_alpha mu_{alpha 2}
 * mb1:		(nv,ns)	buffer matrix of expanded representation of sample category
 */
void pij_gassist_llr3_E1_buffed(const VECTORG* g,const VECTORF* t1,const MATRIXF* t2,size_t nv,VECTORF* llr3,const VECTORF* vb1,VECTORF *vratio,VECTORF* vmean1,MATRIXF *mmean,VECTORF *vcov,VECTORF *vb2,VECTORF *vb3,MATRIXF *mb1);

/* Calculates log likelihood ratio for block of step 3 with buffer provided
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix for A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix for B
 * nv:		Number of possible values for each genotype
 * llr3:	MATRIXF (ng,nt). Log likelihood ratios for test 3.
 * vb1:		VECTORF (ns). Constant buffer vector, must be set to 1 initially.
 * mratio:	MATRIXF (nv,ng). Buffer matrix for categorical ratio
 * mmean1:	MATRIXF (nv,ng). Buffer matrix set for categorical mean
 * mmean2:	MATRIXF[nv] (ng,nt). Buffer matrix set for categorical mean
 * mcov:	MATRIXF (ng,nt). Buffer covariance matrix
 * mmb1:	MATRIXF[nv] (ng,nt). Buffer matrix
 * mmb2:	MATRIXF[nv] (ng,ns). Buffer matrix
 */
void pij_gassist_llr3_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr3,const VECTORF* vb1,MATRIXF *mratio,MATRIXF *mmean1,MATRIXF **mmean2,MATRIXF *mcov,MATRIXF **mmb1,MATRIXF **mmb2);

/* Multithread calculation of log likelihood ratios for 3 steps.
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix of A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix of B
 * llr1:	VECTORF (ng). Log likelihood ratios for test 1.
 * llr2b:	MATRIXF (ng,nt). Log likelihood ratios for test 2 bold.
 * llr2c:	MATRIXF (ng,nt). Log likelihood ratios for test 2 conservative.
 * llr3:	MATRIXF (ng,nt). Log likelihood ratios for test 3.
 * nv:		Number of possible values for each genotype
 * Return:	0 on success.
 */
int pij_gassist_llr(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2b,MATRIXF* llr2c,MATRIXF* llr3,size_t nv);

#ifdef __cplusplus
}
#endif
#endif
