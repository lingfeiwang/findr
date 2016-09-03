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
//void pij_gassist_llr2c_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,const VECTORF* vb1);

//int pij_gassist_llr2c(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2);

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
//void pij_gassist_llr2b_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,MATRIXF** mmb2,MATRIXF* mb1,MATRIXF* mb2,MATRIXF* mb3,const VECTORF* vb1,VECTORF* vb2);

/* Calculates the log likelihood ratio of step 2 when the transcripts of those with eQTLs are different with those
 * to be tested against: A<-E->B with A--B v.s. E->A<-B. A and B can have different sets of transcripts here.
 * g:		MATRIXF (ng,ns) of genotype data
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns)	of transcript data for B
 * nv:		number of possible values of g.
 * llr2:	MATRIXF (ng,nt) of output
 * Return:	0 on success
 */
//int pij_gassist_llr2b(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2);

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
//FTYPE pij_gassist_llr3_one_from_stats_buffed(FTYPE fcov,const VECTORF* vratio,const VECTORF* vmean1,const VECTORF* vmean2,VECTORF* vb1);

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
//void pij_gassist_llr3_E1_buffed(const VECTORG* g,const VECTORF* t1,const MATRIXF* t2,size_t nv,VECTORF* llr3,const VECTORF* vb1,VECTORF *vratio,VECTORF* vmean1,MATRIXF *mmean,VECTORF *vcov,VECTORF *vb2,VECTORF *vb3,MATRIXF *mb1);

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
//void pij_gassist_llr3_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr3,const VECTORF* vb1,MATRIXF *mratio,MATRIXF *mmean1,MATRIXF **mmean2,MATRIXF *mcov,MATRIXF **mmb1,MATRIXF **mmb2);

/* Multithread calculation of log likelihood ratios for 5 tests.
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix of A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix of B
 * llr1:	VECTORF (ng). Log likelihood ratios for test 1. Tests E->A v.s. E  A.
 * llr2:	MATRIXF (ng,nt). Log likelihood ratios for test 2. Tests E->B v.s. E  B.
 * llr3:	MATRIXF (ng,nt). Log likelihood ratios for test 3. Tests E->A->B v.s. E->A->B with E->B.
 * llr4:	MATRIXF (ng,nt). Log likelihood ratios for test 4. Tests E->A->B with E->B v.s. E->A  B.
 * llr5:	MATRIXF (ng,nt). Log likelihood ratios for test 5. Tests E->A->B with E->B v.s. A<-E->B.
 * nv:		Number of possible values for each genotype
 * Return:	0 on success.
 */
int pij_gassist_llr(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5,size_t nv);



























#ifdef __cplusplus
}
#endif
#endif
