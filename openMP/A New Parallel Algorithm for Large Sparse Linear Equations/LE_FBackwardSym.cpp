#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "LE_SymSprsMatFunc.h"
#include "LE_SymSprsMatDef.h"

#define NUM_THREADS 34
#define UPBOUND 4096

extern int dispatch[NUM_THREADS+1][UPBOUND];
extern int predispatch[NUM_THREADS+1][UPBOUND];
extern int mappingf[UPBOUND];
extern int mappingb[UPBOUND];
extern double x_[UPBOUND];
extern int *rs_u_;
extern int *i_u;
extern double *u_u_;
//////////////////////////////////////////////////////////////////////
// 函 数 名:          //LE_FBackwardSym
// 描    述:          //对称矩阵前推回代方法求解方程组
//
// 输入参数:          // U阵结构及U阵值，右端项b
// 输出参数:          // 右端项x，维数为pU的维数（解向量）
// 返 回 值:          // 无
// 其    他:          // 不会影响U阵中矩阵的任何信息,created by xdc 2014/6/16
//////////////////////////////////////////////////////////////////////
void LE_FBackwardSym(SprsUMatRealStru *pFU,double b[],double x[])
{
	int i,j,k,l,m;
	int ks,ke;
	int iDim;
	int *rs_u,*j_u;
	double *d_u,*u_u;
	double xc;

	d_u=pFU->d_u;
	u_u=pFU->u_u;
	rs_u=pFU->uMax.rs_u;
	j_u=pFU->uMax.j_u;
	iDim=pFU->uMax.iDim;

	//////////////////////////////////////////////////
	//						//
	//	paralleled by wangdong 			//
	//	wangd1991@foxmail.com			//
	//						//
	//////////////////////////////////////////////////

#pragma omp parallel for private(m,i,ks,ke,k,j,l,xc) num_threads(NUM_THREADS)
	for(l=1;l<=NUM_THREADS;l++)
	{
		for(m=0; m<dispatch[l][0]; m++)
		{
			i=dispatch[l][m+1];
			ks = rs_u_[i];
			ke = rs_u_[i+1];
			x_[i] = b[mappingf[i]];
			for(k = ks; k < ke; k ++)
			{
				j = i_u[k];
				x_[i] -= u_u_[k]*x_[j];
			}
		}
	}
#pragma omp parallel for private(m,i,ks,ke,k,j,l,xc) num_threads(NUM_THREADS)
	for(l=1;l<=NUM_THREADS;l++)
	{
		for(m=0; m<predispatch[l][0]; m++)
		{
			i=predispatch[l][m+1];
			ks = rs_u_[i];
			ke = rs_u_[i+1];
			xc=b[mappingf[i]];
			for(k = ks; k < ke; k ++)
			{
				j = i_u[k];
				xc -= u_u_[k]*x_[j];
			}
			x_[i]=xc;
		}
	}

	for(m=0; m<dispatch[0][0]; m++)
	{
		i=dispatch[0][m+1];
		ks = rs_u_[i];
		ke = rs_u_[i+1];
		xc=b[mappingf[i]];
		for(k = ks; k < ke; k ++)
		{
			j = i_u[k];
			xc -= u_u_[k]*x_[j];
		}
		x_[i]=xc;
	}


	for(m=dispatch[0][0]-1;m>=0;m--)
	{
		i=dispatch[0][m+1];
		ks = rs_u[i];
		ke = rs_u[i+1];
		xc=x_[i]*d_u[i];
		for(k=ke-1; k>=ks; k--)
		{
			j = j_u[k];
			xc -= u_u[k]*x_[j];
		}
		x_[i]=xc;
		x[mappingf[i]] = xc;
	}
#pragma omp parallel for private(m,i,ks,ke,k,j,l,xc) num_threads(NUM_THREADS)
	for(l=1;l<=NUM_THREADS;l++)
	{
		for(m=predispatch[l][0]-1;m>=0;m--)
		{
			i=predispatch[l][m+1];
			ks = rs_u[i];
			ke = rs_u[i+1];
			xc=x_[i]*d_u[i];
			for(k=ke-1; k>=ks; k--)
			{
				j = j_u[k];
				xc -= u_u[k]*x_[j];
			}
			x_[i]=xc;
			x[mappingf[i]] = xc;
		}
	}

#pragma omp parallel for private(m,i,ks,ke,k,j,l,xc) num_threads(NUM_THREADS)
	for(l=1;l<=NUM_THREADS;l++)
	{
		for(m=dispatch[l][0]-1;m>=0;m--)
		{
			i=dispatch[l][m+1];
			ks = rs_u[i];
			ke = rs_u[i+1];
			xc=x_[i]*d_u[i];
			for(k=ke-1; k>=ks; k--)
			{
				j = j_u[k];
				xc -= u_u[k]*x_[j];
			}
			x_[i]=xc;
			x[mappingf[i]] = xc;
		}
	}
}
