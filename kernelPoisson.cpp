#include "hip/hip_runtime.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <assert.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(hipError_t code, const char *file, int line, bool abort=true) {

   if (code != hipSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }

}

__device__ int gpumerge(int a, int b, bool test){

   return (test)? a: b;
}

__device__ void gpuInvertMatrix3x3(__restrict__ double *G, __restrict__ double *GI, double detG, int elem){

   double s = 1.0/detG;
   int start = 9*elem;

     GI[start] =  s * (G[start+4]*G[start+8] - G[start+5]*G[start+7]);
   GI[start+1] = -s * (G[start+1]*G[start+8] - G[start+2]*G[start+7]);
   GI[start+2] =  s * (G[start+1]*G[start+5] - G[start+2]*G[start+4]);

   GI[start+3] = -s * (G[start+3]*G[start+8] - G[start+5]*G[start+6]);
   GI[start+4] =  s * (G[start]  *G[start+8] - G[start+2]*G[start+6]);
   GI[start+5] = -s * (G[start]  *G[start+5] - G[start+2]*G[start+3]);

   GI[start+6] =  s * (G[start+3]*G[start+7] - G[start+4]*G[start+6]);
   GI[start+7] = -s * (G[start]  *G[start+7] - G[start+1]*G[start+6]);
   GI[start+8] =  s * (G[start]  *G[start+4] - G[start+1]*G[start+3]);

}

__device__ double gpumyElementMetric(int nd, __restrict__ double *xx, __restrict__ double *yy, __restrict__ double *zz,
	                  	     int dim, __restrict__ double *dLBasisdx, __restrict__ double *LtoGMap, 
				     int elem, int sizexx1, int sizexx2, int l, 
				     __restrict__ double *G, __restrict__ double *GI, __restrict__ double *dx, 
				     __restrict__ double *Metric) {

   int start = 9*elem;
   int n = (sizexx2 < nd) ? sizexx2 : nd;
   
   for (int i=0; i<dim; i++) {
      double sumx =0.0, sumy =0.0, sumz =0.0;
      for (int j=0; j<n; j++) {
         sumx += xx[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
         sumy += yy[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
         sumz += zz[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
      }
      dx[start+i*dim]   = sumx;
      dx[start+i*dim+1] = sumy;
      dx[start+i*dim+2] = sumz;
   }

   for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
         double sum = 0.0;
         for (int k=0; k<dim; k++) 
            sum += dx[start+i*dim+k]*dx[start+j*dim+k];
         G[start+j*dim+i] = sum;
      }
   }

   double DetG =   G[start] * ( G[start+4]*G[start+8] - G[start+7]*G[start+5] ) + 
	         G[start+3] * ( G[start+7]*G[start+2] - G[start+1]*G[start+8] ) + 
                 G[start+6] * ( G[start+1]*G[start+5] - G[start+4]*G[start+2] );

   gpuInvertMatrix3x3(G, GI, DetG, elem);
   for (int i=0; i<9; i++) 
      Metric[start+i] = GI[start+i];
   DetG = sqrt(DetG);

   for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
         double sum = 0.0;
         for (int k=0; k<dim; k++) 
            sum += dx[start+k*dim+i]*Metric[start+j*dim+k];
         LtoGMap[start+j*dim+i] = sum;
      }
   }

   return DetG;

}

__device__ void gpuModuleLocalMatrixVecSO(int n, int nd, int nb, __restrict__ double *x, __restrict__ double *y, 
		                          __restrict__ double *z, int dim, __restrict__ double *refbasis, 
					  __restrict__ double *refdBasisdx, __restrict__ double *ips, int ngp, 
					  int elem, int *l2g, __restrict__ double *values, int *cols, 
					  int *rows, __restrict__ double *rhs, int sizexx1, int sizexx2, 
					  __restrict__ double *dBasisdx, __restrict__ double *DetJ, __restrict__ double *STIFF, 
					  __restrict__ double *FORCE, __restrict__ double *DiffCoeff, 
					  __restrict__ double *SourceCoeff, __restrict__ double *LtoGMap, int max_active, 
					  __restrict__ double *G, __restrict__ double *GI, __restrict__ double *dx, __restrict__ double *Metric) {

   int size = dim*ngp*nd;
   for (int i=0; i<size; i++) 
      dBasisdx[elem*size+i] = 0.0;

   for (int l=0; l<ngp; l++) {
      DetJ[elem*ngp+l] = gpumyElementMetric(nd, x, y, z, dim, refdBasisdx, LtoGMap, elem, 
		                                  sizexx1, sizexx2, l, G, GI, dx, Metric);
      for (int m=0; m<nd; m++)
         for (int k=0; k<dim; k++)
            for (int j=0; j<dim; j++)
	       dBasisdx[elem*size+nd*ngp*j+ngp*m+l] += refdBasisdx[dim*nd*l+nd*k+m]*
		                                       LtoGMap[elem*dim*dim+dim*k+j];
   }

   for (int i=0; i<nd*nd; i++) 
      STIFF[elem*nd*nd+i] = 0.0;
   for (int i=0; i<ngp; i++) 
      DiffCoeff[elem*ngp+i] = 1.0;
   for (int j=0; j<nd; j++)
      for (int i=0; i<nd; i++)
         for (int k=0; k<dim; k++)
            for (int l=0; l<ngp; l++)
	       STIFF[elem*nd*nd+j*nd+i] += dBasisdx[elem*size+nd*ngp*k+ngp*i+l]*
		                           dBasisdx[elem*size+nd*ngp*k+ngp*j+l]*
					   DiffCoeff[elem*ngp+l]*DetJ[elem*ngp+l]*ips[l];

   for (int i=0; i<ngp; i++) 
      SourceCoeff[elem*ngp+i] = 1.0;
   for (int i=0; i<nd; i++) 
      FORCE[elem*nd+i] = 0.0;
   for (int i=0; i<nd; i++) 
      for (int l=0; l<ngp; l++)
         FORCE[elem*nd+i] += refbasis[i*nd+l]*SourceCoeff[elem*ngp+l]*
		             DetJ[elem*ngp+l]*ips[l];

   for (int i=0; i<nd; i++) {
      for (int j=0; j<nd; j++) {
        int colind = 0;
	for (int k=rows[l2g[i*sizexx1+elem]-1]; k<rows[l2g[i*sizexx1+elem]]; k++)
	  colind += gpumerge(k, 0, (l2g[j*sizexx1+elem] == cols[k-1]));
	values[colind-1] += STIFF[elem*nd*nd+j*nd+i];
      }
      rhs[l2g[i*sizexx1+elem]-1] += FORCE[elem*nd+i];
   }
}
__global__ void kernel_ModuleLocalMatrixVecSO(int active, int *elemdofs, __restrict__ double *x, __restrict__ double *y,
		                              __restrict__ double *z, int dim, __restrict__ double *refbasis, 
					      __restrict__ double *refdBasisdx, __restrict__ double *ips, int ngp, int *l2g, 
					      __restrict__ double *values, int *cols, int *rows, __restrict__ double *rhs, 
					      int sizexx1, int sizexx2, __restrict__ double *dBasisdx, 
					      __restrict__ double *DetJ, __restrict__ double *STIFF, __restrict__ double *FORCE,
                                              __restrict__ double *DiffCoeff, __restrict__ double *SourceCoeff, 
					      __restrict__ double *LtoGMap, __restrict__ double *G, __restrict__ double *GI, __restrict__ double *dx, 
					      __restrict__ double *Metric, int max_active) {
   
   int elem = blockIdx.x*blockDim.x + threadIdx.x;
   if (elem < active) {
      int n = elemdofs[elem];
      int nb = elemdofs[2*active+elem];
      int nd = nb + elemdofs[active+elem];
      gpuModuleLocalMatrixVecSO(n, nd, nb, x, y, z, dim, refbasis, refdBasisdx, ips, ngp, elem, l2g, 
		                values, cols, rows, rhs, sizexx1, sizexx2, dBasisdx, DetJ, STIFF, 
				FORCE, DiffCoeff, SourceCoeff, LtoGMap, max_active, G, GI, dx, Metric);
   }

}
extern "C" {
   void launchModuleLocalMatrixVecSO(int active, int *elemdofs, double *x, double *y, double *z, int dim, 
		                     double *refbasis, double *refdBasisdx, double *ips, int ngp, int *l2g,
				     double *values, int *cols, int *rows, double *rhs, int sizexx1, 
				     int sizexx2, double *dBasisdx, double *DetJ, double *STIFF, 
				     double *FORCE, double *DiffCoeff, double *SourceCoeff, double *LtoGMap,
				     double *G, double *GI, double *dx, double *Metric, int thread_per_block, 
				     int max_active) {
      
      dim3 nblocks=dim3((active+thread_per_block)/thread_per_block,1,1);
      dim3 nthreads=dim3(thread_per_block,1,1);

      kernel_ModuleLocalMatrixVecSO<<<nblocks,nthreads,0,0>>> (active, elemdofs, x, y, z, dim, refbasis, 
		                                               refdBasisdx, ips, ngp, l2g, values, cols, rows, 
							       rhs, sizexx1, sizexx2, dBasisdx, DetJ, STIFF, 
							       FORCE, DiffCoeff, SourceCoeff, LtoGMap, G, GI, 
							       dx, Metric, max_active);
      hipStreamSynchronize(0);
   }
}
