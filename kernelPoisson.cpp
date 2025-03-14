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

__device__ int gpumerge(int a, int b, int k){

   return ((a == b)? k : 0);

}

__device__ void gpuInvertMatrix3x3(double *G, double *GI, double detG, int thread){

   double s = 1.0/detG;
   int start = 9*thread;

   GI[start]   =  s * (G[start+4]*G[start+8] - G[start+7]*G[start+5]);
   GI[start+1] = -s * (G[start+1]*G[start+8] - G[start+7]*G[start+2]);
   GI[start+2] =  s * (G[start+1]*G[start+5] - G[start+4]*G[start+2]);
   GI[start+3] = -s * (G[start+3]*G[start+8] - G[start+6]*G[start+5]);
   GI[start+4] =  s * (G[start]*G[start+8]   - G[start+6]*G[start+2]);
   GI[start+5] = -s * (G[start]*G[start+5]   - G[start+3]*G[start+2]);
   GI[start+6] =  s * (G[start+3]*G[start+7] - G[start+6]*G[start+4]);
   GI[start+7] = -s * (G[start]*G[start+7]   - G[start+6]*G[start+1]);
   GI[start+8] =  s * (G[start]*G[start+4]   - G[start+3]*G[start+1]);

}

__device__ double gpumyElementMetric(int nDOFs, double *xx, double *yy, double *zz,
	                  	     int dim, double *dLBasisdx, double *LtoGMap, 
				     int elem, int sizexx1, int sizexx2, int thread, 
				     int l, int max_nd, double *G, double *GI, 
				     double *dx, double *Metric) {

   int start = 9*thread;
   int n = (sizexx2 < nDOFs) ? sizexx2 : nDOFs;
   
   for (int i=0; i<dim; i++) {
      double sum =0.0;
      for (int j=0; j<n; j++) 
         sum += xx[sizexx1*j+elem]*dLBasisdx[l*dim*max_nd+i*dim+j];
      dx[start+i*dim] = sum;
   }
   for (int i=0; i<dim; i++) {
      double sum =0.0;
      for (int j=0; j<n; j++) 
         sum += yy[sizexx1*j+elem]*dLBasisdx[l*dim*max_nd+i*dim+j];
      dx[start+i*dim+1] = sum;
   }
   for (int i=0; i<dim; i++) {
      double sum =0.0;
      for (int j=0; j<n; j++) 
	 sum += zz[sizexx1*j+elem]*dLBasisdx[l*dim*max_nd+i*dim+j];
      dx[start+i*dim+2] = sum;
   }
   for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
         double sum = 0.0;
         for (int k=0; k<dim; k++) 
            sum += dx[start+i*dim+k]*dx[start+j*dim+k];
         G[start+j*dim+i] = sum;
      }
   }
   double DetG = G[start]   * ( G[start+4]*G[start+8] - G[start+5]*G[start+7] ) + 
	         G[start+1] * ( G[start+5]*G[start+6] - G[start+3]*G[start+8] ) + 
                 G[start+2] * ( G[start+3]*G[start+7] - G[start+4]*G[start+6] );
   gpuInvertMatrix3x3(G, GI, DetG, thread);
   for (int i=0; i<9; i++) 
      Metric[start*i] = GI[start+i];
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

__device__ void gpuModuleLocalMatrixVecSO(int n, int nd, int nb, double *x, double *y, 
		                          double *z, int dim, double *refbasis, 
					  double *refdBasisdx, double *ips, int ngp, 
					  int elem, int *l2g, double *values, int *cols, 
					  int *rows, double *rhs, int sizexx1, int sizexx2, 
					  double *dBasisdx, double *DetJ, double *STIFF, 
					  double *FORCE, double *DiffCoeff, 
					  double *SourceCoeff, double *LtoGMap, int thread, 
					  int max_nd, double *G, double *GI, double *dx,
					  double *Metric) {

   for (int i=0; i<dim*ngp*nd; i++) 
      dBasisdx[thread*dim*ngp*max_nd+i] = 0.0;
   for (int l=0; l<ngp; l++) {
      DetJ[thread*ngp+l] = gpumyElementMetric(nd, x, y, z, dim, refdBasisdx, LtoGMap, elem, 
		                              sizexx1, sizexx2, thread, l, max_nd, G, GI, 
					      dx, Metric);
      for (int m=1; m<nd; m++)
         for (int k=1; k<dim; k++)
            for (int j=1; j<dim; j++)
	       dBasisdx[thread*dim*ngp*nd+nd*ngp*j+ngp*m+l] += refdBasisdx[dim*nd*l+nd*k+m]*
		                                               LtoGMap[thread*dim*dim+dim*k+j];
   }
   for (int i=0; i<nd*nd; i++) 
      STIFF[thread*nd*nd+i]=0.0;
   for (int i=0; i<ngp; i++) 
      DiffCoeff[thread*ngp+i]=1.0;
   for (int j=0; j<nd; j++)
      for (int i=0; i<nd; i++)
         for (int k=0; k<dim; k++)
            for (int l=0; l<ngp; l++)
	       STIFF[thread*nd*nd+j*nd+i] += dBasisdx[thread*ngp*nd*dim+nd*ngp*k+ngp*i+l]*
		                             dBasisdx[thread*ngp*nd*dim+nd*ngp*k+ngp*j+l]*
					     DiffCoeff[thread*ngp+l]*DetJ[thread*ngp+l]*
					     ips[l];
   for (int i=0; i<ngp; i++) 
      SourceCoeff[thread*ngp+i] = 1.0;
   for (int i=0; i<nd; i++) 
      FORCE[thread*nd+i] = 0.0;
   for (int i=0; i<nd; i++) 
      for (int l=0; l<ngp; l++)
         FORCE[thread*nd+i] += refbasis[i*ngp+l]*SourceCoeff[thread*ngp+l]*DetJ[thread*ngp+l]*
		               ips[l];
   for (int i=0; i<nd; i++) {
      for (int j=0; j<nd; j++) {
        int colind = 0;
	for (int k=rows[l2g[i*sizexx1+elem]]-1; k<rows[l2g[i*sizexx1+elem]]-2; k++)
	  colind += gpumerge(0, cols[k], k);
	values[colind] += STIFF[thread*nd*nd+j*nd+i];
      }
      rhs[i*nd+elem] += FORCE[thread*n+i];
   }
   
}
__global__ void kernel_ModuleLocalMatrixVecSO(int active, int *elemdofs, double *x, double *y,
		                              double *z, int dim, double *refbasis, 
					      double *refdBasisdx, double *ips, int ngp, int *l2g, 
					      double *values, int *cols, int *rows, double *rhs, 
					      int sizexx1, int sizexx2, double *dBasisdx, 
					      double *DetJ, double *STIFF, double *FORCE,
                                              double *DiffCoeff, double *SourceCoeff, 
					      double *LtoGMap, int max_nd, double *G, double *GI, 
					      double *dx, double *Metric) {
   
   int thread = threadIdx.x;
   int elem = thread+blockIdx.x*blockDim.x;
   int n = elemdofs[elem];
   int nb = elemdofs[2*active+elem];
   int nd = nb + elemdofs[active+elem];

   if (elem < active) 
      gpuModuleLocalMatrixVecSO(n, nd, nb, x, y, z, dim, refbasis, refdBasisdx, ips, ngp, elem, l2g, 
		                values, cols, rows, rhs, sizexx1, sizexx2, dBasisdx, DetJ, STIFF, 
				FORCE, DiffCoeff, SourceCoeff, LtoGMap, thread, max_nd, G, GI, dx, 
				Metric);

}
extern "C" 
   {
   void launchModuleLocalMatrixVecSO(int active, int *elemdofs, double *x, double *y, double *z, int dim, 
		                     double *refbasis, double *refdBasisdx, double *ips, int ngp, int *l2g,
				     double *values, int *cols, int *rows, double *rhs, int sizexx1, 
				     int sizexx2, double *dBasisdx, double *DetJ, double *STIFF, 
				     double *FORCE, double *DiffCoeff, double *SourceCoeff, double *LtoGMap,
				     int max_nd, double *G, double *GI, double *dx, double *Metric, 
				     int thread_per_block) {
      
      dim3 nblocks=dim3((active+thread_per_block)/thread_per_block,1,1);
      dim3 nthreads=dim3(thread_per_block,1,1);

      kernel_ModuleLocalMatrixVecSO<<<nblocks,nthreads,0,0>>> (active, elemdofs, x, y, z, dim, refbasis, 
		                                               refdBasisdx, ips, ngp, l2g, values, cols, rows, 
							       rhs, sizexx1, sizexx2, dBasisdx, DetJ, STIFF, 
							       FORCE, DiffCoeff, SourceCoeff, LtoGMap, max_nd, 
							       G, GI, dx, Metric);
   }
}
