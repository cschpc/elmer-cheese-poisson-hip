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

   return (test)? a : b;
}

__device__ void gpuInvertMatrix3x3(double * __restrict__ G, double * __restrict__ GI, double detG){

   GI[0] =  (G[4]*G[8] - G[5]*G[7])/detG;
   GI[1] = -(G[1]*G[8] - G[2]*G[7])/detG;
   GI[2] =  (G[1]*G[5] - G[2]*G[4])/detG;
   GI[3] = -(G[3]*G[8] - G[5]*G[6])/detG;
   GI[4] =  (G[0]*G[8] - G[2]*G[6])/detG;
   GI[5] = -(G[0]*G[5] - G[2]*G[3])/detG;
   GI[6] =  (G[3]*G[7] - G[4]*G[6])/detG;
   GI[7] = -(G[0]*G[7] - G[1]*G[6])/detG;
   GI[8] =  (G[0]*G[4] - G[1]*G[3])/detG;

}

__device__ double gpumyElementMetric(int nd, double * __restrict__ xx, double * __restrict__ yy, 
		                     double * __restrict__ zz, int dim, double * __restrict__ dLBasisdx, 
				     double * __restrict__ LtoGMap, int elem, int sizexx1, int sizexx2, 
				     int l) {

   double G[9], GI[9], dx[9];
   int n = (sizexx2 < nd) ? sizexx2 : nd;
   
   for (int i=0; i<dim; i++) {
      double sumx =0.0, sumy =0.0, sumz =0.0;
      for (int j=0; j<n; j++) {
         sumx += xx[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
         sumy += yy[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
         sumz += zz[sizexx1*j+elem]*dLBasisdx[l*dim*nd+i*nd+j];
      }
      dx[i*dim]   = sumx;
      dx[i*dim+1] = sumy;
      dx[i*dim+2] = sumz;
   }

   for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
         double sumx = 0.0;
         for (int k=0; k<dim; k++) 
            sumx += dx[i*dim+k]*dx[j*dim+k];
         G[j*dim+i] = sumx;
      }
   }

   double DetG = G[0] * ( G[4]*G[8] - G[7]*G[5] ) + G[3] * ( G[7]*G[2] - G[1]*G[8] ) + 
                 G[6] * ( G[1]*G[5] - G[4]*G[2] );

   gpuInvertMatrix3x3(G, GI, DetG);
   DetG = sqrt(DetG);

   for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
         double sumx = 0.0;
         for (int k=0; k<dim; k++) 
            sumx += dx[k*dim+i]*GI[j*dim+k];
         LtoGMap[j*dim+i] = sumx;
      }
   }

   return DetG;

}

__device__ void gpuModuleLocalMatrixVecSO(int n, int nd, int nb, double * __restrict__ x, double * __restrict__ y, 
		                          double * __restrict__ z, int dim, double * __restrict__ refbasis, 
					  double * __restrict__ refdBasisdx, double * __restrict__ ips, int ngp, 
					  int elem, int * __restrict__ l2g, double * __restrict__ values, 
					  int * __restrict__ cols, int * __restrict__ rows, double * __restrict__ rhs, 
					  int sizexx1, int sizexx2, double * __restrict__ dBasisdx, 
					  double * __restrict__ DetJ, double * __restrict__ STIFF, 
					  double * __restrict__ FORCE, double * __restrict__ DiffCoeff, 
					  double * __restrict__ SourceCoeff) {

   double LtoGMap[9];
   int size = dim*ngp*nd;
   for (int i=0; i<size; i++) 
      dBasisdx[elem*size+i] = 0.0;

   for (int l=0; l<ngp; l++) {
      DetJ[elem*ngp+l] = gpumyElementMetric(nd, x, y, z, dim, refdBasisdx, LtoGMap, elem, sizexx1, sizexx2, l);
      for (int m=0; m<nd; m++)
         for (int k=0; k<dim; k++)
            for (int j=0; j<dim; j++)
	       dBasisdx[elem*size+nd*ngp*j+ngp*m+l] += refdBasisdx[dim*nd*l+nd*k+m]* LtoGMap[dim*k+j];
   }

   for (int i=0; i<ngp; i++) 
      DiffCoeff[elem*ngp+i] = 1.0;
   for (int j=0; j<nd; j++)
      for (int i=0; i<nd; i++) {
	 STIFF[elem*nd*nd+j*nd+i] = 0.0;
         for (int k=0; k<dim; k++) {
            double sum = 0.0;
            for (int l=0; l<ngp; l++)
	       sum += dBasisdx[elem*size+nd*ngp*k+ngp*i+l]* dBasisdx[elem*size+nd*ngp*k+ngp*j+l]* DiffCoeff[elem*ngp+l]*
		      DetJ[elem*ngp+l]*ips[l];
	    STIFF[elem*nd*nd+j*nd+i] += sum;
	 }   
      }	 

   for (int i=0; i<ngp; i++) 
      SourceCoeff[elem*ngp+i] = 1.0;
   for (int i=0; i<nd; i++) {
      double sum = 0.;
      for (int l=0; l<ngp; l++)
         sum += refbasis[i*nd+l]*SourceCoeff[elem*ngp+l]* DetJ[elem*ngp+l]*ips[l];
      FORCE[elem*nd+i] = sum;
   }

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
__global__ void __launch_bounds__(512, 1) kernel_ModuleLocalMatrixVecSO(int active, int * __restrict__ elemdofs, double * __restrict__ x, 
		                                                        double * __restrict__ y, double * __restrict__ z, 
									int dim, double * __restrict__ refbasis, double * __restrict__ refdBasisdx, 
									double * __restrict__ ips, int ngp, int * __restrict__ l2g, 
					                                double * __restrict__ values, int * __restrict__ cols, int * __restrict__ rows, 
									double * __restrict__ rhs, int sizexx1, int sizexx2, 
									double * __restrict__ dBasisdx, double * __restrict__ DetJ, 
									double * __restrict__ STIFF, double * __restrict__ FORCE, 
									double * __restrict__ DiffCoeff, double * __restrict__ SourceCoeff) {
   
   int elem = blockIdx.x*blockDim.x + threadIdx.x;
   if (elem < active) {
      int n = elemdofs[elem];
      int nb = elemdofs[2*active+elem];
      int nd = nb + elemdofs[active+elem];
      gpuModuleLocalMatrixVecSO(n, nd, nb, x, y, z, dim, refbasis, refdBasisdx, ips, ngp, elem, l2g, values, cols, rows, rhs, sizexx1, sizexx2, 
		                dBasisdx, DetJ, STIFF, FORCE, DiffCoeff, SourceCoeff);
   }

}
extern "C" {
   void launchModuleLocalMatrixVecSO(int active, int * __restrict__ elemdofs, double * __restrict__ x, double * __restrict__ y, 
		                     double * __restrict__ z, int dim, double * __restrict__ refbasis, double * __restrict__ refdBasisdx, 
				     double * __restrict__ ips, int ngp, int * __restrict__ l2g, double * __restrict__ values, 
				     int * __restrict__ cols, int * __restrict__ rows, double * __restrict__ rhs, int sizexx1, int sizexx2, 
				     double * __restrict__ dBasisdx, double * __restrict__ DetJ, double * __restrict__ STIFF, 
				     double * __restrict__ FORCE, double * __restrict__ DiffCoeff, double * __restrict__ SourceCoeff, 
				     int thread_per_block) {
      
      dim3 nblocks=dim3((active+thread_per_block)/thread_per_block,1,1);
      dim3 nthreads=dim3(thread_per_block,1,1);

      kernel_ModuleLocalMatrixVecSO<<<nblocks,nthreads,0,0>>> (active, elemdofs, x, y, z, dim, refbasis, refdBasisdx, ips, ngp, l2g, 
		                                               values, cols, rows, rhs, sizexx1, sizexx2, dBasisdx, DetJ, STIFF, 
							       FORCE, DiffCoeff, SourceCoeff);
   }
}
