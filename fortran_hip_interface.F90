module fortran_hip
   use hipfort
   use hipfort_check
  
   implicit none

   interface
      subroutine launchModuleLocalMatrixVecSO(active, elemdofs, x, y, z, dim, refbasis,   &
                                              refdBasisdx, ips, ngp, l2g, values, cols,   &
                                              rows, rhs, sizexx1, sizexx2,  dBasisdx,     &
                                              DetJ, STIFF, FORCE, DiffCoeff, SourceCoeff, &
                                              th_per_block)                               &
                                              bind(c, name="launchModuleLocalMatrixVecSO")


         use hipfort_types
         use iso_c_binding
         
         implicit none

         type(c_ptr),    value :: elemdofs, x, y, z, refbasis, refdBasisdx, ips, l2g, values
         type(c_ptr),    value :: cols, rows, rhs, dBasisdx, DetJ, STIFF, FORCE, DiffCoeff
         type(c_ptr),    value :: SourceCoeff
         integer(c_int), value :: active, dim, ngp, sizexx1, sizexx2, th_per_block

      end subroutine
   end interface

end module 
