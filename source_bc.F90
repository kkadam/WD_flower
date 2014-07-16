!************************************************************************
!*  
!*  SOURCE_BC
!*
!************************************************************************
subroutine source_bc
implicit none
include 'runhydro.h'
!************************************************************************
!*
!
!*
!************************************************************************
!*
!*  Subroutine Arguments


!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: etot
common /energy/ etot

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

logical :: iam_on_top, iam_on_bottom, iam_on_axis,           &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                &
             in_neighbor, out_neighbor, root,                &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                  &
           in_neighbor, out_neighbor, root,                  &
           REAL_SIZE, INT_SIZE, numprocs
#endif 
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,           &
                        iam_on_bottom, iam_on_axis,          &
                        iam_on_edge, down_neighbor,          &
                        up_neighbor, in_neighbor,            &
                        out_neighbor, root, column_num,      &
                        row_num, pe_grid, iam_root,          &
                        REAL_SIZE, INT_SIZE


!*
!************************************************************************
!*
!*   Local Variables

real :: sum

real :: emin

integer :: J, K, L

!*
!************************************************************************
!  initialize the local variables
emin = taumin**gamma

do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if ( etot(J,K,L) < emin ) then
            etot(J,K,L) = emin
         endif
      enddo
   enddo
enddo

if( iam_on_axis ) then

   do L = philwb, phiupb
      do K = zlwb, zupb
         s(rlwb,K,L) = 0.0
      enddo
   enddo

   do K = zlwb, zupb

      sum = 0.0
      do L = philwb, phiupb
         sum = sum + etot(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         etot(rlwb,K,L) = sum
      enddo

      sum = 0.0
      do L = philwb, phiupb
         sum = sum + a(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         a(rlwb,K,L) = sum
      enddo

      sum = 0.0
      do L = philwb, phiupb
         sum= sum + t(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         t(rlwb,K,L) = sum 
      enddo
   enddo

   if( isym == 3 ) then
      do L = philwb, phiupb
         do K = zlwb, zupb
            etot(rlwb-1,K,L)  = etot(rlwb,K,L)
            s(rlwb-1,K,L) = - s(rlwb+1,K,L)
            a(rlwb-1,K,L) =   a(rlwb,K,L)
            t(rlwb-1,K,L) =   t(rlwb,K,L)
         enddo
      enddo
   else
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            etot(rlwb-1,K,L)               = etot(rlwb,K,L+numphi_by_two)
            etot(rlwb-1,K,L+numphi_by_two) = etot(rlwb,K,L)
         enddo
      enddo
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            s(rlwb-1,K,L)               = - s(rlwb+1,K,L+numphi_by_two)
            s(rlwb-1,K,L+numphi_by_two) = - s(rlwb+1,K,L)
         enddo
      enddo
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            t(rlwb-1,K,L)               = t(rlwb,K,L+numphi_by_two)
            t(rlwb-1,K,L+numphi_by_two) = t(rlwb,K,L)
         enddo
      enddo
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            a(rlwb-1,K,L)               = a(rlwb,K,L+numphi_by_two)
            a(rlwb-1,K,L+numphi_by_two) = a(rlwb,K,L)
         enddo
      enddo
   endif
endif

if( iam_on_bottom ) then
   if( isym == 1 ) then
      if( boundary_condition(1) == 1 ) then
         ! wall boundary condition
         do L = philwb, phiupb
            do J = rlwb, rupb
               etot(J,zlwb-1,L)  = etot(J,zlwb,L)
               s(J,zlwb-1,L) =  s(J,zlwb,L)
               t(J,zlwb-1,L) = -t(J,zlwb+1,L)
               t(J,zlwb,L)   = 0.0
               a(J,zlwb-1,L)  =  a(J,zlwb,L)
            enddo
         enddo
      else if( boundary_condition(1) == 2 ) then
         ! free boundary condition
         do L = philwb, phiupb
            do J = rlwb, rupb
               etot(J,zlwb-1,L)  = etot(J,zlwb,L)
               s(J,zlwb-1,L) = s(J,zlwb,L)
               a(J,zlwb-1,L) = a(J,zlwb,L)
               t(J,zlwb-1,L) = t(J,zlwb,L)
            enddo
         enddo
      else  
         ! dirichlet bc case
         do L = philwb, phiupb
            do J = rlwb, rupb
               s(J,zlwb-1,L) = 0.0
               a(J,zlwb-1,L) = 0.0
               t(J,zlwb-1,L) = 0.0
               t(J,zlwb,L) = 0.5*(abs(t(J,zlwb,L)) -    &
                                      t(J,zlwb,L) )
           enddo
         enddo
      endif
   else
      ! if isym /= 1 enforce wall boundary condition at
      ! bottom of the grid
      do L = philwb, phiupb
         do J = rlwb, rupb
            etot(J,zlwb-1,L)  =  etot(J,zlwb,L)
            s(J,zlwb-1,L) =   s(J,zlwb,L)
            a(J,zlwb-1,L) =   a(J,zlwb,L)
            t(J,zlwb-1,L) = - t(J,zlwb+1,L)
            t(J,zlwb,L)   =   0.0
         enddo
      enddo
   endif
endif

if( iam_on_edge ) then
   if( boundary_condition(2) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            etot(rupb+1,K,L)  =  etot(rupb,K,L)
            s(rupb+1,K,L) = - s(rupb-1,K,L)
            s(rupb,K,L)   =   0.0
            a(rupb+1,K,L) =   a(rupb,K,L)
            t(rupb+1,K,L) =   t(rupb,K,L)
         enddo
      enddo
   else if( boundary_condition(2) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            etot(rupb+1,K,L)  =  etot(rupb,K,L)
            s(rupb+1,K,L) = s(rupb,K,L)
            t(rupb+1,K,L) = t(rupb,K,L)
            a(rupb+1,K,L) = a(rupb,K,L)
         enddo
      enddo
   else
      ! dirichlet boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            s(rupb+1,K,L) = 0.0
            s(rupb,K,L)   = 0.5*(abs(s(rupb,K,L)) +     &
                                     s(rupb,K,L))
            a(rupb+1,K,L) = 0.0
            t(rupb+1,K,L) = 0.0
         enddo
      enddo
   endif
endif

if( iam_on_top ) then
   if( boundary_condition(3) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            etot(J,zupb+1,L)  =  etot(J,zupb,L)
            s(J,zupb+1,L) =   s(J,zupb,L)
            a(J,zupb+1,L) =   a(J,zupb,L)
            t(J,zupb+1,L) = - t(J,zupb-1,L)
            t(J,zupb,L)   = 0.0
         enddo
      enddo
   else if( boundary_condition(3) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            etot(J,zupb+1,L)  =  etot(J,zupb,L)
            s(J,zupb+1,L) = s(J,zupb,L)
            t(J,zupb+1,L) = t(J,zupb,L)
            a(J,zupb+1,L) = a(J,zupb,L)
         enddo
      enddo
   else
      ! dirichlet boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            s(J,zupb+1,L) = 0.0
            a(J,zupb+1,L) = 0.0
            t(J,zupb+1,L) = 0.0
            t(J,zupb,L)   = 0.5*(abs(t(J,zupb,L)) +        &
                                     t(J,zupb,L))
         enddo
      enddo
   endif
endif

return
end subroutine source_bc
