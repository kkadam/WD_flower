!**************************************************************
!*
!*  STATE
!*
!**************************************************************
subroutine state
implicit none
include 'runhydro.h'
!**************************************************************
!*
!    Using an ideal gas equation of state, assign the
!    pressure as a function of the density and internal
!    energy per unit mass
!
!      p  =  (gamma - 1) tau ** gamma
!
!    where tau = (rho * eps) ** 1 / gamma
!*
!**************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: etot
common /energy/ etot

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                    &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                         &
             in_neighbor, out_neighbor, root,                         &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                           &
           in_neighbor, out_neighbor, root,                           &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                    &
                        iam_on_bottom, iam_on_axis,                   &
                        iam_on_edge, down_neighbor,                   &
                        up_neighbor, in_neighbor,                     &
                        out_neighbor, root, column_num,               &
                        row_num, pe_grid, iam_root,                   &
                        REAL_SIZE, INT_SIZE

real :: a_hat, b_hat, eight_a_b_hat
common /whitedwarf/ a_hat, b_hat, eight_a_b_hat 

real, dimension(numr_dd,numz_dd,numphi) :: ppoly, pdeg
common /ztte/ ppoly,pdeg
!*
!**************************************************************
!*
!*   Local variables

real :: vel_r
real :: vel_z
real :: vel_p
real :: int_energy, vel_sq
real :: emin
REAL :: x_sq, x, arc_sinh,temp
REAL :: p_den_cutoff,pmin 
REAL, DIMENSION(numr_dd,numz_dd,numphi) :: udeg

integer :: J, K, L

!*
!**************************************************************

x = (1.0e-10/b_hat)**(1.0/3.0)
x_sq = x*x
arc_sinh = log(x+SQRT(1.0+x_sq))
pmin = a_hat*(x*(2.0*x_sq-3.0)*SQRT(1.0+x_sq)+3.0*arc_sinh)
emin= (a_hat*8*(x**3*SQRT(x_sq+1.0d0)-1.0d0)-pmin)!/rho(J,K,L)
p_den_cutoff = 1.0d-8

!ORG p = (gamma - 1.0) * tau**gamma

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
        x = (rho(J,K,L)/b_hat)**(1.0/3.0)
        x_sq = x*x
        arc_sinh = log(x+SQRT(1.0+x_sq))
        pdeg(J,K,L) = a_hat*(x*(2.0*x_sq-3.0)*SQRT(1.0+x_sq)+3.0*arc_sinh)
        udeg(J,K,L) = (a_hat*8*x**3*(SQRT(x_sq+1.0d0)-1.0d0)-pdeg(J,K,L))!/rho(J,K,L)
        IF (pdeg(J,K,L).LT.pmin) pdeg(J,K,L) = pmin
        IF (udeg(J,K,L).LT.emin) udeg(J,K,L) = emin
     enddo
   enddo
enddo



do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         vel_r = 0.5 * ( u(J,K,L) + u(J+1,K,L) )
         vel_z = 0.5 * ( w(J,K,L) + w(J,K+1,L) )
         vel_p = 0.5 * ( jn(J,K,L) + jn(J,K,L+1) ) * rhfinv(J)
         vel_sq = vel_r * vel_r + vel_z * vel_z + vel_p * vel_p
         int_energy = etot(J,K,L) - 0.5 * rho(J,K,L) * vel_sq - udeg(J,K,L)
         if ( int_energy <= 0.0d0 ) then
            int_energy = 0.0d0
         endif
         ppoly(J,K,L) = (gamma - 1.0) * int_energy*(1.0d0-p_den_cutoff/rho(J,K,L))
         IF (ppoly(J,K,L).LT.0.0d0) ppoly(J,K,L) = 0.0d0
         p(J,K,L) = ppoly(J,K,L) + pdeg(J,K,L)
!         p(J,K,L) = pdeg(J,K,L)
      enddo
   enddo
enddo

!L = phiupb
do K = zlwb, zupb
   do J = rlwb, rupb
      vel_r = 0.5 * ( u(J,K,phiupb) + u(J+1,K,phiupb) )
      vel_z = 0.5 * ( w(J,K,phiupb) + w(J,K+1,phiupb) )
      vel_p = 0.5 * ( jn(J,K,phiupb) + jn(J,K,philwb) ) * rhfinv(J)
      vel_sq = vel_r * vel_r + vel_z * vel_z + vel_p * vel_p
      int_energy = etot(J,K,phiupb) - 0.5 * rho(J,K,phiupb) * vel_sq - udeg(J,K,phiupb)
      if ( int_energy <= 0.0d0 ) then
         int_energy = 0.0d0
      endif
      ppoly(J,K,phiupb) = (gamma - 1.0) * int_energy*(1.0d0-p_den_cutoff/rho(J,K,phiupb))
      IF (ppoly(J,K,phiupb).LT.0.0d0) ppoly(J,K,phiupb) = 0.0d0
      p(J,K,phiupb) = ppoly(J,K,phiupb) + pdeg(J,K,phiupb)
!      p(J,K,phiupb) = pdeg(J,K,phiupb)
   enddo
enddo

if ( iam_on_bottom ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
        p(J,zlwb-1,L) = p(J,zlwb,L)
      enddo
   enddo
endif

if ( iam_on_top ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
         p(J,zupb+1,L) = p(J,zupb,L)
      enddo
   enddo
endif

if ( iam_on_edge ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         p(rupb+1,K,L) = p(rupb,K,L)
      enddo
   enddo
endif

call comm(p)

return
end subroutine state
