!***********************************************************************
!*
!*  SCFIN_FRAC
!*
!***********************************************************************
subroutine scfin
implicit none
include 'runhydro.h'
!***********************************************************************
!*
!  scfin reads in the three dimensional arrays for the density,
!  anglular momentum density and radial momentum density.  Then
!  we set the internal energy density array from the mass density
!  and the appropriate polytropic constant for each star.
!
!  Assume that an external program generates the initial conditions
!  for all data on all processors.  This means that guard cell
!  values are assumed to be valid in the input file and appropriate
!  boundary conditions are enforced as well.
!
!  >>NOTE<< for isym = 1, the input files have to extend the entire
!           range of z values.  Communication to replicate data
!           from upper half space to lower half space would be
!           too big a pain in the neck
!
!  5/17/2000 modified to set up the fluid's mass fraction array 
!  for an initial model.  The fluid is all star1 if it is left
!  of x = 0.1 initially... hard coded and stupid.
!*
!***********************************************************************
!*
!*  Global Variable

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: etot
common /energy/ etot

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi 

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr_dd,numz_dd,numphi) :: ppoly, pdeg
common /ztte/ ppoly,pdeg

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

real :: a_hat, b_hat, eight_a_b_hat
common /whitedwarf/ a_hat, b_hat, eight_a_b_hat 

!*
!***********************************************************************
!*
!*  Local Variable

real, dimension(numr_dd,numz_dd,numphi) :: eps
real, dimension(numr_dd,numz_dd,numphi) :: udeg 

real :: temp 

real :: r_velocity
real :: z_velocity
real :: phi_velocity
real :: vel_sq
 
integer :: record_length, three_numphi_by_four

integer :: J, K, L

REAL :: x_sq, x, arc_sinh

!*
!***********************************************************************
!  initialize the local variables
eps = 0.0
temp = 0.0
three_numphi_by_four = 3 * numphi_by_four
inquire(iolength=record_length) rho

open(unit=21, file='input/density', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)

open(unit=22, file='input/ang_mom', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)

open(unit=23, file='input/rad_mom', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)
       
read(21,rec=iam+1) rho

read(22,rec=iam+1) a

read(23,rec=iam+1) s

close(21)
close(22)
close(23)

! 10/08/2003 : Modified x to -0.2
! set up the mass fraction array
!
! if fluid is of density greater then 10 times the vacuum level 
! and is to the left of x = -0.2 then call it fluid of type 1
!
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         frac1(J,K,L) = 0.0
         frac2(J,K,L) = 0.0
         if( rho(J,K,L) > 10.0 * densmin ) then
            if( rhf(J) * cos_cc(L) >= 0.15 ) then
               frac1(J,K,L) = 1.0
            else
               frac2(J,K,L) = 1.0
            endif
         endif 
      enddo
   enddo
enddo


! given the density generate the internal energy per unit mass
!
!  p  =  kappa * rho ** gamma  =  (gamma - 1) * rho * eps
!
!  so that
! 
! eps = kappa * rho**(1/n) / (gamma - 1)
!
temp = gamma - 1.0

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         if( rhf(J) * cos_cc(L) >= 0.15 ) then
           eps(J,K,L) = (kappa2 * rho(J,K,L)**(1.0/pin))/temp
         else
           eps(J,K,L) = (kappa1 * rho(J,K,L)**(1.0/pin))/temp
         endif
      enddo
   enddo
enddo

! now convert eps to be the entropy tracer, tau
! 
!  tau  =  (eps * rho) ** (1/gamma)
!
! and tau is related to the entropy of the fluid as follows
!
!  s  =  c_p * ln( tau / rho ) + constant
!
!  for adiabatic flow tau has no source term in the eulerian
!  fluid equations (just like rho).  this is not true if
!  the effects of viscosity (artificial or real) are included. 
tau = (eps*rho)**gammainv


do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
        x = (rho(J,K,L)/b_hat)**(1.0/3.0)
        x_sq = x*x
        arc_sinh = log(x+SQRT(1.0+x_sq))
        pdeg(J,K,L) = a_hat*(x*(2.0*x_sq-3.0)*SQRT(1.0+x_sq)+3.0*arc_sinh)
        udeg(J,K,L) = (a_hat*8*x**3*(SQRT(x_sq+1)-1)-pdeg(J,K,L))!/rho(J,K,L)
        IF (pdeg(J,K,L).LT.0) pdeg(J,K,L) = 0.0
     enddo
   enddo
enddo



do L = philwb, phiupb-1
   do K = zlwb-1, zupb
      do J = rlwb-1, rupb
         r_velocity = rhfinv(J) * ( s(J,K,L)   * ( rhf(J) - 0.25 * dr ) + &
                                    s(J+1,K,L) * ( rhf(J) + 0.25 * dr ) )
         r_velocity = r_velocity / rho(J,K,L)
         z_velocity = 0.5 * ( t(J,K,L) + t(J,K+1,L) )
         z_velocity = z_velocity / rho(J,K,L)
         phi_velocity = 0.5 * ( a(J,K,L) + a(J,K,L+1) )
         phi_velocity = phi_velocity / rho(J,K,L) * rhfinv(J)
         vel_sq = r_velocity*r_velocity + z_velocity*z_velocity + phi_velocity*phi_velocity
         etot(J,K,L) = 0.5*( rho(J,K,L) * vel_sq ) + udeg(J,K,L) !rho(J,K,L) * eps(J,K,L)
      enddo
   enddo
enddo
L = phiupb
do K = zlwb-1, zupb
   do J = rlwb-1, rupb
      r_velocity = rhfinv(J) * ( s(J,K,phiupb)   * ( rhf(J) - 0.25 * dr ) + &
                                 s(J+1,K,phiupb) * ( rhf(J) + 0.25 * dr ) )
      r_velocity = r_velocity / rho(J,K,phiupb)
      z_velocity = 0.5 * ( t(J,K,phiupb) + t(J,K+1,phiupb) )
      z_velocity = z_velocity / rho(J,K,phiupb)
      phi_velocity = 0.5 * ( a(J,K,phiupb) + a(J,K,philwb) )
      phi_velocity = phi_velocity / rho(J,K,phiupb) * rhfinv(J)
      vel_sq = r_velocity*r_velocity + z_velocity*z_velocity + phi_velocity*phi_velocity
      etot(J,K,phiupb) = 0.5*( rho(J,K,phiupb) * vel_sq ) + udeg(J,K,phiupb) !rho(J,K,phiupb) * eps(J,K,phiupb)
   enddo
enddo

call comm(etot)

return
end subroutine scfin
