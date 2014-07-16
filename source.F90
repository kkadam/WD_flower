!************************************************************************
!*  
!*  SOURCE
!*
!************************************************************************
subroutine source(source_time) 
implicit none
include 'runhydro.h'
!************************************************************************
!*
!   source applies Lagrangian source terms, which for adiabatic flow,
!   only appear in the three equations of motion.  therefore
!   one has to update the three conserved momenta.
!
!   s  =  s  -  dt ( <rho> dPhieff/dR  - <a>^2/(<rho>R^3) - 2 Omegaframe<a>/R )
!                          ^                   ^                   ^
!                          |                   |                   |
!                   gravity, pressure and   curvature of        Coriolis 
!                      centrifugal       cylindrical coords
!
!   t  =  t  -  dt <rho> dPhieff/dZ
!                          ^
!                          |
!                    gravity and pressure
!                       
!   a  =  a  -  dt ( rho dPhieff/dphi  +  2 Omaegaframe s R)
!                          ^                      ^ 
!                          |                      |
!                   gravity and pressure       Coriolis
!
!  where,
!
!  Phieff = (pin + 1) p / rho + Phi  - 0.5 Omegaframe^2 R^2
!
!         =           H     +   Phi  - 0.5 Omegaframe^2 R^2
!
!         =  the effective hydrodynamic potential 
!
!         =  the Bernoulli function
!
!   this is a little tricky here, to center the sourcing of s and a
!   in time need to take different actions depending on whether this
!   is the first or second call to source in a timestep cycle.
!   
!   -->  On the first call: 
!
!        s(n + 1/2S)  -  s(n)
!        --------------------  =  f(rho(n), phieff(n), tempa)
!              0.5 dt(n)
!
!        a(n + 1/2S)  -  a(n)
!        --------------------  =  g(rho(n), phieff(n), temps)
!              0.5 dt(n)
!
!   -->  On the second call: 
!
!        temps  -  s(n + 1/2S + A)
!        -------------------------  =  f(rho(n+1), phieff(n+1), a(n + 1/2S + A))
!                 0.5 dt(n)
!
!        tempa  -  a(n + 1/2S + A)
!        -------------------------  =  g(rho(n+1), phieff(n+1), s(n + 1/2S + A))
!                 0.5 dt(n)
!
!        s(n+1)  -  s(n + 1/2S + A)
!        --------------------------  =  f(rho(n+1), phieff(n+1), tempa)
!                 0.5 dt(n)
!
!        a(n+1)  -  a(n + 1/2S + A)
!        --------------------------  =  g(rho(n+1), phieff(n+1), temps)
!                 0.5 dt(n)
!
!    the first two assignment statements in the second call
!    are actually performed by the surbroutine 'make_source_temp', this 
!    allows the same code to perform the nescessary operations for both
!    the first and second calls to source in a timestep cycle
!
!  04/05/2003 : Changed H (= p/rho) to p
!
!  05/24/2003 : Modified for com motion (only com acceleration in the radial 
!               direction).
!
!  06/09/2003 : Modified the radial and azimuthal equations to include
!               the radial and azimuthal component of the acceleration
!               of the com. Refer to subroutine com_vel_accn.F for 
!               details on how the radial and azimuthal components of
!               aceleration of the com are calculated.
!
!*
!************************************************************************
!*
!*  Subroutine Arguments

integer :: source_time

!*
!************************************************************************
!*
!*   Global variables

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: etot
common /energy/ etot

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega


!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: temp, rhjkl, phieff

real :: p_mr, p_pr
real :: p_mz, p_pz
real :: p_mp, p_pp
real :: rho_mr, rho_pr
real :: rho_mz, rho_pz
real :: rho_mp, rho_pp
real :: divpv
real :: gradphi_r
real :: gradphi_z
real :: gradphi_p
integer :: J, K, L

!*
!************************************************************************
!  initialize the local variables

! 04/05/2003 : Modified (H = p/rho) to p
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         phieff(J,K,L) = pot(J,K,L) - 0.5*omega_frame*            &
                         omega_frame*rhf(J)*rhf(J)

!         phieff(J,K,L) = (pin+1.0)*p(J,K,L)/rho(J,K,L) +          &
!                         pot(J,K,L) - 0.5*omega_frame*            &
!                         omega_frame*rhf(J)*rhf(J)
      enddo
   enddo
enddo

! source terms for  total energy equation
do L = philwb+1, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb

         p_mr = 0.5 * ( p(J,K,L) + p(J-1,K,L) )
         p_pr = 0.5 * ( p(J,K,L) + p(J+1,K,L) )
         p_mz = 0.5 * ( p(J,K,L) + p(J,K-1,L) )
         p_pz = 0.5 * ( p(J,K,L) + p(J,K+1,L) )
         p_mp = 0.5 * ( p(J,K,L) + p(J,K,L-1) )
         p_pp = 0.5 * ( p(J,K,L) + p(J,K,L+1) )

         rho_mr = 0.5 * rinv(J) * ( rho(J,K,L)   * ( r(J) + 0.25 * dr ) + &
                                    rho(J-1,K,L) * ( r(J) - 0.25 * dr ) )
         rho_pr = 0.5 * rinv(J+1) * ( rho(J+1,K,L) * ( r(J+1) + 0.25 * dr ) + &
                                      rho(J,K,L)   * ( r(J+1) - 0.25 * dr ) )
         rho_mz = 0.5 * ( rho(J,K,L)   + rho(J,K-1,L) )
         rho_pz = 0.5 * ( rho(J,K+1,L) + rho(J,K,L) )
         rho_mp = 0.5 * ( rho(J,K,L)   + rho(J,K,L-1) )
         rho_pp = 0.5 * ( rho(J,K,L+1) + rho(J,K,L) )

         if ( r(J) > 0.0 ) then
            divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,L) / rho_pr - &
                                          r(J)   * p_mr * s(J,K,L)   / rho_mr )
         else
            divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,L) / rho_pr )
         endif

         divpv = divpv + rhfinv(J) * dphiinv * ( p_pp * rhfinv(J) * a(J,K,L+1) / rho_pp - &
                                                 p_mp * rhfinv(J) * a(J,K,L)   / rho_mp )

         divpv = divpv + dzinv * ( p_pz * t(J,K+1,L) / rho_pz - p_mz * t(J,K,L) / rho_mz )

         gradphi_r = 0.5 * drinv * rhfinv(J) * ( s(J+1,K,L) * (rhf(J) + 0.25*dr) + &
                                                 s(J,K,L)   * (rhf(J) - 0.25*dr) )
         gradphi_r = 0.5 * gradphi_r * ( phieff(J+1,K,L) - phieff(J-1,K,L) )

         gradphi_z = 0.5 * ( t(J,K+1,L) + t(J,K,L) )
         gradphi_z = 0.5 * dzinv * gradphi_z * ( phieff(J,K+1,L) - phieff(J,K-1,L) )

         gradphi_p = 0.5 * ( a(J,K,L+1) + a(J,K,L) ) * rhfinv(J)
         gradphi_p = 0.5 * dphiinv * rhfinv(J) * gradphi_p * ( phieff(J,K,L+1) - phieff(J,K,L-1) )

         etot(J,K,L) = etot(J,K,L) - dt*( divpv + gradphi_r + gradphi_z + gradphi_p )
      enddo
   enddo
enddo

!L = philwb
do K = zlwb, zupb
   do J = rlwb, rupb
      p_mr = 0.5 * ( p(J,K,philwb) + p(J-1,K,philwb) )
      p_pr = 0.5 * ( p(J,K,philwb) + p(J+1,K,philwb) )
      p_mz = 0.5 * ( p(J,K,philwb) + p(J,K-1,philwb) )
      p_pz = 0.5 * ( p(J,K,philwb) + p(J,K+1,philwb) )
      p_mp = 0.5 * ( p(J,K,philwb) + p(J,K,phiupb) )
      p_pp = 0.5 * ( p(J,K,philwb) + p(J,K,philwb+1) )

      rho_mr = 0.5 * rinv(J) * ( rho(J,K,philwb)   * ( r(J) + 0.25 * dr ) + &
                                 rho(J-1,K,philwb) * ( r(J) - 0.25 * dr ) )
      rho_pr = 0.5 * rinv(J+1) * ( rho(J+1,K,philwb) * ( r(J+1) + 0.25 * dr ) + &
                                   rho(J,K,philwb)   * ( r(J+1) - 0.25 * dr ) )
      rho_mz = 0.5 * ( rho(J,K,philwb)   + rho(J,K-1,philwb) )
      rho_pz = 0.5 * ( rho(J,K+1,philwb) + rho(J,K,philwb) )
      rho_mp = 0.5 * ( rho(J,K,philwb)   + rho(J,K,phiupb) )
      rho_pp = 0.5 * ( rho(J,K,philwb+1) + rho(J,K,philwb) )

      if ( r(J) > 0.0 ) then
         divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,philwb) / rho_pr - &
                                       r(J)   * p_mr * s(J,K,philwb)   / rho_mr )
      else
         divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,philwb) / rho_pr )
      endif

      divpv = divpv + rhfinv(J) * dphiinv * ( p_pp * rhfinv(J) * a(J,K,philwb+1) / rho_pp - &
                                              p_mp * rhfinv(J) * a(J,K,philwb)   / rho_mp )

      divpv = divpv + dzinv * ( p_pz * t(J,K+1,philwb) / rho_pz - p_mz * t(J,K,philwb) / rho_mz )

      gradphi_r = 0.5 * drinv * rhfinv(J) * ( s(J+1,K,philwb) * (rhf(J) + 0.25*dr) + &
                                              s(J,K,philwb)   * (rhf(J) - 0.25*dr) )
      gradphi_r = 0.5 * gradphi_r * ( phieff(J+1,K,philwb) - phieff(J-1,K,philwb) )

      gradphi_z = 0.5 * ( t(J,K+1,philwb) + t(J,K,philwb) )
      gradphi_z = 0.5 * dzinv * gradphi_z * ( phieff(J,K+1,philwb) - phieff(J,K-1,philwb) )

      gradphi_p = 0.5 * ( a(J,K,philwb+1) + a(J,K,philwb) ) * rhfinv(J)
      gradphi_p = 0.5 * dphiinv * rhfinv(J) * gradphi_p * ( phieff(J,K,philwb+1) - phieff(J,K,phiupb) )

      etot(J,K,philwb) = etot(J,K,philwb) - dt*( divpv + gradphi_r + gradphi_z + gradphi_p )
   enddo
enddo

!L = phiupb
do K = zlwb, zupb
   do J = rlwb, rupb
      p_mr = 0.5 * ( p(J,K,phiupb) + p(J-1,K,phiupb) )
      p_pr = 0.5 * ( p(J,K,phiupb) + p(J+1,K,phiupb) )
      p_mz = 0.5 * ( p(J,K,phiupb) + p(J,K-1,phiupb) )
      p_pz = 0.5 * ( p(J,K,phiupb) + p(J,K+1,phiupb) )
      p_mp = 0.5 * ( p(J,K,phiupb) + p(J,K,phiupb-1) )
      p_pp = 0.5 * ( p(J,K,phiupb) + p(J,K,philwb) )

      rho_mr = 0.5 * rinv(J) * ( rho(J,K,phiupb)   * ( r(J) + 0.25 * dr ) + &
                                 rho(J-1,K,phiupb) * ( r(J) - 0.25 * dr ) )
      rho_pr = 0.5 * rinv(J+1) * ( rho(J+1,K,phiupb) * ( r(J+1) + 0.25 * dr ) + &
                                   rho(J,K,phiupb)   * ( r(J+1) - 0.25 * dr ) )
      rho_mz = 0.5 * ( rho(J,K,phiupb)   + rho(J,K-1,phiupb) )
      rho_pz = 0.5 * ( rho(J,K+1,phiupb) + rho(J,K,phiupb) )
      rho_mp = 0.5 * ( rho(J,K,phiupb)   + rho(J,K,phiupb-1) )
      rho_pp = 0.5 * ( rho(J,K,philwb)   + rho(J,K,phiupb) )

      if ( r(J) > 0.0 ) then
         divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,phiupb) / rho_pr - &
                                       r(J)   * p_mr * s(J,K,phiupb)   / rho_mr )
      else
         divpv = rhfinv(J) * drinv * ( r(J+1) * p_pr * s(J+1,K,phiupb) / rho_pr )
      endif

      divpv = divpv + rhfinv(J) * dphiinv * ( p_pp * rhfinv(J) * a(J,K,philwb) / rho_pp - &
                                              p_mp * rhfinv(J) * a(J,K,phiupb) / rho_mp )

      divpv = divpv + dzinv * ( p_pz * t(J,K+1,phiupb) / rho_pz - p_mz * t(J,K,phiupb) / rho_mz )

      gradphi_r = 0.5 * drinv * rhfinv(J) * ( s(J+1,K,phiupb) * (rhf(J) + 0.25*dr) + &
                                              s(J,K,phiupb)   * (rhf(J) - 0.25*dr) )
      gradphi_r = 0.5 * gradphi_r * ( phieff(J+1,K,phiupb) - phieff(J-1,K,phiupb) )

      gradphi_z = 0.5 * ( t(J,K+1,phiupb) + t(J,K,phiupb) )
      gradphi_z = 0.5 * dzinv * gradphi_z * ( phieff(J,K+1,phiupb) - phieff(J,K-1,phiupb) )

      gradphi_p = 0.5 * ( a(J,K,philwb) + a(J,K,phiupb) ) * rhfinv(J)
      gradphi_p = 0.5 * dphiinv * rhfinv(J) * gradphi_p * ( phieff(J,K,philwb) - phieff(J,K,phiupb-1) )

      etot(J,K,phiupb) = etot(J,K,phiupb) - dt*( divpv + gradphi_r + gradphi_z + gradphi_p )
   enddo
enddo

!  center rho at the vertical face centers
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * ( rho(J,K,L) + rho(J,K-1,L) )
      enddo
   enddo
enddo

! 06/17/2003 :
! 04/05/2003 :
!  source t, the vertical momentum density
if (tstep .le. 1000) then
!if (tstep .le. 100000000) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            t(J,K,L) = t(J,K,L) - dt*dzinv*                       &
                       ( (p(J,K,L)-p(J,K-1,L)) +                  &
                         rhjkl(J,K,L)*                            &
                         (phieff(J,K,L) - phieff(J,K-1,L)) )
         enddo
      enddo
   enddo
else
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            t(J,K,L) = t(J,K,L) - dt * (dzinv*                    &
                       ( (p(J,K,L)-p(J,K-1,L)) +                  &
                         rhjkl(J,K,L)*                            &
                         (phieff(J,K,L) - phieff(J,K-1,L)) )      &
                        + rhjkl(J,K,L)*cylin_a_com(3,L) )
         enddo
      enddo
   enddo
endif

!  source s the radial momentum density and a the angular
!  momentum density

if (source_time == 2) call make_source_temp

!  center rho at the radial face centers
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * rinv(J) *                         &
                        ( rho(J,K,L)*(r(J) + 0.25*dr) +         &
                          rho(J-1,K,L)*(r(J) - 0.25*dr) )
      enddo
   enddo
enddo

! center tempa at the radial face centers
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = 0.25 * rinv(J) *                         &
                      ( (tempa(J,K,L) + tempa(J,K,L+1)) *       &
                       (r(J) + 0.25*dr) +                       &
                         (tempa(J-1,K,L) + tempa(J-1,K,L+1)) *  &
                       (r(J) - 0.25*dr) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      temp(J,K,phiupb) = 0.25 * rinv(J) *                                &
                        ( (tempa(J,K,phiupb) + tempa(J,K,philwb)) *      &
                         (r(J) + 0.25*dr) +                              &
                          (tempa(J-1,K,phiupb) + tempa(J-1,K,philwb)) *  &
                         (r(J) - 0.25*dr) )
   enddo
enddo

! 06/17/2003 :
! 06/09/2003 :
! 05/24/2003 : Modified for com motion (only com acceleration in the radial 
!              direction).
! 04/05/2003 : Changed H (= p/rho) to p

! apply source terms to S
if (tstep .le. 1000) then
!if (tstep .le. 100000000) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( rhjkl(J,K,L) > 0.0 ) then
               s(J,K,L) = s(J,K,L) - dt*(                          &
                    drinv*( (p(J,K,L)-p(J-1,K,L)) +                &
                    rhjkl(J,K,L)*(phieff(J,K,L)-phieff(J-1,K,L)) ) &
                    - temp(J,K,L)*temp(J,K,L)*rinv(J)*             &
                      rinv(J)*rinv(J)/rhjkl(J,K,L)                 &
                    - 2.0*omega_frame*temp(J,K,L)*rinv(J) )
            else
               s(J,K,L) = 0.0
            endif
         enddo
      enddo
   enddo
else
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( rhjkl(J,K,L) > 0.0 ) then
               s(J,K,L) = s(J,K,L) - dt*(                          &
                    drinv*( (p(J,K,L)-p(J-1,K,L)) +                &
                    rhjkl(J,K,L)*(phieff(J,K,L)-phieff(J-1,K,L)) ) &
                    - temp(J,K,L)*temp(J,K,L)*rinv(J)*             &
                      rinv(J)*rinv(J)/rhjkl(J,K,L)                 &
                    - 2.0*omega_frame*temp(J,K,L)*rinv(J)          &
                    + rhjkl(J,K,L)* (cylin_a_com(1,L)              &
                      - omega_frame*omega_frame*R_com              &
                      - 2*omega_frame*cylin_v_com(2,L)) )
            else
               s(J,K,L) = 0.0
            endif
         enddo
      enddo
   enddo
endif

! center rho at the azimuthal face centers
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * (rho(J,K,L) + rho(J,K,L-1))
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      rhjkl(J,K,philwb) = 0.5 * (rho(J,K,philwb) + rho(J,K,phiupb))
   enddo
enddo

! center temps at the azimuhtal face centers
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = 0.25 * rhfinv(J) *                        &
                      ( (temps(J+1,K,L) + temps(J+1,K,L-1))*     &
                       (rhf(J) + 0.25*dr) +                      &
                        (temps(J,K,L) + temps(J,K,L-1))*         &
                       (rhf(J) - 0.25*dr) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      temp(J,K,philwb) = 0.25 * rhfinv(J) *                                &
                        ( (temps(J+1,K,philwb) + temps(J+1,K,phiupb))*     &
                         (rhf(J) + 0.25*dr) +                              &
                          (temps(J,K,philwb) + temps(J,K,phiupb))*         &
                         (rhf(J) - 0.25*dr) )
   enddo
enddo

! 06/17/2003 :
! 06/09/2003 : Modified to include azimuthal component of the com 
!              acceleration (cylin_a_com(2))
! 04/05/2003 : Changed H (= p/rho) to p
! apply source terms to A

if (tstep .le. 1000) then
!if (tstep .le. 10000000) then
    do L = philwb+1, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             a(J,K,L) = a(J,K,L) - dt*(                              &
                    dphiinv*( (p(J,K,L)-p(J,K,L-1)) +                &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J,K,L-1)) ) &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,L) )
          enddo
       enddo
    enddo
    do K = zlwb, zupb
       do J = rlwb, rupb
          a(J,K,philwb) = a(J,K,philwb) - dt*(                      &
                      dphiinv*( (p(J,K,philwb)-p(J,K,phiupb)) +     &  
                      rhjkl(J,K,philwb) *                           &
                      (phieff(J,K,philwb) - phieff(J,K,phiupb)) )   &
                      + 2.0*omega_frame*rhf(J)*temp(J,K,philwb) )
       enddo
    enddo
else
    do L = philwb+1, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             a(J,K,L) = a(J,K,L) - dt*(                               &
                    dphiinv*( (p(J,K,L)-p(J,K,L-1)) +                 &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J,K,L-1)) )  &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,L)              &
                    + rhjkl(J,K,L)*rhf(J)* (cylin_a_com(2,L)          &
                       + 2*omega_frame*cylin_v_com(1,L)) )
          enddo
       enddo
    enddo
    do K = zlwb, zupb
       do J = rlwb, rupb
          a(J,K,philwb) = a(J,K,philwb) - dt*(                        &
                    dphiinv*( (p(J,K,philwb)-p(J,K,phiupb)) +         &
                    rhjkl(J,K,philwb) *                               &
                    (phieff(J,K,philwb) - phieff(J,K,phiupb)) )       &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,philwb)         &
                    + rhjkl(J,K,philwb)*rhf(J)*                       &
                      (cylin_a_com(2,philwb)                          &
                      + 2*omega_frame*cylin_v_com(1,philwb)) )
       enddo
    enddo
endif
    
call source_bc

!  share these updated values at processor edges with neighbors
call comm(s)
call comm(t)
call comm(a)
call comm(etot)

return
end subroutine source
