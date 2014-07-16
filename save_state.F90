!*************************************************************
!*
!*  SAVE_STATE
!*
!*************************************************************
subroutine save_state
implicit none
include 'runhydro.h'
!*************************************************************
!*
!  save copies s, t, a, and rho into temp arrays so they
!  are not changed by the first call to advect in a 
!  timestep cycle.  The first call to advect in a cycle
!  is only used to get time centered velocities to do a 
!  full timestep of advection.
!* 
!*************************************************************
!*
!*   Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: etot
common /energy/ etot

real, dimension(numr_dd,numz_dd,numphi) :: s1, t1, a1, rho1, etot1
common /save/ s1, t1, a1, rho1, etot1

!*
!*************************************************************
!*
!*  Local variables

integer :: J, K, L

!*
!*************************************************************

do L = 1, numphi
   do K = 1, numz_dd
      do J = 1, numr_dd
         s1(J,K,L) = s(J,K,L)
         t1(J,K,L) = t(J,K,L)
         a1(J,K,L) = a(J,K,L)
         rho1(J,K,L) = rho(J,K,L)
         etot1(J,K,L) = etot(J,K,L)
      enddo
   enddo
enddo
     
return
end subroutine save_state
