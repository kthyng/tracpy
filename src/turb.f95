subroutine turbuflux(ia,ja,ka,rr,dtmin,ah,imt,jmt,km,uflux,vflux,wflux,ff,do3d,upr)
!====================================================================
! computes the paramterised turbulent velocities u' and v' into upr
!
!  Input:
!    ia,ja,ka       : original position in integers
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    dtmin          : time step based on the interpolation step between model output times.
!                   : sets a limit on the time step that a drifter can go. (seconds)
!    ah             : horizontal diffusion in m^2/s. Input parameter.
!                   : For -turb and -diffusion flags
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!    wflux          : w velocity (vertical) flux field, two time steps [kxt]
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    do3d           : Flag to set whether to use 3d velocities or not
!
!  Output:
!    upr            : parameterized turbulent velocities u', v', and w'
!                     optional because only used if using turb flag for diffusion
!                     size [6,2]. The 2nd dimension is for two time steps.
!                     The 1st dimension is: [u'_ia,u'_ia-1,v'_ja,v'_ja-1,w'_ka,w'_ka-1]
!
!  Other parameters used in function:
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!    im,jm          : generic index for grid index -1 in i and j directions 
!    n              : Index for looping through times (1 and 2)
!    uv             : array for holding time-interpolated horizontal fluxes.
!                   : 4 of the 12 spots are used (for u and v forward and backward in 
!                     space). The other 8 spots aren't used but are there for now.
!    localW         : holds the vertical flux below the drifter grid cell
!    qran           : array of random numbers that gets multipled with other things.
!                   : Only spots 1 and 3 of the 12 are used.
!    amp            : amplitude of turbulence correction
!====================================================================

  
implicit none
 
integer,        intent(in)                                      :: ia,ja,ka, ff, imt,jmt,km
integer,        intent(in)                                      :: do3d
real*8,         intent(in)                                      :: rr, dtmin,ah
real(kind=8),   intent(in),     dimension(imt-1,jmt,km,2)       :: uflux
real(kind=8),   intent(in),     dimension(imt,jmt-1,km,2)       :: vflux
real(kind=8),   intent(in),     dimension(0:km,2)               :: wflux
real*8,         intent(out),    dimension(6,2)                  :: upr  
integer                                                         :: im,jm,n
integer                                                         :: nsm=1,nsp=2
real*8                                                          :: uv(12),rg,localW !,en
real*8                                                          :: qran(12),amp

!amp=ah/sqrt(dt)
!amp=ah/dtmin
!amp=ah/sqrt(dtmin)
amp=ah/( dtmin**(1./3.) )
!print *,ah,dt,amp
!stop 4967
!  call random_number(qran) ! generates a random number between 0 and 1
  
! random generated numbers between 0 and 1
!do n=1,12
!qran(n)=rand()
CALL RANDOM_NUMBER (qran)
!enddo
	
!  qran=2.*qran-1. ! === Max. amplitude of turbulence with random numbers between -1 and 1
                   ! === (varies with the same aplitude as the mean vel)
  !qran=1.*qran-0.5  ! Reduced amplitude of turb.
!  qran=4.*qran-2. ! ===  amplitude of turbulence with random numbers between -2 and 2
!  qran=8.*qran-4. ! ===  amplitude of turbulence with random numbers between -2 and 2 only works with   (upr)(:,1)=upr(:,2)
!  qran=30.*qran-15. ! ===  amplitude of turbulence with random numbers between -2 and 2 only works with   upr(:,1)=upr(:,2)
qran=amp*qran-0.5*amp ! ===  amplitude of turbulence with random numbers between -2 and 2 only works with   upr(:,1)=upr(:,2)
    
rg=1.d0-rr

im=ia-1
! KMT commented this. Isn't this a periodic condition? Why is this hardwired in and
! why isn't there one for the j direction too?
! Is there a need for another sort of boundary condition here, in place of periodic?
!   if(im.eq.0) im=IMT 
jm=ja-1
  
! time interpolated velocities
uv(1)=(rg*uflux(ia,ja,ka,nsp)+rr*uflux(ia,ja,ka,nsm))*ff ! western u
uv(2)=(rg*uflux(im,ja,ka,nsp)+rr*uflux(im,ja,ka,nsm))*ff ! eastern u
uv(3)=(rg*vflux(ia,ja,ka,nsp)+rr*vflux(ia,ja,ka,nsm))*ff ! northern v
uv(4)=(rg*vflux(ia,jm,ka,nsp)+rr*vflux(ia,jm,ka,nsm))*ff ! southern v

!   upr(:,1)=upr(:,2) ! store u' from previous time iteration step (this makes it unstable!)
   
! multiply the time interpolated velocities by random numbers
   
! 4 different random velocities at the 4 walls
!   upr(:,2)=uv(:)*rand(:) 
! or same random velocities at the eastern and western as well as nothern and southern

! print *,'qran(1)=',qran(1),' qran(3)=',qran(3) ! Want qran's to be between -1 and 1
upr(1,2)=uv(1)*qran(1)
upr(2,2)=uv(2)*qran(1)
upr(3,2)=uv(3)*qran(3)
upr(4,2)=uv(4)*qran(3)

!  upr(:,1)=upr(:,2) !  impose same velocities for t-1 and t (this makes it unstable. but why? K.Doos)
  
if(do3d==1) then

    ! === Calculates the w' at the top of the box from 
    ! === the divergence of u' and v' in order
    ! === to respect the continuity equation even 
    ! === for the turbulent velocities
    do n=1,2
     
     !Detta ser ut som en bugg   ! t2
     !  upr(5,n) = w(ka-1) - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     !  upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     
     ! === The vertical velocity is calculated from 
     ! === the divergence of the horizontal turbulence !v2
     ! === The turbulent is evenly spread between the 
     ! === top and bottom of box if not a bottom box

! #ifdef full_wflux
!     localW=wflux(ia,ja,ka-1,nsm)
! #else
        localW=wflux(ka-1,nsm)
! #endif
    
        if(localW.eq.0.d0) then
            upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
            upr(6,n) = 0.d0
        else
            upr(5,n) = - 0.5d0 * ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
            upr(6,n) = - upr(5,n)
        endif

    enddo
endif
 
 
 return
end subroutine turbuflux
!_______________________________________________________________________
