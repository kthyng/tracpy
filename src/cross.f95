subroutine cross(ijk,ia,ja,ka,r0,sp,sn,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb,upr)
  
!====================================================================
! subroutine to compute time (sp,sn) when trajectory 
! crosses face of box (ia,ja,ka) 
! two crossings are considered for each direction:  
! east and west for longitudinal directions, etc.  
! 
!  Input:
!
!    ijk            : considered direction (1=zonal,2=meridional,3=vertical)
!    ia,ja,ka       : original position in integers
!    r0             : original non-dimensional position in the ijk-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!    wflux          : w velocity (vertical) flux field, two time steps [kxt]
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    do3d           : Flag to set whether to use 3d velocities or not
!    doturb         : Flag to set whether or not to use turb/diff and which kind if so
!    upr            : parameterized turbulent velocities u', v', and w'
!                     optional because only used if using turb flag for diffusion
!                     size [6,2]. The 2nd dimension is for two time steps.
!                     The 1st dimension is: [u'_ia,u'_ia-1,v'_ja,v'_ja-1,w'_ka,w'_ka-1]
!
!  Output:
!
!    sp,sn          : crossing time (positive and negative directions) to reach the 
!                     grid box wall (in units of s/m3)
!
!  Other parameters used in function:
!    ba             : linear interpolation of the transport (Eqn 1.6 in user manual)
!    uu             : time-interpolated flux at ia/ja/ka (depending on ijk)
!    um             : time-interpolated flux at ia-1/ja-1/ka-1 (depending on ijk)
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!====================================================================

implicit none

integer,            intent(in)                                      :: ijk,ia,ja,ka,ff,imt,jmt,km
integer,            intent(in)                                      :: do3d, doturb
real(kind=8),       intent(in)                                      :: r0,rr
real(kind=8),       intent(in),     dimension(imt-1,jmt,km,2)       :: uflux
real(kind=8),       intent(in),     dimension(imt,jmt-1,km,2)       :: vflux
real(kind=8),       intent(in),     dimension(0:km,2)               :: wflux
real(kind=8),       intent(out)                                     :: sp,sn
real*8, optional,   intent(out),    dimension(6,2)                  :: upr  
real(kind=8)                                                        :: ba,uu,um,rg
integer                                                             :: ii,im,nsm=1,nsp=2
real(kind=8),       PARAMETER                                       :: UNDEF=1.d20

rg=1.d0-rr

if(ijk.eq.1) then
    ii=ia
    im=ia-1
    ! KMT: This statement occurs many times throughout the code and I don't
    ! understand its purpose if it is supposed to be a hardwired boundary condition.
    ! I am not sure this statement would ever actually be used, unless maybe the
    ! drifter makes it to the very edge of the domain in x, on the west side, and
    ! even then, the uflux is not set up to have values for imt.
    ! Is there another boundary condition that should be taking this place?
!     if(im.eq.0) im=imt
    uu=(rg*uflux(ia,ja,ka,nsp)+rr*uflux(ia,ja,ka,nsm))*ff ! this is interpolation between time fields
    um=(rg*uflux(im,ja,ka,nsp)+rr*uflux(im,ja,ka,nsm))*ff

!     print *,'x: before adding in turb vels'
!     print *,'uu=',uu,' um=',um

    if(doturb==1) then   
        if(r0.ne.dble(ii)) then
            uu=uu+upr(1,2)  
        else
            uu=uu+upr(1,1)  ! add u' from previous iterative time step if on box wall
        endif
        if(r0.ne.dble(im)) then
            um=um+upr(2,2)
        else
            um=um+upr(2,1)  ! add u' from previous iterative time step if on box wall
        endif
    endif

!     print *,'x: after adding in turb vels'
!     print *,'uu=',uu,' um=',um

else if(ijk.eq.2) then
    ii=ja
    uu=(rg*vflux(ia,ja  ,ka,nsp)+rr*vflux(ia,ja  ,ka,nsm))*ff
    um=(rg*vflux(ia,ja-1,ka,nsp)+rr*vflux(ia,ja-1,ka,nsm))*ff

!     print *,'y: before adding in turb vels'
!     print *,'uu=',uu,' um=',um

    if(doturb==1) then   
        if(r0.ne.dble(ja  )) then
            uu=uu+upr(3,2)  
        else
            uu=uu+upr(3,1)  ! add u' from previous iterative time step if on box wall
        endif
        if(r0.ne.dble(ja-1)) then
            um=um+upr(4,2)
        else
            um=um+upr(4,1)  ! add u' from previous iterative time step if on box wall
        endif
    endif

!     print *,'y: after adding in turb vels'
!     print *,'uu=',uu,' um=',um

else if(ijk.eq.3) then
    ii=ka
! #ifdef full_wflux
!       uu=wflux(ia ,ja ,ka   ,nsm)
!       um=wflux(ia ,ja ,ka-1 ,nsm)
! #else
    uu=rg*wflux(ka  ,nsp)+rr*wflux(ka  ,nsm)
    um=rg*wflux(ka-1,nsp)+rr*wflux(ka-1,nsm)
! #endif

    if(do3d==0 .and. doturb==1) then
        if(r0.ne.dble(ka  )) then
            uu=uu+upr(5,2)  
        else
            uu=uu+upr(5,1)  ! add u' from previous iterative time step if on box wall
        endif
        if(r0.ne.dble(ka-1)) then
            uu=uu+upr(6,2)  
        else
            uu=uu+upr(6,1)  ! add u' from previous iterative time step if on box wall
        endif
    endif
endif

! east, north or upward crossing
if(uu.gt.0.d0 .and. r0.ne.dble(ii)) then
    if(um.ne.uu) then
        ba=(r0+dble(-ii+1)) * (uu-um) + um ! linear interpolation of the transport (Eqn 1.6)
        if(ba.gt.0.d0) then
            sp=( dlog(ba) - dlog(uu) )/(um-uu)  ! -or-  sp=-1.d0/(um-uu)*( dlog(uu) - dlog(ba) )
        else
            sp=UNDEF
        endif
    else
        sp=(dble(ii)-r0)/uu
    endif
else
    sp=UNDEF
endif

if(sp.le.0.d0) sp=UNDEF

! west, south or downward crossing
if(um.lt.0.d0 .and. r0.ne.dble(ii-1)) then
    if(um.ne.uu) then
        ba=-((r0-dble(ii))*(uu-um)+uu) 
        if(ba.gt.0.d0) then
            sn=( dlog(ba) - dlog(-um)  )/(um-uu) ! -or- sn=-1.d0/(um-uu)*( dlog(-um) - dlog(ba)  )
        else
            sn=UNDEF
        endif
    else
        sn=(dble(ii-1)-r0)/uu
    endif
else
    sn=UNDEF
endif

if(sn.le.0.d0) sn=UNDEF

return
end subroutine cross