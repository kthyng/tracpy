subroutine pos_orgn(ijk,ia,ja,ka,r0,r1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr)
#ifndef timeanalyt 
!====================================================================
! computes new position (r0 --> r1) of trajectory after time ds
! the new coordinate is still on one of the faces of box at ia,ja,ka
!
!  Input:
!
!    ijk            : considered direction (1=zonal,2=meridional,3=vertical)
!    ia,ja,ka       : original position in integers
!    r0             : original non-dimensional position in the ijk-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!    ds             : crossing time to reach the grid box wall (units=s/m3)
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!    wflux          : w velocity (vertical) flux field, two time steps [kxt]
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    upr            : parameterized turbulent velocities u', v', and w'
!                     optional because only used if using turb flag for diffusion
!                     size [6,2]. The 2nd dimension is for two time steps.
!                     The 1st dimension is: [u'_ia,u'_ia-1,v'_ja,v'_ja-1,w'_ka,w'_ka-1]
!
!  Output:
!    
!    r1             : the new position (coordinate)
!
!  Other parameters used in function:
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    uu             : time-interpolated flux at ia/ja/ka (depending on ijk)
!    um             : time-interpolated flux at ia-1/ja-1/ka-1 (depending on ijk)
!    ii             : generic index for grid index for whichever direction, ijk
!    im             : generic index for grid index -1 for whichever direction, ijk. 
!                     Is only used in the i direction for whatever reason.
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!====================================================================

implicit none

integer,            intent(in)                                  :: ijk,ia,ja,ka,ff,imt,jmt,km
real(kind=8),       intent(in)                                  :: r0,rr,ds
real(kind=8),       intent(in),     dimension(imt-1,jmt,km,2)   :: uflux
real(kind=8),       intent(in),     dimension(imt,jmt-1,km,2)   :: vflux
real(kind=8),       intent(in),     dimension(0:km,2)           :: wflux
real*8, optional,   intent(in),     dimension(6,2)              :: upr  
real(kind=8),       intent(out)                                 :: r1
real*8                                                          :: rg,uu,um
integer                                                         :: ii,im,nsm=1,nsp=2

rg=1.d0-rr

#ifdef twodim   
    if(ijk.eq.3) then
        r1=r0
        return
    endif
#endif


if(ijk.eq.1) then
    ii=ia
    im=ia-1
    if(im.eq.0) im=imt
    uu=(rg*uflux(ia,ja,ka,nsp)+rr*uflux(ia,ja,ka,nsm))*ff
    um=(rg*uflux(im,ja,ka,nsp)+rr*uflux(im,ja,ka,nsm))*ff

#ifdef turb    
    if(r0.ne.dble(ii)) then
        uu=uu+upr(1,2)  
    else
        uu=uu+upr(1,1)  
        ! add u' from previous iterative time step if on box wall
    endif
    if(r0.ne.dble(im)) then
        um=um+upr(2,2)
    else
        um=um+upr(2,1)  
        ! add u' from previous iterative time step if on box wall
    endif
#endif

else if(ijk.eq.2) then
    ii=ja
    uu=(rg*vflux(ia,ja  ,ka,nsp)+rr*vflux(ia,ja  ,ka,nsm))*ff
    um=(rg*vflux(ia,ja-1,ka,nsp)+rr*vflux(ia,ja-1,ka,nsm))*ff

#ifdef turb    
    if(r0.ne.dble(ja  )) then
        uu=uu+upr(3,2)  
    else
        uu=uu+upr(3,1)  
        ! add u' from previous iterative time step if on box wall
    endif
    if(r0.ne.dble(ja-1)) then
        um=um+upr(4,2)
    else
        um=um+upr(4,1)  
        ! add u' from previous iterative time step if on box wall
    endif
#endif

elseif(ijk.eq.3) then
    ii=ka
! #ifdef full_wflux
!      uu=wflux(ia ,ja ,ka   ,nsm)
!      um=wflux(ia ,ja ,ka-1 ,nsm)
! #else
    uu=rg*wflux(ka  ,nsp)+rr*wflux(ka  ,nsm)
    um=rg*wflux(ka-1,nsp)+rr*wflux(ka-1,nsm)
! #endif

#ifndef twodim   
#ifdef turb    
    if(r0.ne.dble(ka  )) then
        uu=uu+upr(5,2)  
    else
        uu=uu+upr(5,1)  
        ! add u' from previous iterative time step if on box wall
    endif
    if(r0.ne.dble(ka-1)) then
        uu=uu+upr(6,2)  
    else
        uu=uu+upr(6,1)  
        ! add u' from previous iterative time step if on box wall
    endif
#endif
#endif
endif

!
! note: consider in future to improve the code below for accuracy 
! in case of um-uu = small; also see subroutine cross
if(um.ne.uu) then
    r1= (r0+(-dble(ii-1) + um/(uu-um))) * dexp( (uu-um)*ds ) + dble(ii-1) - um/(uu-um)
else
    r1=r0+uu*ds
endif

! print '(a,i1,a,f6.2,a,f6.2,a,e6.1,a,f4.2)','ijk=',ijk,' r0=',r0,' r1=',r1,' ds=',ds,' rr=',rr
! print '(a,f8.1,a,f8.1)','uu=',uu,' um=',um
! print '(a,f8.1,a,f8.1,a)','vflux(ia,ja,ka,1)=',vflux(ia,ja,ka,1),' vflux(ia ,ja,ka,2)=',vflux(ia ,ja,ka,2)
! print '(a,f7.1,a,f6.1,a,f5.2,a,f5.2)', 'tt=',tt,' dt=',dt,' ts=',ts,' tss=',tss

!if(abs(um/(uu-um)).gt.1.d10) print *,'possible precision problem?',um/(uu-um),uu,um,ijk,ia,ja,ka,r0,r1,ds,rr

return
#endif
end subroutine pos_orgn
