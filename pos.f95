subroutine pos_orgn(ijk,ia,ja,ka,r0,r1,ds,rr,uflux,vflux,ff,IMT,JMT,KM)
#ifndef timeanalyt 
!====================================================================
! computes new position (r0 --> r1) of trajectory after time ds
! the new coordinate is still on one of the faces of box at ia,ja,ka
!
!  Input:
!
!    ijk      : considered direction (i=zonal,2=meridional,3=vertical)
!    ia,ja,ka : original position in integers
!    r0       : original non-dimensional position in the ijk-direction
!                   of particle (fractions of a grid box side in the 
!                                corresponding direction)
!    rr       : time interpolation constant between 0 and 1 
!    sp       : crossing time to reach the grid box wall (units=s/m3)
!
!  Output:
!    
!    r1       : the new position (coordinate)
!====================================================================

IMPLICIT none

real*8 :: rg,ds,uu,um,vv,vm,en
integer :: ii,im,nsm=1,nsp=2
real(kind=8), intent(in) :: r0,rr
real(kind=8), intent(out) :: r1
integer, intent(in) :: ijk,ia,ja,ka,ff,IMT,JMT,KM
real(kind=8),   intent(in),     dimension(IMT,JMT-1,KM,2)         :: uflux
real(kind=8),   intent(in),     dimension(IMT-1,JMT,KM,2)         :: vflux

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
    if(im.eq.0) im=IMT
    uu=(rg*uflux(ia,ja,ka,nsp)+rr*uflux(ia,ja,ka,nsm))*ff
    um=(rg*uflux(im,ja,ka,nsp)+rr*uflux(im,ja,ka,nsm))*ff
!     print *,'ijk=',ijk,' rg=',rg,' rr=',rr,' ff=',ff
!     print *,'uflux(ia,ja,ka,nsp)=',uflux(ia,ja,ka,nsp),' uflux(ia,ja,ka,nsm)=',uflux(ia,ja,ka,nsm)
!     print *,'uflux(ia,ja-1,ka,nsp)=',uflux(ia,ja-1,ka,nsp),' uflux(ia,ja-1,ka,nsm)=',uflux(ia,ja-1,ka,nsm)
! #ifdef turb    
!      if(r0.ne.dble(ii)) then
!         uu=uu+upr(1,2)  
!      else
!         uu=uu+upr(1,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
!      if(r0.ne.dble(im)) then
!         um=um+upr(2,2)
!      else
!         um=um+upr(2,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
! #endif
else if(ijk.eq.2) then
    ii=ja
    uu=(rg*vflux(ia,ja  ,ka,nsp)+rr*vflux(ia,ja  ,ka,nsm))*ff
    um=(rg*vflux(ia,ja-1,ka,nsp)+rr*vflux(ia,ja-1,ka,nsm))*ff
!     print *,'ijk=',ijk,' rg=',rg,' rr=',rr,' ff=',ff
!     print *,'vflux(ia,ja,ka,nsp)=',vflux(ia,ja,ka,nsp),' vflux(ia,ja,ka,nsm)=',vflux(ia,ja,ka,nsm)
!     print *,'vflux(ia,ja-1,ka,nsp)=',vflux(ia,ja-1,ka,nsp),' vflux(ia,ja-1,ka,nsm)=',vflux(ia,ja-1,ka,nsm)
! #ifdef turb    
!      if(r0.ne.dble(ja  )) then
!         uu=uu+upr(3,2)  
!      else
!         uu=uu+upr(3,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
!      if(r0.ne.dble(ja-1)) then
!         um=um+upr(4,2)
!      else
!         um=um+upr(4,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
! #endif
!   elseif(ijk.eq.3) then
!      ii=ka
! #ifdef full_wflux
!      uu=wflux(ia ,ja ,ka   ,nsm)
!      um=wflux(ia ,ja ,ka-1 ,nsm)
! #else
!      uu=rg*wflux(ka  ,nsp)+rr*wflux(ka  ,nsm)
!      um=rg*wflux(ka-1,nsp)+rr*wflux(ka-1,nsm)
! #endif
! #ifndef twodim   
! #ifdef turb    
!      if(r0.ne.dble(ka  )) then
!         uu=uu+upr(5,2)  
!      else
!         uu=uu+upr(5,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
!      if(r0.ne.dble(ka-1)) then
!         uu=uu+upr(6,2)  
!      else
!         uu=uu+upr(6,1)  
!         ! add u' from previous iterative time step if on box wall
!      endif
! #endif
! #endif
endif

!
! note: consider in future to improve the code below for accuracy 
! in case of um-uu = small; also see subroutine cross
if(um.ne.uu) then
    r1= (r0+(-dble(ii-1) + um/(uu-um))) * dexp( (uu-um)*ds ) + dble(ii-1) - um/(uu-um)
! print *,'r0=',r0,' ii=',ii,' um=',um,' uu=',uu,' ds=',ds,' r1=',r1,' dble(ii-1)=',dble(ii-1)
else
    r1=r0+uu*ds
endif
! print *,'ijk=',ijk,' r0=',r0,' ii=',ii,' um=',um,' uu=',uu,' ds=',ds,' r1=',r1
! print *,'rg=',rg,' uflux=(ia,ja,ka,nsp)=',uflux(ia,ja,ka,nsp),' rr=',rr,' uflux(ia,ja,ka,nsm)=',uflux(ia,ja,ka,nsm)
! print *,'ia=',ia,' ja=',ja,' ka=',ka
! print *,' uflux=(im,ja,ka,nsp)=',uflux(im,ja,ka,nsp),' uflux(im,ja,ka,nsm)=',uflux(ia,ja,ka,nsm),' ff=',ff
!if(abs(um/(uu-um)).gt.1.d10) print *,'possible precision problem?',um/(uu-um),uu,um,ijk,ia,ja,ka,r0,r1,ds,rr

return
#endif
end subroutine pos_orgn
