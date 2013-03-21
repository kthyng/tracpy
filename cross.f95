subroutine cross(ijk,ia,ja,ka,r0,sp,sn,rr,uflux,vflux,ff,KM,JMT,IMT)
  
! subroutine to compute time (sp,sn) when trajectory 
! crosses face of box (ia,ja,ka) 
! two crossings are considered for each direction:  
! east and west for longitudinal directions, etc.  
! 
!  Input:
!
!  ijk      : considered direction (i=zonal, 2=meridional, 3=vertical)
!  ia,ja,ka : original position in integers
!  r0       : original non-dimensional position in the 
!             ijk-direction of particle 
!             (fractions of a grid box side in the 
!              corresponding direction)
!  rr       : time interpolation constant between 0 and 1 
!
!
!  Output:
!
!    sp,sn   : crossing time to reach the grid box wall 
!              (in units of s/m3)
  

!     USE mod_param
!     USE mod_vel
!     USE mod_turb
!     USE mod_time
    
IMPLICIT none
integer, intent(in)     :: ijk,ia,ja,ka,ff,IMT,JMT,KM
real(kind=8),   intent(in),     dimension(IMT,JMT-1,KM,2)         :: uflux
real(kind=8),   intent(in),     dimension(IMT-1,JMT,KM,2)         :: vflux
real(kind=8),   intent(in)  :: r0,rr
real(kind=8) :: ba,uu,um,rg,vv,vm
integer  :: ii,im,nsm=1,nsp=2
real(kind=8), PARAMETER                         :: UNDEF=1.d20
real(kind=8),   intent(out)  :: sp,sn

    rg=1.d0-rr

    if(ijk.eq.1) then
        ii=ia
        im=ia-1
        if(im.eq.0) im=IMT
        uu=(rg*uflux(ia,ja,ka,nsp)+rr*uflux(ia,ja,ka,nsm))*ff ! this is interpolation between time fields
        um=(rg*uflux(im,ja,ka,nsp)+rr*uflux(im,ja,ka,nsm))*ff
!         print *,'in cross: uu=',uu,' um=',um
! #ifdef turb   
!       if(r0.ne.dble(ii)) then
!         uu=uu+upr(1,2)  
!       else
!         uu=uu+upr(1,1)  ! add u' from previous iterative time step if on box wall
!       endif
!       if(r0.ne.dble(im)) then
!         um=um+upr(2,2)
!       else
!         um=um+upr(2,1)  ! add u' from previous iterative time step if on box wall
!       endif
! #endif
    else if(ijk.eq.2) then
        ii=ja
        uu=(rg*vflux(ia,ja  ,ka,nsp)+rr*vflux(ia,ja  ,ka,nsm))*ff
        um=(rg*vflux(ia,ja-1,ka,nsp)+rr*vflux(ia,ja-1,ka,nsm))*ff
!         print *,'in cross: vu=',uu,' vm=',um
! #ifdef turb    
!       if(r0.ne.dble(ja  )) then
!         uu=uu+upr(3,2)  
!       else
!         uu=uu+upr(3,1)  ! add u' from previous iterative time step if on box wall
!       endif
!       if(r0.ne.dble(ja-1)) then
!         um=um+upr(4,2)
!       else
!         um=um+upr(4,1)  ! add u' from previous iterative time step if on box wall
!       endif
! #endif
!     elseif(ijk.eq.3) then
!       ii=ka
! #ifdef full_wflux
!       uu=wflux(ia ,ja ,ka   ,nsm)
!       um=wflux(ia ,ja ,ka-1 ,nsm)
! #else
!       uu=rg*wflux(ka  ,nsp)+rr*wflux(ka  ,nsm)
!       um=rg*wflux(ka-1,nsp)+rr*wflux(ka-1,nsm)
! #endif

! #ifndef twodim && turb   
!       if(r0.ne.dble(ka  )) then
!         uu=uu+upr(5,2)  
!       else
!         uu=uu+upr(5,1)  ! add u' from previous iterative time step if on box wall
!       endif
!       if(r0.ne.dble(ka-1)) then
!         uu=uu+upr(6,2)  
!       else
!         uu=uu+upr(6,1)  ! add u' from previous iterative time step if on box wall
!       endif
! #endif
    endif

! east, north or upward crossing
    if(uu.gt.0.d0 .and. r0.ne.dble(ii)) then
!         print *,'uu=',uu,' r0=',r0,' ii=',ii,' um=',um
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
!               print *,'sn=',sn
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

!      print *,'in cross: sp=',sp,' sn=',sn,' r0=',r0,' ii=',ii,' ba=',ba
!     print *,'in cross: ijk=',ijk,' sp=',sp,' sn=',sn,' uu=',uu,' um=',um,' ba=',ba,' r0=',r0

    return
end subroutine cross