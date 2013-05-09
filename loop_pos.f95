subroutine pos(ia,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1,ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dsc,&
                &ff,imt,jmt,km,rr,rb,uflux,vflux,wflux,upr)
!====================================================================
! calculate the new positions of the trajectory with pos_orgn()   
!
!  Input:
!
!    ia,ja,ka       : original position in grid space indices
!    x0,y0,z0       : original non-dimensional position in the i,j,k-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!    ds             : crossing time to reach the grid box wall (units=s/m3)
!    dse,dsw        : crossing times for drifter to reach the east and west grid box wall
!    dsn,dss        : crossing times for drifter to reach the north and south grid box wall
!    dsu,dsd        : crossing times for drifter to reach the up and down grid box wall
!    dsmin          : time step based on the interpolation step between model output times.
!                   : sets a limit on the time step that a drifter can go.
!    dsc            : Not sure what this is right now
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation (for original time 
!                     at the beginning of the loop?).
!    rb             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation (for next time 
!                     at the end of the loop?).
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!    wflux          : w velocity (vertical) flux field, two time steps [kxt]
!    upr            : parameterized turbulent velocities u', v', and w'
!                     optional because only used if using turb flag for diffusion
!                     size [6,2]. The 2nd dimension is for two time steps.
!                     The 1st dimension is: [u'_ia,u'_ia-1,v'_ja,v'_ja-1,w'_ka,w'_ka-1]
!
!  Output:
!    
!    ib,jb,kb       : new position in grid space indices
!    x1,y1,z1       : updated non-dimensional position in the i,j,k-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!
!  Other parameters used in function:
!    rbg            : rbg=1-rg for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    uu             : time-interpolated flux at ia/ja/ka (depending on ijk)
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!    im             : generic index for grid index -1 for whichever direction, ijk. 
!                     Is only used in the i direction for whatever reason.
!====================================================================
    
implicit none

integer,            intent(in)                                  :: ia, ja, ka, imt, jmt, km,ff
real*8,             intent(in)                                  :: x0, y0, z0,ds,dse,dsw,dss,dsn,dsd,dsu,dsmin,dsc,rb, rr
real(kind=8),       intent(in),     dimension(imt-1,jmt,km,2)   :: uflux
real(kind=8),       intent(in),     dimension(imt,jmt-1,km,2)   :: vflux
real(kind=8),       intent(in),     dimension(0:km,2)           :: wflux
real*8, optional,   intent(in),     dimension(6,2)              :: upr  
integer,            intent(out)                                 :: ib, jb, kb
real*8,             intent(out)                                 :: x1, y1, z1
integer                                                         :: nsm=1,nsp=2,im,iam
real(kind=8)                                                    :: uu, rbg
real(kind=8),       parameter                                   :: UNDEF=1.d20

im=ia-1
rbg =1.d0-rb

! === calculate the new positions ===
! === of the trajectory           ===    
!     scrivi=.false.
if(ds==dse) then ! eastward grid-cell exit 
!        scrivi=.false. ! flag for when to write to file I think
    uu=(rbg*uflux(ia,ja,ka,nsp)+rb*uflux(ia ,ja,ka,nsm))*ff
    ! if the drifter is exiting east and the east transport is positive,
    ! bump the east index up by one to keep it greater than x
    ! and change the x1 value to be the value at the west side of the new grid cell
    if(uu.gt.0.d0) then
        ib=ia+1
        if(ib.gt.imt) ib=ib-imt ! imt is a grid parameter
    endif
    x1=dble(ia)
#ifdef turb
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) 
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr)
#else
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) 
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km)
#endif /*turb*/

else if(ds==dsw) then ! westward grid-cell exit
!        scrivi=.false.
    uu=(rbg*uflux(iam,ja,ka,nsp)+rb*uflux(iam,ja,ka,nsm))*ff
    if(uu.lt.0.d0) then
        ib=iam
    endif
    x1=dble(iam)
#ifdef turb
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! meridional position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! vertical position
#else
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! meridional position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! vertical position
#endif /*turb*/
!       scrivi=.true.      

else if(ds==dsn) then ! northward grid-cell exit
!        scrivi=.false.
    uu=(rbg*vflux(ia,ja,ka,nsp)+rb*vflux(ia,ja,ka,nsm))*ff
    if(uu.gt.0.d0) then
        jb=ja+1
    endif
    y1=dble(ja)
#ifdef turb
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! zonal position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! vertical position
#else
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! zonal position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! vertical position
#endif /*turb*/

else if(ds==dss) then ! southward grid-cell exit       
!        scrivi=.false.
    uu=(rbg*vflux(ia,ja-1,ka,nsp)+rb*vflux(ia,ja-1,ka,nsm))*ff
    if(uu.lt.0.d0) then
        jb=ja-1
! #ifndef ifs 
!     if(jb==0) stop 34578
! #endif
    endif
    y1=dble(ja-1)
#ifdef turb
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! zonal position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! vertical position
#else
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! zonal position
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! vertical position
#endif /*turb*/
       
elseif(ds==dsu) then ! upward grid-cell exit
!        scrivi=.false.
    call vertvel(rr,ia,iam,ja,ka,imt,jmt,km,ff,uflux,vflux,wflux)
! #ifdef full_wflux
!        uu=wflux(ia,ja,ka,nsm)
! #else
    uu=rbg*wflux(ka,nsp)+rb*wflux(ka,nsm)
! #endif
    if(uu.gt.0.d0) then
        kb=ka+1
    endif
    z1=dble(ka)
    if(kb==km+1) then    ! prevent "evaporation" and put particle from the surface
        kb=km           
        z1=dble(km)-0.5d0 ! to the middle of the surface layer
    endif
#ifdef turb
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! zonal position
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! meridional position
#else
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! zonal position
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! meridional position
#endif /*turb*/

elseif(ds==dsd) then ! downward grid-cell exit
!        scrivi=.false.
    call vertvel(rr,ia,iam,ja,ka,imt,jmt,km,ff,uflux,vflux,wflux)
       
! #ifdef full_wflux
!        if(wflux(ia,ja,ka-1,nsm).lt.0.d0) kb=ka-1
! #else
    if(rbg*wflux(ka-1,nsp)+rb*wflux(ka-1,nsm).lt.0.d0) kb=ka-1
! #endif              
    z1=dble(ka-1)
#ifdef turb
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! zonal position
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! meridional position
#else
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! zonal position
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! meridional position
#endif /*turb*/

else if( ds==dsc .or. ds==dsmin) then  
       ! shortest time is the time-steping 
!        scrivi=.true.
       ! If there is no spatial solution, 
       ! which should correspond to a convergence zone
        if(dse==UNDEF .and. dsw==UNDEF .and. dsn==UNDEF .and. & 
          dss==UNDEF .and. dsu==UNDEF .and. dsd==UNDEF ) then
          
            ! move if atmosphere, freeze if ocean
            ib=ia ; jb=ja ; kb=ka

!#ifdef ifs
!          call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr) ! zonal crossing 
!          call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr) ! merid. crossing 
!          call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vert. crossing 
!#else
!          x1=x0 ; y1=y0 ; z1=z0 
!          print *,ib,jb,kb,x1,y1,z1
!#endif  
        ! If there is at least one spatial solution 
        ! but the shortest cross time is the time step
        endif
!       else

#ifdef turb
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! zonal crossing 
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! merid. crossing 
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km,upr) ! vert. crossing 
#else
    call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! zonal crossing 
    call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! merid. crossing 
    call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,wflux,ff,imt,jmt,km) ! vert. crossing 
#endif /*turb*/
!       endif
endif
    

end subroutine pos
