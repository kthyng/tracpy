subroutine pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1,ds,dse,dsw,dsn,dss,dsmin,dsc,ff,IMT,JMT,KM,rr,rbg,rb,uflux,vflux)
! === calculate the new positions ===
! === of the trajectory           ===    
!
! Inputs:
!
! Outputs:
!  
    
    IMPLICIT none
    INTEGER                                    :: nsm=1,nsp=2
    REAL(kind=8)                                       :: uu
    INTEGER, intent(in)                        :: ia, iam, ja, ka, IMT, JMT, KM,ff
    integer, intent(out)  :: ib, jb, kb
    REAL*8, INTENT(IN)                         :: x0, y0, z0,ds,dse,dsw,dss,dsn,dsmin,dsc,rbg,rb, rr
    REAL*8, INTENT(OUT)                        :: x1, y1, z1
    REAL(kind=8), PARAMETER                         :: UNDEF=1.d20
real(kind=8),   intent(in),     dimension(IMT,JMT-1,KM,2)         :: uflux
real(kind=8),   intent(in),     dimension(IMT-1,JMT,KM,2)         :: vflux

    ! === calculate the new positions ===
    ! === of the trajectory           ===    
!     scrivi=.false.
    if(ds==dse) then ! eastward grid-cell exit 
!         print *,'dse'
!        scrivi=.false. ! flag for when to write to file I think
        uu=(rbg*uflux(ia,ja,ka,nsp)+rb*uflux(ia ,ja,ka,nsm))*ff
        ! if the drifter is exiting east and the east transport is positive,
        ! bump the east index up by one to keep it greater than x
        ! and change the x1 value to be the value at the west side of the new grid cell
        if(uu.gt.0.d0) then
            ib=ia+1
            if(ib.gt.IMT) ib=ib-IMT ! IMT is a grid parameter
        endif
        x1=dble(ia)
        call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) 
        call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,ff,IMT,JMT,KM)

    else if(ds==dsw) then ! westward grid-cell exit
!         print *,'dsw'
!        scrivi=.false.
        uu=(rbg*uflux(iam,ja,ka,nsp)+rb*uflux(iam,ja,ka,nsm))*ff
        if(uu.lt.0.d0) then
            ib=iam
        endif
        x1=dble(iam)
        call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! meridional position
        call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! vertical position
!       scrivi=.true.      

    else if(ds==dsn) then ! northward grid-cell exit
!         print *,'dsn'
!        scrivi=.false.
        uu=(rbg*vflux(ia,ja,ka,nsp)+rb*vflux(ia,ja,ka,nsm))*ff
        if(uu.gt.0.d0) then
            jb=ja+1
        endif
        y1=dble(ja)
        call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! zonal position
        call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! vertical position

    else if(ds==dss) then ! southward grid-cell exit
!         print *,'dss'
       
!        scrivi=.false.
        uu=(rbg*vflux(ia,ja-1,ka,nsp)+rb*vflux(ia,ja-1,ka,nsm))*ff
        if(uu.lt.0.d0) then
            jb=ja-1
#ifndef ifs 
        if(jb==0) stop 34578
#endif
       endif
       y1=dble(ja-1)
       call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! zonal position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! vertical position
       
!     elseif(ds==dsu) then ! upward grid-cell exit
!        scrivi=.false.
!        call vertvel(rb,ia,iam,ja,ka)
! #ifdef full_wflux
!        uu=wflux(ia,ja,ka,nsm)
! #else
!        uu=rbg*wflux(ka,nsp)+rb*wflux(ka,nsm)
! #endif
!        if(uu.gt.0.d0) then
!           kb=ka+1
!        endif
!        z1=dble(ka)
!        if(kb==KM+1) then    ! prevent "evaporation" and put particle from the surface
!           kb=KM           
!           z1=dble(KM)-0.5d0 ! to the middle of the surface layer
!        endif
!        call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr)
!        call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr)


!     elseif(ds==dsd) then ! downward grid-cell exit
!        scrivi=.false.
!        call vertvel(rb,ia,iam,ja,ka)
       
! #ifdef full_wflux
!        if(wflux(ia,ja,ka-1,nsm).lt.0.d0) kb=ka-1
! #else
!        if(rbg*wflux(ka-1,nsp)+rb*wflux(ka-1,nsm).lt.0.d0) kb=ka-1
! #endif              
!        z1=dble(ka-1)
!        call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr)
!        call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr)

    else if( ds==dsc .or. ds==dsmin) then  
!         print *,'ds=',ds,' dsc=',dsc,' dsmin=',dsmin
!         print *,'other'
!         print *,'dsc=',dsc,' dsmin=',dsmin,' ds=',ds
!             print *, 'in other in pos: ib=',ib,' ia=',ia

       ! shortest time is the time-steping 
!        scrivi=.true.
       ! If there is no spatial solution, 
       ! which should correspond to a convergence zone
!        if(dse==UNDEF .and. dsw==UNDEF .and. dsn==UNDEF .and. & 
!           dss==UNDEF .and. dsu==UNDEF .and. dsd==UNDEF ) then
        if(dse==UNDEF .and. dsw==UNDEF .and. dsn==UNDEF .and. & 
          dss==UNDEF  ) then
          
          ! move if atmosphere, freeze if ocean
            ib=ia ; jb=ja ; kb=ka
!             print *, 'in other in pos: ib=',ib,' ia=',ia
!          print *,'convergence for ',ib,jb,kb,x0,y0,z0
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
! print *,'before: ia=',ia,' x0=',x0,' x1=',x1
        call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! zonal crossing 
        call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! merid. crossing 
        call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr,uflux,vflux,ff,IMT,JMT,KM) ! vert. crossing 
! print *,'after: ia=',ia,' x0=',x0,' x1=',x1
!       endif
    endif
! print *,'x1=',x1,' y1=',y1
    

end subroutine pos
