subroutine vertvel(rr,ia,ja,ka,imt,jmt,km,ff,uflux,vflux,wflux)
!====================================================================
! Calculate the vertical flux based on the uflux and vflux
!
!  Input:
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    ia,ja,ka       : original position in integers
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!
!  Output:
!    wflux          : w velocity (vertical) flux field, two time steps [kxt]
!
!  Other parameters used in function:
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    uu             : time-interpolated flux at ia/ja/ka (depending on ijk)
!    um             : time-interpolated flux at ia-1/ja-1/ka-1 (depending on ijk)
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!    im             : generic index for grid index -1 for whichever direction, ijk. 
!                     Is only used in the i direction for whatever reason.
!    k              : Index for looping through depth layers
!    n              : Index for looping through times (1 and 2)
!
!  Notes:
!    Computes the vertical velocity by integrating the continuity eq. from the bottom 
!    for the nsm and nsp velocity time steps.
!    This is set up to use neither the -full_wflux nor -explicit_w flags currently, 
!    just the default w flux option.
!====================================================================

#ifndef explicit_w
  
implicit none
  
integer,        intent(in)                                      :: ff,ia,ja,ka,imt,jmt,km
real(kind=8),   intent(in)                                      :: rr
real(kind=8),   intent(in),     dimension(imt-1,jmt,km,2)       :: uflux
real(kind=8),   intent(in),     dimension(imt,jmt-1,km,2)       :: vflux
real(kind=8),   intent(out),    dimension(0:km,2)               :: wflux
real(kind=8)                                                    :: uu,um,rg
integer                                                         :: nsm=1,nsp=2,k,n,im,iam
   

rg=1.d0-rr
wflux=0.d0
im=ia-1 

#ifdef twodim
    return
  
! start 3D code
#else
    kloop: do k=1,ka
  ! these only need to be defined if we use full_wflux, which we aren't doing right now
!      uu=rg*uflux(ia ,ja  ,k,nsp)+rr*uflux(ia ,ja  ,k,nsm)
!      um=rg*uflux(iam,ja  ,k,nsp)+rr*uflux(iam,ja  ,k,nsm)
!      vv=rg*vflux(ia ,ja  ,k,nsp)+rr*vflux(ia ,ja  ,k,nsm)
!      vm=rg*vflux(ia ,ja-1,k,nsp)+rr*vflux(ia ,ja-1,k,nsm)

! ! start ifs code
! #if defined ifs
!     do n=nsm,nsp
!      wflux(k,n) = wflux(k-1,n) - ff * &
!      ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n)  &
!      + (dzt(ia,ja,k,nsp)-dzt(ia,ja,k,nsm))*dxdy(ia,ja)/tseas )  ! time change of the mass the in grid box
!     enddo
! #endif
! ! end ifs code

! start ocean code
! #ifdef  full_wflux
!      wflux(ia,ja,k,nsm)=wflux(ia,ja,k-1,nsm) - ff * ( uu - um + vv - vm )
! #else
        do n=nsm,nsp
            wflux(k,n) = wflux(k-1,n) - ff * &
                        ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + &
                        vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n) )
        enddo
! #endif
!end ocean code
    end do kloop

#endif
! end 3D code
  
!#endif
  return
#endif
end subroutine vertvel

 
