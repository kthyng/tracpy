
subroutine vertvel(rr,ia,iam,ja,ka,imt,jmt,km,ff,uflux,vflux,wflux)

! Inputs:
!  rr
!  ia
! etc
!
! Outputs:
!  wflux    Vertical flux

#ifndef explicit_w
  
  ! === Computes the vertical velocity by integrating ===
  ! === the continuity eq. from the bottom            ===
  ! === for the nsm and nsp velocity time steps       ===
    
  implicit none
  
  
integer k,n
integer, intent(in)     :: ff,ia,iam,ja,ka,imt,jmt,km
real(kind=8) :: uu,um,rg,vv,vm
real(kind=8),   intent(in)  :: rr
real(kind=8),   intent(in),     dimension(imt-1,jmt,km,2)         :: uflux
real(kind=8),   intent(in),     dimension(imt,jmt-1,km,2)         :: vflux
! #ifdef  full_wflux
!     ! Not sure if this is right
!     real(kind=8),   intent(out),     dimension(IMT,JMT,KM+1,2)         :: wflux
! #else
    real(kind=8),   intent(out),     dimension(0:km,2)         :: wflux
! #endif
integer                                         :: nsm=1,nsp=2
   
rg=1.d0-rr
wflux=0.d0
  
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
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n) )
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

 