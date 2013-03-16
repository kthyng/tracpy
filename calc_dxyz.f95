subroutine calc_dxyz(ib,jb,kb,rr,KM,kmt,dxyzarray,dxyz,JMT,IMT)

! Inputs:
!  dzt
!  ib
!  jb
!  kb
!  rr
!  KM
!  hs
!  dz
!  kmt
!  dztb
!  dxdy
!
! Outputs:
!  dxyz

implicit none

! #ifdef zgrid3Dt 
!     REAL, DIMENSION(:,:,:,:)       :: dzt
! #elif zgrid3D
!     REAL, DIMENSION(:,:,:)         :: dzt
! #endif /*zgrid3Dt*/
! #ifdef varbottombox 
!     REAL, DIMENSION(:,:,:)         :: dztb
! #endif /*varbottombox*/
integer,      intent(in)                   :: kmt,ib,jb,kb,KM,JMT,IMT
! real(kind=4),   intent(in),     dimension(IMT,JMT,2)            :: hs
real(kind=8),   intent(in),     dimension(IMT-1,JMT-1,KM,2)         :: dxyzarray
! real(kind=8),   intent(in),     dimension(IMT,JMT)              :: dxdy
! real(kind=8),   intent(in),     dimension(KM)                   :: dz
real(kind=8)                                    :: rg
integer                                         :: nsm=1,nsp=2
real(kind=8), intent(in) :: rr
real(kind=8), intent(out)                                    :: dxyz


    rg=1.d0-rr

! T-box volume in m3
! #ifdef zgrid3Dt 
    dxyz=rg*dxyzarray(ib,jb,kb,nsp)+rr*dxyzarray(ib,jb,kb,nsm)
!     dxyz=rg*dzt(ib,jb,kb,nsp)+rr*dzt(ib,jb,kb,nsm)
! #elif  zgrid3D
!     dxyz=dzt(ib,jb,kb)
! #ifdef freesurface
!     if(kb == KM) dxyz=dxyz+rg*hs(ib,jb,nsp)+rr*hs(ib,jb,nsm)
! #endif /*freesurface*/
! #else
!     print *,'kb = ',kb, ' dz = ',dz
!     dxyz=dz(kb)
! #ifdef varbottombox
!     if(kb == KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
! #endif /*varbottombox*/
! #ifdef freesurface
!     if(kb == KM) dxyz=dxyz+rg*hs(ib,jb,nsp)+rr*hs(ib,jb,nsm)
! #endif /*freesurface*/
! #endif /*zgrid3Dt*/
!     dxyz=dxyz*dxdy(ib,jb)
!     if (dxyz<0) then
!        print *,'====================================='
!        print *,'ERROR: Negative box volume           '
!        print *,'-------------------------------------'
!        print *,'dzt  = ', dxyz/dxdy(ib,jb),dz(kb),hs(ib,jb,:)
!        print *,'dxdy = ', dxdy(ib,jb)
!        print *,'ib  = ', ib, ' jb  = ', jb, ' kb  = ', kb 
!        print *,'-------------------------------------'
!        print *,'The run is terminated'
!        print *,'====================================='
! !        errCode = -60
!        stop
!     end if
end subroutine calc_dxyz
