subroutine calc_dxyz(ib,jb,kb,rr,imt,jmt,km,kmt,dzt,dxdy,dxyz)
!====================================================================
! calculate the grid cell volume at the drifter's location
!
!  Input:
!    ib,jb,kb       : new position in grid space indices
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    kmt            : Number of vertical levels in horizontal space [imt,jmt]
!    dzt 			: Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
!    dxdy 			: Area of horizontal cell-walls [imt,jmt]
!
!  Output:
!    dxyz 			: Volume of drifter's grid cell
!
!  Other parameters used in function:
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    nsm=1,nsp=2    : Time index. nsm picks out the earlier bounding time step and 
!                     nsp picks out the later bounding time step for interpolation.
!
!  Notes:
!    KMT: Assuming zgrid3Dt flag which includes free surface already, I think. Other 
!         flags aren't set up in the allocations or in the input variables,
! 	      for example, hs is not currently input since it is already included in 
!         the dzt.
!====================================================================

implicit none

integer,   	    intent(in)                   					:: ib,jb,kb,km,jmt,imt
real(kind=8),   intent(in),     dimension(imt,jmt,km,2)         :: dzt
real(kind=8),   intent(in),     dimension(imt,jmt)              :: dxdy,kmt
real(kind=8), 	intent(in) 										:: rr
real(kind=8), 	intent(out)                                    	:: dxyz
real(kind=8)                                    				:: rg
integer                                         				:: nsm=1,nsp=2


rg=1.d0-rr

! T-box volume in m3
dxyz=rg*dzt(ib,jb,kb,nsp)+rr*dzt(ib,jb,kb,nsm)
dxyz=dxyz*dxdy(ib,jb)

end subroutine calc_dxyz
