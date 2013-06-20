subroutine diffuse(x1, y1, z1, ib, jb, kb, dt,imt,jmt,km,kmt,dxv,dyu,dzt,h,ah,av,do3d,doturb)

!============================================================================
! Add a small displacement to a particle s.t. it is still in the model area.
! 
!  Input:
!    dt             : time step of trajectory iteration in seconds  
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    kmt            : Number of vertical levels in horizontal space [imt,jmt]
!    dxv            : Horizontal grid cell walls areas in x direction [imt,jmt-1]
!    dyu            : Horizontal grid cell walls areas in y direction [imt-1,jmt]
!    dzt            : Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
!    h              : Depths [imt,jmt]
!    ah             : horizontal diffusion in m^2/s. Input parameter.
!                   : For -turb and -diffusion flags
!    av             : vertical diffusion in m^2/s. Input parameter. For -diffusion flag.
!    do3d           : Flag to set whether to use 3d velocities or not
!    doturb         : Flag to set whether or not to use turb/diff and which kind if so
!
!  Input/Output:
!    x1, y1, z1     : Current positions of the particle. Are updated by the subroutine.
!    ib, jb, kb     : Indices for current box for the particle. Are updated by the subroutine.
!
!  Other parameters used in function:
!    tmpX, tmpY, tmpZ   : Temporal position
!    xd, yd, zd         : Displacement
!    itno               : Number of iterations
!    tmpi, tmpj, tmpk   : Temporal box indices
!    tryAgain           : Tells whether to continue displace
!============================================================================
  
implicit none

integer,        intent(in),     dimension(imt,jmt)          :: kmt
integer,        intent(in)                                  :: do3d, doturb
real(kind=8),   intent(in)                                  :: ah,av,dt
real*8,         intent(in),     dimension(imt-1,jmt)        :: dyu
real*8,         intent(in),     dimension(imt,jmt-1)        :: dxv
real*8,         intent(in),     dimension(imt,jmt)          :: h
real(kind=8),   intent(in),     dimension(imt,jmt,km,2)     :: dzt
integer,        intent(in)                                  :: imt,jmt,km           
real(kind=8),   intent(in out)                              :: x1, y1, z1
integer,        intent(in out)                              :: ib,jb,kb
integer                                                     :: itno,tmpi, tmpj, tmpk
real(kind=8)                                                :: xd, yd, zd, tmpX, tmpY, tmpZ
logical                                                     :: tryAgain       

! print *,''
! print *,'start of diffusion'
! print *,'x1=',x1,' y1=',y1,' ib=',ib,' jb=',jb

tryAgain = .FALSE.
  
! Is particle within model area?
! KMT: I do not understand how the vertical checks are correct here, so I am changing them

!   if(ib>=1 .AND. ib<=imt .AND. jb>=1 .AND. jb<=jmt .AND. km+1-kmt(ib,jb)<=kb .AND. kb>=1 ) then
if(ib>1 .AND. ib<imt .AND. jb>1 .AND. jb<jmt .AND. km>=kb .AND. kb>=1 ) then
    tryAgain = .TRUE.
else
    print *,'outside model domain in diffusion',ib,jb,km,kb
!     print *,'outside model domain in diffusion',ib,jb,km+1-kmt(ib,jb),kb
end if
  
itno=0
do while(tryAgain)
    itno=itno+1
    ! find random displacement 
    call displacement(xd, yd, zd, ib, jb, kb, dt,imt,jmt,h,ah,av,dxv,dyu,do3d,doturb)
!     CALL displacement(xd, yd, zd, ib, jb, kb, dt,imt,jmt,kmt,h,ah,av)
!     print *,'xd=',xd,' yd=',yd,' zd=',zd
!     print *,'dzt(ib,jb,:,1)=',dzt(ib,jb,:,1)
    ! Convert displacement from meters to model coordinates
    ! KMT note: following example in master code for orc project
    xd = xd/dxv(ib,jb)
    yd = yd/dyu(ib,jb)
    zd = zd/dzt(ib,jb,kb,1) !KMT: this should be better than dz()
!     print *,'dzt(ib,jb,kb,1)=',dzt(ib,jb,kb,1)

    ! TO DO: improve so that dzt is calculated including free surface and for correct time
!     zd = zd/dz(kb) !should be replaced for bottom box and include ssh (original note)

!     print *,'xd=',xd,' yd=',yd,' zd=',zd

    ! Update position temporarily
    tmpX = x1 + xd
    tmpY = y1 + yd
    tmpZ = z1 + zd
    ! Update box number temporarily
    tmpi = int(tmpX) + 1
    tmpj = int(tmpY) + 1
    tmpk = int(tmpZ) + 1
!     print *,'temp dzt(ib,jb,kb,1)=',dzt(tmpi,tmpj,tmpk,1)

!     ! Check if particle is on an open boundary. Using rco example.
!     if(tmpi==1 .AND. tmpj>=1 .AND. tmpj<=jmt .AND. km+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
!         tryAgain = .FALSE.
!     end if

    ! check that column is deep enough  
    ! KMT: this is the same check as above that looks wrong and I changed
!     if( 1<=tmpi .AND. tmpi<=imt .AND. 1<=tmpj .AND. tmpj<=jmt .AND. km+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
    if(tmpi>1 .AND. tmpi<imt .AND. tmpj>1 .AND. tmpj<jmt .AND. km>=tmpk .AND. tmpk>=1 ) then
!         print *,'km=',km,' km+1-kmt(tmpi,tmpj)=',km+1-kmt(tmpi,tmpj),' tmpk=',tmpk
        tryAgain = .FALSE. 
!         print *,'found new position for drifter that is within model domain, so move on.'
        ! if false then a new position for the particle has been found and we exit the loop
    end if 

    ! KMT: Need to check to see if new position is on masked land.
    ! Can do this by checking dzt at the new position, which is zero when on land.
    ! If we are on land, keep tryAgain = .TRUE. so that this displacement choice is
    ! not used. This can alter the value of tryAgain just given above.
    if(dzt(tmpi,tmpj,tmpk,1)==0.) then
        tryAgain = .TRUE.
    endif

    ! If tryAgain is still true, the particle is outside model area. The
    ! displacement is not saved, but we make a new try to displace.
    
    ! "Infinite loop?"
    if(itno>=100000 .AND. tryAgain) then
        tryAgain = .FALSE.
        write(*,*)"Particle stuck in infinite diffusion loop. No diffusion added.",ib,jb,kb
        stop 34956
        tmpX=x1 ; tmpY=y1 ; tmpZ=z1
        tmpi=ib ; tmpj=jb ; tmpk=kb 
    end if
    
  enddo
  
! Update return position
! if(abs(x1-tmpX).gt.5. .and. abs(x1-tmpX).lt.1000.) print *,x1,tmpX,xd,ib,jb,dxv(ib,jb)
x1 = tmpX
y1 = tmpY
ib = tmpi
jb = tmpj

! print *,'After diffusion'
! print *,'x1=',x1,' y1=',y1,' ib=',ib,' jb=',jb

if(do3d==1) then
    z1 = tmpZ
    kb = tmpk
endif  

if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! make sure the index corresponds to the right grid cell

!print *,'slut',itno,kmt(ib,jb),ib,jb,kb,x1,y1,z1
  
end subroutine diffuse

!===============================================================================
! Calculate a random displacement
! (sqrt(-4ah*dt*log(1-q1))*cos(2PIq2),
!  sqrt(-4ah*dt*log(1-q1))*sin(2PIq2),
!  sqrt(-4av*dt*log(1-q3))*cos(2PIq4))
! where av and ah are set in run.in, dt is the model time step and q1,q2,q3,q4 are random
! numbers between 0 and 1.
!
! Arguments :
! xd, yd, zd : Variables in which the displacement will be stored
! dt: Model time step
!===============================================================================
subroutine displacement(xd, yd, zd, ib, jb, kb, dt,imt,jmt,h,ah,av,dxv,dyu,do3d,doturb) 

implicit none
  
integer,        intent(in)                          :: ib,jb,kb,imt,jmt   ! box indices
integer,        intent(in)                          :: do3d, doturb
real(kind=8),   intent(in)                          :: dt,ah,av
real*8,         intent(in), dimension(imt,jmt)      :: h
real*8,         intent(in), dimension(imt-1,jmt)    :: dyu
real*8,         intent(in), dimension(imt,jmt-1)    :: dxv
real(kind=8),   intent(out)                         :: xd, yd, zd
real(kind=8),   parameter                           :: pi = 3.14159265358979323846
real(kind=8)                                        :: q1, q2, q3, q4, r, elip, xx, yy, hp, hm

real(kind=8)                                        :: grdx, grdy, grad, theta
integer                                             :: ip,im,jp,jm
    
! random generated numbers between 0 and 1
q1 = rand() 
q2 = rand()
q3 = rand()
q4 = rand()

! Horizontal displacements in meters
R = sqrt(-4*ah*dt*log(1-q1))
xd = R * cos(2*PI*q2)
yd = R * sin(2*PI*q2)

! Vertical displacement in meters
if(do3d==1) then
    R = sqrt(-4*av*dt*log(1-q3))
    zd = R*cos(2*PI*q4)
else
    zd = 0.
endif

! print *,' '
! print *,'circle'
! print *,'xd=',xd,' yd=',yd

if(doturb==3) then
    !_______________________________________________________________________________
    ! The diffusion is here set on an ellipse instead of a circle
    ! so that the diffusion is higher along the isobaths and weaker
    ! in the perpendicular direction (up/downhill)

    ! for index in i direction for one cell positive from the current cell
    ip=ib+1
    ! KMT: Is this a periodic BC? Replace with just the edge value so it won't 
    ! elliptical but just at the edge.
    ! I think the check should happen for ip==imt since in the x direction for u
    ! the grid only goes up to imt-1 for the arakawa c grid.
    if(ip.eq.imt) ip = imt-1 ! use ib instead of ib+1, and here ib=imt-1
!     if(ip.eq.imt+1) ip=1

    ! for index in i direction for one cell negative from the current cell
    im=ib-1
    ! KMT: Is this a periodic BC? Replace with just the edge value so it won't 
    ! elliptical but just at the edge.
    if(im.eq.0) im = 1 ! use ib instead of ib-1, and here ib=1
!     if(im.eq.0) im=imt

     ! for index in j direction for one cell positive from the current cell
    jp=jb+1
    ! KMT: In v they don't use period BC. 
    ! I think the check should happen for jp==jmt since in the y direction for v
    ! the grid only goes up to jmt-1 for the arakawa c grid.
    if(jp.eq.jmt) jp = jmt-1 ! use jb instead of jb+1, and here jb=jmt-1
!     if(jp.eq.jmt+1) jp=jmt

    ! for index in j direction for one cell negative from the current cell
    jm=jb-1
    ! KMT: In v they don't use period BC and this check is fine as is.
    if(jm.eq.0) jm = 1

    ! just depth gradient
    ! grdx=float(kmt(ip,jb)-kmt(im,jb)) ! zonal      depth gradient (zw should be used)
    ! grdy=float(kmt(ib,jp)-kmt(ib,jm)) ! meridional depth gradient
    ! KMT: switching the bathymetry/depths array and grid distances in here
    ! since user manual (pg 21) shows grdx=\Delta h/\Delta x, for example
!     grdx=(h(ip,jb)-h(im,jb))/(0.5*(dxv(im,jb)+dxv(ib,jb))+0.5*(dxv(ib,jb)+dxv(ip,jb))) ! zonal depth gradient
    ! Same as easier-to-read previous line but more efficient:
    grdx=(h(ip,jb)-h(im,jb))/(0.5*dxv(im,jb) + 0.5*dxv(ip,jb) + dxv(ib,jb)) ! zonal depth gradient

!     grdy=(h(ib,jp)-h(ib,jm))/(0.5*(dyu(ib,jm)+dyu(ib,jb)) + 0.5*(dyu(ib,jb)+dyu(ib,jp))) ! meridional depth gradient
    ! Same as easier-to-read previous line but more efficient:
    grdy=(h(ib,jp)-h(ib,jm))/(0.5*dyu(ib,jm) + dyu(ib,jb) + 0.5*dyu(ib,jp)) ! meridional depth gradient

    grad=dsqrt(grdx**2+grdy**2)       ! total depth gradient

    ! angle between the eastward direction and 
    ! the constant depth direction (=0 if only meridional slope)
    if(grad.eq.0.d0) then
        ! theta=0.d0 + pi/2.d0
        theta=0.d0 
    else
        ! theta=dasin(grdx/grad) +pi/2.d0 ! varför 90 grader mer? ny2
        theta=dasin(grdx/grad) ! ny1 
    endif

    ! elliptic horizontal distribution of the diffusion
    !elip=dabs(1.d2*grad)+1.d0 ! för 1/h
    !elip=dabs(1.d-2*grad)+1.d0 ! för h
    elip=dabs(grad)+1.d0 ! för h

    !elip=amin1(10.d0,grad) ! gives slightly higher relative dispersion 
    elip=amin1(5.d0,elip) ! 

    !print *,'elip=',elip,grad,theta

    !The circular disk is stretched into an elliptic disk 
    ! with unchanges surface area by
    xx=xd*elip
    yy=yd/elip  

    ! coordinate transformation to put the diffusion on an 
    ! ellipse with the maxium diffusion along the isobaths
    xd= xx*dcos(theta)-yy*dsin(theta)
    yd=-xx*dsin(theta)+yy*dcos(theta)

!     print *,'ellipse'
!     print *,'xd=',xd,' yd=',yd
    ! print *,'theta=',theta,' xx=',xx,' yy=',yy,' grdx=',grdx,' grdy=',grdy

    !print *,'theta=',theta*180./pi,elip,grdx/grad

    !if(jb.gt.400) stop 3096
    !_______________________________________________________________________________
endif

end subroutine displacement

