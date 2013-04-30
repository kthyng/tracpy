SUBROUTINE diffuse(x1, y1, z1, ib, jb, kb, dt,IMT,JMT,KM,kmt,dxv,dyu,dzt,h,Ah,Av)
!KMT had to move inside subroutine
#ifdef diffusion 

!============================================================================
! Add a small displacement to a particle s.t. it is still in the model area.
! 
! Arguments
! x1, y1, z1 : Current position of the particle. Are updated by the subroutine.
! ib, jb, kb : Current box for the particle. Are updated by the subroutine.
! dt : The model time step
!
!============================================================================
  
  implicit none
INTEGER, PARAMETER                :: DP = SELECTED_REAL_KIND(15, 307)
  INTEGER, intent(in)                     :: IMT,JMT,KM           
  INTEGER, intent(in out)                     :: ib,jb,kb           ! Box indices
  INTEGER                     :: itno               ! Number of iterations
  REAL                        :: xd, yd, zd         ! Displacement
  REAL                        :: tmpX, tmpY, tmpZ   ! Temporal position
  INTEGER                     :: tmpi, tmpj, tmpk   ! Temporal box indices
  REAL (KIND=DP), INTENT(in OUT) :: x1, y1, z1         ! Final position
  REAL (KIND=DP), INTENT(IN)  :: dt                 ! Time step
  LOGICAL                     :: tryAgain           ! Tells whether to continue displace
integer,   intent(in),     dimension(IMT,JMT)              :: kmt,h
  real(kind=8), intent(in) :: Ah,Av
integer,   intent(in),     dimension(IMT-1,JMT)              :: dyu
integer,   intent(in),     dimension(IMT,JMT-1)              :: dxv
real(kind=8),   intent(in),     dimension(IMT,JMT,KM,2)         :: dzt

  tryAgain = .FALSE.
  
  ! Is particle within model area?
  ! KMT: I do not understand how the vertical checks are correct here, so I am changing them
!   if(ib>=1 .AND. ib<=IMT .AND. jb>=1 .AND. jb<=JMT .AND. KM+1-kmt(ib,jb)<=kb .AND. kb>=1 ) then
  if(ib>=1 .AND. ib<=IMT .AND. jb>=1 .AND. jb<=JMT .AND. KM>=kb .AND. kb>=1 ) then
    tryAgain = .TRUE.
  else
    print *,'outside model domain in diffusion',ib,jb,KM+1-kmt(ib,jb),kb
  end if
  
  itno=0
  do while(tryAgain)
      itno=itno+1
    ! find random displacement 
    CALL displacement(xd, yd, zd, ib, jb, kb, dt,IMT,JMT,h,Ah,Av)
!     CALL displacement(xd, yd, zd, ib, jb, kb, dt,IMT,JMT,kmt,h,Ah,Av)
!     print *,'xd=',xd,' yd=',yd,' zd=',zd
!     print *,'dzt(ib,jb,:,1)=',dzt(ib,jb,:,1)
    ! Convert displacement from meters to model coordinates
    ! KMT note: following example in master code for orc project
    xd = xd/dxv(ib,jb)
    yd = yd/dyu(ib,jb)
    zd = zd/dzt(ib,jb,kb,1) !KMT: this should be better than dz()
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

    ! Check if particle is on an open boundary. Using rco example.
    if(tmpi==1 .AND. tmpj>=1 .AND. tmpj<=JMT .AND. KM+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
      tryAgain = .FALSE.
    end if

    ! check that column is deep enough  
    ! KMT: this is the same check as above that looks wrong and I changed
!     if( 1<=tmpi .AND. tmpi<=IMT .AND. 1<=tmpj .AND. tmpj<=JMT .AND. KM+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
    if(tmpi>=1 .AND. tmpi<=IMT .AND. tmpj>=1 .AND. tmpj<=JMT .AND. KM>=tmpk .AND. tmpk>=1 ) then
            print *,'KM=',KM,' KM+1-kmt(tmpi,tmpj)=',KM+1-kmt(tmpi,tmpj),' tmpk=',tmpk
            tryAgain = .FALSE. 
            ! if false then a new position for the particle has been found and we exit the loop
    end if 
    
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
#ifndef twodim   
  z1 = tmpZ
  kb = tmpk
#endif  

if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! make sure the index corresponds to the right grid cell

!print *,'slut',itno,kmt(ib,jb),ib,jb,kb,x1,y1,z1
  
#endif 
!KMT had to move inside subroutine (from very bottom of script)

END SUBROUTINE diffuse

!===============================================================================
! Calculate a random displacement
! (sqrt(-4Ah*dt*log(1-q1))*cos(2PIq2),
!  sqrt(-4Ah*dt*log(1-q1))*sin(2PIq2),
!  sqrt(-4Av*dt*log(1-q3))*cos(2PIq4))
! where Av and Ah are set in run.in, dt is the model time step and q1,q2,q3,q4 are random
! numbers between 0 and 1.
!
! Arguments :
! xd, yd, zd : Variables in which the displacement will be stored
! dt: Model time step
!===============================================================================
SUBROUTINE displacement(xd, yd, zd, ib, jb, kb, dt,IMT,JMT,h,Ah,Av) 
  IMPLICIT NONE
  
INTEGER, PARAMETER                :: DP = SELECTED_REAL_KIND(15, 307)
  REAL                  :: q1, q2, q3, q4, R
  REAL, INTENT(OUT)     :: xd, yd, zd
  REAL (KIND=DP), INTENT(IN)  :: dt
  REAL, PARAMETER       :: PI = 3.14159265358979323846
  INTEGER, intent(in)               :: ib,jb,kb,IMT,JMT   ! Box indices
#ifdef anisodiffusion   
  REAL*8                :: Rx, Ry, grdx, grdy, grad, theta, elip, xx, yy, hp, hm
  INTEGER               :: ip,im,jp,jm
#endif
  real(kind=8), intent(in) :: Ah,Av
integer,   intent(in),     dimension(IMT,JMT)              :: h
    
! random generated numbers between 0 and 1
  q1 = rand() 
  q2 = rand()
  q3 = rand()
  q4 = rand()

! Horizontal displacements in meters
  R = sqrt(-4*Ah*dt*log(1-q1))
  xd = R * cos(2*PI*q2)
  yd = R * sin(2*PI*q2)

!   print *,'R=',R,' q1=',q1,' q2=',q2,' Ah=',Ah,' dt=',dt
    
  ! Vertical displacement in meters
#ifndef twodim   
  R = sqrt(-4*Av*dt*log(1-q3))
  zd = R*cos(2*PI*q4)
!   print *,'R=',R,' Av=',Av,' dt=',dt,' q3=',q3,' q4=',q4,' zd=',zd
#else
    zd = 0.
#endif

#ifdef anisodiffusion
!_______________________________________________________________________________
! The diffusion is here set on an ellipse instead of a circle
! so that the diffusion is higher along the isobaths and weaker
! in the perpendicular direction (up/downhill)

ip=ib+1
if(ip.eq.IMT+1) ip=1
im=ib-1
if(im.eq.0) im=IMT
jp=jb+1
if(jp.eq.JMT+1) jp=JMT
jm=jb-1
if(jm.eq.0) jm=1


! just depth gradient
! grdx=float(kmt(ip,jb)-kmt(im,jb)) ! zonal      depth gradient (zw should be used)
! grdy=float(kmt(ib,jp)-kmt(ib,jm)) ! meridional depth gradient
! KMT: switching the bathymetry/depths array in here
grdx=float(h(ip,jb)-h(im,jb)) ! zonal      depth gradient
grdy=float(h(ib,jp)-h(ib,jm)) ! meridional depth gradient

! print *,'kmt(ip,jb)=',kmt(ip,jb),' kmt(im,jb)=',kmt(im,jb),' kmt(ib,jp)=',kmt(ib,jp),' kmt(ib,jm)=',kmt(ib,jm)

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

! print *,'theta=',theta,' xx=',xx,' yy=',yy,' grdx=',grdx,' grdy=',grdy

!print *,'theta=',theta*180./pi,elip,grdx/grad

!if(jb.gt.400) stop 3096
!_______________________________________________________________________________
#endif

END SUBROUTINE displacement

