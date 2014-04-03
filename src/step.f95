SUBROUTINE step(xstart,ystart,zstart,tseas, &
                & uflux,vflux,ff,imt,jmt,km,kmt,dzt,dxdy,dxv,dyu,h, &
                & ntractot,xend,yend,zend,flag,ttend, &
                & iter,ah,av,do3d,doturb, doperiodic, dostream, N, T0, &
                & ut, vt)
! SUBROUTINE step(xstart,ystart,zstart,tseas, &
!                 & uflux,vflux,ff,imt,jmt,km,kmt,dzt,dxdy,dxv,dyu,h, &
!                 & ntractot,xend,yend,zend,iend,jend,kend,flag,ttend, &
!                 & iter,ah,av,do3d,doturb, dostream, N, T0, &
!                 & ut, vt)

!============================================================================
! Loop to step a numerical drifter forward for the time tseas between two 
! model outputs. uflux and vflux are interpolated in time to iter steps 
! between the model outputs, and drifter tracks will be output at 1/iter
! intervals.
!
! This is put together with code pieces from tracmass to be 
! controlled by outside wrapping Python code.
! Kristen M. Thyng Feb 2013
!
!  Input:
!
!    xstart         : Starting x,y,z position of drifters for this set of two
!    ystart           model outputs, in fraction of the grid cell in each
!    zstart           direction. [ntractoc]
!    tseas          : Time between model-output velocity fields that are input.
!                     Also max time for this loop (seconds)
!    uflux          : u velocity (zonal) flux field, two time steps [ixjxkxt]
!    vflux          : v velocity (meridional) flux field, two time steps [ixjxkxt]
!    ff             : time direction. ff=1 forward, ff=-1 backward
!    imt,jmt,km     : grid index sizing constants in (x,y,z), are for 
!                     horizontal and vertical rho grid [scalar]
!    kmt            : Number of vertical levels in horizontal space [imt,jmt]
!    dzt            : Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
!    dxdy           : Horizontal area of cells defined at cell centers [imt,jmt]
!    dxv            : Horizontal grid cell walls areas in x direction [imt,jmt-1]
!    dyu            : Horizontal grid cell walls areas in y direction [imt-1,jmt]
!    h              : Depths [imt,jmt]
!    ntractoc       : Total number of drifters
!    iter           : number of interpolation time steps to do between model outputs, 
!                     was called nsteps in run.py. This controls the max time step 
!                     between drifter steps but not how often the drifter location is sampled.
!    ah             : horizontal diffusion in m^2/s. Input parameter.
!                   : For -turb and -diffusion flags
!    av             : vertical diffusion in m^2/s. Input parameter. For -diffusion flag.
!    do3d           : Flag to set whether to use 3d velocities or not. 
!                   : do3d=0 makes the run 2d and do3d=1 makes the run 3d
!    doturb         : Flag to set whether or not to use turb/diff and which kind if so
!                   : doturb=0 means no turb/diffusion,
!                   : doturb=1 means adding parameterized turbulence
!                   : doturb=2 means adding diffusion on a circle
!                   : doturb=3 means adding diffusion on an ellipse (anisodiffusion)
!    doperiodic     : Whether to use periodic boundary conditions for drifters and, if so, on which walls.
!                   : 0: do not use periodic boundary conditions
!                   : 1: use a periodic boundary condition in the east-west/x/i direction
!                   : 2: use a periodic boundary condition in the north-south/y/j direction
!    dostream       : Either calculate (dostream=1) or don't (dostream=0) the
!                     Lagrangian stream function variables.
!    N              : The number of samplings of the drifter tracks will be N+1 total.
!  Optional inputs:
!    T0             : (optional) Initial volume transport of drifters (m^3/s)
!    ut, vt     : (optional) Array aggregating volume transports as drifters move [imt,jmt]
!
!  Output:
!
!    flag           : set to 1 for a drifter if drifter shouldn't be stepped in the 
!                     future anymore
!    xend           : the new grid fraction position of drifters in x/y/z [ntractot,iter]
!    yend       
!    zend       
!    ttend          : time in seconds relative to the code start for when particles are
!                     output. [ntractot,iter]
!
!  Other parameters used in function:
!
!    istart         : Starting grid index of drifters in x,y,z for this set of two
!    jstart           model outputs, should be integers and be for grid cell wall
!    kstart           just beyond the drifter location. These are inferred 
!                     from x/y/zstart arrays. [ntractoc]
!    rg             : rg=1-rr for time interpolation between time steps. Controls how much
!                   : of later time step is used in interpolation.
!    rr             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation.
!    rb             : time interpolation constant between 0 and 1. Controls how much
!                   : of earlier time step is used in interpolation (for next time 
!                     at the end of the loop?).
!    dsc            : Not sure what this is right now
!    dstep          : Model time step between time interpolation steps between model outputs
!    dtmin          : time step based on the interpolation step between model output times.
!                   : sets a limit on the time step that a drifter can go. (seconds)
!    dxyz           : Volume of drifter's grid cell (m^3)
!    dt             : time step of trajectory iteration in seconds  
!    dsmin          : time step based on the interpolation step between model output times.
!                   : sets a limit on the time step that a drifter can go. (s/m^3)
!                     dtmin/(dx*dy*dz)
!    ds             : crossing time to reach the grid box wall (units=s/m3). ds=dt/(dx*dy*dz)
!    dse,dsw        : crossing times for drifter to reach the east and west grid box wall
!    dsn,dss        : crossing times for drifter to reach the north and south grid box wall
!    dsu,dsd        : crossing times for drifter to reach the up and down grid box wall
!    x0,y0,z0       : original non-dimensional position in the i,j,k-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!    x1,y1,z1       : updated non-dimensional position in the i,j,k-direction
!                     of particle (fractions of a grid box side in the 
!                     corresponding direction)
!    tt             : time in seconds of the trajectory relative to the code start
!    tss            : Counter for iterations between model outputs. Counts up to iter I think.
!    ts             : Time step of time interpolation between model outputs, non-dimensional,
!                     has values between 0 and 1
!    ntrac          : Index in arrays for current drifter in ntracLoop
!    niter          : Index in arrays for current point in drifter loop in niterLoop
!    Ni             : Counter for N for when to write to file
!    ia,ja,ka       : original position in integers
!    iam            : index for grid index ia-1
!    ib,jb,kb       : new position in grid space indices
!    errCode        : Keeps track of errors that occur. Currently not being well utilized.
!    upr            : parameterized turbulent velocities u', v', and w'
!                     optional because only used if using turb flag for diffusion
!                     size [6,2]. The 2nd dimension is for two time steps.
!                     The 1st dimension is: [u'_ia,u'_ia-1,v'_ja,v'_ja-1,w'_ka,w'_ka-1]
!    rwn, rwp       : Interpolation constants for when outputting.
!
!====================================================================

implicit none

integer,    intent(in)                                  :: ff, imt, jmt, km, ntractot, iter
integer,    intent(in)                                  :: do3d, doturb, doperiodic, dostream, N
integer,    intent(in),     dimension(imt,jmt)          :: kmt
real*8,     intent(in),     dimension(imt-1,jmt)        :: dyu
real*8,     intent(in),     dimension(imt,jmt-1)        :: dxv
real*8,     intent(in),     dimension(ntractot)         :: xstart, ystart, zstart
real*8,     intent(in),     dimension(imt-1,jmt,km,2)   :: uflux
real*8,     intent(in),     dimension(imt,jmt-1,km,2)   :: vflux
real*8,     intent(in),     dimension(imt,jmt,km,2)     :: dzt
real*8,     intent(in),     dimension(imt,jmt)          :: dxdy, h
real*8,     intent(in)                                  :: tseas, ah, av

integer,    intent(out),    dimension(ntractot)         :: flag
real*8,     intent(out),    dimension(ntractot,N)       :: xend, yend, zend, ttend
! integer,    intent(out),    dimension(ntractot,N)    :: iend, jend, kend
integer,                    dimension(ntractot)         :: istart, jstart, kstart

real*8,                     dimension(0:km,2)           :: wflux
real*8                                                  :: rr, rg, rb, dsc, &
                                                           dstep, dtmin, dxyz, dt, &
                                                           dsmin, ds, dse, dsw, dsn, &
                                                           dss, dsd, dsu, &
                                                           x0, y0, z0, x1, y1, z1, &
                                                           tt, tss, ts, rwn, rwp
integer                                                 :: ntrac, niter, ia, ja, ka, &
                                                           iam, ib, jb, kb, errCode, Ni
real*8,     parameter                                   :: UNDEF=1.d20

real*8, dimension(6,2)                                    :: upr
! The following are for calculating Lagrangian stream functions
real*8, optional, intent(in), dimension(ntractot) :: T0
real*8, optional, intent(inout), dimension(imt-1,jmt) :: ut
real*8, optional, intent(inout), dimension(imt,jmt-1) :: vt

!f2py intent(out) :: ut, vt 

! Set indices for x/y/z grid locations to be the ceiling of the grid indices
istart = ceiling(xstart)
jstart = ceiling(ystart)
kstart = ceiling(zstart)

! controls the iterative stepping between the two input model outputs
dstep = 1.d0/dble(iter)
dtmin = dstep*tseas 

flag = 0

!=======================================================
!=== Loop over all trajectories and calculate        ===
!=== a new position for this time step.              ===
!=======================================================

ntracLoop: do ntrac=1,ntractot  
!     print *,'ntractot=',ntractot
    ! start track off with no error
!     errCode = 0
    ! Counter for sub-interations for each drifter. In the source, this was read in but I am not sure why.
    niter = 0
    Ni = 1 ! counter for when to write to file
    ! Initialize drifter parameters for this drifter (x0,y0,x0 will be updated from these at the start of the loop)
    x1 = xstart(ntrac)
    y1 = ystart(ntrac)
    z1 = zstart(ntrac)
    tt = 0 ! time of trajectory in seconds
    ts = 0.d0 ! model dataset time step of trajectory
    tss = 0.d0
    ib = istart(ntrac)
    jb = jstart(ntrac)
    kb = kstart(ntrac)
  
    ! ===  start loop for each trajectory ===
    !         scrivi=.true.
    niterLoop: do 
        niter=niter+1 ! iterative step of trajectory
!         print *,'niter=',niter
        rg = ts 
!         rg=dmod(ts,1.d0) ! time interpolation constant between 0 and 1
        rr = 1.d0-rg
!         print *,'ts=',ts,' rg=',rg,' rr=',rr
        ! Update particle indices and locations
        x0 = x1
        y0 = y1
        z0 = z1
        ia = ib
        iam = ia-1
        ja = jb
        ka = kb
!       print '(a,f8.1,a,f8.1,a,f12.2)','uflux(ia,ja,ka,1)=',uflux(ia,ja,ka,1),' uflux(ia ,ja,ka,2)=',uflux(ia ,ja,ka,2),' dxyz=',dxyz
!       print '(a,f8.1,a,f8.1,a,f12.2)','vflux(ia,ja,ka,1)=',vflux(ia,ja,ka,1),' vflux(ia ,ja,ka,2)=',vflux(ia ,ja,ka,2),' dxyz=',dxyz
!       print *,'ia=',ia,' ja=',ja,' ka=',ka
!       print *,'iter=',iter,'ntractot=',ntractot
!         ! start track off with no error (maybe export these later to keep track of)
        errCode = 0

! print *,'ib=',ib,' ia=',ia,' jb=',jb,' ja=',ja

        call calc_dxyz(ib,jb,kb,rr,imt,jmt,km,dzt,dxdy,dxyz)
        ! I am putting the error checks directly in the loop for convenience
        !         call errorCheck('dxyzError',errCode,dxyz,flag)
        ! Check the grid box volume
        if(dxyz == 0.d0) then
            print *,'====================================='
            print *,'ERROR: dxyz is zero'
            print *,'-------------------------------------'
            print *,'ntrac=',ntrac!,' ints=', ints
            print *,'ib=',ib,'jb=',jb,'kb=',kb
            !           print *,'kmt=',kmt(ib,jb)
            !           print *,'dz=',dz(kb)
            !           print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb)
            print *,'-------------------------------------'
            print *,'The trajectory is killed'
            print *,'====================================='
            errCode = -39
            flag(ntrac) = 1 
        end if

!         call errorCheck('coordBoxError' ,errCode)
        ! ===  Check that coordinates belongs to   ===
        ! ===  correct box. Valuable for debugging ===
        if( dble(ib-1).gt.x1 .or. dble(ib).lt.x1 )  then
            print *,'========================================'
            print *,'ERROR: Particle overshoot in i direction'
            print *,'----------------------------------------'
            print *,ib-1,x1,ib,ntrac,ib,jb,kb
            x1=dble(ib-1)+0.5d0
            ib=idint(x1)+1
            print *,'error i',ib-1,x1,ib,ntrac,ib,jb,kb
            print *,y1,z1
            print *,'-------------------------------------'
            print *,'The run is terminated'
            print *,'====================================='             
            errCode = -42
            stop
        else if( dble(jb-1).gt.y1 .or. dble(jb).lt.y1 )  then
            print *,'========================================'
            print *,'ERROR: Particle overshoot in j direction'
            print *,'----------------------------------------'
            print *,'error j',jb-1,y1,jb,ntrac,x1,z1
            print *,'error j',jb-1,y1,jb,ntrac,ib,jb,kb
            print *,'-------------------------------------'
            print *,'The run is terminated'
            print *,'====================================='    
            errCode = -44
            stop
        else if((dble(kb-1).gt.z1.and.kb.ne.km).or. & 
             dble(kb).lt.z1 ) then
            print *,'========================================'
            print *,'ERROR: Particle overshoot in k direction'
            print *,'----------------------------------------'
            print *,'error k',kb-1,z1,kb,ntrac,x1,y1
            print *,'error k',kb-1,z1,kb,ntrac,ib,jb,kb
            print *,'-------------------------------------'
            print *,'The run is terminated'
            print *,'====================================='
            errCode = -46
            stop
        end if

!         call errorCheck('infLoopError'  ,errCode)

        if(niter.gt.30000) then ! break infinite loops
!         if(niter-nrj(ntrac,4).gt.30000) then ! break infinite loops
!             nloop=nloop+1             
            print *,'====================================='
            print *,'Warning: Particle in infinite loop '
            print *,'ntrac:',ntrac
            print *,'niter:',niter!,'nrj:',nrj(ntrac,4)
!             print *,'dxdy:',dxdy(ib,jb),'dxyz:',dxyz
!             print *,'kmt:',kmt(ia-1,ja-1),'dz(k):',dz(ka-1)
            print *,'ia=',ia,' ib=',ib,' ja=',ja,' jb=',jb, & 
               ' ka=',ka,' kb=',kb
            print *,'x1=',x1,' x0=',x0,' y1=',y1,' y0=',y0, & 
               ' z1=',z1,' z0=',z0
            print *,'dse=',dse,' dsw=',dsw,' dsn=',dsn,' dss=',dss,'dsmin=',dsmin
            print *,'-------------------------------------'
            errCode = -48
            flag(ntrac) = 1 ! want to continue drifters if time was the limiting factor here
        end if

        if (errCode.ne.0) print *,'Error code=',errCode

        if (errCode.ne.0) cycle ntracLoop

        !==============================================! 
        ! calculate the 3 crossing times over the box  ! 
        ! choose the shortest time and calculate the   !
        ! new positions                                !
        !                                              !
        !-- solving the differential equations ---     !
        ! note:                                        !
        ! space variables (x,...) are dimensionless    !
        ! time variables (ds,...) are in seconds/m^3   !
        !==============================================!

        dsmin=dtmin/dxyz
!         print *,'dsmin=',dsmin,' dtmin=',dtmin,' dxyz=',dxyz

        ! ! === calculate the turbulent velocities ===
        if(doturb==1) then
            call turbuflux(ia,ja,ka,rr,dtmin,ah,imt,jmt,km,uflux,vflux,wflux,ff,do3d,upr)
        end if

!         ! === calculate the vertical velocity ===
!              print *,''
!            print *,'ntrac=',ntrac
!         print *,'before vertvel ========================================'  
!         print '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)','x0=',x0,' x1=',x1,&
!             ' y0=',y0,' y1=',y1,' z0=',z0,' z1=',z1
!         print '(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)','ia=',ia,' ib=',ib,&
!             ' ja=',ja,' jb=',jb,' ka=',ka,' kb=',kb
!         print *,'imt=',imt,' size(uflux)=',size(uflux,1),' size(vflux)=',size(vflux,1)
        call vertvel(rr,ia,ja,ka,imt,jmt,km,ff,uflux,vflux,do3d,wflux)
!         print *,'returned from vertvel ========================================'
!         print *,'wflux=',wflux

!  print *,'ja=',ja,' jb=',jb,x1,y1
        if(doturb==1) then
            call cross(1,ia,ja,ka,x0,dse,dsw,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb,upr) ! zonal
            call cross(2,ia,ja,ka,y0,dsn,dss,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb,upr) ! meridional
            call cross(3,ia,ja,ka,z0,dsu,dsd,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb,upr) ! vertical
        else 
            call cross(1,ia,ja,ka,x0,dse,dsw,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb) ! zonal
            call cross(2,ia,ja,ka,y0,dsn,dss,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb) ! meridional
            call cross(3,ia,ja,ka,z0,dsu,dsd,rr,uflux,vflux,wflux,ff,km,jmt,imt,do3d,doturb) ! vertical
        endif
!         print *, 'dse=',dse,' dsw=',dsw,' dsn=',dsn,' dss=',dss,' dsmin=',dsmin
          ds=dmin1(dse,dsw,dsn,dss,dsu,dsd,dsmin)
!         print *,'Before calc_time: ds=',ds,' tss=',tss,' ts=',ts,' tt=',tt,' dtmin=',dtmin

!         call errorCheck('dsCrossError', errCode)
     ! === Can not find any path for unknown reasons ===
        if(ds == UNDEF .or.ds == 0.d0)then 
            print *, " "
            print *, " "
            print *,'==================================================='
            print *,'Warning: not find any path for unknown reason '
            print *, " "
!             write (*,FMT='(A, 5E9.2)'),' ds=',ds,dse,dsw,dsn,dss
!             write (*,FMT='(4E9.2)'), dsu,dsd,dsmin,dxyz
!             print *,'---------------------------------------------------'
!             print *,"   ntrac = ",ntrac
!             write (*,'(A7, I10, A7, I10, A7, I10)'), & 
!                 ' ia= ', ia, ' ja= ', ja, ' ka= ', ka
!             write (*,'(A7, I10, A7, I10, A7, I10)'), & 
!                 ' ib= ', ib, ' jb= ', jb, ' kb= ', kb
!             write (*,'(A7, F10.3, A7, F10.3, A7, F10.3)'), & 
!                 ' x0= ', x0, ' y0= ', y0, ' z0= ', z0
!             write (*,'(A7, F10.3, A7, F10.3, A7, F10.3)'), & 
!                 ' x0= ', x0, ' y0= ', y0, ' z0= ', z0
!             write (*,'(A7, I10, A7, I10, A7, I10)'), & 
!                 ' k_inv= ', km+1-kmt(ia,ja), ' kmt= ', kmt(ia,ja), &
!                 'lnd= ', mask(ia,ja)
            print *,'---------------------------------------------------'
            print *,'The trajectory is killed'
            print *,'==================================================='
    !         nerror=nerror+1
    !         nrj(ntrac,6)=1
            flag(ntrac) = 1
            errCode = -56
        end if
        if (errCode.ne.0) cycle ntracLoop

        call calc_time(ds,dsmin,dt,dtmin,tss,tseas,ts,tt,dxyz,dstep,iter,rb,dsc)

!         print *,'After calc_time: ds=',ds,' tss=',tss,' ts=',ts,' tt=',tt,' dt=',dt,' dtmin=',dtmin

        ! === calculate the new positions ===
        ! === of the trajectory           ===  
!         if(ntrac==42) then
!              print *,''
!            print *,'before pos'  
!             print '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)','x0=',x0,' x1=',x1,&
!                 ' y0=',y0,' y1=',y1,' z0=',z0,' z1=',z1
!             print '(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)','ia=',ia,' ib=',ib,&
!                 ' ja=',ja,' jb=',jb,' ka=',ka,' kb=',kb
!         endif

        if(doturb==1) then
            call pos(ia,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1,ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dsc,&
                ff,imt,jmt,km,rr,rb,uflux,vflux,wflux,do3d,doturb,upr)
        else
            call pos(ia,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1,ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dsc,&
                ff,imt,jmt,km,rr,rb,uflux,vflux,wflux,do3d,doturb)
        endif

!         print *,''
!         print *,'x1/x0=',x1/x0,' y1/y0=',y1/y0

!         if(ntrac==42) then
!             print *,'after pos'
!             print '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)','x0=',x0,' x1=',x1,&
!                 ' y0=',y0,' y1=',y1,' z0=',z0,' z1=',z1
!             !         print '(a,f10.2,a,f10.2,a,f10.2,a,e7.1,a,f4.2,a,f4.2,a,f4.2,a,f5.2)','tt=',tt, ' t0=',t0,' ds=',ds,' rr=',rr,' rbg=',rbg,' rb=',rb,' ts=',ts
!             print '(a,f8.1,a,f8.1,a,f12.2)','uflux(ia,ja,ka,1)=',uflux(ia,ja,ka,1),' uflux(ia ,ja,ka,2)=', &
!                 uflux(ia ,ja,ka,2),' dxyz=',dxyz
!             print '(a,f8.1,a,f8.1,a,f12.2)','uflux(ia-1,ja,ka,1)=',uflux(ia-1,ja,ka,1),' uflux(ia-1 ,ja,ka,2)=', &
!                 uflux(ia-1 ,ja,ka,2)
!             print '(a,f8.1,a,f8.1,a)','vflux(ia,ja,ka,1)=',vflux(ia,ja,ka,1),' vflux(ia ,ja,ka,2)=',vflux(ia ,ja,ka,2)
!             print '(a,f8.1,a,f8.1,a)','vflux(ia,ja-1,ka,1)=',vflux(ia,ja-1,ka,1),' vflux(ia ,ja-1,ka,2)=',vflux(ia ,ja-1,ka,2)
!             print '(a,f7.1,a,f6.1,a,f5.2,a,f5.2)', 'tt=',tt,' dt=',dt,' ts=',ts,' tss=',tss
!             print *,''
!         endif
!         === make sure that trajectory ===
        ! === is inside ib,jb,kb box    ===
!         print *,'ib=',ib,' jb=',jb,' x1=',x1,' idint(x1)=',idint(x1)
        if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! index for correct cell?
        if(y1.ne.dble(idint(y1))) jb=idint(y1)+1 ! index for correct cell?
!         print *,'ib=',ib,' jb=',jb

!             print *,'after box check'  
!             print '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)','x0=',x0,' x1=',x1,&
!                 ' y0=',y0,' y1=',y1,' z0=',z0,' z1=',z1
!             print '(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)','ia=',ia,' ib=',ib,&
!                 ' ja=',ja,' jb=',jb,' ka=',ka,' kb=',kb

!         call errorCheck('boundError', errCode)
        if(ia>imt .or. ib>imt .or. ja>jmt .or. jb>jmt &
             .or. ia<1 .or. ib<1 .or. ja<1 .or. jb<1) then
            print *,'====================================='
            print *,'Warning: Trajectory leaving model area'
            print *,'-------------------------------------'
    !           print *,'iaib',ia,ib,ja,jb,ka,kb
    !           print *,'xyz',x0,x1,y0,y1,z0,z1
    !           print *,'ds',dse,dsw,dsn,dss,dsu,dsd
    !           print *,'dsmin=',ds,dsmin,dtmin,dxyz
    !           print *,'tt=',tt,ts
    !           print *,'ntrac=',ntrac
            print *,'-------------------------------------'
            print *,'The trajectory is killed'
            print *,'====================================='
!            nerror=nerror+1
!            boundError = boundError +1
            errCode = -50
!            if (strict==1) stop
            flag(ntrac) = 1

!            nrj(ntrac,6)=1
        endif
        if (errCode.ne.0) cycle ntracLoop

!         call errorCheck('airborneError', errCode)
         ! if trajectory above sea level,
         ! then put back in the middle of shallowest layer (evaporation)
         if( z1.ge.dble(km) ) then
            z1=dble(km)-0.5d0
            kb=km
            errCode = -50
         endif

!         call errorCheck('corrdepthError', errCode)
         ! sets the right level for the corresponding trajectory depth
         if(z1.ne.dble(idint(z1))) then
            kb=idint(z1)+1
            if(kb == km+1) kb=km  ! (should perhaps be removed)
            errCode = -52
         endif

!         call errorCheck('cornerError', errCode)
         if(x1 == dble(idint(x1)) .and. y1 == dble(idint(y1)))  then
            print *,'corner problem'
            !print *,'corner problem',ntrac,x1,x0,y1,y0,ib,jb
            !print *,'ds=',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
            !stop 34957
            ! corner problems may be solved the following way 
            ! but should really not happen at all
            if(ds == dse .or. ds == dsw) then
                if(y1.ne.y0) then
                    y1=y0 ; jb=ja
                else
                    y1=dble(jb)-0.5d0
                endif
            else if(ds == dsn .or. ds == dss) then
                if(y1.ne.y0) then
                    x1=x0 ; ib=ia 
                else
                    x1=dble(ib)-0.5d0
                endif
            else
                x1=dble(ib)-0.5d0
                y1=dble(jb)-0.5d0
            endif
            errCode = -54
         endif


! print *,'ib(after)=',ib,' ia=',ia
! print *,'x1=',x1,' y1=',y1,' z1=',z1
        ! === diffusion, which adds a random position ===
        ! === position to the new trajectory          ===

        if(doturb==2 .or. doturb==3) then
            call diffuse(x1, y1, z1, ib, jb, kb, dt,imt,jmt,km, dxv,dyu,dzt,h,ah,av,do3d,doturb)
        endif
        !         nout=nout+1 ! number of trajectories that have exited the space and time domain

        ! === Optional periodic boundary conditions ===
        ! If activated, a drifter is moved from the appropriate domain edge and 
        ! wrapped around to the other side of the domain.
        ! === end trajectory if outside chosen domain ===
        ! Note that these are combined with checking for drifters having exited the
        ! domain since if using doperiodic in a direction, it cannot exit.
        ! Also note that I don't think drifters can actually have values outside the box [1,imt-1,1,jmt-1]
        if(x1<=1.d0) then ! at west end of domain
            if(doperiodic==1) then ! using periodic boundary conditions
                x1 = dble(imt-1) ! place drifter at the east end of the domain
                ib = idint(x1)+1 ! need to update grid index too
            else ! not using periodic boundary conditions
                x1 = 1.d0
                flag(ntrac) = 1
                exit niterLoop
            endif
        else if(x1>=imt-1) then ! at east end of domain
            if(doperiodic==1) then ! using periodic boundary conditions
                x1 = 1.d0 ! place drifter at the west end of the domain
                ib = idint(x1)+1 ! need to update grid index too
            else ! not using periodic boundary conditions
                x1 = dble(imt-1)
                flag(ntrac) = 1
                exit niterLoop
            endif
        else if(y1<=1.d0) then ! at south end of domain
            if(doperiodic==2) then ! using periodic boundary conditions
                y1 = dble(imt-1) ! place drifter at the north end of the domain
                jb = idint(y1)+1 ! need to update grid index too
            else ! not using periodic boundary conditions
                y1 = 1.d0
                flag(ntrac) = 1
                exit niterLoop
            endif
        else if(y1>=jmt-1) then ! at north end of domain
            if(doperiodic==2) then ! using periodic boundary conditions
                y1 = 1.d0 ! place drifter at the south end of the domain
                jb = idint(y1)+1 ! need to update grid index too
            else ! not using periodic boundary conditions
                y1 = dble(jmt-1)
                flag(ntrac) = 1
                exit niterLoop
            endif
        endif




!  print *,'ja=',ja,' jb=',jb,x1,y1

!             print *,'end of loop'  
!             print '(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)','x0=',x0,' x1=',x1,&
!                 ' y0=',y0,' y1=',y1,' z0=',z0,' z1=',z1
!             print '(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)','ia=',ia,' ib=',ib,&
!                 ' ja=',ja,' jb=',jb,' ka=',ka,' kb=',kb

        ! If drifter is on a grid cell wall, add/subtract its initial volume
        ! transport to the appropriate grid cells
        ! drifter needs to have just made it to the wall to count
!         print *,'x1=', x1, ' dble(ib)=', dble(ib), ' x0=', x0, ' dble(ia)=', dble(ia)
!         print *,'y1=', y1, ' dble(jb)=', dble(jb), ' y0=', y0, ' dble(ja)=', dble(ja)
        if (dostream==1) then
            if(idint(x1)>idint(x0)) then ! moving in positive x direction
!             if(x1==dble(ia) .and. x1>x0) then ! moving in positive x direction
!                 print *, 'moving in positive x direction'
!                 print *, 'ut(ia, jb)=', ut(ia, jb), &
!                             ' ut(ib, jb)=', ut(ib, jb), &
!                             ' T0(ntrac)=', T0(ntrac)
                ut(ia, ja) = ut(ia, ja) + T0(ntrac) ! positive direction exit
!                 ut(ia, jb) = ut(ia, jb) - T0(ntrac) ! leaving cell
!                 ut(ib, jb) = ut(ib, jb) + T0(ntrac) ! entering cell
!                 print *, 'ut(ia, jb)=', ut(ia, jb), &
!                             ' ut(ib, jb)=', ut(ib, jb)

            else if(idint(x1)<idint(x0)) then
!             else if(x1==dble(ia-1) .and. x1<x0) then ! moving in negative x direction
!                 print *, 'moving in negative x direction'
!                 print *, 'ut(ia, jb)=', ut(ia, jb), &
!                             ' ut(ib, jb)=', ut(ib, jb), &
!                             ' T0(ntrac)=', T0(ntrac)
                ut(ia-1, ja) = ut(ia-1, ja) - T0(ntrac) ! negative direction exit
!                 ut(ia, jb) = ut(ia, jb) - T0(ntrac) ! leaving cell
!                 ut(ib, jb) = ut(ib, jb) + T0(ntrac) ! entering cell
!                 print *, 'ut(ia, jb)=', ut(ia, jb), &
!                             ' ut(ib, jb)=', ut(ib, jb)
            else if(idint(y1)>idint(y0)) then
!             else if(y1==dble(ja) .and. y1>y0) then ! moving in positive y direction
!                 print *, 'moving in positive y direction'
!                 print *, 'vt(ib, ja)=', vt(ib, ja), &
!                             ' vt(ib, jb)=', vt(ib, jb), &
!                             ' T0(ntrac)=', T0(ntrac)
                vt(ia, ja) = vt(ia, ja) + T0(ntrac) ! positive direction exit
!                 vt(ib, ja) = vt(ib, ja) - T0(ntrac)
!                 vt(ib, jb) = vt(ib, jb) + T0(ntrac)
!                 print *, 'vt(ib, ja)=', vt(ib, ja), &
!                             ' vt(ib, jb)=', vt(ib, jb)
            else if(idint(y1)<idint(y0)) then
!             else if(y1==dble(ja-1) .and. y1<y0) then ! moving in negative y direction
!                 print *, 'moving in negative y direction'
!                 print *, 'vt(ib, ja)=', vt(ib, ja), &
!                             ' vt(ib, jb)=', vt(ib, jb), &
!                             ' T0(ntrac)=', T0(ntrac)
                vt(ia, ja-1) = vt(ia, ja-1) - T0(ntrac) ! negative direction exit
!                 vt(ib, ja) = vt(ib, ja) - T0(ntrac)
!                 vt(ib, jb) = vt(ib, jb) + T0(ntrac)
!                 print *, 'vt(ib, ja)=', vt(ib, ja), &
!                             ' vt(ib, jb)=', vt(ib, jb)
            end if
        end if

        ! Check for case in which the model iteration has passed a write time
       if(tt/tseas >= dble(Ni)/dble(N)) then

            ! Time interpolation weights
            ! Time at the output time minus the time for x0 over the time window, weight for x0, etc
            rwn = 1 - (((tseas/dble(N))*dble(Ni) - (tt-dt))/dt)
            rwp = 1 - rwn ! weight for x1, etc

            ! Interpolate to find values at the desired output time.
            ! linear combination of previous step and current step to find output
            ! at the output time between
            xend(ntrac,Ni) = rwn*x0 + rwp*x1 
            yend(ntrac,Ni) = rwn*y0 + rwp*y1
            zend(ntrac,Ni) = rwn*z0 + rwp*z1
!             iend(ntrac,Ni) = rwn*ia + rwp*ib
!             jend(ntrac,Ni) = rwn*ja + rwp*jb
!             kend(ntrac,Ni) = rwn*ka + rwp*kb
            ttend(ntrac,Ni) = (rwn*(tt-dt) + rwp*tt)*ff
            Ni = Ni + 1 ! counter for writing
        endif

! This is the normal stopping routine for the loop. I am going to do a shorter one
        ! === stop trajectory if the choosen time or ===
        ! === water mass properties are exceeded     ===
        ! Have already written to save arrays above
        if(tt.ge.tseas) then ! was .gt. in original code
!         if(tt-t0.ge.tseas) then ! was .gt. in original code, also eliminated t0 since was =0
        !               nexit(NEND)=nexit(NEND)+1
!             print *, 'Stopping trajectory due to time'
!             print *, 'tt=',tt,' t0=',t0,' tseas=',tseas
!             print *, 'x1=',x1,' y1=',y1
            flag(ntrac) = 0 ! want to continue drifters if time was the limiting factor here
            exit niterLoop
        endif
         
    end do niterLoop
    !         nout=nout+1 ! ' trajectories exited the space and time domain'

end do ntracLoop

end SUBROUTINE step



