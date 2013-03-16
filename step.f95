SUBROUTINE step(xstart,ystart,zstart,t0,istart,jstart,kstart,tseas,uflux,vflux,ff,IMT,JMT,KM,dxyzarray,ntractot,xend,yend,zend,flag)

!!--------------------
!! Loop to step a numerical drifter forward one dt.
!!
!! This is put together with code pieces from TRACMASS to be 
!! controlled by outside wrapping Python code.
!! Kristen M. Thyng Feb 2012
!! Code from:   loop.f95
!!
!! Options used from TRACMASS:
!!  * NOT regulardt, which I think is for stepping within model outputs, 
!!    since this will be accomplished in the wrapping Python.
!!  * NOT timeanalyt since it is not as robust as the original code.
!!  * I am including options for both means of adding in turbulence:
!!    turb (for turbuflux) and diffusion (or diffuse).
!!    Use only one of them, not both, if any.
!!  * Removed all writing calls.
!!  * NOT rerun flag
!!  * NOT tracer flag
!!  * ifs flag looks like it is for a specific example, ignoring for now
!!  * Vertical flags:
!!      * two dim flag: sets vertical velocity to zero and drifters are essentially drogued USE THIS
!!      * full_wflux flag: Use a full 3D wflux field, and calculate it within the code
!!      * explicit_w: Use a given vertical velocity
!!  * I have been commenting out the interpolation stuff
!!
!! Time variables:
!!   These are stored for between model outputs (whether original or interpolated)
!!   tt     time in seconds of the trajectory relative to the code start
!!   ts     model dataset time step of trajectory, which changes from 0 to 1 in the script since interpolation between model outputs is done in python
!!   tss
!!   dt     time step of trajectory iteration in seconds    
!!   ds     ds=dt/(dx*dy*dz) Units in (s/m**3)
!!   dtmin   was iterative time steping in seconds, now there is only one iteration between input flux fields
!!   dsmin  dtmin/(dx*dy*dz) Units in (s/m**3)
!!   tseas (sec) velocity fields time step = time between velocity fields
!! Interpolation constants:
!!   rr,rb,rbg      Time interpolation constant between 0 and 1 at different parts in the loop
!!   
!!
!! Input:
!!   xstart       :   Original x position of drifters
!!   ystart       :   Original y position of drifters
!!   zstart
!!   t0       :   Time at original position of drifters
!!   tseas    :   Time between velocity fields that are input
!!   uflux    :   u velocity flux field, two time steps [txkxjxi]
!!   vflux
!!   ff      :   Trajectory direction (forward=1, backward=-1
!!   hs        :    Sea surface height in meters [txjxi]?
!!  
!!
!! Output:
!!   xend       : the new position (coordinate), x
!!   yend       : the new position (coordinate), y
!!   zend       : the new position (coordinate), z
!!   flag       : set to 1 for a drifter if drifter shouldn't be stepped in the future anymore
!!  
!!---------------------

implicit none

integer,        intent(in)                                      :: ff, IMT, JMT, KM, ntractot
integer,        intent(out)                                      :: flag
real(kind=8),   intent(in),     dimension(ntractot)             :: xstart, ystart, zstart
real(kind=8),   intent(out),    dimension(ntractot)             :: xend,yend,zend
real(kind=8),   intent(in),     dimension(IMT,JMT-1,KM,2)         :: uflux
real(kind=8),   intent(in),     dimension(IMT-1,JMT,KM,2)         :: vflux
real(kind=8),   intent(in),     dimension(IMT-1,JMT-1,KM,2)         :: dxyzarray
! real(kind=4),   intent(in),     dimension(IMT,JMT,2)            :: hs
! real(kind=8),   intent(in),     dimension(KM)                   :: dz
! real(kind=8),   intent(in),     dimension(IMT,JMT)              :: dxdy
real(kind=8),   intent(in)                                      :: t0, tseas
integer,        intent(in),     dimension(ntractot)             :: istart, jstart, kstart
real(kind=8)                                                    :: rr, rg, ts, timax, &
                                                                    dstep, dtmin, dxyz, &
                                                                    dsmin, ds, dse,dsw,dsn,dss,dt, &
                                                                    x0,y0,z0,x1,y1,z1,tt,tss,rbg,rb, dsc
integer                                                         :: ntrac, niter, iter, ia,ja,ka,&
                                                                    iam,ib,jb,kb,kmt

! Number of drifters input in x0, y0
! ntractot = size(xstart) !Sum of x0 entries
! The following all seem to be for interpolation between output files in addition to the necessary iterations in the file.
!     ! controls the stepping
!     dstep=1.d0/dble(iter)
!     dtmin=dstep*tseas
dstep = 1.d0 ! The input model info is not being subdivided in this code.
dtmin = dstep*tseas 

timax = tseas! now in seconds *(1.d0/(24.d0*3600.d0)) ! maximum length of a trajectory in days. Just set to tseas.
iter = 1 ! one iteration since model output is interpolated before sending in here

kmt = KM ! Why are these separate? I think kmt is used in calc_dxyz as an array if they number of cells changes in x,y
flag = 0

!=======================================================
!=== Loop over all trajectories and calculate        ===
!=== a new position for this time step.              ===
!=======================================================
     
!      call fancyTimer('advection','start')
ntracLoop: do ntrac=1,ntractot  

    print *,'ntrac=',ntrac

    ! Counter for sub-interations for each drifter. In the source, this was read in but I am not sure why.
    niter = 0
    ! model dataset time step of trajectory. Initialize to the incoming time step for all drifters.
    ts = 0 !t1
    ! Initialize drifter parameters for this drifter (x0,y0,x0 will be updated from these at the start of the loop)
    x1 = xstart(ntrac)
    y1 = ystart(ntrac)
    z1 = zstart(ntrac)
    tt = t0 !+ 0.d0 !time of trajectory in seconds... start at zero for this test case
    ts = 0.d0 ! model dataset time step of trajectory... start at zero for this test case
    tss = 0.d0
!         subvol =  trj(ntrac,5)
!         t0     =  trj(ntrac,7)
    ib = istart(ntrac)
    jb = jstart(ntrac)
    kb = kstart(ntrac)
!        niter  =  nrj(ntrac,4)
!         ts     =  dble(nrj(ntrac,5))
!         tss    =  0.d0
  
    ! ===  start loop for each trajectory ===
    !         scrivi=.true.
    niterLoop: do 
        niter=niter+1 ! iterative step of trajectory
        print *,'niter=',niter
        rg=dmod(ts,1.d0) ! time interpolation constant between 0 and 1
        rr=1.d0-rg
!         print *,'ts=',ts,' rg=',rg,' rr=',rr
        ! Update particle indices and locations
        x0 = x1
        y0 = y1
        z0 = z1
        ia = ib
        iam = ia-1
        ja = jb
        ka = kb

        ! Are tracking fields being properly updated between loops?
! print *,'ib(before)=',ib,' ia=',ia
        call calc_dxyz(ib,jb,kb,rr,KM,kmt,dxyzarray,dxyz,JMT,IMT)

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
        ! #ifdef turb
        !   call turbuflux(ia,ja,ka,rr,dt)
        ! #endif /*turb*/

        ! WHAT IS RR VS RBG VS RG? THEY ARE ALL IN VARIOUES SUBROUTINES

        ! === calculate the vertical velocity ===
        !             call vertvel(rr,ia,iam,ja,ka) ! Start with 2D drifters
        ! supposed to be inputting nondimensional distances x0,y0 I think
!  print *,'ja=',ja,' jb=',jb,x1,y1
        call cross(1,ia,ja,ka,x0,dse,dsw,rr,uflux,vflux,ff,KM,JMT,IMT) ! zonal
        call cross(2,ia,ja,ka,y0,dsn,dss,rr,uflux,vflux,ff,KM,JMT,IMT) ! meridional
        !             call cross(3,ia,ja,ka,z0,dsu,dsd,rr) ! vertical ! Start with 2D drifters
        ds=dmin1(dse,dsw,dsn,dss,dsmin)
        print *, 'dse=',dse,' dsw=',dsw,' dsn=',dsn,' dss=',dss,' dsmin=',dsmin
        !             ds=dmin1(dse,dsw,dsn,dss,dsu,dsd,dsmin)
!         print *,'Before calc_time: ds=',ds,' tss=',tss,' ts=',ts,' tt=',tt,' dtmin=',dtmin
        call calc_time(ds,dsmin,dt,dtmin,tss,tseas,ts,tt,dxyz,dstep,iter,rbg,rb,dsc)

!         print *,'After calc_time: ds=',ds,' tss=',tss,' ts=',ts,' tt=',tt,' dt=',dt,' dtmin=',dtmin

        ! === calculate the new positions ===
        ! === of the trajectory           ===    
        call pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1,ds,dse,dsw,dsn,dss,dsmin,dsc,ff,IMT,JMT,KM,rr,rbg,rb,uflux,vflux)

        if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! index for correct cell?
        if(y1.ne.dble(idint(y1))) jb=idint(y1)+1 ! index for correct cell?

! print *,'ib(after)=',ib,' ia=',ia
! print *,'x1=',x1,' y1=',y1,' z1=',z1
        ! === diffusion, which adds a random position ===
        ! === position to the new trajectory          ===

        ! #if defined diffusion     
        !            call diffuse(x1,y1,z1,ib,jb,kb,dt)
        ! #endif
        !         nout=nout+1 ! number of trajectories that have exited the space and time domain

        ! === end trajectory if outside chosen domain ===

        ! For stream functions?
        !            LBTloop: do k=1,LBT
        !               if(float(ienw(k)) <= x1 .and. x1 <= float(iene(k)) .and. &
        !                  float(jens(k)) <= y1 .and. y1 <= float(jenn(k))  ) then
        !                  nexit(k)=nexit(k)+1
        !                  exit niterLoop                                
        !               endif
        !            enddo LBTLOOP
!  print *,'ja=',ja,' jb=',jb,x1,y1

        ! Need to add other conditions to this. Checking to see if drifter has exited domain.
        if(x1<=2.d0 .or. x1>=IMT-2.d0 .or. y1<=2.d0 .or. y1>=JMT-2.d0) then
            print *, 'Stopping trajectory due to domain'
            print *, 'x1=',x1,' y1=',y1
            if(x1<=2.d0) x1=2.d0
            if(x1>=IMT) x1=dble(IMT)
            if(y1<=2.d0) y1=2.d0
            if(y1>=JMT) y1=dble(JMT)
            xend(ntrac) = x1
            yend(ntrac) = y1
            zend(ntrac) = z1
            flag = 1
!             print *,'I should exit'
            exit niterLoop
!             print *,'I did not exit'
        endif

! This is the normal stopping routine for the loop. I am going to do a shorter one
        ! === stop trajectory if the choosen time or ===
        ! === water mass properties are exceeded     ===
        if(tt-t0.ge.timax) then ! was .gt. in original code
        !               nexit(NEND)=nexit(NEND)+1
            print *, 'Stopping trajectory due to time'
            print *, 'tt=',tt,' t0=',t0,' timax=',timax
            print *, 'x1=',x1,' y1=',y1
            xend(ntrac) = x1
            yend(ntrac) = y1
            zend(ntrac) = z1
            flag = 0 ! want to continue drifters if time was the limiting factor here
            exit niterLoop
        endif


!         print *, 'x0[',niter,']=', x0, ' y0[',niter,']=',y0
!         print *, 'x1[',niter,']=', x1, ' y1[',niter,']=',y1
!         print *,'ts=',ts,' rg=',rg,' rr=',rr
!     print *,'tt=',tt,' t0=',t0,' timax=',timax

! ! this is for shortening the loop for troubleshooting
!         if(niter .gt. 0) then
!             xend(ntrac) = x1
!             yend(ntrac) = y1
!             zend(ntrac) = z1
!             exit niterLoop
!         endif



         
    end do niterLoop
    !         nout=nout+1 ! ' trajectories exited the space and time domain'

end do ntracLoop

end SUBROUTINE step



