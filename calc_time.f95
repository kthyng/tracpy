subroutine calc_time(ds,dsmin,dt,dtmin,tss,tseas,ts,tt,dxyz,dstep,iter,rbg,rb,dsc)


! Inputs/output list here isn't correct
! Inputs:
!  ds
!  dsmin    Iterative time stepping in s/m^3
!  dtmin    Iterative time stepping in seconds
!  tseas    Time between model outputs (s)
!  dxyz
!  dstep    Model time step between time interpolation steps between model outputs
!  iter     Number of interpolations to do between model outputs
!
! Outputs:
!  dt       Timestep of trajectory iteration (in seconds)
!  tt       time in seconds of the trajectory relative to the code start
!  ts       Time step of time interpolation between model outputs, non-dimensional, has values between 0 and 1
!  tss      Counter for iterations between model outputs. Counts up to iter I think.
!  rbg
!  rb
!  dsc


    IMPLICIT NONE

real(kind=8), intent(in) :: dsmin,dtmin,tseas,dxyz,dstep
real(kind=8), intent(in out) :: ds,tss,tt,ts
integer, intent(in)   :: iter
real(kind=8), intent(out) :: dt,rbg,rb,dsc

    ! Don't allow particle to move more than between the model outputs
    if(ds == dsmin) then ! transform ds to dt in seconds
        dt=dtmin  ! this makes dt more accurate
!     print *,'ds=',ds,' dsmin=',dsmin
    else
        dt=ds*dxyz 
    endif
    if(dt.lt.0.d0) then
        print *,'dt=',dt
        stop 49673
    endif
    ! === if time step makes the integration ===
    ! === exceed the time when fiedls change ===
    ! This checks if the previous iteration count, tss, plus the new increment of 
    ! iteration is larger than the total number of iterations, in which case 
    ! everything is set to its value at the finish of all the iterations
    if(tss+dt/tseas*dble(iter).ge.dble(iter)) then 
!         print *,'dt=',dt,' tseas=',tseas,' dtmin=',dtmin,' ds=',ds,' dxyz=',dxyz
!         print *,'original if in calc_time'
        dt=dble(idint(ts)+1)*tseas-tt
        tt=dble(idint(ts)+1)*tseas
        ts=dble(idint(ts)+1) ! total number of time steps, but 1 is the biggest if my output looping is outside tracmass
        tss=dble(iter) ! total number of iterations for iteration counter
        ds=dt/dxyz
        dsc=ds
    ! KMT adding a condition to stop at the time interpolation between model outputs
    ! trying to find when the trajectory timing is switching from below to above an
    ! interpolation step in order to stop at time interpolation step and later write
    elseif(dble(idint(tss))<dble(idint(tss+dt/tseas*dble(iter)))) then
!         print *,'new elseif in calc_time'
        tss=dble(idint(tss)+1) ! time interpolation step (should be whole number)
        ts=tss/dble(iter) ! fractional number of time steps, but 1 is the biggest if my output looping is outside tracmass
        dt=ts*tseas-tt
        tt=ts*tseas
        ds=dt/dxyz
        dsc=ds
    else
        tt=tt+dt
        if(dt == dtmin) then ! If the particle is moving the full time of the outputs, step it accordingly
           ts=ts+dstep ! step to the next interpolated step between model outputs
           tss=tss+1.d0 ! add 1 iteration on to iteration counter tss
        else ! Otherwise, step it within the outputs
           ts =ts +dt/tseas ! non-dimensional, add in normalized incremental time step
           tss=tss+dt/tseas*dble(iter) ! add incremental amount of iteration onto iteration counter tss
        !                 tss=tss+dt/dtmin
        endif
    end if
    ! === time interpolation constant ===
    ! KMT change: using the dmod function makes the final rb,rbg values be switched
    ! in value, so rb=1, rbg=0 when it should be the opposite at the end of a model time step
    ! However, I am not sure why this would be wrong here, so I want to ask in the future.
    rbg=ts/1.d0 
!     rbg=dmod(ts,1.d0) 
    rb =1.d0-rbg

end subroutine calc_time