subroutine calc_time(ds,dsmin,dt,dtmin,tss,tseas,ts,tt,dxyz,dstep,iter,rbg,rb,dsc)


! Inputs/output list here isn't correct
! Inputs:
!  ds
!  dsmin
!  dtmin
!  tseas
!  dxyz
!  dstep
!  iter
!
! Outputs:
!  dt
!  tt     
!  ts
!  tss
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
    if(tss+dt/tseas*dble(iter).ge.dble(iter)) then
!         print *,'dt=',dt,' tseas=',tseas,' dtmin=',dtmin,' ds=',ds,' dxyz=',dxyz
        dt=dble(idint(ts)+1)*tseas-tt
        tt=dble(idint(ts)+1)*tseas
        ts=dble(idint(ts)+1)
        tss=dble(iter)
        ds=dt/dxyz
        dsc=ds
    else
        tt=tt+dt
        if(dt == dtmin) then ! If the particle is moving the full time of the outputs, step it accordingly
           ts=ts+dstep
           tss=tss+1.d0
        else ! Otherwise, step it within the outputs
           ts =ts +dt/tseas
           tss=tss+dt/tseas*dble(iter)
        !                 tss=tss+dt/dtmin
        endif
    end if
    ! === time interpolation constant ===
    rbg=dmod(ts,1.d0) 
    rb =1.d0-rbg

end subroutine calc_time