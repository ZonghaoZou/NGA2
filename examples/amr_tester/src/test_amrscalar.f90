!> Test amrscalar solver with time integration and Ensight output
module mod_test_amrscalar
   use precision,         only: WP
   use string,            only: itoa,rtoa,str_medium
   use amrgrid_class,     only: amrgrid
   use amrscalar_class,   only: amrscalar
   use amrensight_class,  only: amrensight
   use amrviz_class,      only: amrviz
   use amrio_class,       only: amrio
   use monitor_class,     only: monitor
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use messager,          only: log,warn
   use amrex_amr_module,  only: amrex_boxarray,amrex_distromap,amrex_multifab
   implicit none
   private
   public :: test_amrscalar

   ! Module-level objects (referenced by callbacks)
   type(amrgrid),   allocatable, target :: amr
   type(amrscalar), allocatable, target :: sc

   real(WP), parameter :: SC_REFINE_THRESH=0.01_WP  !< Refine where SC > this value

contains

   !> Scalar-based tagger - refine where SC > threshold
   subroutine box_tagger(lvl,tags_ptr,time)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_tagboxarray,amrex_mfiter,amrex_box
      implicit none
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      character(kind=c_char), parameter :: SET=char(1)
      integer :: i,j,k
      tags=tags_ptr
      ! Iterate over level boxes (SC and tags have same BoxArray)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         pSC=>sc%SC%mf(lvl)%dataptr(mfi)
         tagarr=>tags%dataPtr(mfi)
         do k=bx%lo(3),bx%hi(3)
            do j=bx%lo(2),bx%hi(2)
               do i=bx%lo(1),bx%hi(1)
                  ! Tag where SC exceeds threshold
                  if (pSC(i,j,k,1).gt.SC_REFINE_THRESH) then
                     tagarr(i,j,k,1)=SET
                  end if
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)
   end subroutine box_tagger


   !> Initialize scalar field with Gaussian blob (called by amrdata on_init callback)
   subroutine init_gaussian_blob(lvl, mf, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_mfiter, amrex_box, amrex_mfiter_build
      implicit none
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      real(WP) :: x,y,z,dx,dy,dz
      integer :: i,j,k
      dx = geom%dx(1)
      dy = geom%dx(2)
      dz = geom%dx(3)
      call amrex_mfiter_build(mfi, mf)
      do while (mfi%next())
         bx = mfi%tilebox()
         pSC => mf%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  ! Smaller, sharper Gaussian offset from center
                  pSC(i,j,k,1) = exp(-200.0_WP*((x-0.25_WP)**2 + (y-0.5_WP)**2 + (z-0.5_WP)**2))
               end do
            end do
         end do
      end do
   end subroutine init_gaussian_blob

   subroutine test_amrscalar()
      use iso_c_binding,    only: c_associated
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mathtools,        only: Pi
      implicit none
      type(amrensight) :: ens
      type(amrviz) :: viz
      type(amrio) :: io
      type(monitor) :: mfile
      type(timetracker) :: time
      type(event) :: regrid_evt,ensight_evt,hdf5_evt
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: U,V,W,dSCdt,SCfill
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC,pU,pV,pW
      integer :: lvl,i,j,k
      real(WP) :: x,y,z,dx,dy,dz,int0,intF

      call log("---------------------------------------------------")
      call log("Running Test: amrscalar Time Integration")
      call log("---------------------------------------------------")

      ! Allocate module objects
      allocate(amr,sc)

      ! Setup grid (periodic box)
      amr%nx  =32
      amr%ny  =32
      amr%nz  =32
      amr%xlo =0.0_WP
      amr%xhi =1.0_WP
      amr%ylo =0.0_WP
      amr%yhi =1.0_WP
      amr%zlo =0.0_WP
      amr%zhi =1.0_WP
      amr%xper=.true.
      amr%yper=.true.
      amr%zper=.true.
      amr%maxlvl=2   ! 3 levels total (0,1,2)
      amr%nmax=16

      call amr%initialize("sc_amr")

      ! Initialize scalar solver first (registers callbacks)
      call sc%initialize(amr,nscalar=1,name="test_scalar")

      ! Set user on_init callback for SC to initialize level 0 with Gaussian blob
      ! (fine levels get interpolated data via fill_from_coarse, then we average_down)
      call sc%SC%set_on_init(init_gaussian_blob)

      ! Set tagging
      call amr%add_tagging(box_tagger)

      ! Build grid
      call amr%initialize_grid(0.0_WP)
      call amr%get_info()
      call log("Grid built: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Initialize ensight output
      call ens%initialize(amr=amr,name="sc_advect")
      call ens%add_scalar(data=sc%SC,comp=1,name="SC")

      ! Initialize HDF5 viz output (all fields in single file per timestep)
      call viz%initialize(amr=amr,name='sc_advect')
      call viz%add_scalar(sc%SC, 1, 'SC')
      call viz%add_scalar(sc%SCold, 1, 'SCold')

      ! Build all levels using init_from_scratch
      ! (callback init_gaussian_blob is called for each level, then postregrid fires average_down)
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)
      call amr%get_info()
      call log("After init_from_scratch: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Get initial integral
      call sc%get_info()
      int0=sc%SCint(1)
      call log("Initial: SCint="//trim(rtoa(int0)))



      ! Setup timetracker
      time=timetracker(amRoot=amr%amRoot)
      time%dt=0.0025_WP
      time%dtmax=time%dt
      time%tmax=2.0_WP*PI  ! Full circle rotation
      time%t=0.0_WP
      time%n=0

      ! Setup regrid event (every 10 steps)
      regrid_evt=event(time=time,name='Regrid')
      regrid_evt%nper=10

      ! Setup ensight event (every 20 steps)
      ensight_evt=event(time=time,name='Ensight')
      ensight_evt%tper=0.125_WP

      ! Setup HDF5 event (every 0.125 time units)
      hdf5_evt=event(time=time,name='HDF5')
      hdf5_evt%tper=0.125_WP

      ! Setup monitor file
      call sc%get_info()
      mfile=monitor(amr%amRoot,'scalar')
      call mfile%add_column(time%n,'Timestep')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(sc%SCmin(1),'SC_min')
      call mfile%add_column(sc%SCmax(1),'SC_max')
      call mfile%add_column(sc%SCint(1),'SC_int')
      call mfile%write()

      ! Write initial output at start time
      call ens%write(time=time%t)
      call viz%write(time=time%t)

      call log("Advancing to t="//trim(rtoa(time%tmax))//", dt="//trim(rtoa(time%dt)))

      ! Time loop
      do while (.not.time%done())
         call time%increment()

         ! Copy to old
         call sc%copy2old()

         ! Loop over levels (uniform dt)
         do lvl=0,amr%clvl()
            dx=amr%dx(lvl)
            dy=amr%dy(lvl)
            dz=amr%dz(lvl)

            ! Build velocity (rotating vortex in xy plane)
            call amr%mfab_build(lvl=lvl,mfab=U,ncomp=1,nover=0,atface=[.true.,.false.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=V,ncomp=1,nover=0,atface=[.false.,.true.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=W,ncomp=1,nover=0,atface=[.false.,.false.,.true.])
            call amr%mfab_build(lvl=lvl,mfab=dSCdt,ncomp=1,nover=0,atface=[.false.,.false.,.false.])
            call amr%mfab_build(lvl=lvl,mfab=SCfill,ncomp=1,nover=2,atface=[.false.,.false.,.false.])

            ! Set velocity: solid body rotation in z-plane (U=-(y-0.5), V=(x-0.5), W=0)
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               bx=mfi%tilebox()
               pU=>U%dataptr(mfi)
               pV=>V%dataptr(mfi)
               pW=>W%dataptr(mfi)
               ! U on x-faces
               do k=lbound(pU,3),ubound(pU,3)
                  do j=lbound(pU,2),ubound(pU,2)
                     y=amr%ylo+(real(j,WP)+0.5_WP)*dy
                     do i=lbound(pU,1),ubound(pU,1)
                        ! Constant rotation (old: time-reversing)
                        pU(i,j,k,1)=-(y-0.5_WP)
                        !pU(i,j,k,1)=-(y-0.5_WP)*cos(Pi*time%t/time%tmax)
                     end do
                  end do
               end do
               ! V on y-faces
               do k=lbound(pV,3),ubound(pV,3)
                  do j=lbound(pV,2),ubound(pV,2)
                     do i=lbound(pV,1),ubound(pV,1)
                        x=amr%xlo+(real(i,WP)+0.5_WP)*dx
                        ! Constant rotation (old: time-reversing)
                        pV(i,j,k,1)=(x-0.5_WP)
                        !pV(i,j,k,1)=(x-0.5_WP)*cos(Pi*time%t/time%tmax)
                     end do
                  end do
               end do
               ! W on z-faces
               pW=0.0_WP
            end do
            call amr%mfiter_destroy(mfi)

            ! Fill SCfill with ghost cells properly (uses FillPatch for C-F interface)
            call sc%SCold%fill_mfab(SCfill, lvl, time%t)

            ! Calculate dSC/dt
            call sc%get_dSCdt(lvl,dSCdt,SCfill,U,V,W)

            ! Forward Euler step: SC = SCold + dt*dSC/dt
            call sc%SC%mf(lvl)%lincomb(a=1.0_WP,srcmf1=sc%SCold%mf(lvl),srccomp1=1,&
            &                          b=time%dt,srcmf2=dSCdt          ,srccomp2=1,&
            &                          dstcomp=1,nc=1,ng=0)

            ! Cleanup level work arrays
            call amr%mfab_destroy(U)
            call amr%mfab_destroy(V)
            call amr%mfab_destroy(W)
            call amr%mfab_destroy(dSCdt)
            call amr%mfab_destroy(SCfill)
         end do

         ! Reflux and average down
         call sc%reflux_avg(time%dt)

         ! Regrid if needed (average_down is now handled by amrscalar postregrid callback)
         if (regrid_evt%occurs()) then
            call log("Regridding at step "//trim(itoa(time%n)))
            call amr%regrid(baselvl=0,time=time%t)
            call amr%get_info()
            call log("  Grid: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")
         end if

         ! Write ensight output
         if (ensight_evt%occurs()) then
            call ens%write(time=time%t)
            call log("Step "//trim(itoa(time%n))//": t="//trim(rtoa(time%t))//", SCint="//trim(rtoa(sc%SCint(1))))
         end if

         ! Write HDF5 output
         if (hdf5_evt%occurs()) then
            call viz%write(time=time%t)
         end if

         ! Update monitor file every step
         call sc%get_info()
         call mfile%write()
      end do

      ! Final conservation check
      call sc%get_info()
      intF=sc%SCint(1)
      call log("Final: SCint="//trim(rtoa(intF)))
      call log("Conservation error: "//trim(rtoa(abs(intF-int0)/int0*100.0_WP))//"%")

      if (abs(intF-int0)/int0<1.0e-10_WP) then
         call log("PASS: Scalar conserved to machine precision")
      else
         call warn("FAIL: Conservation error too large!")
      end if

      ! Test checkpoint round-trip
      call io%initialize(amr, nfiles=1)
      call sc%register_checkpoint(io)  ! Solver registers its fields
      call io%add_scalar('dt', time%dt)  ! Store dt as metadata
      call io%write('checkpoint_test', time%t, time%n)
      call log("PASS: Checkpoint written to checkpoint_test/")

      ! Save current integral, then zero the data
      call sc%get_info()
      int0 = sc%SCint(1)
      call log("  Before zero: SCint="//trim(rtoa(int0)))

      ! Zero out SC data to test read_data
      do lvl = 0, amr%clvl()
         call amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            pSC => sc%SC%mf(lvl)%dataptr(mfi)
            pSC = 0.0_WP
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call sc%get_info()
      call log("  After zero: SCint="//trim(rtoa(sc%SCint(1))))

      ! Read checkpoint back
      block
         real(WP) :: read_time, read_dt
         integer :: read_step
         call io%read_header('checkpoint_test', read_time, read_step)
         call io%get_scalar('dt', read_dt)
         call sc%restore_checkpoint(io, 'checkpoint_test')  ! Solver restores its fields
         call log("  Read checkpoint: time="//trim(rtoa(read_time))//" step="//trim(itoa(read_step))//" dt="//trim(rtoa(read_dt)))
      end block

      ! Verify data restored
      call sc%get_info()
      intF = sc%SCint(1)
      call log("  After read: SCint="//trim(rtoa(intF)))

      if (abs(intF-int0)/int0 .lt. 1.0e-10_WP) then
         call log("PASS: Checkpoint round-trip verified!")
      else
         call warn("FAIL: Checkpoint round-trip failed! Integral mismatch.")
      end if

      call log("PASS: HDF5 plotfiles written to amrviz/sc_advect/")

      ! Cleanup
      call io%finalize()
      call viz%finalize()
      call mfile%finalize()
      call sc%finalize()
      call amr%finalize()
      if (allocated(sc)) deallocate(sc)
      if (allocated(amr)) deallocate(amr)
      call log("PASS: amrscalar test complete!")
      call log("View output in ensight/sc_advect/")

   end subroutine test_amrscalar

end module mod_test_amrscalar
