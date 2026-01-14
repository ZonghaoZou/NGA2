!> AMR Visualization class: HDF5 plotfile output with field registration
!> Outputs Chombo-compatible HDF5 files with all registered fields in single file
!> readable by VisIt (Chombo reader) and ParaView (VisItChomboReader)
module amrviz_class
   use precision,      only: WP
   use string,         only: str_medium
   use iso_c_binding,  only: c_ptr,c_char,c_int,c_double,c_null_char,c_loc
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use amrex_interface, only: amrplotfile_write_hdf5,amrplotfile_read_time
   use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box, &
   &                           amrex_multifab_build,amrex_multifab_destroy
   implicit none
   private

   public :: amrviz

   ! Scalar field registration node
   type :: scl
      type(scl), pointer :: next => null()
      character(len=str_medium) :: name
      type(amrdata), pointer :: ptr => null()
      integer :: comp
   end type scl

   ! Vector field registration node
   type :: vct
      type(vct), pointer :: next => null()
      character(len=str_medium) :: name
      type(amrdata), pointer :: ptrx => null()
      type(amrdata), pointer :: ptry => null()
      type(amrdata), pointer :: ptrz => null()
      integer :: compx, compy, compz
   end type vct

   !> AMR visualization handler
   type :: amrviz
      class(amrgrid), pointer :: amr => null()   !< Pointer to AMR grid
      character(len=str_medium) :: name          !< Subdirectory name for output
      type(scl), pointer :: first_scl => null()  !< Registered scalars
      type(vct), pointer :: first_vct => null()  !< Registered vectors
      integer :: ntime = 0                       !< File counter for output
      real(WP), allocatable :: time(:)           !< Time values for each file
   contains
      procedure :: initialize                    !< Initialize with grid
      procedure :: add_scalar                    !< Register a scalar field
      procedure :: add_vector                    !< Register a vector field
      procedure :: write                         !< Write all registered fields to single HDF5
      procedure :: finalize                      !< Clean up registered field lists
   end type amrviz

contains

   !> Initialize visualization handler with AMR grid
   !> If output directory already exists with files, reads their times to continue series
   subroutine initialize(this, amr, name)
      use filesys, only: makedir, isdir
      use mpi_f08, only: MPI_BCAST, MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrviz), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name

      character(len=str_medium) :: filename, dirname
      integer :: ierr, n
      real(WP) :: file_time

      this%amr => amr
      this%name = trim(adjustl(name))
      this%ntime = 0
      this%first_scl => null()
      this%first_vct => null()

      ! Create output directory: amrviz/<name>/
      dirname = 'amrviz/'//trim(this%name)
      if (this%amr%amRoot) then
         if (.not.isdir('amrviz')) call makedir('amrviz')
         if (.not.isdir(trim(dirname))) call makedir(trim(dirname))
      end if

      ! Check for existing files and read their times (root only, then broadcast)
      if (this%amr%amRoot) then
         ! First pass: count existing files
         n = 0
         do
            n = n + 1
            write(filename,'(a,"/nga2amr.",i6.6,".h5")') trim(dirname), n
            file_time = amrplotfile_read_time(trim(filename)//c_null_char)
            if (file_time .lt. 0.0_WP) exit  ! File doesn't exist
         end do
         this%ntime = n - 1  ! Last existing file number

         ! Second pass: allocate and fill time array
         if (this%ntime .gt. 0) then
            allocate(this%time(this%ntime))
            do n = 1, this%ntime
               write(filename,'(a,"/nga2amr.",i6.6,".h5")') trim(dirname), n
               this%time(n) = amrplotfile_read_time(trim(filename)//c_null_char)
            end do
         end if
      end if

      ! Broadcast ntime and time array to all ranks
      call MPI_BCAST(this%ntime, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
      if (this%ntime .gt. 0) then
         if (.not.this%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time, this%ntime, MPI_REAL_WP, 0, this%amr%comm, ierr)
      end if
   end subroutine initialize


   !> Register a scalar field for output
   subroutine add_scalar(this, data, comp, name)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      type(amrdata), target, intent(in) :: data
      integer, intent(in) :: comp
      character(len=*), intent(in), optional :: name
      type(scl), pointer :: new_scl

      ! Check that the component is meaningful
      if (comp.le.0) call die('[amrviz add_scalar] comp must be at least one')
      if (comp.gt.data%ncomp) call die('[amrviz add_scalar] comp is too large for provided amrdata')

      ! Create new scalar node
      allocate(new_scl)
      if (present(name)) then
         new_scl%name = trim(adjustl(name))
      else
         new_scl%name = trim(adjustl(data%name))
      end if
      new_scl%ptr => data
      new_scl%comp = comp

      ! Insert at front of list
      new_scl%next => this%first_scl
      this%first_scl => new_scl
   end subroutine add_scalar


   !> Register a vector field for output (3 components from 3 amrdata objects)
   subroutine add_vector(this, datax, compx, datay, compy, dataz, compz, name)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      type(amrdata), target, intent(in) :: datax, datay, dataz
      integer, intent(in) :: compx, compy, compz
      character(len=*), intent(in), optional :: name
      type(vct), pointer :: new_vct

      ! Check that the components are meaningful
      if (compx.le.0) call die('[amrviz add_vector] compx must be at least one')
      if (compx.gt.datax%ncomp) call die('[amrviz add_vector] compx is too large for provided datax')
      if (compy.le.0) call die('[amrviz add_vector] compy must be at least one')
      if (compy.gt.datay%ncomp) call die('[amrviz add_vector] compy is too large for provided datay')
      if (compz.le.0) call die('[amrviz add_vector] compz must be at least one')
      if (compz.gt.dataz%ncomp) call die('[amrviz add_vector] compz is too large for provided dataz')

      ! Create new vector node
      allocate(new_vct)
      if (present(name)) then
         new_vct%name = trim(adjustl(name))
      else
         new_vct%name = 'vector'
      end if
      new_vct%ptrx => datax
      new_vct%ptry => datay
      new_vct%ptrz => dataz
      new_vct%compx = compx
      new_vct%compy = compy
      new_vct%compz = compz

      ! Insert at front of list
      new_vct%next => this%first_vct
      this%first_vct => new_vct
   end subroutine add_vector


   !> Write all registered fields to a single HDF5 plotfile
   !> File: amrviz/<name>/nga2amr.<ntime>.h5
   !> Uses time-based insertion: finds correct position based on physical time
   subroutine write(this, time)
      use messager, only: die
      implicit none
      class(amrviz), intent(inout) :: this
      real(WP), intent(in) :: time

      integer :: nlev, lev, nscl, nvct, ncomp, icomp, i, n
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      real(WP), dimension(:), allocatable :: temp_time
      character(len=str_medium) :: ctime, filename
      real(WP) :: rtime

      ! Temporary combined MultiFab per level
      type(amrex_multifab), allocatable :: combined(:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: dst, src

      ! Arrays for HDF5 writer
      type(c_ptr), allocatable :: mf_ptrs(:), geom_ptrs(:)
      integer(c_int), allocatable :: level_steps(:), ref_ratios(:)
      character(len=64), target, allocatable :: varnames(:)
      type(c_ptr), allocatable :: varname_ptrs(:)

      nlev = this%amr%maxlvl + 1

      ! Time-based insertion with tolerance for floating point
      if (this%ntime.eq.0) then
         ! First time stamp
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         ! There are time stamps already, check where to insert
         n=1
         rewind: do i=this%ntime,1,-1
            ! Safe .lt. check with tolerance for floating point drift
            if (this%time(i).lt.time-1.0e-6_WP) then
               n=i+1; exit rewind
            end if
         end do rewind
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if

      ! Count registered fields
      nscl = 0; my_scl => this%first_scl
      do while (associated(my_scl))
         nscl = nscl + 1
         my_scl => my_scl%next
      end do
      nvct = 0; my_vct => this%first_vct
      do while (associated(my_vct))
         nvct = nvct + 1
         my_vct => my_vct%next
      end do

      ncomp = nscl + 3*nvct
      if (ncomp == 0) return  ! Nothing to write

      ! Allocate combined MultiFab for each level
      allocate(combined(0:this%amr%maxlvl))

      ! Build variable names array
      allocate(varnames(ncomp))
      allocate(varname_ptrs(ncomp))
      icomp = 0

      ! Scalars first
      my_scl => this%first_scl
      do while (associated(my_scl))
         icomp = icomp + 1
         varnames(icomp) = trim(my_scl%name)//c_null_char
         varname_ptrs(icomp) = c_loc(varnames(icomp))
         my_scl => my_scl%next
      end do

      ! Then vectors (x,y,z components)
      my_vct => this%first_vct
      do while (associated(my_vct))
         icomp = icomp + 1
         varnames(icomp) = trim(my_vct%name)//'_x'//c_null_char
         varname_ptrs(icomp) = c_loc(varnames(icomp))
         icomp = icomp + 1
         varnames(icomp) = trim(my_vct%name)//'_y'//c_null_char
         varname_ptrs(icomp) = c_loc(varnames(icomp))
         icomp = icomp + 1
         varnames(icomp) = trim(my_vct%name)//'_z'//c_null_char
         varname_ptrs(icomp) = c_loc(varnames(icomp))
         my_vct => my_vct%next
      end do

      ! Build combined MultiFab for each level and copy data
      do lev = 0, this%amr%maxlvl
         ! Get a reference amrdata to use for building the combined mfab
         my_scl => this%first_scl
         if (.not.associated(my_scl)) call die('[amrviz write] No scalars registered')

         ! Build combined MultiFab with same BoxArray/DistroMap as first scalar
         call amrex_multifab_build(combined(lev), my_scl%ptr%mf(lev)%ba, my_scl%ptr%mf(lev)%dm, ncomp, 0)

         ! Copy scalar data into combined MultiFab
         icomp = 0
         my_scl => this%first_scl
         do while (associated(my_scl))
            icomp = icomp + 1
            ! Use mfiter to copy component by component
            call this%amr%mfiter_build(lev, mfi)
            do while (mfi%next())
               bx = mfi%tilebox()
               src => my_scl%ptr%mf(lev)%dataptr(mfi)
               dst => combined(lev)%dataptr(mfi)
               dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                  src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_scl%comp)
            end do
            call this%amr%mfiter_destroy(mfi)
            my_scl => my_scl%next
         end do

         ! Copy vector data into combined MultiFab
         my_vct => this%first_vct
         do while (associated(my_vct))
            ! X component
            icomp = icomp + 1
            call this%amr%mfiter_build(lev, mfi)
            do while (mfi%next())
               bx = mfi%tilebox()
               src => my_vct%ptrx%mf(lev)%dataptr(mfi)
               dst => combined(lev)%dataptr(mfi)
               dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                  src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_vct%compx)
            end do
            call this%amr%mfiter_destroy(mfi)
            ! Y component
            icomp = icomp + 1
            call this%amr%mfiter_build(lev, mfi)
            do while (mfi%next())
               bx = mfi%tilebox()
               src => my_vct%ptry%mf(lev)%dataptr(mfi)
               dst => combined(lev)%dataptr(mfi)
               dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                  src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_vct%compy)
            end do
            call this%amr%mfiter_destroy(mfi)
            ! Z component
            icomp = icomp + 1
            call this%amr%mfiter_build(lev, mfi)
            do while (mfi%next())
               bx = mfi%tilebox()
               src => my_vct%ptrz%mf(lev)%dataptr(mfi)
               dst => combined(lev)%dataptr(mfi)
               dst(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), icomp) = &
                  src(bx%lo(1):bx%hi(1), bx%lo(2):bx%hi(2), bx%lo(3):bx%hi(3), my_vct%compz)
            end do
            call this%amr%mfiter_destroy(mfi)
            my_vct => my_vct%next
         end do
      end do

      ! Generate filename: amrviz/<name>/nga2amr.<ntime>
      filename = 'amrviz/'//trim(this%name)//'/nga2amr.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime

      ! Prepare pointers for HDF5 writer
      allocate(mf_ptrs(0:this%amr%maxlvl))
      allocate(geom_ptrs(0:this%amr%maxlvl))
      allocate(level_steps(0:this%amr%maxlvl))
      allocate(ref_ratios(0:max(this%amr%maxlvl-1,0)))

      do lev = 0, this%amr%maxlvl
         mf_ptrs(lev) = combined(lev)%p
         geom_ptrs(lev) = this%amr%geom(lev)%p
         level_steps(lev) = this%ntime
      end do
      do lev = 0, this%amr%maxlvl-1
         ref_ratios(lev) = this%amr%rref(lev)
      end do

      ! Write HDF5 file
      call amrplotfile_write_hdf5(trim(filename)//c_null_char, nlev, mf_ptrs, &
         varname_ptrs, ncomp, geom_ptrs, real(time,c_double), level_steps, &
         ref_ratios, c_null_char)

      ! Cleanup combined MultiFabs
      do lev = 0, this%amr%maxlvl
         call amrex_multifab_destroy(combined(lev))
      end do
      deallocate(combined)
      deallocate(mf_ptrs, geom_ptrs, level_steps, ref_ratios)
      deallocate(varnames, varname_ptrs)

   end subroutine write


   !> Finalize: clean up registered field lists
   subroutine finalize(this)
      implicit none
      class(amrviz), intent(inout) :: this
      type(scl), pointer :: cur_scl, next_scl
      type(vct), pointer :: cur_vct, next_vct

      ! Clean up scalar list
      cur_scl => this%first_scl
      do while (associated(cur_scl))
         next_scl => cur_scl%next
         deallocate(cur_scl)
         cur_scl => next_scl
      end do
      nullify(this%first_scl)

      ! Clean up vector list
      cur_vct => this%first_vct
      do while (associated(cur_vct))
         next_vct => cur_vct%next
         deallocate(cur_vct)
         cur_vct => next_vct
      end do
      nullify(this%first_vct)

      ! Reset state
      nullify(this%amr)
      this%ntime = 0
      if (allocated(this%time)) deallocate(this%time)
   end subroutine finalize

end module amrviz_class
