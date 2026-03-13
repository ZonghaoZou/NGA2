!> Definition for a postproc class
module postproc_class
   use precision,         only: WP,SP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use timetracker_class, only: timetracker
   use cclabel_class,     only: cclabel
   use lpt_class,         only: lpt
   use tpns_class,        only: tpns
   use sgsmodel_class,    only: sgsmodel
   implicit none
   private
   
   public :: postproc
   
   !> postproc object
   type :: postproc
      !> Input file for the simulation
      type(inputfile) :: input
      !> Config
      type(config) :: cfg
      !> Time info
      type(timetracker) :: time  
      !> Ensight postprocessing
      type(partmesh) :: pmesh    !< Particle mesh for core output
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(cclabel)  :: ccl
      type(lpt)      :: lp       !< Tracking particles
      type(tpns)     :: fs    !< Two-phase flow solver
      type(sgsmodel) :: sgs   !< SGS model for eddy viscosity
      !> Data arrays
      real(WP), dimension(:,:,:), allocatable :: VF
      real(WP), dimension(:), allocatable :: Lb,by,bz
      real(WP), dimension(:,:), allocatable :: U
      !> Add these for the SGS model
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU ! Adjust dimensions to match your solver's signature
      real(WP), dimension(:,:,:), allocatable :: resU
      real(WP), dimension(:,:,:,:), allocatable :: SR
      real(WP), dimension(:,:,:), allocatable :: U_mean, V_mean, W_mean
   contains
      procedure :: analyze
      procedure, private :: read_ensight_scalar
      procedure, private :: read_ensight_vector
      procedure, private :: read_ensight_part
      procedure, private :: extract_core
      procedure, private :: analyze_core
      procedure, private :: extract_dissipation
      procedure, private :: compute_mean_fields
      ! procedure, private :: extract_EPL
      ! procedure, private :: analyze_EPL
      ! procedure, private :: extract_InletVel
      ! procedure, private :: analyze_InletVel
      ! procedure, private :: extract_droplets
   end type postproc

   real(WP), parameter, public :: dl=0.003_WP   ! Liquid outer pipe diameter 
   real(WP), parameter, public :: dg=0.0206_WP   ! Gas pipe diameter ~(inner+outer)/2
   real(WP), parameter, public :: rl=0.0010_WP   ! Liquid pipe inner radius
   real(WP), parameter, public :: rlo=0.0015_WP   ! Liquid pipe outer radius
      
contains
   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary

   !> Function that localizes liquid stream at -x
   function liq_inlet(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: rad
      isIn=.false.
      rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (rad.lt.0.5_WP*dl.and.i.eq.pg%imin) isIn=.true.
   end function liq_inlet

   !> Function that localizes gas stream at -x
   function gas_inlet(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: rad
      isIn=.false.
      rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (rad.ge.0.5_WP*dl.and.rad.lt.0.5_WP*dg.and.i.eq.pg%imin) isIn=.true.
   end function gas_inlet

   !> Function that localizes the top (y+) of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator

   !> Function that localizes the bottom (y-) of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator

   !> Function that localizes the top (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator

   !> Function that localizes the bottom (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator
   
   !> Read a scalar ensight file to an WP array - handle ghost cells as well
   subroutine read_ensight_scalar(this,filename,SC)
      use mpi_f08
      use parallel, only: group,info_mpiio,MPI_REAL_SP
      use string,   only: str_medium
      use messager, only: die
      class(postproc), intent(inout) :: this
      character(len=str_medium), intent(in) :: filename
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: SC
      integer(kind=MPI_OFFSET_KIND) :: disp
      real(SP), dimension(:,:,:), allocatable :: spbuff
      type(MPI_Status):: status
      type(MPI_File) :: ifile
      integer :: ierr,i
      ! Zero out SC
      SC=0.0_WP
      ! Parallel read the file
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
      disp=244
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      call MPI_FILE_READ_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      SC(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),WP)
      call MPI_FILE_CLOSE(ifile,ierr)
      deallocate(spbuff)
      ! Update ghost cells
      call this%cfg%sync(SC)
      if (this%cfg%iproc.eq.1) then
         do i=this%cfg%imino_,this%cfg%imin_-1
            SC(i,:,:)=SC(this%cfg%imin_,:,:)
         end do
      end if
   end subroutine read_ensight_scalar
   
   subroutine read_ensight_vector(this,filename,U,V,W)
      use mpi_f08
      use parallel, only: group,info_mpiio,MPI_REAL_SP
      use string,   only: str_medium
      use messager, only: die
      class(postproc), intent(inout) :: this
      character(len=str_medium), intent(in) :: filename
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: U
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: V
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: W
      integer(kind=MPI_OFFSET_KIND) :: disp
      real(SP), dimension(:,:,:), allocatable :: spbuff
      type(MPI_Status):: status
      type(MPI_File) :: ifile
      integer :: ierr,i
      ! Zero out SC
      U=0.0_WP;V=0.0_WP;W=0.0_WP
      ! Parallel read the file
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
      disp=244
      allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      U(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),WP)
      disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      V(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),WP)
      disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      W(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),WP)
      call MPI_FILE_CLOSE(ifile,ierr)
      deallocate(spbuff)
      ! Update ghost cells
      call this%cfg%sync(U)
      call this%cfg%sync(V)
      call this%cfg%sync(W)
   end subroutine read_ensight_vector

   subroutine read_ensight_part(this,nfile)
      use mpi_f08
      use parallel, only: group,info_mpiio,MPI_REAL_SP
      use string,   only: str_medium
      use messager, only: die
      class(postproc), intent(inout) :: this
      integer, intent(in) :: nfile
      integer(kind=MPI_OFFSET_KIND) :: disp
      character(len=str_medium) :: filename
      type(MPI_Status):: status
      type(MPI_File) :: ifile
      integer :: i, npart,ierr,iunit
      real(SP), dimension(:,:), allocatable :: sppos_buff,spvel_buff
      real(SP), dimension(:), allocatable ::   spradius_buff,spid_buff,tmp_test
      if (this%cfg%amRoot) then
         filename='ensight/atom/part/particle.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call MPI_FILE_OPEN(MPI_COMM_SELF,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
         if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
         disp=240_MPI_OFFSET_KIND
         ! Read in Npart
         call MPI_FILE_READ_AT(ifile, 240_MPI_OFFSET_KIND, npart, 1, MPI_REAL_SP,status,ierr)
         if (npart.gt.0) then
            ! If we have some particles, first read positions
            allocate(sppos_buff(3,npart))
            disp = 240_MPI_OFFSET_KIND+ int(npart+1, MPI_OFFSET_KIND) * int(SP,MPI_OFFSET_KIND)
            call MPI_FILE_READ_AT(ifile, disp, sppos_buff, 3*npart, MPI_REAL_SP,status,ierr)
            call MPI_FILE_CLOSE(ifile,ierr)
            ! Then read velocities
            allocate(spvel_buff(3,npart))
            filename='ensight/atom/part/velocity.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
            call MPI_FILE_OPEN(MPI_COMM_SELF,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
            if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
            call MPI_FILE_READ_AT(ifile, 80_MPI_OFFSET_KIND, spvel_buff, 3*npart, MPI_REAL_SP,status,ierr)
            call MPI_FILE_CLOSE(ifile,ierr)

            ! Then read radius
            allocate(spradius_buff(npart))
            filename='ensight/atom/part/radius.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
            call MPI_FILE_OPEN(MPI_COMM_SELF,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
            if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
            call MPI_FILE_READ_AT(ifile, 80_MPI_OFFSET_KIND, spradius_buff, npart, MPI_REAL_SP,status,ierr)
            call MPI_FILE_CLOSE(ifile,ierr)

            ! Then read ID
            allocate(spid_buff(npart))
            filename='ensight/atom/part/id.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
            call MPI_FILE_OPEN(MPI_COMM_SELF,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
            if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
            call MPI_FILE_READ_AT(ifile, 80_MPI_OFFSET_KIND, spid_buff, npart, MPI_REAL_SP,status,ierr)
            call MPI_FILE_CLOSE(ifile,ierr)

            ! Add all the particles into my lpt object class
            call this%lp%resize(npart)
            do i = 1,npart
               this%lp%np_=this%lp%np_+1
               this%lp%p(this%lp%np_)%id=int(spid_buff(i),WP)
               this%lp%p(this%lp%np_)%d =real(spradius_buff(i)*2.0,WP)
               this%lp%p(this%lp%np_)%vel=real(spvel_buff(:,i),WP)
               this%lp%p(this%lp%np_)%pos=real(sppos_buff(:,i),WP)
               this%lp%p(this%lp%np_)%flag=0                                                                                        
               this%lp%p(this%lp%np_)%dt  =0.0_WP                                                                                  
               this%lp%p(this%lp%np_)%Acol=0.0_WP                                                                                  
               this%lp%p(this%lp%np_)%Tcol=0.0_WP
               this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])     
            end do
         else
            call MPI_FILE_CLOSE(ifile,ierr)
         end if
      end if
      call this%lp%sync()
   end subroutine read_ensight_part

   !> Extract a pmesh skeleton of the liquid core from CCL data
   subroutine extract_core(this,VFtmp,nfile)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: twoPi
      use vfs_class, only: VFlo
      implicit none
      class(postproc), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: VFtmp
      integer, intent(in) :: nfile
      real(WP), dimension(:), allocatable :: dvol,xmax!,mybv,myby,mybz
      real(WP) :: x,y,z
      real(WP) :: mybv,myby,mybz
      integer:: i,j,k,n,m,ierr
      integer, dimension(1) :: idx
      ! Start by building a CCL
      call this%ccl%build(make_label,same_label)
      ! Allocate droplet stats arrays
      allocate(dvol(1:this%ccl%nstruct)); dvol=0.0_WP
      allocate(xmax(1:this%ccl%nstruct)); xmax=-HUGE(x)
      ! First pass to accumulate volume, position, and velocity
      do n=1,this%ccl%nstruct
         ! Loop over cells in structure
         do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
            dvol(n)=dvol(n)+this%cfg%vol(i,j,k)*VFtmp(i,j,k)
            xmax(n)=max(xmax(n),this%cfg%x(i))
         end do 
      end do 
      ! Get the structure with the largest liquid volume
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      idx=maxloc(dvol); n=idx(1); this%Lb(nfile)=xmax(n); mybv=0.0_WP;myby=0.0_WP;mybz=0.0_WP
      ! Calculate x-dependent y and z bary centers
      do m=1,this%ccl%struct(n)%n_
         ! Get cell indices
         i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
         ! Integrate barycenter only for liquid core outside exit
         if (this%cfg%xm(i).ge.0.0_WP) then
            mybv=mybv+VFtmp(i,j,k)*this%cfg%vol(i,j,k)
            myby=myby+VFtmp(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%ym(j)
            mybz=mybz+VFtmp(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%zm(k)
         end if
      end do 
      call MPI_ALLREDUCE(MPI_IN_PLACE,mybv,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,myby,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mybz,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

      if (mybv.gt.0.0_WP) then
         this%by(nfile)=myby/mybv
         this%bz(nfile)=mybz/mybv         
      end if
      contains
      
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (VFtmp(i,j,k).gt.VFlo) then
            make_label=.true.
         else
            make_label=.false.
         end if
      end function make_label
      
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
   end subroutine extract_core
   
   subroutine analyze_core(this,fstart,fend)
      use string,   only: str_medium
      implicit none
      class(postproc), intent(inout) :: this
      integer, intent(in) :: fstart,fend
      character(len=str_medium) :: filename
      integer :: i,j
      if (this%cfg%amRoot) then
         filename="Lb.csv"
         ! Open file dynamically with append mode
         open(unit=10, file=filename, status="replace", action="write")
         do i=fstart,fend
            write(10, '(F24.16, ",", F24.16)') i*1.0_WP,this%Lb(i)
         end do
         close(unit=10)

         filename="yzbary.csv"
         ! Open file dynamically with append mode
         open(unit=10, file=filename, status="replace", action="write")
         do i=fstart,fend
            write(10, '(F24.16, ",", F24.16,",", F24.16)') i*1.0_WP,this%by(i),this%bz(i)
         end do
         close(unit=10)
      end if
   end subroutine analyze_core

   ! !> Extract a liquid droplets from CCL data
   ! subroutine extract_EPL(this,VF,nfile)
   !    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MAX
   !    use parallel,  only: MPI_REAL_WP
   !    implicit none
   !    class(postproc), intent(inout) :: this
   !    integer, intent(in) :: nfile
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: VF

   !    integer :: i,j,k
   !    ! Sum all the VF in z direction
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             this%VF(i,j,nfile)=this%VF(i,j,nfile)+VF(i,j,k)*this%cfg%dz(k)
   !          end do
   !       end do
   !    end do

   ! end subroutine extract_EPL

   ! !> Extract a liquid droplets from CCL data
   ! subroutine analyze_EPL(this,dir,xconst,yconst)
   !    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MAX
   !    use parallel,  only: MPI_REAL_WP
   !    use string,   only: str_medium
   !    implicit none
   !    class(postproc), intent(inout) :: this
   !    character(len=str_medium) :: filename
   !    integer, intent(in) :: dir
   !    real(WP), intent(in) :: xconst,yconst
   !    real(WP), dimension(:), allocatable :: EPL_avg,EPL_std
   !    integer :: i,j,ierr
   !       select case (dir)
   !          ! constant in x
   !          case(1)
   !             ! allocate the list accounting for the whole y range
   !             allocate(EPL_avg(this%cfg%jmin:this%cfg%jmax));EPL_avg=0.0_WP
   !             allocate(EPL_std(this%cfg%jmin:this%cfg%jmax));EPL_std=0.0_WP
   !             do j=this%cfg%jmin_,this%cfg%jmax_
   !                do i=this%cfg%imin_,this%cfg%imax_
   !                   if (this%cfg%xm(i-1).lt.xconst.and.this%cfg%xm(i).ge.xconst) then
   !                      EPL_avg(j)=sum(this%VF(i,j,:))/size(this%VF(i,j,:))
   !                      EPL_std(j)=sqrt(sum((this%VF(i,j,:)-EPL_avg(j))**2)/size(this%VF(i,j,:)))
   !                   end if
   !                end do
   !             end do
   !             call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_avg,size(EPL_avg),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
   !             call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_std,size(EPL_std),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
   !             if (this%cfg%amRoot) then
   !                filename="x.csv"
   !                ! filename="x="
   !                ! write(filename, '(F0.2, ".csv")') xconst
   !                ! Open file dynamically with append mode
   !                open(unit=10, file=filename, status="replace", action="write")
   !                do j=this%cfg%jmin,this%cfg%jmax
   !                   write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%ym(j),EPL_avg(j),EPL_std(j)
   !                end do
   !                close(unit=10)
   !             end if
   !          ! constant in y
   !          case(2)
   !             ! allocate the list accounting for the whole x range
   !             allocate(EPL_avg(this%cfg%imin:this%cfg%imax));EPL_avg=0.0_WP
   !             allocate(EPL_std(this%cfg%imin:this%cfg%imax));EPL_std=0.0_WP
   !             do i=this%cfg%imin_,this%cfg%imax_
   !                do j=this%cfg%jmin_,this%cfg%jmax_
   !                   if (this%cfg%ym(j-1).lt.yconst.and.this%cfg%ym(j).ge.yconst) then
   !                      EPL_avg(i)=sum(this%VF(i,j,:))/size(this%VF(i,j,:))
   !                      EPL_std(i)=sqrt(sum((this%VF(i,j,:)-EPL_avg(i))**2)/size(this%VF(i,j,:)))
   !                      ! print *, "i1",this%cfg%xm(i),EPL_avg(i),EPL_std(i),sum(this%VF(i,j,:))/size(this%VF(i,j,:)),sqrt(sum((this%VF(i,j,:)-EPL_avg(i))**2)/size(this%VF(i,j,:)))
   !                   end if
   !                end do
   !             end do
   !             call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_avg,size(EPL_avg),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
   !             call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_std,size(EPL_std),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
   !             if (this%cfg%amRoot) then
   !                filename="y.csv"
   !                ! filename="y="
   !                ! write(filename, '(F0.2, ".csv")') xconst
   !                ! Open file dynamically with append mode
   !                open(unit=10, file=filename, status="replace", action="write")
   !                do i=this%cfg%imin,this%cfg%imax
   !                   write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%xm(i),EPL_avg(i),EPL_std(i)
   !                end do
   !                close(unit=10)
   !             end if
   !       end select
   ! end subroutine analyze_EPL

!   !> Extract a pmesh skeleton of the liquid core from CCL data
!    subroutine extract_InletVel(this,U,nfile)
!       use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
!       use parallel,  only: MPI_REAL_WP
!       implicit none
!       class(postproc), intent(inout) :: this
!       real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: U
!       real(WP), dimension(:), allocatable :: Utmp
!       real(WP) :: velmesurement
!       integer, intent(in) :: nfile
!       integer:: i,j,k,ierr
!       velmesurement=0.0003_WP
!       allocate(Utmp(this%cfg%jmin:this%cfg%jmax));Utmp=0.0_WP
!       ! Loop through x domain and get x location
!       do k=this%cfg%kmin_,this%cfg%kmax_
!          do i=this%cfg%imin_,this%cfg%imax_
!             if (this%cfg%zm(k).ge.0.0_WP .and. this%cfg%zm(k-1).lt.0.0_WP ) then
!                if (this%cfg%xm(i).ge.velmesurement .and. this%cfg%xm(i-1).lt.velmesurement ) then
!                   do j = this%cfg%jmin_,this%cfg%jmax_
!                      Utmp(j) = U(i,j,k)
!                   end do
!                end if
!             end if
!          end do
!       end do
!       call MPI_ALLREDUCE(MPI_IN_PLACE,Utmp,size(Utmp),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
!       this%U(:,nfile) = Utmp
!    end subroutine extract_InletVel

!    !> Extract a pmesh skeleton of the liquid core from CCL data
!    subroutine analyze_InletVel(this)
!       use mathtools, only: Pi
!       use string,   only: str_medium
!       implicit none
!       class(postproc), intent(inout) :: this
!       character(len=str_medium) :: filename
!       real(WP), dimension(:), allocatable :: Uinlet_avg,Uinlet_std
!       integer :: j
!       real(WP), parameter :: SLPM2SI=1.66667E-5_WP
!       real(WP) :: dg,Qaxial,Uaxial,Aaxial,dl
!       dg=0.01_WP
!       dl=0.003_WP
!       ! call this%input%read('Total flow rate (SLPM)',Qaxial)
!       Qaxial=150.0_WP*SLPM2SI
!       Aaxial=0.25_WP*Pi*(dg**2-dl**2)
!       Uaxial=Qaxial/Aaxial 
!       allocate(Uinlet_avg(this%cfg%jmin:this%cfg%jmax));Uinlet_avg=0.0_WP
!       allocate(Uinlet_std(this%cfg%jmin:this%cfg%jmax));Uinlet_std=0.0_WP
!       do j=this%cfg%jmin,this%cfg%jmax
!          Uinlet_avg(j)=sum(this%U(j,:))/size(this%U(j,:))
!          Uinlet_std(j)=sqrt(sum((this%U(j,:)-Uinlet_avg(j))**2)/size(this%U(j,:)))
!       end do
!       if (this%cfg%amRoot) then
!          filename="Uinlet.csv"
!          ! Open file dynamically with append mode
!          open(unit=10, file=filename, status="replace", action="write")
!          do j=this%cfg%jmin,this%cfg%jmax
!             write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%ym(j)/dg,Uinlet_avg(j)/Uaxial,Uinlet_std(j)/Uaxial
!          end do
!          close(unit=10)
!       end if

!    end subroutine analyze_InletVel

   ! !> Extract a pmesh skeleton of the liquid core from CCL data
   ! subroutine extract_droplets(this)
   !    use mathtools, only: Pi
   !    use string,   only: str_medium
   !    implicit none
   !    class(postproc), intent(inout) :: this
   !    character(len=str_medium) :: filename
   !    real(WP), dimension(:), allocatable :: Uinlet_avg,Uinlet_std
   !    integer :: j
   !    real(WP), parameter :: SLPM2SI=1.66667E-5_WP
   !    real(WP) :: dg,Qaxial,Uaxial,Aaxial,dl
   !    dg=0.01_WP
   !    dl=0.003_WP
   !    Qaxial=2.0_WP*85.7_WP*SLPM2SI
   !    Aaxial=0.25_WP*Pi*(dg**2-dl**2)
   !    Uaxial=Qaxial/Aaxial 
   !    allocate(Uinlet_avg(this%cfg%jmin:this%cfg%jmax));Uinlet_avg=0.0_WP
   !    allocate(Uinlet_std(this%cfg%jmin:this%cfg%jmax));Uinlet_std=0.0_WP
   !    do j=this%cfg%jmin,this%cfg%jmax
   !       Uinlet_avg(j)=sum(this%U(j,:))/size(this%U(j,:))
   !       Uinlet_std(j)=sqrt(sum((this%U(j,:)-Uinlet_avg(j))**2)/size(this%U(j,:)))
   !    end do
   !    if (this%cfg%amRoot) then
   !       filename="Uinlet.csv"
   !       ! Open file dynamically with append mode
   !       open(unit=10, file=filename, status="replace", action="write")
   !       do j=this%cfg%jmin,this%cfg%jmax
   !          write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%ym(j)/dg,Uinlet_avg(j)/Uaxial,Uinlet_std(j)/Uaxial
   !       end do
   !       close(unit=10)
   !    end if

   ! end subroutine extract_droplets

   !> Calculate local dissipation, Vreman SGS viscosity, and area-weighted scales
   !> Calculate local dissipation and area-weighted scales using built-in SGS model
   !> Calculate local dissipation using native solver routines
   subroutine extract_dissipation(this, U, V, W, VF, nfile, fstart)
      use mpi_f08,   only: MPI_ALLREDUCE, MPI_SUM, MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use sgsmodel_class, only: vreman
      implicit none
      class(postproc), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U, V, W, VF
      integer, intent(in) :: nfile,fstart

      real(WP) :: nu_phys, rho_c, sigma, nu_sgs, eps_local, S2
      real(WP) :: gradVF_x, gradVF_y, gradVF_z, mag_gradVF, vol
      real(WP) :: dx, dy, dz, dt
      real(WP), dimension(:), allocatable :: eps_sum_x, area_sum_x, d_H_x
      integer :: i, j, k, ierr

      ! Extract physical properties
      nu_phys = this%fs%visc_g / this%fs%rho_g
      rho_c   = this%fs%rho_g
      sigma   = this%fs%sigma
      dt = 1.1e-6_WP

      !> 1. Pass Ensight velocity to the flow solver
      this%fs%U = U
      this%fs%V = V
      this%fs%W = W
      
      this%resU = this%fs%rho_g 
      call this%fs%get_gradu(this%gradU)
      call this%sgs%get_visc(type=vreman, dt=dt, rho=this%resU, gradu=this%gradU)

      this%fs%U = U-this%U_mean
      this%fs%V = V-this%V_mean
      this%fs%W = W-this%W_mean
      
      call this%fs%get_strainrate(this%SR)

      allocate(eps_sum_x(this%cfg%imin:this%cfg%imax));  eps_sum_x = 0.0_WP
      allocate(area_sum_x(this%cfg%imin:this%cfg%imax)); area_sum_x = 0.0_WP
      allocate(d_H_x(this%cfg%imin:this%cfg%imax));      d_H_x = 0.0_WP

      !> 3. Loop to calculate Dissipation & Area Weights
      do k = this%cfg%kmin_, this%cfg%kmax_
         dz = this%cfg%dz(k)
         do j = this%cfg%jmin_, this%cfg%jmax_
            dy = this%cfg%dy(j)
            do i = this%cfg%imin_, this%cfg%imax_
               dx = this%cfg%dx(i)
               vol = dx * dy * dz

               ! Kinematic SGS Viscosity
               nu_sgs = this%sgs%visc(i,j,k) / rho_c

               ! S2 Contraction exactly from your SR array mapping:
               ! SR(1,2,3) are diagonal (S11, S22, S33)
               ! SR(4,5,6) are off-diagonal (S12, S23, S13)
               S2 = this%SR(1,i,j,k)**2 + this%SR(2,i,j,k)**2 + this%SR(3,i,j,k)**2 + &
                    2.0_WP * (this%SR(4,i,j,k)**2 + this%SR(5,i,j,k)**2 + this%SR(6,i,j,k)**2)

               ! Total Dissipation
               eps_local = 2.0_WP * (nu_phys + nu_sgs) * S2

               ! VOF Gradient for Interfacial Area Weighting
               gradVF_x = (VF(i+1,j,k) - VF(i-1,j,k)) / (2.0_WP * dx)
               gradVF_y = (VF(i,j+1,k) - VF(i,j-1,k)) / (2.0_WP * dy)
               gradVF_z = (VF(i,j,k+1) - VF(i,j,k-1)) / (2.0_WP * dz)
               mag_gradVF = sqrt(gradVF_x**2 + gradVF_y**2 + gradVF_z**2)

               eps_sum_x(i)  = eps_sum_x(i)  + (eps_local * mag_gradVF * vol)
               area_sum_x(i) = area_sum_x(i) + (mag_gradVF * vol)

            end do
         end do
      end do

      !> 4. Parallel Reduction
      call MPI_ALLREDUCE(MPI_IN_PLACE, eps_sum_x,  size(eps_sum_x),  MPI_REAL_WP, MPI_SUM, this%cfg%comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, area_sum_x, size(area_sum_x), MPI_REAL_WP, MPI_SUM, this%cfg%comm, ierr)

      !> 5. Output (Streaming Tidy Data)
      if (this%cfg%amRoot) then
         open(unit=11, file="Hinze_Scale_x.csv", position="append", action="write")
         ! Write header only on the first file
         if (nfile == fstart) write(11, '("nfile, x, eps_avg, d_H")') 
         
         do i = this%cfg%imin, this%cfg%imax
            if (area_sum_x(i) .gt. 0.0_WP) then
               eps_sum_x(i) = eps_sum_x(i) / area_sum_x(i)
               d_H_x(i) = 0.725_WP * (sigma / rho_c)**0.6_WP * eps_sum_x(i)**(-0.4_WP)
               
               ! Add nfile to the output string
               write(11, '(I6, ",", F24.16, ",", F24.16, ",", F24.16)') &
                     nfile, this%cfg%xm(i), eps_sum_x(i), d_H_x(i)
            end if
         end do
         close(unit=11)
      end if
      
      deallocate(eps_sum_x, area_sum_x, d_H_x)
   end subroutine extract_dissipation

   subroutine compute_mean_fields(this, fstart, fend)
      use string,   only: str_medium
      use messager, only: log
      implicit none
      class(postproc), intent(inout) :: this
      integer, intent(in) :: fstart, fend
      
      real(WP), dimension(:,:,:), allocatable :: U_tmp, V_tmp, W_tmp
      character(len=str_medium) :: filename
      integer :: nfile
      real(WP) :: nfiles_total
      
      allocate(U_tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(V_tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(W_tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      call log('--- Starting Pre-Processing: Computing Mean Velocity Fields ---')
      
      do nfile = fstart, fend
         filename='ensight/atom/velocity/velocity.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call this%read_ensight_vector(filename, U_tmp, V_tmp, W_tmp)
         
         this%U_mean = this%U_mean + U_tmp
         this%V_mean = this%V_mean + V_tmp
         this%W_mean = this%W_mean + W_tmp
      end do
      
      nfiles_total = real(fend - fstart + 1, WP)
      this%U_mean = this%U_mean / nfiles_total
      this%V_mean = this%V_mean / nfiles_total
      this%W_mean = this%W_mean / nfiles_total
      
      ! Sync the ghost cells for the mean fields
      call this%cfg%sync(this%U_mean)
      call this%cfg%sync(this%V_mean)
      call this%cfg%sync(this%W_mean)
      
      deallocate(U_tmp, V_tmp, W_tmp)
      call log('--- Mean Velocity Fields Computed Successfully ---')
   end subroutine compute_mean_fields

   !> Analysis of atom simulation
   subroutine analyze(this)
      use parallel, only: amRoot
      use string,   only: str_medium
      use messager, only: log
      implicit none
      class(postproc), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: VFtmp,U,V,W
      character(len=str_medium) :: filename
      integer :: nfile,fstart,fend
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_postproc')
      ! Create the config
      create_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='postproc')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use tpns_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set the flow properties
         call this%input%read('Liquid dynamic viscosity',this%fs%visc_l)
         call this%input%read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call this%input%read('Liquid density',this%fs%rho_l)
         call this%input%read('Gas density'   ,this%fs%rho_g)
         call this%input%read('Surface tension coefficient',this%fs%sigma)
         ! Define gas and liquid inlet boundary conditions
         call this%fs%add_bcond(name='gas_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=gas_inlet)
         call this%fs%add_bcond(name='liq_inlet',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=liq_inlet)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
      end block create_flow_solver

      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs

      call this%input%read('File start',fstart); call this%input%read('File end',fend);
      ! Initialize CCL
      call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      ! Allocate work arrays
      allocate_data: block
         ! allocate(this%VF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,fstart:fend)); this%VF=0.0_WP
         ! allocate(this%Lb(fstart:fend));this%Lb=0.0_WP
         ! allocate(this%by(fstart:fend)); this%by=0.0_WP
         ! allocate(this%bz(fstart:fend)); this%bz=0.0_WP
         ! allocate(this%U (this%cfg%jmin:this%cfg%jmax,fstart:fend)); this%U =0.0_WP
         allocate(VFtmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); VFtmp=0.0_WP
         allocate(U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); U=0.0_WP
         allocate(V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); V=0.0_WP
         allocate(W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); W=0.0_WP
         allocate(this%U_mean(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U_mean=0.0_WP
         allocate(this%V_mean(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V_mean=0.0_WP
         allocate(this%W_mean(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W_mean=0.0_WP
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_data
      


      ! initialze_lpt: block
      !    this%lp=lpt(cfg=this%cfg,name='spray_analyze')
      !    call this%lp%resize(0)
      ! end block initialze_lpt

      ! create_pmesh: block
      !    integer :: i
      !    this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
      !    this%pmesh%varname(1)='radius'
      !    this%pmesh%varname(2)='id'
      !    this%pmesh%vecname(1)='velocity'
      !    call this%lp%update_partmesh(this%pmesh)
      !    do i=1,this%lp%np_
      !       this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
      !       this%pmesh%var(2,i)=this%lp%p(i)%id
      !       this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
      !    end do
      ! end block create_pmesh

      ! testlpt_ensight: block
      !    this%ens_out=ensight(cfg=this%cfg,name='spray_analyze')
      !    call this%ens_out%add_particle('part',this%pmesh)
      ! end block testlpt_ensight
      call this%compute_mean_fields(fstart,fend)

      ! Run on all files available
      do nfile=fstart,fend
         filename='ensight/atom/VOF/VOF.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call log('Postprocessing file '//trim(filename)//'...')
         call this%read_ensight_scalar(filename,VFtmp)
         ! call log('|----> VOF read successfully')
         filename='ensight/atom/velocity/velocity.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call log('Postprocessing file '//trim(filename)//'...')
         call this%read_ensight_vector(filename,U,V,W)
         call this%extract_dissipation(U,V,W,VFtmp,nfile,fstart)
         ! call this%extract_EPL(VF=VFtmp,nfile=nfile)
         ! call log('|----> EPL calculation done')
         ! call this%extract_core(VFtmp=VFtmp,nfile=nfile)
         ! call log('|----> liquid core extracted')
         ! call this%extract_InletVel(U=U,nfile=nfile)
         ! call log('|----> Inlet Vel extracted')
         ! call this%read_ensight_part(nfile=nfile)
      end do
      ! update_pmesh: block
      !    integer :: i
      !    call this%lp%update_partmesh(this%pmesh)
      !    do i=1,this%lp%np_
      !       this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
      !       this%pmesh%var(2,i)=this%lp%p(i)%id
      !       this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
      !    end do
      ! end block update_pmesh 
      ! call this%ens_out%write_data(this%time%t)
      
      ! call this%analyze_EPL(dir=2,xconst=0.0_WP,yconst=0.0_WP)
      ! call this%analyze_core(fstart,fend)
      ! call this%analyze_InletVel()
   end subroutine analyze
   

end module postproc_class