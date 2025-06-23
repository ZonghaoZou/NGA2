!> Definition for a postproc class
module postproc_class
   use precision,         only: WP,SP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use timetracker_class, only: timetracker
   use cclabel_class,     only: cclabel
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
      type(cclabel)     :: ccl
      !> CCL analysis
      ! type(ccl)  :: cc
      !> Data arrays
      real(WP), dimension(:,:,:), allocatable :: VF
      real(WP), dimension(:), allocatable :: Lb
      ! real(WP), dimension(:,:,:), allocatable :: U
      real(WP), dimension(:,:), allocatable :: by,bz
   contains
      procedure :: analyze
      procedure, private :: read_ensight_scalar
      procedure, private :: extract_core
      procedure, private :: analyze_core
      procedure, private :: extract_EPL
      procedure, private :: analyze_EPL
   end type postproc
      
contains
   
   
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
      real(WP), dimension(:), allocatable :: dvol,xmax,mybv,myby,mybz
      real(WP) :: x,y,z
      integer:: i,j,k,n,m,ierr
      integer, dimension(1) :: idx
      ! Start by building a CCL
      call this%ccl%build(make_label,same_label)
      ! Allocate droplet stats arrays
      allocate(dvol(1:this%ccl%nstruct)); dvol=0.0_WP
      allocate(xmax(1:this%ccl%nstruct)); xmax=-HUGE(x)
      allocate(mybv(this%cfg%imin:this%cfg%imax)); mybv=0.0_WP
      allocate(myby(this%cfg%imin:this%cfg%imax)); myby=0.0_WP
      allocate(mybz(this%cfg%imin:this%cfg%imax)); mybz=0.0_WP
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
      idx=maxloc(dvol); n=idx(1); this%Lb(nfile)=xmax(n) 
      ! Calculate x-dependent y and z bary centers
      do m=1,this%ccl%struct(n)%n_
         ! Get cell indices
         i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
         ! Integrate barycenter
         mybv(i)=mybv(i)+VFtmp(i,j,k)*this%cfg%vol(i,j,k)
         myby(i)=myby(i)+VFtmp(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%ym(j)
         mybz(i)=mybz(i)+VFtmp(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%zm(k)
      end do 

      call MPI_ALLREDUCE(MPI_IN_PLACE,mybv,size(mybv),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,myby,size(myby),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mybz,size(mybz),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)

      do i = this%cfg%imin,this%cfg%imax
         if (mybv(i).gt.0.0_WP) then
            this%by(i,nfile)=myby(i)/mybv(i)
            this%bz(i,nfile)=mybz(i)/mybv(i)
         end if
      end do
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
      end if
   end subroutine analyze_core


   !> Extract a liquid droplets from CCL data
   subroutine extract_EPL(this,VF,nfile)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(postproc), intent(inout) :: this
      integer, intent(in) :: nfile
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: VF

      integer :: i,j,k
      ! Sum all the VF in z direction
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%VF(i,j,nfile)=this%VF(i,j,nfile)+VF(i,j,k)*this%cfg%dz(k)
            end do
         end do
      end do

   end subroutine extract_EPL
   

   !> Extract a liquid droplets from CCL data
   subroutine analyze_EPL(this,dir,xconst,yconst)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      use string,   only: str_medium
      implicit none
      class(postproc), intent(inout) :: this
      character(len=str_medium) :: filename
      integer, intent(in) :: dir
      real(WP), intent(in) :: xconst,yconst
      real(WP), dimension(:), allocatable :: EPL_avg,EPL_std
      integer :: i,j,ierr
         select case (dir)
            ! constant in x
            case(1)
               ! allocate the list accounting for the whole y range
               allocate(EPL_avg(this%cfg%jmin:this%cfg%jmax));EPL_avg=0.0_WP
               allocate(EPL_std(this%cfg%jmin:this%cfg%jmax));EPL_std=0.0_WP
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%cfg%imin_,this%cfg%imax_
                     if (this%cfg%xm(i-1).lt.xconst.and.this%cfg%xm(i).ge.xconst) then
                        EPL_avg(j)=sum(this%VF(i,j,:))/size(this%VF(i,j,:))
                        EPL_std(j)=sqrt(sum((this%VF(i,j,:)-EPL_avg(j))**2)/size(this%VF(i,j,:)))
                     end if
                  end do
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_avg,size(EPL_avg),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_std,size(EPL_std),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
               if (this%cfg%amRoot) then
                  filename="x.csv"
                  ! filename="x="
                  ! write(filename, '(F0.2, ".csv")') xconst
                  ! Open file dynamically with append mode
                  open(unit=10, file=filename, status="replace", action="write")
                  do j=this%cfg%jmin,this%cfg%jmax
                     write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%ym(j),EPL_avg(j),EPL_std(j)
                  end do
                  close(unit=10)
               end if
            ! constant in y
            case(2)
               ! allocate the list accounting for the whole x range
               allocate(EPL_avg(this%cfg%imin:this%cfg%imax));EPL_avg=0.0_WP
               allocate(EPL_std(this%cfg%imin:this%cfg%imax));EPL_std=0.0_WP
               do i=this%cfg%imin_,this%cfg%imax_
                  do j=this%cfg%jmin_,this%cfg%jmax_
                     if (this%cfg%ym(j-1).lt.yconst.and.this%cfg%ym(j).ge.yconst) then
                        EPL_avg(i)=sum(this%VF(i,j,:))/size(this%VF(i,j,:))
                        EPL_std(i)=sqrt(sum((this%VF(i,j,:)-EPL_avg(i))**2)/size(this%VF(i,j,:)))
                        ! print *, "i1",this%cfg%xm(i),EPL_avg(i),EPL_std(i),sum(this%VF(i,j,:))/size(this%VF(i,j,:)),sqrt(sum((this%VF(i,j,:)-EPL_avg(i))**2)/size(this%VF(i,j,:)))
                     end if
                  end do
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_avg,size(EPL_avg),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,EPL_std,size(EPL_std),MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
               if (this%cfg%amRoot) then
                  filename="y.csv"
                  ! filename="y="
                  ! write(filename, '(F0.2, ".csv")') xconst
                  ! Open file dynamically with append mode
                  open(unit=10, file=filename, status="replace", action="write")
                  do i=this%cfg%imin,this%cfg%imax
                     write(10, '(F24.16, ",", F24.16, ",", F24.16)') this%cfg%xm(i),EPL_avg(i),EPL_std(i)
                  end do
                  close(unit=10)
               end if
         end select
   end subroutine analyze_EPL

   !> Analysis of atom simulation
   subroutine analyze(this)
      use parallel, only: amRoot
      use string,   only: str_medium
      use messager, only: log
      implicit none
      class(postproc), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: VFtmp
      character(len=str_medium) :: filename
      integer :: nfile,fstart,fend
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_postproc')
      ! input_atom
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
      
      call this%input%read('File start',fstart); call this%input%read('File end',fend);
      ! Initialize CCL
      call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      ! Allocate work arrays
      allocate_data: block
         allocate(this%VF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,fstart:fend)); this%VF=0.0_WP
         allocate(this%Lb(fstart:fend));this%Lb=0.0_WP
         allocate(this%by(this%cfg%imin:this%cfg%imax,fstart:fend)); this%by=0.0_WP
         allocate(this%bz(this%cfg%imin:this%cfg%imax,fstart:fend)); this%bz=0.0_WP
         ! allocate(this%U (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
         ! allocate(this%Lb(fstart:fend)); this%Lb=0.0_WP
         ! allocate(this%bv(this%cfg%imino:this%cfg%imaxo)); this%bv=0.0_WP
         allocate(VFtmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); VFtmp=0.0_WP
      end block allocate_data
      

      
      ! Run on all files available
      do nfile=fstart,fend
         filename='ensight/atom/VOF/VOF.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call log('Postprocessing file '//trim(filename)//'...')
         call this%read_ensight_scalar(filename,VFtmp)
         call log('|----> File read successfully')
         call this%extract_EPL(VF=VFtmp,nfile=nfile)
         call log('|----> EPL calculation done')
         call this%extract_core(VFtmp=VFtmp,nfile=nfile)
         call log('|----> liquid core extracted')
      end do
      
      call this%analyze_EPL(dir=2,xconst=0.0_WP,yconst=0.0_WP)
      call this%analyze_core(fstart,fend)
   end subroutine analyze
   

end module postproc_class