!> Definition for a breakup class that handles VOF-to-Lagrangian transfer
module breakup_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use lpt_class,         only: lpt
   use cclabel_class,     only: cclabel
   use timetracker_class, only: timetracker
   use string,            only: str_medium

   implicit none
   private

   type :: spline_info
      integer :: n_knots_x, n_knots_y, n_knots_z, flag
      real(WP) :: length
      real(WP), dimension(:), allocatable :: t_knots_x, c_coeffs_x
      real(WP), dimension(:), allocatable :: t_knots_y, c_coeffs_y
      real(WP), dimension(:), allocatable :: t_knots_z, c_coeffs_z
   end type spline_info

   public :: breakup

   !> Breakup object
   type :: breakup

      !> CCLabel objects for structure identification
      type(cclabel) :: ccl,ccl_film,ccl_lig

      class(ibconfig), pointer :: cfg
      type(vfs) , pointer :: vf
      type(tpns), pointer :: fs

      !> Drop transfer parameters
      logical  :: use_drop_transfer
      real(WP) :: dmax,dmin,ddel,emax
      real(WP) :: vof_tf_drop,vof_deleted,lper_d
      integer  :: np_drop

      !> Film transfer parameters
      logical  :: use_film_transfer,output_filmstats
      real(WP) :: frp,fmin,fd0,fbvol2dvol,fthld,minfdrop
      real(WP) :: fnumcell,vof_tf_film,fthc_percent
      integer  :: maxdpcell,np_film,film_tracker

      !> Ligament transfer parameters
      logical  :: use_lig_transfer
      real(WP) :: lmin,lmake,lper,lstratio,labsmin
      real(WP) :: dw,size_ratio,vof_tf_lig
      integer  :: np_lig

      !> Buffer layer size for auto-transfer
      integer  :: nlayer=4

      !< IDs for the droplet types
      !< Not due to boundary exit
      !< Direct conversion to droplets                                  : 1
      !< Conversion from film to droplets                               : 6
      !< Conversion from film to droplets (for special cases)           : 3
      !< Conversion from ligament to droplets                           : 11
      !< Due to boundary exit
      !< Direct conversion to droplets                                  : 2
      !< Conversion from film to droplets                               : 7
      !< Conversion from film to droplets (for special cases)           : 4
      !< Conversion from resolved ligament to droplets                  : 10
      !< Conversion from sub-grid ligament to droplets                  : 12

   contains
      procedure :: init            !< Initialize breakup model
      procedure :: transfer_drops  !< Transfer drops to Lagrangian
      procedure :: transfer_films  !< Transfer films to Lagrangian drops
      procedure :: transfer_ligs   !< Transfer ligaments to Lagrangian drops
      procedure :: fit_spline
      procedure :: distribute_on_spline
      ! procedure :: eigensolve_moi  !< Eigensolve moment of inertia
   end type breakup

   !> Hardcode inlet parameters (shared with atom_class)
   real(WP), parameter :: dl=0.003_WP

contains

   !> Initialize breakup model
   subroutine init(this,input,cfg,vf,fs)
      use filesys,  only: makedir,isdir
      use messager,  only: die
      implicit none
      character(len=str_medium) :: filename
      class(breakup), intent(inout) :: this
      class(inputfile), intent(inout) :: input
      class(ibconfig), target, intent(in) :: cfg
      class(vfs), target, intent(in) :: vf
      class(tpns), target, intent(in) :: fs
      integer :: iunit,ierr
      this%cfg=>cfg
      this%vf=>vf
      this%fs=>fs
      ! Read transfer flags
      call input%read('Transfer drops',this%use_drop_transfer,default=.true.)
      call input%read('Transfer films',this%use_film_transfer,default=.true.)
      call input%read('Transfer ligaments',this%use_lig_transfer,default=.true.)
      call input%read('Output film stats',this%output_filmstats,default=.false.)

      ! Initialize drop transfer
      if (this%use_drop_transfer) then
         call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
         this%ddel=0.2_WP*this%cfg%min_meshsize
         this%dmin=1.5_WP*this%cfg%min_meshsize
         this%dmax=1.0e-3_WP
         this%emax=0.8_WP
         this%lper_d=0.7_WP
         this%vof_tf_drop=0.0_WP
         this%vof_deleted=0.0_WP
         this%np_drop=0
      end if

      ! Initialize film transfer
      if (this%use_film_transfer) then
         call this%ccl_film%initialize(pg=this%cfg%pgrid,name='ccl_film')
         this%frp=0.0_WP
         this%fd0=dl
         this%fbvol2dvol=0.25_WP
         this%fmin=2.3e-6_WP
         this%fthld=10.0_WP*this%fmin
         this%minfdrop=0.1e-6_WP
         this%maxdpcell=50000
         this%fnumcell=50.0_WP
         this%fthc_percent=0.25_WP
         this%film_tracker=0
         this%vof_tf_film=0.0_WP
         this%np_film=0
      end if

      ! Initialize ligament transfer
      if (this%use_lig_transfer) then
         call this%ccl_lig%initialize(pg=this%cfg%pgrid,name='ccl_lig')
         this%dw=0.697_WP
         this%size_ratio=0.015_WP
         this%lmin=1.0_WP
         this%labsmin=6e-5_WP
         this%lmake=1.0_WP
         this%lper=0.9_WP
         this%lstratio=1.5_WP
         this%vof_tf_lig=0.0_WP
         this%np_lig=0
      end if

      if (this%cfg%amroot) then
         if (.not.isdir('spray-all')) call makedir('spray-all')     
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
      end if

   end subroutine init


   !> Transfer droplets to Lagrangian representation
   subroutine transfer_drops(this,lp,timet)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi,twoPi
      use messager,  only: die
      use vfs_class, only: VFlo,VFhi
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      class(lpt),         intent(inout) :: lp
      real(WP), intent(in) :: timet
      character(len=str_medium) :: filename
      real(WP), dimension(:)    , allocatable :: dvol,drem
      real(WP), dimension(:,:)  , allocatable :: dpos,dvel
      real(WP), dimension(:,:,:), allocatable :: dmoi
      integer :: n,m,ierr,i,j,k,iunit,l
      real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
      logical :: transfer
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      logical :: drem_active
      real(IRL_double), dimension(1:3) :: a_aligned_Cylinder
      type(spline_info) :: s_info
      real(WP):: Lrim,Vrim,minor_radius,Oh,tmpvol,lnum,lper
      real(WP), dimension(3) :: tmpxvol,tmpL
      integer  :: nmain,nsat,nneigh_moi,ii,jj,kk
      real(WP), dimension(:,:), allocatable :: points


      ! Build CCL
      call this%ccl%build(make_label,same_label)
      ! Allocate stats arrays
      allocate(dvol(1:this%ccl%nstruct        )); dvol=0.0_WP
      allocate(dpos(1:this%ccl%nstruct,1:3    )); dpos=0.0_WP
      allocate(dvel(1:this%ccl%nstruct,1:3    )); dvel=0.0_WP
      allocate(dmoi(1:this%ccl%nstruct,1:3,1:3)); dmoi=0.0_WP
      allocate(drem(1:this%ccl%nstruct        )); drem=0.0_WP

      ! First pass: accumulate volume, position, velocity
      do n=1,this%ccl%nstruct
         do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m)
            j=this%ccl%struct(n)%map(2,m)
            k=this%ccl%struct(n)%map(3,m)
            ! Get cell position, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL
            y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL
            z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL
            ! Accumulate volume, position, and velocity
            dvol(n  )=dvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            dpos(n,:)=dpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
            dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[sum(this%fs%U(i:i+1,j,k)),sum(this%fs%V(i,j:j+1,k)),sum(this%fs%W(i,j,k:k+1))]*0.5_WP
            ! Check if drop touches auto-transfer layer
            if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
            &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
            &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
            &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
            &   k.ge.this%vf%cfg%kmax-this%nlayer) drem(n)=1.0_WP
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dpos,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvel,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,drem,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

      ! Second pass: accumulate moment of inertia
      do n=1,this%ccl%nstruct
         x0=dpos(n,1)/dvol(n); y0=dpos(n,2)/dvol(n); z0=dpos(n,3)/dvol(n)
         do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
            ! Get cell position relative to drop barycenter, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL-x0
            y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL-y0
            z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL-z0
            ! Accumulate moment of inertia
            dmoi(n,1,1)=dmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
            dmoi(n,2,2)=dmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
            dmoi(n,3,3)=dmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
            dmoi(n,1,2)=dmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
            dmoi(n,1,3)=dmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
            dmoi(n,2,3)=dmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)

      ! Third pass: normalize
      do n=1,this%ccl%nstruct
         dpos(n,:)=dpos(n,:)/dvol(n)
         if (this%vf%cfg%xper.and.dpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) dpos(n,1)=dpos(n,1)+this%vf%cfg%xL
         if (this%vf%cfg%yper.and.dpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%vf%cfg%yL
         if (this%vf%cfg%zper.and.dpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%vf%cfg%zL
         dvel(n,:)=dvel(n,:)/dvol(n)
      end do

      ! Zero monitoring
      this%vof_tf_drop=0.0_WP; this%vof_deleted=0.0_WP; this%np_drop=0

      ! Transfer drops
      do n=1,this%ccl%nstruct
         diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)

         if (diam.gt.this%dmax) then
            transfer=.false.
         else if (diam.le.this%ddel) then
            transfer=.false.
            do m=1,this%ccl%struct(n)%n_
               this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            this%vof_deleted=this%vof_deleted+dvol(n)
         else if (diam.le.this%dmin) then
            transfer=.true.
         else
            A=dmoi(n,:,:)
            call eigensolve_moi(A,d)
            d=max(0.0_WP,d)
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
            ! if (lmin.eq.0.0_WP) lmin=lmid
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            ! transfer=(ecc.le.this%emax)
            if (d(3).lt.1.8*d(1)) then
               transfer=.true.
            else
               transfer=.false.
            end if
         end if

         drem_active=.false.
         if (drem(n).gt.0.0_WP) then
            transfer=.true.; drem_active=.true.
         end if

         ! Check if it might be a ligment if so transfer like a ligament, otherwise convert it based on a sphereical drop below
         if (drem_active) then
            checklig : block
               
               nneigh_moi=2; lper=0.0_WP; lnum=1.0_WP*this%ccl%struct(n)%n_
               ! Do ligament percentage check
               do m=1,this%ccl%struct(n)%n_
                  tmpvol=0.0_WP; tmpxvol=0.0_WP; A=0.0_WP
                  ! Get cell indices
                  i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
                  do kk = k-nneigh_moi,k+nneigh_moi
                     do jj = j-nneigh_moi,j+nneigh_moi
                        do ii = i-nneigh_moi,i+nneigh_moi
                           tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                           tmpxvol = tmpxvol + this%vf%Lbary(:,ii,jj,kk)*this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                        end do
                     end do
                  end do
                  ! Second pass to accumulate moment of inertia
                  if (tmpvol.gt.1.0e-14_WP) tmpxvol = tmpxvol/tmpvol
                  do kk = k-nneigh_moi,k+nneigh_moi
                     do jj = j-nneigh_moi,j+nneigh_moi
                        do ii = i-nneigh_moi,i+nneigh_moi
                           ! Location of film node
                           tmpL = this%vf%Lbary(:,ii,jj,kk) - tmpxvol
                           A(1,1)=A(1,1)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(2)**2+tmpL(3)**2)
                           A(2,2)=A(2,2)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(3)**2)
                           A(3,3)=A(3,3)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(2)**2)
                           A(1,2)=A(1,2)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(2)
                           A(1,3)=A(1,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(3)
                           A(2,3)=A(2,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(2)*tmpL(3)   
                        end do
                     end do
                  end do
                  call eigensolve_moi(A,d)
                  d=max(0.0_WP,d)
                  if ((d(3).lt.1.5_WP*d(2)).and.(d(2).gt.1.5_WP*d(1)) ) lper=lper+1.0_WP
               end do

               call MPI_ALLREDUCE(MPI_IN_PLACE,lnum,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,lper,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
               if (lnum .gt. 0.0_WP) then
                  lper=lper/lnum
               else
                  lper=0.0_WP
               end if

               if ((lper.gt.this%lper_d).or. ecc.gt.this%emax) then
                  transfer=.false.
                  call this%fit_spline(n,s_info,this%ccl)
                  if (s_info%flag.lt.1.0_WP) then
                     Lrim = s_info%length
                  else
                     Lrim=0.0_WP
                  end if
                  Vrim=dvol(n)
                  if (Lrim .le. VFlo .or. Vrim .le. VFlo) cycle
                  minor_radius = 0.0_WP
                  l = 0
                  do m=1,this%ccl%struct(n)%n_
                     ! Get cell indices
                     i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
                     if (this%vf%recon_type(i,j,k).eq.1) then
                        a_aligned_Cylinder = getAlignedCylinder(this%vf%liquid_gas_interface(i,j,k))
                        minor_radius = minor_radius + sqrt(a_aligned_Cylinder(1))
                        l=l+1
                     end if
                  end do
                  call MPI_ALLREDUCE(MPI_IN_PLACE,minor_radius,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
                  call MPI_ALLREDUCE(MPI_IN_PLACE,l,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
                  if (minor_radius.ge.VFlo) then
                     if (l .gt. 0) then
                        minor_radius=minor_radius/l
                     else
                        minor_radius=sqrt(Vrim/pi/Lrim)
                     end if
                  else
                     minor_radius=sqrt(Vrim/pi/Lrim)
                  end if
                  if (minor_radius .le. 0.8*sqrt(Vrim/pi/Lrim)) minor_radius = sqrt(Vrim/pi/Lrim)
                  if (minor_radius .le. VFlo) cycle
                  ! Drop size method from Kim & Moin (2020)
                  nmain=floor(this%dw*Lrim/(twoPi*minor_radius))

                  Oh=this%fs%visc_l/sqrt(this%fs%rho_l*minor_radius*this%fs%sigma)
                  this%size_ratio=1.0_WP/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
                  nsat=nmain+1; diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)

                  if (this%cfg%amRoot) then
                     filename='spray-all/droplets'
                     open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
                     if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
                     allocate(points(3,nsat+nmain))
                     if (s_info%flag.eq.0) call this%distribute_on_spline(nsat+nmain,s_info,points)
                     do l=1,nsat+nmain
                        ! Increment particle counter
                        lp%np_=lp%np_+1
                        ! Make room for new drop
                        call lp%resize(lp%np_)
                        ! Add the drop
                        lp%p(lp%np_)%id  =int(10,8)                                                                               
                        if (mod(l,2).eq.1) then
                           lp%p(lp%np_)%d=diam*this%size_ratio                                                                                    
                        else
                           lp%p(lp%np_)%d=diam                                                                                    
                        end if
                        if (s_info%flag.ne.0) then
                           lp%p(lp%np_)%pos = dpos(n,:)
                        else
                           lp%p(lp%np_)%pos = points(:,l)
                        end if
                        lp%p(lp%np_)%vel =dvel(n,:)
                        lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
                        lp%p(lp%np_)%flag=0                                                                                        
                        lp%p(lp%np_)%dt  =0.0_WP                                                                                  
                        lp%p(lp%np_)%Acol=0.0_WP                                                                                  
                        lp%p(lp%np_)%Tcol=0.0_WP
                        write(iunit,*)timet,lp%p(lp%np_)%d,lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3),&
                        &norm2([lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3)]),lp%p(lp%np_)%pos(1),&
                        &lp%p(lp%np_)%pos(2),lp%p(lp%np_)%pos(3),lp%p(lp%np_)%id  
                     end do
                     if (allocated(points)) deallocate(points)
                     close(iunit)
                     ! Increment monitoring variables
                     
                  end if
                  ! empty out the VF
                  do m=1,this%ccl%struct(n)%n_
                     this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
                  end do    

                  lp%np_new=lp%np_new+nmain+nsat
                  lp%vp_new=lp%vp_new+dvol(n)
                  this%np_lig=this%np_lig+nmain+nsat
                  this%vof_tf_lig=this%vof_tf_lig+dvol(n)
               else
                  transfer=.true.
               end if

            end block checklig 
         end if

         if (transfer) then
            ! Root creates new Lagrangian drop
            if (this%vf%cfg%amRoot) then
               lp%np_=lp%np_+1
               call lp%resize(lp%np_)
               if (drem_active) then
                  lp%p(lp%np_)%id  =int(2,8)
               else
                  lp%p(lp%np_)%id  =int(1,8)
               end if
               lp%p(lp%np_)%d   =diam
               lp%p(lp%np_)%pos =dpos(n,:)
               lp%p(lp%np_)%vel =dvel(n,:)
               lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(dpos(n,:),[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               lp%p(lp%np_)%flag=0
               lp%p(lp%np_)%dt  =0.0_WP
               lp%p(lp%np_)%Acol=0.0_WP
               lp%p(lp%np_)%Tcol=0.0_WP

               filename='spray-all/droplets'
               open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
               ! Output diameter, velocity, and position
               write(iunit,*) timet,lp%p(lp%np_)%d,lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3),&
               &norm2([lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3)]),lp%p(lp%np_)%pos(1),&
               &lp%p(lp%np_)%pos(2),lp%p(lp%np_)%pos(3),lp%p(lp%np_)%id  
               ! Close the file
               close(iunit)
            end if

            do m=1,this%ccl%struct(n)%n_
               this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            this%vof_tf_drop=this%vof_tf_drop+dvol(n)
            this%np_drop=this%np_drop+1
            lp%np_new=lp%np_new+1
            lp%vp_new=lp%vp_new+dvol(n)
         end if
      end do

      ! Synchronize
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      call lp%sync()

      deallocate(dvol,dpos,dvel,dmoi,drem)

   contains
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
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
   end subroutine transfer_drops


   !> Transfer films to Lagrangian drops (Jackiw & Ashgriz model)
   subroutine transfer_films(this,lp,timet)
      use irl_fortran_interface
      use vfs_class, only: VFlo,VFhi
      use mathtools, only: Pi,normalize,cross_product
      use random,    only: random_uniform,random_gamma
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use messager,  only: die
      implicit none
      class(breakup), intent(inout) :: this
      class(lpt), intent(inout) :: lp
      real(WP), intent(in) :: timet
      character(len=str_medium) :: filename
      real(WP), dimension(:), allocatable :: fvol,fthc,frem,fcnt,fcurv,fcsa,fthick_,fthick,fivol_,fivol
      real(WP), dimension(:), allocatable :: sort_ke,fcsar_,fcsar,fcurvr,fcurvr_
      real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
      integer, dimension(:), allocatable ::  sort_id,plist,dispels
      integer :: n,nn,m,i,j,k,ii,jj,kk,ncell_,tmp_id,l,totalnewp,np_start,np_old,count,ierr,ind,ip,iunit,rank
      integer :: totalcell,nfthcpercent,dropcounter
      real(WP) :: tmp_ke,Vt,Vl,Vd,alpha,beta,mySA,tmp_thic,dh,t1,t2
      real(WP), dimension(3) :: nref,tref,sref,film_tangent,grad_h,r_nb
      logical :: sampled, frem_active
      ! Build CCL
      call this%ccl_film%build(make_label,same_label)
      if (this%ccl_film%nstruct.ge.1) then
         ! Allocate film stats arrays
         allocate(fvol (1:this%ccl_film%nstruct)); fvol=0.0_WP
         allocate(fthc (1:this%ccl_film%nstruct)); fthc=10.0_WP*this%cfg%min_meshsize
         allocate(frem (1:this%ccl_film%nstruct)); frem=0.0_WP
         allocate(fcnt (1:this%ccl_film%nstruct)); fcnt=0.0_WP
         allocate(fcurv(1:this%ccl_film%nstruct)); fcurv=0.0_WP
         allocate(fcsa (1:this%ccl_film%nstruct)); fcsa=0.0_WP
         allocate(plist  (0:this%vf%cfg%nproc-1))
         allocate(dispels(0:this%vf%cfg%nproc-1))
         ! Get local thickness of the film to determine if film should be convereted
         call this%vf%get_thickness()
         ! First pass to accumulate volume and get minimum thickness
         do n=1,this%ccl_film%nstruct
            fcnt(n)=1.0_WP*this%ccl_film%struct(n)%n_
            ! For each structure n, individually acquire all the thickness value
            call MPI_AllGATHER(this%ccl_film%struct(n)%n_,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%vf%cfg%comm,ierr)
            totalcell=sum(plist)
            allocate(fthick_(1:this%ccl_film%struct(n)%n_))
            allocate(fthick (1:totalcell))
            count=1
            ! Loop over cells in structure
            do m=1,this%ccl_film%struct(n)%n_
               ! Get cell indices
               i=this%ccl_film%struct(n)%map(1,m); j=this%ccl_film%struct(n)%map(2,m); k=this%ccl_film%struct(n)%map(3,m)
               ! Accumulate liquid volume
               fvol(n)=fvol(n)+this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
               ! Record all thickness and bag volume
               fthick_(count)=this%vf%thickness(i,j,k)
               count=count+1
               ! Sum curvatures
               do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).eq.0) cycle
                  mySA=abs(calculateVolume(this%vf%interface_polygon(l,i,j,k)))
                  fcurv(n)=fcurv(n)+abs(this%vf%curv2p(l,i,j,k))*mySA
                  fcsa(n)=fcsa(n)+mySA
               end do
               ! Check if film touches auto burst layer
               if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
               &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
               &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
               &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
               &   k.ge.this%vf%cfg%kmax-this%nlayer) frem(n)=1.0_WP
            end do
            ! Get the thickness of the 15th percentile
            count = 0
            do rank=0,this%vf%cfg%nproc-1
               dispels(rank) = count
               count = count + plist(rank)
            end do
            call MPI_ALLGATHERV(fthick_,this%ccl_film%struct(n)%n_,MPI_REAL_WP,fthick,plist,dispels,MPI_REAL_WP,this%vf%cfg%comm)
            nfthcpercent  = max(1, int(this%fthc_percent*totalcell))
            if (totalcell > 1) call quicksort_real(fthick, 1, totalcell)
            fthc(n)=fthick(nfthcpercent)
            deallocate(fthick_,fthick)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,fvol ,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,frem ,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,fcnt ,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,fcurv,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,fcsa ,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         ! Zero out monitoring variables
         this%vof_tf_film=0.0_WP
         this%np_film=0
         ! Record initial droplets in each processor for future outputing purpose
         np_start=lp%np_; sampled=.false.
         ! Second pass to decide if the film has reached a minimum thickness to burst
         do n=1,this%ccl_film%nstruct
            if (fcurv(n).gt.0.0_WP .and. fcsa(n).gt.0.0_WP) then
               fcurv(n)=fcurv(n)/fcsa(n)
               if (fcurv(n).gt.1.0_WP/this%cfg%min_meshsize) fcurv(n)=1.0_WP/this%cfg%min_meshsize
            else
               ! Set it to be the largest curvature
               fcurv(n)=1.0_WP/this%cfg%min_meshsize
            end if
            ! Min thickness below threshold and film volume greater than a threshold 
            frem_active = .false.
            if (fthc(n).le.this%fmin .and. fcnt(n).gt.this%fnumcell) then
            else if (frem(n).gt.0.0_WP) then
               ! Too close to the end of domain
               frem_active = .true.
            else
               cycle
            end if
            ! output to confirm
            if (this%vf%cfg%amRoot) print *, "This is a thin film with min_thickness", fthc(n), "and this is id:", n ,"vol is:", fvol(n),"cur for the bag: ",fcurv(n),'number of cells',fcnt(n) ,'end of domain: ',frem_active

            if (this%output_filmstats) then
               ! Now output the local thickness and volume to files for a prior testing
               output_film: block
                  use filesys,  only: makedir,isdir
                  call MPI_AllGATHER(this%ccl_film%struct(n)%n_,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%vf%cfg%comm,ierr)
                  totalcell=sum(plist)
                  allocate(fthick_(1:this%ccl_film%struct(n)%n_));allocate(fivol_(1:this%ccl_film%struct(n)%n_))
                  allocate(fthick (1:totalcell)); allocate(fivol (1:totalcell))
                  allocate(fcsar_(1:this%ccl_film%struct(n)%n_));allocate(fcurvr_(1:this%ccl_film%struct(n)%n_))
                  allocate(fcsar (1:totalcell)); allocate(fcurvr (1:totalcell))
                  count=1
                  do m=1,this%ccl_film%struct(n)%n_
                     ! Get cell indices
                     i=this%ccl_film%struct(n)%map(1,m); j=this%ccl_film%struct(n)%map(2,m); k=this%ccl_film%struct(n)%map(3,m)
                     ! Record all thickness and bag volume
                     fthick_(count)=this%vf%thickness(i,j,k)
                     fivol_(count)=this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
                     do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).eq.0) cycle
                        mySA=abs(calculateVolume(this%vf%interface_polygon(l,i,j,k)))
                        fcsar_(count)=fcsar_(count)+mySA
                        fcurvr_(count)=fcurvr_(count)+abs(this%vf%curv2p(l,i,j,k))*mySA                     
                     end do
                     if (fcsar_(count).gt.0.0_WP) then
                        fcurvr_(count)=fcurvr_(count)/fcsar_(count)
                     else
                        fcurvr_(count)=0.0_WP
                     end if
                     count=count+1
                  end do
                  ! Accumulate the film information across all processors
                  count = 0
                  do rank=0,this%vf%cfg%nproc-1
                     dispels(rank) = count
                     count = count + plist(rank)
                  end do
                  call MPI_ALLGATHERV(fthick_,this%ccl_film%struct(n)%n_,MPI_REAL_WP,fthick,plist,dispels,MPI_REAL_WP,this%vf%cfg%comm)
                  call MPI_ALLGATHERV(fivol_,this%ccl_film%struct(n)%n_,MPI_REAL_WP,fivol,plist,dispels,MPI_REAL_WP,this%vf%cfg%comm)
                  call MPI_ALLGATHERV(fcsar_,this%ccl_film%struct(n)%n_,MPI_REAL_WP,fcsar,plist,dispels,MPI_REAL_WP,this%vf%cfg%comm)
                  call MPI_ALLGATHERV(fcurvr_,this%ccl_film%struct(n)%n_,MPI_REAL_WP,fcurvr,plist,dispels,MPI_REAL_WP,this%vf%cfg%comm)
                  if (this%vf%cfg%amRoot)  then
                     if (.not.isdir('film')) call makedir('film')
                     this%film_tracker=this%film_tracker+1
                     write(filename,'("film/film_data_",I0)') this%film_tracker
                     open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',position='append',iostat=ierr)
                     if (ierr.ne.0) call die('[Write film stats] Could not open file: '//trim(filename))
                     do i = 1,totalcell
                        write(iunit,*) fthick(i),fivol(i),fcsar(i),fcurvr(i)
                     end do
                     close(iunit)
                  end if
                  deallocate(fthick,fivol,fthick_,fivol_,fcsar,fcurvr,fcsar_,fcurvr_)
               end block  output_film
            end if
            ! Assume fd0 across the processor based on the total volume of the film (doesn't really matter to the droplet size generation in the end)
            if (.not.frem_active) then
               this%fd0 =(6.0_WP*fvol(n)/this%fbvol2dvol/Pi)**(1.0_WP/3.0_WP)
            else
               this%fd0 =dl
            end if
            np_old=lp%np_
            ! sort cell index based on local film thickness
            if (this%ccl_film%struct(n)%n_.ge.1) then
               allocate(sort_id(1:this%ccl_film%struct(n)%n_)) 
               allocate(sort_ke(1:this%ccl_film%struct(n)%n_))
               ncell_ = this%ccl_film%struct(n)%n_
               do m=1,ncell_
                  i=this%ccl_film%struct(n)%map(1,m)
                  j=this%ccl_film%struct(n)%map(2,m)
                  k=this%ccl_film%struct(n)%map(3,m)
                  sort_id(m)=m 
                  sort_ke(m)=this%vf%thickness(i,j,k)
               end do
               ! sort based on thickness
               if (ncell_ > 1) call quicksort_real(sort_ke, 1, ncell_, sort_id)
               Vt=0.0_WP; Vl=0.0_WP
               do m=1,ncell_
                  i=this%ccl_film%struct(n)%map(1,sort_id(m)); j=this%ccl_film%struct(n)%map(2,sort_id(m)); k=this%ccl_film%struct(n)%map(3,sort_id(m))
                  ! Accumulate 
                  Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
                  dropcounter=0
                  do while (.true.)
                     if (.not.sampled) then
                        call bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP/fcurv(n))
                        Vd = pi/6.0_WP*(max(min(random_gamma(alpha)*beta*this%fd0,2.0_WP*this%frp),this%minfdrop))**3
                        sampled = .true.
                     end if
                     if (Vl.gt.Vd) then
                        nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
                        select case (maxloc(abs(nref),1))
                        case (1)
                           tref=normalize([+nref(2),-nref(1),0.0_WP])
                        case (2)
                           tref=normalize([0.0_WP,+nref(3),-nref(2)])
                        case (3)
                           tref=normalize([-nref(3),0.0_WP,+nref(1)])
                        end select
                        sref=cross_product(nref,tref)
                        grad_h=0.0_WP
                        do kk = -1, 1
                           do jj = -1, 1
                              do ii = -1, 1
                                 if (ii == 0 .and. jj == 0 .and. kk == 0) cycle
                                 if (this%vf%VF(i+ii, j+jj, k+kk).gt. VFlo)then
                                    ! Difference in thickness (positive if neighbor is thicker)
                                    dh = this%vf%thickness(i+ii,j+jj,k+kk) - this%vf%thickness(i,j,k)
                                    r_nb = this%vf%Lbary(:,i+ii,j+jj,k+kk) - this%vf%Lbary(:,i,j,k)
                                    grad_h = grad_h + (dh * r_nb)
                                 end if
                              end do
                           end do
                        end do
                        grad_h = normalize(grad_h)
                        film_tangent = dot_product(grad_h,tref)*tref + dot_product(grad_h,sref)*sref
                        film_tangent = normalize(film_tangent)
                        ! Increment particle counter
                        lp%np_=lp%np_+1; call lp%resize(lp%np_)
                        ! Add the drop
                        if (frem_active) then
                           lp%p(lp%np_)%id  =int(7,8)
                        else                                   
                           lp%p(lp%np_)%id  =int(6,8)
                        end if
                        lp%p(lp%np_)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            
                        lp%p(lp%np_)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref
                        lp%p(lp%np_)%vel =film_tangent!this%fs%cfg%get_velocity(pos=lp%p(lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                        lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])    !< Place the drop in the proper cell for the lp%cfg
                        lp%p(lp%np_)%flag=0                                          
                        lp%p(lp%np_)%dt  =0.0_WP                                     
                        lp%p(lp%np_)%Acol=0.0_WP                                     
                        lp%p(lp%np_)%Tcol=0.0_WP  

                        ! Update tracked volumes
                        Vl=Vl-Vd
                        Vt=Vt+Vd
                        sampled = .false.

                        ! Increment monitoring variables
                        this%vof_tf_film=this%vof_tf_film+Vd
                        this%np_film=this%np_film+1
                        lp%np_new=lp%np_new+1
                        lp%vp_new=lp%vp_new+Vd

                        dropcounter=dropcounter+1
                     else
                        exit
                     end if
                     if (dropcounter.gt.this%maxdpcell) exit
                  end do
                  ! Remove liquid in that cell
                  this%vf%VF(i,j,k)=0.0_WP
               end do


               deallocate(sort_id,sort_ke)
               ! If for some reason a film with 0 liquid volume has been tagged, skip it
               if (Vt.eq.0.0_WP .and. Vl.eq.0.0_WP) cycle
               ! Based on how many particles were created, decide what to do with left-over volume
               if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                  ! Increment particle counter
                  lp%np_=lp%np_+1
                  ! Make room for new drop
                  call lp%resize(lp%np_)
                  ! Add the drop
                  if (frem_active) then
                     lp%p(lp%np_)%id  =int(4,8)
                  else                                   
                     lp%p(lp%np_)%id  =int(3,8)
                  end if                                   
                  lp%p(lp%np_)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            
                  lp%p(lp%np_)%pos =this%vf%Lbary(:,i,j,k)                     
                  lp%p(lp%np_)%vel =0.0_WP!this%fs%cfg%get_velocity(pos=lp%p(lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
                  lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  lp%p(lp%np_)%flag=0                                          
                  lp%p(lp%np_)%dt  =0.0_WP                                     
                  lp%p(lp%np_)%Acol=0.0_WP                                    
                  lp%p(lp%np_)%Tcol=0.0_WP
                  ! Increment monitoring variables
                  lp%np_new=lp%np_new+1
                  this%np_film=this%np_film+1
               else ! Some particles were created, make them all larger
                  do ip=np_old+1,lp%np_
                     lp%p(ip)%d=lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                  end do
               end if
               ! Increment monitoring variables
               this%vof_tf_film=this%vof_tf_film+Vl
               lp%vp_new=lp%vp_new+Vl
            end if

            ! At this point the film has been converted, we therefore need to apply a kick
            applykick : block
               real(WP) :: sumdsq,ckick
               integer, dimension(3) :: atom_ind
               ! Get summed diameter squared
               sumdsq=0.0_WP
               do ip = np_old+1, lp%np_
                  sumdsq=sumdsq+lp%p(ip)%d**2 
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,sumdsq,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
               ! This is dervied from surface energy difference: sigma*fsca(n)-sum pi d^2
               ! The velocity kick magnitude of each droplet gets is sum 1/2 rho_l (1/6) pi d^3 * (ckick^2)/d
               ckick=sqrt(12.0_WP*this%fs%sigma*(max(fcsa(n)-Pi*sumdsq,0.0_WP))/(Pi*this%fs%rho_l*sumdsq))
               do ip = np_old+1, lp%np_
                  ! Get the droplet location on the atomization mesh
                  atom_ind=this%cfg%get_ijk_global(lp%p(ip)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  lp%p(ip)%vel=lp%p(ip)%vel*ckick/sqrt(lp%p(ip)%d) + &
                  & this%fs%cfg%get_velocity(pos=lp%p(ip)%pos,i0=atom_ind(1),j0=atom_ind(2),k0=atom_ind(3),U=this%fs%U,V=this%fs%V,W=this%fs%W)           
               end do
            end block applykick
         end do

         ! Gather the number of newly generated particles from each processor due to film burst
         totalnewp = 0
         ! Get number of particle generated for each processor
         call MPI_AllGATHER(this%np_film,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%vf%cfg%comm,ierr)
         totalnewp= sum(plist)
         ! If there is any particle generated
         if (totalnewp .gt. 0) then
            allocate(pinfo_(1:9,1:this%np_film))
            allocate(pinfo(1:9,1:totalnewp))
            ! Get info
            do ip = np_start+1, lp%np_
               pinfo_(1,ip-np_start)=lp%p(ip)%d
               pinfo_(2:4,ip-np_start)=lp%p(ip)%vel
               pinfo_(5,ip-np_start)=norm2(lp%p(ip)%vel)
               pinfo_(6:8,ip-np_start)=lp%p(ip)%pos
               pinfo_(9,ip-np_start)=lp%p(ip)%id
            end do
            ! Calculate dispels
            count = 0
            do rank=0,this%vf%cfg%nproc-1
               dispels(rank) = count
               count = count + plist(rank)
            end do
            ! Communicate to root
            do i = 1,9
               call MPI_GATHERV(pinfo_(i,:),this%np_film,MPI_REAL_WP,pinfo(i,:),plist,dispels,MPI_REAL_WP,0,this%vf%cfg%comm)
            end do
            !!! Write to droplet list !!!
            if (this%vf%cfg%amRoot)  then
               filename='spray-all/droplets'
               open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
               do i = 1,totalnewp
               write(iunit,*) timet,pinfo(1,i),pinfo(2,i),pinfo(3,i),pinfo(4,i),pinfo(5,i),pinfo(6,i),pinfo(7,i),pinfo(8,i),INT(pinfo(9,i))
               end do
               close(iunit)
            end if
            deallocate(pinfo_,pinfo)
            ! Synchronize VF fields
            call this%vf%cfg%sync(this%vf%VF)
            call this%vf%clean_irl_and_band()
            ! Synchronize particles
            call lp%sync()
            ! Integrate monitoring variables 
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_film,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_film    ,1,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         end if
         deallocate(fvol,fthc,frem,fcnt,fcurv,fcsa,plist,dispels)
      end if 

      contains
         subroutine bag_droplet_gamma(h,R)
            implicit none
            real(WP), intent(in) :: h,R
            real(WP) :: Utc,ac,b,dr,ds,Oh
            real(WP) :: mean, stdev
            real(WP) :: mult_fact,h_use
            if (h.le.this%fthld) then
               mult_fact=2.0_WP
               if (h.le.this%fmin) then
                  h_use=this%fmin
               else
                  h_use=h
               end if
               mean=mult_fact*(this%fmin+h_use+this%fthld)/(3.0_WP*this%fd0)
               stdev=sqrt(sum((mult_fact*[this%fmin,h_use,this%fthld]/this%fd0-mean)**2)/(3.0_WP))
               ! Gamma distribution parameters
               alpha=(mean/stdev)**2
               beta=stdev**2/mean
               ! Effectively remove upperbound
               this%frp=100.0_WP*this%fmin
            else
               ! Retraction speed
               Utc=sqrt(2.0_WP*this%fs%sigma/this%fs%rho_l/h)
               ! Centripetal acceleration
               ac=Utc**2/R
               ! Rim diameter
               b=sqrt(this%fs%sigma/this%fs%rho_l/ac)
               ! RP droplet diameter
               this%frp=1.508_WP*b
               ! Rim Ohnesorge number
               Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b*this%fs%sigma*0.5_WP)
               ! Satellite droplet diameter
               ds=this%frp/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
               ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
               mean=0.25_WP*(h+b+this%frp+ds)/this%fd0
               stdev=sqrt(0.25_WP*sum(([h,b,this%frp,ds]/this%fd0-mean)**2))
               ! Gamma distribution parameters
               alpha=(mean/stdev)**2
               beta=stdev**2/mean
            end if
         end subroutine bag_droplet_gamma
      
         !> Function that identifies cells that need a label
         logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         ! if ((this%vf%VF(i,j,k).gt.VFlo).and.(this%vf%VF(i,j,k).lt.VFhi).and.((this%vf%norm_pos(i,j,k)-this%vf%norm_neg(i,j,k)).lt.0.5_WP).and.((this%vf%norm_pos(i,j,k)+this%vf%norm_neg(i,j,k)).ge.0.925_WP)) then
            if ((this%vf%VF(i,j,k).gt.VFlo).and.(this%vf%VF(i,j,k).lt.VFhi).and.this%vf%thin_sensor(i,j,k).eq.1.0_WP)then
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
      
   end subroutine transfer_films



   subroutine transfer_ligs(this,lp,timet)
      use vfs_class, only: VFlo,VFhi
      use mathtools, only: pi,twoPi
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use messager, only: die
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      class(lpt), intent(inout) :: lp
      real(WP), intent(in) :: timet
      character(len=str_medium) :: filename
      real(WP), dimension(:)    , allocatable :: lvol,lthc,llen,lnum,lper
      real(WP), dimension(:)    , allocatable :: lrem,lSR,xmin,xmax,ymin,ymax,zmin,zmax
      real(WP), dimension(:,:)  , allocatable :: lpos,lvel
      real(WP), dimension(:,:,:), allocatable :: lmoi
      integer :: n,m,ierr,i,j,k,l,ii,jj,kk,iunit,totalnewp,count,ip,rank
      real(WP) :: x,y,z,x0,y0,z0,lmax,lmid,lmin
      real(WP) :: Vt,Vl,Vd,minor_radius,diam,Vrim,Lrim
      real(WP) :: Trp,Lrp,Tsr,SR_tmp,Oh,b
      real(IRL_double), dimension(1:3) :: a_aligned_Cylinder
      real(WP) :: prob
      real(WP), dimension(1:3) :: tangent
      real(WP), dimension(:,:,:,:), allocatable :: SR
      integer  :: nmain,nsat
      real(WP), dimension(:,:,:), allocatable :: thickness
      integer,  dimension(:,:,:), allocatable :: struct_type
      real(WP), dimension(:,:), allocatable :: points
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      logical :: lrem_active

      ! type :: spline_info
      !    integer :: n_knots_x, n_knots_y, n_knots_z, flag
      !    real(WP) :: length
      !    real(WP), dimension(:), allocatable :: t_knots_x, c_coeffs_x
      !    real(WP), dimension(:), allocatable :: t_knots_y, c_coeffs_y
      !    real(WP), dimension(:), allocatable :: t_knots_z, c_coeffs_z
      ! end type spline_info
   
      type(spline_info) :: s_info

      ! Get thickness and local struct_type for global information calculation
      allocate(thickness  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));thickness=0.0_WP
      allocate(struct_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));struct_type=0
      call get_liginfo()
      ! Start by performing a CCL based on ligament criteria
      call this%ccl_lig%build(make_label,same_label)
      if (this%ccl_lig%nstruct.ge.1) then
         ! Allocate ligament stats arrays
         allocate(lvol(1:this%ccl_lig%nstruct        )); lvol=0.0_WP
         allocate(lthc(1:this%ccl_lig%nstruct        )); lthc=HUGE(x)
         allocate(llen(1:this%ccl_lig%nstruct        )); llen=0.0_WP
         allocate(lnum(1:this%ccl_lig%nstruct        )); lnum=0.0_WP
         allocate(lper(1:this%ccl_lig%nstruct        )); lper=0.0_WP
         allocate(lpos(1:this%ccl_lig%nstruct,1:3    )); lpos=0.0_WP
         allocate(lvel(1:this%ccl_lig%nstruct,1:3    )); lvel=0.0_WP
         allocate(lmoi(1:this%ccl_lig%nstruct,1:3,1:3)); lmoi=0.0_WP
         allocate(lrem(1:this%ccl_lig%nstruct        )); lrem=0.0_WP
         allocate(lSR(1:this%ccl_lig%nstruct        )); lSR=-HUGE(x)
         allocate(xmin(1:this%ccl_lig%nstruct),xmax(1:this%ccl_lig%nstruct)); xmin=HUGE(x);xmax=-HUGE(x)
         allocate(ymin(1:this%ccl_lig%nstruct),ymax(1:this%ccl_lig%nstruct)); ymin=HUGE(x);ymax=-HUGE(x)
         allocate(zmin(1:this%ccl_lig%nstruct),zmax(1:this%ccl_lig%nstruct)); zmin=HUGE(x);zmax=-HUGE(x)
         allocate(SR(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));SR=0.0_WP
         call this%fs%get_strainrate(SR)
         call this%cfg%sync(this%vf%SD)
         ! First pass to accumulate volume, position, min thickness and ligament percentage
         do n=1,this%ccl_lig%nstruct
            ! Loop over cells in structure
            lnum(n)=lnum(n)+1.0_WP*this%ccl_lig%struct(n)%n_
            do m=1,this%ccl_lig%struct(n)%n_
               ! Get cell indices
               i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
               ! Get cell position, accounting for periodicity
               x=this%vf%cfg%xm(i)-this%ccl_lig%struct(n)%per(1)*this%vf%cfg%xL
               y=this%vf%cfg%ym(j)-this%ccl_lig%struct(n)%per(2)*this%vf%cfg%yL
               z=this%vf%cfg%zm(k)-this%ccl_lig%struct(n)%per(3)*this%vf%cfg%zL
               ! Accumulate volume and position. Get min thickness and ligament percentage
               lvol(n  )=lvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
               lpos(n,:)=lpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
               lvel(n,:)=lvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[sum(this%fs%U(i:i+1,j,k)),sum(this%fs%V(i,j:j+1,k)),sum(this%fs%W(i,j,k:k+1))]*0.5_WP
               if (struct_type(i,j,k).eq.1) then
                  lthc(n)=min(lthc(n),thickness(i,j,k))
                  lper(n)=lper(n)+1.0_WP
               end if
               ! Check if drop touches auto-transfer layer
               if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
               &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
               &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
               &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
               &   k.ge.this%vf%cfg%kmax-this%nlayer) lrem(n)=1.0_WP
               ! Get the structures's locally largest and smallest x,y,z locations
               if (this%vf%recon_type(i,j,k).eq.2) then
                  do l=1,2
                     if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                        d = calculateCentroid(this%vf%interface_polygon(l,i,j,k))
                        xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
                        ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
                        zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
                     end if
                  end do
               else if (this%vf%recon_type(i,j,k).eq.1) then
                  if (this%vf%VF(i,j,k).ge.VFlo) then
                     d = this%vf%LBARY(:,i,j,k)
                     xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
                     ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
                     zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
                  end if
               end if
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lvol,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lpos,3*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lvel,3*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lrem,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lthc,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lnum,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,lper,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
            
         ! Second pass to accumulate moment of inertia
         do n=1,this%ccl_lig%nstruct
            ! Get ligament barycenter
            x0=lpos(n,1)/lvol(n); y0=lpos(n,2)/lvol(n); z0=lpos(n,3)/lvol(n)
            ! Loop over cells in structure
            do m=1,this%ccl_lig%struct(n)%n_
               ! Get cell indices
               i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
               ! Get cell position relative to drop barycenter, accounting for periodicity
               x=this%vf%cfg%xm(i)-this%ccl_lig%struct(n)%per(1)*this%vf%cfg%xL-x0
               y=this%vf%cfg%ym(j)-this%ccl_lig%struct(n)%per(2)*this%vf%cfg%yL-y0
               z=this%vf%cfg%zm(k)-this%ccl_lig%struct(n)%per(3)*this%vf%cfg%zL-z0
               ! Accumulate moment of inertia
               lmoi(n,1,1)=lmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
               lmoi(n,2,2)=lmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
               lmoi(n,3,3)=lmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
               lmoi(n,1,2)=lmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
               lmoi(n,1,3)=lmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
               lmoi(n,2,3)=lmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,lmoi,9*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)

         ! Third pass to generalize ligament stats
         do n=1,this%ccl_lig%nstruct
            ! Get ligament, accounting for periodicity
            lpos(n,:)=lpos(n,:)/lvol(n)
            if (this%vf%cfg%xper.and.lpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) lpos(n,1)=lpos(n,1)+this%vf%cfg%xL
            if (this%vf%cfg%yper.and.lpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) lpos(n,2)=lpos(n,2)+this%vf%cfg%yL
            if (this%vf%cfg%zper.and.lpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) lpos(n,3)=lpos(n,3)+this%vf%cfg%zL
            ! Get drop velocity
            lvel(n,:)=lvel(n,:)/lvol(n)
            ! Calculate the percentage of ligament structure type
            lper(n)=lper(n)/lnum(n)
            ! Calculate maximum length of the structure
            A=lmoi(n,:,:)
            call eigensolve_moi(A,d)
            d=max(0.0_WP,d)    
            ! Replace with corrected eigenvectors for future ligament droplet placement
            lmoi(n,:,:)=A
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/lvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/lvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/lvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ! Use max of bounding box and MoI-derived lengths as length
            llen(n) = max(sqrt((xmax(n)-xmin(n))**2+(ymax(n)-ymin(n))**2+(zmax(n)-zmin(n))**2),lmax)
            ! With the tangent direction of the ligament, we can evaluate the strain rate of each cell of the ligament
            tangent = lmoi(n,:,1)
            do m=1,this%ccl_lig%struct(n)%n_
               ! Get cell indices
               i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
               SR_tmp =SR(1,i,j,k)*tangent(1)**2        +SR(2,i,j,k)*tangent(2)**2        +SR(3,i,j,k)*tangent(3)**2 + &
             & 2.0_WP*(SR(4,i,j,k)*tangent(1)*tangent(2)+SR(5,i,j,k)*tangent(2)*tangent(3)+SR(6,i,j,k)*tangent(1)*tangent(3))
               lSR(n) = max(lSR(n),abs(SR_tmp))
            end do
         end do
         ! Find the maximum tangential strain rate of each ligament
         call MPI_ALLREDUCE(MPI_IN_PLACE,lSR,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

         ! Zero out monitoring variables
         this%vof_tf_lig=0.0_WP; this%np_lig=0
         ! Perform transfer
         do n=1,this%ccl_lig%nstruct
            call this%fit_spline(n,s_info,this%ccl_lig)
            
            if (s_info%flag.lt.1.0_WP) then
               Lrim = s_info%length
            else
               Lrim=0.0_WP!llen(n)
            end if
            Vrim=lvol(n)
            if (Lrim .le. VFlo .or. Vrim .le. VFlo) cycle
            minor_radius = 0.0_WP
            l = 0.0_WP
            do m=1,this%ccl_lig%struct(n)%n_
               ! Get cell indices
               i=this%ccl_lig%struct(n)%map(1,m)
               j=this%ccl_lig%struct(n)%map(2,m)
               k=this%ccl_lig%struct(n)%map(3,m)
               if (this%vf%recon_type(i,j,k).eq.1) then
                  a_aligned_Cylinder = getAlignedCylinder(this%vf%liquid_gas_interface(i,j,k))
                  minor_radius = minor_radius + sqrt(a_aligned_Cylinder(1))
                  l=l+1.0_WP
               end if
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,minor_radius,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,l,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
            if (minor_radius.ge.VFlo) then
               if (l .gt. 0) then
                  minor_radius=minor_radius/l
               else
                  minor_radius=sqrt(Vrim/pi/Lrim)
               end if
            else
               minor_radius=sqrt(Vrim/pi/Lrim)
            end if
            if (minor_radius .le. 0.8*sqrt(Vrim/pi/Lrim)) minor_radius = sqrt(Vrim/pi/Lrim)
            if (minor_radius .le. VFlo) cycle
            ! Drop size method from Kim & Moin (2020)
            nmain=floor(this%dw*Lrim/(twoPi*minor_radius))
            ! Calculate breakup time scale based on inviscid RP instability analysis
            Trp=2.91258_WP*sqrt(this%fs%rho_l*minor_radius**3/this%fs%sigma)
            ! Calcuate time scale based on maximum local strainrate
            Tsr=1.0_WP/lSR(n)
            ! Only breakup if minimum thickness is reached, sufficient volume of the ligament, enough of local ligament-like structures,
            ! local time scale asscoiated with strain rate is on par or bigger than the RP time scale, and its length is longer than the inviscid most unstable wavelength
            ! lrem_active = .false.
            ! if ((lthc(n).le.this%lmin*this%cfg%min_meshsize).and.(lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.lper).and.(Trp.le.Tsr).and.(nmain.ge.1)) then
            ! else if(lrem(n).gt.0.0_WP.and.(lvol(n).ge.this%cfg%min_meshsize**3)) then
            !    lrem_active=.true.
            ! else
            !    cycle
            ! end if
            
            call random_number(prob)
            prob = prob * lthc(n)/this%cfg%min_meshsize
            lrem_active = .false.
            ! if(lrem(n).gt.0.0_WP.and.(lvol(n).ge.this%cfg%min_meshsize**3)) then
            !    lrem_active=.true.
            ! else if (at_edge_total.or.(.not.((prob.le.0.01_WP).and.(Lrim.gt.4*minor_radius).and.(Vrim.ge.0.01*this%cfg%min_meshsize**3).and.(lnum(n).ge.10).and.(minor_radius.ge.0.01*this%cfg%min_meshsize).and.(nmain.ge.1).and.(lthc(n).le.2*this%cfg%min_meshsize).and.(lper(n).gt.0.9)))) then
            !    cycle
            ! end if
            ! if ((.not.at_edge_total).or.(.not.((prob.le.0.01_WP).and.(Lrim.gt.4*minor_radius).and.(Vrim.ge.0.01*this%cfg%min_meshsize**3).and.(lnum(n).ge.10).and.(minor_radius.ge.0.01*this%cfg%min_meshsize).and.(nmain.ge.1).and.(lthc(n).le.2*this%cfg%min_meshsize).and.(lper(n).gt.0.9)))) cycle!.and.(Trp.le.dt)!.and.(Trp.le.Tsr).and.(lig_timers_new(n,1).ge.Trp)
            ! if ((lthc(n).le.this%labsmin).and.(lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.lper).and.(nmain.ge.1)) then

            if (((lthc(n).le.this%labsmin).and.(lper(n).ge.this%lper)).or.((lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.this%lper).and.(Trp.le.Tsr).and.(nmain.ge.1))) then
            else if(lrem(n).gt.0.0_WP.and.(lvol(n).ge.this%cfg%min_meshsize**3)) then
               lrem_active=.true.
            else
               cycle
            end if
            
            ! if (this%vf%cfg%amRoot) print *, "This is the min_thickness", lthc(n), ",lig percentage:", lper(n),"max length:",llen(n),&
            ! & "how many cells",lnum(n), "vol:",lvol(n),"nmain", nmain, "Trp:", Trp, "Tsr:", Tsr, "Trp/Tsr", Trp/Tsr,"and id:", n,'End of domain:',lrem_active
            if (this%cfg%amRoot) print *, "This is the min_thickness", lthc(n), "radius:", minor_radius, "lig percentage:", lper(n),"max length:",Lrim,&
                     & "how many cells",lnum(n), "vol:",Lvol(n),"nmain", nmain, "and id:", n
            Oh=this%fs%visc_l/sqrt(this%fs%rho_l*minor_radius*this%fs%sigma)
            this%size_ratio=1.0_WP/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
            nsat=nmain+1
            diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)

            ! Only the main processor is in charge of creating droplets
            if (this%cfg%amRoot) then
               ! Lrp = twoPi*minor_radius/this%dw
               filename='spray-all/droplets'
               open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
               allocate(points(3,nsat+nmain))
               if (s_info%flag.eq.0) call this%distribute_on_spline(nsat+nmain,s_info,points)
               do l=1,nsat+nmain
                  ! Increment particle counter
                  lp%np_=lp%np_+1
                  ! Make room for new drop
                  call lp%resize(lp%np_)
                  ! Add the drop
                  if (lrem_active) then
                     lp%p(lp%np_)%id  =int(12,8)                                                                               
                  else
                     lp%p(lp%np_)%id  =int(11,8)                                                                               
                  end if
                  if (mod(l,2).eq.1) then
                     lp%p(lp%np_)%d=diam*this%size_ratio                                                                                    
                  else
                     lp%p(lp%np_)%d=diam                                                                                    
                  end if
                  if (s_info%flag.ne.0) then
                     lp%p(lp%np_)%pos = lpos(n,:)
                  else
                     lp%p(lp%np_)%pos = points(:,l)
                  end if
                  lp%p(lp%np_)%vel =lvel(n,:)
                  lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
                  lp%p(lp%np_)%flag=0                                                                                        
                  lp%p(lp%np_)%dt  =0.0_WP                                                                                  
                  lp%p(lp%np_)%Acol=0.0_WP                                                                                  
                  lp%p(lp%np_)%Tcol=0.0_WP
                  write(iunit,*)timet,lp%p(lp%np_)%d,lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3),&
                  &norm2([lp%p(lp%np_)%vel(1),lp%p(lp%np_)%vel(2),lp%p(lp%np_)%vel(3)]),lp%p(lp%np_)%pos(1),&
                  &lp%p(lp%np_)%pos(2),lp%p(lp%np_)%pos(3),lp%p(lp%np_)%id  
               end do
               if (allocated(points)) deallocate(points)
               close(iunit)
               ! Increment monitoring variables
               lp%np_new=lp%np_new+nmain+nsat
               lp%vp_new=lp%vp_new+lvol(n)
               this%np_lig=this%np_lig+nmain+nsat
               this%vof_tf_lig=this%vof_tf_lig+lvol(n)
            end if
            ! empty out the VF
            do m=1,this%ccl_lig%struct(n)%n_
               i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
               this%vf%VF(i,j,k)=0.0_WP
            end do    
         end do
         ! Synchronize VF fields
         call this%vf%cfg%sync(this%vf%VF)
         call this%vf%clean_irl_and_band()
         ! Synchronize particles
         call lp%sync()

         ! Integrate monitoring variables 
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_lig,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_lig    ,1,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
      end if
      contains 
         ! Calculate thickness and struct_type based on moment of inertia
         subroutine get_liginfo()
            implicit none 
            real(WP) :: tmpvol,tmparea
            real(WP), dimension(1:3) :: tmpxvol, tmpL
            integer :: nneigh_moi, nneigh_thickness
            nneigh_moi=2; nneigh_thickness=3
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     if (this%vf%VF(i,j,k) .le. VFlo) cycle
                     ! calculate thickness
                     tmpvol=0.0_WP; tmparea=0.0_WP
                     do kk = k-nneigh_thickness,k+nneigh_thickness
                        do jj = j-nneigh_thickness,j+nneigh_thickness
                           do ii = i-nneigh_thickness,i+nneigh_thickness
                              tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                              if (this%vf%recon_type(ii,jj,kk).eq.2.or.this%vf%recon_type(ii,jj,kk).eq.3.or.this%vf%recon_type(ii,jj,kk).eq.0) then
                                 tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)
                              else
                                 tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)*(2.0/sqrt(pi))
                              end if
                           end do
                        end do
                     end do
                     ! Calculate thickness
                     if (this%vf%VF(i,j,k).le.VFlo) then
                        thickness(i,j,k) = 0.0_WP
                     else if (tmparea .gt. 0.0_WP) then    
                        thickness(i,j,k) = 2.0_WP*tmpvol/(tmparea+tiny(1.0_WP))
                     else
                        thickness(i,j,k) = 3.5_WP*this%cfg%min_meshsize
                     end if

                     if (thickness(i,j,k).lt.this%lmake*this%cfg%min_meshsize)then
                        ! Calculate moi
                        tmpvol=0.0_WP; tmpxvol=0.0_WP; A=0.0_WP
                        ! First pass to accumulate volume, surface area, and position
                        do kk = k-nneigh_moi,k+nneigh_moi
                           do jj = j-nneigh_moi,j+nneigh_moi
                              do ii = i-nneigh_moi,i+nneigh_moi
                                 tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                                 tmpxvol = tmpxvol + this%vf%Lbary(:,ii,jj,kk)*this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                              end do
                           end do
                        end do
                        ! Second pass to accumulate moment of inertia
                        if (tmpvol.gt.1.0e-14_WP) tmpxvol = tmpxvol/tmpvol
                        do kk = k-nneigh_moi,k+nneigh_moi
                           do jj = j-nneigh_moi,j+nneigh_moi
                              do ii = i-nneigh_moi,i+nneigh_moi
                                 ! Location of film node
                                 tmpL = this%vf%Lbary(:,ii,jj,kk) - tmpxvol
                                 A(1,1)=A(1,1)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(2)**2+tmpL(3)**2)
                                 A(2,2)=A(2,2)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(3)**2)
                                 A(3,3)=A(3,3)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(2)**2)
                                 A(1,2)=A(1,2)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(2)
                                 A(1,3)=A(1,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(3)
                                 A(2,3)=A(2,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(2)*tmpL(3)   
                              end do
                           end do
                        end do
                        call eigensolve_moi(A,d)
                        d=max(0.0_WP,d)
                        if ((d(3).lt.1.5_WP*d(2)).and.(d(2).gt.1.5_WP*d(1)) ) then
                           struct_type(i,j,k) = 1
                        else if ( (d(3) .gt. 1.5_WP * d(2)) .and. (d(2) .lt. 1.5_WP * d(1)) ) then
                           struct_type(i,j,k) = 2
                        end if
                     end if
                  end do 
               end do 
            end do
            call this%vf%cfg%sync(thickness)
            call this%vf%cfg%sync(struct_type)
         end subroutine get_liginfo

         !> Function that identifies cells that need a label
         logical function make_label(i,j,k)
            implicit none
            integer, intent(in) :: i,j,k
            if ((this%vf%VF(i,j,k).gt.VFlo).and.(thickness(i,j,k).lt.this%lmake*this%cfg%min_meshsize))then
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

   end subroutine transfer_ligs

   subroutine eigensolve_moi(A,d)
      use mathtools, only: pi
      implicit none
      real(WP), dimension(3,3), intent(inout) :: A
      real(WP), dimension(3), intent(inout) :: d
      real(WP), dimension(3,3) :: e, B
      real(WP) :: x1,x2,phi
      real(WP), dimension(3) :: v1, v2, v3, v
      real(WP) :: norm1, norm2, norm3, maxnorm, invnorm
      integer :: i
      ! Eigensolve moment of inertia
      A(2,1) = A(1,2)
      A(3,1) = A(1,3)
      A(3,2) = A(2,3)
      x1 = A(1,1)**2+A(2,2)**2+A(3,3)**2-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+3*(A(1,2)**2+A(1,3)**2+A(2,3)**2)
      x2 = -(2*A(1,1)-A(2,2)-A(3,3))*(2*A(2,2)-A(1,1)-A(3,3))*(2*A(3,3)-A(1,1)-A(2,2))+9.0_WP*((2*A(3,3)-A(1,1)-A(2,2))*A(1,2)**2+(2*A(2,2)-A(1,1)-A(3,3))*A(1,3)**2+(2*A(1,1)-A(2,2)-A(3,3))*A(2,3)**2)-54.0_WP*A(1,2)*A(1,3)*A(2,3)
      phi = atan2(sqrt(max(0.0_WP, 4*x1**3 - x2**2)), x2)
      d(1) = (A(1,1)+A(2,2)+A(3,3)-2*sqrt(x1)*cos(phi/3.0_WP))/3.0_WP
      d(2) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(x1)*cos((phi+pi)/3.0_WP))/3.0_WP
      d(3) = (A(1,1)+A(2,2)+A(3,3)+2*sqrt(x1)*cos((phi-pi)/3.0_WP))/3.0_WP

      ! Calculate eigenvectors by taking cross products of columns of (A - lambda I)
      do i = 1, 3
         B(1,1)=A(1,1)-d(i); B(2,2)=A(2,2)-d(i); B(3,3)=A(3,3)-d(i)
         B(1,2)=A(1,2); B(1,3)=A(1,3); B(2,3)=A(2,3)
         
         v1(1) = B(2,2)*B(3,3) - B(2,3)**2
         v1(2) = B(2,3)*B(1,3) - B(1,2)*B(3,3)
         v1(3) = B(1,2)*B(2,3) - B(2,2)*B(1,3)
         
         v2(1) = v1(2)
         v2(2) = B(1,1)*B(3,3) - B(1,3)**2
         v2(3) = B(1,2)*B(1,3) - B(1,1)*B(2,3)
         
         v3(1) = v1(3)
         v3(2) = v2(3)
         v3(3) = B(1,1)*B(2,2) - B(1,2)**2

         norm1 = v1(1)**2 + v1(2)**2 + v1(3)**2
         norm2 = v2(1)**2 + v2(2)**2 + v2(3)**2
         norm3 = v3(1)**2 + v3(2)**2 + v3(3)**2
         
         maxnorm = max(norm1, norm2, norm3)
         
         if (maxnorm .gt. 1.0e-30_WP) then
            if (maxnorm == norm1) then
               v = v1
            else if (maxnorm == norm2) then
               v = v2
            else
               v = v3
            end if
            invnorm = 1.0_WP / sqrt(maxnorm)
            e(1,i) = v(1)*invnorm
            e(2,i) = v(2)*invnorm
            e(3,i) = v(3)*invnorm
         else
            ! Degenerate Case: Matrix is identity relative to numerical precision
            e(:,i) = 0.0_WP
            e(i,i) = 1.0_WP
         end if
      end do
      
      ! Guarantee orthogonality even for degenerate eigenvalues via Gram-Schmidt
      e(:,2) = e(:,2) - sum(e(:,1)*e(:,2))*e(:,1)
      norm2 = sum(e(:,2)**2)
      if (norm2 .gt. 1.0e-12_WP) then
         e(:,2) = e(:,2) / sqrt(norm2)
      else
         if (abs(e(1,1)) < 0.9_WP) then
            e(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         else
            e(:,2) = [0.0_WP, 1.0_WP, 0.0_WP]
         end if
         e(:,2) = e(:,2) - sum(e(:,1)*e(:,2))*e(:,1)
         e(:,2) = e(:,2) / sqrt(sum(e(:,2)**2))
      end if
      
      ! Third vector via cross product for consistent right-hand orientation
      e(1,3) = e(2,1)*e(3,2) - e(3,1)*e(2,2)
      e(2,3) = e(3,1)*e(1,2) - e(1,1)*e(3,2)
      e(3,3) = e(1,1)*e(2,2) - e(2,1)*e(1,2)

      A = e
   end subroutine eigensolve_moi

   recursive subroutine quicksort_real(arr, first, last, ids)
      implicit none
      real(WP), dimension(:), intent(inout) :: arr
      integer, intent(in) :: first, last
      integer, dimension(:), intent(inout), optional :: ids
      real(WP) :: pivot, temp_arr
      integer :: temp_id, i, j

      if (first >= last) return
      pivot = arr(first + (last - first) / 2)
      i = first
      j = last
      do while (i <= j)
         do while (arr(i) < pivot)
            i = i + 1
         end do
         do while (arr(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp_arr = arr(i); arr(i) = arr(j); arr(j) = temp_arr
            if (present(ids)) then
               temp_id = ids(i); ids(i) = ids(j); ids(j) = temp_id
            end if
            i = i + 1
            j = j - 1
         end if
      end do
      if (first < j) then
         if (present(ids)) then
            call quicksort_real(arr, first, j, ids)
         else
            call quicksort_real(arr, first, j)
         end if
      end if
      if (i < last) then
         if (present(ids)) then
            call quicksort_real(arr, i, last, ids)
         else
            call quicksort_real(arr, i, last)
         end if
      end if
   end subroutine quicksort_real


   subroutine fit_spline(this,n,s_info,ccl)
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use fitpack_core, only: curfit, splev
      use vfs_class, only: VFlo,VFhi
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      integer, intent(in) :: n
      type(spline_info), intent(out) :: s_info
      type(cclabel), intent(inout) :: ccl
      integer :: m, local_point_count, num_procs, total_points, i, j, k, ier, nest_max, lwrk, unique_count,ierr
      real(WP), dimension(3) :: end_p
      real(WP), dimension(:,:), allocatable :: local_points
      real(WP), dimension(:,:), allocatable :: points, unique_points, sorted_points
      integer, dimension(:), allocatable :: recv_counts
      integer, dimension(:), allocatable :: displacements
      real(WP), dimension(:), allocatable :: t_param
      real(WP) :: s, fp_x, fp_y, fp_z
      real(WP), dimension(:), allocatable :: weights
      logical, dimension(:), allocatable :: is_used
      real(WP) :: min_dist_sq, dist_sq
      real(WP), dimension(3) :: last_sorted_point
      integer :: start_point_idx, best_idx

      real(WP), dimension(:), allocatable :: wrk
      integer, dimension(:), allocatable :: iwrk
      real(WP) :: tolerance

      integer :: num_eval
      real(WP), dimension(:), allocatable :: t_eval, x_eval, y_eval, z_eval
      real(WP) :: step, t_current
      integer :: e_flag
      integer, allocatable :: cluster_id(:)
      integer, allocatable :: cluster_counts(:)
      character(len=20) :: my_string

      s_info%length = 0.0_WP
      s_info%flag = 0

      allocate(local_points(3, ccl%struct(n)%n_)); local_points = 0.0_WP
      write (my_string, '(i0)') n
      do m=1,ccl%struct(n)%n_
         i=ccl%struct(n)%map(1,m)
         j=ccl%struct(n)%map(2,m)
         k=ccl%struct(n)%map(3,m)

         if (this%vf%recon_type(i,j,k).eq.1) then
            local_points(:,m) = getDatum(this%vf%liquid_gas_interface(i,j,k))
            local_points(1,m) = local_points(1,m)-ccl%struct(n)%per(1)*this%cfg%xL
            local_points(2,m) = local_points(2,m)-ccl%struct(n)%per(2)*this%cfg%yL
            local_points(3,m) = local_points(3,m)-ccl%struct(n)%per(3)*this%cfg%zL
         else
            local_points(:,m) = this%vf%LBARY(:,i,j,k)
            local_points(1,m) = local_points(1,m)-ccl%struct(n)%per(1)*this%cfg%xL
            local_points(2,m) = local_points(2,m)-ccl%struct(n)%per(2)*this%cfg%yL
            local_points(3,m) = local_points(3,m)-ccl%struct(n)%per(3)*this%cfg%zL
         end if
      end do

      num_procs = this%cfg%nproc

      allocate(recv_counts(num_procs))
      recv_counts = 0

      call MPI_GATHER(ccl%struct(n)%n_, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, this%cfg%comm, ierr)

      if (this%cfg%amRoot) then
         total_points = sum(recv_counts)
         allocate(points(3, total_points))
         allocate(is_used(total_points))
         allocate(displacements(num_procs))
         if (total_points.gt.0) then
            displacements(1) = 0
            do i=2,num_procs
               displacements(i) = displacements(i-1) + recv_counts(i-1)
            end do
            is_used = .false.
         end if
         call MPI_GATHERV(local_points, ccl%struct(n)%n_*3, MPI_REAL_WP, points, recv_counts*3, displacements*3, MPI_REAL_WP, 0, this%cfg%comm, ierr)
      else
         call MPI_GATHERV(local_points, ccl%struct(n)%n_*3, MPI_REAL_WP, local_points, recv_counts, recv_counts, MPI_REAL_WP, 0, this%cfg%comm, ierr)
      end if

      if (this%cfg%amRoot) then 
         tolerance = (4*this%cfg%min_meshsize)**2
         if (total_points .gt. 4) then
            s_info%flag = 0
            end_p = points(:,1)
            start_point_idx = 1

            allocate(cluster_id(total_points))
            cluster_id = 0
            unique_count = 0

            do m = 1, total_points
               if (cluster_id(m) .eq. 0) then
                  unique_count = unique_count + 1
                  cluster_id(m) = unique_count
                  do i = m + 1, total_points
                     if (cluster_id(i) .eq. 0) then
                        if (sum((points(:, m) - points(:, i))**2) .le. tolerance) then
                           cluster_id(i) = unique_count
                        end if
                     end if
                  end do
               end if
            end do
            if (unique_count .gt. 0) then
               allocate(unique_points(3, unique_count), cluster_counts(unique_count))
               unique_points = 0.0_WP
               cluster_counts = 0
               do m = 1, total_points
                  unique_points(:, cluster_id(m)) = unique_points(:, cluster_id(m)) + points(:, m)
                  cluster_counts(cluster_id(m)) = cluster_counts(cluster_id(m)) + 1
               end do
               do m = 1, unique_count
                  if (cluster_counts(m) .gt. 0) then
                     unique_points(:, m) = unique_points(:, m) / real(cluster_counts(m))
                  end if
               end do
            end if
            deallocate(points)
            allocate(points(3, unique_count))
            allocate(sorted_points(3, total_points))
            points = unique_points
            total_points = unique_count
            deallocate(cluster_id, unique_points, cluster_counts)

            if (total_points .gt. 4) then
               call order_points(total_points, points, unique_count, sorted_points)
               total_points = unique_count
               deallocate(points)
               allocate(points(3, total_points))
               points = sorted_points
               deallocate(sorted_points)
               allocate(t_param(total_points))
               allocate(weights(total_points))
               t_param(1) = 0.0_WP
               do m = 2, total_points
                  t_param(m) = t_param(m-1)+sqrt(sqrt((points(1,m)-points(1,m-1))**2+(points(2,m)-points(2,m-1))**2+(points(3,m)-points(3,m-1))**2))
               end do
               if (t_param(total_points) .gt. VFlo) then
                  t_param = t_param/t_param(total_points)
               else
                  s_info%flag = 1
               end if
               s = real(total_points, WP) * (this%cfg%min_meshsize)**2
               weights = 1.0_WP
               k = 3
               nest_max = max(total_points+k+1, 2*k+3)
               lwrk = total_points * (k + 1) + nest_max * (7 + 3 * k)
   
               if (allocated(s_info%t_knots_x)) deallocate(s_info%t_knots_x)
               if (allocated(s_info%c_coeffs_x)) deallocate(s_info%c_coeffs_x)
               if (allocated(s_info%t_knots_y)) deallocate(s_info%t_knots_y)
               if (allocated(s_info%c_coeffs_y)) deallocate(s_info%c_coeffs_y)
               if (allocated(s_info%t_knots_z)) deallocate(s_info%t_knots_z)
               if (allocated(s_info%c_coeffs_z)) deallocate(s_info%c_coeffs_z)
   
               allocate(wrk(lwrk))
               allocate(iwrk(nest_max))
               allocate(s_info%t_knots_x(nest_max))
               allocate(s_info%t_knots_y(nest_max))
               allocate(s_info%t_knots_z(nest_max))
               allocate(s_info%c_coeffs_x(nest_max))
               allocate(s_info%c_coeffs_y(nest_max))
               allocate(s_info%c_coeffs_z(nest_max))
   
               call curfit(iopt=0, m=total_points, x=t_param, y=points(1,:), w=weights, &
               xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
               n=s_info%n_knots_x, t=s_info%t_knots_x, c=s_info%c_coeffs_x, fp=fp_x, &
               wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
               if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for X: ", ier
   
               call curfit(iopt=0, m=total_points, x=t_param, y=points(2,:), w=weights, &
                     xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
                     n=s_info%n_knots_y, t=s_info%t_knots_y, c=s_info%c_coeffs_y, fp=fp_y, &
                     wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
               if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for Y: ", ier
   
               call curfit(iopt=0, m=total_points, x=t_param, y=points(3,:), w=weights, &
                     xb=t_param(1), xe=t_param(total_points), k=k, s=s, nest=nest_max, &
                     n=s_info%n_knots_z, t=s_info%t_knots_z, c=s_info%c_coeffs_z, fp=fp_z, &
                     wrk=wrk, lwrk=lwrk, iwrk=iwrk, ier=ier)
               if (ier .gt. 0) s_info%flag = 1!print *, "Error in CURFIT for Z: ", ier
   
               if (s_info%flag.eq.0) then
                  num_eval = total_points*10
                  e_flag = 0
                  allocate(t_eval(num_eval))
                  allocate(x_eval(num_eval))
                  allocate(y_eval(num_eval))
                  allocate(z_eval(num_eval))
                  
                  step = t_param(total_points) / real(num_eval - 1, WP)
                  do m = 1, num_eval
                     t_eval(m) = real(m-1, WP) * step
                  end do
                  
                  call splev(s_info%t_knots_x, s_info%n_knots_x, s_info%c_coeffs_x, k, t_eval, x_eval, num_eval, e_flag, ier)
                  ! if (ier .ne. 0) print *, "Error in SPLEV for X: ", ier
                  
                  call splev(s_info%t_knots_y, s_info%n_knots_y, s_info%c_coeffs_y, k, t_eval, y_eval, num_eval, e_flag, ier)
                  ! if (ier .ne. 0) print *, "Error in SPLEV for Y: ", ier
                  
                  call splev(s_info%t_knots_z, s_info%n_knots_z, s_info%c_coeffs_z, k, t_eval, z_eval, num_eval, e_flag, ier)
                  ! if (ier .ne. 0) print *, "Error in SPLEV for Z: ", ier
   
                  s_info%length = 0.0_WP
                  do i = 2, num_eval
                     s_info%length = s_info%length + sqrt((x_eval(i)-x_eval(i-1))**2+(y_eval(i)-y_eval(i-1))**2+(z_eval(i)-z_eval(i-1))**2)
                  end do
               end if
            else
               s_info%flag = 1
            end if
         else
            s_info%flag = 1
         end if
      end if
      call MPI_BCAST(s_info%flag, 1, MPI_INTEGER, 0, this%cfg%comm, ierr)
      call MPI_BCAST(s_info%length, 1, MPI_REAL_WP, 0, this%cfg%comm, ierr)
      if (allocated(local_points)) deallocate(local_points)
      if (allocated(recv_counts)) deallocate(recv_counts)
      if (allocated(points)) deallocate(points)
      if (allocated(is_used)) deallocate(is_used)
      if (allocated(displacements)) deallocate(displacements)
      if (allocated(t_param)) deallocate(t_param)
      if (allocated(weights)) deallocate(weights)
      if (allocated(wrk)) deallocate(wrk)
      if (allocated(iwrk)) deallocate(iwrk)
      if (allocated(t_eval)) deallocate(t_eval)
      if (allocated(x_eval)) deallocate(x_eval)
      if (allocated(y_eval)) deallocate(y_eval)
      if (allocated(z_eval)) deallocate(z_eval)
   end subroutine fit_spline

   subroutine distribute_on_spline(this,n_part,s_info,points)
      use fitpack_core, only: splev
      implicit none
      class(breakup), intent(inout) :: this
      integer, intent(in) :: n_part
      type(spline_info), intent(in) :: s_info
      real(WP), dimension(:,:), intent(out) :: points
      integer :: k, ier, n_map, idx, i, m
      real(WP), dimension(:), allocatable :: t_map, x_map, y_map, z_map, dist_map
      real(WP) :: step, total_len, target_dist, step_map, frac
      integer :: e_flag

      if (this%cfg%amRoot) then 
         e_flag = 0
         k = 3
         n_map = 2000

         allocate(t_map(n_map), x_map(n_map), y_map(n_map), z_map(n_map), dist_map(n_map))

         step_map = 1.0_WP / real(n_map - 1, WP)
         do i = 1, n_map
            t_map(i) = real(i-1, WP) * step_map
         end do

         call splev(s_info%t_knots_x, s_info%n_knots_x, s_info%c_coeffs_x, k, t_map, x_map, n_map, e_flag, ier)
         ! if (ier .ne. 0) print *, "Error in SPLEV for X: ", ier
         call splev(s_info%t_knots_y, s_info%n_knots_y, s_info%c_coeffs_y, k, t_map, y_map, n_map, e_flag, ier)
         ! if (ier .ne. 0) print *, "Error in SPLEV for Y: ", ier
         call splev(s_info%t_knots_z, s_info%n_knots_z, s_info%c_coeffs_z, k, t_map, z_map, n_map, e_flag, ier)
         ! if (ier .ne. 0) print *, "Error in SPLEV for Z: ", ier
         dist_map(1) = 0.0_WP
         do i = 2, n_map
            dist_map(i) = dist_map(i-1) + sqrt((x_map(i)-x_map(i-1))**2 + (y_map(i)-y_map(i-1))**2 + (z_map(i)-z_map(i-1))**2)
         end do
         total_len = dist_map(n_map)
         step = total_len / real(n_part + 1, WP)
         
         do m = 1, n_part
            target_dist = real(m, WP) * step
            idx = 1
            do i = 1, n_map-1
               if (dist_map(i+1) .ge. target_dist) then
                  idx = i
                  exit
               end if
            end do
            
            if (abs(dist_map(idx+1) - dist_map(idx)) .gt. 1.0e-12_WP) then
               frac = (target_dist - dist_map(idx)) / (dist_map(idx+1) - dist_map(idx))
            else
               frac = 0.0_WP
            end if
            
            points(1,m) = x_map(idx) + frac * (x_map(idx+1) - x_map(idx))
            points(2,m) = y_map(idx) + frac * (y_map(idx+1) - y_map(idx))
            points(3,m) = z_map(idx) + frac * (z_map(idx+1) - z_map(idx))

            if (this%cfg%xper.and.points(1,m).lt.this%cfg%x(this%cfg%imin)) points(1,m)=points(1,m)+this%cfg%xL
            if (this%cfg%yper.and.points(2,m).lt.this%cfg%y(this%cfg%jmin)) points(2,m)=points(2,m)+this%cfg%yL
            if (this%cfg%zper.and.points(3,m).lt.this%cfg%z(this%cfg%kmin)) points(3,m)=points(3,m)+this%cfg%zL
         end do

         if (allocated(t_map)) deallocate(t_map, x_map, y_map, z_map, dist_map)
      end if
   end subroutine distribute_on_spline

   subroutine order_points(total_points, points, total_ordered_points, ordered_points)
      integer, intent(inout) :: total_points
      real(WP), intent(in) :: points(3, total_points)
      integer, intent(out) :: total_ordered_points
      real(WP), allocatable, intent(out) :: ordered_points(:,:)
      
      logical, allocatable :: in_tree(:)
      real(WP), allocatable :: min_tree_dist(:)
      integer, allocatable :: parent(:)
      integer :: i, j, new_node
      real(WP) :: dist, current_dist
      
      integer, allocatable :: degree(:), node_offset(:), adjacency_list(:), local_offset(:)
      
      integer, allocatable :: queue(:), path(:)
      logical, allocatable :: visited(:)
      integer :: q_head, q_tail, node_a, node_b

      if (total_points .le. 1) then
         total_ordered_points = total_points
         allocate(ordered_points(3, max(1, total_points)))
         if (total_points .eq. 1) ordered_points(:,1) = points(:,1)
         return
      end if

      allocate(parent(total_points), in_tree(total_points), min_tree_dist(total_points))
      in_tree = .false.; parent = 0; min_tree_dist = huge(1.0_WP)
      min_tree_dist(1) = 0.0_WP
      
      do i = 1, total_points
         current_dist = huge(1.0_WP); new_node = -1
         do j = 1, total_points
            if (.not. in_tree(j) .and. min_tree_dist(j) .lt. current_dist) then
               current_dist = min_tree_dist(j); new_node = j
            end if
         end do
         
         if (new_node .eq. -1) exit
         in_tree(new_node) = .true.
         
         do j = 1, total_points
            if (.not. in_tree(j)) then
               dist = sum((points(:, new_node) - points(:, j))**2)
               if (dist .lt. min_tree_dist(j)) then
                  min_tree_dist(j) = dist; parent(j) = new_node
               end if
            end if
         end do
      end do
      deallocate(in_tree, min_tree_dist)

      allocate(degree(total_points), node_offset(total_points+1), local_offset(total_points), adjacency_list(2*total_points - 2))
      degree = 0
      
      do i = 2, total_points
         j = parent(i)
         if (j .gt. 0) then
            degree(i) = degree(i) + 1
            degree(j) = degree(j) + 1
         end if
      end do
      
      node_offset(1) = 1
      do i = 1, total_points
         node_offset(i+1) = node_offset(i) + degree(i)
         local_offset(i) = node_offset(i)
      end do
      
      do i = 2, total_points
         j = parent(i)
         if (j .gt. 0) then
            adjacency_list(local_offset(i)) = j; local_offset(i) = local_offset(i) + 1
            adjacency_list(local_offset(j)) = i; local_offset(j) = local_offset(j) + 1
         end if
      end do
      deallocate(degree, local_offset, parent)

      allocate(queue(total_points), visited(total_points), path(total_points))
      visited = .false.
      
      q_head = 1; q_tail = 1; queue(q_head) = 1; visited(1) = .true.
      node_a = 1
      
      do while (q_head .le. q_tail)
         new_node = queue(q_head); q_head = q_head + 1
         node_a = new_node 
         do i = node_offset(new_node), node_offset(new_node+1) - 1
            j = adjacency_list(i)
            if (.not. visited(j)) then
               visited(j) = .true.
               q_tail = q_tail + 1; queue(q_tail) = j
            end if
         end do
      end do

      visited = .false.; path = 0
      q_head = 1; q_tail = 1; queue(q_head) = node_a; visited(node_a) = .true.
      node_b = node_a
      
      do while (q_head .le. q_tail)
         new_node = queue(q_head); q_head = q_head + 1
         node_b = new_node 
         do i = node_offset(new_node), node_offset(new_node+1) - 1
            j = adjacency_list(i)
            if (.not. visited(j)) then
               visited(j) = .true.; path(j) = new_node
               q_tail = q_tail + 1; queue(q_tail) = j
            end if
         end do
      end do

      total_ordered_points = 0; new_node = node_b
      do while (new_node .ne. 0)
         total_ordered_points = total_ordered_points + 1
         new_node = path(new_node)
      end do

      allocate(ordered_points(3, total_ordered_points))
      new_node = node_b; i = 1
      do while (new_node .ne. 0)
         ordered_points(:, i) = points(:, new_node)
         i = i + 1
         new_node = path(new_node)
      end do

      deallocate(node_offset, adjacency_list, queue, visited, path)
   end subroutine order_points
end module breakup_class
