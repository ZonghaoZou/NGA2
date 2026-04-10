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
   implicit none
   private

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
      real(WP) :: vof_tf_drop,vof_deleted
      integer  :: np_drop

      !> Film transfer parameters
      logical  :: use_film_transfer,output_filmstats
      real(WP) :: frp,fmin,fd0,fbvol2dvol,fthld,minfdrop
      real(WP) :: fnumcell,vof_tf_film,fthc_percent
      integer  :: maxdpcell,np_film,film_tracker

      !> Ligament transfer parameters
      logical  :: use_lig_transfer
      real(WP) :: lmin,lmake,lper,lstratio
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
      !< Conversion from ligament to droplets                           : 12

   contains
      procedure :: init            !< Initialize breakup model
      procedure :: transfer_drops  !< Transfer drops to Lagrangian
      procedure :: transfer_films  !< Transfer films to Lagrangian drops
      procedure :: transfer_ligs   !< Transfer ligaments to Lagrangian drops
      ! procedure :: eigensolve_moi  !< Eigensolve moment of inertia
   end type breakup

   !> Hardcode inlet parameters (shared with atom_class)
   real(WP), parameter :: dl=0.003_WP

contains

   !> Initialize breakup model
   subroutine init(this,input,cfg,vf,fs)
      implicit none
      class(breakup), intent(inout) :: this
      class(inputfile), intent(inout) :: input
      class(ibconfig), target, intent(in) :: cfg
      class(vfs), target, intent(in) :: vf
      class(tpns), target, intent(in) :: fs
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
         this%emax=0.75_WP
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
         this%fthc_percent=0.15_WP
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
         this%lmake=1.5_WP
         this%lper=0.9_WP
         this%lstratio=1.5_WP
         this%vof_tf_lig=0.0_WP
         this%np_lig=0
      end if

   end subroutine init


   !> Transfer droplets to Lagrangian representation
   subroutine transfer_drops(this,lp)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi
      implicit none
      class(breakup), intent(inout) :: this
      class(lpt),         intent(inout) :: lp
      real(WP), dimension(:)    , allocatable :: dvol,drem
      real(WP), dimension(:,:)  , allocatable :: dpos,dvel
      real(WP), dimension(:,:,:), allocatable :: dmoi
      integer :: n,m,ierr,i,j,k
      real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
      logical :: transfer
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      logical :: drem_active
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
            if (lmin.eq.0.0_WP) lmin=lmid
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            transfer=(ecc.le.this%emax)
         end if

         drem_active=.false.
         if (drem(n).gt.0.0_WP) then
            transfer=.true.; drem_active=.true.
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
      call this%vf%sync_interface()
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
   subroutine transfer_films(this,lp)
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
      character(len=str_medium) :: filename
      real(WP), dimension(:), allocatable :: fvol,fthc,frem,fcnt,fcurv,fcsa,fthick_,fthick,fivol_,fivol
      real(WP), dimension(:), allocatable :: sort_ke,fcsar_,fcsar,fcurvr,fcurvr_
      integer, dimension(:), allocatable ::  sort_id,plist,dispels
      integer :: n,nn,m,i,j,k,ii,jj,kk,ncell_,tmp_id,l,totalnewp,np_start,np_old,count,ierr,ind,ip,iunit,rank
      integer :: totalcell,nfthcpercent,dropcounter
      real(WP)  :: tmp_ke,Vt,Vl,Vd,alpha,beta,mySA,tmp_thic
      real(WP), dimension(3) :: nref,tref,sref
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
               Vt=0.0_WP; Vl=0.0_WP; np_old=lp%np_
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
                        lp%p(lp%np_)%vel =this%fs%cfg%get_velocity(pos=lp%p(lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
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
                  lp%p(lp%np_)%vel =this%fs%cfg%get_velocity(pos=lp%p(lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
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
         end do
         ! Gather the number of newly generated particles from each processor due to film burst
         totalnewp = 0
         ! Get number of particle generated for each processor
         call MPI_AllGATHER(this%np_film,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%vf%cfg%comm,ierr)
         totalnewp= sum(plist)
         ! If there is any particle generated
         if (totalnewp .gt. 0) then
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


   subroutine transfer_ligs(this,lp)
      use vfs_class, only: VFlo,VFhi
      use mathtools, only: pi,twoPi
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use messager, only: die
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      class(lpt), intent(inout) :: lp
      real(WP), dimension(:)    , allocatable :: lvol,lthc,llen,lnum,lper
      real(WP), dimension(:)    , allocatable :: lrem,lSR,xmin,xmax,ymin,ymax,zmin,zmax
      real(WP), dimension(:,:)  , allocatable :: lpos,lvel
      real(WP), dimension(:,:,:), allocatable :: lmoi
      integer :: n,m,ierr,i,j,k,l,ii,jj,kk,iunit,totalnewp,count,ip,rank
      real(WP) :: x,y,z,x0,y0,z0,lmax,lmid,lmin
      real(WP) :: Vt,Vl,Vd,minor_radius,diam,Vrim,Lrim
      real(WP) :: Trp,Lrp,Tsr,SR_tmp,Oh,b
      real(WP), dimension(1:3) :: tangent
      real(WP), dimension(:,:,:,:), allocatable :: SR
      integer  :: nmain,nsat
      real(WP), dimension(:,:,:), allocatable :: thickness
      integer,  dimension(:,:,:), allocatable :: struct_type
      real(WP), dimension(:,:), allocatable :: points
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      logical :: lrem_active

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
               lthc(n)=min(lthc(n),thickness(i,j,k))
               if (struct_type(i,j,k).eq.1) lper(n)=lper(n)+1.0_WP
               ! Check if drop touches auto-transfer layer
               if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
               &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
               &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
               &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
               &   k.ge.this%vf%cfg%kmax-this%nlayer) lrem(n)=1.0_WP
               ! Get the structures's locally largest and smallest x,y,z locations
               do l=1,2
                  if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                     d = calculateCentroid(this%vf%interface_polygon(l,i,j,k))
                     xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
                     ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
                     zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
                  end if
               end do
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
            ! Assume a cylinder ligament
            Lrim=llen(n)
            Vrim=lvol(n)
            minor_radius=sqrt(Vrim/pi/Lrim)
            ! Drop size method from Kim & Moin (2020)
            nmain=floor(this%dw*Lrim/(twoPi*minor_radius))
            ! Calculate breakup time scale based on inviscid RP instability analysis
            Trp=2.91258_WP*sqrt(this%fs%rho_l*minor_radius**3/this%fs%sigma)
            ! Calcuate time scale based on maximum local strainrate
            Tsr=1.0_WP/lSR(n)
            ! Only breakup if minimum thickness is reached, sufficient volume of the ligament, enough of local ligament-like structures,
            ! local time scale asscoiated with strain rate is on par or bigger than the RP time scale, and its length is longer than the inviscid most unstable wavelength
            lrem_active = .false.
            if ((lthc(n).le.this%lmin*this%cfg%min_meshsize).and.(lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.this%lper).and.(Trp.le.Tsr).and.(nmain.ge.1)) then
            else if(lrem(n).gt.0.0_WP.and.(lvol(n).ge.this%cfg%min_meshsize**3)) then
               lrem_active=.true.
            else
               cycle
            end if
            
            if (this%vf%cfg%amRoot) print *, "This is the min_thickness", lthc(n), ",lig percentage:", lper(n),"max length:",llen(n),&
            & "how many cells",lnum(n), "vol:",lvol(n),"nmain", nmain, "Trp:", Trp, "Tsr:", Tsr, "Trp/Tsr", Trp/Tsr,"and id:", n,'End of domain:',lrem_active
            
            Oh=this%fs%visc_l/sqrt(this%fs%rho_l*minor_radius*this%fs%sigma)
            this%size_ratio=1.0_WP/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
            nsat=nmain+1
            diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)

            ! Only the main processor is in charge of creating droplets
            if (this%cfg%amRoot) then
               Lrp = twoPi*minor_radius/this%dw
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
                  if (llen(n).eq.0.0_WP) then
                     lp%p(lp%np_)%pos = lpos(n,:)
                  else
                     lp%p(lp%np_)%pos =lpos(n,:)+0.5_WP*Lrp*(l-(nmain+1))*lmoi(n,:,1)
                  end if
                  lp%p(lp%np_)%vel =lvel(n,:)
                  lp%p(lp%np_)%ind =lp%cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
                  lp%p(lp%np_)%flag=0                                                                                        
                  lp%p(lp%np_)%dt  =0.0_WP                                                                                  
                  lp%p(lp%np_)%Acol=0.0_WP                                                                                  
                  lp%p(lp%np_)%Tcol=0.0_WP
               end do
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
                     ! calculate thickness
                     tmpvol=0.0_WP; tmparea=0.0_WP
                     do kk = k-nneigh_thickness,k+nneigh_thickness
                        do jj = j-nneigh_thickness,j+nneigh_thickness
                           do ii = i-nneigh_thickness,i+nneigh_thickness
                              tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
                              tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)
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

                     if ((this%vf%VF(i,j,k).gt.VFlo).and.(thickness(i,j,k).lt.this%lmake*this%cfg%min_meshsize))then
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
                        tmpxvol = tmpxvol/tmpvol
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
                        if (d(3).gt.(this%lstratio*d(1))) struct_type(i,j,k)=struct_type(i,j,k)+1
                        if (d(3).gt.(this%lstratio*d(2))) struct_type(i,j,k)=struct_type(i,j,k)+1
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

end module breakup_class
