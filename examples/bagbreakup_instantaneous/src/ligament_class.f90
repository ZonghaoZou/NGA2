!> Definition for a ligament atomization class
module ligament_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use lpt_class,         only: lpt
   use cclabel_class,     only: cclabel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   use pardata_class,     only: pardata
   implicit none
   private
   
   public :: ligament
   
   !> Ligament object
   type :: ligament
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      !type(ddadi)       :: vs    !< DDADI solver for velocity
      type(timetracker) :: time  !< Time info
      type(cclabel)     :: ccl,ccl_film,ccl_lig   !< CCLabel to transfer droplets

      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      real(WP), dimension(:,:,:), allocatable :: thickness,struct_type
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer  :: nlayer=4                 !< Size of buffer layer for VOF removal
      
      !> Timing info
      type(monitor) :: timefile !< Timing monitoring
      type(timer)   :: tstep    !< Timer for step
      type(timer)   :: tvel     !< Timer for velocity
      type(timer)   :: tpres    !< Timer for pressure
      type(timer)   :: tvof     !< Timer for VOF
      type(timer)   :: ttrans   !< Timer for droplet transfer
      type(timer)   :: tftrans   !< Timer for film transfer
      type(timer)   :: tltrans   !< Timer for ligament transfer
      
      !> Drop transfer modeling
      logical :: use_drop_transfer !< Do we use droplet transfer
      logical :: use_film_transfer !< Do we use film transfer
      logical :: use_lig_transfer  !< Do we use ligament transfer
      type(lpt)      :: lp         !< Lagrangian particle tracking
      type(monitor)  :: pfile      !< Particle monitoring
      type(partmesh) :: pmesh      !< Particle mesh for lpt
      real(WP) :: dmax             !< Maximum diameter for transfer
      real(WP) :: dmin             !< Minimum diameter below which transfer is automatic
      real(WP) :: ddel             !< Minimum diameter below which structure is directly deleted
      real(WP) :: emax             !< Maximum eccentricity for transfer
      real(WP) :: vof_tf_drop      !< Integral of VOF transfered by conversion to droplet
      real(WP) :: vof_deleted      !< Integral of VOF deleted
      integer  :: np_drop

      real(WP) :: frp
      real(WP) :: fmin
      real(WP) :: fd0
      real(WP) :: fnumcell
      real(WP) :: vof_tf_film       
      integer  :: np_film

      real(WP) :: lmin
      real(WP) :: lmake
      real(WP) :: lper
      real(WP) :: lstratio
      real(WP) :: dw
      real(WP) :: ldmin
      real(WP) :: size_ratio
      real(WP) :: vof_tf_lig
      integer  :: np_lig
      !> Provide a pardata and an event tracker for saving restarts
      type(event)   :: save_evt
      type(pardata) :: df
      logical :: restarted
      
   contains
      procedure :: init             !< Initialize nozzle simulation
      procedure :: step             !< Advance nozzle simulation by one time step
      procedure :: final            !< Finalize nozzle simulation
      procedure :: transfer_drops   !< Transfer drops to a Lagrangian representation
      procedure :: transfer_films   !< Transfer films to Lagrangian drops based on Jackiw and Ashgriz's model
      procedure :: transfer_ligs    !< Transfer ligaments to Lagrangian drops based on Kim and Moin's model
   end type ligament
   
   
contains
   
    !> Transfer droplet to Lagrangian representation
    subroutine transfer_drops(this)
        use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
        use parallel,  only: MPI_REAL_WP
        use mathtools, only: pi
        use messager,  only: die
        class(ligament), intent(inout) :: this
        real(WP), dimension(:)    , allocatable :: dvol
        real(WP), dimension(:,:)  , allocatable :: dpos
        real(WP), dimension(:,:)  , allocatable :: dvel
        real(WP), dimension(:,:,:), allocatable :: dmoi
        real(WP), dimension(:)    , allocatable :: drem
        integer :: n,m,ierr,i,j,k,iunit,np_start
        real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
        character(len=str_medium) :: filename
        logical :: transfer
        ! Moment of inertia calculation using lapack
        real(WP), dimension(:), allocatable, save :: work !< Saved!
        integer, save :: lwork                            !< Saved!
        real(WP), dimension(1) :: lwork_query
        real(WP), dimension(3) :: d
        real(WP), dimension(3,3) :: A
        integer :: info
        
        ! Query optimal work array size
        if (.not.allocated(work)) then
        call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
        lwork=int(lwork_query(1)); allocate(work(lwork))
        end if
        
        ! Start by performing a CCL
        call this%ccl%build(make_label,same_label)
        
        ! Allocate droplet stats arrays
        allocate(dvol(1:this%ccl%nstruct        )); dvol=0.0_WP
        allocate(dpos(1:this%ccl%nstruct,1:3    )); dpos=0.0_WP
        allocate(dvel(1:this%ccl%nstruct,1:3    )); dvel=0.0_WP
        allocate(dmoi(1:this%ccl%nstruct,1:3,1:3)); dmoi=0.0_WP
        allocate(drem(1:this%ccl%nstruct        )); drem=0.0_WP
        
        ! First pass to accumulate volume, position, and velocity
        do n=1,this%ccl%nstruct
        ! Loop over cells in structure
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
            dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
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
        
        ! Second pass to accumulate moment of inertia
        do n=1,this%ccl%nstruct
        ! Get drop barycenter
        x0=dpos(n,1)/dvol(n)
        y0=dpos(n,2)/dvol(n)
        z0=dpos(n,3)/dvol(n)
        ! Loop over cells in structure
        do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m)
            j=this%ccl%struct(n)%map(2,m)
            k=this%ccl%struct(n)%map(3,m)
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
        
        ! Third pass to generate normalized drop stats
        do n=1,this%ccl%nstruct
        ! Get drop barycenter, accounting for periodicity
        dpos(n,:)=dpos(n,:)/dvol(n)
        if (this%vf%cfg%xper.and.dpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) dpos(n,1)=dpos(n,1)+this%vf%cfg%xL
        if (this%vf%cfg%yper.and.dpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%vf%cfg%yL
        if (this%vf%cfg%zper.and.dpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%vf%cfg%zL
        ! Get drop velocity
        dvel(n,:)=dvel(n,:)/dvol(n)
        end do
        
        ! Zero out monitoring variables
        this%vof_tf_drop=0.0_WP
        this%vof_deleted=0.0_WP
        this%np_drop=0

        ! Transfer drops based on our criteria
        do n=1,this%ccl%nstruct
        
        ! Compute diameter
        diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)
        
        ! Decide whether to transfer based on diameter
        if (diam.gt.this%dmax) then
            ! Too big to transfer
            transfer=.false.
        else if (diam.le.this%ddel) then
            ! Too small to track, delete immediately
            transfer=.false.
            ! Zero out VF in the structure
            do m=1,this%ccl%struct(n)%n_
                this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            ! Increment monitoring variables
            this%vof_deleted=this%vof_deleted+dvol(n)
        else if (diam.gt.this%ddel.and.diam.le.this%dmin) then
            ! Small enough to transfer automatically
            transfer=.true.
        else
            ! In between, check eccentricity from moment of inertia tensor
            A=dmoi(n,:,:)
            call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
            d=max(0.0_WP,d)                             !< Get rid of very small negative values (due to machine accuracy)
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            if (ecc.gt.this%emax) then
                ! Too eccentric to transfer yet
                transfer=.false.
            else
                ! Spherical enough to transfer
                transfer=.true.
            end if
        end if
        
        ! Force transfer if drop touches auto-transfer layer
        if (drem(n).gt.0.0_WP) transfer=.true.
        
        ! Perform transfer
        if (transfer) then
            
            ! Root creates a new Lagrangian drop
            if (this%vf%cfg%amRoot) then
                np_start=this%lp%np_
                ! Increment particle counter
                this%lp%np_=this%lp%np_+1
                ! Make room for new drop
                call this%lp%resize(this%lp%np_)
                ! Add the drop
                this%lp%p(this%lp%np_)%id  =int(1,8)
                this%lp%p(this%lp%np_)%d   =diam
                this%lp%p(this%lp%np_)%pos =dpos(n,:)
                this%lp%p(this%lp%np_)%vel =dvel(n,:)
                this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(dpos(n,:),[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])
                this%lp%p(this%lp%np_)%flag=0
                this%lp%p(this%lp%np_)%dt  =0.0_WP
                this%lp%p(this%lp%np_)%Acol=0.0_WP
                this%lp%p(this%lp%np_)%Tcol=0.0_WP

                !!! Write to droplet list !!!
                ! Open the file
                filename='spray-all/droplets'
                open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
                if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
                ! Output diameter, velocity, and position
                write(iunit,*) this%lp%p(this%lp%np_)%d,this%lp%p(this%lp%np_)%vel(1),this%lp%p(this%lp%np_)%vel(2),this%lp%p(this%lp%np_)%vel(3),&
                &norm2([this%lp%p(this%lp%np_)%vel(1),this%lp%p(this%lp%np_)%vel(2),this%lp%p(this%lp%np_)%vel(3)]),this%lp%p(this%lp%np_)%pos(1),&
                &this%lp%p(this%lp%np_)%pos(2),this%lp%p(this%lp%np_)%pos(3),this%lp%p(this%lp%np_)%id  
                ! Close the file
                close(iunit)

            end if
            
            ! Zero out VF in the structure
            do m=1,this%ccl%struct(n)%n_
                this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            
            ! Increment monitoring variables
            this%vof_tf_drop=this%vof_tf_drop+dvol(n)
            this%np_drop=this%np_drop+1
            this%lp%np_new=this%lp%np_new+1
            this%lp%vp_new=this%lp%vp_new+dvol(n)

        end if
        
        end do
        
        ! Synchronize VF fields
        call this%vf%sync_interface()
        call this%vf%clean_irl_and_band()
        
        ! Synchronize particles
        call this%lp%sync()
        
        ! Deallocate all but work array
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
   

    subroutine transfer_films(this)
        use irl_fortran_interface
        use messager,  only: die
        use vfs_class, only: VFlo,VFhi
        use mathtools, only: Pi,normalize,cross_product
        use random,    only: random_uniform,random_gamma
        use mpi_f08
        use parallel,  only: MPI_REAL_WP
        implicit none
        class(ligament), intent(inout) :: this
        real(WP), dimension(:), allocatable :: fvol
        real(WP), dimension(:), allocatable :: fthc
        real(WP), dimension(:), allocatable :: frem
        character(len=str_medium) :: filename
        real(WP), dimension(:), allocatable :: sort_ke
        integer, dimension(:), allocatable ::  sort_id,plist,dispels
        real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
        integer :: n,nn,m,i,j,k,ii,jj,kk,ncell_,tmp_id,l,totalnewp,np_start,np_old,count,ierr,ind,ip,iunit,rank
        real(WP)  :: tmp_ke,curv_sum,ncurv,Vt,Vl,Vd,alpha,beta
        real(WP), dimension(3) :: nref,tref,sref
        logical :: sampled

        ! Start by performing a CCL based on film criteria
        call this%ccl_film%build(make_label,same_label)
        if (this%ccl_film%nstruct.ge.1) then
        ! Allocate film stats arrays
        allocate(fvol(1:this%ccl_film%nstruct)); fvol=0.0_WP
        allocate(fthc(1:this%ccl_film%nstruct)); fthc=HUGE(alpha)!5.0_WP*this%cfg%min_meshsize
        allocate(frem(1:this%ccl_film%nstruct)); frem=0.0_WP

        ! Get local thickness of the film to determine if film should be convereted
        call this%vf%get_thickness()

        ! First pass to accumulate volume and get minimum thickness
        do n=1,this%ccl_film%nstruct
        ! Loop over cells in structure
        do m=1,this%ccl_film%struct(n)%n_
           ! Get cell indices
           i=this%ccl_film%struct(n)%map(1,m)
           j=this%ccl_film%struct(n)%map(2,m) 
           k=this%ccl_film%struct(n)%map(3,m)
           ! Accumulate volume
           fvol(n)=fvol(n)+this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
           ! Get minimum thickness
           fthc(n)=min(fthc(n),this%vf%thickness(i,j,k))
           ! Check if film touches auto burst layer
           if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
           &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
           &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
           &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
           &   k.ge.this%vf%cfg%kmax-this%nlayer) frem(n)=1.0_WP
        end do
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,fvol,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,fthc,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,frem,1*this%ccl_film%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

        ! Zero out monitoring variables
        this%vof_tf_film=0.0_WP
        this%np_film=0
        ! Record initial droplets in each processor for future outputing purpose
        np_start=this%lp%np_; sampled=.false.
        ! Second pass to decide if the film has reached a minimum thickness to burst
        do n=1,this%ccl_film%nstruct
        ! Min thickness below threshold and film volume greater than a threshold 
        if (fthc(n).le.this%fmin .and. fvol(n).gt.this%fnumcell*this%fmin*(this%vf%cfg%min_meshsize**2)) then
        ! Too close to the end of domain
        else if (frem(n).gt.0.0_WP) then
        else
            cycle
        end if
        ! output to confirm
        if (this%vf%cfg%amRoot) print *, "This is a thin film with min_thickness", fthc(n), "and this is id:", n ,"vol is:", fvol(n)
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
            do ii = 1, ncell_-1
                do jj = 1, ncell_-ii
                   if (sort_ke(jj).gt.sort_ke(jj+1)) then
                      ! Swap the values
                      tmp_ke = sort_ke(jj)
                      sort_ke(jj) = sort_ke(jj+1)
                      sort_ke(jj+1) = tmp_ke
                      ! Swap the corresponding IDs
                      tmp_id = sort_id(jj)
                      sort_id(jj) = sort_id(jj+1)
                      sort_id(jj+1) = tmp_id
                   end if
                end do
             end do
            Vt=0.0_WP; Vl=0.0_WP
            np_old=this%lp%np_
            do m=1,ncell_
               i=this%ccl_film%struct(n)%map(1,sort_id(m))
               j=this%ccl_film%struct(n)%map(2,sort_id(m))
               k=this%ccl_film%struct(n)%map(3,sort_id(m))
               ! Accumulate 
               Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               if (.not.sampled) then
                  ! Get droplet information based on localized curvature
                  curv_sum=0.0_WP; ncurv=0.0_WP
                  do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                        curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                        ncurv=ncurv+1.0_WP
                     end if
                  end do
                  call bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum)
                  Vd = pi/6.0_WP*(min(random_gamma(alpha)*beta*this%fd0,2.0_WP*this%frp))**3
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
                  this%lp%np_=this%lp%np_+1
                  ! Make room for new drop
                  call this%lp%resize(this%lp%np_)
                  ! Add the drop
                  this%lp%p(this%lp%np_)%id  =int(6,8)                                   
                  this%lp%p(this%lp%np_)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            
                  this%lp%p(this%lp%np_)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref
                  this%lp%p(this%lp%np_)%vel =this%cfg%get_velocity(pos=this%lp%p(this%lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                  this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
                  this%lp%p(this%lp%np_)%flag=0                                          
                  this%lp%p(this%lp%np_)%dt  =0.0_WP                                     
                  this%lp%p(this%lp%np_)%Acol=0.0_WP                                     
                  this%lp%p(this%lp%np_)%Tcol=0.0_WP  
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
                  sampled = .false.

                  ! Increment monitoring variables
                  this%vof_tf_film=this%vof_tf_film+Vd
                  this%np_film=this%np_film+1
                  this%lp%np_new=this%lp%np_new+1
                  this%lp%vp_new=this%lp%vp_new+Vd
               end if
               ! Remove liquid in that cell
               this%vf%VF(i,j,k)=0.0_WP
            end do
  
            deallocate(sort_id,sort_ke)
            ! If for some reason a film with 0 liquid volume has been tagged, skip it
            if (Vt.eq.0.0_WP .and. Vl.eq.0.0_WP) cycle
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Increment particle counter
               this%lp%np_=this%lp%np_+1
               ! Make room for new drop
               call this%lp%resize(this%lp%np_)
               ! Add the drop
               this%lp%p(this%lp%np_)%id  =int(3,8)                                   
               this%lp%p(this%lp%np_)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            
               this%lp%p(this%lp%np_)%pos =this%vf%Lbary(:,i,j,k)                     
               this%lp%p(this%lp%np_)%vel =this%cfg%get_velocity(pos=this%lp%p(this%lp%np_)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
               this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(this%lp%np_)%flag=0                                          
               this%lp%p(this%lp%np_)%dt  =0.0_WP                                     
               this%lp%p(this%lp%np_)%Acol =0.0_WP                                    
               this%lp%p(this%lp%np_)%Tcol =0.0_WP                                    
               ! Increment monitoring variables
               this%lp%np_new=this%lp%np_new+1
               this%np_film=this%np_film+1
            else ! Some particles were created, make them all larger
               do ip=np_old+1,this%lp%np_
                  this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if
            ! Increment monitoring variables
            this%vof_tf_film=this%vof_tf_film+Vl
            this%lp%vp_new=this%lp%vp_new+Vl
         end if
         end do
        ! Gather the number of newly generated particles from each processor due to film burst
        totalnewp = 0
        allocate(plist(0:this%vf%cfg%nproc-1))
        ! Get number of particle generated for each processor
        call MPI_AllGATHER(this%np_film,1,MPI_INTEGER,plist,1,MPI_INTEGER,this%vf%cfg%comm,ierr)
        totalnewp= sum(plist)
        ! If there is any particle generated
        if (totalnewp .gt. 0) then
            allocate(pinfo_(1:9,1:this%np_film))
            allocate(pinfo(1:9,1:totalnewp))
            allocate(dispels(0:this%vf%cfg%nproc-1))
            ! Get info
            do ip = np_start+1, this%lp%np_
                pinfo_(1,ip-np_start)=this%lp%p(ip)%d
                pinfo_(2:4,ip-np_start)=this%lp%p(ip)%vel
                pinfo_(5,ip-np_start)=norm2(this%lp%p(ip)%vel)
                pinfo_(6:8,ip-np_start)=this%lp%p(ip)%pos
                pinfo_(9,ip-np_start)=this%lp%p(ip)%id
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
                write(iunit,'(f24.16,1x,f24.16,1x,f24.16,1x,f24.16,1x,f24.16,f24.16,1x,f24.16,1x,f24.16,1x,I2)')pinfo(1,i),pinfo(2,i),pinfo(3,i)&
                &,pinfo(4,i),pinfo(5,i),pinfo(6,i),pinfo(7,i),pinfo(8,i),INT(pinfo(9,i))
                end do
                close(iunit)
            end if
            ! Synchronize VF fields
            call this%vf%cfg%sync(this%vf%VF)
            call this%vf%clean_irl_and_band()
            ! Synchronize particles
            call this%lp%sync()
            ! Integrate monitoring variables 
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_film,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_film    ,1,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
        end if
        deallocate(fvol,fthc,frem)
        end if 
        contains
           subroutine bag_droplet_gamma(h,R)
            implicit none
            real(WP), intent(in) :: h,R
            real(WP) :: d0,Utc,ac,b,dr,ds,Oh
            real(WP) :: mean, stdev
            ! assert h,R != 0
            ! Drop diameter
            d0=1.0_WP
            ! Retraction speed
            Utc=sqrt(2.0_WP*this%fs%sigma/this%fs%rho_l/h)
            ! Centripetal acceleration
            ac=Utc**2/R
            ! Rim diameter
            b=sqrt(this%fs%sigma/this%fs%rho_l/ac)
            ! RP droplet diameter
            this%frp=1.89_WP*b
            ! Rim Ohnesorge number
            Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b*this%fs%sigma)
            ! Satellite droplet diameter
            ds=this%frp/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
            ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
            mean=0.25_WP*(h+b+this%frp+ds)/d0
            stdev=sqrt(0.25_WP*sum(([h,b,this%frp,ds]/d0-mean)**2))
            ! Gamma distribution parameters
            alpha=(mean/stdev)**2
            beta=stdev**2/mean
          end subroutine bag_droplet_gamma
        
           !> Function that identifies cells that need a label
           logical function make_label(i,j,k)
           implicit none
           integer, intent(in) :: i,j,k
           if ((this%vf%VF(i,j,k).gt.VFlo).and.(this%vf%VF(i,j,k).lt.VFhi).and.((this%vf%norm_pos(i,j,k)-this%vf%norm_neg(i,j,k)).lt.0.5_WP).and.((this%vf%norm_pos(i,j,k)+this%vf%norm_neg(i,j,k)).ge.0.925_WP)) then
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

     subroutine transfer_ligs(this)
      use vfs_class, only: VFlo,VFhi
      use mathtools, only: pi,twoPi
      use mpi_f08
      use parallel,  only: MPI_REAL_WP
      use messager, only: die
      use irl_fortran_interface
      implicit none
      class(ligament), intent(inout) :: this
      real(WP), dimension(:)    , allocatable :: lvol
      real(WP), dimension(:)    , allocatable :: lthc
      real(WP), dimension(:)    , allocatable :: llen
      real(WP), dimension(:)    , allocatable :: lnum
      real(WP), dimension(:)    , allocatable :: lper
      real(WP), dimension(:,:)  , allocatable :: lpos
      real(WP), dimension(:,:)  , allocatable :: lvel
      real(WP), dimension(:,:,:), allocatable :: lmoi
      real(WP), dimension(:)    , allocatable :: lrem
      real(WP), dimension(:)    , allocatable :: lSR
      real(WP), dimension(:)    , allocatable :: xmin,xmax,ymin,ymax,zmin,zmax
      integer :: n,m,ierr,i,j,k,l,ii,jj,kk,iunit,totalnewp,np_start,np_old,count,ip,rank
      real(WP) :: x,y,z,x0,y0,z0,lmax,lmid,lmin
      character(len=str_medium) :: filename
      integer, dimension(:), allocatable ::  plist,dispels
      real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
      real(WP) :: Vt,Vl,Vd,minor_radius,diam,Vrim,Lrim
      real(WP) :: Oh,Trp,Lrp,Tsr,SR_tmp
      real(WP), dimension(1:3) :: tangent
      real(WP), dimension(:,:,:,:), allocatable :: SR
      integer  :: nmain,nsat
      real(WP), dimension(:,:,:), allocatable :: thickness
      integer,  dimension(:,:,:), allocatable :: struct_type
      ! Moment of inertia calculation using lapack
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                            !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info

      ! Query optimal work array size
      if (.not.allocated(work)) then
      call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
      lwork=int(lwork_query(1)); allocate(work(lwork))
      end if

      ! Get thickness and local struct_type for global information calculation
      allocate(thickness  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));thickness=0.0_WP
      allocate(struct_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));struct_type=0
      call get_liginfo()

      this%thickness=thickness
      this%struct_type=struct_type*1.0_WP
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
          i=this%ccl_lig%struct(n)%map(1,m)
          j=this%ccl_lig%struct(n)%map(2,m)
          k=this%ccl_lig%struct(n)%map(3,m)
          ! Get cell position, accounting for periodicity
          x=this%vf%cfg%xm(i)-this%ccl_lig%struct(n)%per(1)*this%vf%cfg%xL
          y=this%vf%cfg%ym(j)-this%ccl_lig%struct(n)%per(2)*this%vf%cfg%yL
          z=this%vf%cfg%zm(k)-this%ccl_lig%struct(n)%per(3)*this%vf%cfg%zL
          ! Accumulate volume and position. Get min thickness and ligament percentage
          lvol(n  )=lvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
          lpos(n,:)=lpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
          lvel(n,:)=lvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
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
      x0=lpos(n,1)/lvol(n)
      y0=lpos(n,2)/lvol(n)
      z0=lpos(n,3)/lvol(n)
      ! Loop over cells in structure
      do m=1,this%ccl_lig%struct(n)%n_
          ! Get cell indices
          i=this%ccl_lig%struct(n)%map(1,m)
          j=this%ccl_lig%struct(n)%map(2,m)
          k=this%ccl_lig%struct(n)%map(3,m)
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
      call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
      d=max(0.0_WP,d)    
      ! Replace with corrected eigenvectors for future ligament droplet placement
      lmoi(n,:,:)=A
      ! Get characteristic lengths of drop
      lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/lvol(n))
      lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/lvol(n))
      lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/lvol(n))
      if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
      ! Use max of bounding box and MoI-derived lengths as length
      !hypot(hypot(xmax(n)-xmin(n),ymax(n)-ymin(n))**2,zmax(n)-zmin(n))
      llen(n) = max(sqrt((xmax(n)-xmin(n))**2+(ymax(n)-ymin(n))**2+(zmax(n)-zmin(n))**2),lmax)

      ! With the tangent direction of the ligament, we can evaluate the strain rate of each cell of the ligament
      tangent = lmoi(n,:,1)
      do m=1,this%ccl_lig%struct(n)%n_
         ! Get cell indices
         i=this%ccl_lig%struct(n)%map(1,m)
         j=this%ccl_lig%struct(n)%map(2,m)
         k=this%ccl_lig%struct(n)%map(3,m)
         SR_tmp =SR(1,i,j,k)*tangent(1)**2        +SR(2,i,j,k)*tangent(2)**2        +SR(3,i,j,k)*tangent(3)**2 + &
       & 2.0_WP*(SR(4,i,j,k)*tangent(1)*tangent(2)+SR(5,i,j,k)*tangent(2)*tangent(3)+SR(6,i,j,k)*tangent(1)*tangent(3))
         lSR(n) = max(lSR(n),abs(SR_tmp))
      end do
      end do
      ! Find the maximum tangential strain rate of each ligament
      call MPI_ALLREDUCE(MPI_IN_PLACE,lSR,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

      ! Zero out monitoring variables
      this%vof_tf_lig=0.0_WP
      this%np_lig=0
      ! Record initial droplets in each processor for future outputing purpose
      np_start=this%lp%np_
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
      if ((lthc(n).le.this%lmin*this%cfg%min_meshsize).and.(lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.this%lper).and.(Trp.le.Tsr).and.(nmain.ge.1)) then
      else if(lrem(n).gt.0.0_WP) then
      else
         cycle
      end if
      
      if (this%vf%cfg%amRoot) print *, "This is the min_thickness", lthc(n), ",lig percentage:", lper(n),"max length:",llen(n),&
      & "how many cells",lnum(n), "vol:",lvol(n),"nmain", nmain, "Trp:", Trp, "Tsr:", Tsr, "Trp/Tsr", Trp/Tsr,"and id:", n
      
      nsat=nmain+1
      diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)

      ! Only the main processor is in charge of creating droplets
      if (this%cfg%amRoot) then
         Lrp = twoPi*minor_radius/this%dw
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
         do l=1,nsat+nmain
            ! Increment particle counter
            this%lp%np_=this%lp%np_+1
            ! Make room for new drop
            call this%lp%resize(this%lp%np_)
            ! Add the drop
            this%lp%p(this%lp%np_)%id  =int(11,8)                                                                               
            if (mod(l,2).eq.1) then
               this%lp%p(this%lp%np_)%d=diam*this%size_ratio                                                                                    
            else
               this%lp%p(this%lp%np_)%d=diam                                                                                    
            end if
            this%lp%p(this%lp%np_)%pos =lpos(n,:)+0.5_WP*Lrp*(l-(nmain+1))*lmoi(n,:,1)
            this%lp%p(this%lp%np_)%vel =lvel(n,:)
            this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])     
            this%lp%p(this%lp%np_)%flag=0                                                                                        
            this%lp%p(this%lp%np_)%dt  =0.0_WP                                                                                  
            this%lp%p(this%lp%np_)%Acol=0.0_WP                                                                                  
            this%lp%p(this%lp%np_)%Tcol=0.0_WP
            ! Output diameter, velocity, and position
            write(iunit,*) this%lp%p(this%lp%np_)%d,this%lp%p(this%lp%np_)%vel(1),this%lp%p(this%lp%np_)%vel(2),this%lp%p(this%lp%np_)%vel(3),&
            &norm2([this%lp%p(this%lp%np_)%vel(1),this%lp%p(this%lp%np_)%vel(2),this%lp%p(this%lp%np_)%vel(3)]),this%lp%p(this%lp%np_)%pos(1),&
            &this%lp%p(this%lp%np_)%pos(2),this%lp%p(this%lp%np_)%pos(3),this%lp%p(this%lp%np_)%id  
         end do
         ! Close the file
         close(iunit)
         ! Increment monitoring variables
         this%lp%np_new=this%lp%np_new+nmain+nsat
         this%np_lig=this%np_lig+nmain+nsat
         this%vof_tf_lig=this%vof_tf_lig+lvol(n)
         this%lp%vp_new=this%lp%vp_new+lvol(n)
      end if
         ! empty out the VF
         do m=1,this%ccl_lig%struct(n)%n_
            i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
            this%vf%VF(i,j,k)=0.0_WP
         end do    
      ! end if
      end do
      ! Synchronize VF fields
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      ! Synchronize particles
      call this%lp%sync()
      ! Integrate monitoring variables 
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_lig,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_lig    ,1,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
      ! end if
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
                  ! Calculate local struct type
                  call dsyev('V','U',3,A,3,d,work,lwork,info)
                  d=max(0.0_WP,d)
                  if (d(3).gt.(this%lstratio*d(1))) struct_type(i,j,k)=struct_type(i,j,k)+1
                  if (d(3).gt.(this%lstratio*d(2))) struct_type(i,j,k)=struct_type(i,j,k)+1
               end do 
            end do 
         end do
         call this%vf%cfg%sync(thickness)
         call this%vf%cfg%sync(struct_type)
      end subroutine

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
   

   !> Initialization of ligament simulation
   subroutine init(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='ligament.input')
      end block read_input
      
      ! Create the ligament mesh
      create_config: block
         use parallel,    only: group
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xlig
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1)); call this%input%read('X ligament',xlig)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xlig
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='Ligament')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%thickness  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%struct_type  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: remap,VFlo,VFhi,plicnet,r2pnet,r2p
         use mms_geom,  only: cube_refine_vol
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2p,transport_method=remap,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         this%vf%flotsam_thld=0.0_WP
         this%vf%maxcurv_times_mesh=1.0_WP
         ! Initialize the interface to a drop/ligament
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_droplet,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a multiphase flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         use tpns_class,      only: dirichlet,clipped_neumann,bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call this%input%read('Density ratio',this%fs%rho_l)
         call this%input%read('Reynolds number',this%fs%visc_g); this%fs%visc_g=1.0_WP/this%fs%visc_g
         call this%input%read('Viscosity ratio',this%fs%visc_l); this%fs%visc_l=this%fs%visc_g*this%fs%visc_l
         call this%input%read('Weber number',this%fs%sigma); this%fs%sigma=1.0_WP/this%fs%sigma
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=1.0_WP
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block create_flow_solver
      
       ! Prepare Lagrangian drop model
      prepare_transfer: block
         use messager,  only: die
         use filesys,  only: makedir,isdir
         integer :: ierr,iunit
         character(len=str_medium) :: filename
         ! Is transfer used?
         call this%input%read('Transfer drops',this%use_drop_transfer,default=.true.)
         call this%input%read('Transfer films',this%use_film_transfer,default=.true.)
         call this%input%read('Transfer ligaments',this%use_lig_transfer,default=.true.)
         ! Only initialize transfer model if used
         if (this%use_drop_transfer) then
            ! Create CCL
            call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
            ! Set parameters for transfer
            this%ddel=0.2_WP*this%cfg%min_meshsize
            this%dmin=1.5_WP*this%cfg%min_meshsize
            this%dmax=7.0e-1_WP!1.0e-3_WP
            this%emax=0.75_WP
            ! Zero out monitoring variables
            this%vof_tf_drop=0.0_WP
            this%np_drop=0
         end if

         if (this%use_film_transfer) then
            ! Create CCL FILM
            call this%ccl_film%initialize(pg=this%cfg%pgrid,name='ccl_film')
            this%frp=0.0_WP
            this%fd0 =1.0_WP
            this%fmin=1.0e-3_WP
            this%fnumcell=100.0_WP
            ! Zero out monitoring variables
            this%vof_tf_film=0.0_WP
            this%np_film=0
         end if

         if (this%use_lig_transfer) then
            ! Create CCL LIG
            call this%ccl_lig%initialize(pg=this%cfg%pgrid,name='ccl_lig')
            this%ldmin=1.0e-2_WP
            this%dw =0.697_WP
            this%size_ratio=0.015_WP!0.707_WP 
            this%lmin=1.0_WP
            this%lmake=1.5_WP
            this%lper=0.9_WP
            this%lstratio=1.5_WP
            ! Zero out monitoring variables
            this%vof_tf_lig=0.0_WP
            this%np_lig=0
         end if

         if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
            ! Create lpt solver
            this%lp=lpt(cfg=this%cfg,name='spray')
            this%lp%rho=this%fs%rho_l
            this%lp%gravity=this%fs%gravity
            this%lp%filter_width=3.5_WP*this%cfg%min_meshsize
            call this%lp%resize(0)

            if (this%lp%cfg%amroot) then
                if (.not.isdir('spray-all')) call makedir('spray-all')
                filename='spray-all/droplets'
                open(newunit=iunit,file=trim(filename),form='formatted',status='replace',access='stream',iostat=ierr)
                if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
                ! Write the header
                write(iunit,*) 'Diameter ','U ','V ','W ','Total velocity ','X ','Y ','Z ','origin','id'
                ! Close the file
                close(iunit)         
             end if

         end if
      end block prepare_transfer
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         integer :: i,j,k
         logical :: partfile_exists
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(timestamp))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
            ! Finally, handle particle I/O
            if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
                ! Check if particle file exists
                inquire(file='restart/part_'//trim(timestamp),exist=partfile_exists)
                ! If so, read it
                if (partfile_exists) call this%lp%read(filename='restart/part_'//trim(timestamp))
             end if
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,np,nplane
         this%smesh=surfmesh(nvar=8,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='ccl_film'
         this%smesh%varname(4)='norm_abs'
         this%smesh%varname(5)='norm_sig'
         this%smesh%varname(6)='ccl_lig'
         this%smesh%varname(7)='thickness_unfilt'
         this%smesh%varname(8)='struct_type'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Calculate thickness
         call this%vf%get_thickness()
         ! Populate nplane and thickness variables
         this%smesh%var(1,:)=1.0_WP
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                        this%smesh%var(3,np)=real(this%ccl_film%id(i,j,k),WP)
                        this%smesh%var(5,np)=this%vf%norm_pos(i,j,k)-this%vf%norm_neg(i,j,k)
                        this%smesh%var(4,np)=this%vf%norm_pos(i,j,k)+this%vf%norm_neg(i,j,k)
                        this%smesh%var(6,np)=real(this%ccl_lig%id(i,j,k),WP)
                        this%smesh%var(7,np)=this%thickness(i,j,k)
                        this%smesh%var(8,np)=this%struct_type(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      ! Create partmesh object for particle output
      if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
        create_pmesh: block
           integer :: i
           this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
           this%pmesh%varname(1)='radius'
           this%pmesh%varname(2)='id'
           this%pmesh%vecname(1)='velocity'
           call this%lp%update_partmesh(this%pmesh)
           do i=1,this%lp%np_
              this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
              this%pmesh%var(2,i)=this%lp%p(i)%id
              this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
           end do
        end block create_pmesh
     end if
      
      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='ligament')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_surface('plic',this%smesh)
         if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) call this%ens_out%add_particle('part',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_atom')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vof_deleted,'VOF deleted')
         call this%mfile%add_column(this%vof_tf_drop,'VOF tf drop')
         call this%mfile%add_column(this%vof_tf_film,'VOF tf film')
         call this%mfile%add_column(this%vof_tf_lig ,'VOF tf lig')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%flotsam_error,'Flotsam error')
         call this%mfile%add_column(this%vf%thinstruct_error,'Film error')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_atom')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
         ! Create particle monitor
         if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
            call this%lp%get_max()
            this%pfile=monitor(amroot=this%lp%cfg%amRoot,name='particles')
            call this%pfile%add_column(this%time%n,'Timestep number')
            call this%pfile%add_column(this%time%t,'Time')
            call this%pfile%add_column(this%lp%np,'Particle number')
            call this%pfile%add_column(this%lp%vp_tot,'Particle volume')
            call this%pfile%add_column(this%lp%np_new,'Npart new')
            call this%pfile%add_column(this%np_drop,'Npart new drop')
            call this%pfile%add_column(this%np_film,'Npart new film')
            call this%pfile%add_column(this%np_lig, 'Npart new lig')
            call this%pfile%add_column(this%lp%vp_new,'Vpart new')
            call this%pfile%add_column(this%lp%np_out,'Npart removed')
            call this%pfile%add_column(this%lp%vp_out,'Vpart removed')
            call this%pfile%add_column(this%lp%Umin,'Particle Umin')
            call this%pfile%add_column(this%lp%Umax,'Particle Umax')
            call this%pfile%add_column(this%lp%Vmin,'Particle Vmin')
            call this%pfile%add_column(this%lp%Vmax,'Particle Vmax')
            call this%pfile%add_column(this%lp%Wmin,'Particle Wmin')
            call this%pfile%add_column(this%lp%Wmax,'Particle Wmax')
            call this%pfile%add_column(this%lp%dmin,'Particle dmin')
            call this%pfile%add_column(this%lp%dmax,'Particle dmax')
            call this%pfile%write()
         end if
      end block create_monitor
      
      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep =timer(comm=this%cfg%comm,name='Timestep')
         this%tvof  =timer(comm=this%cfg%comm,name='VOFsolve')
         this%tvel  =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres =timer(comm=this%cfg%comm,name='Pressure')
         this%ttrans=timer(comm=this%cfg%comm,name='Transferdrop')
         this%tftrans=timer(comm=this%cfg%comm,name='Transferfilm')
         this%tltrans=timer(comm=this%cfg%comm,name='Transferlig')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvof%time  ,trim(this%tvof%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
         call this%timefile%add_column(this%ttrans%time,trim(this%ttrans%name))
         call this%timefile%add_column(this%tftrans%time,trim(this%tftrans%name))
         call this%timefile%add_column(this%tltrans%time,trim(this%tltrans%name))
      end block create_timing
      
      
   contains
      
      
      !> Function that localizes the x- boundary
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      
      !> Function that localizes the x+ boundary
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that defines a level set function for a droplet
      function levelset_droplet(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
      end function levelset_droplet
      
      
      !> Function that defines a level set function for a ligament
      function levelset_ligament(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2)
      end function levelset_ligament
      

   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Reset all timers and start timestep timer
      call this%tstep%reset()
      call this%tvof%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%ttrans%reset()
      call this%tftrans%reset()
      call this%tltrans%reset()
      call this%tstep%start()
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Advance lagrangian droplets
      if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
        this%resU=this%fs%rho_g
        this%resV=this%fs%visc_g
        call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)
     end if

      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
      
      ! VOF solver step
      call this%tvof%start() ! Start VOF timer
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%tvof%stop() ! Stop VOF timer
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Start velocity timer
         call this%tvel%start()
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop velocity timer and start pressure timer
         call this%tvel%stop()
         call this%tpres%start()
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         !call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop pressure timer
         call this%tpres%stop()
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! attempt transfter
      attempt_transfer : block
         ! Zero out monitoring variables
         this%lp%np_new=0
         this%lp%vp_new=0.0_WP
         call this%ttrans%start() ! Start transfer timer
         if (this%use_drop_transfer) call this%transfer_drops()
         call this%ttrans%stop() ! Stop transfer timer
         call this%tftrans%start() ! Start burst timer
         if (this%use_film_transfer) call this%transfer_films()
         call this%tftrans%stop() ! Stop burst timer
         call this%tltrans%start() ! Start burst timer
         if (this%use_lig_transfer) call this%transfer_ligs()
         call this%tltrans%stop() ! Stop burst timer
      end block attempt_transfer
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
      end block remove_vof
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surface mesh
         update_smesh: block
            use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
            integer :: i,j,k,np,nplane
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=1.0_WP
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                           this%smesh%var(3,np)=real(this%ccl_film%id(i,j,k),WP)
                           this%smesh%var(5,np)=this%vf%norm_pos(i,j,k)-this%vf%norm_neg(i,j,k)
                           this%smesh%var(4,np)=this%vf%norm_pos(i,j,k)+this%vf%norm_neg(i,j,k)
                           this%smesh%var(6,np)=real(this%ccl_lig%id(i,j,k),WP)
                           this%smesh%var(7,np)=this%thickness(i,j,k)
                           this%smesh%var(8,np)=this%struct_type(i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Update particle mesh object
         if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
            update_pmesh: block
               integer :: i
               call this%lp%update_partmesh(this%pmesh)
               do i=1,this%lp%np_
                  this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
                  this%pmesh%var(2,i)=this%lp%p(i)%id
                  this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
               end do
            end block update_pmesh 
         end if
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Stop timestep timer
      call this%tstep%stop()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) then
        call this%lp%get_max()
        call this%pfile%write()
     end if
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Finally, handle particle I/O
            if (this%use_drop_transfer.or.this%use_film_transfer.or.this%use_lig_transfer) call this%lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
   end subroutine step
   
   
   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
end module ligament_class