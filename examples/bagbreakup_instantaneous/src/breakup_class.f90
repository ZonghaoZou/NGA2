!> Break-up model class: uses cclabel to identify thin structures and break them up into droplets
module breakup_class
   use precision,     only: WP
   use config_class,  only: config
   use cclabel_class, only: cclabel
   use lpt_class,     only: lpt
   use vfs_class,     only: vfs,VFlo
   use tpns_class,    only: tpns
   use string,  only: str_medium
   use irl_fortran_interface
   implicit none
   private

   ! Expose type/methods
   public :: breakup

   !> Break-up model object
   type :: breakup

      !> Three ccl objects are needed
      type(cclabel):: ccl,ccl_film

      !> Pointers to lpt, vfs, and tpns
      class(vfs),  pointer :: vf
      class(tpns), pointer :: fs
      class(lpt),  pointer :: lp

      !> Break-up model parameters
      real(WP) :: min_filmthickness      =1.0e-3_WP
      real(WP) :: max_eccentricity       =8.0e-1_WP
      real(WP) :: d_threshold            =6.0e-1_WP
      real(WP) :: vol_convert            =0.0_WP
      real(WP) :: d_film_rp
      real(WP) :: numfilmcell            = 20.0_WP
      real(WP), dimension(:,:,:), allocatable :: struct_type            !< Tmp struct_type for output purposes
   contains
      procedure :: initialize
      procedure :: spray_statistics_setup
      procedure :: attempt_breakup
      procedure :: transfer_detached_struct
      procedure :: breakup_film_instantaneous
      procedure :: bag_droplet_gamma
      procedure :: get_cclstats
   end type breakup

   integer, dimension(:,:,:), allocatable :: tmpstruct_type
   real(WP), dimension(:,:,:), allocatable :: tmpVF,tmpthin_sensor

contains

   !> Function that identifies cells that need a label to min thickness region
   logical function make_label_film(i,j,k)
      implicit none
      integer, intent(in) :: i,j,k
      integer :: ii,jj,kk,ncell,sum_f
      ncell = 1; sum_f = 0
      ! count number of local films
      do kk = k-ncell,k+ncell
         do jj = j-ncell,j+ncell
            do ii = i-ncell,i+ncell
               if(tmpstruct_type(ii,jj,kk).eq.2) then
                  sum_f = sum_f + 1
               end if
            end do
         end do
      end do
      if ((tmpVF(i,j,k).gt.VFlo).and.((tmpstruct_type(i,j,k).eq.2).or.(tmpstruct_type(i,j,k).eq.1.and.sum_f.ge.3)).and.tmpthin_sensor(i,j,k).eq.1.0_WP) then
         make_label_film=.true.
      else
         make_label_film=.false.
      end if
   end function make_label_film

   !> Function that identifies cells that need a label to min thickness region 
   logical function make_label(i,j,k)
      implicit none
      integer, intent(in) :: i,j,k
      if(tmpVF(i,j,k).gt.VFlo)then
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

   !> Initialize the breakup model
   subroutine initialize(this,vf,fs,lp)
      implicit none
      class(breakup), intent(inout) :: this
      class(vfs),  target, intent(in) :: vf
      class(tpns), target, intent(in) :: fs
      class(lpt),  target, intent(in) :: lp
      ! Store pointers to our solvers
      this%vf=>vf
      this%fs=>fs
      this%lp=>lp
      ! Create a connected-component labeling object
      call this%ccl%initialize(pg=this%vf%cfg%pgrid,name='ccl')
      call this%ccl_film%initialize(pg=this%vf%cfg%pgrid,name='ccl_film')

      allocate(this%struct_type(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_));this%struct_type=0.0_WP
      allocate(tmpVF(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_))
      allocate(tmpstruct_type(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_))
      allocate(tmpthin_sensor(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_))
      call this%spray_statistics_setup()
   end subroutine initialize

   !> Setup spray statistics folders
   subroutine spray_statistics_setup(this)
      use messager, only: die
      implicit none
      class(breakup), intent(inout) :: this
      character(len=str_medium) :: filename
      integer :: iunit,ierr
      ! Create directory
      if (this%lp%cfg%amroot) then
         call execute_command_line('mkdir -p spray-all')
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
         ! Write the header
         write(iunit,*) 'Diameter ','U ','V ','W ','Total velocity ','X ','Y ','Z ','origin '
         ! Close the file
         close(iunit)         
      end if
   end subroutine spray_statistics_setup

   subroutine attempt_breakup(this)
      implicit none
      class(breakup), intent(inout) :: this
      ! set up tmp variables for ccl detection
      call this%vf%get_localstructtype(tmpstruct_type); this%struct_type = tmpstruct_type
      tmpVF=this%vf%VF; tmpthin_sensor=this%vf%thin_sensor
      ! First build an all encompasing ccl to remove detached liquid structures
      call this%ccl%build(make_label,same_label)
      if (this%ccl%nstruct .ge.1) call this%transfer_detached_struct()
      ! Lastly final all the existing films and perform instantenous breakup based on minium local thickness
      call this%ccl_film%build(make_label_film,same_label)
      if (this%ccl_film%nstruct .ge.1) call this%breakup_film_instantaneous()
   end subroutine attempt_breakup

   subroutine transfer_detached_struct(this)
      use mathtools, only: pi
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      use messager, only: die
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      integer  :: n,nn,i,j,k,ierr,np,ip,m,iunit
      character(len=str_medium) :: filename
      ! Stats of the ccl objects
      real(WP), dimension(:), allocatable :: x,y,z,u,v,w,vol
      real(WP), dimension(:,:), allocatable :: lengths      
      ! Varaibles determing transfer
      real(WP) :: diam,lmin,lmax,eccentricity,myint,integral
      logical :: autotransfer
      ! stats of the ligaments that will be used for modeling breakup
      allocate(x(1:this%ccl%nstruct),y(1:this%ccl%nstruct),z(1:this%ccl%nstruct));x=0.0_WP;y=0.0_WP;z=0.0_WP
      allocate(u(1:this%ccl%nstruct),v(1:this%ccl%nstruct),w(1:this%ccl%nstruct));u=0.0_WP;v=0.0_WP;w=0.0_WP
      allocate(lengths(1:this%ccl%nstruct,1:3));lengths=0.0_WP
      allocate(vol(1:this%ccl%nstruct));vol=0.0_WP;
      myint =0.0_WP; integral =0.0_WP

      call this%get_cclstats(this%ccl,x,y,z,u,v,w,vol,lengths)
      ! Loops over global list of structures and remove detached structures
      do n=1,this%ccl%nstruct
         
         diam=(6.0_WP*vol(n)/pi)**(1.0_WP/3.0_WP)
         autotransfer=.false.
         ! Test if structure is at end of domain
         if (x(n).gt.this%vf%cfg%x(this%vf%cfg%imax-10)) autotransfer=.true.
         if (.not.autotransfer) then
            ! Test if sphericity is compatible with transfer
            lmin=lengths(n,3)
            if (lmin.eq.0.0_WP) lmin=lengths(n,2) ! Handle 2D case
            lmax=lengths(n,1)
            eccentricity=sqrt(1.0_WP-lmin**2/(lmax**2+tiny(1.0_WP)))

            if (eccentricity.gt.this%max_eccentricity) cycle
            if ((diam.eq.0.0_WP).or.(diam.gt.this%d_threshold)) cycle
         end if
         
         ! Create drop from available liquid volume - only one root does that
         if (this%vf%cfg%amRoot) then
            ! Make room for new drop
            np=this%lp%np_+1; call this%lp%resize(np)
            ! Add the drop
            this%lp%p(np)%id  =int(1,8)                                                                                 
            this%lp%p(np)%dt  =0.0_WP                                                                                   
            this%lp%p(np)%Acol =0.0_WP                                                                                  
            this%lp%p(np)%Tcol =0.0_WP                                                                                  
            this%lp%p(np)%d   =diam                                                                                     
            this%lp%p(np)%pos =[x(n),y(n),z(n)] 
            this%lp%p(np)%vel =[u(n),v(n),w(n)] 
            this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])               
            this%lp%p(np)%flag=0                                                                                        
            ! Increment particle counter
            this%lp%np_=np

            !!! Write to droplet list !!!
            ! Open the file
            filename='spray-all/droplets'
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameter, velocity, and position
            write(iunit,*) this%lp%p(np)%d,this%lp%p(np)%vel(1),this%lp%p(np)%vel(2),this%lp%p(np)%vel(3),norm2([this%lp%p(np)%vel(1),this%lp%p(np)%vel(2),this%lp%p(np)%vel(3)]),this%lp%p(np)%pos(1),this%lp%p(np)%pos(2),this%lp%p(np)%pos(3),'detached',this%lp%p(np)%id  
            ! Close the file
            close(iunit)
         end if

         ! Find local structs with matching id
         do nn=1,this%ccl%struct(n)%n_
            i=this%ccl%struct(n)%map(1,nn);j=this%ccl%struct(n)%map(2,nn);k=this%ccl%struct(n)%map(3,nn)
            myint = myint + this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         vol(n) = 0.0_WP
      end do
         
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      call MPI_ALLREDUCE(myint,integral,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr) ! total converted liquid volume this timestep
      this%vol_convert = this%vol_convert + integral 
      call this%lp%sync()
      deallocate(x,y,z,u,v,w,vol,lengths)
   end subroutine transfer_detached_struct

   subroutine breakup_film_instantaneous(this)
      use messager,  only: die
      use mathtools, only: Pi,normalize,cross_product
      use random,    only: random_uniform,random_gamma
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(breakup), intent(inout) :: this
      character(len=str_medium) :: filename
      real(WP), dimension(:), allocatable :: min_thickness,min_thickness_,thickness_list
      integer, dimension(:), allocatable :: id
      real(WP), dimension(:), allocatable :: vol_,vol
      integer :: n,nn,i,j,k,l,np,np_old,np_start,total_id,rank,iunit,ierr,ip,ii,jj,temp_id
      real(WP)  :: d0,curv_sum,ncurv,alpha,beta,Vt,Vl,Vd,myint,integral,temp_thickness
      real(WP), dimension(3) :: nref,tref,sref
      logical :: hassampled
      allocate(min_thickness(1:this%ccl_film%nstruct),min_thickness_(1:this%ccl_film%nstruct))
      min_thickness = this%vf%cfg%min_meshsize; min_thickness_=5.0_WP*this%vf%cfg%min_meshsize
      allocate(vol_(1:this%ccl_film%nstruct),vol(1:this%ccl_film%nstruct));vol_=0.0_WP;vol=0.0_WP
      call this%vf%get_thickness()
      do n=1,this%ccl_film%nstruct
         do nn=1,this%ccl_film%struct(n)%n_
            i=this%ccl_film%struct(n)%map(1,nn); j=this%ccl_film%struct(n)%map(2,nn); k=this%ccl_film%struct(n)%map(3,nn)
            vol_(n) = vol_(n) + this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            min_thickness_(n) = min(min_thickness_(n),this%vf%thickness(i,j,k))
         end do
      end do
      call MPI_ALLREDUCE(vol_,vol,this%ccl_film%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(min_thickness_,min_thickness,this%ccl_film%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)

      d0=1.0_WP; np_start=this%lp%np_; hassampled=.false.
      do n=1,this%ccl_film%nstruct
         ! Breakup film instantaneously if the minimal thickness is below the criterion and the film is large enough
         if (min_thickness(n).gt.this%min_filmthickness) cycle
         if (vol(n).lt.this%numfilmcell*this%min_filmthickness*(this%vf%cfg%min_meshsize**2)) cycle
         if (this%vf%cfg%amRoot) print *, "This is a thin film with min_thickness", min_thickness(n), "and this is id:", n ,"vol is:", vol(n)
         if (this%ccl_film%struct(n)%n_ .ge.1) then
            allocate(id(1:this%ccl_film%struct(n)%n_)) 
            allocate(thickness_list(1:this%ccl_film%struct(n)%n_))
            total_id = this%ccl_film%struct(n)%n_
            do nn=1,total_id
               i=this%ccl_film%struct(n)%map(1,nn); j=this%ccl_film%struct(n)%map(2,nn); k=this%ccl_film%struct(n)%map(3,nn)
               id(nn) = nn; thickness_list(nn) = this%vf%thickness(i,j,k)
            end do
            call sort_by_thickness()

            Vt=0.0_WP; Vl=0.0_WP; np_old=this%lp%np_
            do nn=1,total_id
               i=this%ccl_film%struct(n)%map(1,id(nn)); j=this%ccl_film%struct(n)%map(2,id(nn)); k=this%ccl_film%struct(n)%map(3,id(nn))
               ! Accumulate 
               Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               ! print *, "Vl is: ", Vl, "Current id is: ", id(nn), "total id: ", total_id, "The current cell volume",this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k) ,"The current thickness: ", this%vf%thickness(i,j,k)
               if (.not.hassampled) then
                  ! Get a localized cell curvature
                  curv_sum=0.0_WP; ncurv=0.0_WP
                  do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                        curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                        ncurv=ncurv+1.0_WP
                     end if
                  end do
                  call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
                  Vd = pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,2.0_WP*this%d_film_rp))**3
                  hassampled = .true.
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
                  
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(6,8)                                   !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  this%lp%p(np)%Acol=0.0_WP                                     !< Give zero collision force
                  this%lp%p(np)%Tcol=0.0_WP                                     !< Give zero collision force
                  this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref
                  ! this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-1.0_WP*this%vf%cfg%meshsize(i,j,k),0.0_WP*this%vf%cfg%meshsize(i,j,k))*this%vf%edge_normal(:,i,j,k) ! shedd the droplet backwards
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
                  this%lp%p(np)%flag=0                                          !< Activate it
                  this%lp%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
                  hassampled = .false.
               end if
               ! Remove liquid in that cell
               myint = myint + this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               this%vf%VF(i,j,k)=0.0_WP
            end do

            deallocate(id,thickness_list)
            ! If for some reason a film with 0 liquid volume has been tagged, skip it
            if (Vt.eq.0.0_WP .and. Vl.eq.0.0_WP) cycle
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Add one last drop for remaining liquid volume
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(3,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                     !< Place the drop at the liquid barycenter
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
            else ! Some particles were created, make them all larger
               do ip=np_old+1,this%lp%np_
                  this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if

         end if
      end do
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      call MPI_ALLREDUCE(myint,integral,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)    
      this%vol_convert = this%vol_convert + integral 

      filename='spray-all/droplets'
      do rank=0,this%vf%cfg%nproc-1
         if (rank.eq.this%vf%cfg%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameters and velocities
            do i=np_start+1,this%lp%np_
               write(iunit,*) this%lp%p(i)%d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),norm2([this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]),this%lp%p(i)%pos(1),this%lp%p(i)%pos(2),this%lp%p(i)%pos(3),'film',this%lp%p(np)%id  
            end do
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%vf%cfg%comm,ierr)
      end do

      call this%lp%sync()
      deallocate(min_thickness,min_thickness_,vol_,vol)

      contains
         !> Sort film indices by increasing thickness
         subroutine sort_by_thickness()
            implicit none
            do ii = 1, total_id-1
               do jj = 1, total_id-ii
                  if (thickness_list(jj) > thickness_list(jj+1)) then
                     ! Swap the values
                     temp_thickness = thickness_list(jj)
                     thickness_list(jj) = thickness_list(jj+1)
                     thickness_list(jj+1) = temp_thickness
          
                     ! Swap the corresponding IDs
                     temp_id = id(jj)
                     id(jj) = id(jj+1)
                     id(jj+1) = temp_id
                  end if
               end do
            end do
         end subroutine sort_by_thickness
      
   end subroutine breakup_film_instantaneous

   !> Generate a Gamma distribution for bag droplet formation
   !> where the number PDF is in the form
   !> p_n(x=d;alpha,beta)=x**(alpha-1)*exp(-x/beta)/beta**alpha/gamma(alpha)
   !> Adapted from Jackiw and Ashgriz 2022, JFM
   subroutine bag_droplet_gamma(this,h,R,alpha,beta)
      implicit none
      class(breakup), intent(inout) :: this
      real(WP), intent(in) :: h,R
      real(WP), intent(out) :: alpha,beta
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
      this%d_film_rp=1.89_WP*b
      ! Rim Ohnesorge number
      Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b*this%fs%sigma)
      ! Satellite droplet diameter
      ds=this%d_film_rp/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
      ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
      mean=0.25_WP*(h+b+this%d_film_rp+ds)/d0
      stdev=sqrt(0.25_WP*sum(([h,b,this%d_film_rp,ds]/d0-mean)**2))
      ! Gamma distribution parameters
      alpha=(mean/stdev)**2
      beta=stdev**2/mean
   end subroutine bag_droplet_gamma

   subroutine get_cclstats(this,ccl,x,y,z,u,v,w,vol,lengths)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_INTEGER
      use parallel,  only: MPI_REAL_WP
      use irl_fortran_interface
      implicit none
      class(breakup), intent(inout) :: this
      class(cclabel), intent(in)::ccl
      real(WP), dimension(1:), intent(inout) :: x,y,z,u,v,w,vol
      real(WP), dimension(1:,1:), intent(inout) :: lengths
      integer :: n,nn,nnn,i,j,k,ii,jj,kk,ierr,np,ip,m,iunit,rank,per_x,per_y,per_z
      ! Allocate variables to get stats
      real(WP), dimension(:), allocatable :: vol_,x_vol_,y_vol_,z_vol_
      real(WP), dimension(:), allocatable :: u_vol_,v_vol_,w_vol_
      real(WP), dimension(:), allocatable :: x_min_,x_min,x_max_,x_max
      real(WP), dimension(:), allocatable :: y_min_,y_min,y_max_,y_max
      real(WP), dimension(:), allocatable :: z_min_,z_min,z_max_,z_max
      real(WP), dimension(:,:,:), allocatable :: Imom_,Imom
      real(WP) :: xtmp,ytmp,ztmp
      ! Moment of inertia variable
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer , parameter :: order = 3
      integer  :: lwork,info

      allocate(vol_(1:ccl%nstruct));vol_=0.0_WP
      allocate(x_min(1:ccl%nstruct),x_min_(1:ccl%nstruct));x_min=0.0_WP;x_min_= 10000.0_WP
      allocate(x_max(1:ccl%nstruct),x_max_(1:ccl%nstruct));x_max=0.0_WP;x_max_=-10000.0_WP
      allocate(y_min(1:ccl%nstruct),y_min_(1:ccl%nstruct));y_min=0.0_WP;y_min_= 10000.0_WP
      allocate(y_max(1:ccl%nstruct),y_max_(1:ccl%nstruct));y_max=0.0_WP;y_max_=-10000.0_WP
      allocate(z_min(1:ccl%nstruct),z_min_(1:ccl%nstruct));z_min=0.0_WP;z_min_= 10000.0_WP
      allocate(z_max(1:ccl%nstruct),z_max_(1:ccl%nstruct));z_max=0.0_WP;z_max_=-10000.0_WP
      allocate(x_vol_(1:ccl%nstruct),y_vol_(1:ccl%nstruct),z_vol_(1:ccl%nstruct));x_vol_=0.0_WP;y_vol_=0.0_WP;z_vol_=0.0_WP
      allocate(u_vol_(1:ccl%nstruct),v_vol_(1:ccl%nstruct),w_vol_(1:ccl%nstruct));u_vol_=0.0_WP;v_vol_=0.0_WP;w_vol_=0.0_WP
      allocate(Imom(1:ccl%nstruct,3,3),Imom_(1:ccl%nstruct,3,3));Imom=0.0_WP;Imom_=0.0_WP
      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      do n=1,ccl%nstruct
         ! Periodicity
         per_x = ccl%struct(n)%per(1); per_y = ccl%struct(n)%per(2); per_z = ccl%struct(n)%per(3)
         ! get number of local cells
         do nn=1,ccl%struct(n)%n_
            i=ccl%struct(n)%map(1,nn); j=ccl%struct(n)%map(2,nn); k=ccl%struct(n)%map(3,nn)
            ! Location of struct node
            xtmp = this%vf%cfg%xm(i)-per_x*this%vf%cfg%xL
            ytmp = this%vf%cfg%ym(j)-per_y*this%vf%cfg%yL
            ztmp = this%vf%cfg%zm(k)-per_z*this%vf%cfg%zL
            ! Volume
            vol_(n) = vol_(n) + this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            ! Center of gravity
            x_vol_(n) = x_vol_(n) + xtmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            y_vol_(n) = y_vol_(n) + ytmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            z_vol_(n) = z_vol_(n) + ztmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            ! Average gas velocity inside struct
            u_vol_(n) = u_vol_(n) + this%fs%U(i,j,k)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            v_vol_(n) = v_vol_(n) + this%fs%V(i,j,k)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            w_vol_(n) = w_vol_(n) + this%fs%W(i,j,k)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
         end do
      end do
      ! Sum parallel stats
      call MPI_ALLREDUCE(vol_,vol,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(x_vol_,x,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_vol_,y,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_vol_,z,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(u_vol_,u,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(v_vol_,v,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(w_vol_,w,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      do n=1,ccl%nstruct
         ! Periodicity
         per_x = ccl%struct(n)%per(1); per_y = ccl%struct(n)%per(2); per_z = ccl%struct(n)%per(3)
         do nn=1,ccl%struct(n)%n_
            ! Indices of struct node
            i=ccl%struct(n)%map(1,nn); j=ccl%struct(n)%map(2,nn); k=ccl%struct(n)%map(3,nn)
            xtmp = this%vf%cfg%xm(i)-per_x*this%vf%cfg%xL-x(n)/vol(n)
            ytmp = this%vf%cfg%ym(j)-per_y*this%vf%cfg%yL-y(n)/vol(n)
            ztmp = this%vf%cfg%zm(k)-per_z*this%vf%cfg%zL-z(n)/vol(n)
            ! Moment of Inertia
            Imom_(n,1,1) = Imom_(n,1,1) + (ytmp**2 + ztmp**2)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            Imom_(n,2,2) = Imom_(n,2,2) + (xtmp**2 + ztmp**2)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            Imom_(n,3,3) = Imom_(n,3,3) + (xtmp**2 + ytmp**2)*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            Imom_(n,1,2) = Imom_(n,1,2) - xtmp*ytmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            Imom_(n,1,3) = Imom_(n,1,3) - xtmp*ztmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            Imom_(n,2,3) = Imom_(n,2,3) - ytmp*ztmp*this%vf%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            do nnn=1,2
               if (getNumberOfVertices(this%vf%interface_polygon(nnn,i,j,k)).gt.0) then
                  d = calculateCentroid(this%vf%interface_polygon(nnn,i,j,k))
                  x_min_(n) = min(x_min_(n),d(1)); x_max_(n) = max(x_max_(n),d(1))
                  y_min_(n) = min(y_min_(n),d(2)); y_max_(n) = max(y_max_(n),d(2))
                  z_min_(n) = min(z_min_(n),d(3)); z_max_(n) = max(z_max_(n),d(3))
               end if
            end do
         end do 
      end do

      ! Sum parallel stat on Imom
      do i=1,3
         do j=1,3
            call MPI_ALLREDUCE(Imom_(:,i,j),Imom(:,i,j),ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         end do
      end do
      ! Get extents
      call MPI_ALLREDUCE(x_min_,x_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(x_max_,x_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_min_,y_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_max_,y_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_min_,z_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_max_,z_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      ! Store data
      do n=1,ccl%nstruct
         ! Center of gravity
         x(n) = x(n)/vol(n); y(n) = y(n)/vol(n); z(n) = z(n)/vol(n)
         ! Periodicity: transport back inside domain if needed
         if (x(n).lt.this%vf%cfg%x(this%vf%cfg%imin)) x(n) = x(n)+this%vf%cfg%xL
         if (y(n).lt.this%vf%cfg%y(this%vf%cfg%jmin)) y(n) = y(n)+this%vf%cfg%yL
         if (z(n).lt.this%vf%cfg%z(this%vf%cfg%kmin)) z(n) = z(n)+this%vf%cfg%zL
         u(n)=u(n)/vol(n); v(n)=v(n)/vol(n); w(n)=w(n)/vol(n)
         ! Eigenvalues/eigenvectors of moments of inertia tensor
         A = Imom(n,:,:); nnn = 3
         ! On exit, A contains eigenvectors, and d contains eigenvalues in ascending order
         call dsyev('V','U',nnn,A,nnn,d,work,lwork,info)
         ! Get rid of very small negative values (due to machine accuracy)
         d = max(0.0_WP,d)
         ! Store characteristic lengths
         lengths(n,1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol(n))
         lengths(n,2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol(n))
         lengths(n,3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol(n))
         ! Zero out length in 3rd dimension if 2D
         if (this%vf%cfg%nx.eq.1.or.this%vf%cfg%ny.eq.1.or.this%vf%cfg%nz.eq.1) lengths(n,3)=0.0_WP
      end do
      ! Deallocate arrays
      deallocate(vol_,x_vol_,y_vol_,z_vol_,u_vol_,v_vol_,w_vol_,Imom_,Imom)
      deallocate(x_min_,y_min_,z_min_,x_max_,y_max_,z_max_)
      deallocate(x_min,y_min,z_min,x_max,y_max,z_max)
   end subroutine get_cclstats
end module
