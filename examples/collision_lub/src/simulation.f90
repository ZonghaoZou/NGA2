!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use cclabel_class,     only: cclabel
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   !type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   type(cclabel)     :: ccl

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final,get_gasP,apply_gasP,get_thickness!,solveUs!,getGP,solveUs
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   ! real(WP), dimension(:,:,:), allocatable :: Us,Vs,Ws
   ! real(WP), dimension(:,:,:), allocatable :: Usold,Vsold,Wsold
   real(WP), dimension(:,:,:), allocatable :: Pg,Pd
   real(WP), dimension(:,:,:), allocatable :: radialU,verticalU
   real(WP), dimension(:,:,:), allocatable :: thickness_old,thickness_new
   ! real(WP), dimension(:,:,:), allocatable :: rho
   ! real(WP), dimension(:,:,:), allocatable :: alpha_x,alpha_y,alpha_z
   
   !> Problem definition
   real(WP), dimension(3) :: center1,center2,vel1,vel2
   real(WP) :: radius1,radius2
   ! real(WP), parameter :: HamakerC=-19.4318  ! Written in log form!3.7e-20_WP
   real(WP), parameter :: HamakerC=3.7e-20_WP  ! Written in log form!
   
contains
   
   subroutine get_gasP
      use irl_fortran_interface
      use mpi_f08!,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE,MPI_MIN
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi
      use messager,  only: die
      implicit none
      ! Parameters for moment of inertia
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                            !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info
      integer :: total_cell,ierr,rank,count,i,j,k,m,n,l,rind,ni
      real(WP) :: x,y,z,x0,y0,z0,lambdaAir,Kn,myvol,u1,u2,u3,xr,yr,rmag,signmeasure,voltmp,ptmp,pmin_,pmin
      real(WP), dimension(3) :: mybary,myvel
      real(WP), dimension(4) :: plane
      real(WP), dimension(5) :: tmprow
      integer, dimension(:), allocatable:: plist,dispels
      real(WP), dimension(:)    , allocatable :: dgvol,rdhdt_int!,dr0,dP0,
      real(WP), dimension(:,:)  , allocatable :: dgpos,gpinfo_,gpinfo 
      real(WP), dimension(:,:,:), allocatable :: dmoi
      ! Output to learn what's going on
      radialU=0.0_WP;verticalU=0.0_WP
      ! First get the global minium gas pressure as my base reference
      pmin_=huge(1.0_WP);pmin=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (cfg%VF(i,j,k).gt.0.0_WP .and. vf%VF(i,j,k).eq.0.0_WP) pmin_=min(pmin_,fs%P(i,j,k) )
            end do
         end do
      end do
      call MPI_ALLREDUCE(pmin_,pmin,1,MPI_REAL_WP,MPI_MIN,cfg%comm,ierr)
      Pg=pmin
      lambdaAir= 69e-9_WP
      ! Query optimal work array size
      if (.not.allocated(work)) then
         call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
         lwork=int(lwork_query(1)); allocate(work(lwork))
      end if

      ! Build ccl to get the thin gas region for caculating gas pressure
      call ccl%build(make_label,same_label)

      ! Allocate fields for calculation
      allocate(dgvol(1:ccl%nstruct        )); dgvol=0.0_WP
      allocate(dgpos(1:ccl%nstruct,1:3    )); dgpos=0.0_WP
      allocate(dmoi(1:ccl%nstruct,1:3,1:3)); dmoi=0.0_WP
      ! allocate(dr0(1:ccl%nstruct        )); dr0=0.0_WP
      ! allocate(dP0(1:ccl%nstruct        )); dP0=huge(1.0_WP)
      allocate(plist(0:cfg%nproc-1))
      allocate(dispels(0:cfg%nproc-1))
      ! First pass to accumulate position for moment of inertia
      do n=1,ccl%nstruct
         ! Loop over cells in structure
         do m=1,ccl%struct(n)%n_
            ! Get cell indices
            i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
            ! Get cell position, accounting for periodicity
            x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL
            y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL
            z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL

            ! Accumulate volume and position
            dgvol(n  )=dgvol(n  )+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))
            dgpos(n,:)=dgpos(n,:)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*[x,y,z]
        end do 
      end do 
      call MPI_ALLREDUCE(MPI_IN_PLACE,dgvol,1*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dgpos,3*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)

      ! Second pass to accumulate moment of inertia
      do n=1,ccl%nstruct
         ! Get the region gas barycenter
         x0=dgpos(n,1)/dgvol(n)
         y0=dgpos(n,2)/dgvol(n)
         z0=dgpos(n,3)/dgvol(n)
         ! Loop over cells in structure
         do m=1,ccl%struct(n)%n_
             ! Get cell indices
             i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
             ! Get cell position relative to drop barycenter, accounting for periodicity
             x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL-x0
             y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL-y0
             z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL-z0
             ! Accumulate moment of inertia
             dmoi(n,2,2)=dmoi(n,2,2)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(z**2+x**2)
             dmoi(n,3,3)=dmoi(n,3,3)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x**2+y**2)
             dmoi(n,1,1)=dmoi(n,1,1)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(y**2+z**2)
             dmoi(n,1,2)=dmoi(n,1,2)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x*y)
             dmoi(n,1,3)=dmoi(n,1,3)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x*z)
             dmoi(n,2,3)=dmoi(n,2,3)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(y*z)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      do n=1,ccl%nstruct
         ! Get the moi directions
         A=dmoi(n,:,:)
         call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
         dmoi(n,:,:)=A ! dmoi(n,:,2) and dmoi(n,:,3) are the two principle axes marking the tagential plane
         ! allocate the fields needed for getting the gas pressure
         ! 1 is r, 2 is h, 3 is deltP, 4 is ur, 5 is uz
         allocate(gpinfo_(1:ccl%struct(n)%n_,1:5));gpinfo_=0.0_WP
         ! Fill out all the info of the required parameters
         do m=1,ccl%struct(n)%n_
            ! Get cell indices
            i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
            ! Get cell position relative to drop barycenter, accounting for periodicity
            x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL
            y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL
            z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL
            ! Get the r magnitude
            xr=dot_product([x,y,z],dmoi(n,:,1)); yr=dot_product([x,y,z],dmoi(n,:,2))
            gpinfo_(m,1)=sqrt(xr**2+yr**2)
            gpinfo_(m,2)=thickness_new(i,j,k)
            Kn = lambdaAir/thickness_new(i,j,k)
            ! gpinfo_(m,2)=vf%thickness(i,j,k)
            ! Kn = lambdaAir/vf%thickness(i,j,k)
            gpinfo_(m,3)=1.0_WP+6.88_WP*Kn+6.0_WP*Kn*LOG(1.0_WP+2.76_WP*Kn+0.127_WP*Kn**2)/pi
            ! ! Record pressure as P_0 that has the largest r
            ! if (gpinfo_(m,1).gt.dr0(n)) then
            !    dr0(n) = gpinfo_(m,1)
            !    dP0(n) = fs%P(i,j,k)
            ! end if
            voltmp=0.0_WP
            ! For each interface
            do ni=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
               if (getNumberOfVertices(vf%interface_polygon(ni,i,j,k)).ne.0) then
                  ! add the surface area of the polygon
                  myvol=abs(calculateVolume(vf%interface_polygon(ni,i,j,k)))
                  ! Get the barycetner of the polygon
                  mybary=calculateCentroid(vf%interface_polygon(ni,i,j,k))
                  ! get my velocity based on my polygon barycenter
                  myvel=cfg%get_velocity(mybary,i,j,k,fs%U,fs%V,fs%W)
                  u1=dot_product(myvel,dmoi(n,:,1)); u2=dot_product(myvel,dmoi(n,:,2)); u3=dot_product(myvel,dmoi(n,:,3))
                  gpinfo_(m,4)=gpinfo_(m,4)+myvol*(sqrt(u1**2+u2**2))
                  ! Using the current interface as a measure for film expansion direction
                  ! plane=getPlane(vf%liquid_gas_interface(i,j,k),ni-1)
                  ! signmeasure=dot_product(plane(1:3),myvel)
                  ! ! If the velocity is aligned with outward normal, the gas film thickness decreases, therefore negative
                  ! if (signmeasure.gt.0) then
                  !    gpinfo_(m,5)=gpinfo_(m,5)-myvol*(abs(u3))
                  ! else
                  !    gpinfo_(m,5)=gpinfo_(m,5)+myvol*(abs(u3))
                  ! end if
                  ! gpinfo_(m,5)=(vf%thickness(i,j,k)-thickness_old(i,j,k))/time%dt
                  
                  ! adding the volume of each polygon
                  voltmp=voltmp+myvol
               end if
            end do
            gpinfo_(m,5)=(thickness_new(i,j,k)-thickness_old(i,j,k))/time%dt
            ! surface area weighted average of the velocities
            gpinfo_(m,4)=gpinfo_(m,4)/voltmp
            ! gpinfo_(m,5)=2.0_WP*gpinfo_(m,5)/voltmp
            radialU(i,j,k)=gpinfo_(m,4)
            verticalU(i,j,k)=gpinfo_(m,5)
         end do
         ! call MPI_ALLREDUCE(MPI_IN_PLACE,dr0,1*ccl%nstruct,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         ! call MPI_ALLREDUCE(MPI_IN_PLACE,dP0,1*ccl%nstruct,MPI_REAL_WP,MPI_MIN,cfg%comm,ierr)
         call MPI_AllGATHER(ccl%struct(n)%n_,1,MPI_INTEGER,plist,1,MPI_INTEGER,cfg%comm,ierr)
         total_cell=sum(plist)
         if (total_cell.gt.0) then
            allocate(gpinfo(1:total_cell,1:5));gpinfo=0.0_WP
            allocate(rdhdt_int(1:total_cell));rdhdt_int=0.0_WP
            ! Calculate dispels
            count = 0
            do rank=0,cfg%nproc-1
               dispels(rank) = count
               count = count + plist(rank)
            end do
            ! Communicate to root
            do i = 1,5
               call MPI_ALLGATHERV(gpinfo_(:,i),ccl%struct(n)%n_,MPI_REAL_WP,gpinfo(:,i),plist,dispels,MPI_REAL_WP,cfg%comm)
            end do
            ! Now I want to bubble sort them in ascending order
            do i = 1, total_cell-1
               do j = i+1, total_cell
                  if (gpinfo(j,1) < gpinfo(i,1)) then
                     tmprow = gpinfo(i,:)     ! swap whole rows
                     gpinfo(i,:)  = gpinfo(j,:)
                     gpinfo(j,:)  = tmprow
                  end if
               end do
            end do

            do l=2,total_cell
               do m=1,l-1
                  ! integral of r*dh/dt
                  rdhdt_int(l)=rdhdt_int(l)+0.5_WP*(gpinfo(m,1)*gpinfo(m,5)+gpinfo(m+1,1)*gpinfo(m+1,5))*(gpinfo(m+1,1)-gpinfo(m,1))
               end do 
            end do

            ! Pg =Pg+dP0(n)
            ! Based on the gpinfo, we can now loop through each cell in the processor to get the corresponding gas pressure
            do m=1,ccl%struct(n)%n_
               ! Get cell indices
               i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
               ! Get cell position relative to drop barycenter, accounting for periodicity
               x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL
               y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL
               z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL
               ! Get the r magnitude
               xr=dot_product([x,y,z],dmoi(n,:,1)); yr=dot_product([x,y,z],dmoi(n,:,2))
               rmag=sqrt(xr**2+yr**2)
               do l =1, total_cell
                  if (gpinfo(l,1)>=rmag) then
                     rind=l
                     exit
                  end if 
               end do 
               
               ptmp=0.0_WP
               do l=rind,total_cell-1
                  ! 1/2 * (r(l+1)-r(l))*(f(l+1)+f(l))
                  ! f(l) = -12 mu_g * 1/deltaP  *(ur/h^2 + r dhdt/2h^3)          ! 1 is r, 2 is h, 3 is deltP, 4 is ur, 5 is uz
                  ! ptmp=ptmp-6.0_WP*fs%visc_g*((gpinfo(l,4)/(gpinfo(l,2)**2)+gpinfo(l,1)*gpinfo(l,5)/(2.0_WP*gpinfo(l,2)**3))/gpinfo(l,3)+&
                  ! & (gpinfo(l+1,4)/(gpinfo(l+1,2)**2)+gpinfo(l+1,1)*gpinfo(l+1,5)/(2.0_WP*gpinfo(l+1,2)**3))/gpinfo(l+1,3))*(gpinfo(l+1,1)-gpinfo(l,1))
                  ! f(l) = -12 mu_g * 1/deltaP  *(ur/h^2 + int_(r dhdt)/rh^3)    
                  ptmp=ptmp-6.0_WP*fs%visc_g*((gpinfo(l,4)/(gpinfo(l,2)**2)+rdhdt_int(l)/(gpinfo(l,1)*gpinfo(l,2)**3))/gpinfo(l,3)+&
                  & (gpinfo(l+1,4)/(gpinfo(l+1,2)**2)+rdhdt_int(l+1)/(gpinfo(l+1,1)*gpinfo(l+1,2)**3))/gpinfo(l+1,3))*(gpinfo(l+1,1)-gpinfo(l,1))
               end do
               Pg(i,j,k)=ptmp
               ! Pd(i,j,k)=-10**(HamakerC-log10(6*pi*thickness_new(i,j,k)**3))
               Pd(i,j,k)=-HamakerC/(6*pi*thickness_new(i,j,k)**3)
            end do
            deallocate(gpinfo,rdhdt_int)
         end if

         deallocate(gpinfo_)
      end do
      call cfg%sync(Pg)
      call cfg%sync(verticalU)
      call cfg%sync(radialU)
      contains
         !> Function that identifies cells that need a label
         logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         ! if (vf%VF(i,j,k).gt.0.0_WP) then
         if (vf%thin_sensor(i,j,k).eq.2.0_WP) then
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

   end subroutine get_gasP
   
   subroutine apply_gasP
      implicit none
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: Pgradx,Pgrady,Pgradz
      allocate(Pgradx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Pgrady(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Pgradz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               Pgradx(i,j,k)=sum(fs%divu_x(:,i,j,k)*Pg(i-1:i,j,k))
               Pgrady(i,j,k)=sum(fs%divv_y(:,i,j,k)*Pg(i,j-1:j,k))
               Pgradz(i,j,k)=sum(fs%divw_z(:,i,j,k)*Pg(i,j,k-1:k))
            end do
         end do
      end do
      ! Sync it
      call cfg%sync(Pgradx)
      call cfg%sync(Pgrady)
      call cfg%sync(Pgradz)
      resU=resU+Pgradx
      resV=resV+Pgrady
      resW=resW+Pgradz
   end subroutine apply_gasP


   subroutine get_thickness(thickness_in)
      use vfs_class, only: VFlo,VFhi
      implicit none 
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(out) :: thickness_in
      real(WP) :: tmplvol,tmpgvol,tmparea
      real(WP), dimension(1:3) :: tmpxvol, tmpL
      integer :: nneigh_thickness,i,j,k,ii,jj,kk
      nneigh_thickness=3
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ! calculate thickness
               tmplvol=0.0_WP; tmpgvol=0.0_WP;tmparea=0.0_WP
               do kk = k-nneigh_thickness,k+nneigh_thickness
                  do jj = j-nneigh_thickness,j+nneigh_thickness
                     do ii = i-nneigh_thickness,i+nneigh_thickness
                        ! tmpvol = tmpvol + vf%VF(ii,jj,kk)*cfg%vol(i,j,k)
                        ! tmparea = tmparea + vf%SD(ii,jj,kk)*cfg%vol(i,j,k)
                        tmplvol=tmplvol+(       vf%VF(ii,jj,kk))*cfg%vol(ii,jj,kk)
                        tmpgvol=tmpgvol+(1.0_WP-vf%VF(ii,jj,kk))*cfg%vol(ii,jj,kk)
                        tmparea=tmparea+        vf%SD(ii,jj,kk) *cfg%vol(ii,jj,kk)
                     end do
                  end do
               end do
               ! Calculate thickness
               if (vf%VF(i,j,k).ge.VFhi) then
                  thickness_in(i,j,k) = 0.0_WP
               else if (tmparea .gt. 0.0_WP) then    
                  thickness_in(i,j,k) = 2.0_WP*tmpgvol/(tmparea+tiny(1.0_WP))
               else
                  thickness_in(i,j,k) = 3.5_WP*cfg%min_meshsize
               end if
            end do 
         end do 
      end do
      call cfg%sync(thickness_in)
   end subroutine get_thickness

   ! A subroutine that solves the slip velocity field based on current info
   ! subroutine solveUs
   !    use vfs_class, only: VFlo,VFhi
   !    implicit none
   !    integer :: i,j,k,ii,jj,kk
   !    do k=cfg%kmino_  ,cfg%kmaxo_; do j=cfg%jmino_  ,cfg%jmaxo_; do i=cfg%imino_+1,cfg%imaxo_
   !       alpha_x(i,j,k)=sum(fs%itpr_x(:,i,j,k)*vf%VF(i-1:i,j,k))
   !    end do; end do; end do
   !    do k=cfg%kmino_  ,cfg%kmaxo_; do j=cfg%jmino_+1,cfg%jmaxo_; do i=cfg%imino_  ,cfg%imaxo_
   !       alpha_y(i,j,k)=sum(fs%itpr_y(:,i,j,k)*vf%VF(i,j-1:j,k))
   !    end do; end do; end do
   !    do k=cfg%kmino_+1,cfg%kmaxo_; do j=cfg%jmino_  ,cfg%jmaxo_; do i=cfg%imino_  ,cfg%imaxo_
   !       alpha_z(i,j,k)=sum(fs%itpr_z(:,i,j,k)*vf%VF(i,j,k-1:k))
   !    end do; end do; end do
   !    ! Handle non-periodic borders
   !    if (.not.cfg%xper.and.cfg%iproc.eq.1) alpha_x(cfg%imino,:,:)=vf%VF(cfg%imino,:,:)
   !    if (.not.cfg%yper.and.cfg%jproc.eq.1) alpha_y(:,cfg%jmino,:)=vf%VF(:,cfg%jmino,:)
   !    if (.not.cfg%zper.and.cfg%kproc.eq.1) alpha_z(:,:,cfg%kmino)=vf%VF(:,:,cfg%kmino)
   !    ! Synchronize boundaries
   !    call cfg%sync(alpha_x)
   !    call cfg%sync(alpha_y)
   !    call cfg%sync(alpha_z)
   !    ! Get a density in case I need to use it
   !    ! rho=vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
   !    get_dmomdt : block
   !       real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
   !       real(WP), dimension(:,:,:), allocatable :: drhoUdt,drhoVdt,drhoWdt
   !       ! Allocate flux arrays
   !       allocate(FX(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       allocate(FY(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       allocate(FZ(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       allocate(drhoUdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       allocate(drhoVdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       allocate(drhoWdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
   !       ! Zero out drhoUVW/dt arrays
   !       drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
   !       do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
   !          do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
   !             do ii=fs%cfg%imin_,fs%cfg%imax_+1
   !                ! Fluxes on x-face
   !                i=ii-1; j=jj-1; k=kk-1
   !                FX(i,j,k)=-sum(fs%itpu_x(:,i,j,k)*Us(i:i+1,j,k))*sum(fs%itpu_x(:,i,j,k)*fs%U(i:i+1,j,k)) &
   !                &         -sum(fs%itpu_x(:,i,j,k)*fs%U(i:i+1,j,k))*sum(fs%itpu_x(:,i,j,k)*Us(i:i+1,j,k)) 
   !                ! Fluxes on y-face
   !                i=ii; j=jj; k=kk
   !                FY(i,j,k)=-sum(fs%itpv_x(:,i,j,k)*Vs(i-1:i,j,k))*sum(fs%itpu_y(:,i,j,k)*fs%U(i,j-1:j,k)) &
   !                &         -sum(fs%itpv_x(:,i,j,k)*fs%V(i-1:i,j,k))*sum(fs%itpu_y(:,i,j,k)*Us(i,j-1:j,k))  
   !                ! Fluxes on z-face
   !                i=ii; j=jj; k=kk
   !                FZ(i,j,k)=-sum(fs%itpw_x(:,i,j,k)*Ws(i-1:i,j,k))*sum(fs%itpu_z(:,i,j,k)*fs%U(i,j,k-1:k)) &
   !                &         -sum(fs%itpw_x(:,i,j,k)*fs%W(i-1:i,j,k))*sum(fs%itpu_z(:,i,j,k)*Us(i,j,k-1:k)) 
   !             end do
   !          end do
   !       end do
   !       ! Time derivative of rhoU
   !       do k=fs%cfg%kmin_,fs%cfg%kmax_
   !          do j=fs%cfg%jmin_,fs%cfg%jmax_
   !             do i=fs%cfg%imin_,fs%cfg%imax_
   !                drhoUdt(i,j,k)=sum(fs%divu_x(:,i,j,k)*FX(i-1:i,j,k))+&
   !                &              sum(fs%divu_y(:,i,j,k)*FY(i,j:j+1,k))+&
   !                &              sum(fs%divu_z(:,i,j,k)*FZ(i,j,k:k+1))
   !                ! Makes sure the interpolated VF and (1-VF) are not 0
   !                if (alpha_x(i,j,k).gt.VFlo .and. alpha_x(i,j,k) .lt. VFhi) then
   !                   drhoUdt(i,j,k)=drhoUdt(i,j,k)+sum(fs%divu_x(:,i,j,k)*Pg(i-1:i,j,k)*(1.0_WP-vf%VF(i-1:i,j,k)))/(alpha_x(i,j,k)*fs%rho_g)
   !                else
   !                   drhoUdt(i,j,k)=-Usold(i,j,k)/time%dt
   !                end if
   !                   ! drhoUdt(i,j,k)= drhoUdt(i,j,k)-sum(fs%divu_x(:,i,j,k)*fs%P(i-1:i,j,k))/(sum(fs%itpr_x(:,i,j,k)*vf%VF(i-1:i,j,k))*fs%rho_l)+&
   !                   ! &              sum(fs%divu_x(:,i,j,k)*Pg(i-1:i,j,k)*(1.0_WP-vf%VF(i-1:i,j,k)))*&
   !                   ! &              sum(fs%itpr_x(:,i,j,k)*rho(i-1:i,j,k))/(fs%rho_l*fs%rho_g*sum(fs%itpr_x(:,i,j,k)*vf%VF(i-1:i,j,k))*sum(fs%itpr_x(:,i,j,k)*(1.0_WP-vf%VF(i-1:i,j,k))))
   !             end do
   !          end do
   !       end do
   !       ! Sync it
   !       call fs%cfg%sync(drhoUdt)
   !       ! Flux of rhoV
   !       do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
   !          do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
   !             do ii=fs%cfg%imin_,fs%cfg%imax_+1
   !                ! Fluxes on x-face
   !                i=ii; j=jj; k=kk
   !                FX(i,j,k)=-sum(fs%itpu_y(:,i,j,k)*Us(i,j-1:j,k))*sum(fs%itpv_x(:,i,j,k)*fs%V(i-1:i,j,k)) &
   !                &         -sum(fs%itpu_y(:,i,j,k)*fs%U(i,j-1:j,k))*sum(fs%itpv_x(:,i,j,k)*Vs(i-1:i,j,k)) 
   !                ! Fluxes on y-face
   !                i=ii-1; j=jj-1; k=kk-1
   !                FY(i,j,k)=-sum(fs%itpv_y(:,i,j,k)*Vs(i,j:j+1,k))*sum(fs%itpv_y(:,i,j,k)*fs%V(i,j:j+1,k)) &
   !                &         -sum(fs%itpv_y(:,i,j,k)*fs%V(i,j:j+1,k))*sum(fs%itpv_y(:,i,j,k)*Vs(i,j:j+1,k)) 
   !                ! Fluxes on z-face
   !                i=ii; j=jj; k=kk
   !                FZ(i,j,k)=-sum(fs%itpw_y(:,i,j,k)*Ws(i,j-1:j,k))*sum(fs%itpv_z(:,i,j,k)*fs%V(i,j,k-1:k)) &
   !                &         -sum(fs%itpw_y(:,i,j,k)*fs%W(i,j-1:j,k))*sum(fs%itpv_z(:,i,j,k)*Vs(i,j,k-1:k)) 
   !             end do
   !          end do
   !       end do
   !       ! Time derivative of rhoV
   !       do k=fs%cfg%kmin_,fs%cfg%kmax_
   !          do j=fs%cfg%jmin_,fs%cfg%jmax_
   !             do i=fs%cfg%imin_,fs%cfg%imax_
   !                drhoVdt(i,j,k)=sum(fs%divv_x(:,i,j,k)*FX(i:i+1,j,k))+&
   !                &              sum(fs%divv_y(:,i,j,k)*FY(i,j-1:j,k))+&
   !                &              sum(fs%divv_z(:,i,j,k)*FZ(i,j,k:k+1))
   !                ! Makes sure the interpolated VF and (1-VF) are not 0
   !                if (alpha_y(i,j,k).gt.VFlo .and. alpha_y(i,j,k).lt.VFhi) then
   !                   drhoVdt(i,j,k)=drhoVdt(i,j,k)+sum(fs%divv_y(:,i,j,k)*Pg(i,j-1:j,k)*(1.0_WP-vf%VF(i,j-1:j,k)))/(alpha_y(i,j,k)*fs%rho_g)
   !                else
   !                   drhoVdt(i,j,k)=-Vsold(i,j,k)/time%dt
   !                end if
   !                   ! drhoVdt(i,j,k)= drhoVdt(i,j,k)-sum(fs%divv_y(:,i,j,k)*fs%P(i,j-1:j,k))/(sum(fs%itpr_y(:,i,j,k)*vf%VF(i,j-1:j,k))*fs%rho_l)+&
   !                   ! &              sum(fs%divv_y(:,i,j,k)*Pg(i,j-1:j,k)*(1.0_WP-vf%VF(i,j-1:j,k)))*&
   !                   ! &              sum(fs%itpr_y(:,i,j,k)*rho(i,j-1:j,k))/(fs%rho_l*fs%rho_g*sum(fs%itpr_y(:,i,j,k)*vf%VF(i,j-1:j,k))*sum(fs%itpr_y(:,i,j,k)*(1.0_WP-vf%VF(i,j-1:j,k))))
   !             end do
   !          end do
   !       end do
   !       ! Sync it
   !       call fs%cfg%sync(drhoVdt)

   !       ! Flux of rhoW
   !       do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
   !          do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
   !             do ii=fs%cfg%imin_,fs%cfg%imax_+1
   !                ! Fluxes on x-face
   !                i=ii; j=jj; k=kk
   !                FX(i,j,k)=-sum(fs%itpu_z(:,i,j,k)*Us(i,j,k-1:k))*sum(fs%itpw_x(:,i,j,k)*fs%W(i-1:i,j,k)) &
   !                &         -sum(fs%itpu_z(:,i,j,k)*fs%U(i,j,k-1:k))*sum(fs%itpw_x(:,i,j,k)*Ws(i-1:i,j,k)) 
   !                ! Fluxes on y-face
   !                i=ii; j=jj; k=kk
   !                FY(i,j,k)=-sum(fs%itpv_z(:,i,j,k)*Vs(i,j,k-1:k))*sum(fs%itpw_y(:,i,j,k)*fs%W(i,j-1:j,k)) &
   !                &         -sum(fs%itpv_z(:,i,j,k)*fs%V(i,j,k-1:k))*sum(fs%itpw_y(:,i,j,k)*Ws(i,j-1:j,k)) 
   !                ! Fluxes on z-face
   !                i=ii-1; j=jj-1; k=kk-1
   !                FZ(i,j,k)=-sum(fs%itpw_z(:,i,j,k)*Ws(i,j,k:k+1))*sum(fs%itpw_z(:,i,j,k)*fs%W(i,j,k:k+1)) &
   !                &         -sum(fs%itpw_z(:,i,j,k)*fs%W(i,j,k:k+1))*sum(fs%itpw_z(:,i,j,k)*Ws(i,j,k:k+1)) 

   !             end do
   !          end do
   !       end do
   !       ! Time derivative of rhoW
   !       do k=fs%cfg%kmin_,fs%cfg%kmax_
   !          do j=fs%cfg%jmin_,fs%cfg%jmax_
   !             do i=fs%cfg%imin_,fs%cfg%imax_
   !                drhoWdt(i,j,k)=sum(fs%divw_x(:,i,j,k)*FX(i:i+1,j,k))+&
   !                &              sum(fs%divw_y(:,i,j,k)*FY(i,j:j+1,k))+&
   !                &              sum(fs%divw_z(:,i,j,k)*FZ(i,j,k-1:k))
   !                ! Makes sure the interpolated VF and (1-VF) are not 0
   !                if (alpha_z(i,j,k).gt.VFlo .and. alpha_z(i,j,k).lt.VFhi) then
   !                   drhoWdt(i,j,k)=drhoWdt(i,j,k)+sum(fs%divw_z(:,i,j,k)*Pg(i,j,k-1:k)*(1.0_WP-vf%VF(i,j,k-1:k)))/(alpha_z(i,j,k)*fs%rho_g)
   !                else
   !                   drhoWdt(i,j,k)=-Wsold(i,j,k)/time%dt
   !                end if
   !                   ! drhoWdt(i,j,k)= drhoWdt(i,j,k)-sum(fs%divw_z(:,i,j,k)*fs%P(i,j,k-1:k))/(sum(fs%itpr_z(:,i,j,k)*vf%VF(i,j,k-1:k))*fs%rho_l)+&
   !                   ! &              sum(fs%divw_z(:,i,j,k)*Pg(i,j,k-1:k)*(1.0_WP-vf%VF(i,j,k-1:k)))*&
   !                   ! &              sum(fs%itpr_z(:,i,j,k)*rho(i,j,k-1:k))/(fs%rho_l*fs%rho_g*sum(fs%itpr_z(:,i,j,k)*vf%VF(i,j,k-1:k))*sum(fs%itpr_z(:,i,j,k)*(1.0_WP-vf%VF(i,j,k-1:k))))
   !             end do
   !          end do
   !       end do
   !       ! Sync it
   !       call fs%cfg%sync(drhoWdt)
   !       resU=drhoUdt; resV=drhoVdt; resW=drhoWdt
   !    end block get_dmomdt

   !    ! Only focused on current interfacial cells
   !    Us=Usold+resU*time%dt
   !    Vs=Vsold+resV*time%dt
   !    Ws=Wsold+resW*time%dt

   !    call fs%update_laplacian_slip(alpha_x,alpha_y,alpha_z,VFlo)
   !    call fs%get_div_slip(Us,Vs,Ws)
   !    fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
   !    fs%psolv%sol=0.0_WP
   !    call fs%psolv%solve()
   !    call fs%shift_p(fs%psolv%sol)
   !    ! Correct velocity
   !    call fs%get_pgrad_slip(fs%psolv%sol,resU,resV,resW)

   !    Us=Us-time%dt*resU
   !    Vs=Vs-time%dt*resV
   !    Ws=Ws-time%dt*resW
   !    ! do k=cfg%kmino_,cfg%kmaxo_
   !    !    do j=cfg%jmino_,cfg%jmaxo_
   !    !       do i=cfg%imino_, cfg%imaxo_
   !    !          if (alpha_x(i,j,k).eq. 0.0_WP)  then
   !    !             Us(i,j,k)=0.0_WP
   !    !          else
   !    !             Us(i,j,k)=Us(i,j,k)-time%dt*resU(i,j,k)/(alpha_x(i,j,k)*fs%rho_l)
   !    !          end if

   !    !          if (alpha_y(i,j,k).eq.0.0_WP) then
   !    !             Vs(i,j,k)=0.0_WP
   !    !          else
   !    !             Vs(i,j,k)=Vs(i,j,k)-time%dt*resV(i,j,k)/(alpha_y(i,j,k)*fs%rho_l)   
   !    !          end if

   !    !          if (alpha_z(i,j,k).eq.0.0_WP) then
   !    !             Ws(i,j,k)=0.0_WP
   !    !          else
   !    !             Ws(i,j,k)=Ws(i,j,k)-time%dt*resW(i,j,k)/(alpha_z(i,j,k)*fs%rho_l)
   !    !          end if
   !    !       end do 
   !    !    end do
   !    ! end do
   !    ! Us=0.0_WP 
   !    ! Vs=0.0_WP 
   !    ! Ws=0.0_WP 

   !    Usold=Us
   !    Vsold=Vs
   !    Wsold=Ws
   ! end subroutine solveUs
   !> Function that defines a level set function for colliding drops problem
   function levelset_colliding_drops(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,G1,G2
      ! Create droplet 1
      G1=radius1-sqrt(sum((xyz-center1)**2))
      ! Create droplet 1
      G2=radius2-sqrt(sum((xyz-center2)**2))
      ! Combine
      G=max(G1,G2)
   end function levelset_colliding_drops
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Ui=0.0_WP
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vi=0.0_WP
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wi=0.0_WP
         ! allocate(Us  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Us=0.0_WP 
         ! allocate(Vs  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vs=0.0_WP 
         ! allocate(Ws  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Ws=0.0_WP 
         ! allocate(Usold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Usold=0.0_WP 
         ! allocate(Vsold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vsold=0.0_WP 
         ! allocate(Wsold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wsold=0.0_WP 
         ! allocate(rho  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Pg  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Pg=0.0_WP
         allocate(Pd  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Pd=0.0_WP
         allocate(radialU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));radialU=0.0_WP
         allocate(verticalU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));verticalU=0.0_WP

         allocate(thickness_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));thickness_old=0.0_WP
         allocate(thickness_new(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));thickness_new=0.0_WP
         ! allocate(alpha_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(alpha_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(alpha_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: plicnet,elvira,VFhi,VFlo,r2p,r2pnet
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with r2p reconstruction
         ! call vf%initialize(cfg=cfg,reconstruction_method=r2pnet,name='VOF')
         ! call vf%initialize(cfg=cfg,reconstruction_method=r2p,name='VOF')
         call vf%initialize(cfg=cfg,reconstruction_method=r2pnet,name='VOF')
         vf%thin_thld_max=1.5_WP
         vf%twoplane_thld2=0.8_WP
         ! Initialize two droplets
         call param_read('Droplet 1 diameter',radius1); radius1=0.5_WP*radius1
         call param_read('Droplet 1 position',center1)
         call param_read('Droplet 2 diameter',radius2); radius2=0.5_WP*radius2
         call param_read('Droplet 2 position',center2)
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[cfg%x(i+si),cfg%y(j+sj),cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_colliding_drops,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([cfg%xm(i),cfg%ym(j),cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[cfg%xm(i),cfg%ym(j),cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[cfg%xm(i),cfg%ym(j),cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
         ! Initialize a density field
         ! rho=vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
         ! rho=1.0_WP
      end block create_and_initialize_vof
      
      initialize_ccl: block
         call ccl%initialize(pg=cfg%pgrid,name='ccl')
      end block initialize_ccl
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=20
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         !vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)!,implicit_solver=vs)
         ! Initial droplet velocity
         call param_read('Droplet 1 velocity',vel1)
         call param_read('Droplet 2 velocity',vel2)
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  ! U velocity
                  xyz=[fs%cfg%x(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%U(i,j,k)=vel1(1)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%U(i,j,k)=vel2(1)
                  ! V velocity
                  xyz=[fs%cfg%xm(i),fs%cfg%y(j),fs%cfg%zm(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%V(i,j,k)=vel1(2)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%V(i,j,k)=vel2(2)
                  ! V velocity
                  xyz=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%z(k)]
                  if (radius1-sqrt(sum((xyz-center1)**2)).ge.0.0_WP) fs%W(i,j,k)=vel1(3)
                  if (radius2-sqrt(sum((xyz-center2)**2)).ge.0.0_WP) fs%W(i,j,k)=vel2(3)
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Set slip vel to 0
         ! Us=0.0_WP; Vs=0.0_WP; Ws=0.0_WP
         ! Usold=0.0_WP; Vsold=0.0_WP; Wsold=0.0_WP
      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=7,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='radialU'
         smesh%varname(4)='verticalU'
         smesh%varname(5)='gasP'
         smesh%varname(6)='thickness_old'
         smesh%varname(7)='thickness_new'
         ! smesh%varname(3)='edge_sensor'
         ! smesh%varname(4)='thin_sensor'
         ! smesh%varname(5)='thickness'
         ! smesh%varname(6)='norm_abs'
         ! smesh%varname(7)='norm_sig'
         ! smesh%varname(8)='ccl'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                        smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                        smesh%var(3,np)=radialU(i,j,k)
                        smesh%var(4,np)=verticalU(i,j,k)
                        smesh%var(5,np)=Pg(i,j,k)
                        smesh%var(6,np)=thickness_old(i,j,k)
                        smesh%var(7,np)=thickness_new(i,j,k)
                        ! smesh%var(3,np)=vf%edge_sensor(i,j,k)
                        ! smesh%var(4,np)=vf%thin_sensor(i,j,k)
                        ! smesh%var(5,np)=vf%thickness  (i,j,k)
                        ! smesh%var(6,np)=vf%norm_pos(i,j,k)+vf%norm_neg(i,j,k)
                        ! smesh%var(7,np)=vf%norm_pos(i,j,k)-vf%norm_neg(i,j,k)
                        ! smesh%var(8,np)=real(ccl%id(i,j,k),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='CollidingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         ! call ens_out%add_vector('slipvel',resU,resV,resW)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('Gpressure',Pg)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('radialU',radialU)
         call ens_out%add_scalar('verticalU',verticalU)
         call ens_out%add_surface('vofplic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(vf%flotsam_error,'Flotsam error')
         call mfile%add_column(vf%thinstruct_error,'Film error')
         call mfile%add_column(vf%SDint,'SD integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - fs mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc,harmonic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! thickness_old=vf%thickness
         call get_thickness(thickness_old)
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! Update the velocity to be the corrected liquid velocity
         ! resU=fs%U+(Us*(1.0_WP-alpha_x)*fs%rho_g)/(alpha_x*fs%rho_l+(1.0_WP-alpha_x)*fs%rho_g)
         ! resV=fs%V+(Vs*(1.0_WP-alpha_y)*fs%rho_g)/(alpha_y*fs%rho_l+(1.0_WP-alpha_y)*fs%rho_g)
         ! resW=fs%W+(Ws*(1.0_WP-alpha_z)*fs%rho_g)/(alpha_z*fs%rho_l+(1.0_WP-alpha_z)*fs%rho_g)
         ! resU=fs%U-(Us*alpha_x*fs%rho_l)/(alpha_x*fs%rho_l+(1.0_WP-alpha_x)*fs%rho_g)
         ! resV=fs%V-(Vs*alpha_y*fs%rho_l)/(alpha_y*fs%rho_l+(1.0_WP-alpha_y)*fs%rho_g)
         ! resW=fs%W-(Ws*alpha_z*fs%rho_l)/(alpha_z*fs%rho_l+(1.0_WP-alpha_z)*fs%rho_g)
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         ! call vf%advance(dt=time%dt,U=resU,V=resV,W=resW)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
         call get_thickness(thickness_new)
         ! Get the gas pressure
         call get_gasP()

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! call apply_gasP()
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            ! call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho_W
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            call fs%add_surface_tension_jump_thin(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         ! Solve the slip velocity
         ! call solveUs()
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! !> Calculate the interpolated velocity, including overlap and ghosts
         ! interp_vel: block         
         !    integer :: i,j,k
         !    ! Calculate as far as possible each component
         !    do k=cfg%kmino_,cfg%kmaxo_
         !       do j=cfg%jmino_,cfg%jmaxo_
         !          do i=cfg%imino_,cfg%imaxo_-1
         !             resU(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Us(i:i+1,j,k))
         !          end do
         !       end do
         !    end do
         !    do k=cfg%kmino_,cfg%kmaxo_
         !       do j=cfg%jmino_,cfg%jmaxo_-1
         !          do i=cfg%imino_,cfg%imaxo_
         !             resV(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vs(i,j:j+1,k))
         !          end do
         !       end do
         !    end do
         !    do k=cfg%kmino_,cfg%kmaxo_-1
         !       do j=cfg%jmino_,cfg%jmaxo_
         !          do i=cfg%imino_,cfg%imaxo_
         !             resW(i,j,k)=sum(fs%itpw_z(:,i,j,k)*Ws(i,j,k:k+1))
         !          end do
         !       end do
         !    end do
         !    ! Add last layer in each direction
         !    if (.not.cfg%xper.and.cfg%iproc.eq.cfg%npx) resU(cfg%imaxo,:,:)=Us(cfg%imaxo,:,:)
         !    if (.not.cfg%yper.and.cfg%jproc.eq.cfg%npy) resV(:,cfg%jmaxo,:)=Vs(:,cfg%jmaxo,:)
         !    if (.not.cfg%zper.and.cfg%kproc.eq.cfg%npz) resW(:,:,cfg%kmaxo)=Ws(:,:,cfg%kmaxo)
         !    ! Sync it
         !    call cfg%sync(Us)
         !    call cfg%sync(Vs)
         !    call cfg%sync(Ws)
         ! end block interp_vel

         ! Output to ensight
         if (ens_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Also populate nplane variable
               smesh%var(1,:)=1.0_WP
               np=0
               do k=cfg%kmin_,cfg%kmax_
                  do j=cfg%jmin_,cfg%jmax_
                     do i=cfg%imin_,cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                              smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                              smesh%var(3,np)=radialU(i,j,k)
                              smesh%var(4,np)=verticalU(i,j,k)
                              smesh%var(5,np)=Pg(i,j,k)
                              smesh%var(6,np)=thickness_old(i,j,k)
                              smesh%var(7,np)=thickness_new(i,j,k)
                              ! smesh%var(3,np)=vf%edge_sensor(i,j,k)
                              ! smesh%var(4,np)=vf%thin_sensor(i,j,k)
                              ! smesh%var(5,np)=vf%thickness  (i,j,k)
                              ! smesh%var(6,np)=vf%norm_pos(i,j,k)+vf%norm_neg(i,j,k)
                              ! smesh%var(7,np)=vf%norm_pos(i,j,k)-vf%norm_neg(i,j,k)
                              ! smesh%var(8,np)=real(ccl%id(i,j,k),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            ! Perform ensight output
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
