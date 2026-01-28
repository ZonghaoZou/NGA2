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
   use lpt_class,         only: lpt
   use partmesh_class,    only: partmesh
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
   
   type(cclabel)  :: ccl
   type(lpt)      :: lp         !< Lagrangian particle for estimating flattented radius
   type(partmesh) :: pmesh      !< Particle mesh for showing the particle

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final,get_gasVel,get_thickness,record_thickness
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Us,Vs,Ws
   real(WP), dimension(:,:,:), allocatable :: thickness_old,thickness_new
   integer,  dimension(:,:,:), allocatable :: region_indicator,mask_IB
   
   !> Problem definition
   real(WP), dimension(3) :: center1,center2,vel1,vel2
   real(WP) :: radius1,radius2
   real(WP) :: HamakerC,lambdaAir
   real(WP) :: anew,aold,amax,minThickness,thickthd1,thickthd2,mulfact,thicknew,minThicknessold
   real(WP) :: meshaligned_angle,uzvel
   integer :: curr_case,mesh_ratio
   logical :: activated,breakup,default
   contains
! This is based on Zhang and Law's theoretical gas pressure derivation
subroutine get_gasVel
   use irl_fortran_interface
   use mpi_f08
   use parallel,  only: MPI_REAL_WP
   use mathtools, only: pi,normalize
   use messager,  only: die
   use vfs_class, only: VFlo,VFhi
   use string,    only: str_medium
   implicit none
   ! Parameters for moment of inertia
   real(WP), dimension(:), allocatable, save :: work !< Saved!
   integer, save :: lwork                            !< Saved!
   real(WP), dimension(1) :: lwork_query
   real(WP), dimension(3) :: d
   real(WP), dimension(3,3) :: A
   integer :: info
   integer :: ierr,i,j,k,m,n,iunit
   real(WP) :: x,y,z,myvol,a_cell,Kn,xr,yr,zr
   real(WP) :: x0,y0,z0,dhdt,thickold,DeltaKn,dlogadt
   character(len=str_medium) :: filename
   real(WP), dimension(3) :: mybary,t1,t2,t3
   real(WP), dimension(:)    , allocatable :: fgvol
   real(WP), dimension(:)    , allocatable :: fct
   real(WP), dimension(:,:)  , allocatable :: fgpos
   real(WP), dimension(:,:,:), allocatable :: fmoi
   
   ! Set indicator to 0
   anew=0.0_WP;x0=0.0_WP;y0=0.0_WP;z0=0.0_WP;a_cell=0.0_WP
   mask_IB=1;thickold=0.0_WP;thicknew=0.0_WP;dhdt=0.0_WP
   vf%indicator=0
   ! Query optimal work array size
   if (.not.allocated(work)) then
      call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
      lwork=int(lwork_query(1)); allocate(work(lwork))
   end if
   ! Build ccl to get the thin gas region for caculating gas pressure
   call ccl%build(make_label,same_label)

   ! Allocate fields for calculation
   allocate(fgvol(1:ccl%nstruct        )); fgvol=0.0_WP
   allocate(fct  (1:ccl%nstruct        )); fct=0.0_WP
   allocate(fgpos(1:ccl%nstruct,1:3    )); fgpos=0.0_WP
   allocate(fmoi (1:ccl%nstruct,1:3,1:3)); fmoi =0.0_WP

   ! First pass to accumulate position for moment of inertia
   do n=1,ccl%nstruct
      ! Loop over cells in structure
      do m=1,ccl%struct(n)%n_
         ! Get cell indices
         i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
         x=vf%cfg%xm(i)-ccl%struct(n)%per(1)*cfg%xL!x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL! x=vf%cfg%xm(i)-ccl%struct(n)%per(1)*cfg%xL
         y=vf%cfg%ym(j)-ccl%struct(n)%per(2)*cfg%yL!y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL! y=vf%cfg%ym(j)-ccl%struct(n)%per(2)*cfg%yL
         z=vf%cfg%zm(k)-ccl%struct(n)%per(3)*cfg%zL!z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL! z=vf%cfg%zm(k)-ccl%struct(n)%per(3)*cfg%zL
         ! Accumulate volume and position
         if (thickness_new(i,j,k).le.0.7_WP*cfg%min_meshsize) then
            fgvol(n  )=fgvol(n  )+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))
            fgpos(n,:)=fgpos(n,:)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*[x,y,z]
            thickold=thickold+thickness_old(i,j,k)
            thicknew=thicknew+thickness_new(i,j,k)
            fct(n)=fct(n)+1.0_WP
            mask_IB(i,j,k)=0
            vf%indicator(i,j,k)=1
         end if 
     end do 
   end do
   call cfg%sync(mask_IB)
   call MPI_ALLREDUCE(MPI_IN_PLACE,fgvol   ,1*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,fgpos   ,3*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,fct     ,1*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,thickold,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,thicknew,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   ! Second pass to accumulate moment of inertia
   do n=1,ccl%nstruct
      if (fct(n).eq.0.0_WP) cycle
      ! Get the region gas barycenter
      x0=fgpos(n,1)/fgvol(n)
      y0=fgpos(n,2)/fgvol(n)
      z0=fgpos(n,3)/fgvol(n)
      thicknew=thicknew/fct(n)
      thickold=thickold/fct(n)
      dhdt=(thicknew-thickold)/time%dt
      ! Loop over cells in structure
      do m=1,ccl%struct(n)%n_
          ! Get cell indices
          i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
          x=vf%cfg%xm(i)-ccl%struct(n)%per(1)*cfg%xL-x0!x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL-x0!  x=vf%cfg%xm(i)-ccl%struct(n)%per(1)*cfg%xL-x0
          y=vf%cfg%ym(j)-ccl%struct(n)%per(2)*cfg%yL-y0!y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL-y0!  y=vf%cfg%ym(j)-ccl%struct(n)%per(2)*cfg%yL-y0
          z=vf%cfg%zm(k)-ccl%struct(n)%per(3)*cfg%zL-z0!z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL-z0!  z=vf%cfg%zm(k)-ccl%struct(n)%per(3)*cfg%zL-z0
          ! Accumulate moment of inertia
          fmoi(n,2,2)=fmoi(n,2,2)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(z**2+x**2)
          fmoi(n,3,3)=fmoi(n,3,3)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x**2+y**2)
          fmoi(n,1,1)=fmoi(n,1,1)+cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(y**2+z**2)
          fmoi(n,1,2)=fmoi(n,1,2)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x*y)
          fmoi(n,1,3)=fmoi(n,1,3)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(x*z)
          fmoi(n,2,3)=fmoi(n,2,3)-cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(y*z)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,fmoi,9*ccl%nstruct,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   ! Get all the moment of inertia
   do n=1,ccl%nstruct
      if (fct(n).eq.0.0_WP) cycle
      x0=fgpos(n,1)/fgvol(n)
      y0=fgpos(n,2)/fgvol(n)
      z0=fgpos(n,3)/fgvol(n)
      ! Get the moi directions
      A=fmoi(n,:,:)
      call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
      fmoi(n,:,:)=A ! fmoi(n,:,1) and fmoi(n,:,2) are the two principle axes marking the tagential plane
      t1=fmoi(n,:,1); t2=fmoi(n,:,2); t3=fmoi(n,:,3)
      do m=1,ccl%struct(n)%n_
         ! Get cell indices
         i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
         ! Get cell position relative to drop barycenter, accounting for periodicity
         x=vf%cfg%xm(i)-ccl%struct(n)%per(1)*cfg%xL-x0!x=vf%Gbary(1,i,j,k)-ccl%struct(n)%per(1)*cfg%xL-x0!
         y=vf%cfg%ym(j)-ccl%struct(n)%per(2)*cfg%yL-y0!y=vf%Gbary(2,i,j,k)-ccl%struct(n)%per(2)*cfg%yL-y0!
         z=vf%cfg%zm(k)-ccl%struct(n)%per(3)*cfg%zL-z0!z=vf%Gbary(3,i,j,k)-ccl%struct(n)%per(3)*cfg%zL-z0!
         xr=dot_product([x,y,z],fmoi(n,:,1)); yr=dot_product([x,y,z],fmoi(n,:,2))
         a_cell=max(sqrt(xr**2+yr**2),a_cell)
      end do
   end do
   call MPI_ALLREDUCE(MPI_IN_PLACE,a_cell,1*ccl%nstruct,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
   do n=1,ccl%nstruct
      if (fct(n).eq.0.0_WP) cycle
      ! If this is first time detected a thin gas region use the cell based estimation
      if (.not. activated .and. ccl%nstruct.ge.1)  then
         ! Estimation of the flattened radius
         activated=.true.
         anew=a_cell
         aold=anew
         if (cfg%amRoot) then
            lp%np_=lp%np_+1
            call lp%resize(lp%np_)
            lp%p(lp%np_)%id=int(1,8)
            lp%p(lp%np_)%d=1.0e-7_WP 
            lp%p(lp%np_)%pos=anew*fmoi(n,:,1)+[x0,y0,z0]
            lp%p(lp%np_)%vel =0.0_WP
            lp%p(lp%np_)%ind =cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
            lp%p(lp%np_)%flag=0                                                                                        
            lp%p(lp%np_)%dt  =0.0_WP                                                                                  
            lp%p(lp%np_)%Acol=0.0_WP                                                                                  
            lp%p(lp%np_)%Tcol=0.0_WP
         end if
         ! Move the particle to the correct processor
         call lp%sync()
      else if (activated) then
         ! If already activated store the old value
         ! Advance the particle based on its local velocity and get the new radius 
         do i = 1, lp%np_
            lp%p(lp%np_)%vel = cfg%get_velocity(pos=lp%p(lp%np_)%pos,i0=lp%p(lp%np_)%ind(1),j0=lp%p(lp%np_)%ind(2),k0=lp%p(lp%np_)%ind(3),U=fs%U,V=fs%V,W=fs%W)
            lp%p(lp%np_)%pos = lp%p(lp%np_)%pos + time%dt*lp%p(lp%np_)%vel
            lp%p(lp%np_)%ind =cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
            xr=dot_product(lp%p(lp%np_)%pos,fmoi(n,:,1)); yr=dot_product(lp%p(lp%np_)%pos,fmoi(n,:,2))
            anew=sqrt(xr**2+yr**2)
            if (abs(anew-a_cell).gt.0.1_WP*a_cell) then
               anew=a_cell
               lp%p(lp%np_)%pos = anew*fmoi(n,:,1)+[x0,y0,z0]
               lp%p(lp%np_)%ind =cfg%get_ijk_global(lp%p(lp%np_)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])     
            end if
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,anew,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         ! Move the particle to the correct processor
         call lp%sync()
      end if
   end do
   ! Get global parameters
   ! Kn=lambdaAir/thicknew
   thicknew=minThickness
   thickold=minThicknessold
   dlogadt=(log(anew)-log(aold))/time%dt
   dhdt=(thicknew-thickold)/time%dt
   aold=anew
   ! thicknew=minThickness
   Kn=lambdaAir/thicknew
   if (Kn .ge.1) then
      DeltaKn = 8.7583_WP*Kn**1.1551
   else
      DeltaKn = 1.0_WP + 6.0966_WP*Kn + 0.965_WP*Kn**2 +0.6967_WP*Kn**3
   end if
   dlogadt=0.0_WP
   ! if (cfg%amRoot) print *, "t1: ", t1, "t2: ",t2,"t3: ",t3
   ! DeltaKn=1.0_WP
   ! dlogadt=0.0_WP
   ! if (cfg%amRoot) then
   !    ! write(filename,'("thickness_p",I0,".csv")') 
   !    filename="output.txt"
   !    open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',position='append',iostat=ierr)
   !    write(iunit,*) time%t,"dlogadt", dlogadt, "dhdt",dhdt,"hnew",thicknew,"hold",thickold,"anew", anew, "aold", aold, "anewm",anew/cfg%min_meshsize,"aoldm", aold/cfg%min_meshsize,"dt", time%dt
   !    close(iunit)
   ! end if
   ! if (cfg%amRoot) print *, "dlogadt", dlogadt, "dhdt",dhdt,"hnew",thicknew,"hold",thickold,"anew", anew, "aold", aold, "anewm",anew/cfg%min_meshsize,"aoldm", aold/cfg%min_meshsize,"dt", time%dt
   ! We can now directly prescribe my gas lubrication velocity in r and z directions
   getUs: block
      real(WP) :: Ur,Uz,rmag
      real(WP), dimension(3) :: er,ef
      integer, dimension(:,:,:), allocatable :: mask_IB_x,mask_IB_y,mask_IB_z,tmp_mask
      allocate(tmp_mask(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_));tmp_mask=mask_IB
      ! Extend mask_IB for one more layer
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (mask_IB(i,j,k).eq.0) then
                  tmp_mask(i-1,j,k)=0; tmp_mask(i+1,j,k)=0
                  tmp_mask(i,j-1,k)=0; tmp_mask(i,j+1,k)=0
                  tmp_mask(i,j,k-1)=0; tmp_mask(i,j,k+1)=0
               end if
            end do 
         end do 
      end do 
      call cfg%sync(tmp_mask); mask_IB=tmp_mask
      allocate(mask_IB_x(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_));mask_IB_x=1
      allocate(mask_IB_y(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_));mask_IB_y=1
      allocate(mask_IB_z(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_));mask_IB_z=1
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (sum(mask_IB(i-1:i,j,k)).lt.2) mask_IB_x(i,j,k)=0
               if (sum(mask_IB(i,j-1:j,k)).lt.2) mask_IB_y(i,j,k)=0
               if (sum(mask_IB(i,j,k-1:k)).lt.2) mask_IB_z(i,j,k)=0
            end do 
         end do 
      end do 
      ! Handle non-periodic borders
      if (.not.cfg%xper.and.cfg%iproc.eq.1) mask_IB_x(cfg%imino,:,:)=mask_IB(cfg%imino,:,:)
      if (.not.cfg%yper.and.cfg%jproc.eq.1) mask_IB_y(:,cfg%jmino,:)=mask_IB(:,cfg%jmino,:)
      if (.not.cfg%zper.and.cfg%kproc.eq.1) mask_IB_z(:,:,cfg%kmino)=mask_IB(:,:,cfg%kmino)
      call cfg%sync(mask_IB_x)
      call cfg%sync(mask_IB_y)
      call cfg%sync(mask_IB_z)

      uzvel=mulfact*(2*(5.0_WP-4.0_WP*DeltaKn)*dlogadt*thicknew+5.0_WP*dhdt)/(16*DeltaKn)
      Us=0.0_WP;Vs=0.0_WP;Ws=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (mask_IB_x(i,j,k).eq.0) then
                  x=cfg%x(i)-x0; y=cfg%ym(j)-y0; z=cfg%zm(k)-z0
                  xr=dot_product([x,y,z],t1); yr=dot_product([x,y,z],t2); zr=dot_product([x,y,z],t3)
                  rmag=sqrt(xr**2+yr**2); er=normalize(xr*t1+yr*t2); ef=[1.0_WP,0.0_WP,0.0_WP]
                  ! Ur=2.0_WP*rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Ur=rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Uz=mulfact*(2*(5.0_WP-4.0_WP*DeltaKn)*dlogadt*thicknew+5.0_WP*dhdt)/(16*DeltaKn)
                  ! Ur=mulfact*-rmag*(dlogadt*(DeltaKn-1.0_WP)+DeltaKn*dhdt/(2.0_WP*thicknew))
                  ! Uz=mulfact*(2.0_WP*dlogadt*thicknew*(5*DeltaKn-4.0_WP)+5*DeltaKn*dhdt)/16.0_WP
                  ! if (vf%VF(i-1,j,k).gt.VFlo.and.vf%VF(i-1,j,k).lt.VFhi .and. vf%VF(i,j,k).gt.VFlo.and.vf%VF(i,j,k).lt.VFhi) then
                  ! if (abs(zr))
                  !    ! Could add another constraint on seeing if the face is contained in gas, but skip that for now
                     ! Us(i,j,k)=Ur*dot_product(er,ef)
                  ! else
                     Ur=0.0_WP
                     Us(i,j,k)=Ur*dot_product(er,ef)+sign(1.0_WP,zr)*dot_product(t3,ef)*Uz!-fs%U(i,j,k)
                  ! end if
                  ! Us(i,j,k)=sign(1.0_WP,dot_product([x,y,z],t3))*dot_product(t3,ef)*Uz
               end if

               if (mask_IB_y(i,j,k).eq.0) then
                  x=cfg%xm(i)-x0; y=cfg%y(j)-y0; z=cfg%zm(k)-z0
                  xr=dot_product([x,y,z],t1); yr=dot_product([x,y,z],t2); zr=dot_product([x,y,z],t3)
                  rmag=sqrt(xr**2+yr**2); er=normalize(xr*t1+yr*t2); ef=[0.0_WP,1.0_WP,0.0_WP]
                  ! Ur=2.0_WP*rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Ur=rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Uz=mulfact*(2*(5.0_WP-4.0_WP*DeltaKn)*dlogadt*thicknew+5.0_WP*dhdt)/(16*DeltaKn)
                  ! if (vf%VF(i,j-1,k).gt.VFlo.and.vf%VF(i,j-1,k).lt.VFhi .and. vf%VF(i,j,k).gt.VFlo.and.vf%VF(i,j,k).lt.VFhi) then
                  ! !    ! Could add another constraint on seeing if the face is contained in gas, but skip that for now
                     ! Vs(i,j,k)=Ur*dot_product(er,ef)
                  ! else
                  Ur=0.0_WP
                     Vs(i,j,k)=Ur*dot_product(er,ef)+sign(1.0_WP,zr)*dot_product(t3,ef)*Uz!-fs%V(i,j,k)
                  ! end if
                  ! Vs(i,j,k)=sign(1.0_WP,dot_product([x,y,z],t3))*dot_product(t3,ef)*Uz
               end if

               if (mask_IB_z(i,j,k).eq.0) then
                  x=cfg%xm(i)-x0; y=cfg%ym(j)-y0; z=cfg%z(k)-z0
                  xr=dot_product([x,y,z],t1); yr=dot_product([x,y,z],t2); zr=dot_product([x,y,z],t3)
                  rmag=sqrt(xr**2+yr**2); er=normalize(xr*t1+yr*t2); ef=[0.0_WP,0.0_WP,1.0_WP]
                  ! Ur=2.0_WP*rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Ur=rmag*(2.0_WP*(DeltaKn-1.0_WP)*dlogadt*thicknew-dhdt)/(2.0_WP*DeltaKn*thicknew)
                  Uz=mulfact*(2*(5.0_WP-4.0_WP*DeltaKn)*dlogadt*thicknew+5.0_WP*dhdt)/(16*DeltaKn)
                  ! if (vf%VF(i,j,k-1).gt.VFlo.and.vf%VF(i,j,k-1).lt.VFhi .and. vf%VF(i,j,k).gt.VFlo.and.vf%VF(i,j,k).lt.VFhi) then
                  ! !    ! Could add another constraint on seeing if the face is contained in gas, but skip that for now
                     ! Ws(i,j,k)=Ur*dot_product(er,ef)
                  ! else
                  Ur=0.0_WP
                     Ws(i,j,k)=Ur*dot_product(er,ef)+sign(1.0_WP,zr)*dot_product(t3,ef)*Uz!-fs%W(i,j,k)
                  ! end if
                  ! Ws(i,j,k)=sign(1.0_WP,dot_product([x,y,z],t3))*dot_product(t3,ef)*Uz
               end if
            end do 
         end do 
      end do 
      deallocate(mask_IB_x,mask_IB_y,mask_IB_z)
   end block getUs
   call cfg%sync(Us)
   call cfg%sync(Vs)
   call cfg%sync(Ws)
   contains
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
      implicit none
      integer, intent(in) :: i,j,k
      ! if (vf%VF(i,j,k).gt.0.0_WP) then
      if (vf%thin_sensor(i,j,k).eq.2.0_WP) then! .and.thickness_new(i,j,k).le.1.1*cfg%min_meshsize) then
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

end subroutine get_gasVel

   subroutine get_thickness(thickness_in)
      use mpi_f08
      use vfs_class, only: VFlo,VFhi
      use parallel,  only: MPI_REAL_WP
      implicit none 
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(out) :: thickness_in
      real(WP) :: tmplvol,tmpgvol,tmparea,area
      real(WP), dimension(1:3) :: tmpxvol, tmpL
      real(WP), dimension(:,:,:), allocatable :: tmp
      integer :: nneigh_thickness,i,j,k,ii,jj,kk,ierr
      nneigh_thickness=2
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ! calculate thickness
               tmplvol=0.0_WP; tmpgvol=0.0_WP;tmparea=0.0_WP
               do kk = k-nneigh_thickness,k+nneigh_thickness
                  do jj = j-nneigh_thickness,j+nneigh_thickness
                     do ii = i-nneigh_thickness,i+nneigh_thickness
                        tmplvol=tmplvol+(       vf%VF(ii,jj,kk))*cfg%vol(ii,jj,kk)
                        tmpgvol=tmpgvol+(1.0_WP-vf%VF(ii,jj,kk))*cfg%vol(ii,jj,kk)
                        tmparea=tmparea+        vf%SD(ii,jj,kk) *cfg%vol(ii,jj,kk)
                     end do
                  end do
               end do
               ! Calculate thickness
               if (vf%VF(i,j,k).ge.VFhi) then
                  thickness_in(i,j,k) = 3.5_WP*cfg%min_meshsize
               else if (tmparea .gt. 0.0_WP) then    
                  thickness_in(i,j,k) = min(2.0_WP*tmpgvol/(tmparea+tiny(1.0_WP)),3.5_WP*cfg%min_meshsize)
               else
                  thickness_in(i,j,k) = 3.5_WP*cfg%min_meshsize
               end if
            end do 
         end do 
      end do
      call cfg%sync(thickness_in)

      ! allocate(tmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp=3.5_WP*cfg%min_meshsize
      ! do k=cfg%kmin_,cfg%kmax_
      !    do j=cfg%jmin_,cfg%jmax_
      !       do i=cfg%imin_,cfg%imax_
      !          ! Skip wall/bcond/full cells
      !          if (vf%mask(i,j,k).ne.0) cycle
      !          if (vf%VF(i,j,k).lt.VFlo.or.vf%VF(i,j,k).gt.VFhi) cycle
      !          ! Surface-average thickness
      !          area=0.0_WP;tmp(i,j,k)=0.0_WP
      !          do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
      !             area      =area      +vf%SD(ii,jj,kk)*cfg%vol(ii,jj,kk)
      !             tmp(i,j,k)=tmp(i,j,k)+vf%SD(ii,jj,kk)*cfg%vol(ii,jj,kk)*thickness_in(ii,jj,kk)
      !          end do; end do; end do
      !          if (area.gt.0.0_WP) tmp(i,j,k)=tmp(i,j,k)/area
      !       end do
      !    end do
      ! end do
      ! call cfg%sync(tmp)
      ! thickness_in=tmp
      ! deallocate(tmp)


      ! thickness_in=vf%thickness
      minThicknessold=minThickness
      minThickness=3.5_WP*cfg%min_meshsize
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               minThickness=min(minThickness,thickness_in(i,j,k))
            end do 
         end do 
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,minThickness,1,MPI_REAL_WP,MPI_MIN,cfg%comm,ierr)


   end subroutine get_thickness

   subroutine record_thickness()
      use param, only: param_read
      use string,    only: str_medium
      implicit none
      real(WP):: We
      integer :: iunit,ival,ierr,angle
      character(len=str_medium) :: filename
      if (cfg%amRoot) then
         call param_read('Weber number',We)
         ival = nint(We*100)
         call param_read('Mesh aligned angle',angle)
         ! build filename as thickness_pXX.csv
         if (default) then
            write(filename,'("thickness_de_p",I0,"_",I0,"_","dist_",I0,".csv")') ival,angle,mesh_ratio
         else
            write(filename,'("thickness_p",I0,"_",I0,"_","dist_",I0,".csv")') ival,angle,mesh_ratio
            ! write(filename,'("thickness_p",I0,".csv")') ival
            ! write(filename,'("thickness_p",I0,".csv")') ival
         end if
         open(newunit=iunit,file=trim(filename),form='formatted',status='unknown',access='stream',position='append',iostat=ierr)
         write(iunit,*) time%t, minThickness,uzvel,uzvel*time%dt
         close(iunit)
      end if
   end subroutine record_thickness

   subroutine attempt_breakup
      use vfs_class, only : VFhi,VFlo
      implicit none
      integer :: n,i,j,k,m
      call ccl%build(make_label,same_label)
      do n=1,ccl%nstruct
         if ((anew.gt. 1.0_WP .and. minThickness .lt. thickthd1) .or. (anew.lt. 1.0_WP .and. minThickness .lt. thickthd2).or. breakup) then
               do m=1,ccl%struct(n)%n_
                  i=ccl%struct(n)%map(1,m); j=ccl%struct(n)%map(2,m); k=ccl%struct(n)%map(3,m)
                  if (vf%VF(i,j,k).gt.0.5_WP) vf%VF(i,j,k)=1.0_WP
                  ! if (thickness_new(i,j,k).le.cfg%min_meshsize) vf%VF(i,j,k)=1.0_WP
               end do
               breakup=.true.
               exit
         end if
      end do
      if (breakup) then 
         Us=0.0_WP
         Vs=0.0_WP   
         Ws=0.0_WP
      end if
      call cfg%sync(vf%VF)
      contains
         !> Function that identifies cells that need a label
         logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         ! if (vf%VF(i,j,k).gt.0.0_WP) then
         if (vf%thin_sensor(i,j,k).eq.2.0_WP .and. vf%VF(i,j,k).gt.VFlo .and. vf%VF(i,j,k).lt. VFhi) then! .and.thickness_new(i,j,k).le.1.1*cfg%min_meshsize) then
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
   end subroutine attempt_breakup
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
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));resU=0.0_WP
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));resV=0.0_WP
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));resW=0.0_WP
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Ui=0.0_WP
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vi=0.0_WP
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Wi=0.0_WP
         allocate(Us  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Us=0.0_WP 
         allocate(Vs  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Vs=0.0_WP 
         allocate(Ws  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));Ws=0.0_WP 
         allocate(mask_IB  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));mask_IB=0
         allocate(thickness_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));thickness_old=0.0_WP
         allocate(thickness_new(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));thickness_new=0.0_WP
         activated=.false.; breakup=.false.
         amax=0.0_WP;anew=0.0_WP;aold=0.0_WP;uzvel=0.0_WP
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
         use mathtools, only: pi
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: plicnet,elvira,VFhi,VFlo,r2p,r2pnet
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area,Lx,init_dist,ip
         integer:: nx
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with r2p reconstruction
         ! call vf%initialize(cfg=cfg,reconstruction_method=r2pnet,name='VOF')
         ! call vf%initialize(cfg=cfg,reconstruction_method=r2p,name='VOF')
         call vf%initialize(cfg=cfg,reconstruction_method=r2pnet,name='VOF')
         vf%thin_thld_max=1.5_WP
         vf%twoplane_thld2=0.8_WP
         vf%thin_thld_min=0.0_WP
         ! Initialize two droplets
         call param_read('Lx',Lx); call param_read('nx',nx)
         call param_read('Initial Location',init_dist); call param_read('Impact Parameter',ip)
         call param_read('Mesh aligned angle',meshaligned_angle)
         call param_read('Meshratio',mesh_ratio)
         meshaligned_angle=meshaligned_angle/180*pi
         radius1=1.0_WP;radius2=1.0_WP
         center1=init_dist*[-cos(meshaligned_angle),-sin(meshaligned_angle),0.0_WP]+ip*[sin(meshaligned_angle),-cos(meshaligned_angle),0.0_WP]
         ! Case1
         ! if (curr_case.eq.1) then
         !    center2=(init_dist+3.0_WP*Lx/nx)*[cos(meshaligned_angle),sin(meshaligned_angle),0.0_WP]+ip*[-sin(meshaligned_angle),cos(meshaligned_angle),0.0_WP]
         ! else if (curr_case .eq.2) then
         ! ! Case2
         !    center2=(init_dist+2.0_WP*Lx/nx)*[cos(meshaligned_angle),sin(meshaligned_angle),0.0_WP]+ip*[-sin(meshaligned_angle),cos(meshaligned_angle),0.0_WP]
         ! end if
         center2=(init_dist+(2.0_WP+mesh_ratio*0.1_WP)*Lx/nx)*[cos(meshaligned_angle),sin(meshaligned_angle),0.0_WP]+ip*[-sin(meshaligned_angle),cos(meshaligned_angle),0.0_WP]

         ! center1=[-init_dist,-ip,0.0_WP]; center2=[init_dist+Lx/nx,+ip,0.0_WP]
         ! center1=[-1.1_WP/sqrt(2.0_WP),-1.1_WP/sqrt(2.0_WP),0.0_WP]; center2=[1.1_WP/sqrt(2.0_WP),1.1_WP/sqrt(2.0_WP),0.0_WP]
         call param_read('Threshold 1', thickthd1); call param_read('Threshold 2', thickthd2)
         call param_read('default',default)
         call param_read('mulfact',mulfact)
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
         ! Perform interface sensing
         if (vf%two_planes) call vf%sense_interface()
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

      initialize_lpt: block
         lp=lpt(cfg=cfg,name='spray')
         call lp%resize(0)
      end block initialize_lpt

      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=2,nvec=1,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%varname(2)='id'
         pmesh%vecname(1)='velocity'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(1,i)=0.5_WP*lp%p(i)%d
            pmesh%var(2,i)=lp%p(i)%id
            pmesh%vec(:,1,i)=lp%p(i)%vel
         end do
      end block create_pmesh

     
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         integer :: i,j,k
         real(WP) :: Re,We,r,m
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')

         ! Read in adimensional parameters
         call param_read('Reynolds number',Re)   ! Re=rho_l*2U*2R0/mu_l=4/mu_l
         call param_read('Weber number',We)      ! We=2R0rho_l(2U)^2/sigma=8/sigma
         call param_read('Density ratio',r)      ! r=rho_l/rho_g=1/rho_g
         call param_read('Viscosity ratio',m)    ! m=mu_l/mu_g
         call param_read('Hamaker Constant',HamakerC)    ! m=mu_l/mu_g
         call param_read('Mean Free Path',lambdaAir)    
         fs%rho_l=1.0_WP 
         fs%rho_g=fs%rho_l/r
         fs%visc_l=(Re/4.0_WP)**(-1.0_WP)
         fs%visc_g=fs%visc_l/m
         fs%sigma=(We/8.0_WP)**(-1.0_WP)
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
         ! vel1=[+1.0_WP,0.0_WP,0.0_WP]
         ! vel2=[-1.0_WP,0.0_WP,0.0_WP]
         ! vel1=[+1.0_WP/sqrt(2.0_WP),+1.0_WP/sqrt(2.0_WP),0.0_WP]
         ! vel2=[-1.0_WP/sqrt(2.0_WP),-1.0_WP/sqrt(2.0_WP),0.0_WP]
         vel1=[+cos(meshaligned_angle),+sin(meshaligned_angle),0.0_WP]
         vel2=[-cos(meshaligned_angle),-sin(meshaligned_angle),0.0_WP]
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
      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=5,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='thin_sensor'
         smesh%varname(4)='thickness_new'
         smesh%varname(5)='thickness_old'
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
                        smesh%var(3,np)=vf%thin_sensor(i,j,k)
                        smesh%var(4,np)=thickness_new(i,j,k)
                        smesh%var(5,np)=thickness_old(i,j,k)
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
         call ens_out%add_vector('slipvel',resU,resV,resW)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('Indicator',mask_IB)
         call ens_out%add_surface('vofplic',smesh)
         call ens_out%add_particle('part',pmesh)
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
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         fs%U=fs%U+Us
         fs%V=fs%V+Vs
         fs%W=fs%W+Ws
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         fs%U=fs%Uold
         fs%V=fs%Vold
         fs%W=fs%Wold
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! ! Get the two thickness to integrate gas pressure field
         thickness_old=thickness_new
         call get_thickness(thickness_new)
         ! if (.not. breakup) then
         if (.not. default) then
            call get_gasVel()
            resU=fs%U+Us
            resV=fs%V+Vs
            resW=fs%W+Ws
            call vf%advance(dt=time%dt,U=resU,V=resV,W=resW)
            ! call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         else
         ! VOF solver ste
            ! Remember old VOF
            vf%VFold=vf%VF
            call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         end if
         
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
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
            call fs%add_surface_tension_jump_twoVF(dt=time%dt,div=fs%div,vf=vf)
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
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         call record_thickness()
         ! call attempt_breakup()
         ! !> Calculate the interpolated velocity, including overlap and ghosts
         interp_vel: block         
            integer :: i,j,k
            ! Calculate as far as possible each component
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_-1
                     resU(i,j,k)=sum(fs%itpu_x(:,i,j,k)*Us(i:i+1,j,k))
                  end do
               end do
            end do
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_-1
                  do i=cfg%imino_,cfg%imaxo_
                     resV(i,j,k)=sum(fs%itpv_y(:,i,j,k)*Vs(i,j:j+1,k))
                  end do
               end do
            end do
            do k=cfg%kmino_,cfg%kmaxo_-1
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     resW(i,j,k)=sum(fs%itpw_z(:,i,j,k)*Ws(i,j,k:k+1))
                  end do
               end do
            end do
            ! Add last layer in each direction
            if (.not.cfg%xper.and.cfg%iproc.eq.cfg%npx) resU(cfg%imaxo,:,:)=Us(cfg%imaxo,:,:)
            if (.not.cfg%yper.and.cfg%jproc.eq.cfg%npy) resV(:,cfg%jmaxo,:)=Vs(:,cfg%jmaxo,:)
            if (.not.cfg%zper.and.cfg%kproc.eq.cfg%npz) resW(:,:,cfg%kmaxo)=Ws(:,:,cfg%kmaxo)
            ! Sync it
            call cfg%sync(resU)
            call cfg%sync(resV)
            call cfg%sync(resW)
         end block interp_vel

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
                              smesh%var(3,np)=vf%thin_sensor(i,j,k)
                              smesh%var(4,np)=thickness_new(i,j,k)
                              smesh%var(5,np)=thickness_old(i,j,k)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh

            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=0.5_WP*lp%p(i)%d
                  pmesh%var(2,i)=lp%p(i)%id
                  pmesh%vec(:,1,i)=lp%p(i)%vel
               end do
            end block update_pmesh 

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