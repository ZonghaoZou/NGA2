!> Flux register class wrapping AMReX FluxRegister
!> Designed to be managed by amrgrid Registry
module amrflux_class
   use iso_c_binding,    only: c_ptr,c_associated
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_fluxregister,amrex_fluxregister_build,&
   &                           amrex_fluxregister_destroy,amrex_multifab,&
   &                           amrex_boxarray,amrex_distromap,amrex_spacedim,&
   &                           amrex_real,amrex_geometry
   implicit none
   private

   public :: amrflux

   !> Flux register object wrapping array of FluxRegisters (one per fine level)
   type :: amrflux
      ! Underlying AMReX objects (array 1:maxlvl, no level 0)
      type(amrex_fluxregister), dimension(:), allocatable :: fr
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRFLUX'
      integer :: ncomp=1   !< Number of components
   contains
      procedure :: initialize    !< Initialize amrflux with parameters
      procedure :: build_level
      procedure :: destroy_level
      procedure :: destroy
      ! Wrapper methods that forward to underlying fluxregister
      procedure :: crseinit
      procedure :: fineadd
      procedure :: reflux
      procedure :: setval
   end type amrflux

contains

   !> Initialize amrflux with parameters
   !> Allocates the fr array
   subroutine initialize(this, maxlvl, name, ncomp)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: maxlvl              !< Max level (from amr%maxlvl)
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      ! Set metadata
      this%name = name
      this%ncomp = ncomp
      ! Allocate fr array (flux registers only exist for levels >= 1)
      if (.not.allocated(this%fr)) allocate(this%fr(1:maxlvl))
   end subroutine initialize


   !> Build flux register at a specific fine level
   !> Requires BoxArray and DistroMap of the fine level
   subroutine build_level(this,lvl,ba,dm,ref_ratio)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      integer, intent(in) :: ref_ratio
      ! Flux registers only exist for levels >= 1
      if (lvl.lt.1) return
      if (.not.allocated(this%fr)) return
      if (lvl.gt.ubound(this%fr,1)) return
      ! Destroy if already built
      call this%destroy_level(lvl)
      ! Build new register
      call amrex_fluxregister_build(this%fr(lvl),ba,dm,ref_ratio,lvl,this%ncomp)
   end subroutine build_level

   !> Destroy flux register at a specific level
   subroutine destroy_level(this,lvl)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      if (.not.allocated(this%fr)) return
      if (lvl.lt.1.or.lvl.gt.ubound(this%fr,1)) return
      if (c_associated(this%fr(lvl)%p)) call amrex_fluxregister_destroy(this%fr(lvl))
   end subroutine destroy_level

   !> Destroy all flux registers
   subroutine destroy(this)
      class(amrflux), intent(inout) :: this
      integer :: i
      if (allocated(this%fr)) then
         do i=1,ubound(this%fr,1)
            if (c_associated(this%fr(i)%p)) call amrex_fluxregister_destroy(this%fr(i))
         end do
         deallocate(this%fr)
      end if
   end subroutine destroy

   !> Initialize coarse-side fluxes (resets the register)
   subroutine crseinit(this,lvl,fluxes,scale)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
      real(amrex_real), intent(in) :: scale
      if (.not.allocated(this%fr)) return
      if (lvl.lt.1.or.lvl.gt.ubound(this%fr,1)) return
      if (c_associated(this%fr(lvl)%p)) call this%fr(lvl)%crseinit(fluxes,scale)
   end subroutine crseinit

   !> Add fine-side fluxes
   subroutine fineadd(this,lvl,fluxes,scale)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
      real(amrex_real), intent(in) :: scale
      if (.not.allocated(this%fr)) return
      if (lvl.lt.1.or.lvl.gt.ubound(this%fr,1)) return
      if (c_associated(this%fr(lvl)%p)) call this%fr(lvl)%fineadd(fluxes,scale)
   end subroutine fineadd

   !> Apply correction to coarse-level data
   subroutine reflux(this,lvl,mf,scale,geom)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      real(amrex_real), intent(in) :: scale
      type(amrex_geometry), intent(in) :: geom
      if (.not.allocated(this%fr)) return
      if (lvl.lt.1.or.lvl.gt.ubound(this%fr,1)) return
      if (c_associated(this%fr(lvl)%p)) call this%fr(lvl)%reflux(mf,scale)
   end subroutine reflux

   !> Set all values in the register
   subroutine setval(this,lvl,val)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      real(amrex_real), intent(in) :: val
      if (.not.allocated(this%fr)) return
      if (lvl.lt.1.or.lvl.gt.ubound(this%fr,1)) return
      if (c_associated(this%fr(lvl)%p)) call this%fr(lvl)%setval(val)
   end subroutine setval

end module amrflux_class
