!> Abstract base class for AMR solvers
!> Provides context for coupled BC callbacks, regrid handling, and checkpoint registration
module amrsolver_class
   use iso_c_binding, only: c_ptr
   use precision,     only: WP
   use string,        only: str_medium
   use amrio_class,   only: amrio
   implicit none
   private

   public :: amrsolver

   !> Abstract solver base - extended by concrete solvers (amrscalar, mpcomp, etc.)
   type, abstract :: amrsolver
      ! Name for logging/debugging
      character(len=str_medium) :: name='UNNAMED_SOLVER'
   contains
      ! Deferred procedures that concrete solvers must implement
      procedure(fillbc_interface), deferred :: fillbc              !< Apply physical BCs
      procedure(regrid_interface), deferred :: on_regrid           !< Handle regrid event
      procedure(register_checkpoint_interface), deferred :: register_checkpoint !< Register fields for checkpoint
      procedure(restore_checkpoint_interface), deferred :: restore_checkpoint   !< Restore fields from checkpoint
   end type amrsolver

   !> Interface for fillbc - receives solver context, can access all fields
   abstract interface
      subroutine fillbc_interface(this,mf_ptr,geom_ptr,time,scomp,ncomp)
         import :: amrsolver,c_ptr,WP
         class(amrsolver), intent(inout) :: this
         type(c_ptr), value :: mf_ptr     !< MultiFab pointer
         type(c_ptr), value :: geom_ptr   !< Geometry pointer
         real(WP), value :: time
         integer, value :: scomp,ncomp
      end subroutine fillbc_interface

      subroutine regrid_interface(this,lvl,time)
         import :: amrsolver,WP
         class(amrsolver), intent(inout) :: this
         integer, value :: lvl
         real(WP), value :: time
      end subroutine regrid_interface

      subroutine register_checkpoint_interface(this,io)
         import :: amrsolver,amrio
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
      end subroutine register_checkpoint_interface

      subroutine restore_checkpoint_interface(this,io,dirname)
         import :: amrsolver,amrio
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
         character(len=*), intent(in) :: dirname
      end subroutine restore_checkpoint_interface
   end interface

end module amrsolver_class
