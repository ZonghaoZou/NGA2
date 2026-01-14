!> Test for solver callback pattern
module mod_test_registry
   use iso_c_binding
   use amrgrid_class
   use amrdata_class
   use messager, only: log,warn
   use precision, only: WP
   use amrex_amr_module, only: amrex_boxarray,amrex_distromap
   implicit none
   private

   public :: test_registry

   ! Module-level data for callback access
   type(amrdata), pointer :: velocity_ptr=>null()

contains

   !> Callback for init level
   subroutine on_init_level(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: data
      call c_f_pointer(ctx,data)
      call data%define(lvl,ba,dm)
   end subroutine on_init_level

   !> Callback for clear level
   subroutine on_clear_level(ctx,lvl)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrdata), pointer :: data
      call c_f_pointer(ctx,data)
      call data%clear_level(lvl)
   end subroutine on_clear_level

   subroutine test_registry()
      type(amrgrid) :: amr
      type(amrdata), target :: velocity

      call log("---------------------------------------------------")
      call log("Running Test: Solver Callback Pattern")
      call log("---------------------------------------------------")

      ! 1. Initialize AMR Config
      amr%nx=32
      amr%ny=32
      amr%nz=32
      amr%xlo=0.0_WP
      amr%xhi=1.0_WP
      amr%ylo=0.0_WP
      amr%yhi=1.0_WP
      amr%zlo=0.0_WP
      amr%zhi=1.0_WP
      amr%maxlvl=0
      amr%nmax=32

      call amr%initialize("TestGrid")
      call log("AMR Config Initialized.")

      ! 2. Configure amrdata
      velocity%name="velocity"
      velocity%ncomp=3
      velocity%ng=1
      allocate(velocity%mf(0:amr%maxlvl))
      velocity_ptr=>velocity
      call log("Data 'velocity' configured (ncomp=3, ng=1).")

      ! 3. Register callbacks (pass data object as context)
      call amr%add_on_init(on_init_level,c_loc(velocity))
      call amr%add_on_clear(on_clear_level,c_loc(velocity))
      call log("Callbacks registered.")

      ! 4. Initialize Grid (Triggers callbacks)
      call amr%initialize_grid(0.0d0)
      call log("Grid Initialized.")

      ! 5. Verify Allocation
      if (allocated(velocity%mf)) then
         if (c_associated(velocity%mf(0)%p)) then
            call log("PASS: Data allocated via callback!")
         else
            call warn("FAIL: Data MF pointer is null!")
         end if
      else
         call warn("FAIL: Data MF array not allocated!")
      end if

      ! Cleanup
      velocity_ptr=>null()
      call amr%finalize()
      call log("PASS: Callback pattern test complete.")

   end subroutine test_registry

end module mod_test_registry
