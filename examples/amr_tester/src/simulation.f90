!> Simulation hooks for AMR Tester
module simulation
   use mod_test_registry,      only: test_registry
   use mod_test_visualization, only: test_visualization
   use mod_test_amrscalar,     only: test_amrscalar
   use mod_test_multidata,     only: test_multidata
   use mod_test_amrmg,         only: test_amrmg
   use mod_test_amrabeclap,    only: test_amrabeclap
   use messager, only: log
   implicit none

contains

   !> Initialization hook
   subroutine simulation_init()
      call log("Simulation Init: AMR Tester")
   end subroutine simulation_init

   !> Clean run hook
   subroutine simulation_run()
      call log("Simulation Run: Executing Tests")
      ! call test_multidata()  ! DISABLED FOR DEBUGGING
      ! call test_registry()   ! DISABLED FOR DEBUGGING
      ! call test_visualization() ! DISABLED FOR DEBUGGING
      ! call test_amrscalar()  ! DISABLED FOR DEBUGGING
      ! call test_amrmg()      ! Constant coeff Poisson
      call test_amrabeclap()   ! Variable coeff solver
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
   end subroutine simulation_final

end module simulation
