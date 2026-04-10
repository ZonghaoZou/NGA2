!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use nozzle_class,      only: nozzle
   implicit none
   private
   
   !> Injector simulation
   type(nozzle) :: injector

   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      call injector%init()
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      do while (.not. injector%time%done()) 
         call injector%step()
      end do 
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none


      ! Finalize injector simulation
      call injector%final()
         
   end subroutine simulation_final
   

end module simulation
