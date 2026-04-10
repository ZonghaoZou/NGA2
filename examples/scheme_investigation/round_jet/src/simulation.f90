!> Various definitions and tools for running an NGA2 simulation
module simulation
   use roundjet_class,     only: roundjet
   use coupler_class,      only: coupler
   implicit none
   private
   
   !> Round jet simulation
   type(roundjet) :: jet

   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none

      call jet%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Jet drives overall time integration
      do while (.not.jet%time%done())
         
         ! coupling: block
         !    use tpns_class, only: bcond
         !    integer :: n,i,j,k
         !    type(bcond), pointer :: mybc
         !    ! Exchange data using coupler
         !    if (in_pipe_group) call cpl%push(pipe%fs%U,loc='x'); call cpl%transfer(); call cpl%pull(jet%resU,loc='x')
         !    if (in_pipe_group) call cpl%push(pipe%fs%V,loc='y'); call cpl%transfer(); call cpl%pull(jet%resV,loc='y')
         !    if (in_pipe_group) call cpl%push(pipe%fs%W,loc='z'); call cpl%transfer(); call cpl%pull(jet%resW,loc='z')
         !    ! Apply time-varying Dirichlet conditions
         !    call jet%fs%get_bcond('inflow',mybc)
         !    do n=1,mybc%itr%no_
         !       i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !       jet%fs%U(i  ,j,k)=jet%resU(i  ,j,k)
         !       jet%fs%V(i-1,j,k)=jet%resV(i-1,j,k)
         !       jet%fs%W(i-1,j,k)=jet%resW(i-1,j,k)
         !    end do
         ! end block coupling
         
         ! Advance jet simulation
         call jet%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize jet simulation
      call jet%final()
      
   end subroutine simulation_final
   
   
end module simulation
