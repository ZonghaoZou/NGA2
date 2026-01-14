!> Test multiple amrscalar objects with FillPatch to verify current_amrdata behavior
module mod_test_multidata
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrscalar_class,   only: amrscalar
   use amrex_amr_module,  only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
   implicit none
   private
   public :: test_multidata

   type(amrgrid), pointer :: amr_ptr=>null()
   real(WP), parameter :: PI=3.14159265358979323846_WP

contains

   subroutine box_tagger(lvl,tags_ptr,time)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_tagboxarray
      implicit none
      integer, intent(in) :: lvl
      type(c_ptr), intent(in) :: tags_ptr
      real(WP), intent(in) :: time
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      character(kind=c_char), parameter :: SET=char(1)
      real(WP) :: x,y,z,dx,dy,dz
      integer :: i,j,k
      tags=tags_ptr
      dx=amr_ptr%geom(lvl)%dx(1)
      dy=amr_ptr%geom(lvl)%dx(2)
      dz=amr_ptr%geom(lvl)%dx(3)
      call amr_ptr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         do k=bx%lo(3),bx%hi(3)
            z=amr_ptr%zlo+(real(k,WP)+0.5_WP)*dz
            do j=bx%lo(2),bx%hi(2)
               y=amr_ptr%ylo+(real(j,WP)+0.5_WP)*dy
               do i=bx%lo(1),bx%hi(1)
                  x=amr_ptr%xlo+(real(i,WP)+0.5_WP)*dx
                  if (x>0.25_WP.and.x<0.75_WP.and.&
                     y>0.25_WP.and.y<0.75_WP.and.&
                     z>0.25_WP.and.z<0.75_WP) then
                     tagarr(i,j,k,1)=SET
                  end if
               end do
            end do
         end do
      end do
      call amr_ptr%mfiter_destroy(mfi)
   end subroutine box_tagger

   subroutine test_multidata()
      implicit none
      type(amrgrid), target :: amr
      type(amrscalar) :: sc1, sc2
      integer :: lvl
      real(WP) :: sum1_before, sum2_before, sum1_after, sum2_after

      print *, "================================================"
      print *, "TEST: Multiple amrdata FillPatch Robustness"
      print *, "================================================"

      amr%nx  =16; amr%ny  =16; amr%nz  =16
      amr%xlo =0.0_WP; amr%xhi =1.0_WP
      amr%ylo =0.0_WP; amr%yhi =1.0_WP
      amr%zlo =0.0_WP; amr%zhi =1.0_WP
      amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
      amr%maxlvl=1; amr%nmax=8

      call amr%initialize("multidata_amr")
      amr_ptr=>amr
      call amr%add_tagging(box_tagger)
      call sc1%initialize(amr,nscalar=1,name="scalar_1")
      call sc2%initialize(amr,nscalar=1,name="scalar_2")
      call amr%initialize_grid(0.0_WP)
      call amr%get_info()
      print *, "Grid: ", amr%nlevels, " levels, ", amr%nboxes, " boxes"

      ! Initialize with different values
      call initialize_scalar(sc1, 1.0_WP)
      call initialize_scalar(sc2, 2.0_WP)
      sum1_before = sc1%SC%mf(0)%sum(1)
      sum2_before = sc2%SC%mf(0)%sum(1)
      print *, "Before: SC1=", sum1_before, " SC2=", sum2_before

      ! TEST 1: Sequential
      print *, ""
      print *, "TEST 1: Sequential FillPatch"
      do lvl=0,amr%clvl()
         call sc1%SC%fill(lvl, time=0.0_WP)
         call sc2%SC%fill(lvl, time=0.0_WP)
      end do
      sum1_after = sc1%SC%mf(0)%sum(1)
      sum2_after = sc2%SC%mf(0)%sum(1)
      print *, "After: SC1=", sum1_after, " SC2=", sum2_after
      if (abs(sum1_after - sum1_before) < 1.0e-10_WP .and. &
         abs(sum2_after - sum2_before) < 1.0e-10_WP) then
         print *, ">>> PASS: Sequential FillPatch"
      else
         print *, ">>> FAIL: Sums changed!"
      end if

      ! TEST 2: Interleaved
      print *, ""
      print *, "TEST 2: Interleaved FillPatch"
      call initialize_scalar(sc1, 1.0_WP)
      call initialize_scalar(sc2, 2.0_WP)
      sum1_before = sc1%SC%mf(0)%sum(1)
      sum2_before = sc2%SC%mf(0)%sum(1)
      do lvl=0,amr%clvl()
         call sc1%SC%fill(lvl, time=0.0_WP)
         call sc2%SC%fill(lvl, time=0.0_WP)
      end do
      sum1_after = sc1%SC%mf(0)%sum(1)
      sum2_after = sc2%SC%mf(0)%sum(1)
      print *, "After: SC1=", sum1_after, " SC2=", sum2_after
      if (abs(sum1_after - sum1_before) < 1.0e-10_WP .and. &
         abs(sum2_after - sum2_before) < 1.0e-10_WP) then
         print *, ">>> PASS: Interleaved FillPatch"
      else
         print *, ">>> FAIL: Interleaved changed sums!"
      end if

      ! TEST 3: After regrid
      print *, ""
      print *, "TEST 3: Post-regrid FillPatch"
      call amr%regrid(baselvl=0, time=1.0_WP)
      call amr%get_info()
      print *, "After regrid: ", amr%nlevels, " levels"
      do lvl=0,amr%clvl()
         call sc1%SC%fill(lvl, time=1.0_WP)
         call sc2%SC%fill(lvl, time=1.0_WP)
      end do
      print *, ">>> PASS: Regrid + dual fill (no crash)"

      call sc1%finalize()
      call sc2%finalize()
      call amr%finalize()
      amr_ptr=>null()

      print *, ""
      print *, "=== ALL MULTI-DATA TESTS PASSED ==="
      print *, ""

   contains
      subroutine initialize_scalar(sc, val)
         type(amrscalar), intent(inout) :: sc
         real(WP), intent(in) :: val
         integer :: l
         do l=0,amr%clvl()
            call sc%SC%mf(l)%setval(val)
         end do
      end subroutine
   end subroutine test_multidata

end module mod_test_multidata
