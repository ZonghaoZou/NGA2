module orr_sommerfeld_kh
    use precision, only: WP
    implicit none
    ! --- State Variables from GEVP Solve ---
    integer :: N_cheb
    real(WP) :: alpha_wave, Ug_inf, Ul_inf, del_g, del_l
    real(WP) :: Lg_dom, Ll_dom
    real(WP) :: stretch_g_nd, stretch_l_nd  ! ND stretches for chebev bounds
    complex(WP) :: c_phase                  ! Eigenvalue (dimensional phase speed)
    complex(WP), allocatable :: V_kh(:)     ! Full eigenvector
    complex(WP), allocatable :: phi_l_kh(:), phi_g_kh(:)  ! Split eigenvector
    real(WP), allocatable :: Tp_kh(:,:)     ! Inverse Chebyshev transform
    real(WP), allocatable :: D1_l_kh(:,:), D1_g_kh(:,:)  ! 1st derivative matrices (ND)
    logical :: printed_params = .false.

    ! --- Precomputed Chebyshev coefficients for fast eval_kh ---
    ! These are computed once in setup_kh, then eval_kh uses only Clenshaw recurrence
    real(WP), allocatable :: cheb_dphi_r_l(:), cheb_dphi_i_l(:)  ! D1*phi liquid: real/imag
    real(WP), allocatable :: cheb_dphi_r_g(:), cheb_dphi_i_g(:)  ! D1*phi gas: real/imag
    real(WP), allocatable :: cheb_phi_r_l(:), cheb_phi_i_l(:)    ! phi liquid: real/imag
    real(WP), allocatable :: cheb_phi_r_g(:), cheb_phi_i_g(:)    ! phi gas: real/imag
    real(WP) :: phi_r_interface, phi_i_interface  ! gas eigenfunction at y=0
    real(WP) :: G_coeff1, G_coeff2  ! precomputed level set coefficients

contains

    !> Orr-Sommerfeld GEVP solver for two-phase Kelvin-Helmholtz instability
    !! Ported from the old kelvin_helmholtz.f90 (Boomkamp 1997 formulation)
    !! Eigenvalue is omega = alpha*c (angular frequency)
    !! Output c_eval = omega/alpha (dimensional phase speed)
    subroutine setup_kh(N, alpha_in, Ug_inf_in, Ul_inf_in, del_g_in, del_l_in, &
                        rho_g_in, rho_l_in, mu_g_in, mu_l_in, sigma_in,        &
                        Lg_in, Ll_in, c_eval)
        use mathtools, only: Pi
        implicit none
        ! Inputs
        integer,  intent(in)  :: N
        real(WP), intent(in)  :: alpha_in
        real(WP), intent(in)  :: Ug_inf_in, Ul_inf_in
        real(WP), intent(in)  :: del_g_in, del_l_in
        real(WP), intent(in)  :: rho_g_in, rho_l_in
        real(WP), intent(in)  :: mu_g_in, mu_l_in
        real(WP), intent(in)  :: sigma_in
        real(WP), intent(in)  :: Lg_in, Ll_in
        complex(WP), intent(out) :: c_eval

        ! Internal grid size (1-indexed, nos = N+1)
        integer :: nos

        ! Chebyshev grids
        real(WP), allocatable :: zos(:), yos_l(:), yos_g(:)

        ! Differentiation matrices
        real(WP), allocatable :: Id(:,:)
        real(WP), allocatable :: Tmat(:,:), Gos(:,:), Tp(:,:), Gtmp(:,:)
        real(WP), allocatable :: D1_g(:,:), D2_g(:,:), D3_g(:,:), D4_g(:,:)
        real(WP), allocatable :: D1_l(:,:), D2_l(:,:), D3_l(:,:), D4_l(:,:)

        ! Base flow
        real(WP), allocatable :: Ubase_l(:), Up_l(:), Upp_l(:)
        real(WP), allocatable :: Ubase_g(:), Up_g(:), Upp_g(:)

        ! GEVP matrices (full 2*nos and reduced 2*(nos-1))
        complex(WP), allocatable :: A(:,:), B(:,:)
        complex(WP), allocatable :: Ain(:,:), Bin(:,:)

        ! LAPACK
        integer :: ierr, lwork
        complex(WP), allocatable :: work(:)
        real(WP),    allocatable :: rwork(:)
        complex(WP), allocatable :: eigval_a(:), eigval_b(:)
        complex(WP), allocatable :: eigvec_l(:,:), eigvec_r(:,:)

        ! Flow parameters
        real(WP) :: m_visc, r_dens, Re_l, S_st, F_grav
        real(WP) :: stretch_g, stretch_l

        ! Complex constants
        complex(WP), parameter :: zero = (0.0_WP, 0.0_WP)
        complex(WP), parameter :: ii   = (0.0_WP, 1.0_WP)

        ! Various
        integer  :: i, j, jj, mode_index
        integer  :: iliq_min, iliq_max, igas_min, igas_max, igas, iliq
        real(WP) :: tmp, a2, a4, alpha_nd
        complex(WP) :: omega, omega_

        ! =============================================
        ! Setup parameters -- NON-DIMENSIONALIZE by U_g and delta_g
        ! This makes the GEVP coefficient iα·Re_l correct
        ! (old code assumes U_g·δ_g = 1, which holds in ND space)
        ! =============================================
        nos = N + 1

        ! Non-dimensional parameters
        m_visc  = mu_g_in / mu_l_in
        r_dens  = rho_g_in / rho_l_in
        Re_l    = rho_l_in * Ug_inf_in * del_g_in / mu_l_in
        if (sigma_in > 1.0e-12_WP) then
            S_st = sigma_in / (rho_l_in * Ug_inf_in**2 * del_g_in)
        else
            S_st = 0.0_WP
        end if
        F_grav = 0.0_WP

        ! Non-dimensional stretch (in units of delta_g)
        stretch_g = Lg_in / del_g_in
        stretch_l = Ll_in / del_g_in

        ! Non-dimensional wavenumber and powers
        alpha_nd = alpha_in * del_g_in
        a2 = alpha_nd**2
        a4 = alpha_nd**4

        ! Store module-level state (keep dimensional for external use)
        N_cheb     = N
        alpha_wave = alpha_in
        Ug_inf     = Ug_inf_in
        Ul_inf     = Ul_inf_in
        del_g      = del_g_in
        del_l      = del_l_in
        Lg_dom     = Lg_in

        ! One-time diagnostic print
        ! if (.not.printed_params) then
        !     printed_params = .true.
        !     print '(A)',       '=== Orr-Sommerfeld Parameters ==='
        !     print '(A,ES12.4)', '  rho_g        = ', rho_g_in
        !     print '(A,ES12.4)', '  rho_l        = ', rho_l_in
        !     print '(A,ES12.4)', '  mu_g         = ', mu_g_in
        !     print '(A,ES12.4)', '  mu_l         = ', mu_l_in
        !     print '(A,ES12.4)', '  sigma        = ', sigma_in
        !     print '(A,ES12.4)', '  U_g          = ', Ug_inf_in
        !     print '(A,ES12.4)', '  U_l          = ', Ul_inf_in
        !     print '(A,ES12.4)', '  delta_g      = ', del_g_in
        !     print '(A,ES12.4)', '  delta_l      = ', del_l_in
        !     print '(A,ES12.4)', '  m  (mu_g/mu_l)  = ', m_visc
        !     print '(A,ES12.4)', '  r  (rho_g/rho_l)= ', r_dens
        !     print '(A,ES12.4)', '  Re_l         = ', Re_l
        !     print '(A,ES12.4)', '  Re_g         = ', rho_g_in*Ug_inf_in*del_g_in/mu_g_in
        !     print '(A,ES12.4)', '  S  (sigma)   = ', S_st
        !     print '(A,ES12.4)', '  F  (gravity) = ', F_grav
        !     print '(A,ES12.4)', '  alpha_nd     = ', alpha_in*del_g_in
        !     print '(A,ES12.4)', '  stretch_g_nd = ', stretch_g
        !     print '(A,ES12.4)', '  stretch_l_nd = ', stretch_l
        !     print '(A,I6)',     '  N_cheb       = ', N
        !     print '(A)',       '================================='
        ! end if
        Ll_dom     = Ll_in

        ! =============================================
        ! Create Chebyshev grids (NON-DIMENSIONAL: y/delta_g)
        !   Liquid: yos_l from -stretch_l to 0
        !   Gas:    yos_g from  0 to stretch_g
        ! =============================================
        allocate(zos(nos), yos_l(nos), yos_g(nos))
        do j = 1, nos
            zos(j) = cos(Pi * real(j-1, WP) / real(nos-1, WP))
        end do
        do j = 1, nos
            yos_l(j) = stretch_l * (-zos(j) - 1.0_WP) * 0.5_WP
            yos_g(j) = stretch_g * (-zos(j) + 1.0_WP) * 0.5_WP
        end do

        ! =============================================
        ! Build D matrices via Chebyshev T·G·Tp approach
        ! D matrices are d/d(y/delta_g) -- NON-DIMENSIONAL
        ! =============================================
        allocate(Id(nos,nos), Tmat(nos,nos), Gos(nos,nos), Tp(nos,nos), Gtmp(nos,nos))
        allocate(D1_g(nos,nos), D2_g(nos,nos), D3_g(nos,nos), D4_g(nos,nos))
        allocate(D1_l(nos,nos), D2_l(nos,nos), D3_l(nos,nos), D4_l(nos,nos))

        Id = 0.0_WP
        D1_l = 0.0_WP; D1_g = 0.0_WP

        do i = 1, nos
            Id(i,i) = 1.0_WP
            do j = 1, nos
                Tmat(i,j) = cos(Pi * real((i-1)*(j-1), WP) / real(nos-1, WP))
                Tp(i,j) = 2.0_WP * cos(Pi * real((i-1)*(j-1), WP) / real(nos-1, WP)) &
                         / real(nos-1, WP)
                if (i == 1 .or. i == nos) Tp(i,j) = 0.5_WP * Tp(i,j)
                if (j == 1 .or. j == nos) Tp(i,j) = 0.5_WP * Tp(i,j)
                if (i >= j .or. mod(i+j, 2) == 0) then
                    Gos(i,j) = 0.0_WP
                else
                    if (i == 1 .or. i == nos) then
                        Gos(i,j) = real(j-1, WP)
                    else
                        Gos(i,j) = 2.0_WP * real(j-1, WP)
                    end if
                end if
            end do
        end do

        ! D1 in z-space: D1 = T * G * Tp
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, Tmat, nos, Gos, nos, 0.0_WP, Gtmp, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, Gtmp, nos, Tp,  nos, 0.0_WP, D1_g, nos)
        D1_l = D1_g

        ! Map to ND coordinates: d/d(y_nd) = (dz/dy_nd) * d/dz = (-2/stretch_nd) * d/dz
        do i = 1, nos
            D1_l(i,:) = D1_l(i,:) / (-stretch_l * 0.5_WP)
            D1_g(i,:) = D1_g(i,:) / (-stretch_g * 0.5_WP)
        end do

        ! Higher order by matrix multiplication
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D1_g, nos, D1_g, nos, 0.0_WP, D2_g, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D1_l, nos, D1_l, nos, 0.0_WP, D2_l, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D1_g, nos, D2_g, nos, 0.0_WP, D3_g, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D1_l, nos, D2_l, nos, 0.0_WP, D3_l, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D2_g, nos, D2_g, nos, 0.0_WP, D4_g, nos)
        call dgemm('N', 'N', nos, nos, nos, 1.0_WP, D2_l, nos, D2_l, nos, 0.0_WP, D4_l, nos)

        ! =============================================
        ! Create base flow (NON-DIMENSIONAL: U/U_g)
        ! y_nd = y/delta_g, delta_l_nd = delta_l/delta_g
        ! =============================================
        allocate(Ubase_l(nos), Up_l(nos), Upp_l(nos))
        allocate(Ubase_g(nos), Up_g(nos), Upp_g(nos))
        block
            real(WP) :: delta_l_nd, U_l_nd
            delta_l_nd = del_l_in / del_g_in
            U_l_nd     = Ul_inf_in / Ug_inf_in
            do j = 1, nos
                Ubase_l(j) = U_l_nd * erf(yos_l(j) / delta_l_nd)
                Ubase_g(j) = 1.0_WP * erf(yos_g(j) / 1.0_WP)
                Up_l(j)    = (2.0_WP * U_l_nd) / &
                             (sqrt(Pi) * delta_l_nd * exp(yos_l(j)**2 / delta_l_nd**2))
                Up_g(j)    = (2.0_WP * 1.0_WP) / &
                             (sqrt(Pi) * 1.0_WP * exp(yos_g(j)**2 / 1.0_WP**2))
                Upp_l(j)   = -1.0_WP * (4.0_WP * U_l_nd * yos_l(j)) / &
                             (sqrt(Pi) * delta_l_nd**3 * exp(yos_l(j)**2 / delta_l_nd**2))
                Upp_g(j)   = -1.0_WP * (4.0_WP * 1.0_WP * yos_g(j)) / &
                             (sqrt(Pi) * 1.0_WP**3 * exp(yos_g(j)**2 / 1.0_WP**2))
            end do
        end block

        ! =============================================
        ! Assemble GEVP: A * phi = omega * B * phi
        ! Eigenvalue omega = alpha * c (angular frequency)
        ! Ref: Boomkamp 1997, Chebyshev Collocation Method
        ! =============================================
        allocate(A(2*nos, 2*nos), B(2*nos, 2*nos))
        A = zero; B = zero

        iliq_min = 1;       iliq_max = nos
        igas_min = nos + 1; igas_max = 2 * nos

        ! --- Interior: Liquid Orr-Sommerfeld (j=3..nos-2) ---
        do j = 3, nos - 2
            A(j, iliq_min:iliq_max) = &
                  D4_l(j,:) - 2.0_WP*a2*D2_l(j,:) + a4*Id(j,:) &
                + ii*( alpha_nd*Re_l*Upp_l(j)*Id(j,:)            &
                     - alpha_nd*Re_l*Ubase_l(j)*(D2_l(j,:) - a2*Id(j,:)))
            B(j, iliq_min:iliq_max) = &
                  ii*(-alpha_nd*Re_l*(D2_l(j,:) - a2*Id(j,:)))
        end do

        ! --- Interior: Gas Orr-Sommerfeld (j=3..nos-2) ---
        do j = 3, nos - 2
            A(j+nos, igas_min:igas_max) = &
                  D4_g(j,:) - 2.0_WP*a2*D2_g(j,:) + a4*Id(j,:)            &
                + ii*( alpha_nd*Re_l*(r_dens/m_visc)*Upp_g(j)*Id(j,:)      &
                     - alpha_nd*Re_l*(r_dens/m_visc)*Ubase_g(j)            &
                       *(D2_g(j,:) - a2*Id(j,:)))
            B(j+nos, igas_min:igas_max) = &
                  ii*(-alpha_nd*Re_l*(r_dens/m_visc)*(D2_g(j,:) - a2*Id(j,:)))
        end do

        ! --- Far-field BCs ---
        ! phi=0 at boundaries: enforced by matrix reduction (rows/cols 1 and 2*nos stripped)
        ! phi'=0 at boundaries:
        igas = nos; iliq = 1
        A(2,        iliq_min:iliq_max) = D1_l(iliq,:)
        A(2*nos-1,  igas_min:igas_max) = D1_g(igas,:)

        ! --- Interface conditions (liquid: j=nos, gas: j=1) ---
        igas = 1; iliq = nos

        ! 1. Continuity of streamfunction: phi_l - phi_g = 0
        A(nos-1, iliq_min:iliq_max) =  Id(iliq,:)
        A(nos-1, igas_min:igas_max) = -Id(igas,:)

        ! 2. Tangential velocity: U_l'*phi_l - U_g'*phi_g = omega*(phi_g' - phi_l')
        A(nos,   iliq_min:iliq_max) =  Up_l(iliq)*Id(iliq,:)
        A(nos,   igas_min:igas_max) = -Up_g(igas)*Id(igas,:)
        B(nos,   iliq_min:iliq_max) = -D1_l(iliq,:)
        B(nos,   igas_min:igas_max) =  D1_g(igas,:)

        ! 3. Tangential stress: (D2_l + a2)*phi_l - m*(D2_g + a2)*phi_g = 0
        A(nos+1, iliq_min:iliq_max) =        (D2_l(iliq,:) + a2*Id(iliq,:))
        A(nos+1, igas_min:igas_max) = -m_visc*(D2_g(igas,:) + a2*Id(igas,:))

        ! 4. Normal stress (Boomkamp Eq.)
        A(nos+2, iliq_min:iliq_max) = &
             (D3_l(iliq,:) - 3.0_WP*a2*D1_l(iliq,:))                   &
           + ii*alpha_nd*Re_l*Up_l(iliq)*Id(iliq,:)                     &
           + ii*alpha_nd*Re_l*(F_grav + a2*S_st)                        &
             *D1_l(iliq,:)/(Up_l(iliq) - Up_g(igas))
        A(nos+2, igas_min:igas_max) = -m_visc*                          &
             (D3_g(igas,:) - 3.0_WP*a2*D1_g(igas,:))                   &
           - ii*r_dens*alpha_nd*Re_l*Up_g(igas)*Id(igas,:)              &
           - ii*alpha_nd*Re_l*(F_grav + a2*S_st)                        &
             *D1_g(igas,:)/(Up_l(iliq) - Up_g(igas))
        B(nos+2, iliq_min:iliq_max) = -ii*alpha_nd*Re_l*D1_l(iliq,:)
        B(nos+2, igas_min:igas_max) =  ii*r_dens*alpha_nd*Re_l*D1_g(igas,:)

        ! --- Row scaling for numerical conditioning ---
        do i = 1, 2*nos
            tmp = max(maxval(abs( real(A(i,:)))), maxval(abs( real(B(i,:)))), &
                      maxval(abs(aimag(A(i,:)))), maxval(abs(aimag(B(i,:)))))
            if (tmp == 0.0_WP) tmp = 1.0_WP
            A(i,:) = A(i,:) / tmp
            B(i,:) = B(i,:) / tmp
        end do

        ! =============================================
        ! Reduce system: strip row/col 1 and 2*nos (phi=0 BC)
        ! =============================================
        allocate(Ain(2*(nos-1), 2*(nos-1)))
        allocate(Bin(2*(nos-1), 2*(nos-1)))
        Ain = A(2:2*nos-1, 2:2*nos-1)
        Bin = B(2:2*nos-1, 2:2*nos-1)

        ! =============================================
        ! Solve GEVP with LAPACK ZGGEV
        ! Use workspace query to get optimal workspace
        ! =============================================
        allocate(rwork(8*2*(nos-1)))
        allocate(eigval_a(2*(nos-1)))
        allocate(eigval_b(2*(nos-1)))
        allocate(eigvec_l(1, 2*(nos-1)))
        allocate(eigvec_r(2*(nos-1), 2*(nos-1)))

        ! Workspace query
        allocate(work(1))
        lwork = -1
        call ZGGEV('N', 'V', 2*(nos-1), Ain, 2*(nos-1), Bin, 2*(nos-1), &
                   eigval_a, eigval_b, eigvec_l, 1, eigvec_r, 2*(nos-1), &
                   work, lwork, rwork, ierr)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))

        ! Actual solve
        call ZGGEV('N', 'V', 2*(nos-1), Ain, 2*(nos-1), Bin, 2*(nos-1), &
                   eigval_a, eigval_b, eigvec_l, 1, eigvec_r, 2*(nos-1), &
                   work, lwork, rwork, ierr)
        if (ierr /= 0) then
            print *, "Error: LAPACK ZGGEV failed with info =", ierr
            c_eval = (0.0_WP, 0.0_WP)
            ! Cleanup and return
            deallocate(zos, yos_l, yos_g, Id, Tmat, Gos, Tp, Gtmp)
            deallocate(D1_g, D2_g, D3_g, D4_g, D1_l, D2_l, D3_l, D4_l)
            deallocate(Ubase_l, Up_l, Upp_l, Ubase_g, Up_g, Upp_g)
            deallocate(A, B, Ain, Bin, work, rwork, eigval_a, eigval_b)
            deallocate(eigvec_l, eigvec_r)
            return
        end if
        ! =============================================
        ! Find the most unstable mode
        ! Filter: c_i > 0 (unstable), c_i < 1 (physical bound),
        !         |c_r| < 2 (reject spurious Chebyshev modes)
        ! In ND space, physical c_r ∈ [-U_l/U_g, 1]
        ! =============================================
        mode_index = -1
        omega = -ii*huge(1.0_WP)

        do jj = 1, 2*(nos-1)
            if (abs(eigval_b(jj)) /= 0.0_WP) then
                omega_ = eigval_a(jj) / (eigval_b(jj) + epsilon(abs(eigval_b(jj))))
                if (abs(omega_) > 1.0e-10_WP) then
                    ! Physical mode filters (non-dimensional)
                    if (aimag(omega_) > aimag(omega) .and. &
                        aimag(omega_) < 1.0_WP .and. &
                        abs(real(omega_)) < 2.0_WP) then
                        mode_index = jj
                        omega = omega_
                    end if
                end if
            end if
        end do

        ! The GEVP eigenvalue is the ND complex phase speed c/U_g
        ! Convert back to dimensional: c = c_nd * U_g
        if (mode_index == -1) then
            print *, "orr_sommerfeld: Could not find a correct mode"
            c_eval = (0.0_WP, 0.0_WP)
        else
            c_eval = omega * Ug_inf_in
        end if

        ! Store eigenvalue (dimensional phase speed)
        c_phase = c_eval

        ! Store eigenvector split into liquid/gas parts
        if (allocated(phi_l_kh)) deallocate(phi_l_kh)
        if (allocated(phi_g_kh)) deallocate(phi_g_kh)
        allocate(phi_l_kh(nos), phi_g_kh(nos))
        phi_l_kh = (0.0_WP, 0.0_WP)
        phi_g_kh = (0.0_WP, 0.0_WP)
        if (mode_index > 0) then
            phi_l_kh(2:nos) = eigvec_r(1:nos-1, mode_index)
            phi_g_kh(1:nos-1) = eigvec_r(nos:2*(nos-1), mode_index)
        end if

        ! Store Tp and D1 matrices at module level (ND) for eval_kh
        if (allocated(Tp_kh)) deallocate(Tp_kh)
        if (allocated(D1_l_kh)) deallocate(D1_l_kh)
        if (allocated(D1_g_kh)) deallocate(D1_g_kh)
        allocate(Tp_kh(nos,nos), D1_l_kh(nos,nos), D1_g_kh(nos,nos))
        Tp_kh = Tp
        D1_l_kh = D1_l
        D1_g_kh = D1_g

        ! Store ND stretches for chebev bounds
        stretch_g_nd = stretch_g
        stretch_l_nd = stretch_l

        ! =============================================
        ! Precompute Chebyshev coefficients for eval_kh
        ! This avoids matmul and allocation on every call
        ! =============================================
        block
            complex(WP) :: dphi_tmp(nos)
            real(WP) :: wreal_tmp(nos), wimag_tmp(nos)

            ! Allocate coefficient arrays
            if (allocated(cheb_dphi_r_l)) deallocate(cheb_dphi_r_l)
            if (allocated(cheb_dphi_i_l)) deallocate(cheb_dphi_i_l)
            if (allocated(cheb_dphi_r_g)) deallocate(cheb_dphi_r_g)
            if (allocated(cheb_dphi_i_g)) deallocate(cheb_dphi_i_g)
            if (allocated(cheb_phi_r_l))  deallocate(cheb_phi_r_l)
            if (allocated(cheb_phi_i_l))  deallocate(cheb_phi_i_l)
            if (allocated(cheb_phi_r_g))  deallocate(cheb_phi_r_g)
            if (allocated(cheb_phi_i_g))  deallocate(cheb_phi_i_g)
            allocate(cheb_dphi_r_l(nos), cheb_dphi_i_l(nos))
            allocate(cheb_dphi_r_g(nos), cheb_dphi_i_g(nos))
            allocate(cheb_phi_r_l(nos),  cheb_phi_i_l(nos))
            allocate(cheb_phi_r_g(nos),  cheb_phi_i_g(nos))

            ! Liquid: D1*phi_l (for u-perturbation)
            dphi_tmp = matmul(D1_l, phi_l_kh)
            cheb_dphi_r_l = matmul(Tp, real(dphi_tmp))
            cheb_dphi_i_l = matmul(Tp, aimag(dphi_tmp))

            ! Liquid: phi_l (for v-perturbation)
            cheb_phi_r_l = matmul(Tp, real(phi_l_kh))
            cheb_phi_i_l = matmul(Tp, aimag(phi_l_kh))

            ! Gas: D1*phi_g
            dphi_tmp = matmul(D1_g, phi_g_kh)
            cheb_dphi_r_g = matmul(Tp, real(dphi_tmp))
            cheb_dphi_i_g = matmul(Tp, aimag(dphi_tmp))

            ! Gas: phi_g
            cheb_phi_r_g = matmul(Tp, real(phi_g_kh))
            cheb_phi_i_g = matmul(Tp, aimag(phi_g_kh))

            ! Interface values for level set (gas side at y_nd=0)
            phi_r_interface = chebev(0.0_WP, stretch_g, cheb_phi_r_g, nos, 0.0_WP)
            phi_i_interface = chebev(0.0_WP, stretch_g, cheb_phi_i_g, nos, 0.0_WP)

            ! Precompute level set coefficient
            block
                real(WP) :: c_r, c_i, denom
                c_r = real(c_phase) / Ug_inf
                c_i = aimag(c_phase) / Ug_inf
                denom = alpha_nd**2 * (c_r**2 + c_i**2)
                if (abs(denom) > 1.0e-30_WP) then
                    G_coeff1 = alpha_nd**2 / denom * c_i
                    G_coeff2 = alpha_nd**2 / denom * c_r
                else
                    G_coeff1 = 0.0_WP
                    G_coeff2 = 0.0_WP
                end if
            end block
        end block

        ! Store eigenvector for backward compat
        if (allocated(V_kh)) deallocate(V_kh)
        allocate(V_kh(2*nos))
        V_kh = zero
        if (mode_index > 0) then
            V_kh(2:2*nos-1) = eigvec_r(:, mode_index)
        end if

        ! Cleanup (keep Tp_kh, D1_kh, cheb_* stored at module level)
        deallocate(zos, yos_l, yos_g)
        deallocate(Id, Tmat, Gos, Tp, Gtmp)
        deallocate(D1_g, D2_g, D3_g, D4_g, D1_l, D2_l, D3_l, D4_l)
        deallocate(Ubase_l, Up_l, Upp_l, Ubase_g, Up_g, Upp_g)
        deallocate(A, B, Ain, Bin)
        deallocate(work, rwork, eigval_a, eigval_b, eigvec_l, eigvec_r)

    end subroutine setup_kh


    !> Evaluate KH perturbation at a physical point (x,y)
    !! Returns: u_out = base_flow + u_perturbation
    !!          v_out = v_perturbation
    !!          G_out = signed distance level set
    !! Uses precomputed Chebyshev coefficients — NO allocations, NO matmuls
    subroutine eval_kh(x_in, y_in, eps_in, u_out, v_out, G_out)
        implicit none
        real(WP), intent(in)  :: x_in, y_in, eps_in
        real(WP), intent(out) :: u_out, v_out, G_out

        integer :: nos
        real(WP) :: y_nd, alpha_dim
        real(WP) :: phi_r, phi_i, dphi_r, dphi_i
        real(WP) :: cos_ax, sin_ax

        nos = N_cheb + 1
        alpha_dim = alpha_wave
        y_nd = y_in / del_g

        ! Precompute trig (used for all quantities)
        cos_ax = cos(alpha_dim * x_in)
        sin_ax = sin(alpha_dim * x_in)

        if (y_in <= 0.0_WP) then
            if (y_nd < -stretch_l_nd) then
                dphi_r = 0.0_WP
                dphi_i = 0.0_WP
                phi_r  = 0.0_WP
                phi_i  = 0.0_WP
            else
                ! --- Liquid side ---
                dphi_r = chebev(-stretch_l_nd, 0.0_WP, cheb_dphi_r_l, nos, y_nd)
                dphi_i = chebev(-stretch_l_nd, 0.0_WP, cheb_dphi_i_l, nos, y_nd)
                phi_r  = chebev(-stretch_l_nd, 0.0_WP, cheb_phi_r_l,  nos, y_nd)
                phi_i  = chebev(-stretch_l_nd, 0.0_WP, cheb_phi_i_l,  nos, y_nd)
            end if
            ! Base flow (dimensional)
            u_out = Ul_inf * erf(y_in / del_l)
        else
            if (y_nd > stretch_g_nd) then
                dphi_r = 0.0_WP
                dphi_i = 0.0_WP
                phi_r  = 0.0_WP
                phi_i  = 0.0_WP
            else
                ! --- Gas side ---
                dphi_r = chebev(0.0_WP, stretch_g_nd, cheb_dphi_r_g, nos, y_nd)
                dphi_i = chebev(0.0_WP, stretch_g_nd, cheb_dphi_i_g, nos, y_nd)
                phi_r  = chebev(0.0_WP, stretch_g_nd, cheb_phi_r_g,  nos, y_nd)
                phi_i  = chebev(0.0_WP, stretch_g_nd, cheb_phi_i_g,  nos, y_nd)
            end if
            ! Base flow (dimensional)
            u_out = Ug_inf * erf(y_in / del_g)
        end if

        ! u-perturbation: u' = eps * (dphi_r cos(ax) - dphi_i sin(ax))
        u_out = u_out + eps_in * (dphi_r * cos_ax - dphi_i * sin_ax)

        ! v-perturbation: v' = eps * delta_g * alpha * (phi_i cos(ax) + phi_r sin(ax))
        ! From stream function: psi_dim = eps*Ug*delta_g*phi_nd*exp(i*alpha*x)
        ! v = -dpsi/dx => Re[v] = eps*Ug*delta_g*alpha*(phi_i*cos(ax) + phi_r*sin(ax))
        ! eps_in = factor*Ug already, so:
        v_out = eps_in * del_g * alpha_dim * &
                (phi_i * cos_ax + phi_r * sin_ax)

        ! Level set from paper Eq 71 (precomputed coefficients)
        ! G = y + eps*delta_g * [G_coeff1*(phi_r_int*cos + phi_i_int*sin)
        !                      - G_coeff2*(phi_r_int*cos - phi_i_int*sin)]
        G_out = y_in + eps_in * del_g * ( &
            G_coeff1 * (phi_r_interface * cos_ax + phi_i_interface * sin_ax) &
          - G_coeff2 * (phi_r_interface * cos_ax - phi_i_interface * sin_ax) )

    end subroutine eval_kh


    !> Clenshaw recurrence for evaluating a Chebyshev series
    !! a, b are the physical domain bounds, y is the evaluation point
    pure function chebev(a, b, c, n, y) result(val)
        implicit none
        real(WP), intent(in) :: a, b, y
        integer,  intent(in) :: n
        real(WP), intent(in) :: c(n)
        real(WP) :: val
        real(WP) :: d, dd, sv, y2, ym
        integer  :: j

        ! Map y to [-1,1]
        ym = (2.0_WP * y - a - b) / (b - a)
        y2 = 2.0_WP * ym

        ! Clenshaw recurrence (backward)
        d = 0.0_WP; dd = 0.0_WP
        do j = n, 2, -1
            sv = d
            d = y2 * d - dd + c(j)
            dd = sv
        end do
        val = ym * d - dd + c(1)
    end function chebev

end module orr_sommerfeld_kh
