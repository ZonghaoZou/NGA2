module orr_sommerfeld_kh
    use precision,         only: WP
    implicit none
    ! --- State Variables from GEVP Solve ---
    integer :: N_cheb
    real(WP) :: alpha_wave, Ug_inf, Ul_inf, del_g, del_l, lambda_amp
    real(WP) :: Lg_dom, Ll_dom, dzg_dy, dzl_dy
    complex(WP), allocatable :: V_kh(:)
    real(WP), allocatable :: z_cheb(:), D1_cheb(:,:)

contains

    subroutine setup_kh(N, alpha_in, Ug_inf_in, Ul_inf_in, del_g_in, del_l_in, rho_g_in, rho_l_in, mu_g_in, mu_l_in, sigma_in, Lg_in, Ll_in, c_eval)
        use mathtools, only: Pi
        implicit none
        ! Inputs
        integer, intent(in) :: N                 
        real(WP), intent(in) :: alpha_in            
        real(WP), intent(in) :: Ug_inf_in, Ul_inf_in   
        real(WP), intent(in) :: del_g_in, del_l_in     
        real(WP), intent(in) :: rho_g_in, rho_l_in     
        real(WP), intent(in) :: mu_g_in, mu_l_in
        real(WP), intent(in) :: sigma_in            
        real(WP), intent(in) :: Lg_in, Ll_in               
        complex(WP), intent(out) :: c_eval
        
        ! Local Variables
        real(WP) :: m_visc, r_dens, Re_l, We_l
        real(WP), allocatable :: D2(:,:), D3(:,:), D4(:,:)
        real(WP), allocatable :: D1g(:,:), D2g(:,:), D3g(:,:), D4g(:,:)
        real(WP), allocatable :: D1l(:,:), D2l(:,:), D3l(:,:), D4l(:,:)
        complex(WP), allocatable :: A(:,:), B(:,:), W(:), Vl(:,:), Vr(:,:)
        complex(WP) :: c_max, phi_g
        integer :: info, i, j, lwork, max_idx,k
        complex(WP), allocatable :: work(:)
        real(WP), allocatable :: rwork(:)
        complex(WP), allocatable :: alpha_eig(:), beta_eig(:)
        
        real(WP) :: a2, a4
        real(WP), allocatable :: yy_g(:), Ug_arr(:), Ugg_arr(:)
        real(WP), allocatable :: yy_l(:), Ul_arr(:), Ull_arr(:)
        complex(WP) :: fac_g, fac_l, visc_fac_l, visc_fac_g
        real(WP) :: Ug_prime_int, Ul_prime_int, ST_fac,eye
        complex(WP), parameter :: im = (0.0_WP, 1.0_WP)
        real(WP) :: cr, ci, cr_min, cr_max
        ! Store parameters
        N_cheb = N
        ! alpha_wave = alpha
        ! Ug_inf = Ug_inf_in
        ! Ul_inf = Ul_inf_in
        ! del_g = del_g_in
        ! del_l = del_l_in
        ! Lg_dom = Lg
        ! Ll_dom = Ll
        ! lambda_amp = lambda
        

        ! 1. Dimensionless Groups (Calculated using pure physical inputs)
        m_visc = mu_g_in / mu_l_in
        r_dens = rho_g_in / rho_l_in
        Re_l   = (rho_l_in * Ug_inf_in * del_g_in) / mu_l_in
        
        ! Protect against zero surface tension dividing by zero
        if (sigma_in > 1.0e-12_WP) then
            We_l = (rho_l_in * (Ug_inf_in**2) * del_g_in) / sigma_in
        else
            We_l = 1.0e99_WP 
        end if

        ! 2. Unified Length Scales (Everything divided by gas thickness)
        Lg_dom = Lg_in / del_g_in
        Ll_dom = Ll_in / del_g_in    ! Corrected to del_g_in
        del_l  = del_l_in / del_g_in
        del_g  = 1.0_WP
        
        ! 3. Unified Velocity Scales (Everything divided by gas velocity)
        Ul_inf = Ul_inf_in / Ug_inf_in
        Ug_inf = 1.0_WP

        ! 4. Wavenumbers
        alpha_wave = alpha_in * del_g_in
        a2 = alpha_wave**2
        a4 = alpha_wave**4
        ! ! 1. Non-dimensional parameters
        ! m_visc = mu_g / mu_l
        ! r_dens = rho_g / rho_l
        ! Re_l = (rho_l * Ug_inf * del_g) / mu_l
        ! We_l = (rho_l * Ug_inf**2 * del_g) / sigma
        ! Lg_dom=Lg_dom/del_g_in
        ! Ll_dom=Ll_dom/del_g_in
        ! alpha_wave=alpha*del_g_in
        ! del_g=1.0_WP
        ! del_l=del_l/del_g_in
        ! Ug_inf = 1.0_WP
        ! Ul_inf=Ul_inf/Ug_inf_in
        ! a2 = alpha_wave**2
        ! a4 = alpha_wave**4
        
        ! 2. Setup Chebyshev Collocation
        if (allocated(z_cheb)) deallocate(z_cheb, D1_cheb, V_kh)
        allocate(z_cheb(0:N), D1_cheb(0:N,0:N), D2(0:N,0:N), D3(0:N,0:N), D4(0:N,0:N), V_kh(0:2*N+1))
        call chebyshev_differentiation(N, z_cheb, D1_cheb, D2, D3, D4)
        
        ! Map D matrices to physical gas and liquid domains
        allocate(D1g(0:N,0:N), D2g(0:N,0:N), D3g(0:N,0:N), D4g(0:N,0:N))
        allocate(D1l(0:N,0:N), D2l(0:N,0:N), D3l(0:N,0:N), D4l(0:N,0:N))
        
        dzg_dy = -2.0_WP / Lg_dom
        D1g = D1_cheb * dzg_dy
        D2g = D2 * (dzg_dy**2)
        D3g = D3 * (dzg_dy**3)
        D4g = D4 * (dzg_dy**4)
        
        dzl_dy = -2.0_WP / Ll_dom
        D1l = D1_cheb * dzl_dy
        D2l = D2 * (dzl_dy**2)
        D3l = D3 * (dzl_dy**3)
        D4l = D4 * (dzl_dy**4)
        
        ! Base flow arrays
        allocate(yy_g(0:N), Ug_arr(0:N), Ugg_arr(0:N))
        allocate(yy_l(0:N), Ul_arr(0:N), Ull_arr(0:N))
        do i = 0, N
            ! 1. Map Chebyshev coordinates to physical domains
            yy_g(i) = (1.0_WP - z_cheb(i)) * Lg_dom / 2.0_WP
            yy_l(i) = (-1.0_WP - z_cheb(i)) * Ll_dom / 2.0_WP
            
            ! 2. Gas phase base flow and its second derivative
            Ug_arr(i)  = Ug_inf * erf(yy_g(i) / del_g)
            Ugg_arr(i) = Ug_inf * (-4.0_WP * yy_g(i) / (del_g**3 * sqrt(Pi))) * exp(-(yy_g(i)/del_g)**2)
            
            ! 3. Liquid phase base flow and its second derivative
            Ul_arr(i)  = Ul_inf * erf(yy_l(i) / del_l)
            Ull_arr(i) = Ul_inf * (-4.0_WP * yy_l(i) / (del_l**3 * sqrt(Pi))) * exp(-(yy_l(i)/del_l)**2)
        end do
        ! Evaluate the analytical first derivatives exactly at the interface (y=0)
        Ug_prime_int = (2.0_WP * Ug_inf) / (del_g * sqrt(Pi))
        Ul_prime_int = (2.0_WP * Ul_inf) / (del_l * sqrt(Pi))

        ! 3. Construct the GEVP Matrices: A * X = c * B * X
        allocate(A(0:2*N+1, 0:2*N+1), B(0:2*N+1, 0:2*N+1))
        A=0.0_WP; B=0.0_WP

        fac_g = m_visc / (r_dens * im * alpha_wave * Re_l)
        fac_l = 1.0_WP / (im * alpha_wave * Re_l)

        k = 0

        ! --- Gas Outer Boundaries (y = Lg, index N) ---
        A(k, 0:N-1) = 0.0_WP; A(k, N) = 1.0_WP; B(k, :) = 0.0_WP; k = k + 1 ! phi_g = 0
        A(k, 0:N) = D1g(N, :); B(k, :) = 0.0_WP; k = k + 1              ! phi_g' = 0

        ! --- Gas Interior Equations (Indices 2 to N-2) ---
        do i = 2, N-2
            do j = 0, N
                ! Scalar terms only apply to the diagonal
                if (i == j) then
                    eye = 1.0_WP
                else
                    eye = 0.0_WP
                end if
                A(k, j) = Ug_arr(i)*(D2g(i,j) - a2*eye) - Ugg_arr(i)*eye - &
                          fac_g*(D4g(i,j) - 2.0_WP*a2*D2g(i,j) + a4*eye)
                B(k, j) = D2g(i,j) - a2*eye
            end do
            k = k + 1
        end do

        ! --- Liquid Interior Equations (Indices 2 to N-2) ---
        do i = 2, N-2
            do j = 0, N
                if (i == j) then
                    eye = 1.0_WP
                else
                    eye = 0.0_WP
                end if
                A(k, j+N+1) = Ul_arr(i)*(D2l(i,j) - a2*eye) - Ull_arr(i)*eye - &
                              fac_l*(D4l(i,j) - 2.0_WP*a2*D2l(i,j) + a4*eye)
                B(k, j+N+1) = D2l(i,j) - a2*eye
            end do
            k = k + 1
        end do

        ! --- Liquid Outer Boundaries (y = -Ll, index 0 for liquid, mapped to N+1 global) ---
        A(k, N+2:2*N+1) = 0.0_WP; A(k, N+1) = 1.0_WP; B(k, :) = 0.0_WP; k = k + 1 ! phi_l = 0
        A(k, N+1:2*N+1) = D1l(0, :); B(k, :) = 0.0_WP; k = k + 1                ! phi_l' = 0

        ! --- Interface Coupling Conditions (y=0 -> Gas index 0, Liquid index N) ---

        ! 1. Kinematic Condition (phi_g = phi_l)
        A(k, 0) = 1.0_WP; A(k, 2*N+1) = -1.0_WP; B(k, :) = 0.0_WP; k = k + 1

        ! 2. Tangential Velocity (Eq 67b: Ug'*phi_g - Ul'*phi_l = c*(-phi_g' + phi_l'))
        A(k, :) = 0.0_WP; B(k, :) = 0.0_WP   ! Clear the row first
        do j = 0, N
            ! A matrix gets the base flow derivative terms (only at specific interface nodes)
            if (j == 0) A(k, j)       = Ug_prime_int
            if (j == N) A(k, j+N+1)   = -Ul_prime_int
            
            ! B matrix gets the spatial derivative terms with the correct inverted signs
            B(k, j)       = -D1g(0, j)
            B(k, j+N+1)   =  D1l(N, j)
        end do
        k = k + 1

        ! 3. Normal Stress (Linearized Eq 68 with corrected sign)
        ! Note: The text notes this linearization is only valid when Ul_prime != Ug_prime
        if (abs(Ul_prime_int - Ug_prime_int) > 1.0e-12_WP .and. abs(We_l) > 1.0e-12_WP) then
            ST_fac = -(a2) / (r_dens * We_l * (Ul_prime_int - Ug_prime_int))
        else
            ST_fac = 0.0_WP
        end if
        visc_fac_l = 1.0_WP / (im * alpha_wave * r_dens * Re_l)
        visc_fac_g = m_visc / (im * alpha_wave * r_dens * Re_l)

        A(k, :) = 0.0_WP; B(k, :) = 0.0_WP   ! Clear the row

        do j = 0, N
            ! --- Gas Terms (Columns 0 to N) ---
            ! Derivative components
            A(k, j) = -ST_fac * D1g(0, j) & 
                      + visc_fac_g * (D3g(0, j) - 3.0_WP * a2 * D1g(0, j))
            
            ! Scalar component (Ug_prime * phi_g) only exists at the node j=0
            if (j == 0) A(k, j) = A(k, j) + Ug_prime_int
            
            ! B matrix gas component
            B(k, j) = -D1g(0, j)

            ! --- Liquid Terms (Columns N+1 to 2N+1) ---
            ! Derivative components
            A(k, j+N+1) = ST_fac * D1l(N, j) & 
                        - visc_fac_l * (D3l(N, j) - 3.0_WP * a2 * D1l(N, j))
            
            ! Scalar component (-1/r * Ul_prime * phi_l) only exists at node j=N
            if (j == N) A(k, j+N+1) = A(k, j+N+1) - (1.0_WP / r_dens) * Ul_prime_int
            
            ! B matrix liquid component
            B(k, j+N+1) = (1.0_WP / r_dens) * D1l(N, j)
        end do
        k = k + 1

        ! 4. Tangential Stress (Eq 67d rearranged)
        A(k, :) = 0.0_WP; B(k, :) = 0.0_WP
        
        do j = 0, N
            ! --- Gas Terms (Columns 0 to N) ---
            ! m * phi_g''
            A(k, j) = m_visc * D2g(0, j)
            
            ! m * a2 * phi_g (Scalar term, only at node j=0)
            if (j == 0) A(k, j) = A(k, j) + m_visc * a2
            
            ! --- Liquid Terms (Columns N+1 to 2N+1) ---
            ! - phi_l''
            A(k, j+N+1) = -D2l(N, j)
            
            ! - a2 * phi_l (Scalar term, only at node j=N)
            if (j == N) A(k, j+N+1) = A(k, j+N+1) - a2
        end do
        k = k + 1
        

        ! 4. Solve GEVP using LAPACK (ZGGEV)
        allocate(alpha_eig(0:2*N+1), beta_eig(0:2*N+1), Vl(1,1), Vr(0:2*N+1, 0:2*N+1))
        
        lwork = Maximum(1, 4*(2*N+2))
        allocate(work(lwork), rwork(8*(2*N+2)))
        
        ! Note: LDVL is now perfectly matched to Vl(1,1)
        call zggev('N', 'V', 2*N+2, A, 2*N+2, B, 2*N+2, alpha_eig, beta_eig, &
                   Vl, 1, Vr, 2*N+2, work, lwork, rwork, info)
                   
        if (info /= 0) then
            print *, "Error: LAPACK ZGGEV failed with info =", info
            return
        end if
        
        ! 5. Find the Most Unstable Physical Mode
        allocate(W(0:2*N+1))
        W = (0.0_WP, 0.0_WP)
        
        max_idx = 0
        c_max = (0.0_WP, -1.0e10_WP)
        
        ! Define Howard's Semicircle Bounds based on non-dimensional velocities
        ! We add a 10% buffer (0.1) just to account for slight numerical diffusion/viscous slip
        ! Corrected bounds
        ! Replace your wide cr_min / cr_max bounds with this:
        ! The KH wave travels near the interface velocity (0.0)
        cr_min = -0.3_WP
        cr_max =  0.3_WP
        
        do i = 0, 2*N+1
            if (abs(beta_eig(i)) > 1.0e-12_WP) then
                W(i) = alpha_eig(i) / beta_eig(i)
                cr = real(W(i), WP)
                ci = aimag(W(i))
                
                ! The tighter filter isolates the interfacial mode
                if (cr >= cr_min .and. cr <= cr_max .and. ci < 10.0_WP) then
                    if (ci > aimag(c_max)) then
                        c_max = W(i)
                        max_idx = i
                    end if
                end if
            end if
        end do
        ! Normalize eigenvector so phi_g(0) is exactly 1.0 + 0.0i
        phi_g = Vr(0, max_idx)
        if (abs(phi_g) > 1.0e-10_WP) then
            Vr(:, max_idx) = Vr(:, max_idx) / phi_g
        end if

        ! Save the most unstable mode to the module-level state variable
        V_kh(:) = Vr(:, max_idx)
        
        ! Output the c_evlautae
        c_eval = c_max * Ug_inf_in

        deallocate(D2, D3, D4, A, B, W, Vl, Vr, work, rwork)
        deallocate(D1g, D2g, D3g, D4g, D1l, D2l, D3l, D4l)
        deallocate(alpha_eig, beta_eig, yy_g, Ug_arr, Ugg_arr, yy_l, Ul_arr, Ull_arr)
        
    end subroutine setup_kh

    ! subroutine eval_kh(x, y, u_init, v_init, G_init)
    !     real(WP), intent(in) :: x, y
    !     real(WP), intent(out) :: u_init, v_init, G_init
        
    !     real(WP) :: zg, zl, Ug_base, Ul_base
    !     complex(WP) :: phi, phi_prime, phi_g
    !     real(WP) :: a_x
        
    !     a_x = alpha_wave * x
        
    !     ! Evaluate interface perturbation at y=0 (z=1)
    !     call eval_cheb(N_cheb, V_kh(0:N_cheb), 1.0_WP, D1_cheb, phi, phi_prime)
    !     phi_g = phi
    !     G_init = y - lambda_amp * (real(phi)*cos(a_x) - aimag(phi)*sin(a_x))
        
    !     if (y > 0.0_WP) then
    !         if (y > Lg_dom) then
    !            zg = -1.0_WP
    !         else
    !            zg = -2.0_WP * y / Lg_dom + 1.0_WP
    !         end if
    !         call eval_cheb(N_cheb, V_kh(0:N_cheb), zg, D1_cheb, phi, phi_prime)
    !         phi_prime = phi_prime * dzg_dy
    !         Ug_base = Ug_inf * erf(y / del_g)
            
    !         u_init = Ug_base + lambda_amp * (real(phi_prime) * cos(a_x) - aimag(phi_prime) * sin(a_x))
    !         v_init = lambda_amp * alpha_wave * (real(phi) * sin(a_x) + aimag(phi) * cos(a_x))
            
    !     else
    !         if (y < -Ll_dom) then
    !            zl = 1.0_WP
    !         else
    !            zl = -2.0_WP * y / Ll_dom - 1.0_WP
    !         end if
    !         call eval_cheb(N_cheb, V_kh(N_cheb+1:2*N_cheb+1), zl, D1_cheb, phi, phi_prime)
    !         phi_prime = phi_prime * dzl_dy
    !         Ul_base = Ul_inf * erf(y / del_l)
            
    !         u_init = Ul_base + lambda_amp * (real(phi_prime) * cos(a_x) - aimag(phi_prime) * sin(a_x))
    !         v_init = lambda_amp * alpha_wave * (real(phi) * sin(a_x) + aimag(phi) * cos(a_x))
    !     end if
        
    ! end subroutine eval_kh

    ! ! Chebyshev evaluation
    ! subroutine eval_cheb(N, coeffs, z_val, D1, phi, phi_prime)
    !     integer, intent(in) :: N
    !     complex(WP), intent(in) :: coeffs(0:N)
    !     real(WP), intent(in) :: z_val
    !     real(WP), intent(in) :: D1(0:N,0:N)
    !     complex(WP), intent(out) :: phi, phi_prime
        
    !     integer :: j
    !     real(WP) :: x_j, w_j, diff, num_r, den_r
    !     complex(WP) :: p_num, p_den, phi_p_num, phi_p_den
        
    !     ! Check if we are exactly on a node to avoid division by zero
    !     do j = 0, N
    !         x_j = cos(dble(j)*acos(-1.0_WP)/dble(N))
    !         if (abs(z_val - x_j) < 1.0d-12) then
    !             phi = coeffs(j)
    !             phi_prime = sum(D1(j,:) * coeffs)
    !             return
    !         end if
    !     end do

    !     ! Barycentric interpolation for phi
    !     p_num = (0.0_WP, 0.0_WP)
    !     p_den = (0.0_WP, 0.0_WP)
        
    !     do j = 0, N
    !         x_j = cos(dble(j)*acos(-1.0_WP)/dble(N))
    !         w_j = (-1.0_WP)**j
    !         if (j == 0 .or. j == N) w_j = 0.5_WP * w_j
            
    !         diff = z_val - x_j
    !         p_num = p_num + (w_j / diff) * coeffs(j)
    !         p_den = p_den + (w_j / diff)
    !     end do
    !     phi = p_num / p_den
        
    !     ! For phi_prime, use the exact spectral derivative coefficients
    !     ! Evaluate the barycentric interpolation of the *derivative* field
    !     phi_p_num = (0.0_WP, 0.0_WP)
    !     phi_p_den = (0.0_WP, 0.0_WP)
        
    !     do j = 0, N
    !         x_j = cos(dble(j)*acos(-1.0_WP)/dble(N))
    !         w_j = (-1.0_WP)**j
    !         if (j == 0 .or. j == N) w_j = 0.5_WP * w_j
            
    !         diff = z_val - x_j
    !         ! Multiply D1 row j with coeffs to get the derivative at node j
    !         phi_p_num = phi_p_num + (w_j / diff) * sum(D1(j,:) * coeffs)
    !         phi_p_den = phi_p_den + (w_j / diff)
    !     end do
    !     phi_prime = phi_p_num / phi_p_den

    ! end subroutine eval_cheb

        ! Chebyshev D matrices (Stable Recursion Algorithm)
    subroutine chebyshev_differentiation(N, z, D1, D2, D3, D4)
        use mathtools, only: Pi
        implicit none
        integer, intent(in) :: N
        real(WP), intent(out) :: z(0:N), D1(0:N,0:N), D2(0:N,0:N), D3(0:N,0:N), D4(0:N,0:N)
        integer :: i, j, L
        real(WP) :: c_i, c_j, sum_diag
        real(WP) :: D_prev(0:N,0:N), D_curr(0:N,0:N)
        
        do i = 0, N
            z(i) = cos(real(i, kind=WP)*Pi/real(N, kind=WP))
        end do
        
        D1 = 0.0_WP
        
        ! 1. Compute D1 exactly as before
        do i = 0, N
            c_i = 1.0_WP; if (i == 0 .or. i == N) c_i = 2.0_WP
            do j = 0, N
                c_j = 1.0_WP; if (j == 0 .or. j == N) c_j = 2.0_WP
                if (i /= j) then
                    D1(i,j) = (c_i/c_j) * ((-1.0_WP)**(i+j)) / (z(i)-z(j))
                end if
            end do
        end do
        
        do i = 1, N-1
            D1(i,i) = -z(i) / (2.0_WP*(1.0_WP - z(i)**2))
        end do
        
        D1(0,0) = real(2*N**2 + 1, kind=WP) / 6.0_WP
        D1(N,N) = -D1(0,0)
        
        ! 2. Compute D2, D3, D4 using the Stable Recursion
        D_prev = D1
        do L = 2, 4
            D_curr = 0.0_WP
            do i = 0, N
                c_i = 1.0_WP; if (i == 0 .or. i == N) c_i = 2.0_WP
                sum_diag = 0.0_WP
                do j = 0, N
                    if (i /= j) then
                        c_j = 1.0_WP; if (j == 0 .or. j == N) c_j = 2.0_WP
                        D_curr(i,j) = (real(L, kind=WP) / (z(i)-z(j))) * &
                                      ( (c_i/c_j) * ((-1.0_WP)**(i+j)) * D_prev(i,i) - D_prev(i,j) )
                        sum_diag = sum_diag + D_curr(i,j)
                    end if
                end do
                ! The diagonal is strictly the negative sum of the off-diagonals
                D_curr(i,i) = -sum_diag
            end do
            
            ! Truncate numerical noise on symmetries (optional but very clean)
            do i = 0, N
               do j = 0, N
                  if (abs(D_curr(i,j)) < 1.0e-14_WP) D_curr(i,j) = 0.0_WP
               end do
            end do
            
            ! Store the matrix
            if (L == 2) D2 = D_curr
            if (L == 3) D3 = D_curr
            if (L == 4) D4 = D_curr
            D_prev = D_curr
        end do
        
    end subroutine chebyshev_differentiation


    
    integer function Maximum(a, b)
        integer, intent(in) :: a, b
        if (a > b) then
            Maximum = a
        else
            Maximum = b
        end if
    end function Maximum

end module orr_sommerfeld_kh
