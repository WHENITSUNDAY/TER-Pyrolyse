program couplage

    use mod_constantes
    use mod_algebre
    use mod_schema

    implicit none

    real(PR), dimension(3)  :: k, k_new
    real(PR)                :: t, t0, tf, dt, dx, T_0, T_init, T_Nx, khi, mv_bois, lambda_phalf, lambda_mhalf
    real(PR), dimension(:), allocatable :: T_new, T_n, rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp
    integer                 :: type_bois, Nx, Nt, n, i
    
    

    type_bois = 1
    tf = 1._PR
    t0 = 0._PR
    T_0 = 400._PR
    T_init = 300._PR
    T_Nx = 400._PR
    Nx = 10
    Nt = 10

    dx = 1._PR/Nx
    dt = 1._PR/Nt

    khi = humidite(type_bois)
    mv_bois = densite_bois(type_bois)

    allocate(rho_b(0:Nx), rho_c(0:Nx), rho_g(0:Nx), rho_l(0:Nx), rho_v(0:Nx), T_n(0:Nx), T_new(0:Nx), rhoCp(0:Nx), lambda_h(0:Nx))

    rho_b = mv_bois*(1._PR-khi)
    rho_c = 0
    rho_g = 0
    rho_l = mv_bois*khi
    rho_v = 0

    dt = 0.1_PR
    do n = 0, Nt

        do i = 0, Nx

            call arrhenius(k, T_n(i))
            call arrhenius(k_new, T_new(i))
            call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k, k_new, dt)
            
            rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)
            
        end do

        T_new(0) = T_0
        T_new(Nx) = T_Nx

        do i = 1, Nx-1
            call calcul_diffusivity(lambda_mhalf, lambda_phalf, rho_b, rho_c, i)
            T_new(i) = T_n(i) + dt/(dx**2 * rhoCp) * (lambda(i)*T_n(i+1) - (lambda_h + lambda_mhalf) * T_n(i) + lambda_mhalf * T_n(i-1))  

        end do

        t = t + dt
    end do 

    contains


        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius

        subroutine calcul_diffusivity(lambda_mhalf, lambda_phalf, rho_b, rho_c, i)

            real(PR), intent(out) :: lambda_mhalf, lambda_phalf 
            real(PR), dimension(:), intent(in) :: rho_b, rho_c
            real(PR)                :: eta, eta_m, eta_p, lambda, lambda_p, lambda_m
            integer                 :: i

            eta = rho_b(i)/(rho_b(i)+rho_c(i))
            eta_p = rho_b(i+1)/(rho_b(i+1)+rho_c(i+1))
            eta_m = rho_b(i-1)/(rho_b(i-1)+rho_c(i-1))

            lambda = eta*0.105_PR + (1-eta)*(0.166_PR + 0.369_PR*khi)
            lambda_p = eta_p*0.105_PR + (1-eta_p)*(0.166_PR + 0.369_PR*khi)
            lambda_m = eta_p*0.105_PR + (1-eta_m)*(0.166_PR + 0.369_PR*khi)

            lambda_mhalf = (2 * lambda_m * lambda)/(lambda_m + lambda)
            lambda_phalf = (2 * lambda_p * lambda)/(lambda_p + lambda)

        end subroutine calcul_diffusivity

        ! subroutine calcul_rhoCp(rhoCp)

        ! end subroutine calcul_rhoCp
end program