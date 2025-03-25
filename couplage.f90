program couplage

    use mod_constantes
    use mod_algebre
    use mod_schema

    implicit none

    real(PR), dimension(3)  :: k, k_new
    real(PR)                :: t, t0, tf, dt, dx, T_0, T_init, T_Nx, khi, mv_bois, lambda_phalf, lambda_mhalf, eta, L
    real(PR), dimension(:), allocatable :: T_new, T_n, rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda
    integer                 :: type_bois, Nx, Nt, n, i
    character(len=20)           :: nc
    
    

    type_bois = 1
    tf = 1._PR
    t0 = 0._PR
    T_0 = 800._PR
    T_init = 300._PR
    T_Nx = 800._PR
    eta = 0._PR
    L = 0.1_PR


    Nx = 100
    dx = L/Nx
    

    khi = humidite(type_bois)
    mv_bois = densite_bois(type_bois)

    allocate(rho_b(0:Nx), rho_c(0:Nx), rho_g(0:Nx), rho_l(0:Nx), rho_v(0:Nx), T_n(0:Nx), T_new(0:Nx), rhoCp(0:Nx), lambda(0:Nx))

    rho_b = mv_bois*(1._PR-khi)
    rho_c = 0._PR
    rho_g = 0._PR
    rho_l = mv_bois*khi
    rho_v = 0._PR

    Nt = 1000
    dt = (tf-t0)/Nt

    T_n = T_init
    T_new = T_init
    T_n(0) = T_0
    T_n(Nx) = T_Nx
    T_new(0) = T_0
    T_new(Nx) = T_Nx

    do n = 0, Nt

        print *, t, T_new(5)
        

        if ( mod(n, 100) == 0) then
                
            write(nc,*) n
            open(unit = 100, file='data_temp/temperature_t'//trim(adjustl(nc))//'.dat')
        
        end if
        
        do i = 0, Nx

            rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)

            eta = rho_b(i)/(rho_b(i)+rho_c(i))
            lambda(i) = eta*0.105_PR + (1-eta)*(0.166_PR + 0.369_PR*khi)
        end do

        
        write(100,*), dx*0, T_new(0)
        do i = 1, Nx-1
            
            lambda_mhalf = (2 * lambda(i-1)* lambda(i))/(lambda(i+1) + lambda(i))
            lambda_phalf = (2 * lambda(i+1) * lambda(i))/(lambda(i+1) + lambda(i))
            T_new(i) = T_n(i) + dt/(dx**2 * rhoCp(i)) * &
            (lambda_phalf*T_n(i+1) - (lambda_phalf + lambda_mhalf) * T_n(i) + lambda_mhalf * T_n(i-1))  
            
            if ( mod(n, 100) == 0) then
                
                write(100,*), dx*i, T_new(i)
            
            end if
        end do

        write(100,*) dx*Nx, T_new(Nx)

        
        do i = 0, Nx

            call arrhenius(k, T_n(i))
            call arrhenius(k_new, T_new(i))
            call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k, k_new, dt)

        end do

        T_n = T_new
        t = t + dt

        close(100)
    end do 

    deallocate(rho_b, rho_c, rho_g, rho_l, rho_v, T_n, T_new, rhoCp, lambda)
    contains


        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius

        ! subroutine calcul_rhoCp(rhoCp)

        ! end subroutine calcul_rhoCp
end program