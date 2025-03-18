program couplage

    use mod_constantes
    use mod_algebre
    use mod_schema

    implicit none

    real(PR), dimension(3)  :: k
    real(PR)                :: t0, tf, dt, T_0, T_i, T_Nx, khi, rhob
    real(PR), dimension(:), allocatable :: T_new, T_old, rhoCp
    real(PR), dimension(5)  :: rho
    integer                 :: type_bois, Nx, Nt
 
    type_bois = 1
    tf = 1._PR
    t0 = 0._PR
    T_0 = 400._PR
    T_i = 300._PR
    T_Nx = 400._PR
    Nx = 10
    Nt = 10

    dx = 1._PR/Nx
    dt = 1._PR/Nt

    khi = humidite(type_bois)
    rhob = densite_bois(type_bois)

    rho = 0
            rho(1) = rhob*(1._PR-khi)
            rho(4) = rhob*khi
    dt = 0.1_PR

    allocate(rhoCp(0:Nx), T_old(0:Nx), T_new(0:Nx))
    call arrhenius(k, Temp)


        do i = 1, Nx

            call arrhenius(k_new, T_old(i))

            call CK2_step(rho, k, k_new, dt)
            t = t + dt
            k = k_new
            
            
        end do

        do i = 1, Nx

        end do

    contains

        subroutine initialize_rho(rho, type_bois)

            real(PR), dimension(:), intent(out) :: rho
            integer, intent(in)                 :: type_bois

            

        end subroutine initialize_rho


        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius
end program