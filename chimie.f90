program chimie

    use constantes
    use schema
    
    implicit none 

    real(PR), dimension(3) :: k
    real(PR) :: tf, dt, t, Temp
    real(PR), dimension(5) :: rho
    integer :: N, i, bois

    Temp = 800
    k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)
    tf = 5
    dt = min(0.9_PR/(k(1)+k(2)), 0.9_PR/(k(3)))
    
    print *, dt
    N = INT(tf/dt) + 1
    t = 0  

    open(unit=100, file = "parametres.dat", action="read")

    close(100)
    
    call initialize_rho(rho, bois)

    open(unit=101, file = "densite.dat")

    do i = 1, N

        call RK4_step(rho, k, f, dt)
        print *, rho
        write(101,*) rho, dt*(i-1)

    end do
    close(101)


    contains


        function f(rho, k)
            real(PR), dimension(:), intent(in) :: rho, k
            real(PR), dimension(5) :: f

            f(1) = -(k(1)+k(2))*rho(1) 
            f(2) = k(1)*rho(1) 
            f(3) = k(2)*rho(1) 
            f(4) = -k(3)*rho(4) 
            f(5) = k(3)*rho(4) 

        end function f


        subroutine initialize_rho(rho, bois, densite_bois, humidite)

            real(PR), dimension(:) :: rho, densite_bois, humidite
            integer :: type_bois

            rho = 0
            rho(1) = densite_bois(type_bois)*(1-humidite(type_bois))
            rho(4) = densite_bois(type_bois)*humidite(type_bois)
  
        end subroutine initialize_rho
    
end program chimie
