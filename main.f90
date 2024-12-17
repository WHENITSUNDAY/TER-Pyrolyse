program main

    use constantes
    use schema
    
    implicit none 

    integer :: N, i
    real(PR) :: tf, dt, t

    real(PR) :: Temp
    real(PR), dimension(5) :: rho
    real(PR), dimension(3) :: k
    character(len=10) :: bois

    Temp = 800
    rho = 0 
    k = (/  (A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)
    tf = 5
    dt = min(0.9_PR/(k(1)+k(2)), 0.9_PR/(k(3)))
    
    print *, dt
    N = INT(tf/dt) + 1
    t = 0  

    print *, k

    print *, "Sélectionnez l'essence de bois"
    print *, "Chêne - Bouleau - Tremble - Epicea - Pin blanc - Meleze"
    read *, bois
    
    call initialize_rho(rho, bois, densite_bois, humidite)

    open(unit=100, file = "densite.dat")

    do i = 1, N

        call RK4_step(rho, k, f, dt)
        print *, rho
        write(100,*) rho, dt*(i-1)

    end do
    close(100)


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
            character(len=10) :: bois 
            select case (bois)

                case ("Chêne")
                    rho(1) = densite_bois(1)*(1-humidite(1))
                    rho(4) = densite_bois(1)*humidite(1)
                case ("Bouleau")
                    rho(1) = densite_bois(2)*(1-humidite(2))
                    rho(4) = densite_bois(2)*humidite(2)
                case ("Tremble")
                    rho(1) = densite_bois(3)*(1-humidite(3))
                    rho(4) = densite_bois(3)*humidite(3)
                case ("Epicea")
                    rho(1) = densite_bois(4)*(1-humidite(4))
                    rho(4) = densite_bois(4)*humidite(4)
                case ("Pin blanc")
                    rho(1) = densite_bois(5)*(1-humidite(5))
                    rho(4) = densite_bois(5)*humidite(5)
                case ("Meleze")
                    rho(1) = densite_bois(6)*(1-humidite(6))
                    rho(4) = densite_bois(6)*humidite(6)
        
                case default
        
                    print *, "Veuillez entrer une essence de bois correcte" 
            
            end select
        
        end subroutine initialize_rho
    
end program main
