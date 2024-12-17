program main

    use constantes
    use schema
    
    implicit none 

    integer, parameter :: PR = 8
    
    integer :: N, i
    real(PR) :: tf, dt, t

    real(PR) :: T, R, lambda, D
    real(PR), dimension(5) :: rho
    real(PR), dimension(6) :: densite_bois, humidite
    real(PR), dimension(3) :: A, E, k
    character(len=10) :: bois

    !Paramètres annexes

    densite_bois = (/888, 740, 582, 469, 360, 469/)
    humidite = (/0.149_PR , 0.153_PR, 0.144_PR, 0.141_PR, 0.146_PR, 0.142_PR/) 

    T = 700
    R = 8.3144261_PR
    A = (/7.38e5_PR, 1.44e4_PR, 5.13e10_PR/)
    E = (/106.5e3_PR, 88.6e3_PR, 88e3_PR/)

    
    
    !Définiton des paramètres de la pyrolyse
    rho = 0 

    k = (/  (A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)
    !k = 0.5

    !Définition des variables de euler explicite
    tf = 10
    dt = min(0.9_PR/(k(1)+k(2)), 0.9_PR/(k(3)))
    
    print *, dt
    N = tf/dt
    t = 0  

    print *, k

    print *, "Sélectionnez l'essence de bois"
    print *, "Chêne - Bouleau - Tremble - Epicea - Pin blanc - Meleze"
    read *, bois
    
    call initialize_rho(rho, bois, densite_bois, humidite)

    open(unit=1,  status="replace", file = "densite.dat")


    close(1)

    contains


        function f(t, rho, k)
            real(PR), dimension(:), intent(in) :: rho, k
            real(PR) :: t
            real(PR), dimension(5), intent(out) :: f

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
