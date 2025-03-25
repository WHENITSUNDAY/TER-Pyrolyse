program chimie

    use constantes
    use schema
    
    implicit none 

    real(PR), dimension(3)  :: k
    real(PR)                :: t, tf, dt, Temp, Tempf
    real(PR), dimension(5)  :: rho
    integer                 :: type_bois
 
    call read_parameters(t, tf, Temp, Tempf, type_bois)
    call initialize_rho(rho, type_bois)
    call arrhenius(k, Temp)

    dt = 0.1_PR
    print *, "dt :", dt
    print *, "k :", k

    ! call create_data_CK2(rho, k, t, tf, dt, Temp, Tempf)
    ! call create_data_Euler(rho, k, t, tf, dt, Temp, Tempf)
    ! call error_Euler(rho, k, t, tf, dt, Temp, Tempf)
    ! call error_CK2(rho, k, t, tf, dt, Temp, Tempf)

    contains

        subroutine create_data_Euler(rho, k, t, tf, dt, Temp, Tempf)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt, Temp, Tempf
            real(PR)                                :: t1, t2, alpha
            integer                                 :: N, i

            open(unit=101, file = "data/densite_Euler.dat")
            
            N = INT(tf/dt) !Calcul du nb d'opérations

            alpha = (Tempf-Temp)/tf !Le coefficient directeur de T(t)

            call cpu_time(t1)
            do i = 1, N

                write(101,*) rho, t, Temp
                Temp = Temp + dt*alpha

                call arrhenius(k, Temp)
                call Euler_step(rho, k, dt)
                t = t + dt
            
            end do

            call cpu_time(t2)
            print *, "Temps d'execution EI :", t2 - t1

            close(101)

        end subroutine create_data_Euler


        subroutine create_data_CK2(rho, k, t, tf, dt, Temp, Tempf)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt, Temp, Tempf
            real(PR)                                :: t1, t2, alpha
            real(PR), dimension(3)                  :: k_new
            integer                                 :: N, i

            open(unit=101, file = "data/densite_CK2.dat")
            
            N = INT(tf/dt) !Calcul du nb d'opérations

            alpha = (Tempf-Temp)/tf !Le coefficient directeur de T(t)

            k_new = k
            call cpu_time(t1)
            do i = 1, N

                write(101,*) rho, t, Temp

                
                Temp = Temp + dt*alpha
                call arrhenius(k_new, Temp)

                call CK2_step(rho, k, k_new, dt)
                t = t + dt
                k = k_new
            
            end do

            call cpu_time(t2)
            print *, "Temps d'execution CK2 :", t2 - t1

            close(101)

        end subroutine create_data_CK2

        subroutine error_Euler(rho, k, t, tf, dt, Temp, Tempf)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt, Temp, Tempf
            real(PR)                                :: t1, t2, alpha, Temp10
            real(PR), dimension(3)                  :: k_new, k10_new, k10
            real(PR), dimension(5)                  :: rho_E1, rho_E10, rho_CK
            integer                                 :: N, i, j

            open(unit=101, file = "data/error_Euler.dat")
            
            N = INT(tf/dt) !Calcul du nb d'opérations

            alpha = (Tempf-Temp)/tf !Le coefficient directeur de T(t)

            k_new = k
            k10_new = k
            k10 = k

            rho_E1 = rho
            rho_E10 = rho
            rho_CK = rho

            Temp10 = Temp

            do i = 1, N

                write(101,*) abs(rho_E1 - rho_E10), t, Temp

                Temp = Temp + dt*alpha
                call arrhenius(k_new, Temp)

                do j = 1, 10

                    Temp10 = Temp10 + dt*alpha/10._PR

                    call arrhenius(k10_new, Temp10)
                    call Euler_step(rho_E10, k10_new, dt/10._PR)
                    
                end do
                
                call Euler_step(rho_E1, k_new, dt)       
                t = t + dt
            
            end do

            call cpu_time(t2)
            print *, "Temps d'execution Erreur EI :", t2 - t1

            close(101)

        end subroutine error_Euler

        subroutine error_CK2(rho, k, t, tf, dt, Temp, Tempf)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt, Temp, Tempf
            real(PR)                                :: t1, t2, alpha, Temp10
            real(PR), dimension(3)                  :: k_new, k10_new, k10
            real(PR), dimension(5)                  :: rho_CK, rho_CK10
            integer                                 :: N, i, j

            open(unit=102, file = "data/error_CK2.dat")
            
            N = INT(tf/dt) !Calcul du nb d'opérations

            alpha = (Tempf-Temp)/tf !Le coefficient directeur de T(t)

            k_new = k
            k10_new = k
            k10 = k

            rho_CK = rho
            rho_CK10 = rho

            Temp10 = Temp

            do i = 1, N

                write(102,*) abs(rho_CK - rho_CK10), t, Temp

                Temp = Temp + dt*alpha
                call arrhenius(k_new, Temp)

                do j = 1, 10

                    Temp10 = Temp10 + dt*alpha/10._PR

                    call arrhenius(k10_new, Temp10)
                    call CK2_step(rho_CK10, k10, k10_new, dt/10._PR)
                    

                    k10 = k10_new
                end do

                call CK2_step(rho_CK, k, k_new, dt)
               
                t = t + dt
                k = k_new
            
            end do

            call cpu_time(t2)
            print *, "Temps d'execution Erreur CK2 :", t2 - t1

            close(102)

        end subroutine error_CK2


        subroutine read_parameters(t, tf, Temp, Tempf, type_bois)
            real(PR), intent(out)           :: t, tf, Temp, Tempf
            integer, intent(out)            :: type_bois 


            open(unit=100, file="parametres.dat", action="read")

            read(100,*) t
            read(100,*) tf
            read(100,*) Temp
            read(100,*) Tempf
            read(100,*) type_bois

            close(100)
        
        end subroutine read_parameters


        subroutine initialize_rho(rho, type_bois)

            real(PR), dimension(:), intent(out) :: rho
            integer, intent(in)                 :: type_bois

            rho = 0
            rho(1) = densite_bois(type_bois)*(1._PR-humidite(type_bois))
            rho(4) = densite_bois(type_bois)*humidite(type_bois)
  
        end subroutine initialize_rho


        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius
    
end program chimie
