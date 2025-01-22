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

    dt = 0.1 
    print *, "dt :", dt
    print *, "k :", k

    !call create_data_CK2(rho, k, t, tf, dt, "densite_CK2.dat")
    call create_data_Euler(rho, k, t, tf, dt, Temp, Tempf)

    contains

        subroutine create_data_Euler(rho, k, t, tf, dt, Temp, Tempf)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt, Temp, Tempf
            real(PR)                                :: t1, t2
            integer                                 :: N, i
            real(PR)                                :: alpha

            open(unit=101, file = "densite_Euler.dat")
            
            N = INT(tf/dt) !Calcul du nb d'opérations

            alpha = (Tempf-Temp)/tf !Le coefficient directeur de T(t)

            call cpu_time(t1)
            do i = 1, N

                write(101,*) rho, t, Temp
                print *, rho, t, Temp

                Temp = Temp + dt*alpha
                call arrhenius(k, Temp)
                call Euler_step(rho, k, dt)
                t = t + dt
            
            end do

            call cpu_time(t2)
            print *, t2 - t1

            close(101)

        end subroutine create_data_Euler


        subroutine create_data_CK2(rho, k, t, tf, dt, file_path)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt
            real(PR)                                :: t1, t2
            real(PR), dimension(3)                  :: k_new
            character(len=*), intent(in)            :: file_path
            integer                                 :: N, i

            open(unit=101, file = file_path)
            
            k_new = k
            N = INT(tf/dt) !Calcul du nb d'opérations
            print *, N
            call cpu_time(t1)
            do i = 1, N

                Temp = Temp + dt
                call arrhenius(k_new, Temp)
                write(101,*) rho, t, Temp
                print *, rho, t, Temp
                call CK2_step(rho, k, k_new, dt)
                t = t + dt
                k = k_new
            
            end do

            call cpu_time(t2)
            print *, t2 - t1

            close(101)

        end subroutine create_data_CK2



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
            rho(1) = densite_bois(type_bois)*(1-humidite(type_bois))
            rho(4) = densite_bois(type_bois)*humidite(type_bois)
  
        end subroutine initialize_rho


        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius
    
end program chimie
