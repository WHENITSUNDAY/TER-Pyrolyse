program chimie

    use constantes
    use schema
    
    implicit none 

    real(PR), dimension(3)  :: k
    real(PR)                :: t, tf, dt, Temp
    real(PR), dimension(5)  :: rho
    integer                 :: type_bois
 
    call read_parameters(t, tf, Temp, type_bois)
    call initialize_rho(rho, type_bois)
    call initialize_arrhenius(k, Temp)

<<<<<<< HEAD
    dt = min(0.9_PR/(k(1)+k(2)), 0.9_PR/(k(3)))
    print *, dt

    !call create_data(rho, k, t, tf, dt, Euler_step, "densite_Euler.dat")
    call create_data(rho, k, t, tf, dt, Heun_step, "densite_Heun.dat")
    !call create_data(rho, k, t, tf, dt, RK4_step, "densite_RK4.dat")
=======
    dt = 0.009_PR/(k(1)+k(2))
    print *, "dt :", dt
    print *, "k :", k
    !call create_data_RK4(rho, k, t, tf, dt, "densite_RK4.dat")
    call create_data_Heun(rho, k, t, tf, dt, "densite_Heun.dat")
>>>>>>> 73beb1fde9f1fc190369466cf732721abdc60403

    contains

        function f(rho, k)
            real(PR), dimension(:), intent(in)  :: rho, k
            real(PR), dimension(5)              :: f

            f(1) = -(k(1)+k(2))*rho(1) 
            f(2) = k(1)*rho(1) 
            f(3) = k(2)*rho(1) 
            f(4) = -k(3)*rho(4) 
            f(5) = k(3)*rho(4) 

        end function f


        subroutine create_data(rho, k, t, tf, dt, schema, file_path)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt
            real(PR)                                :: t1, t2
            real(PR), dimension(3)                  :: k_new
            character(len=*), intent(in)            :: file_path
            integer                                 :: N, i
            
            interface

                subroutine schema(X, Y, f, dt )
                    real(8), dimension(:), intent(inout)   :: X, Y
                    real(8), intent(in)                    :: dt
                    real(8), dimension(5)                  :: P, C

                    interface

                        function f(X, Y)
                            real(8), dimension(:), intent(in)   :: X, Y
                            real(8), dimension(5)               :: f
                        end function f 
                    
                    end interface

                end subroutine
            
            end interface

            open(unit=101, file = file_path)
            
            N = INT(tf/dt)
            print *, N
            call cpu_time(t1)
            do i = 1, N

                call initialize_arrhenius(k_new, Temp)
                write(101,*) rho, t, T
                !print *, rho
<<<<<<< HEAD
                call schema(rho, k, f, dt)
=======
                call Heun_step(rho, k, k_new, f, dt)
>>>>>>> 73beb1fde9f1fc190369466cf732721abdc60403
                t = t + dt
            
            end do

            call cpu_time(t2)
            print *, t2 - t1

            close(101)

        end subroutine create_data


        subroutine read_parameters(t, tf, Temp, type_bois)
            real(PR), intent(out)           :: t, tf, Temp
            integer, intent(out)            :: type_bois 


            open(unit=100, file="parametres.dat", action="read")

            read(100,*) t
            read(100,*) tf
            read(100,*) Temp
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


        subroutine initialize_arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine initialize_arrhenius
    
end program chimie
