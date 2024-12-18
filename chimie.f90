program chimie

    use constantes
    use schema
    
    implicit none 

    real(PR), dimension(3)  :: k
    real(PR)                :: t, tf, dt, Temp
    real(PR), dimension(5)  :: rho
    integer                 :: type_bois
 
    call read_parameters(t, tf, Temp, type_bois, "parametres.dat")
    call initialize_rho(rho, type_bois)
    call initialize_arrhenius(k, Temp)

    dt = min(0.9_PR/(k(1)+k(2)), 0.9_PR/(k(3)))
    
    call create_data_RK4(rho, k, t, tf, dt, "densite_RK4.dat")

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

        
        subroutine create_data_RK4(rho, k, t, tf, dt, file_path)

            real(PR), dimension(:), intent(inout)   :: rho, k
            real(PR), intent(inout)                 :: t, tf, dt
            character(len=*), intent(in)            :: file_path
            
            open(unit=101, file = file_path)

            do while (t < tf)

                write(101,*) rho, t
                print *, rho
                call RK4_step(rho, k, f, dt)
                t = t + dt
            
            end do

            close(101)

        end subroutine create_data_RK4


        subroutine read_parameters(t, tf, Temp, type_bois, file_path)
            real(PR), intent(out)           :: t, tf, Temp
            integer, intent(out)            :: type_bois 
            character(len=*), intent(in)    :: file_path


            open(unit=100, file=file_path, action="read")

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
