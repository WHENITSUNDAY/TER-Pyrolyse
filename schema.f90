module schema
    
    use constantes

    implicit none

    contains

        subroutine Euler_step(rho, k_new, dt)

            real(PR), dimension(:), intent(inout)   :: rho, k_new
            real(PR), intent(in)                    :: dt


            rho(1) = rho(1)/(1 + dt * (k_new(1) + k_new(2)))
            rho(2) = rho(2) + dt * k_new(1) * rho(1)
            rho(3) = rho(3) + dt * k_new(2) * rho(1)
            rho(4) = rho(4)/(1 + dt * k_new(3))
            rho(5) = rho(5) + dt * k_new(3) * rho(4)

        end subroutine Euler_step


        subroutine CK2_step(rho, k, k_new, dt)

            real(PR), dimension(:), intent(inout)   :: rho, k, k_new
            real(PR), intent(in)                    :: dt
            real(PR)                                :: rhob, rhol


            rhob = rho(1)
            rhol = rho(4)
            rho(1) = rho(1) * ((2 - dt * (k(1) + k(2)))/(2 + dt*(k_new(1) + k_new(2))))
            rho(2) = rho(2) + dt/2 * (k_new(1) * rho(1) + k(1) * rhob)
            rho(3) = rho(3) + dt/2 * (k_new(2) * rho(1) + k(2) * rhob)
            rho(4) = rho(4) * ((1 - dt * k(3)/(1 + dt* k_new(3))))
            rho(5) = rho(5) + dt/2 * (k_new(3) * rho(4) + k(3) * rhol)

        end subroutine CK2_step

end module schema