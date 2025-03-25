module mod_schema
    
    use mod_constantes

    implicit none

    contains

        subroutine Euler_step(rho_b, rho_c, rho_g, rho_l, rho_v, k, dt)

            real(PR), dimension(:), intent(inout)   :: k
            real(PR), intent(inout)                 :: rho_b, rho_c, rho_g, rho_l, rho_v
            real(PR), intent(in)                    :: dt


            rho_b = rho_b/(1._PR + dt * (k(1) + k(2)))
            rho_c = rho_c + dt * k(1) * rho_b
            rho_g = rho_g + dt * k(2) * rho_b
            rho_l = rho_l/(1._PR + dt * k(3))
            rho_v = rho_v + dt * k(3) * rho_l

        end subroutine Euler_step


        subroutine CK2_step(rho_b, rho_c, rho_g, rho_l, rho_v, k, k_new, dt)

            real(PR), dimension(:), intent(inout)   :: k, k_new
            real(PR), intent(inout)                 :: rho_b, rho_c, rho_g, rho_l, rho_v
            real(PR), intent(in)                    :: dt
            real(PR)                                :: rho_b_old, rho_l_old


            rho_b_old = rho_b
            rho_l_old = rho_l

            rho_b = rho_b * ((1._PR - dt/2._PR * (k(1) + k(2)))/(1._PR + dt/2._PR *(k_new(1) + k_new(2))))
            rho_c = rho_c + dt/2._PR * (k_new(1) * rho_b + k(1) * rho_b_old)
            rho_g = rho_g + dt/2._PR * (k_new(2) * rho_b + k(2) * rho_b_old)
            rho_l = rho_l * ((1._PR - dt/2._PR * k(3))/(1._PR + dt/2._PR * k_new(3)))
            rho_v = rho_v + dt/2._PR * (k_new(3) * rho_l + k(3) * rho_l_old)

        end subroutine CK2_step

end module mod_schema