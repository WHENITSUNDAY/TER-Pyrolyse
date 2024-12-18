module schema
    
    use constantes

    implicit none

    contains

        subroutine RK4_step(X, Y, f, dt)

            real(PR), dimension(:), intent(inout)   :: X, Y
            real(PR), intent(in)                    :: dt
            real(PR), dimension(5)                  :: rk1, rk2, rk3, rk4 !Coefficient de Runge Kutta associé au vecteur f

            interface

                function f(X, Y)
                    real(8), dimension(:), intent(in)   :: X, Y !Il ne faut pas confondre k d'Arrhenius et les k de RK4, ici on appelle Y les coefficients chimiques
                    real(8), dimension(5)               :: f
                end function f 
            
            end interface

            rk1 = f(X, Y)
            rk2 = f(X + dt*rk1/2, Y)
            rk3 = f(X + dt*rk2/2, Y)
            rk4 = f(X + dt*rk3, Y)

            X = X + dt*(rk1 + 2 * rk2 + 2 * rk3 + rk4)/6

        end subroutine RK4_step


end module schema