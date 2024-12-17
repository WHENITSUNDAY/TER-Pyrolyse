module schema
    
    use constantes

    implicit none

    contains

        subroutine RK4_step(X, f, t, dt)

            real(PR), dimension(:), intent(inout) :: X
            real(PR), intent(in) :: dt
            integer :: i
            real(PR), dimension(4,5) :: rk !Coefficient de Runge Kutta associé au vecteur f

            interface

                function f(t ,X, Y)
                    real(PR) :: t
                    real(PR), dimension(:), intent(in) :: X, Y !Il ne faut pas confondre k d'Arrhenius et les k de RK4, ici on appelle Y les coefficients chimiques
                end function f 
            
            end interface

            do i = 1, size(X)

                k1 = f(t, X)
                k2 = f(t + dt/2, X + dt*rk1/2)
                k3 = f(t + dt/2, X + dt*rk2/2)
                k4 = f(t + dt, X + dt*rk3)

                X(i) = X(i) + dt*(k1 + 2 * k2 + 2 * k3 + k4)/6

            end do 

            N = INT(T/dt) 
        end subroutine RK4

end module schema