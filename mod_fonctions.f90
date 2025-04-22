module mod_fonctions

    use mod_param

    implicit none

    contains

        function Tinit(x, ci)
            
            real(PR) :: Tinit, x
            integer :: ci

            select case (ci)
            
                case (1)

                    Tinit = 300._PR

                case (2)

                    Tinit = exp(-(x-5._PR)**2) 
            end select

        end function Tinit


        function Tg(t, cl)

            real(PR) :: Tg, t
            integer :: cl

            select case (cl)
                
                case (1)

                    Tg = 800._PR

                case (2)

                    Tg = 0._PR
                
                case (3)

                    Tg = 0._PR
            end select

        end function Tg

        subroutine arrhenius(k, Temp)

            real(PR), dimension(:), intent(out) :: k
            real(PR), intent(in)                :: Temp
            integer                             :: i

            k = (/(A(i)*exp(-E(i)/(R*Temp)), i = 1, 3)/)

        end subroutine arrhenius

end module mod_fonctions