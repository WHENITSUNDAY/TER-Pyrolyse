module mod_algebre

    use mod_param

    implicit none
    
    contains

        subroutine remplissage_A(A, imax, lambda, sign, rhoCp, dx, dt)
            real(PR), dimension(0:imax+1) :: lambda, rhoCp
            real(PR), dimension(1:imax, 1:imax) :: A
            real(PR) :: lambda_ph, lambda_mh, alpha, dx, dt, sign
            integer :: i, imax

            A = 0.0_PR

            lambda_mh = 2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0))

            do i = 1, imax
                lambda_ph = 2.0_PR * lambda(i) * lambda(i+1) / (lambda(i) + lambda(i+1))
                !lambda_mh = 2.0_PR * lambda(i) * lambda(i-1) / (lambda(i) + lambda(i-1))
                alpha = sign*dt/(rhoCp(i) * dx**2)

                if (i == 1) then

                    A(i, i) = 1._PR - alpha*(lambda_mh + lambda_ph)
                    A(i, i+1) = alpha * lambda_ph

                else if (i == imax) then

                    A(i, i) = 1._PR - alpha * lambda_mh 
                    A(i, i-1) = alpha * lambda_mh

                else
                    A(i, i) = 1._PR - alpha * (lambda_ph + lambda_mh)
                    A(i, i+1) = alpha * lambda_ph
                    A(i, i-1) = alpha * lambda_mh
                end if

                lambda_mh = lambda_ph
                
            end do
        
        end subroutine remplissage_A

        subroutine lu_tridiagonal(n, A, b, x)
            implicit none
            integer, intent(in) :: n
            real(PR), intent(in) :: A(n, n)
            real(PR), intent(in) :: b(n)
            real(PR), intent(out) :: x(n)
        
            real(PR), dimension(n) :: d, y
            real(PR), dimension(n-1) :: l, u
            integer :: i

            do i = 1, n
                d(i) = A(i, i)
            end do

            do i = 1, n - 1
                u(i) = A(i, i+1)
                l(i) = A(i+1, i) / d(i)
                d(i+1) = d(i+1) - l(i) * u(i)
            end do

            y(1) = b(1)
            do i = 2, n
                y(i) = b(i) - l(i-1) * y(i-1)
            end do

            x(n) = y(n) / d(n)
            do i = n-1, 1, -1
                x(i) = (y(i) - u(i) * x(i+1)) / d(i)
            end do

        end subroutine lu_tridiagonal


        subroutine print_mat(A)
            
            real(PR), dimension(:,:)    :: A
            integer                     :: i

            print *, "Matrice :"
            do i = 1, size(A,1)
                
                print *, A(i,:)
            
            end do

        end subroutine print_mat


        subroutine lu_decomp(A, M)
            real(PR), dimension(:,:), intent(in)        :: A
            real(PR), dimension(:,:), intent(inout)     :: M
            integer                                     :: n, j, k

            n = size(A,1)

            M = A
            
            do k = 1, n - 1
                do j = k+1, n
                    
                    M(j,k) = M(j,k)/M(k,k)
                    M(j,k+1:n) = M(j,k+1:n) - M(j,k)*M(k,k+1:n)

                end do
            end do
        
        end subroutine lu_decomp


        subroutine lu_res(M, b, x)
            real(PR), dimension(:,:), intent(in)    :: M
            real(PR), dimension(:), intent(in)      :: b
            real(PR), dimension(:), intent(out)     :: x
            real(PR), dimension(size(x))            :: y
            integer                                 :: i, n
        
            n = size(M,1)

            y = 0._PR

            do i = 1, n
                y(i) = b(i) - sum(M(i,1:i-1) * y(1:i-1)) 
            end do

            do i = n, 1, -1
                x(i) = y(i) - sum(M(i,i+1:n) * x(i+1:n))
                x(i) = x(i) / M(i,i)
            end do
        
        end subroutine lu_res


        subroutine chol_decomp(A, L)
            real(PR), dimension(:,:), intent(in)        :: A
            real(PR), dimension(:,:), intent(inout)     :: L
            integer                                     :: n, i, j

            
            n = size(A,1)
            L = 0._PR

            do i = 1, n
                L(i,i) = SQRT(A(i,i) - sum(L(i,1:i-1) ** 2))

                do j = i+1, n
                    L(j,i) = (A(j,i) - sum(L(i,1:i-1) * L(j,1:i-1))) / L(i,i)
                end do
            end do

        end subroutine chol_decomp

        subroutine chol_res(L, b, x)
            real(PR), dimension(:,:), intent(in)        :: L
            real(PR), dimension(:), intent(in)          :: b
            real(PR), dimension(:), intent(out)         :: x
            real(PR), dimension(size(x))                :: y
            real(PR), dimension(size(L,1),size(L,2))    :: Lt
            integer                                     :: i, n
        
            n = size(L, 1)
            Lt = TRANSPOSE(L)

            y = 0._PR

            do i = 1, n
                y(i) = (b(i) - sum(L(i,1:i-1) * y(1:i-1))) / L(i,i)
            end do

            do i = n, 1, -1
                x(i) = (y(i) - sum(Lt(i,i+1:n) * x(i+1:n))) / Lt(i,i)
            end do
        
        end subroutine chol_res



end module mod_algebre