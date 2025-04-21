module mod_algebre
    
    use mod_constantes
    
    implicit none

    contains 

        subroutine remplissage_A(A, Nx, lambda, rhoCp, dx, dt)
            real(PR), dimension(0:Nx+1) :: lambda, rhoCp
            real(PR), dimension(1:Nx, 1:Nx) :: A
            real(PR) :: lambda_ph, lambda_mh, alpha, dx, dt
            integer :: i, Nx

            A = 0.0_PR

            lambda_mh = 2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0))

            do i = 1, Nx
                lambda_ph = 2.0_PR * lambda(i) * lambda(i+1) / (lambda(i) + lambda(i+1))

                alpha = dt/(rhoCp(i) * dx**2)

                if (i == 1) then

                    A(i, i) = -alpha * (lambda_mh + lambda_ph)
                    A(i, i+1) = alpha * lambda_ph

                else if (i == Nx) then

                    A(i, i) = alpha * lambda_mh
                    A(i, i-1) = alpha * lambda_mh

                else
                    A(i, i) = -alpha * (lambda_ph + lambda_mh)
                    A(i, i+1) = alpha * lambda_ph
                    A(i, i-1) = alpha * lambda_mh
                end if

                lambda_mh = lambda_ph
                
            end do
        
        end subroutine remplissage_A

        subroutine print_mat(A)
            
            real(PR), dimension(:,:)    :: A
            integer                     :: i

            do i = 1, size(A,1)
                
                print *, A(i,:)
            
            end do

        end subroutine print_mat


        subroutine lu_decomp(A, M)
            real(PR), dimension(:,:), intent(in)        :: A
            real(PR), dimension(:,:), intent(inout)     :: M
            integer                                     :: n, j, k

            n = size(A,1)

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
            integer                                 :: i, j, n

            n = size(A,1)

            y = 0.0_PR

            do i = 1, n
                y(i) = b(i)
                do j = 1, i-1
                    y(i) = y(i) + M(i,j)*y(j)
                end do
            end do

            do i = 1, n
                x(i) = y(i)
                do j = i+1, n
                    x(i) = x(i) - (M(i,j)*y(j))/M(i,i)
                end do
            end do

        end subroutine lu_res
        
        

end module mod_algebre