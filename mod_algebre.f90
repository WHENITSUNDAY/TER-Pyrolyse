module mod_algebre
    
    use constantes
    
    implicit none

    contains 

        function A_n(n, dx)

            integer                         :: n, i
            real(PR)                        :: dx
            real(PR), dimension(n,n)        :: A_n

            A_n = 0._PR

            do i = 1, n

                A_n(i,i) = -2._PR

                if (i > 1) then
                    
                    A_n(i-1, i) = 1._PR
                end if

                if (i < n) then
                    
                    A_n(i+1, i) = 1._PR
                end if
            end do

            A_n = A_n * (-1._PR/(dx**2))
                    
        end function A_n
        

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
            real(PR)                                    :: s_1, s_2

            n = size(A,1)

            do k = 1, n - 1
                do j = k+1, n
                    
                    M(j,k) = M(j,k)/M(k,k)
                    M(j,k+1:n) = M(j,k+1:n) - M(j,k)*M(k,k+1:n)

                end do
            end do
        
        end subroutine lu_decomp


        subroutine lu_res(M, b, x, check)
            real(PR), dimension(:,:), intent(in)    :: M
            real(PR), dimension(:), intent(in)      :: b
            real(PR), dimension(:), intent(out)     :: x
            real(PR), dimension(size(x))            :: y
            real(PR)                                :: s
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
        
        subroutine chol_decomp(A, L)
            real(PR), dimension(:,:), intent(in)        :: A
            real(PR), dimension(:,:), intent(inout)     :: L
            integer                                     :: n, i, j, k
            real(PR)                                    :: s_1, s_2

            
            n = size(A,1)
            L = 0._PR

            do i = 1, n
                s_1 = 0._PR

                do k = 1, i-1
                    s_1 = s_1 + L(i,k)**2
                end do

                L(i,i) = SQRT(A(i,i) - s_1)

                do j = i+1, n
                    s_2 = 0._PR

                    do k = 1, i-1
                        s_2 = s_2 + L(i,k)*L(j,k)
                    end do

                    L(j,i) = (A(j,i) - s_2)/L(i,i)

                end do
            end do
        
        end subroutine chol_decomp


        subroutine chol_res(L, b, x)
            real(PR), dimension(:,:), intent(in)        :: L
            real(PR), dimension(:), intent(in)          :: b
            real(PR), dimension(:), intent(out)         :: x
            real(PR), dimension(size(x))                :: y
            real(PR), dimension(size(L,1),size(L,2))    :: Lt
            real(PR)                                    :: s_3, s_4
            integer                                     :: i, j, n
        
            n = size(L, 1)
            Lt = TRANSPOSE(L)

            y = 0.0_PR

            do i = 1, n
                s_3 = 0.0_PR

                do j = 1, i-1
                    s_3 = s_3 + L(i,j) * y(j)
                end do

                y(i) = (b(i) - s_3) / L(i,i)
            end do

            do i = n, 1, -1
                s_4 = 0.0_PR

                do j = i+1, n
                    s_4 = s_4 + Lt(i,j) * x(j)
                end do

                x(i) = (y(i) - s_4) / Lt(i,i)

            end do
        
        end subroutine chol_res

end module mod_algebre