program temp

    use cstes

    implicit none

    integer  :: nt, nx
    real(PR) :: temp_im, temp_i, temp_ip, temp_ext ! temp_i-1, temp_i, temp_i+1, temp(x=0)
    real(PR) :: t, dt, t_max, temp_initiale, x, x_max, dx
    real(PR) :: rho, lambda, c_p, d
    integer  :: k
    real(PR), dimension(:), allocatable :: x_b ! tableau sur la longueur du bois
    

    rho    = 1
    lambda = 1
    c_p    =  1
    temp_initiale = 600  ! temp en K

    d = lambda/(c_p*rho)
    temp_ext = 280

    x_max = 1
    dx    = 0.1
    nx    = INT(x_max/dx)+1
    dx    = x_max/nx 
    allocate(x_b(nx))
    do k = 1, size(x_b)
        x_b(k) = temp_initiale
        
    end do
    x_b(1)  = temp_ext
    x_b(nx) = temp_ext
    
    t_max = 0.1
    dt    = 0.0002
    !nt    = INT(t_max/dt)+1
    !dt    = t_max/nt 
    t     = 0 

    !temp_im = 1
    !temp_i  = 1
    !temp_ip = 1
    
    open(unit=1, file='temperature.dat')

    write(1, *) t, x_b(:) 

    do while (t <= t_max)

        x = 0
        k = 2

        do k=2,nx-1
            
            temp_im = x_b(k-1)
            temp_ip = x_b(k+1)
            temp_i = (x_b(k) + dt*d/dx**2*(temp_ip+temp_im))/(1+2*d*dt/dx**2)
            x_b(k) = temp_i

        end do

        !x_b(nx) = x_b(nx-1)

        t = t + dt

        write(1, *) t, x_b(:)

    end do


    
end program temp
