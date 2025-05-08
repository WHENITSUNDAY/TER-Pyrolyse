module mod_schemas

    use mod_param
    use mod_fonctions
    use mod_algebre

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

        subroutine save_temp(schema, ci, cl, L, imax, tmax, type_bois, cfl, freq)

            integer, intent(in)             :: imax, type_bois, freq, schema, ci, cl 
            real(PR), intent(in)            :: cfl, tmax, L
            real(PR), dimension(0:imax+1)   :: T, Tnew , rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda, Qr
            real(PR), dimension(imax, imax) :: A, A1, A2
            real(PR), dimension(imax)       :: b
            real(PR), dimension(3, 0:imax+1):: k, knew
            real(PR)                        :: tn, dt, dx, eta, khi, mv_bois, xi, rhoCp0, lambda0, t1, t2
            integer                         :: i, ct
            character(20)                   :: chl, chbois, chtn, chimax, chschema
            character(len=200)              :: fichier

            
            call cpu_time(t1)

            select case (schema)
                case (1)
                    write(chschema,*) 'EE'
                case (2)
                    write(chschema,*) 'EI'
                case (3)
                    write(chschema,*) 'CK'
            end select

            !write(chci,'(1I1)') ci
            !write(chcl,'(1I1)') cl
            write(chl,'(1F4.1)') L
            write(chbois,*) type_bois 
            write(chimax,'(1I4)') imax
            !write(chcfl,'(1F6.2)') cfl



            !Initialisation Chimie
            khi = humidite(type_bois)
            mv_bois = densite_bois(type_bois)
            rho_b = mv_bois*(1._PR-khi)
            rho_c = 0._PR
            rho_g = 0._PR
            rho_l = mv_bois*khi
            rho_v = 0._PR

            
            !Initialisation Schéma
            tn = 0._PR
            ct = 0
            dx = L/(imax+1)
            
            !Calcul dt
            select case (schema)

                case (1)

                    lambda0 = (0.166_PR + 0.369_PR * khi) 
                    rhoCp0 = rho_b(0)*Cp(1) + rho_l(0)*Cp(4)   
                    dt = cfl * rhoCp0 * dx**2 / (2*lambda0)

                case (2)
                    dt = cfl/tmax !Schema implicite : pas de condition CFL
                    
                case(3)
                    dt = cfl/tmax
            
            end select

            print *, "Pas de temps :", dt

            !Initialisation Température
            T(0) = Tg(0._PR, cl)
            Tnew(0) = Tg(0._PR, cl)

            do i = 1, imax+1
                xi = i*dx
                T(i) = Tinit(xi, ci)
                Tnew(i) = Tinit(xi, ci)
                
            end do

            do i = 0, imax+1
                call arrhenius(k(:,i), T(i))
            end do

            do while (tn < tmax)

                if ( mod(ct, freq) == 0) then

                    write(chtn,'(1F6.1)') tn

                    fichier = 'data/1D/temp_'//trim(adjustl(chschema)) // &
                        & '_L'// trim(adjustl(chl)) &
                        & // '_bois'// trim(adjustl(chbois)) &
                        & // '_tn' // trim(adjustl(chtn)) &
                        & // '_imax' // trim(adjustl(chimax)) //'.dat' !&
                        !& // '_cfl' // trim(adjustl(chcfl)) //'.dat'

                    open(unit = 100, file=fichier)
 
                end if
                
                do i = 0, imax+1

                    rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)
                    eta = rho_c(i) / (rho_b(i) + rho_c(i))
                    lambda(i) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)

                    Qr(i) = k(1,i) * rho_b(i) * (dH(1) + (Cp(2) - Cp(1)) * (T(i) - Tinit(i*dx, ci))) + &
                        & k(2,i) * rho_b(i) * (dH(2) + (Cp(3) - Cp(1)) * (T(i) - Tinit(i*dx, ci))) + &
                        & k(3,i) * rho_l(i) * (dH(3) + (Cp(5) - Cp(4)) * (T(i) - Tinit(i*dx, ci)))
                end do


                select case (schema)

                    case (1) !Euler Explicite
                        
                    case (2) !Euler Implicite
                        
                        call remplissage_A(A, imax, lambda, -1._PR, rhoCp, dx, dt)
                        T(0) = Tg(tn, cl)
                        Tnew(0) = Tg(tn, cl)

                        b = T(1:imax) + dt*Qr(1:imax)/rhoCp(1:imax)

                        !La condition de Dirichlet à gauche impose :
                        b(1) = T(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2) &
                        + dt*Qr(1)/rhoCp(1)

                        call lu_tridiagonal(imax, A, b, Tnew(1:imax))
                        
                        !La condition de Neumann à droite impose :
                        Tnew(imax+1) = Tnew(imax)

                    case (3) !Crank-Nicolson

                        T(0) = Tg(tn, cl)

                        call remplissage_A(A1, imax, lambda, -0.5_PR, rhoCp, dx, dt)
                        call remplissage_A(A2, imax, lambda, 0.5_PR, rhoCp, dx, dt)

                        b = MATMUL(A2, T(1:imax)) + dt*Qr(1:imax)/rhoCp(1:imax)

                        !La condition de Dirichlet à gauche impose :
                        b(1) = b(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2)

                        
                        call lu_tridiagonal(imax, A1, b, Tnew(1:imax))

                        Tnew(imax+1) = Tnew(imax)
                end select 

                do i = 0, imax+1

                    call arrhenius(knew(:,i), Tnew(i))
                    call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k(:,i), knew(:,i), dt)
                    if ( mod(ct, freq) == 0) write(100,*) dx*i, Tnew(i), rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i)
                    
                end do

                T = Tnew
                tn = tn + dt
                k = knew
                ct = ct + 1
                
                if ( mod(ct, freq) == 0) close(100)

            end do 

            call cpu_time(t2)

            print *, "1D, Temps d'exécution :", t2-t1
        end subroutine save_temp


        subroutine save_temp_2D(schema, ci, cl, Lx, Ly, nx, ny, tmax, type_bois, cfl, freq)

            integer, intent(in)             :: nx, ny, type_bois, freq, schema, ci, cl 
            real(PR), intent(in)            :: cfl, tmax, Lx, Ly
            real(PR), dimension(0:nx+1, 0:ny+1) :: T, Tnew, rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda, Qr
            real(PR), dimension(nx*ny)      :: b
            real(PR), dimension(3, 0:nx+1, 0:ny+1) :: k, knew
            real(PR)                        :: tn, dt, dx, dy, eta, khi, mv_bois, rhoCp0, lambda0, t1, t2, err, tol
            real(PR)                        :: alpha_x, alpha_y, gamma, lambda_ph_x, lambda_mh_x, lambda_ph_y, lambda_mh_y
            integer                         :: i, j, ct, ite, max_ite
            character(20)                   :: chLx, chLy, chbois, chtn, chnx, chny, chschema
            character(len=200)              :: fichier

            call cpu_time(t1)

            select case (schema)
                case (1)
                    write(chschema,*) 'EE'
                case (2)
                    write(chschema,*) 'EI'
                case (3)
                    write(chschema,*) 'CK'
            end select

            write(chLx,'(1F4.1)') Lx
            write(chLy,'(1F4.1)') Ly
            write(chbois,*) type_bois 
            write(chnx,'(1I4)') nx
            write(chny,'(1I4)') ny

            ! Initialisation Chimie
            khi = humidite(type_bois)
            mv_bois = densite_bois(type_bois)
            rho_b = mv_bois * (1._PR - khi)
            rho_c = 0._PR
            rho_g = 0._PR
            rho_l = mv_bois * khi
            rho_v = 0._PR

            ! Initialisation Schéma
            tol = 1.e-6_PR
            tn = 0._PR
            ct = 0
            dx = Lx / (nx + 1)
            dy = Ly / (ny + 1)

            ! Calcul dt
            select case (schema)
                case (1)
                    lambda0 = (0.166_PR + 0.369_PR * khi) 
                    rhoCp0 = rho_b(0,0)*Cp(1) + rho_l(0,0)*Cp(4)   
                    dt = cfl * rhoCp0 * min(dx**2, dy**2) / (4 * lambda0)
                case (2)
                    dt = cfl / tmax
                case (3)
                    dt = cfl / tmax
            end select

            print *, "Pas de temps :", dt

            do j = 0, ny+1
                do i = 0, nx+1
                    T(i,j) = Tinit(0._PR, ci)
                    Tnew(i,j) = Tinit(0._PR, ci)
                end do
            end do


            T(0,nx/4:3*nx/4) = Tg(0._PR, cl)
            Tnew(0,nx/4:3*nx/4) = Tg(0._PR, cl)
            
            do j = 0, ny+1
                do i = 0, nx+1
                    call arrhenius(k(:, i, j), T(i, j))
                end do
            end do

            do while (tn < tmax)
                if (mod(ct, freq) == 0) then
                    write(chtn,'(1F6.1)') tn

                    fichier = 'data/2D/temp_'//trim(adjustl(chschema)) // &
                    & '_Lx'// trim(adjustl(chLx)) // &
                    & '_Ly'// trim(adjustl(chLy)) //&
                    & '_bois'// trim(adjustl(chbois)) //&
                    & '_tn' // trim(adjustl(chtn)) //&
                    & '_nx' // trim(adjustl(chnx)) //'.dat'

                    open(unit = 100, file=fichier)
                end if

                do j = 0, ny+1
                    do i = 0, nx+1
                        rhoCp(i, j) = rho_b(i, j)*Cp(1) + rho_c(i, j)*Cp(2) + rho_g(i, j)*Cp(3) + &
                            rho_l(i, j)*Cp(4) + rho_v(i, j)*Cp(5)
                        eta = rho_c(i, j) / (rho_b(i, j) + rho_c(i, j))
                        lambda(i, j) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)

                        Qr(i, j) = k(1, i, j) * rho_b(i, j) * (dH(1) + (Cp(2) - Cp(1)) * (T(i, j) - Tinit(i*dx, ci))) + &
                        k(2, i, j) * rho_b(i, j) * (dH(2) + (Cp(3) - Cp(1)) * (T(i, j) - Tinit(i*dx, ci))) + &
                        k(3, i, j) * rho_l(i, j) * (dH(3) + (Cp(5) - Cp(4)) * (T(i, j) - Tinit(i*dx, ci)))
                    end do
                end do

                select case (schema)
                    case (1) ! Euler Explicite
                    case (2) ! Euler Implicite
                        !Résolution via Gauss-Seidel (sans matrice pentadiagonale)

                        err = 1._PR
                        Tnew = T
                        !Tnew(0,:) = Tg(tn, cl)
                        Tnew(0,nx/4:3*nx/4) = Tg(tn, cl)
                        
                        ite = 0
                        max_ite = 100
                        do while (err > tol .and. ite < max_ite)
                            err = 0._PR
                            do j = 1, ny
                                do i = 1, nx
                                lambda_ph_x = 2.0_PR * lambda(i, j) * lambda(i+1, j) / (lambda(i, j) + lambda(i+1, j))
                                lambda_mh_x = 2.0_PR * lambda(i, j) * lambda(i-1, j) / (lambda(i, j) + lambda(i-1, j))
                                lambda_ph_y = 2.0_PR * lambda(i, j) * lambda(i, j+1) / (lambda(i, j) + lambda(i, j+1))
                                lambda_mh_y = 2.0_PR * lambda(i, j) * lambda(i, j-1) / (lambda(i, j) + lambda(i, j-1))

                                alpha_x = dt / (rhoCp(i, j) * dx**2)
                                alpha_y = dt / (rhoCp(i, j) * dy**2)
                                gamma = 1.0 + alpha_x * (lambda_ph_x + lambda_mh_x) + &
                                        alpha_y * (lambda_ph_y + lambda_mh_y)

                                Tnew(i, j) = (T(i, j) + &
                                        alpha_x * (lambda_ph_x * Tnew(i+1, j) + lambda_mh_x * Tnew(i-1, j)) + &
                                        alpha_y * (lambda_ph_y * Tnew(i, j+1) + lambda_mh_y * Tnew(i, j-1)) &
                                        + dt * Qr(i, j) / rhoCp(i, j))/ gamma

                                err = max(err, abs(T(i,j) - Tnew(i,j)))
                                end do
                            end do
                            ite = ite + 1
                        end do  

                        do i = 1, nx
                            Tnew(i, ny+1) = Tnew(i, ny)
                            Tnew(i, 0) = Tnew(i, 1)
                        end do

                        do j = 0, ny+1
                            Tnew(nx+1, j) = Tnew(nx, j)
                        end do

                        do j = 0, ny/4-1
                            Tnew(0, j) = Tnew(1, j)
                        end do

                        do j = 3*ny/4+1, ny+1 
                            Tnew(0, j) = Tnew(1, j)
                        end do
                    case (3) ! Crank-Nicolson
                end select

                do j = 0, ny+1
                    do i = 0, nx+1

                        call arrhenius(knew(:, i, j), Tnew(i,j))
                        call CK2_step(rho_b(i, j), rho_c(i, j), rho_g(i, j), rho_l(i, j),&
                        & rho_v(i, j), k(:, i, j), knew(:, i, j), dt)
                        if (mod(ct, freq) == 0) write(100,*) dx*i, dy*j, Tnew(i, j), rho_b(i, j), rho_c(i, j) 
                    end do
                end do

                T = Tnew
                tn = tn + dt
                k = knew
                ct = ct + 1

                if (mod(ct, freq) == 0) close(100)

            end do

            call cpu_time(t2)

            print *, "2D, Temps d'exécution :", t2-t1
        end subroutine save_temp_2D

        subroutine save_temp_2D_2W(schema, ci, cl, Lx, Ly, nx, ny, tmax, type_bois_1, type_bois_2, cfl, freq)

            integer, intent(in)             :: nx, ny, type_bois_1, type_bois_2, freq, schema, ci, cl 
            real(PR), intent(in)            :: cfl, tmax, Lx, Ly
            real(PR), dimension(0:nx+1, 0:ny+1) :: T, Tnew, rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda, Qr
            real(PR), dimension(nx*ny)      :: b
            real(PR), dimension(3, 0:nx+1, 0:ny+1) :: k, knew
            real(PR)                        :: tn, dt, dx, dy, eta, khi, mv_bois, rhoCp0, lambda0, t1, t2, err, tol
            real(PR)                        :: alpha_x, alpha_y, gamma, lambda_ph_x, lambda_mh_x, lambda_ph_y, lambda_mh_y
            integer                         :: i, j, ct, ite, max_ite
            character(20)                   :: chLx, chLy, chbois1, chbois2, chtn, chnx, chny, chschema
            character(len=200)              :: fichier

            call cpu_time(t1)

            select case (schema)
            case (1)
                write(chschema,*) 'EE'
            case (2)
                write(chschema,*) 'EI'
            case (3)
                write(chschema,*) 'CK'
            end select

            write(chLx,'(1F4.1)') Lx
            write(chLy,'(1F4.1)') Ly
            write(chbois1,*) type_bois_1
            write(chbois2,*) type_bois_2
            write(chnx,'(1I4)') nx
            write(chny,'(1I4)') ny

            ! Initialisation Chimie
            do j = 0, ny+1
                if (j <= ny/2) then
                    khi = humidite(type_bois_1)
                    mv_bois = densite_bois(type_bois_1)
                    do i = 0, nx+1
                        rho_b(i, j) = mv_bois * (1._PR - khi)
                        rho_c(i, j) = 0._PR
                        rho_g(i, j) = 0._PR
                        rho_l(i, j) = mv_bois * khi
                        rho_v(i, j) = 0._PR
                    end do
                else
                    khi = humidite(type_bois_2)
                    mv_bois = densite_bois(type_bois_2)
                    do i = 0, nx+1
                        rho_b(i, j) = mv_bois * (1._PR - khi)
                        rho_c(i, j) = 0._PR
                        rho_g(i, j) = 0._PR
                        rho_l(i, j) = mv_bois * khi
                        rho_v(i, j) = 0._PR
                    end do
                end if
                
            end do

            ! Initialisation Schéma
            tol = 1.e-6_PR
            tn = 0._PR
            ct = 0
            dx = Lx / (nx + 1)
            dy = Ly / (ny + 1)

            ! Calcul dt
            select case (schema)
                case (1)
                    lambda0 = (0.166_PR + 0.369_PR * khi) 
                    rhoCp0 = rho_b(0,0)*Cp(1) + rho_l(0,0)*Cp(4)   
                    dt = cfl * rhoCp0 * min(dx**2, dy**2) / (4 * lambda0)
                case (2)
                    dt = cfl / tmax
                case (3)
                    dt = cfl / tmax
            end select
            print *, "Pas de temps :", dt

            do j = 0, ny+1
            do i = 0, nx+1
                T(i,j) = Tinit(0._PR, ci)
                Tnew(i,j) = Tinit(0._PR, ci)
            end do
            end do

            T(0,:) = Tg(0._PR, cl)
            Tnew(0,:) = Tg(0._PR, cl)
            
            do j = 0, ny+1
            do i = 0, nx+1
                call arrhenius(k(:, i, j), T(i, j))
            end do
            end do

            do while (tn < tmax)
            if (mod(ct, freq) == 0) then
                write(chtn,'(1F6.1)') tn

                fichier = 'data/2D_2W/temp_'//trim(adjustl(chschema)) // &
                & '_Lx'// trim(adjustl(chLx)) // &
                & '_Ly'// trim(adjustl(chLy)) //&
                & '_bois'// trim(adjustl(chbois1)) //&
                & '_bois'// trim(adjustl(chbois2)) //&
                & '_tn' // trim(adjustl(chtn)) //&
                & '_nx' // trim(adjustl(chnx)) //'.dat'

                open(unit = 100, file=fichier)
            end if

            do j = 0, ny+1
                if (j <= ny/2) then
                khi = humidite(type_bois_1)
                else
                khi = humidite(type_bois_2)
                end if
                do i = 0, nx+1
                rhoCp(i, j) = rho_b(i, j)*Cp(1) + rho_c(i, j)*Cp(2) + rho_g(i, j)*Cp(3) + &
                    rho_l(i, j)*Cp(4) + rho_v(i, j)*Cp(5)
                eta = rho_c(i, j) / (rho_b(i, j) + rho_c(i, j))
                lambda(i, j) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)

                Qr(i, j) = k(1, i, j) * rho_b(i, j) * (dH(1) + (Cp(2) - Cp(1)) * (T(i, j) - Tinit(i*dx, ci))) + &
                k(2, i, j) * rho_b(i, j) * (dH(2) + (Cp(3) - Cp(1)) * (T(i, j) - Tinit(i*dx, ci))) + &
                k(3, i, j) * rho_l(i, j) * (dH(3) + (Cp(5) - Cp(4)) * (T(i, j) - Tinit(i*dx, ci)))
                end do
            end do

            select case (schema)
                case (1) ! Euler Explicite

                case (2) ! Euler Implicite
                !Résolution via Gauss-Seidel (sans matrice pentadiagonale)

                err = 1._PR
                Tnew = T
                Tnew(0,:) = Tg(tn, cl)
                
                ite = 0
                max_ite = 100
                do while (err > tol .and. ite < max_ite)
                    err = 0._PR
                    do j = 1, ny
                    do i = 1, nx
                    lambda_ph_x = 2.0_PR * lambda(i, j) * lambda(i+1, j) / (lambda(i, j) + lambda(i+1, j))
                    lambda_mh_x = 2.0_PR * lambda(i, j) * lambda(i-1, j) / (lambda(i, j) + lambda(i-1, j))
                    lambda_ph_y = 2.0_PR * lambda(i, j) * lambda(i, j+1) / (lambda(i, j) + lambda(i, j+1))
                    lambda_mh_y = 2.0_PR * lambda(i, j) * lambda(i, j-1) / (lambda(i, j) + lambda(i, j-1))

                    alpha_x = dt / (rhoCp(i, j) * dx**2)
                    alpha_y = dt / (rhoCp(i, j) * dy**2)
                    gamma = 1.0 + alpha_x * (lambda_ph_x + lambda_mh_x) + &
                        alpha_y * (lambda_ph_y + lambda_mh_y)

                    Tnew(i, j) = (T(i, j) + &
                        alpha_x * (lambda_ph_x * Tnew(i+1, j) + lambda_mh_x * Tnew(i-1, j)) + &
                        alpha_y * (lambda_ph_y * Tnew(i, j+1) + lambda_mh_y * Tnew(i, j-1))  &
                        + dt * Qr(i, j) / rhoCp(i, j)) / gamma

                    err = max(err, abs(T(i,j) - Tnew(i,j)))
                    end do
                    end do
                    ite = ite + 1
                end do  

                do i = 1, nx
                    Tnew(i, ny+1) = Tnew(i, ny)
                    Tnew(i, 0) = Tnew(i, 1)
                end do

                do j = 0, ny+1
                    Tnew(nx+1, j) = Tnew(nx, j)
                end do

                case (3) ! Crank-Nicolson
                !Même principe avec Gauss Seidel...
            end select

            do j = 0, ny+1
                do i = 0, nx+1

                call arrhenius(knew(:, i, j), Tnew(i,j))
                call CK2_step(rho_b(i, j), rho_c(i, j), rho_g(i, j), rho_l(i, j),&
                & rho_v(i, j), k(:, i, j), knew(:, i, j), dt)
                if (mod(ct, freq) == 0) write(100,*) dx*i, dy*j, Tnew(i, j), rho_b(i, j), rho_c(i, j) 
                end do
            end do

            T = Tnew
            tn = tn + dt
            k = knew
            ct = ct + 1

            if (mod(ct, freq) == 0) close(100)

            end do

            call cpu_time(t2)

            print *, "2D_2W, Temps d'exécution :", t2-t1
        end subroutine save_temp_2D_2W


        subroutine save_vect_temp(schema, T, ci, cl, L, imax, tmax, type_bois, cfl)

            integer, intent(in)             :: imax, type_bois, schema, ci, cl 
            real(PR), intent(in)            :: cfl, tmax, L
            real(PR), dimension(0:imax+1)   :: T, Tnew , rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda, Qr
            real(PR), dimension(imax, imax) :: A, A1, A2
            real(PR), dimension(imax)       :: b
            real(PR), dimension(3, 0:imax+1):: k, knew
            real(PR)                        :: tn, dt, dx, eta, khi, mv_bois, xi
            integer                         :: i, ct

            !Initialisation Chimie
            khi = humidite(type_bois)
            mv_bois = densite_bois(type_bois)
            rho_b = mv_bois*(1._PR-khi)
            rho_c = 0._PR
            rho_g = 0._PR
            rho_l = mv_bois*khi
            rho_v = 0._PR

            
            !Initialisation Schéma
            tn = 0._PR
            ct = 0
            dx = L/(imax+1)
            
            !Calcul dt
            select case (schema)

                case (2)
                    dt = cfl/tmax !Schema implicite : pas de condition CFL
                    
                case(3)
                    dt = cfl/tmax
            
            end select

            T(0) = Tg(0._PR, cl)
            Tnew(0) = Tg(0._PR, cl)

            do i = 1, imax+1
                xi = i*dx
                T(i) = Tinit(xi, ci)
                Tnew(i) = Tinit(xi, ci)
                
            end do

            
            do i = 0, imax+1
                call arrhenius(k(:,i), T(i))
            end do

            do while (tn < tmax)
                do i = 0, imax+1

                    rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)
                    eta = rho_c(i) / (rho_b(i) + rho_c(i))
                    lambda(i) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)

                    Qr(i) = k(1,i) * rho_b(i) * (dH(1) + (Cp(2) - Cp(1)) * (T(i) - Tinit(i*dx, ci))) + &
                        & k(2,i) * rho_b(i) * (dH(2) + (Cp(3) - Cp(1)) * (T(i) - Tinit(i*dx, ci))) + &
                        & k(3,i) * rho_l(i) * (dH(3) + (Cp(5) - Cp(4)) * (T(i) - Tinit(i*dx, ci)))
                end do


                select case (schema)

                    case (2) !Euler Implicite
                        
                        call remplissage_A(A, imax, lambda, -1._PR, rhoCp, dx, dt)
                        T(0) = Tg(tn, cl)
                        Tnew(0) = Tg(tn, cl)

                        b = T(1:imax) + dt*Qr(1:imax)/rhoCp(1:imax)

                        !La condition de Dirichlet à gauche impose :
                        b(1) = T(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2) &
                        + dt*Qr(1)/rhoCp(1)

                        call lu_tridiagonal(imax, A, b, Tnew(1:imax))
                        
                        !La condition de Neumann à droite impose :
                        Tnew(imax+1) = Tnew(imax)

                    case (3) !Crank-Nicolson

                        T(0) = Tg(tn, cl)

                        call remplissage_A(A1, imax, lambda, -0.5_PR, rhoCp, dx, dt)
                        call remplissage_A(A2, imax, lambda, 0.5_PR, rhoCp, dx, dt)

                        b = MATMUL(A2, T(1:imax)) + dt*Qr(1:imax)/rhoCp(1:imax)

                        !La condition de Dirichlet à gauche impose :
                        b(1) = b(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2)

                        
                        call lu_tridiagonal(imax, A1, b, Tnew(1:imax))

                        Tnew(imax+1) = Tnew(imax)
                end select 

                do i = 0, imax+1

                    call arrhenius(knew(:,i), Tnew(i))
                    call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k(:,i), knew(:,i), dt)
                    
                end do


                T = Tnew
                tn = tn + dt
                k = knew
                ct = ct + 1

            end do 

        end subroutine save_vect_temp
end module mod_schemas