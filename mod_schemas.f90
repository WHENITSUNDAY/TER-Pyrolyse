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
            real(PR), dimension(0:imax+1)   :: T, Tnew , rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda
            real(PR), dimension(imax, imax) :: A, A1, A2
            real(PR), dimension(3)          :: k, knew
            real(PR)                        :: tn, dt, dx, eta, khi, mv_bois, xi, rhoCp0, lambda0
            integer                         :: i, ct
            character(20)                   :: chci, chcl, chl, chbois, chtn, chimax, chcfl, chschema
            character(len=200)              :: fichier

            
            select case (schema)
                case (1)
                    write(chschema,*) 'EE'
                case (2)
                    write(chschema,*) 'EI'
                case (3)
                    write(chschema,*) 'CK'
            end select

            write(chci,'(1I1)') ci
            write(chcl,'(1I1)') cl
            write(chl,'(1F4.0)') L
            write(chbois,*) type_bois 
            write(chimax,'(1I4)') imax
            write(chcfl,'(1F5.3)') cfl



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
            ct = 1
            dx = 1._PR/(imax+1)
            
            !Calcul CFL
            lambda0 = (0.166_PR + 0.369_PR * khi) 
            rhoCp0 = rho_b(0)*Cp(1) + rho_l(0)*Cp(4)   
            dt = cfl * rhoCp0 * dx**2 / (2*lambda0)

            print *, dt
            !Initialisation Température

            T(0) = Tg(0._PR, cl)
            Tnew(0) = Tg(0._PR, cl)

            do i = 1, imax+1
                xi = i*dx
                T(i) = Tinit(xi, ci)
                Tnew(i) = Tinit(xi, ci)
            end do

            print *, T
            do while (tn < tmax)

                if ( mod(ct, freq) == 0) then

                    write(chtn,'(1F6.1)') tn

                    fichier = 'data/temp_'//trim(adjustl(chschema)) // &
                        & '_ci'// trim(adjustl(chci)) &
                        & // '_cl' // trim(adjustl(chcl)) &
                        & // '_L'// trim(adjustl(chl)) &
                        & // '_bois'// trim(adjustl(chbois)) &
                        & // '_tn' // trim(adjustl(chtn)) &
                        & // '_imax' // trim(adjustl(chimax)) &
                        & // '_cfl' // trim(adjustl(chcfl)) //'.dat'

                    open(unit = 100, file=fichier)
 
                end if
                
                do i = 0, imax+1

                    rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)
                    eta = rho_c(i) / (rho_b(i) + rho_c(i))
                    lambda(i) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)

                    
                end do


                select case (schema)

                    case (1) !Euler Explicite
                        
                    case (2) !Euler Implicite
                        
                        call remplissage_A(A, imax, lambda, -1._PR, rhoCp, dx, dt)

                        T(0) = Tg(tn, CL)
                        T(1) = T(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2)

                        print *, tn, T(1)
                        call lu_tridiagonal(imax, A, T(1:imax), Tnew(1:imax))
                        
                        Tnew(imax+1) = Tnew(imax)
                    case (3) !Crank-Nicolson

                        T(0) = Tg(tn, CL)
                        print *, dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2)
                        T(1) = T(1) + T(0)*dt*(2.0_PR * lambda(1) * lambda(0) / (lambda(1) + lambda(0)))/(rhoCp(1)*dx**2)

                        print *, T(1)
                        call remplissage_A(A1, imax, lambda, -1._PR, rhoCp, dx, dt)
                        call remplissage_A(A2, imax, lambda, 1._PR, rhoCp, dx, dt)

                        call lu_tridiagonal(imax, A1, MATMUL(A2,T(1:imax)), Tnew(1:imax))
                end select 

                do i = 0, imax+1

                    call arrhenius(k, T(i))
                    call arrhenius(knew, Tnew(i))
                    call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k, knew, dt)

                    

                    if ( mod(ct, freq) == 0) then

                        write(100,*) dx*i, Tnew(i), rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i)

                    end if

                    
                end do

                T = Tnew
                tn = tn + dt
                ct = ct + 1

                if ( mod(ct, freq) == 0) then

                    close(100)

                end if
            end do 
        end subroutine save_temp

end module mod_schemas