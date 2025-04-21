program couplage

    use mod_constantes
    use mod_algebre
    use mod_schema

    implicit none

    real(PR)                :: tf, L
    integer                 :: type_bois, Nx, schema
    
    
    L = 0.1_PR
    Nx = 64
    tf = 500
    type_bois = 1
    schema = 1

    contains


        subroutine save_temp(Nx, tf, type_bois, freq, schema)
            integer, intent(in)             :: Nx, type_bois, freq, schema 
            real(PR), intent(in)            :: tf
            real(PR)                        :: tn, dt, dx, eta, khi, mv_bois
            real(PR), dimension(0:Nx+1)     :: T, Tnew , rho_b, rho_c, rho_g, rho_l, rho_v, rhoCp, lambda
            real(PR), dimension(Nx)         :: diag
            real(PR), dimension(Nx-1)       :: udiag
            real(PR), dimension(Nx, Nx)     :: A, Id
            real(PR), dimension(3)          :: k, knew
            integer                         :: i, ct
            character(len=20)               :: tn_str, sch_str

            

            !Initialisation Schéma

            tn = 0._PR
            ct = 0
            dx = 1._PR/(Nx+1)
            Id = 0._PR

            do i = 1, Nx
                Id(i,i) = 1._PR
            end do
            !Initialisation Température
            T(0) = 800._PR

            do i = 1, Nx+1
                T(i) = 300._PR !K
            end do

            !Initialisation Chimie
            khi = humidite(type_bois)
            mv_bois = densite_bois(type_bois)
            rho_b = mv_bois*(1._PR-khi)
            rho_c = 0._PR
            rho_g = 0._PR
            rho_l = mv_bois*khi
            rho_v = 0._PR
            
            do while (tn < tf)

                if ( mod(ct, freq) == 0) then

                    write(tn_str,'(F10.2)') tn
                    write(sch_str,*) schema
                    open(unit = 100, file='data/temp_t'//trim(adjustl(tn_str))//'_sch'//trim(adjustl(sch_str))//'.dat')
 
                end if
                
                do i = 0, Nx+1
                    rhoCp(i) = rho_b(i)*Cp(1) + rho_c(i)*Cp(2) + rho_g(i)*Cp(3) + rho_l(i)*Cp(4) + rho_v(i)*Cp(5)
                    eta = rho_b(i) / (rho_b(i) + rho_c(i))
                    lambda(i) = eta * 0.105_PR + (1.0_PR - eta) * (0.166_PR + 0.369_PR * khi)
                end do

                call remplissage_A(A, Nx, lambda, rhoCp, dx, dt)

                !Résolution du système linéaire en fonction du schéma
                select case (schema)

                    case (1) !Euler Explicite
                        
                    case (2) !Euler Implicite

                    case (3) !Crank-Nicolson

                end select 


                do i = 0, Nx+1

                    call arrhenius(k, T(i))
                    call arrhenius(knew, Tnew(i))
                    call CK2_step(rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i), k, knew, dt)

                    if ( mod(ct, freq) == 0) then

                        write(101,*) dx*i, Tnew(i), rho_b(i), rho_c(i), rho_g(i), rho_l(i), rho_v(i)

                    end if

                end do

                T = Tnew
                t = t + dt

                ct = ct + 1

                close(100)
            end do 
        end subroutine save_temp
end program