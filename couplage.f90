program couplage

    use mod_param
    use mod_algebre
    use mod_fonctions
    use mod_schemas


    implicit none
    real(PR), dimension(0:513) :: T, T_f
    real(PR) :: cfl, cfl0  
    integer :: i, Nx
    !call save_temp(schema=2, ci=1, cl=1, L=0.2_PR, imax=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)
    !call save_temp(schema=3, ci=1, cl=1, L=0.2_PR, imax=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)

    !call save_temp_2D(schema=2, ci=1, cl=2, Lx=0.2_PR, Ly=0.2_PR, nx=128, ny=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)
    !call save_temp_2D_2W(schema=2, ci=1, cl=1, Lx=0.2_PR, Ly=0.2_PR, nx=128, &
    !ny=128, tmax=200._PR, type_bois_1=1, type_bois_2=5, cfl=2._PR, freq=100)

    
    ! T= 0._PR
    ! T_f = 0._PR
    ! Nx = 512
    ! cfl0 = 0.02_PR
    ! call save_vect_temp(schema=2, T=T, ci=1, cl=1, L=0.2_PR, imax=Nx, tmax=20._PR, type_bois=1, cfl=cfl0)

    ! print *, "Euler Implicite"
    ! do i = 1, 2
    !     T_f = 0._PR
    !     cfl = cfl0/(2**i)
    !     call save_vect_temp(schema=2, T=T_f, ci=1, cl=1, L=0.2_PR, imax=Nx, tmax=20._PR, type_bois=1, cfl=cfl)
    !     call calcul_L2_err(T, T_f, 2, Nx, cfl)
    ! end do
    
    ! print *, '------------------'
    ! print *, "Crank Nicolson"
    ! call save_vect_temp(schema=3, T=T, ci=1, cl=1, L=0.2_PR, imax=Nx, tmax=20._PR, type_bois=1, cfl=cfl0)
    ! do i = 1, 2
    !     T_f = 0._PR
    !     cfl = cfl0/(2**i)
    !     call save_vect_temp(schema=3, T=T_f, ci=1, cl=1, L=0.2_PR, imax=Nx, tmax=20._PR, type_bois=1, cfl=cfl)
    !     call calcul_L2_err(T, T_f, 3, Nx, cfl)
    ! end do

    contains 

        subroutine calcul_L2_err(T, T_f, schema, Nx, cfl)
            real(PR), dimension(0:Nx) :: T, T_f
            real(PR) :: cfl, Err
            integer :: i, schema, Nx

            Err = 0._PR
            do i=0, Nx
                Err = Err + (T(i) - T_f(i))**2
            end do
            Err = sqrt(Err*0.2_PR/(Nx+1))

            print *, 'Delta t : ', cfl/20._PR
            print *, 'Erreur L2 : ', Err
        end subroutine calcul_L2_err

end program couplage