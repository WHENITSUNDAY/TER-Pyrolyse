program couplage

    use mod_param
    use mod_algebre
    use mod_fonctions
    use mod_schemas


    implicit none

    !call save_temp(schema=2, ci=1, cl=1, L=0.2_PR, imax=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)
    !call save_temp(schema=3, ci=1, cl=1, L=0.2_PR, imax=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)

    call save_temp_2D(schema=2, ci=1, cl=2, Lx=0.2_PR, Ly=0.2_PR, nx=128, ny=128, tmax=200._PR, type_bois=1, cfl=2._PR, freq=100)
    !call save_temp_2D_2W(schema=2, ci=1, cl=1, Lx=0.2_PR, Ly=0.2_PR, nx=128, &
    !ny=128, tmax=200._PR, type_bois_1=1, type_bois_2=5, cfl=2._PR, freq=100)
end program couplage