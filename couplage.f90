program couplage

    use mod_param
    use mod_algebre
    use mod_fonctions
    use mod_schemas


    implicit none

    call save_temp(schema=2, ci=1, cl=1, L=0.1_PR, imax=64, tmax=500._PR, type_bois=1, cfl=10._PR, freq=5)


end program couplage