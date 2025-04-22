program couplage

    use mod_param
    use mod_algebre
    use mod_fonctions
    use mod_schemas


    implicit none
 
    call save_temp(schema=2, ci=1, cl=1, L=0.10_PR, imax=64, tmax=500._PR, type_bois=1, cfl=1._PR, freq=20)


end program couplage