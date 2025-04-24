module mod_param

    implicit none

    integer, parameter                  :: PR = 8
    real(PR), parameter                 :: R = 8.3144261_PR !Constantes des gaz parfaits
    real(PR), dimension(6), parameter   :: densite_bois = [888, 740, 582, 469, 360, 469] 
    real(PR), dimension(6), parameter   :: humidite = [0.149_PR , 0.153_PR, 0.144_PR, 0.141_PR, 0.146_PR, 0.142_PR]
    real(PR), dimension(3), parameter   :: A_arr = [7.38e5_PR, 1.44e4_PR, 5.13e10_PR] !Paramètres Arrhenius
    real(PR), dimension(3), parameter   :: E_arr = [106.5e3_PR, 88.6e3_PR, 88e3_PR] !Paramètres Arrhenius
    real(PR), dimension(5), parameter   :: Cp = [1.95_PR, 1.39_PR, 2.4_PR, 4.18_PR, 1.58_PR]
    real(PR), dimension(3), parameter   :: dH = [-420._PR, -420._PR, -2440._PR]
    
end module mod_param