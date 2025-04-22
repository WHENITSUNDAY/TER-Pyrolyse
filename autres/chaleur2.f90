program equa_chal

  use euler_explicite
  use densite
  
  implicit none

  ! Paramètres de l'algorithme
  integer, parameter :: Nx = 1000
  real(pr), parameter :: Tf = 1., T0 = 293., L = 1.
  real(pr), parameter :: dx = L/real(Nx, pr)
  real(pr), parameter :: D = 1.
  real(pr), parameter :: dt = dx**2 / (4._pr*D)
  integer, parameter :: Nt = int(Tf/dt, 4)
  real(pr), dimension(:, :), allocatable :: T, rb, rc, rg, rl, rv
  integer :: x, k

  ! Variables utiles à la lecture du code
  real(pr) :: lambda0, lambda1, rcp, Qr

  ! Allocation
  allocate (T(0:Nx, 0:Nt))
  allocate (rb(0:Nx, 0:Nt))
  allocate (rc(0:Nx, 0:Nt))
  allocate (rg(0:Nx, 0:Nt))
  allocate (rl(0:Nx, 0:Nt))
  allocate (rv(0:Nx, 0:Nt))
  
  ! Initialisation
  do x = 0, Nx
     T(x, 0) = T0
     rb(x, 0) = rb0
     rc(x, 0) = 0
     rg(x, 0) = 0
     rl(x, 0) = rl0
     rv(x, 0) = 0
  end do
  
  ! Calcul
  do x = 1, Nx - 1
     do k = 0, Nt
        lambda0 = calcul_ls (lc, calcul_lb (rl(x, k), rc(x, k), rb(x, k)), rc(x, k), rb(x, k))
        lambda1 = calcul_ls (lc, calcul_lb (rl(x+1, k), rc(x+1, k), rb(x+1, k)), rc(x+1, k), rb(x+1, k))
        rcp = calcul_rcp (rb(x, k), rc(x, k), rl(x, k))
        Qr = calcul_Qr (rb(x, k), rl(x, k), T(x, k))
        
        T(x, k+1) = (dx**2)*T(x, k) + dt*(lambda1-lambda0)*(T(x+1, k)-T(x, k)) + dt*lambda0*(T(x+1, k)+T(x-1, k)) + (dx**2)*dt*Qr
        T(x, k+1) = T(x, k+1) / ((dx**2)*rcp + 2*dt*lambda0)

     end do
  end do
  
  
  ! Désallocation
  deallocate (T)
  deallocate (rb)
  deallocate (rc)
  deallocate (rg)
  deallocate (rl)
  deallocate (rv)

  
contains

  function calcul_k (i, T) result (k)
    ! Déclaration des arguments
    integer, intent(in) :: i
    real(pr), intent(in) :: T
    real(pr) :: k

    ! Calcul de ki
    select case (i)
    case(1)
       k = A1 * exp (-E1/(R*T))
    case(2)
       k = A2 * exp (-E2/(R*T))
    case (3)
       k = A3 * exp (-E3/(R*T))
    case default
       error stop "Y a que 3 phases hein..."
    end select
  end function calcul_k

  function calcul_lb (rl, rc, rb) result (lb)
    ! Déclaration des arguments
    real(pr), intent(in) :: rl, rc, rb
    real(pr) :: lb

    ! Déclaration de variables locales
    real(pr) :: x
    
    ! Calcul de lambda
    x = rl / (rl + rc + rb)
    lb = 0.166 + 0.369 * x
  end function calcul_lb
  
  function calcul_ls (lc, lb, rc, rb) result (ls)
    ! Déclaration des arguments
    real(pr), intent(in) :: lc, lb, rc, rb
    real(pr) :: ls

    ! Déclaration de variables locales
    real(pr) :: eta
    
    ! Calcul de lambda
    eta = rc / (rc + rb)
    ls = eta * lc + (1 - eta) * lb
  end function calcul_ls

  function calcul_Qr (rb, rl, T) result (Qr)
    ! Déclaration des arguments
    real(pr), intent(in) :: rb, rl, T
    real(pr) :: Qr

    ! Déclaration de variables locales
    real(pr) :: k1, k2, k3

    ! Initialisation des variables
    k1 = calcul_k (1, T)
    k2 = calcul_k (2, T)
    k3 = calcul_k (3, T)

    ! Calcul de Qr
    Qr = k1*rb * (dh10 + (Cc - Cb)*(T - T0))
    Qr = Qr + k2 * rb * (dh20 + (Cg - Cb) * (T - T0))
    Qr = Qr + k3 * rl * (dh30 + (Cv - Cl) * (T - T0))
  end function calcul_Qr

  function calcul_rcp (rb, rc, rl) result (rcp)
    ! Déclaration des variables
    real(pr), intent(in) :: rb, rc, rl
    real(pr) :: rcp

    ! Calcul de rho * C
    rcp = rb * Cb + rc * Cc + rl * Cl
  end function calcul_rcp

end program equa_chal
