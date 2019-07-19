module ddpcm_lib
use ddcosmo
implicit none

contains

  subroutine ddpcm(phi, psi, esolv)
  implicit none
  real*8, intent(in) :: phi(ncav,nsph), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  

  end subroutine ddpcm

  subroutine bld_rprec()
  end subroutine bld_rprec
  
  subroutine rx()
  end subroutine rx
  
  subroutine dx()
  end subroutine dx

end module ddpcm_lib
