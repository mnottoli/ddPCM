subroutine simpleddcosmo(phi,psi,esolv)
use ddcosmo
! simple implemantation of ddcosmo driver for faster debug
implicit none
real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
real*8, intent(inout) :: esolv
real*8 :: tol
integer :: isph, n_iter
logical :: ok
external :: lx, ldm1x, hnorm
real*8, allocatable :: rhs(:,:), x(:,:), g(:,:)
logical :: dodiag

! setup
allocate(rhs(nbasis,nsph),x(nbasis,nsph),g(ngrid,nsph))
tol = 10.0d0**(-iconv)

write(6,*) 'here'

! build rhs
call wghpot(phi,g)
do isph = 1, nsph
  call intrhs(isph,g(:,isph),rhs(:,isph))
end do

write(6,*) 'here'

! guess
do isph = 1, nsph
  x(:,isph) = facl(:)*rhs(:,isph)
enddo

! solve ddcosmo
n_iter = 200
dodiag = .false.
call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,x,n_iter, &
    & ok,lx,ldm1x,hnorm)
write(6,*) 'ddcosmo iterations:', n_iter

! compute energy
esolv = pt5*((eps - one)/eps)*sprod(nsph*nbasis,x,psi)

deallocate(rhs,x,g)
end subroutine simpleddcosmo
