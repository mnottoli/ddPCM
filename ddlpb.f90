subroutine ddlpb(phi,psi,gradphi,sigma,esolv)
use ddcosmo
implicit none
logical :: converged = .false.
integer :: iteration = 1
real*8, intent(inout) :: esolv
real*8, intent(inout) :: sigma(nylm,nsph)
real*8, intent(in) :: phi(ncav), gradphi(3,ncav)
real*8, intent(in) :: psi(nylm,nsph)
real*8, allocatable :: xr(:,:), xe(:,:), rhs_1(:,:), rhs_2(:,:)
real*8, allocatable :: g(:,:), f(:,:), g0(:), f0(:)
integer :: isph
integer :: i
external :: lx, ldm1x, hnorm, lstarx
logical :: ok = .false.
real*8 :: tol
integer :: n_iter
allocate(g(ngrid,nsph),f(ngrid,nsph))
allocate(g0(nylm),f0(nylm))
allocate(rhs_1(nylm,nsph),rhs_2(nylm,nsph))
allocate(xr(nylm,nsph),xe(nylm,nsph))
!
! Build the right hand side
!
call wghpot(phi,g)
! TODO: optimize wghpot_f
call wghpot_f(gradphi,f)

do isph = 1, nsph
  call intrhs(isph,g(:,isph),g0)
  call intrhs(isph,f(:,isph),f0)
  rhs_1(:,isph) = g0 + f0
  rhs_2(:,isph) = f0
end do

write(6,*) '#########'
do isph = 1, nsph
  do i = 1, nylm
    write(6,'(2F15.8)') rhs_1(i,isph), rhs_2(i,isph)
  end do
end do

tol = 10.0d0**(-iconv)

do while (.not.converged)
  ! solve ddcosmo step
  ! AXr = RHSr
  write(6,*) ok, n_iter
  call print_ddvector('rhs',rhs_1)
  !call convert_ddcosmo(1,rhs_1)
  n_iter = 200
  call jacobi_diis(nsph*nylm,iprint,ndiis,4,tol,rhs_1,xr,n_iter,&
&                  ok,lx,ldm1x,hnorm)
  call convert_ddcosmo(-1,xr)
  call print_ddvector('xr',xr)
  ! solve ddlpb step
  ! BXe = RHSe

  ! update the RHS

  ! check for convergency
  converged = .true.
  write(6,*) 'Iteration', iteration
  iteration = iteration + 1
end do
esolv = 1.2345d0
write(6,*) 'Done'
return
end subroutine ddlpb

subroutine convert_ddcosmo(direction,vector)
use ddcosmo
implicit none
integer :: isph, l, m, ind
integer, intent(in) :: direction
real*8, intent(inout) :: vector(nylm,nsph)
real*8 :: fac

do isph = 1, nsph
  do l = 0, lmax
    ind = l*l + l + 1
    fac = four*pi/(two*dble(l) + one) 
    if (direction.eq.1) fac = one/fac
    do m = -l, l
      vector(ind + m,isph) = fac*vector(ind + m,isph)
    end do
  end do
end do
return
end subroutine convert_ddcosmo

subroutine print_ddvector(label,vector)
use ddcosmo
implicit none
character(len=*) :: label
real*8 :: vector(nylm,nsph)
integer :: isph, lm

write(6,*) label
do isph = 1, nsph
  do lm = 1, nylm
    write(6,*) vector(lm,isph)
  end do
end do
return
end subroutine print_ddvector


