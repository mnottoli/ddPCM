module newschwarz
use ddcosmo
implicit none
! preconditioner array
real*8, allocatable :: nlprec(:,:,:)
real*8, parameter :: p = 1.0d0

contains

subroutine nddcosmo(phi,psi,esolv)
  ! new main (inside the module for clarity) 
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nylm,nsph), esolv
  real*8, allocatable :: x(:,:), rhs(:,:), scr(:,:)
  integer :: isph
  allocate(x(nylm,nsph),rhs(nylm,nsph),scr(ngrid,nsph))

  ! build the RHS 
  call wghpot(phi, scr)
  do isph = 1, nsph
    call intrhs(isph,scr(:,isph),rhs(:,isph)) 
  end do
  call prtsph('rhs of the ddCOSMO equation',nsph,0,rhs)
  deallocate(scr)

  ! assemble and store the preconditioner
  call build_nlprec()
  stop
  return
end subroutine nddcosmo 

subroutine nlx()
  ! perform new LX multiplication
end subroutine nlx

subroutine apply_nlprec()
  ! apply preconditioner
end subroutine apply_nlprec

subroutine build_nlprec()
  ! build preconditioner
  implicit none
  integer :: isph, its, l1, m1, ind, lm
  real*8 :: fac1, fac2
  integer :: istatus
  real*8, allocatable :: nlprec_bk(:,:,:), res(:,:)
  integer, allocatable :: ipiv(:)
  real*8, allocatable :: work(:)

  ! initialize the preconditioner
  allocate(nlprec(nylm,nylm,nsph))
  ! debug only
  allocate(nlprec_bk(nylm,nylm,nsph),res(nylm,nylm))
  nlprec = zero

  ! allocate stuff for lapack matrix inversion
  allocate(ipiv(nylm),work(nylm))

  ! dense contribution 
  do isph = 1, nsph
    fac1 = four*pi/(p*rsph(isph))
    do its = 1, ngrid
      fac1 = fac1*w(its)*(one - ui(its,isph))
      do l1 = 0, lmax
        ind = l1*l1 + l1 + 1
        do m1 = -l1, l1
          fac1 = fac1*basis(ind+m1,its)*dble(l1)/(two*dble(l1) + one)
          do lm = 1, nylm
            write(6,*) fac1*basis(lm,its)
            nlprec(lm,ind + m1,isph) = nlprec(lm,ind + m1,isph) + &
              & fac1*basis(lm,its)
          end do 
        end do
      end do
    end do
  end do

  ! diagonal contribution
  do isph = 1, nsph
    do l1 = 1, lmax
      fac1 = four*pi/(two*dble(l1) + one)
      do m1 = -l1, l1
        nlprec(ind + m1,ind + m1,isph) = nlprec(ind + m1,ind + m1,isph) &
          & + fac1
      end do
    end do 
  end do

  ! invert 
  nlprec_bk = nlprec
  do isph = 1, nsph
    call printmatrix(8,nlprec(:,:,isph),nylm,nylm)
    call dgetrf(nylm,nylm,nlprec(:,:,isph),nylm,ipiv,istatus)
    if (istatus.ne.0) then
      write(6,*) 'LU failed with code', istatus
      stop
    end if
    call dgetri(nylm,nlprec(:,:,isph),nylm,ipiv,work,nylm,istatus)
    if (istatus.ne.0) then
      write(6,*) 'Inversion failed'
      stop
    end if
  end do

  ! debug
  do isph = 1, nsph
    call dgemm('n','n',nylm,nylm,nylm,one,nlprec(:,:,isph),nylm, &
      & nlprec_bk(:,:,isph),nylm,zero,res,nylm)
    do lm = 1, nylm
      write(6,*) res(lm,lm)
    end do
  end do
  return
end subroutine build_nlprec


subroutine printmatrix(iout,a,m,n)
  implicit none 
  integer :: i, j
  integer, intent(in) :: iout, m, n
  real*8, intent(in) :: a(m,n)
  do i = 1, m
    do j = 1, n
      write(iout,'(F10.5 $)') a(i,j)
    end do 
    write(iout,*)
  end do
  return
end subroutine printmatrix

end module newschwarz
