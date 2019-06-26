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
  real*8, intent(in) :: phi(ncav), psi(nylm,nsph)
  real*8, intent(out) :: esolv
  real*8, allocatable :: x(:,:), rhs(:,:), scr(:,:)
  integer :: isph, n_iter
  real*8 :: tol
  logical :: ok
  external :: hnorm, lx
  allocate(x(nylm,nsph),rhs(nylm,nsph),scr(ngrid,nsph))

  ! build the RHS 
  call wghpot(phi, scr)
  do isph = 1, nsph
    call intrhs(isph,scr(:,isph),rhs(:,isph)) 
  end do
  ! call prtsph('rhs of the ddCOSMO equation',nsph,0,rhs)
  deallocate(scr)

  ! assemble and store the preconditioner
  call build_nlprec()

  ! solve ddcosmo 
  tol     = 10.0d0**(-iconv)
  n_iter  = 200
  call jacobi_diis(nsph*nylm,iprint,ndiis,4,tol,rhs,x,n_iter,ok,nlx, &
    & apply_nlprec,hnorm)
  write(6,*) n_iter 

  ! compute the energy
  esolv = pt5*((eps - one)/eps)*sprod(nsph*nylm,x,psi)
  return
end subroutine nddcosmo 

subroutine nlx(n,x,y)
  ! perform new LX multiplication
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nylm,nsph)
  real*8, intent(inout) :: y(nylm,nsph)
  real*8, allocatable :: pot(:), basloc(:), dbsloc(:,:)
  real*8, allocatable :: vplm(:), vcos(:), vsin(:)
  integer :: isph, jsph, its, l1, m1, ind
  integer :: istatus

  allocate(pot(ngrid),vplm(nylm),basloc(nylm),dbsloc(3,nylm), &
    & vcos(lmax+1),vsin(lmax+1),stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'Allocation failed in nlx.'
    stop
  end if

  do isph = 1, nsph
    call calcnv(isph,pot,x,basloc,dbsloc,vplm,vcos,vsin)
    call intrhs(isph,pot,y(:,isph))
    y(:,isph) = - y(:,isph)
    !call printmatrix(6,y(:,isph),1,nylm)
  end do

  deallocate(pot,basloc,dbsloc,vplm,vcos,vsin)
  return 
end subroutine nlx

subroutine calcnv(isph,pot,sigma,basloc,dbsloc,vplm,vcos,vsin)
  implicit none
  integer, intent(in) :: isph
  real*8, intent(in) :: sigma(nylm,nsph)
  real*8, intent(inout) :: basloc(nylm), dbsloc(3,nylm), vplm(nylm), &
    & vcos(lmax + 1), vsin(lmax + 1), pot(ngrid)
  integer :: its, jsph, ij, l1, m1, ind, icomp, jcomp
  real*8 :: fac1, fac2, fac3, wij, vvij, tij, tt, res
  real*8 :: vij(3), sij(3), j(3,3), sg(3)

  pot = zero
  do its = 1, ngrid
    if (ui(its,isph).lt.one) then
      do ij = inl(isph),inl(isph+1)-1
        jsph = nl(ij)
        ! compute geometrical variables
        vij = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij = vvij/rsph(jsph)
        res = zero
        if (tij.lt.one) then
          sij = vij/vvij
          ! compute ddcosmo wij
          wij = fsw(tij,se,eta)
          if (fi(its,isph).gt.one) then
            wij = wij/fi(its,isph)
          end if 
          ! assemble local basis and gradient
          call dbasis(sij,basloc,dbsloc,vplm,vcos,vsin)
          ! build the jacobian Ji (sij) 
          j = zero 
          j(1,1) = one
          j(2,2) = one
          j(3,3) = one
          do icomp = 1, 3 
            do jcomp = 1, 3
              j(icomp,jcomp) = (j(icomp,jcomp) - &
                & sij(icomp)*sij(jcomp))/vvij
            end do
          end do
          ! TODO: build directly the vector
          sg(1) = grid(1,its)*j(1,1) + grid(2,its)*j(2,1) + grid(3,its)*j(3,1) 
          sg(2) = grid(1,its)*j(1,2) + grid(2,its)*j(2,2) + grid(3,its)*j(3,2) 
          sg(3) = grid(1,its)*j(1,3) + grid(2,its)*j(2,3) + grid(3,its)*j(3,3) 
          ! assemble the rhs contraction
          tt = one
          do l1 = 0, lmax
            fac1 = tt/(two*dble(l1) + one)
            fac2 = fac1*(dble(l1)/(p*rsph(jsph)*tij)* &
              & (sij(1)*grid(1,its) + sij(2)*grid(2,its) + &
              & sij(3)*grid(3,its)) + one)
            ind = l1*l1 + l1 + 1
            do m1 = -l1, l1
              fac3 = fac1*(sg(1)*dbsloc(1,ind+m1) + & 
                & sg(2)*dbsloc(2,ind+m1) + sg(3)*dbsloc(3,ind+m1))/p
              res = res + sigma(ind + m1,jsph)*(fac2*basloc(ind+m1) + &
                & fac3)
            end do
            tt = tt*tij
          end do
          ! accumulate the result
          pot(its) = pot(its) + four*pi*wij*res
        end if
      end do
    end if
  end do
end subroutine calcnv

subroutine apply_nlprec(n,x,y)
  ! apply preconditioner
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nylm,nsph)
  real*8, intent(inout) :: y(nylm,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the stored preconditioner 
  do isph = 1, nsph
    call dgemm('n','n',nylm,1,nylm,one,nlprec(:,:,isph),nylm, &
      & x(:,isph),nylm,zero,y(:,isph),nylm)
  end do 
  end subroutine apply_nlprec

subroutine build_nlprec()
  ! build preconditioner
  implicit none
  integer :: isph, its, l1, m1, ind, lm
  real*8 :: fac1, fac2, fac3
  integer :: istatus
  real*8, allocatable :: nlprec_bk(:,:,:), res(:,:)
  integer, allocatable :: ipiv(:)
  real*8, allocatable :: work(:)

  ! initialize the preconditioner
  allocate(nlprec(nylm,nylm,nsph),stat=istatus)
  nlprec = zero

  ! debug for matrix inversion
  ! allocate(nlprec_bk(nylm,nylm,nsph),res(nylm,nylm))

  ! allocate stuff for lapack matrix inversion
  allocate(ipiv(nylm),work(nylm),stat=istatus)

  ! dense contribution 
  do isph = 1, nsph
    fac1 = four*pi/(p*rsph(isph))
    do its = 1, ngrid
      fac2 = fac1*w(its)*(one - ui(its,isph))
      !write(6,*) fac2, fac1, w(its), ui(its,isph)
      do l1 = 0, lmax
        ind = l1*l1 + l1 + 1
        do m1 = -l1, l1
          fac3 = fac2*basis(ind+m1,its)*dble(l1)/(two*dble(l1) + one)
          !write(6,*) fac3, fac2, basis(ind+m1,its), dble(l1)/(two*dble(l1)+one)
          do lm = 1, nylm
            nlprec(lm,ind + m1,isph) = nlprec(lm,ind + m1,isph) + &
              & fac3*basis(lm,its)
          end do 
        end do
      end do
    end do
  end do
  
  ! diagonal contribution
  do isph = 1, nsph
    do l1 = 0, lmax
      fac1 = four*pi/(two*dble(l1) + one)
      ind = l1*l1 + l1 + 1
      do m1 = -l1, l1
        nlprec(ind + m1,ind + m1,isph) = nlprec(ind + m1,ind + m1,isph) &
          & + fac1
      end do
    end do 
  end do

  ! write(8,*) '#', nylm, 1
  ! call printmatrix(8,nlprec(:,:,1),nylm,nylm)

  ! invert 
  nlprec_bk = nlprec
  do isph = 1, nsph
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
  ! do isph = 1, nsph
  !   call dgemm('n','n',nylm,nylm,nylm,one,nlprec(:,:,isph),nylm, &
  !     & nlprec_bk(:,:,isph),nylm,zero,res,nylm)
  ! end do

  deallocate(ipiv,work)
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
