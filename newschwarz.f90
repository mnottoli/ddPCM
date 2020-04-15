module newschwarz
use ddcosmo
implicit none
! preconditioner array
real*8, allocatable :: nlprec(:,:,:)
real*8 :: p
contains

subroutine nddcosmo(phi,psi,esolv)
  ! new main (inside the module for clarity) 
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(out) :: esolv
  real*8, allocatable :: x(:,:), rhs(:,:), scr(:,:)
  integer :: isph, n_iter
  real*8 :: tol
  logical :: ok
  external :: hnorm, lx
  integer :: cr, c1, c2, c3
  integer :: i
  allocate(x(nbasis,nsph),rhs(nbasis,nsph),scr(ngrid,nsph))

  ! initialize the timer
  call system_clock(count_rate=cr)

  ! build the RHS 
  call wghpot(phi, scr)
  do isph = 1, nsph
    call intrhs(isph,scr(:,isph),rhs(:,isph)) 
  end do
  ! call prtsph('rhs of the ddCOSMO equation',nsph,0,rhs)
  deallocate(scr)


  do i = 1, 4
  p = dble(i)*0.1d0
  ! assemble and store the preconditioner
  call system_clock(count=c1)
  call build_nlprec()
  call system_clock(count=c2)

  ! solve ddcosmo
  call apply_nlprec(nbasis*nsph,rhs,x) 
  tol = 10.0d0**(-iconv)
  n_iter  = 200
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,x,n_iter,ok,nlx, &
    & apply_nlprec,hnorm)
  !write(6,*) p, n_iter 

  call system_clock(count=c3)
  !write(6,'(A15,F15.8)') 'precond time', dble(c2-c1)/dble(cr)
  !write(6,'(A15,F15.8)') 'solver time', dble(c3-c2)/dble(cr)

  ! compute the energy
  esolv = pt5*((eps - one)/eps)*sprod(nsph*nbasis,x,psi)
  write(6,'(F4.2,I5,3F15.8)') p, n_iter, dble(c3-c2)/dble(cr), dble(c2-c1)/dble(cr), esolv
  end do
  return
end subroutine nddcosmo 

subroutine nlx(n,x,y)
  ! perform new LX multiplication
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8, allocatable :: pot(:), basloc(:), dbsloc(:,:)
  real*8, allocatable :: vplm(:), vcos(:), vsin(:)
  integer :: isph, jsph, its, l1, m1, ind
  integer :: istatus

  allocate(pot(ngrid),vplm(nbasis),basloc(nbasis),dbsloc(3,nbasis), &
    & vcos(lmax+1),vsin(lmax+1),stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'Allocation failed in nlx.'
    stop
  end if

  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph,pot,basloc,dbsloc,vplm,vcos,vsin)
  do isph = 1, nsph
    call calcnv(isph,pot,x,basloc,dbsloc,vplm,vcos,vsin)
    call intrhs(isph,pot,y(:,isph))
    y(:,isph) = - y(:,isph)
    !call printmatrix(6,y(:,isph),1,nbasis)
  end do

  deallocate(pot,basloc,dbsloc,vplm,vcos,vsin)
  return 
end subroutine nlx

subroutine calcnv(isph,pot,sigma,basloc,dbsloc,vplm,vcos,vsin)
  implicit none
  integer, intent(in) :: isph
  real*8, intent(in) :: sigma(nbasis,nsph)
  real*8, intent(inout) :: basloc(nbasis), dbsloc(3,nbasis), vplm(nbasis), &
    & vcos(lmax + 1), vsin(lmax + 1), pot(ngrid)
  integer :: its, jsph, ij, l1, m1, ind, icomp, jcomp
  real*8 :: fac1, fac2, fac3, fac4, wij, vvij, tij, tt, res
  real*8 :: vij(3), sij(3), j(3,3), sg(3), c(3)

  pot = zero
  do its = 1, ngrid
    c = csph(:,isph) + rsph(isph)*grid(:,its)
    if (ui(its,isph).lt.one) then
      do ij = inl(isph), inl(isph+1) - 1
        jsph = nl(ij)
        ! compute geometrical variables
        vij = c - csph(:,jsph)
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

          ! compute the jacobian
          sg(1) = + grid(1,its)*(one - sij(1)*sij(1)) & 
              &   - grid(2,its)*sij(2)*sij(1) &
              &   - grid(3,its)*sij(3)*sij(1)

          sg(2) = - grid(1,its)*sij(1)*sij(2) &
              &   + grid(2,its)*(one - sij(2)*sij(2)) &
              &   - grid(3,its)*sij(3)*sij(2)

          sg(3) = - grid(1,its)*sij(1)*sij(3) & 
              &   - grid(2,its)*sij(2)*sij(3) &
              &   + grid(3,its)*(one - sij(3)*sij(3))
          sg = sg/vvij

          fac4 = dot_product(sij,grid(:,its))/(p*rsph(jsph)*tij)

          tt = one
          do l1 = 0, lmax
            fac1 = tt/(two*dble(l1) + one)
            fac2 = fac1*(dble(l1)*fac4 + one)
            ind = l1*l1 + l1 + 1
            do m1 = -l1, l1
              fac3 = fac1*dot_product(sg,dbsloc(:,ind+m1))/p
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
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the stored preconditioner 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph)
  do isph = 1, nsph
    call dgemm('n','n',nbasis,1,nbasis,one,nlprec(:,:,isph),nbasis, &
      & x(:,isph),nbasis,zero,y(:,isph),nbasis)
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
  allocate(nlprec(nbasis,nbasis,nsph),stat=istatus)
  nlprec = zero

  ! debug for matrix inversion
  ! allocate(nlprec_bk(nbasis,nbasis,nsph),res(nbasis,nbasis))

  ! allocate stuff for lapack matrix inversion
  allocate(ipiv(nbasis),work(nbasis),stat=istatus)

  ! dense contribution 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph,its,l1,m1,ind,lm,fac1,fac2,fac3)
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
          do lm = 1, nbasis
            nlprec(lm,ind + m1,isph) = nlprec(lm,ind + m1,isph) + &
              & fac3*basis(lm,its)
          end do 
        end do
      end do
    end do
  end do
  
  ! diagonal contribution
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph,l1,m1,ind,fac1)
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
  !$omp barrier

  ! write(8,*) '#', nbasis, 1
  ! call printmatrix(8,nlprec(:,:,1),nbasis,nbasis)
  ! debug
  ! nlprec_bk = nlprec

  ! invert 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph,ipiv,work,istatus)
  do isph = 1, nsph
    call dgetrf(nbasis,nbasis,nlprec(:,:,isph),nbasis,ipiv,istatus)
    if (istatus.ne.0) then
      write(6,*) 'LU failed with code', istatus
      stop
    end if
    call dgetri(nbasis,nlprec(:,:,isph),nbasis,ipiv,work,nbasis,istatus)
    if (istatus.ne.0) then
      write(6,*) 'Inversion failed'
      stop
    end if
  end do

  ! debug
  ! do isph = 1, nsph
  !   call dgemm('n','n',nbasis,nbasis,nbasis,one,nlprec(:,:,isph),nbasis, &
  !     & nlprec_bk(:,:,isph),nbasis,zero,res,nbasis)
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
