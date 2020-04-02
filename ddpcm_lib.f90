module ddpcm_lib
use ddcosmo, only: nbasis, nsph, ngrid, ncav, lmax, iconv, iprint, &
    & wghpot, intrhs, prtsph, zero, pt5, one, two, four, pi, basis, &
    & eps, csph, rsph, grid, w, ui, ndiis, sprod, ylmbas
!use ddcosmo
implicit none

real*8, allocatable :: rx_prc(:,:,:)
real*8, allocatable :: rhs(:,:), phieps(:,:), xs(:,:)
real*8, allocatable :: g(:,:)
logical :: dodiag

contains

  subroutine ddpcm(phi, psi, esolv)
  ! main ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  real*8 :: tol
  integer :: isph, n_iter
  logical :: ok
  external :: lx, ldm1x, hnorm
  
  allocate(rx_prc(nbasis,nbasis,nsph))
  allocate(rhs(nbasis,nsph),phieps(nbasis,nsph),xs(nbasis,nsph))
  allocate(g(ngrid,nsph))
  tol = 10.0d0**(-iconv)

  ! build the preconditioner
  call mkprecsvd

  ! build the RHS
  !write(6,*) 'pot', ncav
  !do isph = 1, ncav
  !  write(6,*) phi(isph)
  !end do
  g = zero
  xs = zero
  call wghpot(phi,g)
  do isph = 1, nsph
    call intrhs(isph,g(:,isph),xs(:,isph))
  end do

  ! call prtsph('phi',nsph,0,xs)
  ! call prtsph('psi',nsph,0,psi)

  ! rinf rhs
  dodiag = .true.
  call rinfx(nbasis*nsph,xs,rhs)
  ! call prtsph('rhs',nsph,0,rhs)

  ! solve the ddpcm linear system
  n_iter = 200
  dodiag = .false.
  phieps = xs
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,phieps,n_iter, &
      & ok,rx,apply_rx_prec,hnorm)
  write(6,*) 'ddpcm step iterations:', n_iter
  ! call prtsph('phie',nsph,0,phieps)

  ! solve the ddcosmo linear system
  n_iter = 200
  dodiag = .false.
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  write(6,*) 'ddcosmo step iterations:', n_iter
  ! call prtsph('x',nsph,0,xs)

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

  return
  end subroutine ddpcm

  subroutine mkprecsvd
  ! Assemble the diagonal blocks of the Reps matrix
  ! then invert them to build the preconditioner
  implicit none
  integer :: isph, lm, ind, l1, m1, ind1, its, istatus, lm1
  integer :: lwork
  real*8  :: f, f1
  real*8,  allocatable :: s(:), u(:,:), vt(:,:), work(:)
  real*8, allocatable :: scr(:,:)

  lwork = 6*nbasis
  allocate(work(lwork),s(nbasis),u(nbasis,nbasis), &
    & vt(nbasis,nbasis),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) 'mkprecsvd : allocation failed !'
    stop
  endif
  allocate(scr(nbasis,nbasis))

  rx_prc(:,:,:) = zero

  ! off diagonal part
  do isph = 1, nsph
    do its = 1, ngrid
      f = two*pi*ui(its,isph)*w(its)
      do l1 = 0, lmax
        ind1 = l1*l1 + l1 + 1
        do m1 = -l1, l1
          f1 = f*basis(ind1 + m1,its)/(two*dble(l1) + one)
          do lm = 1, nbasis
            rx_prc(lm,ind1 + m1,isph) = rx_prc(lm,ind1 + m1,isph) + &
                & f1*basis(lm,its)
          end do
        end do
      end do
    end do
  end do

  ! add diagonal
  f = two*pi*(eps + one)/(eps - one)
  do isph = 1, nsph
    do lm = 1, nbasis
      rx_prc(lm,lm,isph) = rx_prc(lm,lm,isph) + f
    end do
  end do

  ! do a pseudoinversion using singular value decomposition
  do isph = 1, nsph
    scr = rx_prc(:,:,isph)
    call dgesvd('A','A',nbasis,nbasis,rx_prc(1,1,isph),nbasis, &
      & s,u,nbasis,vt,nbasis,work,lwork,istatus)
    rx_prc(:,:,isph) = zero
    do lm = 1, nbasis
      f = one/s(lm)
      do lm1 = 1, nbasis
        vt(lm,lm1) = vt(lm,lm1)*f
      end do
    end do
    call dgemm('t','t',nbasis,nbasis,nbasis,one,vt,nbasis,u,nbasis, &
      & one,rx_prc(1,1,isph),nbasis)
    call dgemm('n','n',nbasis,nbasis,nbasis,one,scr,nbasis, &
      & rx_prc(1,1,isph),nbasis,zero,u,nbasis)
    write(6,*) 'isph', isph
    do lm = 1, nbasis
      do lm1 = 1, nbasis
        write(6,'(F12.6$)') u(lm,lm1)
      end do 
      write(6,*)
    end do
  end do
  deallocate(work)
  deallocate(scr)
  endsubroutine mkprecsvd

  subroutine mkprec
  ! Assemble the diagonal blocks of the Reps matrix
  ! then invert them to build the preconditioner
  implicit none
  integer :: isph, lm, ind, l1, m1, ind1, its, istatus
  real*8  :: f, f1
  integer, allocatable :: ipiv(:)
  real*8, allocatable :: work(:)
  real*8, allocatable :: scr(:,:), u(:,:)
  integer :: lm1
  allocate(scr(nbasis,nbasis),u(nbasis,nbasis))
  allocate(ipiv(nbasis),work(nbasis*nbasis),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) 'mkprec : allocation failed !'
    stop
  endif

  rx_prc(:,:,:) = zero

  ! off diagonal part
  do isph = 1, nsph
    do its = 1, ngrid
      f = two*pi*ui(its,isph)*w(its)
      do l1 = 0, lmax
        ind1 = l1*l1 + l1 + 1
        do m1 = -l1, l1
          f1 = f*basis(ind1 + m1,its)/(two*dble(l1) + one)
          do lm = 1, nbasis
            rx_prc(lm,ind1 + m1,isph) = rx_prc(lm,ind1 + m1,isph) + &
                & f1*basis(lm,its)
          end do
        end do
      end do
    end do
  end do

  ! add diagonal
  f = two*pi*(eps + one)/(eps - one)
  do isph = 1, nsph
    do lm = 1, nbasis
      rx_prc(lm,lm,isph) = rx_prc(lm,lm,isph) + f
    end do
  end do

  ! invert the blocks
  do isph = 1, nsph
    scr = rx_prc(:,:,isph)
    call DGETRF(nbasis,nbasis,rx_prc(:,:,isph),nbasis,ipiv,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'LU failed in mkprc'
      stop
    end if
    call DGETRI(nbasis,rx_prc(:,:,isph),nbasis,ipiv,work, &
        & nbasis*nbasis,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'Inversion failed in mkprc'
      stop
    end if
    call dgemm('n','n',nbasis,nbasis,nbasis,one,scr,nbasis, &
      & rx_prc(1,1,isph),nbasis,zero,u,nbasis)
    write(6,*) 'isph', isph
    do lm = 1, nbasis
      do lm1 = 1, nbasis
        write(6,'(F12.6$)') u(lm,lm1)
      end do 
      write(6,*)
    end do
  end do

  deallocate (work,ipiv,stat=istatus)
  if (istatus.ne.0) then
    write(*,*) 'mkprec : deallocation failed !'
    stop
  end if
  return
  endsubroutine mkprec


  subroutine rx(n,x,y)
  ! Computes Y = Reps X =
  ! = (2*pi*(eps + 1)/(eps - 1) - D) X 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rx


  subroutine rinfx(n,x,y)
  ! Computes Y = Rinf X = 
  ! = (2*pi -  D) X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi
    y = y + fac*x
  end if
  end subroutine rinfx

  
  subroutine dx(n,x,y)
  ! Computes Y = D X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
  real*8 :: c(3), vij(3), sij(3)
  real*8 :: vvij, tij, fourpi, tt, f, f1
  integer :: its, isph, jsph, l, m, ind, lm, istatus
  
  allocate(vts(ngrid),vplm(nbasis),basloc(nbasis),vcos(lmax + 1), &
      & vsin(lmax + 1),stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: allocation failed !'
    stop
  end if
  y = zero
  fourpi = four*pi

  !$omp parallel do default(none) schedule(dynamic) &
  !$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vij, &
  !$omp vvij,tij,sij,tt,l,ind,f,m,vts,c) &
  !$omp shared(nsph,ngrid,ui,csph,rsph,grid, &
  !$omp lmax,fourpi,dodiag,x,y,basis)
  do isph = 1, nsph
    ! compute the "potential" from the other spheres
    ! at the exposed Lebedv points of the i-th sphere 
    vts = zero
    do its = 1, ngrid
      if (ui(its,isph).gt.zero) then
        c = csph(:,isph) + rsph(isph)*grid(:,its)
        do jsph = 1, nsph
          if (jsph.ne.isph) then
            ! build the geometrical variables
            vij = c - csph(:,jsph)
            vvij = sqrt(dot_product(vij,vij))
            tij = vvij/rsph(jsph)
            sij = vij/vvij 
            ! build the local basis
            call ylmbas(sij,basloc,vplm,vcos,vsin)
            ! with all the required stuff, finally compute
            ! the "potential" at the point 
            tt = one/tij 
            do l = 0, lmax
              ind = l*l + l + 1
              f = fourpi*dble(l)/(two*dble(l) + one)*tt
              do m = -l, l
                vts(its) = vts(its) + f*x(ind + m,jsph) * &
                    & basloc(ind + m)
              end do
              tt = tt/tij
            end do
          else if (dodiag) then
            ! add the diagonal contribution
            do l = 0, lmax
              ind = l*l + l + 1
              f = (two*dble(l) + one)/fourpi
              do m = -l, l
                vts(its) = vts(its) - pt5*x(ind + m,isph) * &
                    & basis(ind + m,its)/f
              end do
            end do
          end if 
        end do
        vts(its) = ui(its,isph)*vts(its) 
      end if
    end do
    ! now integrate the potential to get its modal representation
    call intrhs(isph,vts,y(:,isph))
  end do

  deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: deallocation failed !'
    stop
  end if
  end subroutine dx

  subroutine apply_rx_prec(n,x,y)
  ! apply the block diagonal preconditioner
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the stored preconditioner 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph)
  do isph = 1, nsph
    call dgemm('n','n',nbasis,1,nbasis,one,rx_prc(:,:,isph),nbasis, &
      & x(:,isph),nbasis,zero,y(:,isph),nbasis)
  end do
  end subroutine apply_rx_prec

end module ddpcm_lib
