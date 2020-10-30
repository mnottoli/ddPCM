module ddpcm_lib
use ddcosmo, only: nbasis, nsph, ngrid, ncav, lmax, iconv, iprint, &
    & wghpot, intrhs, prtsph, zero, pt5, one, two, four, pi, basis, &
    & eps, csph, rsph, grid, w, ui, ndiis, sprod, ylmbas, facl, ddmkxi, &
    & ptcart
!use ddcosmo
implicit none

real*8, allocatable :: rx_prc(:,:,:)
real*8, allocatable :: rhs(:,:), phieps(:,:), xs(:,:)
real*8, allocatable :: s(:,:), y(:,:)
real*8, allocatable :: g(:,:)
logical :: dodiag

contains

  subroutine ddpcm_init()
  ! initialize ddpcm module by allocating the preconditioner and
  ! various arrays, then build the preconditioner
  implicit none
  integer :: istatus
  allocate(rx_prc(nbasis,nbasis,nsph),s(nbasis,nsph),y(nbasis,nsph), &
    & xs(nbasis,nsph),phieps(nbasis,nsph),stat=istatus)
  if (istatus.ne.0) write(6,*) 'ddpcm allocation failed'
  call mkprec
  end subroutine ddpcm_init

  subroutine ddpcm_finalize()
  ! deallocate various arrays
  implicit none
  integer :: istatus
  deallocate(rx_prc,s,y,xs,phieps,stat=istatus)
  if (istatus.ne.0) write(6,*) 'ddpcm deallocation failed'
  end subroutine ddpcm_finalize

  subroutine ddpcm(do_adjoint, phi, psi, esolv)
  ! main ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  logical, intent(in) :: do_adjoint
  real*8 :: tol, fac
  integer :: isph, n_iter
  logical :: ok
  external :: lx, ldm1x, lstarx, hnorm
  
  allocate(rhs(nbasis,nsph),g(ngrid,nsph))
  tol = 10.0d0**(-iconv)

  ! build RHS
  g = zero
  xs = zero
  call wghpot(phi,g)
  do isph = 1, nsph
    call intrhs(isph,g(:,isph),xs(:,isph))
  end do

  ! rinf rhs
  dodiag = .true.
  call rinfx(nbasis*nsph,xs,rhs)

  ! solve the ddpcm linear system
  n_iter = 200
  dodiag = .false.
  phieps = xs
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,phieps,n_iter, &
      & ok,rx,apply_rx_prec,hnorm)
  write(6,*) 'ddpcm step iterations:', n_iter

  ! solve the ddcosmo linear system
  n_iter = 200
  dodiag = .false.
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  write(6,*) 'ddcosmo step iterations:', n_iter

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

  if (do_adjoint) then

    ! solve ddcosmo adjoint system
    n_iter = 200
    dodiag = .false.
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,psi,s,n_iter, &
      & ok,lstarx,ldm1x,hnorm)
    call prtsph('S',nsph,0,s)

    ! solve ddpcm adjoint system
    n_iter = 200
    dodiag = .false.
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,s,y,n_iter, &
      & ok,rstarx,apply_rstarx_prec,hnorm)
    call prtsph('Y',nsph,0,y)

    ! recover effect of Rinf^*
    fac = two*pi*(one - (eps + one)/(eps - one))
    y = s + fac*y
    call prtsph('adjoint ddpcm solution',nsph,0,y)
  end if
  end subroutine ddpcm


  subroutine mkprec
  ! Assemble the diagonal blocks of the Reps matrix
  ! then invert them to build the preconditioner
  implicit none
  integer :: isph, lm, ind, l1, m1, ind1, its, istatus
  real*8  :: f, f1
  integer, allocatable :: ipiv(:)
  real*8,  allocatable :: work(:)

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

  subroutine rstarx(n,x,y)
  ! Computes Y = Reps X =
  ! = (2*pi*(eps + 1)/(eps - 1) - D^*) X 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dstarx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rstarx

  subroutine dstarx(n,x,y)
  ! Computes Y = D^* X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
  real*8 :: c(3), vji(3), sji(3)
  real*8 :: vvji, tji, fourpi, tt, f, f1
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
  !$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vji, &
  !$omp vvji,tji,sji,tt,l,ind,f,m,vts,c) &
  !$omp shared(nsph,ngrid,ui,csph,rsph,grid,facl, &
  !$omp lmax,fourpi,dodiag,x,y,basis,w,nbasis)
  do isph = 1, nsph
    do jsph = 1, nsph
      if (jsph.ne.isph) then
        do its = 1, ngrid
          if (ui(its,jsph).gt.zero) then
            ! build the geometrical variables
            vji = csph(:,jsph) + rsph(jsph)*grid(:,its) - csph(:,isph)
            vvji = sqrt(dot_product(vji,vji))
            tji = vvji/rsph(isph)
            sji = vji/vvji
            ! build the local basis
            call ylmbas(sji,basloc,vplm,vcos,vsin)
            tt = w(its)*ui(its,jsph)*dot_product(basis(:,its),x(:,jsph))/tji
            do l = 0, lmax
              ind = l*l + l + 1
              f = dble(l)*tt/facl(ind)
              do m = -l, l
                y(ind+m,isph) = y(ind+m,isph) + f*basloc(ind+m)
              end do
              tt = tt/tji
            end do
          end if
        end do
      else if (dodiag) then
        do its = 1, ngrid
          f = pt5*w(its)*ui(its,jsph)*dot_product(basis(:,its),x(:,jsph))
          do ind = 1, nbasis
            y(ind,isph) = y(ind,isph) - f*basis(ind,its)/facl(ind)
          end do
        end do
      end if
    end do
  end do

  deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: deallocation failed !'
    stop
  end if
  end subroutine dstarx

  subroutine apply_rstarx_prec(n,x,y)
  ! apply the block diagonal preconditioner
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the transposed preconditioner 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph)
  do isph = 1, nsph
    call dgemm('t','n',nbasis,1,nbasis,one,rx_prc(:,:,isph),nbasis, &
      & x(:,isph),nbasis,zero,y(:,isph),nbasis)
  end do
  end subroutine apply_rstarx_prec

end module ddpcm_lib
