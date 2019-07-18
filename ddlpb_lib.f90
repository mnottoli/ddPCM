module ddlpb_lib
use ddcosmo
implicit none
logical :: firstiniter, firstoutiter, firstcosmo, firsthsp
logical :: do_diag
integer :: matAB
! just hardcode these
integer, parameter :: lmax0 = 6
integer, parameter :: nbasis0 = 49
real*8, parameter :: epsp = 1.0d0
! as they are in Chaoyu's code 
real*8, allocatable :: wij(:,:), bas_sij(:,:,:), &
    & fac_cosmo_ij(:,:,:), fac_hsp_ij(:,:,:)
real*8, allocatable :: coefvec(:,:,:), Pchi(:,:,:), &
    & Qmat(:,:,:), coefY(:,:,:)
real*8, allocatable :: SI_ri(:,:), DI_ri(:,:), SK_ri(:,:), &
    & DK_ri(:,:), termimat(:,:)
real*8  :: kappa, rp, tol_gmres, n_iter_gmres

contains

  subroutine ddlpb(phi,q,psi,gradphi,sigma,esolv)
  ! main ddlpb
  implicit none
  logical :: converged = .false.
  integer :: iteration = 1
  real*8, intent(inout) :: esolv
  real*8 :: inc, old_esolv
  real*8, intent(inout) :: sigma(nylm,nsph)
  real*8, intent(in) :: phi(ncav), gradphi(3,ncav)
  real*8, intent(in) :: psi(nylm,nsph), q(nsph)
  real*8, allocatable :: xr(:,:), xe(:,:), rhs_r(:,:), rhs_e(:,:), &
      & rhs_r_init(:,:), rhs_e_init(:,:)
  real*8, allocatable :: g(:,:), f(:,:), g0(:), f0(:)
  integer :: isph
  integer :: i
  external :: lx, ldm1x, hnorm, lstarx
  logical :: ok = .false.
  real*8 :: tol
  integer :: n_iter
  integer :: its

  call ddlpb_init

  allocate(g(ngrid,nsph),f(ngrid,nsph))
  allocate(g0(nylm),f0(nylm))
  allocate(rhs_r(nylm,nsph),rhs_e(nylm,nsph))
  allocate(rhs_r_init(nylm,nsph),rhs_e_init(nylm,nsph))
  allocate(xr(nylm,nsph),xe(nylm,nsph))

  !do i = 1, ncav
  !  write(6,'(3F15.8)') phi(i), gradphi(:,i)
  !end do
  !stop

  ! Build the right hand side
  ! do i = 1, ncav
  !   write(6,'(4F20.10)') phi(i), gradphi(:,i)
  ! end do

  call wghpot(phi,g)
  call wghpot_f(gradphi,f)

  ! do isph = 1, nsph
  !   do i = 1, ngrid
  !     write(6,'(2F20.10)') g(i,isph), f(i,isph)
  !   end do
  ! end do
  
  do isph = 1, nsph
    call intrhs(isph,g(:,isph),g0)
    call intrhs(isph,f(:,isph),f0)
    ! rhs 
    rhs_r_init(:,isph) = g0 + f0
    rhs_e_init(:,isph) = f0
  end do

  rhs_r = rhs_r_init
  rhs_e = rhs_e_init

  tol = 10.0d0**(-iconv)

  ! TODO: remove this
  firstcosmo = .true.
  firsthsp = .true.
  firstiniter = .true.
  firstoutiter = .true.

  do while (.not.converged)

    ! solve the ddcosmo step
    ! A X_r = RHS_r 
    n_iter = 200
    call jacobi_diis(nsph*nylm,iprint,ndiis,4,tol,rhs_r,xr,n_iter,&
      & ok,lx,ldm1x,hnorm)
    call convert_ddcosmo(1,xr)
    ! call print_ddvector('xr',xr)
  
    ! solve ddlpb step
    ! B X_e = RHS_e
    call lpb_hsp(rhs_e,xe)
    ! call print_ddvector('xe',xe)
  
    ! update the RHS
    ! / RHS_r \ = / g + f \ - / c1 c2 \ / X_r \
    ! \ RHS_e /   \ f     /   \ c1 c2 / \ X_e /
    call update_rhs(rhs_r_init,rhs_e_init,rhs_r,rhs_e,xr,xe)
    ! call print_ddvector('rhs_r',rhs_r)
    ! call print_ddvector('rhs_e',rhs_e)

    ! compute energy
    !esolv = pt5*sprod(nsph*nylm,xr,psi) ???
    esolv = zero
    do isph = 1, nsph
      esolv = esolv + pt5*q(isph)*Xr(1,isph)*(one/(two*sqrt(pi)))
    end do

    ! check for convergency
    inc = abs(esolv - old_esolv)/abs(esolv)
    old_esolv = esolv
    if ((iteration.gt.1) .and. (inc.lt.tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) iteration, esolv, inc
    iteration = iteration + 1

    ! to be removed
    firstcosmo = .false.
    firsthsp = .false.
    firstiniter = .false.
    firstoutiter = .false.
  end do

  return
  end subroutine ddlpb


  subroutine ddlpb_init
  use bessel
  implicit none
  integer :: istatus, isph, NM
  allocate(SI_ri(0:lmax,nsph),DI_ri(0:lmax,nsph),SK_ri(0:lmax,nsph), &
      & DK_ri(0:lmax,nsph),termimat(0:lmax,nsph),stat=istatus)
  if (istatus.ne.0) then
    write(*,*)'ddinit : [1] allocation failed !'
    stop
  end if
  do isph = 1, nsph
    call SPHI_bessel(lmax,rsph(isph)*kappa,NM,SI_ri(:,isph), &
        & DI_ri(:,isph))
    call SPHK_bessel(lmax,rsph(isph)*kappa,NM,SK_ri(:,isph), &
        & DK_ri(:,isph))
  end do
  return
  end subroutine ddlpb_init


  subroutine wghpot_f( gradphi, f )
  use bessel
  implicit none
  real*8, dimension(3, ncav),       intent(in)  :: gradphi
  real*8, dimension(ngrid,nsph), intent(out) :: f

  integer :: isph, ig, ic, ind, ind0, jg, l, m, jsph
  real*8 :: nderphi, sumSijn, rijn, coef_Ylm, sumSijn_pre, termi, &
      & termk, term
  real*8, dimension(3) :: sijn, vij
  real*8, allocatable :: SK_rijn(:), DK_rijn(:)

  integer :: l0, m0, NM, kep, istatus
  real*8, dimension(nylm, nsph) :: c0
  real*8, dimension(0:lmax, nsph) :: coef_bessel
  real*8, allocatable :: vplm(:), basloc(:), vcos(:), vsin(:)

  ! initialize
  allocate(vplm(nylm),basloc(nylm),vcos(lmax+1),vsin(lmax+1))
  allocate(SK_rijn(0:lmax0),DK_rijn(0:lmax0))
  ic = 0 ; f(:,:)=0.d0

  ! compute c0
  do isph = 1, nsph
    do ig = 1, ngrid
      if ( ui(ig,isph).ne.zero ) then
        ic = ic + 1
        nderphi = dot_product( gradphi(:,ic), grid(:,ig) )
        c0(:, isph) = c0(:,isph) + w(ig)*ui(ig,isph)*nderphi*basis(:,ig)
      end if
    end do
  end do

  allocate (coefY(ncav,nbasis0,nsph), stat = istatus )
  ! memuse = memuse + ncav*nbasis0*nsph
  ! memmax = max(memmax,memuse)
  if ( istatus .ne. 0 ) then
      write(*,*)'wghpot_f : [1] allocation failed!'
      stop
  end if

  ! compute coef_bessel  
  ! write(6,*) 'lmax', lmax0, 'kappa', kappa
  ! write(6,*) tol_inf, tol_zero 
  do jsph = 1, nsph
    do l0 = 0, lmax0
      if (max(DI_ri(l0,jsph), SI_ri(l0,jsph)).gt.tol_inf) then
        termi = kappa
      else if (min(DI_ri(l0,jsph), SI_ri(l0,jsph)).lt.tol_zero) then
        termi = l0/rsph(jsph) + &
            & (l0 + one)*(kappa**2*rsph(jsph))/((two*l0 + one) * &
            & (two*l0 + three))
      else
        termi = DI_ri(l0,jsph)/SI_ri(l0,jsph)*kappa
      end if
      !write(*,*) SI_ri(l0,jsph), termi

      if (SK_ri(l0,jsph).gt.tol_inf) then
        termk = - (l0 + one)/rsph(jsph) - &
            & l0*(kappa**2*rsph(jsph))/((two*l0 - one)*(two*l0 + one))
      else if (SK_ri(l0,jsph).lt.tol_zero) then
        termk = -kappa
      else
        termk = DK_ri(l0,jsph)/SK_ri(l0,jsph)*kappa
      end if

      !write(*,*) SK_ri(l0,jsph), termk
      coef_bessel(l0,jsph) = one/(termi - termk)
      !write(*,*) DI_ri(l0,jsph), SI_ri(l0,jsph), coef_bessel(l0,jsph)
      !write(*,*) (min(-DK_ri(l0,jsph), SK_ri(l0,jsph)).lt.tol_zero), &
      !    & DK_ri(l0,jsph), termk
    end do
  end do

  kep = 0
  do isph = 1, nsph
    do ig = 1, ngrid
      if (ui(ig,isph).gt.zero) then
        kep = kep + 1
        sumSijn = zero
        do jsph = 1, nsph
          sumSijn_pre = sumSijn
          vij  = csph(:,isph) + rsph(isph)*grid(:,ig) - csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sijn = vij/rijn
          
          call SPHK_bessel(lmax0,rijn*kappa,NM,SK_rijn,DK_rijn)
          call ylmbas(sijn,basloc,vplm,vcos,vsin)

          do l0 = 0,lmax0
            if (SK_ri(l0,jsph).gt.tol_inf) then
              term = (rsph(jsph)/rijn)**(l0+1)
            else if (SK_ri(l0,jsph).lt.tol_zero) then
              term = (rsph(jsph)/rijn)*exp(-kappa*(rijn-rsph(jsph)))
            else
              term = SK_rijn(l0)/SK_ri(l0,jsph)
            end if
            coef_Ylm =  coef_bessel(l0,jsph)*term
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              sumSijn = sumSijn + c0(ind0,jsph)*coef_Ylm*basloc(ind0)
              coefY(kep,ind0,jsph) = coef_Ylm*basloc(ind0)
            end do
          end do
        end do

        !write(6,*) sumSijn, epsp, eps, ui(ig,isph)
        f(ig,isph) = -(epsp/eps)*ui(ig,isph) * sumSijn
      end if
    end do
  end do 

  deallocate( vplm, basloc, vcos, vsin, SK_rijn, DK_rijn  )
  return
  end subroutine wghpot_f


  subroutine matabx( n, x, y )
  implicit none 
  integer, intent(in) :: n
  real*8, dimension(nylm,nsph), intent(in) :: x
  real*8, dimension(nylm,nsph), intent(inout) :: y
  integer :: isph, istatus
  real*8, allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
  integer :: i
  ! allocate workspaces
  allocate( pot(ngrid), vplm(nylm), basloc(nylm), vcos(lmax+1), &
      & vsin(lmax+1) , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'Bx: allocation failed !'
    stop
  endif

  y = zero
  do isph = 1, nsph
    call calcv2_lpb( isph, pot, x, basloc, vplm, vcos, vsin )
    call intrhs( isph, pot, y(:,isph) )
    ! action of off-diagonal blocks
    y(:,isph) = - y(:,isph)
    ! add action of diagonal block
    y(:,isph) = y(:,isph) + x(:,isph)
  end do

  deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'matABx: allocation failed !'
    stop
  endif
  end subroutine matabx


  subroutine convert_ddcosmo(direction,vector)
  implicit none
  integer :: isph, l, m, ind
  integer, intent(in) :: direction
  real*8, intent(inout) :: vector(nylm,nsph)
  real*8 :: fac
  
  do isph = 1, nsph
    do l = 0, lmax
      ind = l*l + l + 1
      fac = four*pi/(two*dble(l) + one) 
      if (direction.eq.-1) fac = one/fac
      do m = -l, l
        vector(ind + m,isph) = fac*vector(ind + m,isph)
      end do
    end do
  end do
  return
  end subroutine convert_ddcosmo


  subroutine print_ddvector(label,vector)
  implicit none
  character(len=*) :: label
  real*8 :: vector(nylm,nsph)
  integer :: isph, lm
  
  write(6,*) label
  do isph = 1, nsph
    do lm = 1, nylm
      write(6,'(F15.8)') vector(lm,isph)
    end do
  end do
  return
  end subroutine print_ddvector

  
  subroutine lpb_hsp(rhs,Xe)
  implicit none
  real*8, dimension(nylm,nsph), intent(in) :: rhs
  real*8, dimension(nylm,nsph), intent(inout) :: xe
  integer :: isph, istatus, n_iter, info, c1, c2, cr
  real*8 :: tol, r_norm
  real*8, allocatable :: work(:,:)
  integer, parameter  :: gmm = 20, gmj = 25
  
  allocate(work(nsph*nylm,0:2*gmj + gmm + 2 - 1),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) ' LPB-HSP: failed allocation for GMRES'
    stop
  endif
   
  work = zero
  xe = rhs
  
  matAB = 1
  tol_gmres = 1.0d-8
  call gmresr(.false.,nsph*nylm,gmj,gmm,rhs,Xe,work,tol_gmres,'rel', &
      & n_iter_gmres,r_norm,matABx,info)

  deallocate(work)
  endsubroutine lpb_hsp


  subroutine calcv2_lpb ( isph, pot, x, basloc, vplm, vcos, vsin )
  integer, intent(in) :: isph
  real*8, dimension(nylm,nsph), intent(in) :: x
  real*8, dimension(ngrid), intent(inout) :: pot
  real*8, dimension(nylm), intent(inout) :: basloc
  real*8, dimension(nylm), intent(inout) :: vplm
  real*8, dimension(lmax+1), intent(inout) :: vcos
  real*8, dimension(lmax+1), intent(inout) :: vsin
  real*8, dimension(nylm) :: fac_cosmo, fac_hsp
  integer :: its, ij, jsph
  real*8 :: vij(3), sij(3)
  real*8 :: vvij, tij, xij, oij

  pot = zero
  do its = 1, ngrid
    if (ui(its,isph).lt.one) then
      do ij = inl(isph), inl(isph+1)-1
        jsph = nl(ij)

        ! compute geometrical variables
        vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij  = vvij/rsph(jsph) 

        if ( tij.lt.one ) then
          sij = vij/vvij
          call ylmbas(sij,basloc,vplm,vcos,vsin)
          call inthsp(vvij,rsph(jsph),jsph,basloc,fac_hsp)
          xij = fsw(tij,se,eta)
          if (fi(its,isph).gt.one) then
            oij = xij/fi(its,isph)
          else
            oij = xij
          end if
          pot(its) = pot(its) + oij*dot_product(fac_hsp,x(:,jsph))
        end if
      end do
    end if
  end do
  endsubroutine calcv2_lpb


  subroutine inthsp(rijn, ri, isph, basloc, fac_hsp)
  use bessel
  implicit none
  integer, intent(in) :: isph
  real*8, intent(in) :: rijn, ri
  real*8, dimension(nylm), intent(in) :: basloc
  real*8, dimension(nylm), intent(inout) :: fac_hsp
  real*8, dimension(0:lmax) :: SI_rijn, DI_rijn
  integer :: l, m, ind, NM

  SI_rijn = 0
  DI_rijn = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call SPHI_bessel(lmax,rijn*kappa,NM,SI_rijn,DI_rijn)
  
  do l = 0, lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      if ((SI_rijn(l).lt.zero) .or. (SI_ri(l,isph).lt.tol_zero) &
          & .or. (SI_rijn(l)/SI_ri(l,isph).gt.(rijn/ri)**l)) then
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if ( SI_ri(l,isph).gt.tol_inf) then
        fac_hsp(ind) = zero
      else
        fac_hsp(ind) = SI_rijn(l)/SI_ri(l,isph)*basloc(ind)
      end if
    end do
  end do
  endsubroutine inthsp


  subroutine intcosmo(tij, basloc, fac_cosmo)
  implicit none
  real*8,  intent(in) :: tij
  real*8, dimension(nylm), intent(in) :: basloc
  real*8, dimension(nylm), intent(inout) :: fac_cosmo
  integer :: l, m, ind
  do l = 0, lmax
    do  m = -l, l
        ind = l*l + l + 1 + m
        fac_cosmo(ind) = tij**l*basloc(ind)
    end do
  end do
  end subroutine intcosmo


  subroutine update_rhs(rhs_cosmo_init,rhs_hsp_init,rhs_cosmo, & 
      & rhs_hsp,Xr,Xe)
  use ddcosmo
  use bessel
  implicit none
  real*8, dimension(nylm,nsph), intent(in) :: rhs_cosmo_init, &
      & rhs_hsp_init
  real*8, dimension(nylm,nsph), intent(inout) :: rhs_cosmo, rhs_hsp
  real*8, dimension(nylm,nsph) :: rhs_plus
  real*8, dimension(nylm,nsph), intent(in) :: Xr, Xe
  integer :: isph, jsph, ig, kep, ind, l1,m1, ind1, ind0, count, istatus
  real*8, dimension(3) :: vij
  real*8, dimension(nylm,nsph) :: diff_re
  real*8, dimension(nbasis0,nsph) :: diff0
  real*8, dimension(nylm,nylm,nsph) :: smat
  real*8, dimension(ncav) :: diff_ep
  real*8 :: Qval, rijn, val
  integer :: c0, cr, c_qmat, c_init, c_ep0, c_ep1 !, nbasis_appro
      
  if (firstoutiter) then
    allocate(coefvec(ngrid,nylm,nsph),Pchi(nylm,nbasis0,nsph), &
        & stat = istatus)
    if (istatus.ne.0) then
      write(*,*)'update_rhs : [1] allocation failed!'
      stop
    end if
  end if

  ! compute P_chi matrix
  ! TODO: probably has to be declared somewhere
  ! and i have to recover mkpmat 
  if (firstoutiter) then      
    do jsph = 1, nsph
      call mkpmat(jsph, Pchi(:,:,jsph))
    end do    
  end if 
      
  ! precompute coefvec of Qmat, Cost: linear scaling
  ! TODO: remove all the precomputations
  ! or, if they are really necessary, do them in a separate subroutine
  if (firstoutiter) then 
    do isph = 1,nsph
      do ig = 1,ngrid
        if (ui(ig, isph) .gt. 0) then
          do ind  = 1, nylm 
            coefvec(ig,ind,isph) = w(ig)*ui(ig,isph)*basis(ind,ig)
          end do
        end if
      end do
    end do
  end if
      
  ! precompute termimat
  ! TODO: same as before
  if (firstoutiter) then
    do jsph = 1, nsph
      do l1 = 0, lmax
        if (max(DI_ri(l1,jsph),SI_ri(l1,jsph)).gt.tol_inf) then
          termimat(l1,jsph) = kappa
        else if (min(DI_ri(l1,jsph),SI_ri(l1,jsph)).lt.tol_zero) then
          termimat(l1,jsph) = l1/rsph(jsph) + &
              & (l1 + 1)*(kappa**2*rsph(jsph))/((two*l1 + one) * &
              & (two*l1 + three))
        else
          termimat(l1,jsph) = DI_ri(l1,jsph)/SI_ri(l1,jsph)*kappa
        end if
      end do
    end do
  end if

  if (firstoutiter) then
    if (iprint .gt. 0) then  
      write(*,999) dble(c_init-c0)/dble(cr)
 999  format('Time of initializing Pchi, coefvec, termi: ',f8.3, &
          & ' secs.')
    end if
  end if

  ! diff_re = epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
  diff_re = zero 
  do jsph = 1, nsph
    do l1 = 0, lmax
      do m1 = -l1,l1
        ind1 = l1**2 + l1 + m1 + 1
        diff_re(ind1,jsph) = ((epsp/eps)*(l1/rsph(jsph)) * &
            & xr(ind1,jsph) - termimat(l1,jsph)*Xe(ind1,jsph))
      end do
    end do
  end do

  ! diff0 = Pchi * diff_er, linear scaling
  ! TODO: probably doing PX on the fly is better 
  diff0 = zero 
  do jsph = 1, nsph
    do ind0 = 1, nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & Pchi(:,ind0, jsph))
    end do
  end do

  ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
  diff_ep = zero
  do kep = 1, ncav
    val = zero
    do jsph = 1, nsph 
      do ind0 = 1, nbasis0
        val = val + diff0(ind0,jsph)*coefY(kep,ind0,jsph)
      end do
    end do
    diff_ep(kep) = val 
  end do

  rhs_plus = zero
  kep = 0
  do isph = 1, nsph
    do ig = 1, ngrid
      if (ui(ig,isph).gt.zero) then 
        kep = kep + 1
        do ind = 1, nylm
          rhs_plus(ind,isph) = rhs_plus(ind,isph) + &
              & coefvec(ig,ind,isph)*diff_ep(kep)
        end do
      end if
    end do
  end do

  rhs_cosmo = rhs_cosmo_init - rhs_plus
  rhs_hsp = rhs_hsp_init - rhs_plus 

  return
  end subroutine update_rhs  

  subroutine mkpmat( isph, pmat )
  implicit none
  integer,  intent(in) :: isph
  real*8, dimension(nylm, (lmax0+1)**2), intent(inout) :: pmat
  integer :: l, m, ind, l0, m0, ind0, its, nbasis0
  real*8  :: f, f0

  pmat(:,:) = zero
  do its = 1,ngrid
    if (ui(its,isph).ne.0) then
      do l = 0, lmax
        ind = l*l + l + 1
        do m = -l,l
          f = w(its) * basis(ind+m,its) * ui(its,isph)
          do l0 = 0, lmax0
            ind0 = l0*l0 + l0 + 1
            do m0 = -l0, l0
              f0 = basis(ind0+m0,its)
              pmat(ind+m,ind0+m0) = pmat(ind+m,ind0+m0) + f * f0
            end do
          end do
        end do
      end do
    end if
  end do
  endsubroutine mkpmat

end module
