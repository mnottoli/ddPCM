module ddlpb_lib
use ddcosmo
implicit none
logical :: firstiniter, firstoutiter, firstcosmo, firsthsp
logical :: do_diag
integer :: matAB

contains


subroutine ddlpb(phi,psi,gradphi,sigma,esolv)
  ! main ddlpb
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

  ! Build the right hand side
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
  firstcosmo = .true.
  firsthsp = .true.
  firstiniter = .true.
  firstoutiter = .true.

  do while (.not.converged)
    ! solve the ddcosmo step
    ! AXr = RHSr
    !call print_ddvector('rhs',rhs_1)
    n_iter = 200
    call jacobi_diis(nsph*nylm,iprint,ndiis,4,tol,rhs_1,xr,n_iter,&
      & ok,lx,ldm1x,hnorm)
    call convert_ddcosmo(1,xr)
    !call print_ddvector('xr',xr)
  
    ! solve ddlpb step
    ! BXe = RHSe
    call lpb_hsp(rhs_2,xe)
    call print_ddvector('xe',xe)
  
    ! update the RHS
    call update_rhs(f0 + g0,f0,rhs_1,rhs_2,xr,xe)
    call print_ddvector('rhs_1',rhs_1)
    call print_ddvector('rhs_2',rhs_2)
     
    stop
    ! check for convergency
    converged = .true.
    write(6,*) 'Iteration', iteration
    iteration = iteration + 1
  end do
  esolv = 1.2345d0
  write(6,*) 'Done'
  return
end subroutine ddlpb

subroutine matabx( n, x, y )
!      
! compute BX
! TODO: put in dd_operators
      implicit none 
      integer,                                         intent(in)      :: n
      real*8,  dimension(nylm,nsph), intent(in)      :: x
      real*8,  dimension(nylm,nsph), intent(inout) :: y
!
      integer             :: isph, istatus
      real*8, allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
      integer :: i
!
!     allocate workspaces
      allocate( pot(ngrid), vplm(nylm), basloc(nylm), vcos(lmax+1), &
                vsin(lmax+1) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'Bx: allocation failed !'
        stop
      endif
!
!     initialize
      y = zero
!$omp parallel do default(shared) private(isph,pot,basloc,vplm,vcos,vsin) &
!$omp schedule(dynamic)
!
!     loop over spheres
      do isph = 1, nsph
        call calcv2_lpb( isph, pot, x, basloc, vplm, vcos, vsin )
        write(6,*) 'I am sphere', isph, do_diag
        do i = 1, ngrid
          write(6,*) pot(i) 
        end do
        call intrhs( isph, pot, y(:,isph) )
!       action of off-diagonal blocks
        y(:,isph) = - y(:,isph)
!       add action of diagonal block
        y(:,isph) = y(:,isph) + x(:,isph)
      enddo
!
!     deallocate workspaces
      deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'matABx: allocation failed !'
        stop
      endif
!
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
integer             :: isph, istatus, n_iter, info, c1, c2, cr
real*8              :: tol, r_norm
real*8, allocatable :: work(:,:)
integer, parameter  :: gmm = 20, gmj = 25
!external            :: matABx

allocate( work(nsph*nylm, 0:2*gmj+gmm+2 -1) , stat=istatus )
if ( istatus.ne.0 ) then
    write(*,*) ' LPB-HSP: failed allocation for GMRES'
    stop
endif
work = zero
 
! initial guess
xe = rhs

matAB = 1 
firsthsp = .true.
call gmresr((iprint.gt.1), nsph*nylm, gmj, gmm, rhs, Xe, work, tol_gmres, 'rel', n_iter_gmres, r_norm, matABx, info )
endsubroutine lpb_hsp

subroutine calcv2_lpb ( isph, pot, x, basloc, vplm, vcos, vsin )
!
      integer,                        intent(in)    :: isph
      real*8, dimension(nylm,nsph), intent(in)    :: x
      real*8, dimension(ngrid),       intent(inout) :: pot
      real*8, dimension(nylm),      intent(inout) :: basloc
      real*8, dimension(nylm),      intent(inout) :: vplm
      real*8, dimension(lmax+1),      intent(inout) :: vcos
      real*8, dimension(lmax+1),      intent(inout) :: vsin
!
      real*8, dimension(nylm) :: fac_cosmo, fac_hsp
      
      integer :: its, ij, jsph
      real*8  :: vij(3), sij(3)
      real*8  :: vvij, tij, xij, oij
!
!------------------------------------------------------------------------
!     initialize
      pot(:) = zero
!     loop over grid points
      do its = 1, ngrid
!       contribution from integration point present
        if ( ui(its,isph).lt.one ) then
!         loop over neighbors of i-sphere
          do ij = inl(isph), inl(isph+1)-1
!           neighbor is j-sphere
            jsph = nl(ij)
!           compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
            vvij = sqrt( dot_product( vij, vij ) )
            tij  = vvij / rsph(jsph) 
!           compute s_n^ij = ( r_i + \rho_i s_n - r_j ) / | ... |
            sij = vij / vvij
!           point is INSIDE j-sphere
!           ------------------------
            if ( tij.lt.one ) then
!             compute Y_l^m( s_n^ij )
              call ylmbas( sij, basloc, vplm, vcos, vsin )
              call inthsp( vvij, rsph(jsph), jsph, basloc, fac_hsp)

              xij = fsw( tij, se, eta )
              if ( fi(its,isph).gt.one ) then
                oij = xij / fi(its,isph)
              else
                oij = xij
              end if
              pot(its) = pot(its) + oij * dot_product(fac_hsp, x(:,jsph))
              !write(6,'(3I8,3F15.8)') its, isph, jsph, pot(its), oij, dot_product(fac_hsp,x(:,jsph))
            end if
          end do
        end if
      end do
endsubroutine calcv2_lpb

subroutine inthsp(rijn, ri, isph, basloc, fac_hsp)
      use bessel
      implicit none
      integer,  intent(in) :: isph
      real*8,  intent(in) :: rijn, ri
      real*8, dimension(nylm), intent(in) :: basloc
      real*8, dimension(nylm), intent(inout) :: fac_hsp
      real*8, dimension(0:lmax) :: SI_rijn, DI_rijn
      
      integer :: l, m, ind, NM
            
       SI_rijn = 0
       DI_rijn = 0
       fac_hsp = 0
      !   Computation of modified spherical Bessel function values      
      call SPHI_bessel (lmax, rijn*kappa, NM, SI_rijn, DI_rijn)
      
      do l = 0, lmax
        do  m = -l, l
            ind = l*l + l + 1 + m
            if ( (SI_rijn(l) .lt. 0) .or. (SI_ri(l,isph).lt. tol_zero) .or. (SI_rijn(l)/SI_ri(l,isph) .gt. (rijn/ri)**l) ) then
                fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
                !write(*,*) isph, ind, SI_ri(l,isph), fac_hsp(ind)
                !stop
            else if ( SI_ri(l,isph).gt. tol_inf ) then
                fac_hsp(ind) = 0
            else
                fac_hsp(ind)  = SI_rijn(l)/SI_ri(l,isph)*basloc(ind)
                
                    !write(*,*) 'Bessel warning [1]: unaccurate bessel function in inthsp.'
!                      write(*,*) isph, ind, fac_hsp(ind)
                    !write(*,*) ind, rijn*kappa, ri*kappa, SI_rijn(l), SI_ri(l,isph)
!                      write(*,*) SI_rijn
!                      stop
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
endsubroutine intcosmo

subroutine update_rhs( rhs_cosmo_init, rhs_hsp_init, rhs_cosmo, rhs_hsp, Xr, Xe)
      use ddcosmo
      use bessel
      implicit none
      real*8, dimension(nylm,nsph), intent(in) :: rhs_cosmo_init, rhs_hsp_init
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
      integer :: c0, cr, c_qmat, c_init, c_ep0, c_ep1!nbasis_appro
      
!       real*8 :: dist_far
!       integer :: nbasis_far
!       logical :: use_appro
!       
!       dist_far = 64
!       nbasis_far = 1
!       use_appro = .false.
      
      if (firstoutiter) then
            allocate( coefvec(ngrid, nbasis, nsph), Pchi(nbasis,nbasis0,nsph), stat = istatus)
            memuse = memuse + ngrid*nbasis*nsph + nbasis*nbasis0*nsph
            memmax = max(memmax,memuse)
            if ( istatus .ne. 0 ) then
                write(*,*)'update_rhs : [1] allocation failed!'
                stop
            end if
       end if

      call system_clock(count_rate=cr)
      call system_clock(count=c0)
      
!     compute P_chi matrix
      if (firstoutiter) then      
          do jsph = 1, nsph
            call mkpmat( jsph, Pchi(:,:,jsph) )
          end do    
      end if 
      
!     compute Qmat, if Chi_e is zero, set Qmat to 0
!      if (firstoutiter) then  
!         !$omp parallel do default(shared) private(isph,kep,jsph,ind,ind0,Qval, Qmat)
!         do isph = 1,nsph
!             ! not buried
!             if ( (iep(isph+1)-iep(isph)) .gt. 0 ) then 
!                 ! loop over exposed grid points
!                 do kep = iep(isph), iep(isph+1)-1 
!                     ig = ep(kep) ! kep represents (isph, ig)
!                     do jsph = 1,nsph
!                         if ( (iep(jsph+1)-iep(jsph)) .gt. 0 ) then
!                             vij = csph(:,isph)-csph(:,jsph)
!                             if ( ((dot_product(vij,vij)) .gt. dist_far) .and. (use_appro) ) then ! far way atom pair 
!                                 ! like fast multipole method 
!                                 nbasis_appro  = nbasis_far 
!                             else 
!                                 nbasis_appro = nbasis
!                             end if 
!                             do ind1 = 1, nbasis_appro ! ind0 = 1,  zero contribution!
!                                 Qmat(kep, ind1, jsph) = dot_product( Pchi(ind1,:, jsph), coefY(kep,:,jsph) )
!                             end do
!                         end if
!                     end do
!                 end do
!             end if
!         end do
!         !$omp end parallel do
!        
!        ! save Qmat        
!         open (unit=100,file="Qmat.txt",action="write",status="replace")
!         do kep = 1, ncav
!               write (100,'(*(f10.6))') Qmat(kep,:,:)
!         end do
!         close (100)
!
!      end if

!     precompute coefvec of Qmat, Cost: linear scaling
      if (firstoutiter) then 
        do isph = 1,nsph
            do ig = 1,ngrid
                if (ui(ig, isph) .gt. 0) then
                    do ind  = 1, nbasis
                      coefvec(ig, ind, isph) = w(ig)*ui(ig, isph)*basis(ind,ig)
                    end do
                end if
            end do
        end do
      end if
      
!     precompute termimat      
      if (firstoutiter) then
        do jsph = 1,nsph
            do l1 = 0,lmax
                ! termi
                if ( max(DI_ri(l1,jsph), SI_ri(l1,jsph)) .gt. tol_inf) then
                    termimat(l1,jsph) = kappa
                else if ( min(DI_ri(l1,jsph), SI_ri(l1,jsph)) .lt. tol_zero) then
                    termimat(l1,jsph) = l1/rsph(jsph)+ (l1+1)*(kappa**2*rsph(jsph))/( (2*l1+1)*(2*l1+3) )
                else
                    termimat(l1,jsph) = DI_ri(l1,jsph)/SI_ri(l1,jsph)*kappa
                end if
            end do
        end do
      end if
      
      if (firstoutiter) then
        call system_clock(count=c_init)
        if (iprint .gt. 0) then  
            write(*,999) dble(c_init-c0)/dble(cr)
 999     format('Time of initializing Pchi, coefvec, termi: ',f8.3,' secs.')   
        end if
      end if
      
      ! diff_re = eps1/eps2*l1/ri*Xr - i'(ri)/i(ri)*Xe,
      diff_re = 0
      do jsph = 1,nsph
        do l1 = 0,lmax
            do m1 = -l1,l1
                ind1 = l1**2+l1+m1+1
                diff_re(ind1,jsph) = ( eps1/eps2*l1/rsph(jsph)*Xr(ind1,jsph) - termimat(l1,jsph)*Xe(ind1,jsph) )
             end do
        end do
      end do
      
      ! diff0 = Pchi * diff_er, linear scaling
      diff0 = 0
      do jsph = 1,nsph
          do ind0 = 1,nbasis0
            diff0(ind0, jsph) = dot_product(diff_re(:,jsph), Pchi(:,ind0, jsph))
          end do
      end do
      
      ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
      call system_clock(count=c_ep0)
      diff_ep = 0
      do kep = 1,ncav
        val = 0
        do jsph = 1,nsph
            do ind0 = 1,nbasis0
            val = val + diff0(ind0,jsph)*coefY(kep,ind0,jsph)
            end do
        end do
        diff_ep(kep) = val !sum( diff0(:,:)*coefY(kep,:,:) ) ! of dimesion (nbasis0, nsph) , accelaration?
      end do
      call system_clock(count=c_ep1)
      
      write(*,1000) dble(c_ep1-c_ep0)/dble(cr)
 1000     format('Time of updating diff_ep: ',f8.3,' secs.')   
 
!           ! Save coefY
!         open (unit=100,file="coefY.txt",action="write",status="replace")
!         do kep = 1, ncav
!               write (100,'(*(f10.6))') coefY(kep,:,:)
!         end do
!         close (100)
!         stop
              
!     diff_ep = 0
!       do kep = 1,ncav
!           do jsph = 1,nsph
!             do ind1 = 1,nbasis
!                 diff_ep(kep) = diff_ep(kep) + diff_re(ind1,jsph)*Qmat(kep,ind1,jsph)
!             end do
!           end do
!       end do
      
!write(*,*) 'diff ', diff_re(1,1), Xe(1,1)
          
!     compute rhs_plus = - C1 Xr - C2 Xe,   COST: linear scaling M*nbasis*Nleb
      rhs_plus = 0
      do isph = 1,nsph
         if ( (iep(isph+1)-iep(isph)) .gt. 0 ) then ! not buried
            do ind = 1,nbasis
                do kep = iep(isph), iep(isph+1)-1
                    ig = ep(kep) 
                    rhs_plus(ind,isph) = rhs_plus(ind,isph) + coefvec(ig, ind, isph)*diff_ep(kep)
                end do
            end do
         end if
      end do
      
!     update rhs_cosmo, rhs_hsp
      rhs_cosmo = rhs_cosmo_init - rhs_plus
      rhs_hsp      = rhs_hsp_init - rhs_plus 
      
      
endsubroutine update_rhs  



end module
