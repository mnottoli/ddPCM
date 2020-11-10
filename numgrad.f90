program numgrad
use ddcosmo
use ddpcm_lib, only: ddpcm, ddpcm_init, ddpcm_finalize, ddpcm_forces, ddpcm_zeta
! 
!      888      888  .d8888b.   .d88888b.   .d8888b.  888b     d888  .d88888b.  
!      888      888 d88P  Y88b d88P" "Y88b d88P  Y88b 8888b   d8888 d88P" "Y88b 
!      888      888 888    888 888     888 Y88b.      88888b.d88888 888     888 
!  .d88888  .d88888 888        888     888  "Y888b.   888Y88888P888 888     888 
! d88" 888 d88" 888 888        888     888     "Y88b. 888 Y888P 888 888     888 
! 888  888 888  888 888    888 888     888       "888 888  Y8P  888 888     888 
! Y88b 888 Y88b 888 Y88b  d88P Y88b. .d88P Y88b  d88P 888   "   888 Y88b. .d88P 
!  "Y88888  "Y88888  "Y8888P"   "Y88888P"   "Y8888P"  888       888  "Y88888P"  
!                                                                              
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Paolo Gatto        !
!  Eric Cancès, Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and          !
!  Benedetta Mennucci.                                                         !
!                             ALL RIGHT RESERVED.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! A modular implementation of COSMO using a domain decomposition linear scaling
! strategy.
!
! This code is governed by the LGPL license and abiding by the rules of 
! distribution of free software.
! This program is distributed in the hope that it will be useful, but  
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Lesser General Public License for more details.
!
! Users of this code are asked to include the following references in their
! publications:
!
! [1] E. Cancès, Y. Maday, B. Stamm
!     "Domain decomposition for implicit solvation models"
!     J. Chem. Phys. 139, 054111 (2013)
!
! [2] F. Lipparini, B. Stamm, E. Cancès, Y. Maday, B. Mennucci
!     "Fast Domain Decomposition Algorithm for Continuum Solvation Models: 
!      Energy and First Derivatives"
!     J. Chem. Theory Comput. 9, 3637–3648 (2013)
!
! Also, include one of the three following reference depending on whether you
! use this code in conjunction with a QM [3], Semiempirical [4] or Classical [5]
! description of the solute:
!
! [3] F. Lipparini, G. Scalmani, L. Lagardère, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum, classical, and hybrid QM/MM calculations in solution: General 
!      implementation of the ddCOSMO linear scaling strategy"
!     J. Chem. Phys. 141, 184108 (2014)
!     (for quantum mechanical models)
!
! [4] F. Lipparini, L. Lagardère, G. Scalmani, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum Calculations in Solution for Large to Very Large Molecules: 
!      A New Linear Scaling QM/Continuum Approach"
!     J. Phys. Chem. Lett. 5, 953-958 (2014)
!     (for semiempirical models)
!
! [5] F. Lipparini, L. Lagardère, C. Raynaud, B. Stamm, E. Cancès, B. Mennucci
!     M. Schnieders, P. Ren, Y. Maday, J.-P. Piquemal
!     "Polarizable Molecular Dynamics in a Polarizable Continuum Solvent"
!     J. Chem. Theory Comput. 11, 623-634 (2015)
!     (for classical models, including polarizable force fields
!     
! The users of this code should also include the appropriate reference to the
! COSMO model. This distribution includes the routines to generate lebedev
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the routines
! are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015
!            Paolo Gatto,       December 2017
!            Filippo Lipparini, March 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
implicit none
!
integer :: i, ii, isph, ig, n
real*8  :: tobohr, esolv, xx(1)
real*8, parameter :: toang=0.52917721092d0, tokcal=627.509469d0, step=0.001
real*8, allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
real*8, allocatable :: phi(:), psi(:,:)
real*8, allocatable :: sigma(:,:), s(:,:)
real*8, allocatable :: numfx(:,:), fx(:,:), zeta(:), ef(:,:)
logical :: do_adjoint
real*8 :: ep, em

! read parameters from input file
open (unit=100,file='Input.txt',form='formatted',access='sequential')
read(100,*) iprint      ! printing flag
read(100,*) nproc       ! number of openmp threads
read(100,*) lmax        ! max angular momentum of spherical harmonics basis
read(100,*) ngrid       ! number of lebedev points
read(100,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
read(100,*) igrad       ! whether to compute (1) or not (0) forces
read(100,*) eps         ! dielectric constant of the solvent
read(100,*) eta         ! regularization parameter
read(100,*) n           ! number of atoms

! read geometry
allocate (x(n),y(n),z(n),rvdw(n),charge(n))
do i = 1, n
  read(100,*) charge(i), x(i), y(i), z(i), rvdw(i)
end do
tobohr = 1.0d0/toang
x    = x*tobohr
y    = y*tobohr
z    = z*tobohr
rvdw = rvdw*tobohr
close (100)

allocate(fx(3,n),numfx(3,n))

! analytical forces
call ddpcm_solver(x,y,z,charge,rvdw,n,esolv,.true.,fx)

! numerical forces
do i = 1, n
  x(i) = x(i) + step
  call ddpcm_solver(x,y,z,charge,rvdw,n,ep,.false.,xx)
  x(i) = x(i) - 2.0d0*step
  call ddpcm_solver(x,y,z,charge,rvdw,n,em,.false.,xx)
  x(i) = x(i) + step
  numfx(1,i) = (ep - em)/(2.0d0*step)

  y(i) = y(i) + step
  call ddpcm_solver(x,y,z,charge,rvdw,n,ep,.false.,xx)
  y(i) = y(i) - 2.0d0*step
  call ddpcm_solver(x,y,z,charge,rvdw,n,em,.false.,xx)
  y(i) = y(i) + step
  numfx(2,i) = (ep - em)/(2.0d0*step)

  z(i) = z(i) + step
  call ddpcm_solver(x,y,z,charge,rvdw,n,ep,.false.,xx)
  z(i) = z(i) - 2.0d0*step
  call ddpcm_solver(x,y,z,charge,rvdw,n,em,.false.,xx)
  z(i) = z(i) + step
  numfx(3,i) = (ep - em)/(2.0d0*step)
end do

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, n
  write(6,'(6F20.10)') fx(1,i), fx(2,i), fx(3,i), numfx(1,i), numfx(2,i), numfx(3,i)
end do

deallocate(fx,numfx,x,y,z,rvdw,charge)
contains 

  subroutine ddpcm_solver(x,y,z,charge,rvdw,n,esolv,do_forces,forces)
  ! just a wrapper routine
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(n), y(n), z(n), charge(n), rvdw(n)
  logical, intent(in) :: do_forces
  real*8, intent(out) :: forces(3,n)
  real*8, intent(out) :: esolv
  real*8, allocatable :: phi(:), psi(:,:), ef(:,:), zeta(:)
  real*8, allocatable :: fele(:,:)

  call ddinit(n,x,y,z,rvdw)
  allocate (phi(ncav),psi(nbasis,n))
  call mkrhs(n,charge,x,y,z,ncav,ccav,phi,nbasis,psi)

  ! solve ddpcm
  call ddpcm_init()
  call ddpcm(do_forces,phi,psi,esolv)

  ! compute analytical forces
  if (do_forces) then
    allocate(fele(3,n))

    ! first the geometrical contribution to the forces
    call ddpcm_forces(phi,fx)
    
    ! form the "zeta" intermediate 
    allocate (zeta(ncav))
    call ddpcm_zeta(zeta)
  
    ! compute solute electric field
    allocate(ef(3,max(n,ncav)))
    call efld(n,charge,csph,ncav,ccav,ef)
  
    fele = zero
    ! compute rhs contribution to the forces
    ii = 0
    do isph = 1, nsph
      do ig = 1, ngrid
        if (ui(ig,isph).gt.zero) then 
          ii = ii + 1
          fele(:,isph) = fele(:,isph) - zeta(ii)*ef(:,ii)
        end if
      end do
    end do
    call efld(ncav,zeta,ccav,n,csph,ef)
    do isph = 1, nsph
      fele(:,isph) = fele(:,isph) - ef(:,isph)*charge(isph)
    end do

    write(6,*) 'Electrostatic forces'
    do isph = 1, nsph
      write(6,'(1x,i5,3f16.8)') isph, fele(:,isph)
    end do

    fx = fx + fele
    
    deallocate (zeta,ef,fele)
  end if
  deallocate(phi,psi)
  call memfree
  call ddpcm_finalize()
  end subroutine ddpcm_solver

end program
