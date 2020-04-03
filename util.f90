  subroutine prtmat(iout,a,m,n)
  implicit none
  integer, intent(in) :: iout, m, n
  real*8, intent(in) :: a(m,n)
  integer i, j
  do i = 1, m
    do j = 1, n
      write(iout,'(F12.6$)') a(i,j)
    end do
    write(iout,*)
  end do
  end subroutine prtmat

  subroutine cavpdb(string,crd,q,natoms,cav,ncav)
  implicit none
  character*(*), intent(in) :: string
  integer, intent(in) :: natoms, ncav
  real*8, intent(inout) :: crd(3,natoms), q(natoms), cav(3,ncav)
  real*8, parameter :: toang=0.52917721092d0
  integer :: i, j, k

  ! convert to angstrom
  crd = crd*toang
  cav = cav*toang

  ! pdb format
  10 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)

  open(unit=100,file=string)
  k = 1
  do i = 1, natoms
    write(100,10) 'ATOM  ',k,'C',' ',' AT','A',1,' ',crd(1,i),crd(2,i), &
      & crd(3,i),q(i),1.0d0,' ',' '
    k = k + 1
  end do
  do i = 1, ncav
    write(100,10) 'ATOM  ',k,'He',' ','CAV','A',1,' ',cav(1,i),cav(2,i), &
      & cav(3,i),1.0d0,1.0d0,' ',' '
    k = k + 1
  end do
  close(100)

  ! back to bohr
  crd = crd/toang
  cav = cav/toang
  end subroutine cavpdb
