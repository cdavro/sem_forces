program matrix
 real*8,allocatable :: mat(:,:)
 real*8,dimension(3,3) :: smat
 character(len=8),allocatable :: atomlist (:)
 real*8,allocatable :: atomcoord(:,:)

 real*4 :: xvector
 real*4 :: yvector
 real*4 :: zvector

 integer :: atm1
 integer :: atm2
 integer :: smx
 integer :: smy
 integer :: nbatoms
 integer :: matsize

 integer ( kind = 4 ), parameter :: n = 3
 real ( kind = 8 ) d(n)
 real ( kind = 8 ) error_frobenius
 integer ( kind = 4 ) it_max
 integer ( kind = 4 ) it_num
 integer ( kind = 4 ) rot_num
 real ( kind = 8 ) v(n,n)
 
! print *, "Number of atmoms ?"
! read(*,*) nbatoms

 nbatoms=7
 matsize=nbatoms*3

 open(unit=50,file="matrix")
  allocate(atomlist(matsize))
  read(50,*) atomlist
  print *, "Atomlist : "
  print *
  print *, atomlist
!  print *, "Atom number 1 ?"
!  read(*,*) atm1
!  print *, "Atom number 2 ?"
!  read(*,*) atm2
  atm1 = 1
  atm2 = 2
  smx=((atm1*3)-3)
  smy=((atm2*3)-3)
  allocate(mat(matsize,matsize))
  read(50,*)mat
 close(50)
 do i=1,3
  do j=1,3
   smat(i,j) = mat((smx+i),(smy+j))*627.509
  end do
 end do
 smat=-1*transpose(smat)

 call r8mat_print ( n, n, smat, '  Forces matrix between Atom 1 and Atom 2:' )
  it_max = 1000
 call jacobi_eigenvalue ( n, smat, it_max, v, d, it_num, rot_num )
 write(*,"(a)") &
  " ", &
  " ", &
  "Eigen Values:", &
  " "
 do i=1,3
  write(*, "(E14.5E2,a)", advance="no") &
   d(i), &
   " "
 end do
 write (*,"(a)", advance="yes") &
  " ", &
  " ", &
  "Eigen Vectors:", &
  " "
 do i=1,3
  do j=1,3
   write(*, "(E14.6E2,a)" , advance="no") &
     v(i,j), &
     "  "
  end do
   print*
   end do
call r8mat_is_eigen_right ( n, n, smat, v, d, error_frobenius )
  write ( *, '(a)' ) ''
  write ( *, '(a,E14.6E2)' ) &
   'Frobenius norm error in eigensystem A*V-D*V = ', &  
   error_frobenius 

 open(unit=51,file="coord")
  allocate(atomcoord(5,nbatoms))
  read(51,*) atomcoord
 close(51)
 atomcoord = transpose(atomcoord)
 xvector = atomcoord(atm2,3) - atomcoord(atm1,3)
 yvector = atomcoord(atm2,4) - atomcoord(atm1,4)
 zvector = atomcoord(atm2,5) - atomcoord(atm1,5)
 write (*, "(a)") &
 " ", &
 "Direction vector between Atom 1 and Atom 2:", &
 " "
 write (*, "(E14.6E2)" ,advance="NO") xvector
 write (*, "(E14.6E2)" ,advance="NO") yvector
 write (*, "(E14.6E2)") zvector
end program matrix

