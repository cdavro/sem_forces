program matrix
	real,allocatable :: mat(:,:)
	real*8,dimension(3,3) :: smat
	character(len=8),allocatable :: title (:)
	
	integer :: smx
	integer :: smy
	integer :: atoms
        
	integer ( kind = 4 ), parameter :: n = 3
	real ( kind = 8 ) d(n)
	real ( kind = 8 ) error_frobenius
	integer ( kind = 4 ) it_max
	integer ( kind = 4 ) it_num
	integer ( kind = 4 ) rot_num
	real ( kind = 8 ) v(n,n)
	
!	print *, "Number of atmoms ?"
!	read(*,*) atoms
!	print *, "Atoms number 1 ?"
!	read(*,*) smx
!	print *, "Atoms number 2 ?"
!	read(*,*) smy

	atoms=7*3
	smx=((1*3)-3)
	smy=((2*3)-3)

	open(unit=50,file="matrix")
	allocate(title(atoms))
	read(50,*) title
	allocate(mat(atoms,atoms))
	read(50,*)mat
	close(50)
	do i=1,3
	 do j=1,3
	  smat(i,j) = mat((smx+i),(smy+j))*627.509
	 end do
        end do
	smat=-1*transpose(smat)

 call r8mat_print ( n, n, smat, '  Input matrix A:' )
  it_max = 1000
 call jacobi_eigenvalue ( n, smat, it_max, v, d, it_num, rot_num )
!  write ( *, '(a)' ) ''
!  write ( *, '(a,i3)' ) '  Number of iterations = ', it_num
!  write ( *, '(a,i3)' ) '  Number of rotations  = ', rot_num
!
  call r8vec_print ( n, d, '  Eigenvalues D:' )
!
  call r8mat_print ( n, n, v, '  Eigenvector matrix V:' )
!!
!!  Compute eigentest.
!!
  call r8mat_is_eigen_right ( n, n, smat, v, d, error_frobenius )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  Frobenius norm error in eigensystem A*V-D*V = ', &
    error_frobenius
  return
end program matrix

