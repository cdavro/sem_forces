program matrix
implicit none
 real*8,allocatable :: mat(:,:), atomcoord (:,:)
 real*8,dimension(3,3) :: smat, egvec
 real*8,dimension(3) ::  egval

 character(len=4),allocatable :: crdlist (:)
 character(:),allocatable :: atmlc
 character(len=4),allocatable :: atl (:)
 
 integer :: i, j, n = 1
 real*4 :: xvector, yvector, zvector
 integer*4 :: atm1, atm2, smx, smy, nbatoms, matsize

 real*8 :: k1, k2, k3, kf
 
 integer*4 :: o = 3, it_max=1000, it_num, rot_num
 real*8 :: error_frobenius
 
 write (*,"(a)") "Please input the number of atoms ?"
 read(*,*) nbatoms

 matsize=nbatoms*3

 open(unit=50,file="matrix")
  allocate(crdlist(matsize))
  read(50,*) crdlist
  allocate(atl(nbatoms))
  do i=1,nbatoms
   atmlc = crdlist(n)
   atmlc = trim(atmlc(2:4))
   atl (i) = atmlc
   n=n+3
  end do
  write (*,"(a)") "", "This is the list of atoms:"
  do i=1,nbatoms
   write (*,"(a)", advance="no") atl(i)
  end do
  write (*,"(a)") "", "", "Select atom number one:"
  read(*,*) atm1
  write (*,"(a)") "", "Select atom number two:"
  read(*,*) atm2
  smx=((atm1*3)-3)
  smy=((atm2*3)-3)
  allocate(mat(matsize,matsize))
  read(50,*)mat
 close(50)
 do i=1,3
  do j=1,3
   smat(i,j) = mat((smx+i),(smy+j))
  end do
 end do
 
 smat=-1*transpose(smat)
 
 write (*,"(a)") ""
 write (*,"(5a)") "-1 * Cartesian Force Matrix (kcal/mol) between atom ", trim(atl(atm1)), " and atom ", trim(atl(atm2)), " :"
 write (*,"(a)") ""
 do i=1,3
  do j=1,3
    write(*, "(E14.6E2,a)" , advance="no") smat(i,j), ""
  end do
  write (*,"(a)") ""
 end do

 call jacobi_eigenvalue ( 3, smat, it_max, egvec, egval, it_num, rot_num )
    write ( *, '(a,i4)' ) '  Number of iterations = ', it_num
  write ( *, '(a,i4)' ) '  Number of rotations  = ', rot_num

 write(*,"(a)") "", "Eigen Values (kcal/mol):", ""
 do i=1,3
  write(*, "(E14.5E2,a)", advance="no") egval(i), ""
 end do
 write (*,"(a)") "", "", "Eigen Vectors:", ""
 do i=1,3
  do j=1,3
   write(*, "(E14.6E2,a)" , advance="no") &
     egvec(i,j), &
     "  "
  end do
  write (*,"(a)") ""  
 end do

call r8mat_is_eigen_right ( o, o, smat, egvec, egval, error_frobenius )
  write ( *, '(a)' ) ""
  write ( *, '(a,E14.6E2)' ) 'Frobenius norm error in eigensystem A*V-D*V =', error_frobenius 

 open(unit=51,file="coord")
  allocate(atomcoord(5,nbatoms))
  read(51,*) atomcoord
 close(51)
 atomcoord = transpose(atomcoord)
 xvector = atomcoord(atm2,3) - atomcoord(atm1,3)
 yvector = atomcoord(atm2,4) - atomcoord(atm1,4)
 zvector = atomcoord(atm2,5) - atomcoord(atm1,5)
 write (*,"(a)") ""
 write (*,"(5a)") "Direction vector between atom ", trim(atl(atm1)), " and atom ", trim(atl(atm2)), " :"
 write (*,"(a)") ""
 write (*,"(E14.6E2)" ,advance="NO") xvector
 write (*,"(E14.6E2)" ,advance="NO") yvector
 write (*,"(E14.6E2)") zvector
 k1 = egval(1) * abs( (egvec(1,1) * xvector) + (egvec(1,2) * yvector) + (egvec(1,3) * zvector)) 
 k2 = egval(2) * abs( (egvec(2,1) * xvector) + (egvec(2,2) * yvector) + (egvec(2,3) * zvector)) 
 k3 = egval(3) * abs( (egvec(3,1) * xvector) + (egvec(3,2) * yvector) + (egvec(3,3) * zvector)) 

 kf=(k1+k2+k3)
 write (*,"(a)") ""
 write (*,"(5a,F10.4,a)") "Approximate force constant between atom ", trim(atl(atm1)), " and atom ", trim(atl(atm2)), &
 " : ", kf*2240.874995, " kcal/mol"
 deallocate (mat)
 deallocate (atl)
end program matrix

