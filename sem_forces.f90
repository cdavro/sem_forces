program matrix
implicit none
 
 real,allocatable :: atmat(:,:), atmcrd (:,:)
 
 integer,allocatable :: atmlnb (:), atmlan (:)
 
 real,dimension(3,3) :: satmat, egvec
 real,dimension(3) ::  egval

 character(len=2),allocatable :: atmlna (:)
 character(len=4), allocatable :: atmlnbc (:)

  
 integer :: o = 3, it_max=1000, it_num, rot_num
 real :: error_frobenius

 real :: xvector, yvector, zvector
 integer :: atm1, atm2, smx, smy, nbatm, matsize
 real :: k1, k2, k3, kf

 character :: revar
 integer :: i, j, n = 1
  
 write (*,"(a)") "Please input the number of atoms ?"
 read(*,*) nbatm

 matsize=nbatm*3
 
 open(unit=50,file="toc",status='old', access='sequential', form='formatted', action='read')
  allocate(atmcrd(5,nbatm))
  allocate(atmlnb(matsize))
  allocate(atmlan(matsize))
  allocate(atmlna(matsize))
  allocate(atmat(matsize,matsize))
  read(50,*) revar ! advance one line
   read(50,*) atmcrd ! fill the coordinates matrix (ATMLNB,ATMLAN,X,Y,Z)
  read(50,*) revar !advance one line
    read(50,*) atmlnb ! fill the atom number list (ATM1,ATM1,ATM1,ATM2,ATM2,ATM2,..)
    read(50,*) atmlan ! fill the atomic number list (same)
    read(50,*) atmlna ! fill the atom name list (same)
  read(50,*) revar ! advance one line
    read(50,*) atmat
 close(50)
 print*, atmlnb
 allocate(atmlnbc(matsize))
 do i=1,matsize
 print*, atmlnb(i) 
 write(atmlnbc(i),*) atmlnb(i) 
! atmlnbc(i)=atmlnb(i)

 enddo
! write(atmlnbc,*) atmlnb
 print*, "hola"
 write (*,"(a)") "", "This is the list of atoms:"
 do i=1,nbatm
   write (*,"(2a)", advance="no") trim(atmlna(i*3-2)), trim(atmlnbc((i*3-2)))
 end do
 write (*,"(a)") "", "", "Select atom number one:"
 read(*,*) atm1
 write (*,"(a)") "", "Select atom number two:"
 read(*,*) atm2
 smx=((atm1*3)-3)
 smy=((atm2*3)-3)
 do i=1,3
  do j=1,3
   satmat(i,j) = atmat((smx+i),(smy+j))
  end do
 end do

 satmat=-1*transpose(satmat)
 print*, atmlnb 
 write (*,"(a)") ""
 write (*,"(7a)") "-1 * Cartesian Force Matrix (kcal/mol) between atom ", atmlna(atm1*3-2), atmlnb(atm1*3-2), &
  " and atom ", atmlna(atm2*3-2), atmlnb(atm2*3-2), " :"
 write (*,"(a)") ""
 do i=1,3
  do j=1,3
    write(*, "(E14.6E2,a)" , advance="no") satmat(i,j), ""
  end do
  write (*,"(a)") ""
 end do

 call jacobi_eigenvalue ( 3, satmat, it_max, egvec, egval, it_num, rot_num )
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

call r8mat_is_eigen_right ( o, o, satmat, egvec, egval, error_frobenius )
  write ( *, '(a)' ) ""
  write ( *, '(a,E14.6E2)' ) 'Frobenius norm error in eigensystem A*V-D*V =', error_frobenius 

 atmcrd = transpose(atmcrd)
 xvector = atmcrd(atm2,3) - atmcrd(atm1,3)
 yvector = atmcrd(atm2,4) - atmcrd(atm1,4)
 zvector = atmcrd(atm2,5) - atmcrd(atm1,5)
 write (*,"(a)") ""
 write (*,"(a)") ""
 write (*,"(E14.6E2)" ,advance="NO") xvector
 write (*,"(E14.6E2)" ,advance="NO") yvector
 write (*,"(E14.6E2)") zvector
 k1 = egval(1) * abs( (egvec(1,1) * xvector) + (egvec(1,2) * yvector) + (egvec(1,3) * zvector)) 
 k2 = egval(2) * abs( (egvec(2,1) * xvector) + (egvec(2,2) * yvector) + (egvec(2,3) * zvector)) 
 k3 = egval(3) * abs( (egvec(3,1) * xvector) + (egvec(3,2) * yvector) + (egvec(3,3) * zvector)) 

 kf=(k1+k2+k3)
 write (*,"(a)") ""
! write (*,"(5a,F10.4,a)") "Approximate force constant between atom ", trim(atl(atm1)), " and atom ", trim(atl(atm2)), &
! " : ", kf*2240.874995, " kcal/mol"
end program matrix

