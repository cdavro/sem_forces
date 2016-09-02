program matrix
 implicit none
 integer :: N=3 
 integer :: LDA=3,LDVL=3,LDVR=3
 integer,parameter :: LWMAX=1000
 integer :: INFO, LWORK 
 real*8,dimension(3,3) :: AM, VR, VL
 real*8,dimension(3) ::  WR, WI
 real*8,allocatable :: WORK (:)
 
 real,allocatable :: atmat(:,:), atmcrd (:,:)
 integer,allocatable :: atmlnb (:), atmlan (:)
 character(len=2),allocatable :: atmlna (:) 
 real,dimension(3) :: dvector
 real :: dist
 real :: xvector, yvector, zvector
 integer :: atm1, atm2, smx, smy, nbatm, matsize
 real :: k1, k2, k3, kf

 character :: revar
 integer :: i, j

 external DGEEV
 
 write (*,"(a)") "Please input the number of atoms ?"
 read(*,*) nbatm
 write (*,"(a)") ""
 
 matsize=nbatm*3
!READING FILE 
  open(unit=50,file="data",status='old', access='sequential', form='formatted', action='read')
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
 write (*,"(a)") "", "This is the list of atoms:"
 do i=1,nbatm
   write (*,"(a,I0,a)" ,advance="NO") trim(atmlna(i*3-2)), atmlnb((i*3-2)), " "
 end do
 write (*,"(a)") "", "", "Select atom number one:"
 read(*,*) atm1
 write (*,"(a)") "", "Select atom number two:"
 read(*,*) atm2
 smx=((atm1*3)-3)
 smy=((atm2*3)-3)
 do i=1,3
  do j=1,3
   AM(i,j) = atmat((smx+i),(smy+j))
  end do
 end do
!END READ FILE
!PRINT FORCE MATRIX
 AM=-1*transpose(AM)
 write (*,"(a)") ""
write (*,"(2a,I0,2a,I0,a)") "Interatomic Cartesian Force Matrix (a.u.) between atom ", trim(atmlna(atm1*3-2)), atmlnb(atm1*3-2), &
  " and atom ", trim(atmlna(atm2*3-2)), atmlnb(atm2*3-2), " :"
 write (*,"(a)") ""
 do i=1,3
  do j=1,3
    write(*, "(E14.6E2,a)" , advance="no") (AM(i,j)), ""
  end do
  write (*,"(a)") ""
 end do
 !END PRINT FORCE MATRIX
!CALL FOR INTEL MKL 
 LWORK = -1
 allocate(WORK(LWORK))
 call DGEEV ('V', 'V', N, AM, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
 LWORK = MIN (LWMAX, INT(WORK(1)))
 deallocate(WORK)
 allocate(WORK(LWORK)) 
 call DGEEV ('V', 'V', N, AM, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
 IF( INFO.GT.0 ) THEN
  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
   STOP
  END IF
 deallocate(WORK)
!END CALL FOR INTEL MKL
 write(*,"(a)") ""
 write(*,"(a)") "", "Eigen Values (real) (a.u.):", ""
 do i=1,3
  write(*, "(E14.5E2,a)", advance="no") (WR(i)), ""
 end do
 write(*,"(a)") ""
 write(*,"(a)") "", "Eigen Values (comp) (a.u.):", ""
 do i=1,3
  write(*, "(E14.5E2,a)", advance="no") (WI(i)), ""
 end do
 write (*,"(a)") "", "", "Eigen Vectors:", ""
 do i=1,3
  do j=1,3
   write(*, "(E14.6E2,a)" , advance="no") &
     VR(i,j), &
     "  "
  end do
  write (*,"(a)") ""  
 end do
 xvector = (atmcrd(3,atm2) - atmcrd(3,atm1))
 yvector = (atmcrd(4,atm2) - atmcrd(4,atm1))
 zvector = (atmcrd(5,atm2) - atmcrd(5,atm1))
 dist = SQRT(abs(xvector**2)+abs(yvector**2)+abs(zvector**2)) 
 dvector(1) = xvector
 dvector(2) = yvector
 dvector(3) = zvector
 write (*,"(a)") ""
 write (*,"(a)") ""
 write(*,"(a,F6.2,a)", advance='NO') "Distance (Å): ", dist*0.529177249, ""
 write (*,"(a)") ""
 k1=WR(1) * abs(dot_product(dvector,VR(:,1)))
 k2=WR(2) * abs(dot_product(dvector,VR(:,2)))
 k3=WR(3) * abs(dot_product(dvector,VR(:,3)))
 kf=(k1+k2+k3)
 write (*,"(a)") ""
 write (*,"(2a,I0,2a,I0,a,F10.2,a)") "Approximate force constant between atom " , trim(atmlna(atm1*3-2)), atmlnb(atm1*3-2), &
       " and atom ", trim(atmlna(atm2*3-2)), atmlnb(atm2*3-2), " : ", kf, " a.u."
 write (*,"(a)") ""
end program matrix
