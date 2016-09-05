program sem_forces
 implicit none
 integer, parameter :: dp=kind(0.0d0)         ! doublz precision
 integer, parameter :: sp=kind(0.0)           ! single precision
 
 real(dp), parameter           :: BA = 0.529177249 ! Bohr to Angstrom
 real(dp), parameter           :: HKC = 627.5095   ! Hartree to kcal/mol
 integer, parameter            :: N = 3 
 integer, parameter            :: LDA = N, LDVL = N, LDVR = N
 integer, parameter            :: LWMAX = 1000
 integer                       :: INFO, LWORK 
 real(dp), dimension(3,3)      :: AM, VR, VL
 real(dp), dimension(3)        :: WR, WI, VL1, VL2, VL3
 real(dp), allocatable         :: WORK (:)
 
 real(dp), allocatable         :: atmat(:,:), atmcrd (:,:)
 integer, allocatable          :: atmlnb (:), atmlan (:)
 character(len=2), allocatable :: atmlna (:) 
 real(dp), dimension(3)        :: dvector
 real(dp)                      :: dist, kf
 integer                       :: atm1, atm2, smx, smy, nbatm, matsize
 integer                       :: i, j

 external DGEEV
 
!READ FILE 
 open(unit=50, file="data", status='old', access='sequential', form='formatted', action='read')
   read(50, *) nbatm 
   matsize = nbatm * 3
  allocate(atmcrd(5,nbatm))
  allocate(atmlnb(matsize))
  allocate(atmlan(matsize))
  allocate(atmlna(matsize))
  allocate(atmat(matsize,matsize))
   read(50, *) ! advance one line
   read(50, *) atmcrd ! fill the coordinates matrix (ATMLNB,ATMLAN,X,Y,Z)
   read(50, *) !advance one line
   read(50, *) atmlnb ! fill the atom number list (ATM1,ATM1,ATM1,ATM2,ATM2,ATM2,..)
   read(50, *) atmlan ! fill the atomic number list (same)
   read(50, *) atmlna ! fill the atom name list (same)
   read(50, *) ! advance one line
   read(50, *) atmat
 close(50)
!END READ FILE
!SELECT ATOMS AND FILL INTERATOMIC FORCE MATRIX
  write (*, " (a) ") "This is the list of atoms:"
  do i = 1,nbatm
   write (*, " (a,I0,a) " ,advance="NO") trim(atmlna(i*3-2)), atmlnb((i*3-2)), " "
  end do
  write (*, " (a) ")
  write (*, " (a) ")
  write (*, " (a) ") "Select atom number one:"
   read(*, *) atm1
  write (*, " (a) ")
  write (*, " (a) ") "Select atom number two:"
   read(*, *) atm2
  smx=((atm1*3)-3)
  smy=((atm2*3)-3)
  do i=1,3
   do j=1,3
    AM(i,j) = atmat( (smx+i), (smy+j) )
   end do
  end do
!PRINT FORCE MATRIX
  AM = -1*AM
 write (*, " (a) ")
 write (*, " (2a,I0,2a,I0,a) ") "Interatomic Cartesian Force Matrix (a.u.) between atom ", &
                                trim(atmlna(atm1*3-2)), atmlnb(atm1*3-2), " and atom ", & 
                                trim(atmlna(atm2*3-2)), atmlnb(atm2*3-2), " :"
 write (*, " (a) ")
 do j=1,3
  do i=1,3
   write(*, " (E14.6E2,a) " , advance="no") (AM(i,j)), ""
   end do
   write (*, " (a) ")
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
 if( INFO.GT.0 ) then
  write(*, *) 'The algorithm failed to compute eigenvalues.'
   stop
  end if
 deallocate(WORK)
!END CALL FOR INTEL MKL
!PRINT Eigen Values & Vectors
 write(*, " (a) ")
 write(*, " (a) ") "Eigen Values (real) (a.u.):"
 write(*, " (a) ")
 do i=1,3
  write(*, "(E14.5E2,a)", advance="no") (WR(i)), ""
 end do
 write(*, " (a) ")
 write(*, " (a) ")
 write(*, " (a) ") "Eigen Values (comp) (a.u.):"
 write(*, " (a) ")
 do i=1,3
  write(*, " (E14.5E2,a) ", advance="no") (WI(i)), ""
 end do
 write(*, " (a) ")
 write(*, " (a) ")
 write(*, " (a) ") "Right Eigen Vectors:"
 write(*, " (a) ")
 do j=1,3
  do i=1,3
   write(*, " (E14.6E2,a) ", advance="no") VR(i,j), ""
  end do
  write(*, " (a) ")
 end do
 write(*, " (a) ")
 write(*, " (a) ") "Left Eigen Vectors:"
 write(*, " (a) ")
 do j=1,3
  do i=1,3
   write(*, " (E14.6E2,a) ", advance="no") VL(i,j), ""
  end do
  write(*, " (a) ")
 end do
!END PRINT Eigen Values & Vectors
!CALCULATE distance vector, distance and unit vector
 do i=1,3
 dvector(i) = (atmcrd(i+2,atm2) - atmcrd(i+2,atm1)) ! distance vector
 end do
 dist = SQRT(abs(dvector(1)**2)+abs(dvector(2)**2)+abs(dvector(3)**2)) ! distance 
 do i=1,3
 dvector(i) = dvector(i)/abs(dist) ! distance unit vector
 end do
!END CALCULATE distance vector, distance and unit vector
!PRINT DISTNCE
 write(*, " (a) ")
 write(*, " (a,F6.2,a) ", advance='NO') "Interatomic distance: ", (dist * BA), " Å"
 write(*, " (a) ")
!END PRINT DISTANCE
!CALCULATE RFC AND PRINT
 kf = 0.0_dp
 do i=1,3
  kf = kf + WR(i) * abs(dot_product(VR(:,i), dvector))
 end do 
 write(*, " (a) ")
 write(*, " (2a,I0,2a,I0,a,F7.1,a) ") "Force constant (k(r-r0)^2) using right EV between atom " ,&
                                      trim(atmlna(atm1*3-2)), atmlnb(atm1*3-2), " and atom ",& 
                                      trim(atmlna(atm2*3-2)), atmlnb(atm2*3-2), " : ",&
                                      kf*(HKC/(BA**2))*0.5, " kcal*mol^-1*Â^-2"
!END CALCULATE RFC AND PRINT
!CALCULATE LFC AND PRINT
 kf = 0.0_dp
 do i=1,3
  kf = kf + WR(i) * abs(dot_product(VL(:,i), dvector))
 end do
 write(*, " (a) ")
 write(*, " (2a,I0,2a,I0,a,F7.1,a) ") "Force constant (k(r-r0)^2) using left EV between atom " ,&
                                      trim(atmlna(atm1*3-2)), atmlnb(atm1*3-2), " and atom ",& 
                                      trim(atmlna(atm2*3-2)), atmlnb(atm2*3-2), " : ",&
                                      kf*(HKC/(BA**2))*0.5, " kcal*mol^-1*Å^-2."
!END CALCULATE LFC AND PRINT
 write(*, " (a) ")
  deallocate(atmcrd)
  deallocate(atmlnb)
  deallocate(atmlan)
  deallocate(atmlna)
  deallocate(atmat)
 stop
end program sem_forces
