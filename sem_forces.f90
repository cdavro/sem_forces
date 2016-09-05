! Seminario method for bond force constant (length and angle)
! Reference: Seminario, J. M. Int. J. Quantum Chem. 1996, 60, 1271.
! Author: Rolf David
! Date: 05/09/2016
! Lastest modification: 05/09/2016
! Version: 1.1.0

program sem_forces
 implicit none
 integer, parameter :: dp=kind(0.0d0)         ! double precision
 integer, parameter :: sp=kind(0.0)           ! single precision
 
 real(dp), parameter           :: BA = 0.529177249 ! Bohr to Angstrom
 real(dp), parameter           :: HKC = 627.5095   ! Hartree to kcal/mol
 integer, parameter            :: N = 3 
 integer, parameter            :: LDA = N, LDVL = N, LDVR = N
 integer, parameter            :: LWMAX = 1000
 integer                       :: INFO, LWORK 
 real(dp), dimension(3,3)      :: AM_AB, AM_BC, VRAB, VLAB, VRBC, VLBC
 real(dp), dimension(3)        :: WRAB, WIAB, WRBC, WIBC
 real(dp), allocatable         :: WORK (:)
 
 real(dp), allocatable         :: atmat(:,:), atmcrd (:,:)
 integer, allocatable          :: atmlnb (:), atmlan (:)
 character(len=2), allocatable :: atmlna (:) 
 real(dp), dimension(3)        :: vecAB, vecBC, vecNABC, vecPA, vecPC ! vectors
 real(dp)                      :: distAB, distBC, distNABC ! distance
 real(dp)                      :: kRAB, kRBC, kLAB, kLBC, kavgAB, kavgBC, kRABC, kLABC, kavgABC ! force constant
 integer                       :: atmA, atmB, atmC, smA, smB, smC, nbatm, matsize
 integer                       :: i, j

 external DGEEV
 
!Read file 
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

!Ask user to select atoms
 write (*, " (a) ") "This is the list of atoms:"
 do i = 1,nbatm
  write (*, " (a,I0,a) " ,advance="NO") trim(atmlna(i*3-2)), atmlnb((i*3-2)), " "
 end do
 write (*, " (a) ")
 write (*, " (a) ")
 write (*, " (a) ") "Select the first atom:"
  read(*, *) atmA
 write (*, " (a) ")
 write (*, " (a) ") "Select the second atom:"
  read(*, *) atmB
 write (*, " (a) ")
 write (*, " (a) ") "Select the third atom:"
 write (*, " (a) ") "Put 0 if you want the bond lenght force constant"
  read(*, *) atmC 
  
!Whatever the case load AB interatomic force constant matrix
 smA=((atmA * 3)-3)
 smB=((atmB * 3)-3)
 do i=1,3
  do j=1,3
   AM_AB(i,j) = atmat( (smA+i), (smB+j) )
  end do
 end do
 AM_AB = -1 * AM_AB

!Test if bond length or bond angle ?
 if (atmC == 0) then !bond

!Print force matrix for bond length
  write (*, " (a) ")
  write (*, " (2a,I0,2a,I0,a) ") "Interatomic Cartesian Force Matrix (a.u.) between atom ", &
                                 trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), " and atom ", & 
                                 trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), " :"
  write (*, " (a) ")
  do j=1,3
   do i=1,3
    write(*, " (E14.6E2,a) " , advance="no") (AM_AB(i,j)), ""
   end do
   write (*, " (a) ")
  end do

!Call for intel MKL or LAPACK
  LWORK = -1
  allocate(WORK(LWORK))
  call DGEEV ('V', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
  LWORK = MIN (LWMAX, INT(WORK(1)))
  deallocate(WORK)
  allocate(WORK(LWORK)) 
  call DGEEV ('V', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
  deallocate(WORK)
  if( INFO > 0 ) then
   write(*, *) 'The algorithm failed to compute eigenvalues.'
   stop
  end if

!Print Eigen Values & Vectors
  write(*, " (a) ")
  write(*, " (a) ") "Eigen Values (real) (a.u.):"
  write(*, " (a) ")
  do i=1,3
   write(*, "(E14.5E2,a)", advance="no") (WRAB(i)), ""
  end do
  write(*, " (a) ")
  write(*, " (a) ")
  write(*, " (a) ") "Eigen Values (comp) (a.u.):"
  write(*, " (a) ")
  do i=1,3
   write(*, " (E14.5E2,a) ", advance="no") (WIAB(i)), ""
  end do
  write(*, " (a) ")
  write(*, " (a) ")
  write(*, " (a) ") "Right Eigen Vectors:"
  write(*, " (a) ")
  do j=1,3
   do i=1,3
    write(*, " (E14.6E2,a) ", advance="no") VRAB(i,j), ""
   end do
   write(*, " (a) ")
  end do
  write(*, " (a) ")
  write(*, " (a) ") "Left Eigen Vectors:"
  write(*, " (a) ")
  do j=1,3
   do i=1,3
    write(*, " (E14.6E2,a) ", advance="no") VLAB(i,j), ""
   end do
   write(*, " (a) ")
  end do

!Calculate distance vector, its norm and normalize
  do i=1,3
   vecAB(i) = (atmcrd(i+2,atmB) - atmcrd(i+2,atmA)) ! distance vector
  end do
  distAB = SQRT(abs(vecAB(1)**2)+abs(vecAB(2)**2)+abs(vecAB(3)**2)) ! norm
  do i=1,3
   vecAB(i) = vecAB(i)/abs(distAB) ! distance unit vector
  end do

!Print the distance
  write(*, " (a) ")
  write(*, " (a,F6.2,a) ", advance='NO') "Interatomic distance: ", (distAB * BA), " Å"
  write(*, " (a) ")
  write(*, " (a) ")

!Calculate the RFC and print it
  kRAB = 0.0_dp
  do i=1,3
   kRAB = kRAB + WRAB(i) * abs(dot_product(VRAB(:,i), vecAB))
  end do 
  write(*, " (2a,I0,2a,I0,a,F7.1,a) ") "Bond length force constant (k(r-r0)^2) using right EV for bond ", &
                                       trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
                                       trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), " is equal to : ", &
                                       kRAB*(HKC/(BA**2))*0.5, " kcal*mol^-1*Â^-2"

!Calculate the LFC and print it
  kLAB = 0.0_dp
  do i=1,3
   kLAB = kLAB + WRAB(i) * abs(dot_product(VLAB(:,i), vecAB))
  end do
  write(*, " (2a,I0,2a,I0,a,F7.1,a) ") "Bond length force constant (k(r-r0)^2) using  left EV for bond " , &
                                       trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", & 
                                       trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), " is equal to : ", &
                                       kLAB*(HKC/(BA**2))*0.5, " kcal*mol^-1*Å^-2."

!Calculate the average of the two and print it 
  kavgAB = (kRAB + kLAB) * 0.5
  write(*, " (a) ")
  write(*, " (2a,I0,2a,I0,a,F7.1,a) ") "Averaged bond lenght force constant (k(r-r0)^2) for bond " ,&
                                       trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|",&
                                       trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), " is equal to : ",&
                                       kavgAB*(HKC/(BA**2))*0.5, " kcal*mol^-1*Å^-2."
  write(*, " (a) ")

!Deallocate and exit 
  deallocate(atmcrd)
  deallocate(atmlnb)
  deallocate(atmlan)
  deallocate(atmlna)
  deallocate(atmat)
  stop

!If a bond angle force constant is requested
 else
  smC=((atmC * 3)-3)
  do i=1,3
   do j=1,3
    AM_BC(i,j) = atmat( (smB+i), (smC+j) )
   end do
  end do
  AM_BC = -1*AM_BC

!Call for intel MKL or LAPACK
  LWORK = -1
  allocate(WORK(LWORK))
  call DGEEV ('V', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
  LWORK = MIN (LWMAX, INT(WORK(1)))
  deallocate(WORK)
  allocate(WORK(LWORK)) 
  call DGEEV ('V', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
  deallocate(WORK)
  if( INFO > 0 ) then
   write(*, *) 'The algorithm failed to compute eigenvalues.'
   stop
  end if
  LWORK = -1
  allocate(WORK(LWORK))
  call DGEEV ('V', 'V', N, AM_BC, LDA, WRBC , WIBC, VLBC, LDVL, VRBC, LDVR, WORK, LWORK, INFO)
  LWORK = MIN (LWMAX, INT(WORK(1)))
  deallocate(WORK)
  allocate(WORK(LWORK)) 
  call DGEEV ('V', 'V', N, AM_BC, LDA, WRBC, WIBC, VLBC, LDVL, VRBC, LDVR, WORK, LWORK, INFO)
  deallocate(WORK)
  if( INFO > 0 ) then
   write(*, *) 'The algorithm failed to compute eigenvalues.'
   stop
  end if
  
!Calculate both distance vectors, their norms and normalize
  do i=1,3
   vecAB(i) = (atmcrd(i+2,atmB) - atmcrd(i+2,atmA)) ! distance vector
   vecBC(i) = (atmcrd(i+2,atmC) - atmcrd(i+2,atmB)) !
  end do
  distAB = SQRT(abs(vecAB(1)**2)+abs(vecAB(2)**2)+abs(vecAB(3)**2)) ! norm
  distBC = SQRT(abs(vecBC(1)**2)+abs(vecBC(2)**2)+abs(vecBC(3)**2)) ! 
  do i=1,3
   vecAB(i) = vecAB(i)/abs(distAB) ! unit vector
   vecBC(i) = vecBC(i)/abs(distBC) !
  end do
 
!Calculate the vector perpendicular to the plane ABC, its norm and normalize
  vecNABC(1)=((vecAB(2) * vecBC(3)) - (vecAB(3) * vecBC(2))) ! vector 
  vecNABC(2)=((vecAB(3) * vecBC(1)) - (vecAB(1) * vecBC(3))) !
  vecNABC(3)=((vecAB(1) * vecBC(2)) - (vecAB(2) * vecBC(1))) !
  distNABC=SQRT(abs(vecNABC(1)**2)+abs(vecNABC(2)**2)+abs(vecNABC(3)**2)) ! norm
  do i=1,3
   vecNABC(i)=vecNABC(i)/abs(distNABC) ! normalize
  enddo
 
!Calculate unit vector perpendicaulr to AB and BC on the ABC plane
  vecPA(1)=((vecNABC(2) * vecAB(3)) - (vecNABC(3) * vecAB(2))) ! First one
  vecPA(2)=((vecNABC(3) * vecAB(1)) - (vecNABC(1) * vecAB(3))) !
  vecPA(3)=((vecNABC(1) * vecAB(2)) - (vecNABC(2) * vecAB(1))) !
  vecPC(1)=((vecNABC(2) * vecBC(3)) - (vecNABC(3) * vecBC(2))) ! Second
  vecPC(2)=((vecNABC(3) * vecBC(1)) - (vecNABC(1) * vecBC(3))) !
  vecPC(3)=((vecNABC(1) * vecBC(2)) - (vecNABC(2) * vecBC(1))) !

!Calculate the force contributions right then left  
  do i=1,3
   kRAB = kRAB + WRAB(i) * abs(dot_product(VRAB(:,i), vecPA))
   kRBC = kRBC + WRBC(i) * abs(dot_product(VRBC(:,i), vecPC))
  end do
  do i=1,3
   kLAB = kLAB + WRAB(i) * abs(dot_product(VLAB(:,i), vecPA))
   kLBC = kLBC + WRBC(i) * abs(dot_product(VLBC(:,i), vecPC))
  end do

!Calculate the bond angle force constant, right, left and average
  kRABC = 1 / ((distAB**2)*(kRAB)) + 1 / ((distBC**2)*(kRBC))
  kRABC = 1 / kRABC
  kLABC = 1 / ((distAB**2)*(kLAB)) + 1 / ((distBC**2)*(kLBC))
  kLABC = 1 / kLABC 
  kavgABC = (kRABC + kLABC) * 0.5

!Calculate angle 
!TBD

!Print the angle
  write(*, " (a) ")
  write(*, " (a) ", advance='NO') "Angle equal to: ", " °"
  write(*, " (a) ")
 
!Print all  
  write(*, " (a) ")
  write(*, " (2a,I0,2a,I0,2a,I0,a,F7.1,a) ") "Bond angle force constant (k(r-r0)^2) using right EV for angle ", &
                                             trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
                                             trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), "|", &
                                             trim(atmlna(atmC*3-2)), atmlnb(atmC*3-2), " is equal to : ", &
                                             (kRABC*(HKC)*0.5), " kcal*mol^-1*rad^-2"
  write(*, " (2a,I0,2a,I0,2a,I0,a,F7.1,a) ") "Bond angle force constant (k(r-r0)^2) using  left EV for angle ", &
                                             trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
                                             trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), "|", &
                                             trim(atmlna(atmC*3-2)), atmlnb(atmC*3-2), " is equal to : ", &
                                             (kLABC*(HKC)*0.5), " kcal*mol^-1*rad^-2"
  write(*, " (a) ")
  write(*, " (2a,I0,2a,I0,2a,I0,a,F7.1,a) ") "Averaged bond angle force constant (k(r-r0)^2) for angle ", &
                                             trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
                                             trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), "|", &
                                             trim(atmlna(atmC*3-2)), atmlnb(atmC*3-2), " is equal to : ", &
                                             (kavgABC*(HKC)*0.5), " kcal*mol^-1*rad^-2"
  write(*, " (a) ")

!Deallocate and exit 
  deallocate(atmcrd)
  deallocate(atmlnb)
  deallocate(atmlan)
  deallocate(atmlna)
  deallocate(atmat)
  stop

!Exit 
 end if
 stop
end program sem_forces

