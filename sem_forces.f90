! Seminario method for force constant (bonds, angles and torsions)
! Reference: Seminario, J. M. Int. J. Quantum Chem. 1996, 60, 1271.
! Original Author: Rolf David
! License: MIT
! Date: 05/09/2016
! Lastest modification: 23/09/2019
! Version: 2.0.1

! Mathmodule for the cross product between two vectors
module mathmodule
    implicit none
contains
    function cross(a,b)
        implicit none
        integer, parameter :: dp=kind(0.0d0)         ! double precision
        integer, parameter :: sp=kind(0.0)           ! single precision
        real(dp), dimension(3), intent(in) :: a,b
        real(dp), dimension(3)             :: cross
        cross(1)=a(2) * b(3) - a(3) * b(2)
        cross(2)=a(3) * b(1) - a(1) * b(3)
        cross(3)=a(1) * b(2) - a(2) * b(1)
    end function cross
end module mathmodule

! Main program
program sem_forces
    use mathmodule
    implicit none
    character(len=64) :: filename

    integer, parameter :: dp=kind(0.0d0)         ! double precision
    integer, parameter :: sp=kind(0.0)           ! single precision

    real(dp), parameter           :: BA = 0.529177249      ! Bohr to Angstrom
    real(dp), parameter           :: HKC = 627.509469      ! Hartree to kcal/mol
    real(dp), parameter           :: pi = 3.1415926535897932384

    ! Variables for DGEEV
    integer, parameter            :: N = 3
    integer, parameter            :: LDA = N, LDVL = N, LDVR = N
    integer, parameter            :: LWMAX = 1000
    integer                       :: INFO, LWORK
    real(dp), dimension(3,3)      :: AM_AB, AM_CB, AM_DC, VRAB, VLAB, VRCB, VLCB, VRDC, VLDC, AM_BA, VRBA, VLBA
    real(dp), dimension(3)        :: WRAB, WIAB, WRCB, WICB, WRDC, WIDC, WRBA, WIBA
    real(dp), allocatable         :: WORK (:)

    ! Variables
    real(dp), allocatable         :: atmat(:,:), atmcrd (:,:)
    integer, allocatable          :: atmlnb (:), atmlan (:)
    character(len=2), allocatable :: atmlna (:)
    real(dp), dimension(3)        :: vecAB, vecBA, vecCB, vecNABC, vecPA, vecPC, vecBC, vecDC, vecNBCD, vecCD                        ! vectors
    real(dp)                      :: distAB, distBA, distCB, distNABC, distBC, distDC, distNBCD, distCD                                     ! distance/norm
    real(dp)                      :: angleABC, angleABCD, angleABCDsign                                                                           ! angle
    real(dp)                      :: kAB, kBA, kABavg, kCB, kABC, kBCD, kABC2, kBCD2, kABCD                                                     ! FC     
    integer                       :: atmA, atmB, atmC, atmD, smA, smB, smC, smD, nbatm, matsize
    integer                       :: i, j

    external DGEEV

    ! Read file
    CALL get_command_argument(1,filename)
    open(unit=50, file=filename, status='old', access='sequential', form='formatted', action='read')
    read(50, *) nbatm
    matsize = nbatm * 3
    allocate(atmcrd(nbatm,5))
    allocate(atmlnb(matsize))
    allocate(atmlan(matsize))
    allocate(atmlna(matsize))
    allocate(atmat(matsize,matsize))
    read(50, *) ! advance one line
    do i=1,nbatm
        read(50, *) atmcrd(i,:) ! fill the coordinates matrix (ATMLNB,ATMLAN,X,Y,Z)
    enddo
    read(50, *) ! advance one line
    read(50, *) atmlnb ! fill the atom number list (ATM1,ATM1,ATM1,ATM2,ATM2,ATM2,..)
    read(50, *) atmlan ! fill the atomic number list (same)
    read(50, *) atmlna ! fill the atom name list (same)
    read(50, *) ! advance one line
    read(50, *) atmat
    close(50)
    ! Ask user to select atoms
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
    write (*, " (a) ") "Put 0 if not needed (bond FC), put the atom number for bond angle FC"
    read(*, *) atmC
    write (*, " (a) ")
    write (*, " (a) ") "Select the fourth atom:"
    write (*, " (a) ") "Put 0 if not needed (bond FC, bond angle FC), put the atom number for dihedral angle FC"
    read(*, *) atmD

    ! Whatever the case load AB interatomic force constant matrix
    smA=((atmA * 3)-3)
    smB=((atmB * 3)-3)
    do i=1,3
        do j=1,3
            AM_AB(i,j) = atmat( (smA+i), (smB+j) )
        end do
    end do
    AM_AB = -1.0 * AM_AB
    ! And the BA interatomic FC matrix also
    do i=1,3
        do j=1,3
            AM_BA(i,j) = atmat( (smB+i), (smA+j) )
        end do
    end do
    AM_BA = -1.0 * AM_BA
!----------------------------------- BOND FORCE CONSTANTS ----------------------------------! 
    if (atmC .eq. 0 .and. atmD .eq. 0) then ! bond only

        ! Print force matrix for bond length
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

        ! Call for intel MKL or LAPACK
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_BA, LDA, WRBA, WIBA, VLBA, LDVL, VRBA, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_BA, LDA, WRBA, WIBA, VLBA, LDVL, VRBA, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if

        ! Calculate distance vector, its norm and normalize
        do i=1,3
            vecAB(i) = (atmcrd(atmB,i+2) - atmcrd(atmA,i+2)) ! AB vector
        end do
        do i=1,3
            vecBA(i) = (atmcrd(atmA,i+2) - atmcrd(atmB,i+2)) ! BA vector
        end do

        distAB = norm2(vecAB) ! Norm of AB
        distBA = norm2(vecBA) ! Norm of BA

        do i=1,3
            vecAB(i) = vecAB(i)/abs(distAB) ! AB unit vector
        end do
        do i=1,3
            vecBA(i) = vecBA(i)/abs(distBA) ! BA unit vector
        end do

        ! Print the distance
        write(*, " (a) ")
        write(*, " (a,F6.2,a,a,F6.2,a) ", advance='NO') "Interatomic distance: ", (distAB * BA), " Å", " | ", (distBA * BA), " Å"
        write(*, " (a) ")
        write(*, " (a) ")

        ! Calculate the bond FC and print it
        kAB = 0.0_dp
        do i=1,3
            kAB = kAB + WRAB(i) * abs(dot_product(VRAB(:,i), vecAB))
        end do
        do i=1,3
            kBA = kBA + WRBA(i) * abs(dot_product(VRBA(:,i), vecBA))
        end do
        write(*, " (a,F8.2,a,F8.2) ", advance='NO') "Bond FC: ", kAB*(HKC/(BA**2)), " | ", kBA*(HKC/(BA**2))
        write(*, " (a) ")
        write(*, " (a) ")
        kABavg = (kAB + kBA)/2
        write(*, " (2a,I0,2a,I0,a,F8.2,a,F6.2,a) ") "Ebond = (1/2)k(r-r0)^2 for bond ", &
            trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "-", &
            trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), " with k =", &
            kABavg*(HKC/(BA**2)), " kcal*mol^-1*Å^-2 and r0 =", &
            (distAB * BA), " Å"

        ! Deallocate and exit
        deallocate(atmcrd)
        deallocate(atmlnb)
        deallocate(atmlan)
        deallocate(atmlna)
        deallocate(atmat)
        stop

!----------------------------------- BOND ANGLE FORCE CONSTANTS ----------------------------------!
    else if (atmC .gt. 0 .and. atmD .eq. 0) then
        smC=((atmC * 3)-3)
        do i=1,3
            do j=1,3
                AM_CB(i,j) = atmat( (smC+i), (smB+j) )
            end do
        end do
        AM_CB = -1*AM_CB

        ! Call for intel MKL or LAPACK
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_CB, LDA, WRCB , WICB, VLCB, LDVL, VRCB, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_CB, LDA, WRCB, WICB, VLCB, LDVL, VRCB, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if

        ! Calculate both distance vectors, their norms and normalize
        do i=1,3
            vecAB(i) = (atmcrd(atmB,i+2) - atmcrd(atmA,i+2)) ! AB vector
            vecCB(i) = (atmcrd(atmB,i+2) - atmcrd(atmC,i+2)) ! CB vector
        end do
        distAB = norm2(vecAB) ! Norm of AB
        distCB = norm2(vecCB) ! Norm of CB
        do i=1,3
            vecAB(i) = vecAB(i)/abs(distAB) ! AB unit vector
            vecCB(i) = vecCB(i)/abs(distCB) ! CB unit vector
        end do

        ! Calculate the vector perpendicular to the plane ABC, its norm and normalize
        vecNABC = cross(vecCB,vecAB) ! NormABC vector
        distNABC=norm2(vecNABC) ! Norm of the NormABC vector
        do i=1,3
            vecNABC(i)=vecNABC(i)/abs(distNABC) ! NormABC unit vector
        enddo

        ! Calculate unit vector perpendiculr to AB and CB on the ABC plane
        vecPA = cross(vecNABC,vecAB)
        vecPC = cross(vecCB,vecNABC)

        ! Calculate the force contributions
        kAB = 0.0_dp
        kCB = 0.0_dp
        do i=1,3
            kAB = kAB + WRAB(i) * abs(dot_product(VRAB(:,i), vecPA))
            kCB = kCB + WRCB(i) * abs(dot_product(VRCB(:,i), vecPC))
        end do

        ! Calculate the bond angle FC
        kABC = 0.0_dp
        kABC = 1 / ((distAB**2)*(kAB)) + 1 / ((distCB**2)*(kCB))
        kABC = 1 / kABC

        ! Calculate angle
        angleABC = acos(dot_product(vecAB,vecCB))*180.0/pi

        ! Print the angle
        write(*, " (a) ")
        write(*, " (a,F7.1,a) ", advance='NO') "Angle equal to: ", angleABC, " °"
        write(*, " (a) ")

        ! Print the bond angle FC
        write(*, " (a) ")
        write(*, " (2a,I0,2a,I0,2a,I0,a,F7.1,a) ") "Bond angle FC ((1/2)k(T-T0)^2) for angle ", &
            trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
            trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), "|", &
            trim(atmlna(atmC*3-2)), atmlnb(atmC*3-2), " is equal to : ", &
            (kABC*(HKC)), " kcal*mol^-1*rad^-2"

        ! Deallocate and exit
        deallocate(atmcrd)
        deallocate(atmlnb)
        deallocate(atmlan)
        deallocate(atmlna)
        deallocate(atmat)
        stop

!----------------------------------- DIHEDRAL ANGLE FORCE CONSTANTS ----------------------------------!
    else if (atmC .gt. 0 .and. atmD .gt. 0) then
        smC=((atmC * 3)-3)
        smD=((atmD * 3)-3)

        do i=1,3
            do j=1,3
                AM_DC(i,j) = atmat( (smD+i), (smC+j) )
            end do
        end do
        AM_DC = -1*AM_DC
        ! Call for intel MKL or LAPACK
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_AB, LDA, WRAB, WIAB, VLAB, LDVL, VRAB, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if
        LWORK = -1
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_DC, LDA, WRDC , WIDC, VLDC, LDVL, VRDC, LDVR, WORK, LWORK, INFO)
        LWORK = MIN (LWMAX, INT(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        call DGEEV ('N', 'V', N, AM_DC, LDA, WRDC , WIDC, VLDC, LDVL, VRDC, LDVR, WORK, LWORK, INFO)
        deallocate(WORK)
        if( INFO > 0 ) then
            write(*, *) 'The algorithm failed to compute eigenvalues.'
            stop
        end if

        ! Calculate all distance vectors, their norms and normalize
        do i=1,3
            vecAB(i) = (atmcrd(atmB,i+2) - atmcrd(atmA,i+2)) ! AB vector
            vecBA(i) = (atmcrd(atmA,i+2) - atmcrd(atmB,i+2)) ! BA vector
            vecCB(i) = (atmcrd(atmB,i+2) - atmcrd(atmC,i+2)) ! CB vector
            vecBC(i) = (atmcrd(atmC,i+2) - atmcrd(atmB,i+2)) ! BC vector
            vecDC(i) = (atmcrd(atmC,i+2) - atmcrd(atmD,i+2)) ! DC vector
            vecCD(i) = (atmcrd(atmD,i+2) - atmcrd(atmC,i+2)) ! CD vector
        end do

        distAB = norm2(vecAB) ! Norm of AB
        distBA = norm2(vecBA) ! Norm of BA
        distCB = norm2(vecCB) ! Norm of CB
        distBC = norm2(vecBC) ! Norm of BC
        distDC = norm2(vecDC) ! Norm of DC
        distCD = norm2(vecCD) ! Norm of CD
        do i=1,3
            vecAB(i) = vecAB(i)/abs(distAB) ! AB unit vector
            vecBA(i) = vecBA(i)/abs(distBA) ! BA unit vector
            vecCB(i) = vecCB(i)/abs(distCB) ! CB unit vector
            vecBC(i) = vecBC(i)/abs(distBC) ! BC unit vector
            vecDC(i) = vecDC(i)/abs(distDC) ! DC unit vector
            vecCD(i) = vecCD(i)/abs(distCD) ! CD unit vector
        end do

        ! Calculate the vector perpendicular to the plane ABC, its norm and normalize
        vecNABC = cross(vecCB,vecAB) ! NormABC vector
        distNABC=norm2(vecNABC) ! Norm of the NormABC vector
        do i=1,3
            vecNABC(i)=vecNABC(i)/abs(distNABC) ! NormABC unit vector
        enddo

        ! Calculate the vector perpendicular to the plane BCD, its norm and normalize
        vecNBCD = cross(vecDC,vecBC) ! NormBCD vector
        distNBCD=norm2(vecNBCD) ! Norm of the NormBCDC vector
        do i=1,3
            vecNBCD(i)=vecNBCD(i)/abs(distNBCD) ! NormBCD unit vector
        enddo

        ! Calculate the force contributions

        kABC = 0.0_dp
        kBCD = 0.0_dp

        do i=1,3
            kABC = kABC + WRAB(i) * abs(dot_product(VRAB(:,i), vecNABC))
            kBCD = kBCD + WRDC(i) * abs(dot_product(VRDC(:,i), vecNBCD))
        enddo

        kABC2 = 0.0_dp
        kABC2 = (distAB**2)*(norm2(cross(vecAB,vecBC)))**2*kABC
        kBCD2 = 0.0_dp
        kBCD2 = (distDC**2)*(norm2(cross(vecBC,vecCD)))**2*kBCD
        kABCD = 1 / kABC2 + 1 / kBCD2
        kABCD = 1 / kABCD

        ! Calculate sign and the angle

        vecNABC = cross(vecCB,vecAB) ! NormABC vector
        distNABC=norm2(vecNABC) ! Norm of the NormABC vector
        do i=1,3
            vecNABC(i)=vecNABC(i)/abs(distNABC) ! NormABC unit vector
        enddo
       
        angleABCDsign = dot_product(vecNABC,vecBA)
        angleABCD = acos((dot_product(cross(vecBA,vecBC),cross(vecCB,vecCD)))&
        /(sin(acos(dot_product(vecAB,vecCB)))*sin(acos(dot_product(vecBC,vecCD)))))*180.0/pi
        angleABCD = sign(angleABCD,angleABCDsign)
        
        ! Print the angle
        write(*, " (a) ")
        write(*, " (a,F7.1,a) ", advance='NO') "Angle equal to: ", angleABCD, " °"
        write(*, " (a) ")

        ! Print the bond angle FC
        write(*, " (a) ")
        write(*, " (2a,I0,2a,I0,2a,I0,2a,I0,a,F7.1,a) ") "Dihedral angle FC ((1/2)k(P-P0)^2) for angle ", &
            trim(atmlna(atmA*3-2)), atmlnb(atmA*3-2), "|", &
            trim(atmlna(atmB*3-2)), atmlnb(atmB*3-2), "|", &
            trim(atmlna(atmC*3-2)), atmlnb(atmC*3-2), "|", &
            trim(atmlna(atmD*3-2)), atmlnb(atmD*3-2), " is equal to : ", &
            (kABCD*(HKC)), " kcal*mol^-1*rad^-2"

    else 
        stop
    end if
    stop
end program sem_forces
