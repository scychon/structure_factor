!****************************************************************
! 
! this module reads topology file and save it to the type topfile
!
!****************************************************************

MODUlE topol


    implicit none
    private
 
!   *** topfile type
!   box     - triclinic pbc box of the configuration
!   NATOMS  - number of atoms in the configuration.
!   pos     - positions read in (3,NATOMS)
!   prec    - precision of the coordinates read in
!   STEP    - step number of configuration.
!   STAT    - status of operation. 0 = good
!   time    - time of the configuration
!   xd      - pointer from libxdrfile.
 
!   Should always call init first. Then call read in a loop and do your
!   calculations. After the loops call close.
 
    type, public :: molecule
      Character(5) :: moltypename
      integer      :: numatom, nummol
      integer, dimension(3) :: molPlaneAtomIdxs
      real*8       :: molarmass, molcharge
      real*8, allocatable :: atommass(:)
      real*8, allocatable :: atomcharge(:)
      integer, allocatable :: atomTypeIdx(:)
      character*4, allocatable :: atomtype(:)
    end type
 
    type, public :: topfile
      type(molecule), allocatable :: moltype(:)
      character*4, allocatable :: atomtype(:)
      integer, allocatable :: pairlist(:,:)
!      type(molecule), pointer, allocatable :: sysmol(:)
      integer      :: numsysmol, numsysatom, nummoltype, numatomtype, nsize_mol, nsize_atom
    contains
      procedure :: init => init_top
    end type

 
 
 
contains
 
    ! our wrappers for the topfile class
    subroutine init_top(top,filename_in)
 
        implicit none
        class(topfile), intent(inout) :: top
        character (len=*), intent(in) :: filename_in
        character (len=206) :: filename
        character(50):: line
        character(5):: atomname,molname,atname1,atname2
        logical :: ex
        integer :: i,j,k,nmoltype,nmol,natom,nmolsys,natomtype
        integer :: nUniqueAtomTypes,nGrPairs,attype1,attype2
        integer,dimension(3) :: idxs
        real*8  :: mass, molmass, charge
        character*4 :: typename
        character*4, allocatable :: tempAtomType(:)
        integer, allocatable :: arrGrPairList(:,:)
        logical :: bNewAtomType, bDrude
 
        inquire(file=trim(filename_in),exist=ex)
        write(6,'(A,A)') 'reading topology file ',trim(filename_in)
        write(6,*) ''
 
        if (ex .eqv. .false.) then
            write(0,*)
            write(0,'(a)') " Error: "//trim(filename_in)//" does not exist."
            write(0,*)
            stop
        end if
 
        ! Set the file name to be read in for C.
        filename = trim(filename_in)
        OPEN(unit=7,file=filename,status='old')
 
        ! Get number of moltypes
        read(7, '(A)') line
        read(7,*) nmoltype
        top % nummoltype = nmoltype

        allocate(top % moltype(nmoltype))

        ! Read atomic masses for each molecule type
        natomtype = 0
        bDrude = .false.
        do i=1, nmoltype
          read(7,'(A)') molname
          read(7,*) natom
          top % moltype(i) % moltypename = molname
          top % moltype(i) % molarmass = 0
          top % moltype(i) % molcharge = 0
          top % moltype(i) % numatom = natom
          allocate(top % moltype(i) % atomtype(natom))
          allocate(top % moltype(i) % atomTypeIdx(natom))
          allocate(top % moltype(i) % atommass(natom))
          allocate(top % moltype(i) % atomcharge(natom))

          do j=1, natom
            read(7,*) atomname, mass, charge
            top % moltype(i) % atomtype(j) = atomname
            top % moltype(i) % atommass(j) = mass
            top % moltype(i) % atomcharge(j) = charge
            top % moltype(i) % molarmass = top % moltype(i) % molarmass + mass
            top % moltype(i) % molcharge = top % moltype(i) % molcharge + charge
            if ( mass .le. 0) bDrude = .true.
          enddo

          natomtype = natomtype + natom
        enddo

        ! store unique atomtypes in the topology
        allocate(tempAtomType(natomtype))
        nUniqueAtomTypes = 0
        do i=1, nmoltype
          do j=1, top % moltype(i) % numatom
            typename = top % moltype(i) % atomtype(j)
            ! check if it's shell atomtype
            if (top % moltype(i) % atommass(j) .le. 0) then
                top % moltype(i) % atomtypeidx(j) = -1
                cycle
            endif
            bNewAtomType = .true.
            ! check whether the atomtype is already defined
            do k=1, nUniqueAtomTypes
              if (typename .eq. tempAtomType(k)) then
                bNewAtomType = .false.
                top % moltype(i) % atomtypeidx(j) = k
                exit
              endif
            enddo
            if (bNewAtomType) then
              nUniqueAtomTypes = nUniqueAtomTypes + 1
              tempAtomType(nUniqueAtomTypes) = typename
              top % moltype(i) % atomtypeidx(j) = nUniqueAtomTypes
            endif
          enddo
        enddo

        allocate(top % atomtype(nUniqueAtomTypes))
        top % atomtype = tempAtomType(:nUniqueAtomTypes)


        ! Read total # of molecules in system
        read(7, '(A)') line
        read(7, *) nmol
        top % numsysmol = nmol
!        allocate(top % sysmol(nmol))
        nmolsys = nmol

        ! Read # of molecules for each moltype
        k = 0
        natom = 0
        do i=1, nmoltype
          read(7,*) molname, nmol
          top % moltype(i) % nummol = nmol
          natom = natom + nmol * top % moltype(i) % numatom
          write(*,*) 'molcharge', top % moltype(i) % molcharge
          write(*,*) 'molarmass', top % moltype(i) % molarmass
          k = k + nmol
!          do j=1, nmol
!            k = k+1
!            top % sysmol(k) => top % moltype(i)
!          enddo
        enddo
        top % numsysatom = natom

        if(k .ne. nmolsys) then
            write(0,*)
            write(0,'(a)') " Error: number of molecules of each types does not sum up to the total number of molecules in system."
            write(0,*)
            stop
        end if

        ! rotational ACF calculation
        ! Read three atom idxs of molecules which compose the molecular plane
        read(7, '(A)') line
        do i=1, nmoltype
          read(7,*) molname, idxs
          top % moltype(i) % molPlaneAtomIdxs = idxs
        enddo
        
        ! g(r) calculation pairs
        ! Read three atom idxs of molecules which compose the molecular plane
        read(7, '(A)') line
        read(7, *) nGrPairs
        allocate(arrGrPairList(nGrPairs,2))
        allocate(top % pairlist(nGrPairs,2))
        do i=1, nGrPairs
          read(7,*) atname1, atname2
          attype1 = 0
          attype2 = 0
          ! check whether the atomtype is already defined
          do k=1, nUniqueAtomTypes
            if (atname1 .eq. top % atomtype(k)) then
              attype1 = k
            endif
            if (atname2 .eq. top % atomtype(k)) then
              attype2 = k
            endif
          enddo
          if ((attype1 .gt. 0) .and. (attype2 .gt. 0)) then
            arrGrPairList(i,:) = (/ attype1, attype2 /)
            write(0,*) "gr pair list ",atname1," ",atname2," ",attype1," ",attype2
          else
            write(0,*)
            write(0,*) " Error: atom types in  ",i,"th gr pair list ",atname1," or ",atname2, " is not defined in the topology "
            write(0,*)
          endif
        enddo
        top % pairlist = arrGrPairList

        close(7)

    end subroutine init_top

!    recursive subroutine read_itp(top, itpfile)
!        implicit none
!        class(topfile), intent(inout) :: top
!        character (len=*), intent(in) :: itpfile
!        character (len=206) :: filename
!        character(50):: line
!        character(5):: atomname,molname
!        logical :: ex
!        integer :: i,j,k,nmoltype,nmol,natom,nmolsys
!        integer :: readstat
!        real*8  :: mass, molmass
!        real*8, allocatable :: temp(:)
! 
!        inquire(file=trim(itpfile),exist=ex)
! 
!        if (ex .eqv. .false.) then
!            write(0,*)
!            write(0,'(a)') " Error: "//trim(itpfile)//" does not exist."
!            write(0,*)
!            stop
!        end if
! 
!        ! Set the file name to be read in for C.
!        filename = trim(itpfile)
!        OPEN(unit=7,file=filename,status='old')
!
!        readstat = 0 
!
!        ! Read until end of file
!        do while (readstat .ge. 0)
!          read(15, '(A)',IOSTAT=readstat) line
!
!          !remove comments
!
!          !separate line into variables
!
!          !include statements
!
!          !define statements
!
!          !section titles
!
!
!          !read each section parameters
!
!
!          !
!          if (readstat .ge. 0) then
!            i=i+1
!            if(readstat .eq. 0) then
!              j=j+1
!              if(size(re_to_e_list) .lt. j) then
!                allocate(temp(size(re_to_e_list)+nsize))
!                temp=0
!                temp(:size(re_to_e_list))=re_to_e_list
!                deallocate(re_to_e_list)
!                allocate(re_to_e_list(size(temp)))
!                re_to_e_list=temp
!                deallocate(temp)
!              endif
!              re_to_e_list(j)=re_to_e
!              if(re_to_e .gt. rmax) rmax=re_to_e
!            else
!              write(*,*) 'something is wrong with %ith data', i
!            endif
!          else
!            write(*,*) 'reading file finished'
!          endif
!        enddo
!
!        ! Get number of moltypes
!        read(7, '(A)') line
!    end subroutine

end module topol
