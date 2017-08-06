!  XDR Fortran Interface xtc Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
 
program calc_xqCF_rotACF
 
    ! 1. Use the xdr interface
    use traj, only: trajfile
    use topol, only: topfile
    use fvector
    use filehandle
    use variables
    use omp_lib
 
    implicit none

    ! 2. Declare a variable of type xtcfile
    type(trajfile) :: trajin
    type(topfile) :: sys

!  --------------------------------------------------------------------------
  integer            :: nsize=10000                 !# of bins 
  integer            :: i,j,k,l,icnt,idx,ixtc,idxmol,idxatom,idxmoltype,idxatomtype           !counters
  integer            :: isize,ifile                  !counters
  integer            :: ctime_0, ctime_1, ctime_rate, ctime_max          !counters_for_timer
  integer            :: narg, cptArg       !#of arg & counter of arg
  integer            :: nframe
  integer            :: nsysatoms, nmolcurr
  integer            :: idframe,ndframe,idt
  integer,dimension(3)   :: idxs
  integer time_array_0(8), time_array_1(8), time_array_t(8)
  real start_time, finish_time
  integer            :: attype1, attype2, maxridx,minridx,natomtype
  integer            :: nattypepairs,ngrpairs,itype2inpair,ixtcGr,iline,ISTAT
  integer            :: numthread,threadid
  logical            :: bmatch


  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: dt,dt0,t1,t2,tstart,tend,halfL, tempr
  real*8             :: mass, molmass, dr, rclustcat, rclustani,charge,molcharge
  real*8             :: murotsq, mutranssq, murottrans, volume
  real*8             :: rotsq,transsq,rottrans,rotval,transval
  real*8             :: debye, eps0, kb, enm, multiple, qelec,mult2,multiConduct
  real*8             :: volavg, muavg
  real*8             :: invvol, rho, sqbinpre, sqbinpre2, sqval, sqval2, form_mol_sum, form_mol_sum2
  real*8             :: pirSinPir, ival, kval, rval
  real*8,dimension(3)    :: sysbox,dist,box,mutot
  real*8,dimension(3)    :: xpos,compos,dxpos   ! temporary position vector for x, center of mass(com)
  real*8,dimension(3)    :: murot, mutrans
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:) :: grsum     !matrix for sum of g(r)
  real*8, allocatable,dimension(:) :: rmat     !matrix for r values in g(r)
  integer, allocatable,dimension(:) :: molTypeIdxs     ! array for molecule type idx of each moleclue ( molecule idx )
  integer, allocatable,dimension(:) :: nmolMolType    ! array for number of molecules in each moltype ( moltype idx )
  integer, allocatable,dimension(:) :: atomTypeIdxs    ! array for atomtype idx of each atom in system ( atomtype idx )
  integer, allocatable,dimension(:,:) :: grPairIdxs    ! array for g(r) pair idx of each atom pairs in system ( atom_i idx, atom_j idx ) returns 0 if no matching g(r) pair is defined
  integer, allocatable,dimension(:,:) :: npartInPair   ! array for number of particles of each atom type in g(r) pair list (gr_pair_idx, atomtype_idx)
  integer, allocatable,dimension(:) :: npartInAtType   ! array for number of particles of each atom type in the system topology (atomtype_idx)

  ! matrices for time stamps
  real*8, allocatable,dimension(:) :: timestamp, timestep     !matrix for timestamps

  ! array of count of frame pairs having time difference dt
  integer,allocatable,dimension(:)   :: nDiffTime   !  (time difference dt in unit of time step)
  ! array of dt in unit of dt0
  integer,allocatable,dimension(:,:)   :: idtmat   !  (ixtc, jxtc)

  ! for conductivity calculation
  ! matrices
  real*8, allocatable,dimension(:,:,:) :: xqcomTraj     ! matrix for q*xyz_com of each molecule at each timestep (frame idx, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:,:) :: xqAtomsTraj     ! matrix for q*xyz of each atom at each timestep (frame idx, atom idx, q*xyz)
  real*8, allocatable,dimension(:) :: xqAtomsCFTime     ! matrix for conduct correlation function of each atom at each dt (time difference dt in unit of time step, atom idx)
  real*8, allocatable,dimension(:) :: xqcomCFTime     ! matrix for conduct correlation function of each molecule at each dt (time difference dt in unit of time step, mol idx, q*xyz)
  real*8, allocatable,dimension(:,:) :: xqcomDiff     ! matrix for difference of xq_xyz of each molecules at each time step combination (ixtc,jxtc,idxmol,3)
  real*8,dimension(3)    :: xqdiff1,xqdiff2


  ! for MSD calculation
  real*8                    :: msd
  real*8,dimension(3)       :: xdiff
  real*8, allocatable,dimension(:,:,:)  :: msdTime      ! translational mean square displacement of each molecule types at each time difference (time difference dt in unit of time step, x;y;z;tot, moltype idx i)

  ! for rotational ACF
  ! matrices for unit normal vector of molecules
  ! c1, c21, c22 forms the plane for imidazolium ring, b or p, f1, f2 forms the plane for bf4 or pf6
  real*8                    :: rotacf
  real*8, dimension(3,3)    :: xPlaneAtom   ! temporary position vector for the three atoms that compose the molecular plane
  real*8, allocatable,dimension(:,:,:) :: unitNormMolTraj     ! unit normal vector of each molecules at each timestep (frame idx, mol idx, vec_xyz)
  real*8, allocatable,dimension(:,:) :: unitNormMolt1, unitNormMolt2     ! unit normal vector of each molecules at each timestep (frame idx, mol idx, vec_xyz)
  !real*8, allocatable,dimension(:,:,:) :: rotacf     ! rotational ACF <ui(t2).ui(t1)> of each molecules at each time step combination (ixtc, jxtc, idxmol)
  real*8, allocatable,dimension(:,:) :: rotACFTime,rotACFTimeP2     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )
  real*8, allocatable,dimension(:) :: rotACFt0     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )

  character(len=30)          :: strfmt
  character(len=256)         :: strout,strtitle
  character(len=3000)         :: str, line
  character(len=256),allocatable,dimension(:) :: aStrArgs
! ----------------------------------------------------------------

  murotsq = 0
  mutranssq = 0
  murottrans = 0
  muavg = 0
  nsysatoms = 0
  volavg = 0
  volume = 0
  nskip = 100      ! If not set in the param file
  nskipgr = 100      ! If not set in the param file
  nbins = 200      ! If not set in the param file
  sqRmax = 0      ! If not set in the param file

  write(*,*) 'get omp num threads'
  !$OMP PARALLEL
  numthread = omp_get_num_threads()
  !$OMP END PARALLEL
  write(6,*) 'Using ',numthread,' number of threads'

  debye = 3.33564E-30
  eps0 = 8.85418781762E-12
  kb = 1.380648813E-23
  enm = 0.020819434
  qelec = 1.60217656535E-19

  multiple = debye**2/(enm**2 * eps0 * kb)
  mult2 = qelec**2 * 1e-18/(eps0*kb)
  multiConduct = qelec**2 / kb

  write(6,*) multiple, mult2,multiConduct

  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
    !loop across options
    allocate(aStrArgs(narg))
    do cptArg=1,narg
      call get_command_argument(cptArg,str)
      aStrArgs(cptArg) = str
    enddo
    !assign in and out files
    call getargs(aStrArgs)
  else
    write(6,*) 'no arguments are given'
    write(6,*) 'usage : calc_cond -f param.dat -o outfile.dat'
    stop
  endif


  ! Read param file to locate strXtcfiles
  call readparam

  ! 3. Read the topology file and initialize atomic mass data for system
  write(6,'(A,A)') 'reading topoly file ',trim(strTopFile)
  write(6,*) ''
  call sys % init(strTopFile)
  call readaff(sys % atomtype)
  natomtype = size(sys % atomtype)

  ! Initialize matrices according to the topology
  nmolsys = sys % numsysmol
  nmoltype = sys % nummoltype
  nsysatoms = sys % numsysatom
  ngrpairs = size(sys % pairlist(:,1))
  allocate(molTypeIdxs(nmolsys))
  allocate(atomTypeIdxs(nsysatoms))
  allocate(nmolMolType(nmoltype))
  allocate(npartInAtType(natomtype))
  idxmol = 0
  idxatom = 0
  npartInAtType = 0
  do i=1, sys % nummoltype
    nmolMolType(i) = sys % moltype(i) % nummol
    do j=1, nmolMolType(i)
      idxmol = idxmol + 1
      molTypeIdxs(idxmol) = i
      do k=1, sys % moltype(i) % numatom
        idxatom = idxatom + 1
        idx = sys % moltype(i) % atomtypeidx(k)
        atomTypeIdxs(idxatom) = idx
        npartInAtType(idx) = npartInAtType(idx) + 1
      enddo
    enddo
  enddo

  allocate(npartInPair(ngrpairs,2))
  allocate(grPairIdxs(nsysatoms,nsysatoms))
  grPairIdxs = 0
  npartInPair = 0
  form_mol = 0
  form_mol2 = 0
  !$OMP PARALLEL DO private (i,k)
  do k = 1,MAX_Q_FORM
    do i=1,natomtype
      form_mol(k) = form_mol(k) + npartInAtType(i) * atomtype_form(i,k)
      form_mol2(k) = form_mol2(k) + npartInAtType(i) * atomtype_form(i,k) ** 2
    enddo
  enddo
  !$OMP END PARALLEL DO

  do k=1,ngrpairs
    do i=1, nsysatoms
      if ((atomTypeIdxs(i) == sys % pairlist(k,1)) .or. &
          (atomTypeIdxs(i) == sys % pairlist(k,2))) then
        if (atomTypeIdxs(i) == sys % pairlist(k,1)) then
            npartInPair(k,1) = npartInPair(k,1) + 1
            if (sys % pairlist(k,1) == sys % pairlist(k,2)) npartInPair(k,2) = npartInPair(k,2) + 1
            itype2inpair = 2
        else
            npartInPair(k,2) = npartInPair(k,2) + 1
            itype2inpair = 1
        endif
        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (grPairIdxs, npartInPair)
        !$OMP DO
        do j=i+1,nsysatoms
          if (atomTypeIdxs(j) == sys % pairlist(k,itype2inpair)) then
            if (itype2inpair == 2) then
              grPairIdxs(i,j) = k
              if (sys % pairlist(k,1) == sys % pairlist(k,2)) grPairIdxs(j,i) = k
            else
              grPairIdxs(j,i) = k
            endif
          endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
      endif
    enddo
    write(6,*) "pair num : ",k, npartInPair(k,1) , npartInPair(k,2)
  enddo

  ! initialize all matrices for first trajectory 
  nattypepairs = natomtype**2
  allocate(grsum(nbins,nattypepairs+1))
  allocate(sqPart(max_q_form,nattypepairs))
  allocate(rmat(nbins))
 

  ! 3. Initialize it with the name of trajin file you want to read in.
  ixtc=1
  ifile=1
    strInFile = trim(strXtcfiles(ifile))
    write(6,'(A,A)') 'reading trajectory file ',trim(strInFile)
    write(6,*) ''
    call trajin % init(strInFile)

    ! 4. Read in each configuration. Everything is stored in the xtcfile
    ! type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for
    !    more details.
    !    You can save the positions in the loop for your calculations in
    !    another array, or 
    !    do your calculations after each read.
 
    call trajin % read

    ! initialize the sin(Qr)/Qr array using the initial box size
    call calcQr(minval((/ (trajin % box(i,i), i=1,3) /)))
    dr_bins = minval((/ (trajin % box(i,i), i=1,3) /))/2/nbins
    minridx = nbins
    if (sqRmax .eq. 0) then
        sqRmax = minridx * dr_bins
    endif

    ! check number of atoms in system
    if(nsysatoms .ne. trajin % NATOMS) then
      write(6,*) 'number of atoms in ',ifile,'th trajectory file does not match other trajectories'
      write(6,*) 'nSysAtoms = ',nsysatoms,'  trajin % NATOMS = ', trajin % NATOMS
      stop
    endif

    OPEN (UNIT=6,FORM='FORMATTED',CARRIAGECONTROL='FORTRAN')
        ! Just an example to show what was read in
        write(6,'(a,f12.6,a,i0)') " Time (ps): ", trajin % time, "  Step: ", trajin % STEP
        write(6,'(a,f12.6,a,i0)') " Precision: ", trajin % prec, "  No. Atoms: ", trajin % NATOMS
        ! This is the same order as found in the GRO format fyi
        write(6,'(9f9.5)')  trajin % box(1,1), trajin % box(2,2), trajin % box(3,3), &
                            trajin % box(1,2), trajin % box(1,3), & 
                            trajin % box(2,1), trajin % box(2,3), &
                            trajin % box(3,1), trajin % box(3,2)
        write(6,*) '' 
        volume = trajin % box(1,1) * trajin % box(2,2) * trajin % box(3,3)
!        write(6,*) 'volume', volume, trajin % box(1,1), trajin % box(2,2), trajin % box(3,3)


        ! check number of atoms in system
        if(nsysatoms .ne. idxatom) then
          write(6,*) 'number of atoms in ',ixtc,'th trajectory does not match other frames'
          write(6,*) 'nSysAtoms = ',nsysatoms,'  idxatom = ', idxatom
          stop
        endif

        !call date_and_time(values=time_array_t)
        !write(6,*) 'time : ', time_array_t
        write(6,100) ixtc,'th frame has finished  ' 
 100    FORMAT('+', I7,A)  

    ! 5. Close the file
    call trajin % close

    volavg = volume        ! average volume in nm^3
    invvol = 1.d0 / volume        ! inverse volume in nm^-3
    rho = nsysatoms * invvol      ! particle density in nm^-3
    volavg = volavg * 1e-27       ! convert the unit to m^3
    multiple = multiple / (temperature * volavg * 3)
    multiConduct = qelec**2 / (6. * kb * temperature * volavg)
    write(6,*) 'temp mult vol', multiConduct, temperature, volavg

    write(6,*) 'start analyzing the data' 

    write(6,*) 'matrix size has adjusted' 
    do k=1,ngrpairs
      write(6,*) "pair num : ",k, npartInPair(k,:)
    enddo

    open(17,file=strGrFile,status='old')
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nattypepairs+1)
    ISTAT = 0
    iline = 0
    do while (ISTAT .eq. 0)
      read(17, '(A)',IOSTAT=ISTAT) line
      if(ISTAT .eq. 0) then
        iline = iline + 1
        read(line,strfmt) rmat(iline), grsum(iline,nattypepairs+1),grsum(iline,1:nattypepairs)
      endif
    enddo
    close(17)
    minridx = iline
    write(6,*) 'gr file is loaded'
    if (minridx*dr_bins .gt. sqRmax) then 
        minridx = sqRmax/dr_bins
    endif
    write(*,*) minridx, 'minridx', sqRmax, dr_bins

    ! average the gr matrice and print the result
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nattypepairs+1)
    sq = 0
    sq2 = 0
    sqlorch = 0
    sqlorch2 = 0
    sqtest = 0
    rmax = dr_bins * minridx
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (sq,sq2,sqPart,sqlorch,sqlorch2,atomtype_form,grsum)
    !$OMP DO
    do k=1, MAX_Q_FORM
      kval = q_grid_form(k)
      do j=1, nattypepairs
        attype1 = (j-1)/natomtype +1
        attype2 = mod((j-1),natomtype)+1
        sqval = invvol * npartInAtType(attype1) * npartInAtType(attype2) * atomtype_form(attype1, k) * atomtype_form(attype2, k)
        if (k==1) then
            write(6,*) attype1, attype2, atomtype_form(attype1,k),atomtype_form(attype2,k)
!            write(6,*) sqval,invvol
        endif
        do i=1, minridx-1
          ival = i + 0.5d0
          rval = dr_bins * ival
          pirSinPir = dsin(pi * rval / rmax ) / (pi * rval / rmax)
          sqbinpre = 4.d0 * pi * rval * dsin( kval * rval) / kval
          sqbinpre2 = sqbinpre * pirSinPir
!          if(k.eq.1 .and. j .eq.1) write(*,*) rval, pirSinPir, sqbinpre
          sqPart(k,j) = sqPart(k,j) + dr_bins * sqval * sqbinpre * (grsum(i,j) - 1.d0)
!          if( (attype1 .eq. 1 .and. attype2 .eq. 9) .or. (attype1 .eq. 9 .and. attype2 .eq. 1))
          sq(k) = sq(k) + dr_bins * sqval * sqbinpre * (grsum(i,j) - 1.d0)
          sqlorch(k) = sqlorch(k) + dr_bins * sqval * sqbinpre2 * (grsum(i,j) - 1.d0)
        enddo
        !sq(k) = sq(k) - sqval * 4.d0 * pi * ( dsin(kval * rmax) - kval * rmax * dcos(kval * rmax) ) / kval**3
        !sqlorch(k) = sqlorch(k) - sqval * 4.d0 * pi * rmax**2 * dsin(kval * rmax) / (pi**2 * kval - kval**3 * rmax**2)
      enddo
      sq2(k) = sq(k) / form_mol2(k)
      sq(k) = sq(k) * nsysatoms / form_mol(k)**2
      sqPart(k,:) = sqPart(k,:) * nsysatoms / form_mol(k)**2

      sqlorch2(k) = sqlorch(k) / form_mol2(k)
      sqlorch(k) = sqlorch(k) * nsysatoms / form_mol(k)**2

      !sqval = sqval / form_mol(k)**2
      !write(17,strfmt) q_grid_form(k), sqval,sqval2*rho
      !write(6,strfmt) q_grid_form(k), sqval,sqval2*rho
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

!    do i=1, minridx-1
!        write(*,*) dr_bins*i, grsum(i,nattypepairs)
!    enddo
    write(strfmt,'("( F12.6, ",I0,"ES15.7 )")') (nattypepairs+4)
    open(17,file=strSqFile2)
    write(17,*) "# k       S(k)       S(k)_form^2       S(k)_Lorch          S(k)_Lorch_form^2"
    write(6,*) "# k       S(k)       S(k)_form^2       S(k)_Lorch          S(k)_Lorch_form^2"
    do k=1, MAX_Q_FORM
      write(17,strfmt) q_grid_form(k), sq(k), sq2(k), sqlorch(k), sqlorch2(k),sqPart(k,:)
!      write(6,*) form_mol(k), form_mol2(k)
      if ( k .eq. 1) write(6,'(F12.3,4ES15.7)') q_grid_form(k), sq(k), sq2(k), sqlorch(k), sqlorch2(k)
    enddo
    close(17)

end program calc_xqCF_rotACF
