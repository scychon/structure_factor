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
  integer            :: nattypepairs,ngrpairs,itype2inpair,ixtcGr
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
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,mutraj,mutrajcat,mutrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:,:) :: comUnwrapTraj   ! unwrapped trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: comMolt1,comMolt2,comMolDiff   ! array for com of each molecules at specific time (mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: boxtraj,tempbox      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: mutrajrot, mutrajtrans      !matrix for murot,mutran
  real*8, allocatable,dimension(:,:) :: sqtraj     !matrix for S(q)
  real*8, allocatable,dimension(:,:,:) :: grframe     !matrix for g(r)
  real*8, allocatable,dimension(:,:) :: grsum,grframesum     !matrix for sum of g(r)
  real*8, allocatable,dimension(:) :: temptime     !matrix for timestamps
  real*8, allocatable,dimension(:) :: Sc     !matrix for q*xyz for each timestep
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
  character(len=256)         :: str, line
  character(len=256),allocatable,dimension(:) :: aStrArgs
! ----------------------------------------------------------------

  murotsq = 0
  mutranssq = 0
  murottrans = 0
  muavg = 0
  nframe = 0
  nsysatoms = 0
  volavg = 0
  volume = 0
  nskip = 100      ! If not set in the param file
  nskipgr = 100      ! If not set in the param file
  nbins = 200      ! If not set in the param file

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
  form_mol_sum = sum(form_mol)
  form_mol_sum2 = sum(form_mol2)

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
  allocate(comMolt1(nmolsys,3))
  !allocate(comMolt2(nmolsys,3))
  allocate(comMolDiff(nmolsys,3))
  allocate(comtraj(nsize,nmolsys,3))
  allocate(comUnwrapTraj(nsize,nmolsys,3))
  allocate(mutraj(nsize,nmolsys,3))
  allocate(unitNormMolTraj(nsize,nmolsys,3))
!  allocate(xqAtomsTraj(nsize,nsysatoms,3))
  allocate(boxtraj(nsize,3))
  !allocate(sqtraj(nsize/nskip,MAX_Q_FORM))
  nattypepairs = natomtype**2
  allocate(grframe(numthread,nbins,nattypepairs+1))
  allocate(grframesum(nbins,nattypepairs+1))
  !allocate(grframe(numthread,nbins,ngrpairs+2))
  !allocate(grframesum(nbins,ngrpairs+2))
  grframesum = 0.d0
  allocate(grsum(nbins,nattypepairs+1))
  !allocate(grsum(nbins,ngrpairs+2))
  grsum = 0.d0
  allocate(timestamp(nsize))
  allocate(timestep(nsize))
  write(6,*) 'comMolt1 & comMolDiff allocated' 
 

  ! 3. Initialize it with the name of trajin file you want to read in.
  ixtc=0
  do ifile=1,nxtcfile
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

    if(ixtc .eq. 0) then
!      nsysatoms = trajin % NATOMS
      ! initialize the sin(Qr)/Qr array using the initial box size
      call calcQr(minval((/ (trajin % box(i,i), i=1,3) /)))
      dr_bins = minval((/ (trajin % box(i,i), i=1,3) /))/2/nbins
      minridx = nbins
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
    do while ( trajin % STAT == 0 )
        ixtc = ixtc + 1
        volume = trajin % box(1,1) * trajin % box(2,2) * trajin % box(3,3)
        volavg = volavg + volume
!        write(6,*) 'volume', volume, trajin % box(1,1), trajin % box(2,2), trajin % box(3,3)

        ! check the size of matrices 
        isize = size(comtraj(:,1,1))
        if(isize .lt. ixtc) then
          write(6,*) 'trajectory is longer than ', nsize
            call date_and_time(values=time_array_t)
            write(6,*) 'time : ', time_array_t
          call expand3D(comtraj,nsize,0,0)
          call expand3D(comUnwrapTraj,nsize,0,0)
          call expand3D(mutraj,nsize,0,0)
          call expand3D(unitNormMolTraj,nsize,0,0)
        !  call expand3D(xqAtomsTraj,nsize,0,0)
          call expand2D(boxtraj,nsize,0)
          !call expand2D(sqtraj,nsize/nskip,0)
            write(6,*) size(timestamp(:))
          call expand1D(timestamp,nsize)
          call expand1D(timestep,nsize)
          nsize = nsize*2
            call date_and_time(values=time_array_t)
            write(6,*) 'time : ', time_array_t
            write(6,*) ''
            write(6,*) size(timestamp(:))
        endif

        ! record current time
        timestamp(ixtc) = trajin % time
        timestep(ixtc) = trajin % STEP


        do k=1,3
          boxtraj(ixtc,k) = trajin % box(k,k)
        enddo

        ! Compute and store the trajectory of the center of mass postion,
        ! the rotational principle vector, and the dipole moment vector
        idxatom = 0
        idxmol = 0
        do i=1, sys % nummoltype
          !$OMP PARALLEL &
          !$OMP   DEFAULT (FIRSTPRIVATE) &
          !$OMP   SHARED (unitNormMolTraj, comtraj, mutraj)
          !$OMP DO
          do j=1, sys % moltype(i) % nummol
            compos = 0
            mutot = 0
            idxs = sys % moltype(i) % molPlaneAtomIdxs

            do k=1, sys % moltype(i) % numatom
              idx = idxatom + (j-1) * sys % moltype(i) % numatom + k
              mass = sys % moltype(i) % atommass(k)
              charge = sys % moltype(i) % atomcharge(k)
              xpos(:) = trajin % pos(:,idx)
              compos(:) = compos(:) + xpos(:)*mass
            !  xqAtomsTraj(ixtc,idx,:) = xpos(:)*charge
              mutot(:) = mutot(:) + xpos(:)*charge
              if (k .eq. idxs(1) ) then
                  xPlaneAtom(1,:) = xpos(:)
              elseif (k .eq. idxs(2) ) then
                  xPlaneAtom(2,:) = xpos(:)
              elseif (k .eq. idxs(3) ) then
                  xPlaneAtom(3,:) = xpos(:)
              end if
            enddo
            unitNormMolTraj(ixtc,idxmol+j,:) = vecUnitNorm( xPlaneAtom(1,:)-xPlaneAtom(2,:), xPlaneAtom(1,:)-xPlaneAtom(3,:) )

            molmass = sys % moltype(i) % molarmass
            compos(:) = compos(:)/molmass
            comtraj(ixtc,idxmol+j,:) = compos(:)
            mutraj(ixtc,idxmol+j,:) = mutot(:)
!            write(6,*) comtraj(ixtc,idxmol+j,:)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
          idxmol = idxmol + sys % moltype(i) % nummol
          idxatom = idxatom + sys % moltype(i) % nummol * sys % moltype(i) % numatom
        enddo 

        ! unwrap the trajectory
        if(ixtc .eq. 1) then
          comUnwrapTraj(ixtc,:,:) = comtraj(ixtc,:,:)
        else if (ixtc .gt. 1) then
          !$OMP PARALLEL &
          !$OMP   DEFAULT (FIRSTPRIVATE) &
          !$OMP   SHARED (comtraj, comMolDiff, comUnwrapTraj)
          !$OMP DO
          do i=1, nmolsys
            comMolDiff(i,:) = comtraj(ixtc,i,:) - comtraj(ixtc-1,i,:)
            comMolDiff(i,1) = comMolDiff(i,1) - NINT(comMolDiff(i,1)/boxtraj(ixtc,1))*boxtraj(ixtc,1)
            comMolDiff(i,2) = comMolDiff(i,2) - NINT(comMolDiff(i,2)/boxtraj(ixtc,2))*boxtraj(ixtc,2)
            comMolDiff(i,3) = comMolDiff(i,3) - NINT(comMolDiff(i,3)/boxtraj(ixtc,3))*boxtraj(ixtc,3)
            comUnwrapTraj(ixtc,i,:) = comUnwrapTraj(ixtc-1,i,:)+comMolDiff(i,:)
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        endif 

        ! check number of atoms in system
        if(nsysatoms .ne. idxatom) then
          write(6,*) 'number of atoms in ',ixtc,'th trajectory does not match other frames'
          write(6,*) 'nSysAtoms = ',nsysatoms,'  idxatom = ', idxatom
          stop
        endif

        !call date_and_time(values=time_array_t)
        !write(6,*) 'time : ', time_array_t

        halfL = minval(boxtraj(ixtc,:))/2
        maxridx = int(halfL/dr_bins)
        if(maxridx .lt. minridx) then
            minridx = maxridx
        endif
        ! Compute and store the g(r)_ij and the S(q)_ij for each snapshot
        if ( mod(ixtc,nskipgr) == 0) then
            !call date_and_time(values=time_array_t)
            !write(6,*) 'time : ', time_array_t
            ixtcGr = ixtc/nskipgr
            grframe = 0.d0
            !sqtraj(ixtcGr,:) = 0
            !write(6,*) 'Starting to generate qr and sq matrices'
            do i=1, nsysatoms-1
                attype1 = atomTypeIdxs(i)
                if (attype1 .le. 0) cycle
              !$OMP PARALLEL &
              !$OMP   DEFAULT (FIRSTPRIVATE) &
              !$OMP   SHARED (grframe,grPairIdxs,trajin)
                threadid = omp_get_thread_num()+1
                !write(6,*) threadid
              !$OMP DO
                do j=i+1, nsysatoms
                    attype2 = atomTypeIdxs(j)
                    if (attype2 .le. 0) cycle
                    dxpos(:) = trajin % pos(:,i) - trajin % pos(:,j)
                    dr = getdr(dxpos, boxtraj(ixtc,:))
                    if(dr .lt. halfL) then
                        idx = int(dr/dr_bins)
                        if (idx .le. nbins) then
                            grframe(threadid,idx,nattypepairs+1) = grframe(threadid,idx,nattypepairs+1) + 2
                            !grframe(threadid,idx,nattypepairs+2) = grframe(threadid,idx,nattypepairs+2) + 2*atomtype_form(attype1,1)*atomtype_form(attype2,1)
                            !write(6,*) 'added ', grframe(threadid,idx,nattypepairs+1)
                            idxatomtype = (attype1-1) * natomtype + attype2
                            grframe(threadid,idx,idxatomtype) = grframe(threadid,idx,idxatomtype) + 1
                            idxatomtype = (attype2-1) * natomtype + attype1
                            grframe(threadid,idx,idxatomtype) = grframe(threadid,idx,idxatomtype) + 1
                            !grframe(threadid,idx,grPairIdxs(i,j)) = grframe(threadid,idx,grPairIdxs(i,j)) + 1
                            !if(grPairIdxs(i,j) .gt. 0) then
                            !    grframe(threadid,idx,grPairIdxs(i,j)) = grframe(threadid,idx,grPairIdxs(i,j)) + 1
                            !endif
                            !if(grPairIdxs(j,i) .gt. 0) then
                            !    grframe(threadid,idx,grPairIdxs(j,i)) = grframe(threadid,idx,grPairIdxs(j,i)) + 1
                            !endif
                            !    do k=1,max_q_form
                            !        sqtraj(ixtcGr,k) = sqtraj(ixtcGr,k) + atomtype_form(attype1,k)*atomtype_form(attype2,k)*sinqr_over_qr(k,idx)*2/nsysatoms
                            !    enddo
                            !endif
                        endif
                    endif
                enddo
              !$OMP END DO
              !$OMP END PARALLEL
            enddo


            grframesum = sum(grframe,dim=1)
            !$OMP PARALLEL &
            !$OMP   DEFAULT (FIRSTPRIVATE) &
            !$OMP   SHARED (grframe,grframesum,grsum,npartInPair)
            !$OMP DO
            do i=1,minridx
                grsum(i,nattypepairs+1) = grsum(i,nattypepairs+1) + grframesum(i, nattypepairs+1) / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * float(nsysatoms * nsysatoms)/volume )
                !grsum(i,nattypepairs+2) = grsum(i,nattypepairs+2) + grframesum(i, nattypepairs+2) / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * float(nsysatoms) )
                do j=1,nattypepairs
                    grsum(i,j) = grsum(i,j) + grframesum(i,j) / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * float(npartInAtType((j-1)/natomtype + 1) * npartInAtType(mod((j-1),natomtype)+1))/volume )
                enddo
                !grsum(i,1:nattypepairs) = grsum(i,1:nattypepairs) + grframesum(i,1:nattypepairs) / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * float(npartInPair(:,1) * npartInPair(:,2))/volume )
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !call date_and_time(values=time_array_t)
            !write(6,*) 'time : ', time_array_t
            !write(6,*) ''
        endif


        !call date_and_time(values=time_array_t)
        !write(6,*) 'time : ', time_array_t
        write(6,100) ixtc,'th frame has finished  ' 
 100    FORMAT('+', I7,A)  

        call trajin % read
 
    end do

    ! 5. Close the file
    call trajin % close
    end do

    nxtc = ixtc
    nframe = nframe + nxtc
    volavg = volavg / nxtc        ! average volume in nm^3
    invvol = 1.d0 / volavg        ! inverse volume in nm^-3
    rho = nsysatoms * invvol      ! particle density in nm^-3
    volavg = volavg * 1e-27       ! convert the unit to m^3
    multiple = multiple / (temperature * volavg * 3)
    multiConduct = qelec**2 / (6. * kb * temperature * volavg)
    write(6,*) 'temp mult vol', multiConduct, temperature, volavg

    write(6,*) 'start analyzing the data' 

    ! check the size of trajectory matrix
    isize = size(comtraj(:,1,1))
    if(isize .gt. nxtc) then
      call shrink3D(comtraj,nxtc,nmolsys,3)
      call shrink3D(comUnwrapTraj,nxtc,nmolsys,3)
      call shrink3D(mutraj,nxtc,nmolsys,3)
      call shrink3D(unitNormMolTraj,nxtc,nmolsys,3)
    !  call shrink3D(xqAtomsTraj,nxtc,nsysatoms,3)
      call shrink2D(boxtraj,nxtc,3)
      !call shrink2D(sqtraj,ixtcGr,max_q_form)
      call shrink1D(timestamp,nxtc)
      call shrink1D(timestep,nxtc)
    endif

    write(6,*) 'matrix size has adjusted' 
    do k=1,ngrpairs
      write(6,*) "pair num : ",k, npartInPair(k,:)
    enddo

    ! average the gr matrice and print the result
    grsum = grsum / float(ixtcGr)
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nattypepairs+1)
    sq = 0
    sq2 = 0
    sqlorch = 0
    sqlorch2 = 0
    sqtest = 0
    rmax = dr_bins * minridx
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (sq,sq2,sqlorch,sqlorch2,atomtype_form,grsum)
    !$OMP DO
    do k=1, MAX_Q_FORM
      kval = q_grid_form(k)
      do j=1, nattypepairs
        attype1 = (j-1)/natomtype +1
        attype2 = mod((j-1),natomtype)+1
        sqval = invvol * npartInAtType(attype1) * npartInAtType(attype2) * atomtype_form(attype1, k) * atomtype_form(attype2, k)
        if (k==1) then
            write(6,*) attype1, attype2, atomtype_form(attype1,k),atomtype_form(attype2,k)
        endif
        do i=1, minridx-1
          ival = i + 0.5d0
          rval = dr_bins * ival
          pirSinPir = dsin(pi * rval / rmax ) / (pi * rval / rmax)
          sqbinpre = 4.d0 * pi * rval * dsin( kval * rval) / kval
          sqbinpre2 = sqbinpre * pirSinPir
          sq(k) = sq(k) + dr_bins * sqval * sqbinpre * (grsum(i,j) - 1.d0)
          sqlorch(k) = sqlorch(k) + dr_bins * sqval * sqbinpre2 * (grsum(i,j) - 1.d0)
        enddo
        !sq(k) = sq(k) - sqval * 4.d0 * pi * ( dsin(kval * rmax) - kval * rmax * dcos(kval * rmax) ) / kval**3
        !sqlorch(k) = sqlorch(k) - sqval * 4.d0 * pi * rmax**2 * dsin(kval * rmax) / (pi**2 * kval - kval**3 * rmax**2)
      enddo
      sq2(k) = sq(k) / form_mol2(k)
      sq(k) = sq(k) * nsysatoms / form_mol(k)**2

      sqlorch2(k) = sqlorch(k) / form_mol2(k)
      sqlorch(k) = sqlorch(k) * nsysatoms / form_mol(k)**2

      !sqval = sqval / form_mol(k)**2
      !write(17,strfmt) q_grid_form(k), sqval,sqval2*rho
      !write(6,strfmt) q_grid_form(k), sqval,sqval2*rho
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    open(17,file=strSqFile)
    write(17,*) "# k       S(k)       S(k)_form^2       S(k)_Lorch          S(k)_Lorch_form^2"
    write(6,*) "# k       S(k)       S(k)_form^2       S(k)_Lorch          S(k)_Lorch_form^2"
    do k=1, MAX_Q_FORM
      write(17,'(F12.3,4ES15.7)') q_grid_form(k), sq(k), sq2(k), sqlorch(k), sqlorch2(k)
      write(6,'(F12.3,4ES15.7)') q_grid_form(k), sq(k), sq2(k), sqlorch(k), sqlorch2(k)
    enddo
    close(17)
    open(17,file=strGrFile)
    write(strfmt,'("( F12.6, ",I0,"ES15.7 )")') (nattypepairs+1)
    do i=1, minridx-1
      write(17,strfmt) dr_bins*float(i), grsum(i,nattypepairs+1),grsum(i,1:nattypepairs)
    enddo
    close(17)
    write(6,*) 'first gr saved'

    ! calculate size of time difference matrices
    ! assume the first two steps have minimum dt
    ! and the trajectories are time-ordered
    dt0 = timestamp(2)-timestamp(1)
    dt = timestamp(nxtc)-timestamp(1)
    ndframe = int(dt/dt0 + 0.00001)

    write(6,*) 'ndframe is set' 
    ! generate and initialize time difference matrices
!    allocate(idtmat(nxtc-1,nxtc))
    write(6,*) 'idtmat allocated' 
    allocate(unitNormMolt1(nmolsys,3))
    allocate(unitNormMolt2(nmolsys,3))
    write(6,*) 'unitNormMolt1 & unitNormMolDiff allocated' 

!    allocate(xqAtomsCFTime(ndframe))
    allocate(xqcomCFTime(ndframe))
    write(6,*) 'xqcomCFTime allocated' 
    allocate(xqcomDiff(nmolsys,3))
    write(6,*) 'xqcomDiff allocated' 
    allocate(msdTime(ndframe,4,nmoltype+1))
    write(6,*) 'msdTime allocated' 
    allocate(rotACFTime(ndframe,nmoltype+1))
    write(6,*) 'rotACFTime allocated' 
    allocate(rotACFTimeP2(ndframe,nmoltype+1))
    write(6,*) 'rotACFTimeP2 allocated' 
 !   allocate(rotacf(nxtc-1,nxtc,nmolsys))
!    allocate(xqAtomsDiff(nxtc-1,nxtc,nsysatoms,3))
    allocate(xqcomTraj(nxtc,nmolsys,3))
    write(6,*) 'xqcomTraj allocated' 
 !   allocate(nACFTime(ndframe))
    allocate(nDiffTime(ndframe))
    write(6,*) 'nDiffTime allocated' 
    allocate(rotACFt0(nmoltype+1))
!    xqAtomsCFTime = 0
    xqcomCFTime = 0
    xqcomDiff = 0
!    xqAtomsDiff = 0
    msdTime = 0
    rotACFTime = 0
    rotACFTimeP2 = 0
    rotACFt0 = 1
 !   rotacf = 0
    xqcomTraj = 0
 !   nACFTime = 0
    nDiffTime = 0

    write(6,*) 'generating comUnwrapTraj and xqcomTraj matrix' 
    idxmol = 0
    do i=1, sys % nummoltype
      molcharge = sys % moltype(i) % molcharge
      !$OMP PARALLEL &
      !$OMP   DEFAULT (FIRSTPRIVATE) &
      !$OMP   SHARED (comtraj,comUnwrapTraj,xqcomTraj)
      !$OMP DO
      do j=1, sys % moltype(i) % nummol
        !xqcomTraj(:,idxmol+j,:) = comtraj(:,idxmol+j,:)*molcharge
        xqcomTraj(:,idxmol+j,:) = comUnwrapTraj(:,idxmol+j,:)*molcharge
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      idxmol = idxmol + sys % moltype(i) % nummol
    enddo


    write(6,*) 'start generating rotacf and xqdiff matrices'
    write(6,*) 'first frame start'
    
    call date_and_time(values=time_array_0)
    write(6,*) 'Start time : ', time_array_0
    do i=1,nxtc-1,nskip
      t1 = timestamp(i)
      unitNormMolt1(:,:) = unitNormMolTraj(i,:,:)
      !$OMP PARALLEL &
      !$OMP   DEFAULT (FIRSTPRIVATE) &
      !$OMP   SHARED (nDiffTime, rotACFTime, rotACFTimeP2,msdTime, &
      !$OMP          xqcomCFTime, unitNormMolTraj,comtraj,comUnwrapTraj,xqcomTraj)
      !$OMP DO
      do j=i+1,nxtc
        t2 = timestamp(j)
        dt = t2 - t1
        idt = int(dt/dt0 + 0.00001)
        if(idt .le. 0) then
          write(6,*) 'idt is less than 0', idt, t1, t2, dt
          cycle
        endif
        !comMolDiff(:,:) = comtraj(i,:,:)- comtraj(j,:,:)
        comMolDiff(:,:) = comUnwrapTraj(i,:,:)- comUnwrapTraj(j,:,:)
        unitNormMolt2(:,:) = unitNormMolTraj(j,:,:)

     !   nACFTime(idt) = nACFTime(idt) + 1
        nDiffTime(idt) = nDiffTime(idt) + 1

        idxmol = 0
        xqcomDiff(:,:) = xqcomTraj(i,:,:)-xqcomTraj(j,:,:)
        xqdiff1 = 0
        do k=1,nmolsys
          idxmoltype = molTypeIdxs(k)
          nmolcurr = nmolMolType(idxmoltype)
          xdiff(:) = comMolDiff(k,:)
          msd = dot_product(xdiff, xdiff)
          msdTime(idt,1:3,idxmoltype) = msdTime(idt,1:3,idxmoltype) + xdiff(:)/nmolcurr
          msdTime(idt,1:3,nmoltype+1) = msdTime(idt,1:3,nmoltype+1) + xdiff(:)/nmolsys
          msdTime(idt,4,idxmoltype) = msdTime(idt,4,idxmoltype) + msd/nmolcurr
          msdTime(idt,4,nmoltype+1) = msdTime(idt,4,nmoltype+1) + msd/nmolsys
          
          rotacf = dot_product(unitNormMolt1(k,:),unitNormMolt2(k,:))
          rotACFTime(idt,idxmoltype) = rotACFTime(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTime(idt,nmoltype+1) = rotACFTime(idt,nmoltype+1) + rotacf/nmolsys
          rotacf = rotacf*rotacf
          rotACFTimeP2(idt,idxmoltype) = rotACFTimeP2(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTimeP2(idt,nmoltype+1) = rotACFTimeP2(idt,nmoltype+1) + rotacf/nmolsys
          xqdiff1(:) = xqdiff1(:) + xqcomDiff(k,:)
        enddo

            xqcomCFTime(idt) = xqcomCFTime(idt) + dot_product(xqdiff1,xqdiff1)
!        do k=1,nmolsys
!          xqdiff1(:) = xqcomDiff(k,:)
!          do l=1,nmolsys
!            xqdiff2(:) = xqcomDiff(l,:)
!          enddo
!        enddo

      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      write(6,100) i,'th frame has finished  ' 
    enddo

    call date_and_time(values=time_array_1)
    write(6,*) 'End time : ', time_array_1

      start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)
    write(6,*) 'Elapsed time : ', finish_time - start_time, ' s'

    write(6,*) 'rotacf and xqdiff matrices generated'
    write(6,*) 'generating correlation function matrices'


!    open(21,file='testout')
!    do i=1,nmolsys
!      xqdiff1(:) = xqcomTraj(1,i,:)-xqcomTraj(11,i,:)
!      write(21,*) i, xqcomTraj(1,i,:),xqcomTraj(11,i,:),dot_product(xqdiff1,xqdiff1)
!      write(21,*) i,unitNormMolTraj(1,i,:),unitNormMolTraj(2,i,:), dot_product(unitNormMolTraj(1,i,:),unitNormMolTraj(2,i,:))
!    enddo

    write(6,*) 'data analysis done'
    

    write(6,*) 'start generating output files'

    call backupfile(strMSDFile)
    call backupfile(strConductFile)
    call backupfile(strRotACFFile)
    call backupfile(strRotACFP2File)
    open(18,file=strMSDFile)
    open(19,file=strConductFile)
    open(20,file=strRotACFFile)
    open(21,file=strRotACFP2File)

!    xqAtomsCFTime = xqAtomsCFTime*multiConduct
    xqcomCFTime = xqcomCFTime*multiConduct


    write(18,*) '# MSD of each molecular species, last trajectory from : ',strInFile
    write(18,*) '# n_molType : ', nmoltype
    write(18,*) '# x1 : time (ps)'
    idx=1
    do i=1, sys % nummoltype
      write(18,*) '# x',idx+1,' : msd_x for ',i,' th molecule - ',sys % moltype(i) % moltypename
      write(18,*) '# x',idx+2,' : msd_y for ',i,' th molecule - ',sys % moltype(i) % moltypename
      write(18,*) '# x',idx+3,' : msd_z for ',i,' th molecule - ',sys % moltype(i) % moltypename
      write(18,*) '# x',idx+4,' : msd_tot for ',i,' th molecule - ',sys % moltype(i) % moltypename
      idx = idx + 4
    enddo
    !write(18,'(A5,3E15.7)') 'tot  ',rotsq,transsq,rottrans
    !write(18,'(A5,3E15.7)') 'tot  ',rotsq,transsq,rottrans
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nmoltype+1)*4
    ! write (0,1) for auto correlation functions
    write(20,strfmt) 0,rotACFt0
    write(21,strfmt) 0,rotACFt0
    do idt=1,ndframe
     ! icnt = nACFTime(idt)
      icnt = nDiffTime(idt)
!      if(xqAtomsCFTime(idt) .ne. 0) then
!        write(18,'(2E15.7)') dt0*idt,xqAtomsCFTime(idt)/icnt
!      endif
      if(msdTime(idt,4,nmoltype+1) .ne. 0) then
        write(18,strfmt) dt0*idt,msdTime(idt,:,:)/icnt
      endif
      if(xqcomCFTime(idt) .ne. 0) then
        write(19,'(F12.3,ES15.7)') dt0*idt,xqcomCFTime(idt)/icnt
      endif

      if(rotACFTime(idt,nmoltype+1) .ne. 0) then
        write(20,strfmt) dt0*idt,rotACFTime(idt,:)/icnt
        rotACFTimeP2(idt,:) = 1.5d0*(rotACFTimeP2(idt,:)/icnt)-0.5d0
        write(21,strfmt) dt0*idt,rotACFTimeP2(idt,:)
      endif
    enddo


    write(6,*) 'finished writing output files'
    close(6)
    close(18)
    close(19)
    close(20)
    close(21)

end program calc_xqCF_rotACF
