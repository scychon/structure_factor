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
    use read_frame
    use omp_lib
 
    implicit none

    ! 2. Declare a variable of type xtcfile
    type(trajfile) :: trajin
    type(topfile) :: sys

!  --------------------------------------------------------------------------
  integer            :: i,j,k,l,icnt,idx,idxcnt,ixtc,ixtc0,idxmol,idxatom,idxmoltype,idxatomtype           !counters
  integer            :: isize,ifile,istat,icurstat                  !counters
  integer            :: ctime_0, ctime_1, ctime_rate, ctime_max          !counters_for_timer
  integer            :: nframe
  integer            :: nsysatoms, nmolcurr, nsyspart
  integer            :: idframe,ndframe,idt
  integer,dimension(3)   :: idxs
  integer time_array_0(8), time_array_1(8), time_array_t(8)
  real start_time, finish_time
  integer            :: attype1, attype2, maxridx,minridx,natomtype
  integer            :: nattypepairs,itype2inpair,ixtcGr,nqtot
  integer            :: numthread,threadid
  logical            :: bmatch


  real*8             :: delr=0.02d0,rmax=3.0d0,boxsize=0.3d0      !
  real*8             :: dt,dt0,t1,t2,tstart,tend,halfL, tempr
  real*8             :: charge,molcharge
  real*8             :: volume
  real*8             :: multiple, mult2,multiConduct
  real*8             :: volavg, muavg
  real*8             :: invvol, rho, sqbinpre, sqbinpre2, sqval, sqval2
  real*8             :: pirSinPir, ival, kval, dkval
  real*8,dimension(3)    :: sysbox,dist,box,dk, dkavg
  real*8,dimension(3)    :: murot, mutrans
  real*8,dimension(300)  :: rbin
  real*8, allocatable,dimension(:,:,:) :: comtraj,comtrajcat,comtrajani,mutraj,mutrajcat,mutrajani,temp   ! trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:,:) :: comUnwrapTraj   ! unwrapped trajectory of com of each molecules (frame idx, mol idx, xyz)
  real*8, allocatable,dimension(:,:) :: comMolDiff   ! array for com of each molecules at specific time (xyz, mol idx)
  real*8, allocatable,dimension(:,:) :: boxtraj      !matrix for # of water molecules along z direction
  real*8, allocatable,dimension(:,:) :: mutrajrot, mutrajtrans      !matrix for murot,mutran
  real*8, allocatable,dimension(:,:) :: sqtraj     !matrix for S(q)
  real*8, allocatable,dimension(:,:) :: xposframe   !matrix for positions in the frame
  real*8, allocatable,dimension(:,:) :: grsum,grout     !matrix for sum of g(r)
  real*8, allocatable,dimension(:) :: sqsum,sqavg     !matrix for sum of s(q)
  integer, allocatable,dimension(:) :: sqcnt     !matrix for sum of s(q)
  real*8, allocatable,dimension(:) :: temptime     !matrix for timestamps
  real*8, allocatable,dimension(:) :: Sc     !matrix for q*xyz for each timestep
  integer, allocatable,dimension(:) :: molTypeIdxs     ! array for molecule type idx of each moleclue ( molecule idx )
  integer, allocatable,dimension(:) :: nmolMolType    ! array for number of molecules in each moltype ( moltype idx )
  integer, allocatable,dimension(:) :: atomTypeIdxs    ! array for atomtype idx of each atom in system ( atomtype idx )
  integer, allocatable, dimension(:,:) :: atomPairIdxs     ! array for atompairIdxs

  integer, allocatable,dimension(:) :: npartInAtType   ! array for number of particles of each atom type in the system topology (atomtype_idx)
  real*8, allocatable,dimension(:) :: chargeAtType   ! array for number of particles of each atom type in the system topology (atomtype_idx)
  real*8, allocatable,dimension(:)   :: chargeAtoms, chargeMols, chargeMolType   !matrix for charges of each atoms / molecules

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
  real*8, allocatable,dimension(:,:) :: xqcomDiff     ! matrix for difference of xq_xyz of each molecules at each time step combination (3, idxmol)
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
  real*8, allocatable,dimension(:,:,:) :: unitNormMolTraj     ! unit normal vector of each molecules at each timestep (vec_xyz, mol idx, frame idx)
  real*8, allocatable,dimension(:,:) :: unitNormMolt1, unitNormMolt2     ! unit normal vector of each molecules at each timestep (vec_xyz, mol idx)
  !real*8, allocatable,dimension(:,:,:) :: rotacf     ! rotational ACF <ui(t2).ui(t1)> of each molecules at each time step combination (ixtc, jxtc, idxmol)
  real*8, allocatable,dimension(:,:) :: rotACFTime,rotACFTimeP2     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )
  real*8, allocatable,dimension(:) :: rotACFt0     ! rotational ACF <ui(t2).ui(t1)> of each molecule types at each time difference (time difference dt in unit of time step, moltype idx i )

  character(len=30)          :: strfmt
  character(len=256)         :: strFilename
! ----------------------------------------------------------------

  nxtc = 0
  nframe = 0
  volavg = 0
  nskip = 100      ! If not set in the param file
  nskipgr = 100      ! If not set in the param file
  nbins = 200      ! If not set in the param file
  nqgrid = 10

  write(*,*) 'get omp num threads'
  !$OMP PARALLEL
  !$OMP SINGLE
  numthread = omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  write(6,*) 'Using ',numthread,' number of threads'

  multiple = DEBYE**2/(ENM**2 * EPS0 * BOLTZMANN_CONST)
  mult2 = ELEMENTARY_CHARGE**2 * 1e-18/(EPS0*BOLTZMANN_CONST)

  write(6,*) multiple, mult2

  !Check if any arguments are found
  call checkargs

  ! Read param file to locate strXtcfiles
  call readparam

  ! 3. Read the topology file and initialize atomic mass data for system
  call sys % init(strTopFile)

  ! Initialize matrices according to the topology
  nmolsys = sys % numsysmol
  nmoltype = sys % nummoltype
  nsysatoms = sys % numsysatom
  natomtype = size(sys % atomtype)
  allocate(molTypeIdxs(nmolsys))
  allocate(atomTypeIdxs(nsysatoms))
  allocate(nmolMolType(nmoltype))
  allocate(npartInAtType(natomtype))
  allocate(chargeAtType(natomtype))
  allocate(xposframe(3,nsysatoms))
  allocate(chargeAtoms(nsysatoms))
  allocate(chargeMols(nmolsys))
  allocate(chargeMolType(nmoltype))

  idxmol = 0
  idxatom = 0
  npartInAtType = 0
  atomTypeIdxs = 0
  chargeAtType = 0
  chargeAtoms = 0
  chargeMols = 0
  do i=1, sys % nummoltype
    nmolMolType(i) = sys % moltype(i) % nummol
    molcharge = sys % moltype(i) % molcharge
    chargeMolType(i) = molcharge
    do j=1, nmolMolType(i)
      idxmol = idxmol + 1
      molTypeIdxs(idxmol) = i
      chargeMols(idxmol) = molcharge
      do k=1, sys % moltype(i) % numatom
        idxatom = idxatom + 1
        idx = sys % moltype(i) % atomtypeidx(k)
        charge = sys % moltype(i) % atomcharge(k)
        atomTypeIdxs(idxatom) = idx
        chargeAtoms(idxatom) = charge
        if ((chargeAtType(idx) .ne. 0) .and. (chargeAtType(idx) .ne. charge)) then
            write(*,*) 'Error : charge of atomtype ',sys % moltype(i) % atomtype(k),' is not unique!'
        endif
        chargeAtType(idx) = charge
        if( idx .gt. 0) npartInAtType(idx) = npartInAtType(idx) + 1
      enddo
    enddo
  enddo

  !call set_atompairs(atomTypeIdxs, atomPairIdxs)

  write(*,*) 'start qform generation'
  form_mol = 0
  form_mol2 = 0
  !!$OMP PARALLEL DO private (i,k)
  do k = 1,MAX_Q_FORM
    do i=1,natomtype
      form_mol(k) = form_mol(k) + npartInAtType(i) * chargeAtType(i)
      form_mol2(k) = form_mol2(k) + npartInAtType(i) * chargeAtType(i) ** 2
    enddo
  enddo
  !!$OMP END PARALLEL DO

  ! initialize all matrices for first trajectory 
  allocate(comMolDiff(3,nmolsys))
  if(bGrCOM) then
      nattypepairs = nmoltype**2
      nsyspart = nmolsys
      write(6,*) 'maxmoltype',nmoltype,maxval(molTypeIdxs),nattypepairs
  else
      nattypepairs = natomtype**2
      nsyspart = sum(npartInAtType)
      !nsyspart = nsysatoms
  endif
  allocate(grsum(nattypepairs+1,nbins))
  grsum = 0.d0
  allocate(sqsum(3*nqgrid(1)*nqgrid(2)*nqgrid(3)))
  sqsum = 0.d0
  dkavg = 0.d0
  allocate(sqPart(max_q_form,nattypepairs))
  allocate(sqPartDirect(max_q_form,nattypepairs))
  write(6,*) 'comMolDiff allocated' 
 

  ! 3. Initialize it with the name of trajin file you want to read in.
  ixtc=0
  ixtcGr=0
  istat=0
  icurstat = 0
  OPEN (UNIT=6,FORM='FORMATTED')
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
 

    write(*,*) 'nframes : ', nXtcFrames(ifile)
    write(*,*) ''
    if (trajin % NFRAMES == 0) then
        if (nXtcFrames(ifile) == 0) then
            write(6,*) 'number of frames is unknown. Start counting them'
            call trajin % count(strInFile, nframe)
            write(6,*) 'number of frames : ',nframe
            write(*,*) ''
        else
            nframe = nXtcFrames(ifile)
        endif
    elseif (nXtcFrames(ifile) .ne. 0) then
        if (trajin % NFRAMES .ne. nXtcFrames(ifile)) then
            write(6,*) 'number of frames in trajfile ', trajin % NFRAMES, ' does not match the parameter file ', nXtcFrames(ifile), '.'
            nframe = min(trajin % NFRAMES, nXtcFrames(ifile))
            write(6,*) 'The smaller value ', nframe, ' is used here.'
        else
            nframe = trajin % NFRAMES
        endif
    else
        nframe = trajin % NFRAMES
    endif

    if(ixtc .eq. 0) then
        write(6,*) 'Starting allocating initial arrays.'
        allocate(comtraj(3,nmolsys,nframe))
        allocate(comUnwrapTraj(3,nmolsys,nframe))
        allocate(mutraj(3,nmolsys,nframe))
        allocate(unitNormMolTraj(3,nmolsys,nframe))
        allocate(boxtraj(3,nframe))
        allocate(timestamp(nframe))
        allocate(timestep(nframe))
        write(6,*) 'Initialization is done !'
          call date_and_time(values=time_array_t)
          write(6,*) 'time : ', time_array_t
          write(6,*) ''
    else
        write(6,*) 'expanding trajectory by ', nframe
        call date_and_time(values=time_array_t)
        write(6,*) 'time : ', time_array_t
        call expand3D(comtraj,0,0,nframe)
        call expand3D(comUnwrapTraj,0,0,nframe)
        call expand3D(mutraj,0,0,nframe)
        call expand3D(unitNormMolTraj,0,0,nframe)
        call expand2D(boxtraj,0,nframe)
          write(6,*) size(timestamp(:))
        call expand1D(timestamp,nframe)
        call expand1D(timestep,nframe)
          call date_and_time(values=time_array_t)
          write(6,*) 'time : ', time_array_t
          write(6,*) ''
          write(6,*) size(timestamp(:))
    endif


    idx = nxtc
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (istat,ixtc,ixtcGr,dkavg,dr_bins,minridx,trajin,volavg) &
    !$OMP   SHARED (unitNormMolTraj, comtraj, mutraj,chargeAtoms) &
    !$OMP   SHARED (comMolDiff, comUnwrapTraj,timestamp,timestep,boxtraj,grsum,sqsum)
    threadid = omp_get_thread_num()+1
    !write(6,*) threadid
    !$OMP DO SCHEDULE(DYNAMIC)
    do i = 1, nframe
        !$OMP CRITICAL
        call trajin % read
        if (trajin % STAT == 0) then
            ixtc = ixtc + 1
            idx = ixtc
            if(i .eq. 1) then
                if(ixtc .eq. 1) then
                    ! initialize the sin(Qr)/Qr array using the initial box size
                    call calcQr(minval((/ (trajin % box(k,k), k=1,3) /)))
                    dr_bins = minval((/ (trajin % box(k,k), k=1,3) /))/(2.d0*nbins)
                    minridx = nbins
                    write(6,'(a,f12.6,a,i0)') " Time (ps): ", trajin % time, "  Step: ", trajin % STEP
                    write(6,'(a,f12.6,a,i0)') " Precision: ", trajin % prec, "  No. Atoms: ", trajin % NATOMS
                    ! This is the same order as found in the GRO format fyi
                    write(6,'(9f9.5)')  trajin % box(1,1), trajin % box(2,2), trajin % box(3,3), &
                                        trajin % box(1,2), trajin % box(1,3), & 
                                        trajin % box(2,1), trajin % box(2,3), &
                                        trajin % box(3,1), trajin % box(3,2)
                endif
                write(6,*) trajin % NFRAMES
                write(6,*) ''
    
                ! check number of atoms in system
                if(nsysatoms .ne. trajin % NATOMS) then
                  write(6,*) 'number of atoms in ',ifile,'th trajectory file does not match other trajectories'
                  write(6,*) 'nSysAtoms = ',nsysatoms,'  trajin % NATOMS = ', trajin % NATOMS
                  stop
                endif
            endif
    
            volume = trajin % box(1,1) * trajin % box(2,2) * trajin % box(3,3)
            volavg = volavg + volume
     
            ! record current time
            timestamp(idx) = trajin % time
            timestep(idx) = trajin % STEP
    
            do k=1,3
              sysbox(k) = trajin % box(k,k)
              boxtraj(k,idx) = sysbox(k)
            enddo
    
            xposframe = trajin % pos
            icurstat = 0
            !$OMP FLUSH (dr_bins)
        else
            istat = 1
            icurstat = 1
            !$OMP FLUSH (istat)
        endif
        !$OMP END CRITICAL

        if (icurstat == 0) then
            ! Compute and store the trajectory of the center of mass postion,
            ! the rotational principle vector, and the dipole moment vector
            call get_comtraj_frame(sys, xposframe, unitNormMolTraj(:,:,idx), comtraj(:,:,idx), mutraj(:,:,idx))
    
            !call date_and_time(values=time_array_t)
            !write(6,*) 'time : ', time_array_t
    
            halfL = minval(sysbox(:))/2
            maxridx = int(halfL/dr_bins)
            if(maxridx .lt. minridx) then
                minridx = maxridx
            endif
            ! Compute and store the g(r)_ij and the S(q)_ij for each snapshot
            if ( mod(idx,nskipgr) == 0) then
                dk = 2.d0*pi/sysbox
                !$OMP CRITICAL
                ixtcGr = ixtcGr + 1
                dkavg = dkavg + dk
                !$OMP END CRITICAL
                if (bGrCOM) then
                    if (bCalcGr) call get_gr_frame_2D(comtraj(:,:,idx),sysbox,molTypeIdxs, grsum)    ! old method
                    if (bCalcSq) call get_sq_frame(comtraj(:,:,idx),dk,chargeMols,sqsum)
                else
                    if (bCalcGr) call get_gr_frame_2D(xposframe,sysbox,atomTypeIdxs, grsum)    ! old method
                    if (bCalcSq) call get_sq_frame(xposframe,dk,chargeAtoms,sqsum)
                endif
                !if (bCalcSq) call get_sq_frame(xposframe,sysbox,sqsum)
                !call get_gr_frame(xposframe,sysbox,atomPairIdxs, grsum)       ! new method (atomtypeidxs initialized)
            endif
    
    
            !call date_and_time(values=time_array_t)
            !write(6,*) 'time : ', time_array_t
            if (mod(idx,20)==0) write(6,100,advance='no') creturn, idx,'th frame has finished  ' 
        endif
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    nxtc = ixtc

    ! 5. Close the file
    call trajin % close
    end do
 100    FORMAT(A, I8, A)  

    nxtc = ixtc
    volavg = volavg / nxtc        ! average volume in nm^3
    invvol = 1.d0 / volavg        ! inverse volume in nm^-3
    rho = nsysatoms * invvol      ! particle density in nm^-3
    volavg = volavg * 1e-27       ! convert the unit to m^3
    multiple = multiple / (temperature * volavg * 3)
    multiConduct = ELEMENTARY_CHARGE**2 / (6. * BOLTZMANN_CONST * temperature * volavg)
    write(6,*) 'nxtc, volavg', nxtc, volavg
    write(6,*) 'temp mult vol', multiConduct, temperature, volavg

    write(6,*) 'start analyzing the data' 
    call date_and_time(values=time_array_t)
    write(6,*) 'time : ', time_array_t

    ! check the size of trajectory matrix
    isize = size(comtraj(1,1,:))
    if(isize .gt. nxtc) then
      call shrink3D(comtraj,3,nmolsys,nxtc)
      call shrink3D(comUnwrapTraj,3,nmolsys,nxtc)
      call shrink3D(mutraj,3,nmolsys,nxtc)
      call shrink3D(unitNormMolTraj,3,nmolsys,nxtc)
      call shrink2D(boxtraj,3,nxtc)
      call shrink1D(timestamp,nxtc)
      call shrink1D(timestep,nxtc)
    endif

    write(6,*) 'matrix size has adjusted' 
    volavg = volavg * 1e27
    mult2 =  ELEMENTARY_CHARGE**2 * 1e9 * nsyspart / (EPS0 * BOLTZMANN_CONST * temperature * volavg)
    write(6,*) 'mult2 : ',mult2

    ! average the gr matrice and print the result
    if (bCalcGr) then
        grsum = grsum / dble(ixtcGr)
        allocate(grout(nattypepairs+1,nbins))
        nbins = minridx
        !call shrink2D(grsum,nattypepairs+1,nbins)
        grout = grsum

        grsum = grsum * volavg
        if(bGrCOM) then
            call norm_grsum_2D(grsum, nmolMolType)
            call calc_sq_from_grsum(grsum, grout, nmolMolType,chargeMolType,invvol)
        else
            call norm_grsum_2D(grsum, npartInAtType)
            call calc_sq_from_grsum(grsum, grout, npartInAtType,chargeAtType,invvol)
        endif

        strFileName = trim(adjustl(strSqFile)) // '_gr'
        call backupfile(strFileName)
        write(strfmt,'("( F12.6, ",I0,"ES15.7 )")') (nattypepairs+4)
        open(17,file=strFileName)
        write(17,*) "# k       S(k)       S(k)_Lorch         "
        write(6,*) "# k       S(k)       S(k)_Lorch          "
        do k=1, MAX_Q_FORM
          write(17,strfmt) k*kunit, sq(k), sq2(k), sqlorch(k), sqlorch(k)*mult2, sqPart(k,:)
          if ( k .eq. 1) write(6,'(F12.3,3ES15.7)') k*kunit, sq(k), sq2(k), sqlorch(k)
        enddo
        close(17)
        call backupfile(strGrFile)
        open(17,file=strGrFile)
        write(strfmt,'("( F12.6, ",I0,"ES15.7 )")') (nattypepairs+2)
        do i=1, minridx-1
          write(17,strfmt) dr_bins*dble(i), grout(nattypepairs+1,i),grsum(nattypepairs+1,i),grsum(1:nattypepairs,i)
        enddo
        close(17)
    endif
    if(bCalcSq) then
        call backupfile(strSqFile)
        sqsum = sqsum / dble(ixtcGr*nsyspart)
        dkavg = dkavg / dble(ixtcGr)
        dkval = sum(dkavg) / dble(9)
        open(17,file=strSqFile)
        write(17,*) "# k       S(k)         S(k)/k^2          S(k)/k^2/lamdab   "
        write(6,*) "# k       S(k)/k^2/lamdab          "
        nqtot = nqgrid(1)*nqgrid(2)*nqgrid(3)
        allocate(sqavg(sum(nqgrid)*10))
        allocate(sqcnt(sum(nqgrid)*10))
        sqavg = 0.d0
        sqcnt = 0
        do idx=1,nqtot
            i=idx/nqgrid(1)
            j=mod(idx,nqgrid(1))/nqgrid(2)
            k=mod(mod(idx,nqgrid(1)),nqgrid(2))
            do l=1,3
                kval = norm((/ dkavg(l)*i, dkavg(mod(l,3)+1)*j, dkavg(mod(l+1,3)+1)*k /))
                idxcnt = int(kval/dkval + 0.5d0)
                if(idxcnt < size(sqavg)) then
                    sqavg(idxcnt) = sqavg(idxcnt) + sqsum(idx)
                    sqcnt(idxcnt) = sqcnt(idxcnt) + 1
                endif
                write(17,'(F12.3,6ES15.7)') kval, sqsum((idx-1)*3+l), sqsum((idx-1)*3+l)/kval/kval, sqsum((idx-1)*3+l)*mult2/kval/kval, kval*kval/mult2/sqsum((idx-1)*3+l)
            enddo
            if ( k .eq. 1) write(6,'(F12.3,ES15.7)') kval, sqsum(idx*3)*mult2/kval/kval
        enddo
        close(17)
        strFileName = trim(adjustl(strSqFile)) // '_avg'
        call backupfile(strFileName)
        open(17,file=strFileName)
        write(17,*) "# k       S(k)         S(k)/k^2          S(k)/k^2/lamdab   "
        write(6,*) "# k       S(k)/k^2/lamdab          "
        do idx=1,size(sqavg)
            if (sqcnt(idx)>0) then
                kval = dkval* dble(idx)
                sqavg(idx) = sqavg(idx)/sqcnt(idx)
                write(17,'(F12.3,5ES15.7)') kval, sqavg(idx), sqavg(idx)/kval/kval, sqavg(idx)*mult2/kval/kval
            endif
        enddo
        close(17)
    endif
    write(6,*) 'first gr saved'

    ! calculate size of time difference matrices
    ! assume the first two steps have minimum dt
    ! and the trajectories are time-ordered
    dt0 = timestamp(2)-timestamp(1)
    dt = timestamp(nxtc)-timestamp(1)
    ndframe = int(dt/dt0 + 0.00001)

    write(6,*) 'dt is set' , dt
    write(6,*) 'ndframe is set' , ndframe
    ! generate and initialize time difference matrices
!    allocate(idtmat(nxtc-1,nxtc))
    write(6,*) 'idtmat allocated' 
    allocate(unitNormMolt1(3,nmolsys))
    allocate(unitNormMolt2(3,nmolsys))
    write(6,*) 'unitNormMolt1 & unitNormMolDiff allocated' 

!    allocate(xqAtomsCFTime(ndframe))
    allocate(xqcomCFTime(ndframe))
    write(6,*) 'xqcomCFTime allocated' 
    allocate(xqcomDiff(3,nmolsys))
    write(6,*) 'xqcomDiff allocated' 
    allocate(msdTime(ndframe,4,nmoltype+1))
    write(6,*) 'msdTime allocated' 
    allocate(rotACFTime(ndframe,nmoltype+1))
    write(6,*) 'rotACFTime allocated' 
    allocate(rotACFTimeP2(ndframe,nmoltype+1))
    write(6,*) 'rotACFTimeP2 allocated' 
 !   allocate(rotacf(nxtc-1,nxtc,nmolsys))
!    allocate(xqAtomsDiff(nxtc-1,nxtc,nsysatoms,3))
    allocate(xqcomTraj(3,nmolsys,nxtc))
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
        ! unwrap the trajectory
    do ixtc=1, nxtc
        if(ixtc .eq. 1) then
          comUnwrapTraj(:,:,ixtc) = comtraj(:,:,ixtc)
        else if (ixtc .gt. 1) then
          sysbox(:) = boxtraj(:,ixtc)
          do i=1, nmolsys
            comMolDiff(:,i) = comtraj(:,i,ixtc) - comtraj(:,i,ixtc-1)
            comMolDiff(1,i) = comMolDiff(1,i) - NINT(comMolDiff(1,i)/sysbox(1))*sysbox(1)
            comMolDiff(2,i) = comMolDiff(2,i) - NINT(comMolDiff(2,i)/sysbox(2))*sysbox(2)
            comMolDiff(3,i) = comMolDiff(3,i) - NINT(comMolDiff(3,i)/sysbox(3))*sysbox(3)
            comUnwrapTraj(:,i,ixtc) = comUnwrapTraj(:,i,ixtc-1)+comMolDiff(:,i)
          enddo
        endif 
    enddo

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
      unitNormMolt1(:,:) = unitNormMolTraj(:,:,i)
      !$OMP PARALLEL &
      !$OMP PRIVATE(t2,dt,idt,comMolDiff,unitNormMolt2,idxmol,xqcomDiff, &
      !$OMP         xqdiff1,k,idxmoltype,nmolcurr,xdiff,msd,rotacf)
      !$OMP DO
      do j=i+1,nxtc
        t2 = timestamp(j)
        dt = t2 - t1
        idt = int(dt/dt0 + 0.00001)
        if(idt .le. 0) then
          write(6,*) 'idt is less than 0', idt, t1, t2, dt
          cycle
        elseif((idt .gt. 100) .and. (mod(idt,10) .ne. 0) ) then
          cycle
        elseif((idt .gt. 1000) .and. (mod(idt,100) .ne. 0) ) then
          cycle
        endif
        comMolDiff(:,:) = comUnwrapTraj(:,:,i)- comUnwrapTraj(:,:,j)
        unitNormMolt2(:,:) = unitNormMolTraj(:,:,j)

        nDiffTime(idt) = nDiffTime(idt) + 1

        idxmol = 0
        xqcomDiff(:,:) = xqcomTraj(:,:,i)-xqcomTraj(:,:,j)
        xqdiff1 = 0
        do k=1,nmolsys
          idxmoltype = molTypeIdxs(k)
          nmolcurr = nmolMolType(idxmoltype)
          xdiff(:) = comMolDiff(:,k)
          msd = dot_product(xdiff, xdiff)
          msdTime(idt,1:3,idxmoltype) = msdTime(idt,1:3,idxmoltype) + xdiff(:)**2/nmolcurr
          msdTime(idt,1:3,nmoltype+1) = msdTime(idt,1:3,nmoltype+1) + xdiff(:)**2/nmolsys
          msdTime(idt,4,idxmoltype) = msdTime(idt,4,idxmoltype) + msd/nmolcurr
          msdTime(idt,4,nmoltype+1) = msdTime(idt,4,nmoltype+1) + msd/nmolsys
          
          rotacf = dot_product(unitNormMolt1(:,k),unitNormMolt2(:,k))
          rotACFTime(idt,idxmoltype) = rotACFTime(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTime(idt,nmoltype+1) = rotACFTime(idt,nmoltype+1) + rotacf/nmolsys
          rotacf = rotacf*rotacf
          rotACFTimeP2(idt,idxmoltype) = rotACFTimeP2(idt,idxmoltype) + rotacf/nmolcurr
          rotACFTimeP2(idt,nmoltype+1) = rotACFTimeP2(idt,nmoltype+1) + rotacf/nmolsys
          xqdiff1(:) = xqdiff1(:) + xqcomDiff(:,k)
        enddo

        xqcomCFTime(idt) = xqcomCFTime(idt) + dot_product(xqdiff1,xqdiff1)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      write(6,100,advance='no') creturn, i,'th frame has finished  ' 
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
    write(strfmt,'("( F12.3, ",I0,"ES15.7 )")') (nmoltype+1)*4
    ! write (0,1) for auto correlation functions
    write(20,strfmt) 0.d0,rotACFt0
    write(21,strfmt) 0.d0,rotACFt0
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

!!!!! deprecated functionalities
!---- block grpairtypes (deprecated)
!---- was calculating gr pairs which were set in topology file
!---- now all the possible pairs of atoms with atomtype > 0 will be calculated
!!! variables
!  integer            :: ngrpairs
!  integer, allocatable,dimension(:,:) :: npartInPair   ! array for number of particles of each atom type in g(r) pair list (gr_pair_idx, atomtype_idx)
!  integer, allocatable,dimension(:,:) :: grPairIdxs    ! array for g(r) pair idx of each atom pairs in system ( atom_i idx, atom_j idx ) returns 0 if no matching g(r) pair is defined
!
!!! during qform generation
!  ngrpairs = size(sys % pairlist(:,1))
!  allocate(npartInPair(ngrpairs,2))
!  allocate(grPairIdxs(nsysatoms,nsysatoms))
!  grPairIdxs = 0
!  npartInPair = 0
!
!!! generate gr pairs listed in the pair list
!  write(*,*) 'start gr generation'
!  do k=1,ngrpairs
!    !$OMP PARALLEL &
!    !$OMP   DEFAULT (FIRSTPRIVATE) &
!    !$OMP   SHARED (grPairIdxs, npartInPair)
!    !$OMP DO
!    do i=1, nsysatoms
!      if ((atomTypeIdxs(i) == sys % pairlist(k,1)) .or. &
!          (atomTypeIdxs(i) == sys % pairlist(k,2))) then
!        if (atomTypeIdxs(i) == sys % pairlist(k,1)) then
!            !$OMP CRITICAL
!            npartInPair(k,1) = npartInPair(k,1) + 1
!            if (sys % pairlist(k,1) == sys % pairlist(k,2)) npartInPair(k,2) = npartInPair(k,2) + 1
!            !$OMP END CRITICAL
!            itype2inpair = 2
!        else
!            !$OMP CRITICAL
!            npartInPair(k,2) = npartInPair(k,2) + 1
!            !$OMP END CRITICAL
!            itype2inpair = 1
!        endif
!        do j=i+1,nsysatoms
!          if (atomTypeIdxs(j) == sys % pairlist(k,itype2inpair)) then
!            if (itype2inpair == 2) then
!              grPairIdxs(i,j) = k
!              if (sys % pairlist(k,1) == sys % pairlist(k,2)) grPairIdxs(j,i) = k
!            else
!              grPairIdxs(j,i) = k
!            endif
!          endif
!        enddo
!      endif
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!    
!    write(6,*) "pair num : ",k, npartInPair(k,1) , npartInPair(k,2)
!  enddo
!
!!! after matrix size has adjusted
!    do k=1,ngrpairs
!      write(6,*) "pair num : ",k, npartInPair(k,:)
!    enddo
