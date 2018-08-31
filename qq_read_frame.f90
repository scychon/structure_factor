module read_frame
    implicit none

contains
  subroutine get_comtraj_frame(sys, xposframe, unitNormMolFrame, comframe, muframe)
    use topol, only: topfile
    use fvector
    use omp_lib

    implicit none
    type(topfile), intent(in) :: sys
    real*8, intent(in) :: xposframe(:,:)   !matrix for positions in the frame (xyz, atom idx)
    real*8, intent(inout) :: unitNormMolFrame(:,:), comframe(:,:), muframe(:,:)   ! frame of com of each molecules (vec_xyz, mol idx)

    integer :: i,j,k,l
    integer :: idxatom, idxmol, idx
    real*8  :: mass, charge, molmass
    integer,dimension(3)   :: idxs
    real*8,dimension(3)    :: xpos,compos   ! temporary position vector for x, center of mass(com)
    real*8,dimension(3)    :: mutot
    real*8, dimension(3,3)    :: xPlaneAtom   ! temporary position vector for the three atoms that compose the molecular plane (vec_xyz, vec idx)

    ! Compute and store the trajectory of the center of mass postion,
    ! the rotational principle vector, and the dipole moment vector
    idxatom = 0
    idxmol = 0
    do i=1, sys % nummoltype
      do j=1, sys % moltype(i) % nummol
        compos = 0
        mutot = 0
        idxs = sys % moltype(i) % molPlaneAtomIdxs

        do k=1, sys % moltype(i) % numatom
          idx = idxatom + (j-1) * sys % moltype(i) % numatom + k
          mass = sys % moltype(i) % atommass(k)
          charge = sys % moltype(i) % atomcharge(k)
          xpos(:) = xposframe(:,idx)
          compos(:) = compos(:) + xpos(:)*mass
          mutot(:) = mutot(:) + xpos(:)*charge
          if (k .eq. idxs(1) ) then
              xPlaneAtom(1,:) = xpos(:)
          elseif (k .eq. idxs(2) ) then
              xPlaneAtom(2,:) = xpos(:)
          elseif (k .eq. idxs(3) ) then
              xPlaneAtom(3,:) = xpos(:)
          end if
        enddo
        unitNormMolFrame(:,idxmol+j) = vecUnitNorm( xPlaneAtom(1,:)-xPlaneAtom(2,:), xPlaneAtom(1,:)-xPlaneAtom(3,:) )

        molmass = sys % moltype(i) % molarmass
        compos(:) = compos(:)/molmass
        comframe(:,idxmol+j) = compos(:)
        muframe(:,idxmol+j) = mutot(:)
      enddo
      idxmol = idxmol + sys % moltype(i) % nummol
      idxatom = idxatom + sys % moltype(i) % nummol * sys % moltype(i) % numatom
    enddo 

    ! check number of atoms in system
    if(sys % numsysatom .ne. idxatom) then
        write(6,*) 'number of atoms does not match other frames'
        write(6,*) 'nSysAtoms = ',sys % numsysatom,'  idxatom = ', idxatom
        stop
    endif
  end subroutine get_comtraj_frame

! initialize the nsysatoms*nsysatoms atompairs atomtype array 
  subroutine set_atompairs(atomTypeIdxs, atomPairIdxs)
    use fvector
    use variables
    use omp_lib

    implicit none
    integer, intent(in) :: atomTypeIdxs(:) ! array for atomtype idx of each atom in system ( atomtype idx )
    integer, allocatable, intent(inout) :: atomPairIdxs(:,:)

    integer :: i,j,k
    integer :: idxatomtype
    integer :: attype1, attype2, natomtype, nattypepairs, nsysatoms
    integer*8 :: idx,ntot

    nsysatoms = size(atomTypeIdxs)
    natomtype = maxval(atomTypeIdxs(:))
    nattypepairs = natomtype*natomtype
    ntot = (nsysatoms*int(nsysatoms-1,8))/2
    write(6,*) 'starting assigning atompairIdxs', nsysatoms, natomtype, nattypepairs, ntot

    if(allocated(atomPairIdxs)) deallocate(atomPairIdxs)
    allocate(atomPairIdxs(4,ntot))
    write(6,*) 'atomPairIdxs allocated ', size(atomPairIdxs(1,:)), (nsysatoms*(nsysatoms-1))/2
    idx = 0
    do i=1, nsysatoms-1
        attype1 = atomTypeIdxs(i)
        if (attype1 .le. 0) cycle
        do j=i+1, nsysatoms
            attype2 = atomTypeIdxs(j)
            if (attype2 .le. 0) cycle
            idx = idx + 1
            atomPairIdxs(1,idx) = i
            atomPairIdxs(2,idx) = j
            atomPairIdxs(3,idx) = (attype1-1) * natomtype + attype2
            atomPairIdxs(4,idx) = (attype2-1) * natomtype + attype1
        enddo
    enddo
    write(6,*) 'atompairIdxs', size(atompairIdxs(1,:))
    call shrink2D(atomPairIdxs,int(4,8),idx)
    write(6,*) 'atompairIdxs', size(atompairIdxs(1,:))
    write(6,*) ''
  end subroutine set_atompairs

! Compute and store the g(r)_ij for each snapshot
  subroutine get_gr_frame(xposframe, sysbox, atomPairIdxs, grsum)
    use fvector
    use variables
    use omp_lib

    implicit none
    real*8, intent(in) :: xposframe(:,:)   !matrix for positions in the frame
    real*8, intent(in) :: sysbox(:)
    integer, intent(in) :: atomPairIdxs(:,:) ! array for atomtype idx of each atom in system ( atomtype idx )
    real*8, intent(inout) :: grsum(:,:)   !matrix for positions in the frame

    integer :: i,j,k,idx
    integer :: idxatomtype, idxAtomPair, natomPairs
    integer :: nsysatoms, nattypepairs
    real*8  :: dr, halfL
    real*8,dimension(3)    :: dxpos   ! temporary position difference vector for dx
    real*8, allocatable, dimension(:,:) :: grframe   !matrix for positions in the frame

    nsysatoms = size(xposframe(1,:))
    nattypepairs = size(grsum(:,1))-1
    natomPairs = size(atomPairIdxs(1,:))

    allocate(grframe(nattypepairs+1,nbins))
    grframe = 0.d0
    halfL = minval(sysbox(:))/2
    do idxAtomPair=1, natomPairs
        i = atomPairIdxs(1,idxAtomPair)
        j = atomPairIdxs(2,idxAtomPair)
        dxpos(:) = xposframe(:,i) - xposframe(:,j)
        dr = getdr(dxpos, sysbox)
        if(dr .lt. halfL) then
            idx = int(dr/dr_bins)
            if (idx .le. nbins) then
                grframe(nattypepairs+1,idx) = grframe(nattypepairs+1,idx) + 2
                idxatomtype = atomPairIdxs(3,idxAtomPair)
                grframe(idxatomtype,idx) = grframe(idxatomtype,idx) + 1
                idxatomtype = atomPairIdxs(4,idxAtomPair)
                grframe(idxatomtype,idx) = grframe(idxatomtype,idx) + 1
            endif
        endif
    enddo

    !OMP CRITICAL
    do i=1,nbins
        grsum(nattypepairs+1,i) = grsum(nattypepairs+1,i) + grframe(nattypepairs+1,i) 
        do j=1,nattypepairs
            grsum(j,i) = grsum(j,i) + grframe(j,i) 
        enddo
    enddo
    !OMP END CRITICAL
  end subroutine get_gr_frame

! Compute and store the g(r)_ij for each snapshot (old method)
  subroutine get_gr_frame_2D(xposframe, sysbox, atomTypeIdxs, grsum)
    use fvector
    use variables
    use omp_lib

    implicit none
    real*8, intent(in) :: xposframe(:,:)   !matrix for positions in the frame
    real*8, intent(in) :: sysbox(:)
    integer, intent(in) :: atomTypeIdxs(:) ! array for atomtype idx of each atom in system ( atomtype idx )
    real*8, intent(inout) :: grsum(:,:)   !matrix for positions in the frame

    integer :: i,j,k
    integer :: idxatomtype, idx
    integer :: attype1, attype2, natomtype
    integer :: nsysatoms, nattypepairs
    real*8  :: dr, halfL
    real*8,dimension(3)    :: xposi, dxpos   ! temporary position difference vector for dx
    real*8, allocatable, dimension(:,:) :: grframe   !matrix for positions in the frame

    nsysatoms = size(xposframe(1,:))
    nattypepairs = size(grsum(:,1))-1
    natomtype = maxval(atomTypeIdxs(:))

    allocate(grframe(nattypepairs+1,nbins))
    grframe = 0.d0
    halfL = minval(sysbox(:))/2
    do i=1, nsysatoms-1
        attype1 = atomTypeIdxs(i)
        if (attype1 .le. 0) cycle
        xposi(:) = xposframe(:,i)
        do j=i+1, nsysatoms
            attype2 = atomTypeIdxs(j)
            if (attype2 .le. 0) cycle
            dxpos(:) = xposi(:) - xposframe(:,j)
            dr = getdr(dxpos, sysbox)
            if(dr .lt. halfL) then
                idx = int(dr/dr_bins)
                if (idx .le. nbins) then
                    grframe(nattypepairs+1,idx) = grframe(nattypepairs+1,idx) + 2
                    idxatomtype = (attype1-1) * natomtype + attype2
                    grframe(idxatomtype,idx) = grframe(idxatomtype,idx) + 1
                    idxatomtype = (attype2-1) * natomtype + attype1
                    grframe(idxatomtype,idx) = grframe(idxatomtype,idx) + 1
                endif
            endif
        enddo
    enddo

    !OMP CRITICAL
    do i=1,nbins
        grsum(nattypepairs+1,i) = grsum(nattypepairs+1,i) + grframe(nattypepairs+1,i) 
        do j=1,nattypepairs
            grsum(j,i) = grsum(j,i) + grframe(j,i) 
        enddo
    enddo
    !OMP END CRITICAL
  end subroutine get_gr_frame_2D

! Compute and store the g(r)_ij for each snapshot (old method)
  subroutine norm_grsum_2D(grsum, npartInType)
    use fvector
    use variables
    use omp_lib

    implicit none
    integer, intent(in)    :: npartInType(:) ! array for number of particles per type (atom / molecule)
    real*8,  intent(inout) :: grsum(:,:)   !matrix for positions in the frame

    integer :: i,j,k
    integer :: idxatomtype, idx
    integer :: attype1, attype2, ntype, ntypepairs, nsyspart

    nsyspart = sum(npartInType)
    ntype = size(npartInType)
    ntypepairs = ntype**2

    do i=1,ntype
        write(6,*), 'ntype : ',npartInType(i),nsyspart
    enddo
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (grsum,npartInType)
    !$OMP DO
    do i=1,nbins
        grsum(ntypepairs+1,i) = grsum(ntypepairs+1,i)*1.d0 / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * dble(nsyspart * nsyspart) )
        do j=1,ntypepairs
            grsum(j,i) = grsum(j,i)*1.d0 / ( (4.d0/3.d0)*pi* ((i+1)**3 - i**3)*dr_bins**3 * dble(npartInType((j-1)/ntype + 1) * npartInType(mod((j-1),ntype)+1)) )
        enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine norm_grsum_2D

  subroutine calc_sq_from_grsum(grsum,grout,npartInType,scaleParticles,invvol)
    use fvector
    use variables
    use omp_lib

    implicit none
    integer, intent(in)    :: npartInType(:) ! array for number of particles per type (atom / molecule)
    real*8,  intent(in)    :: grsum(:,:), grout(:,:), scaleParticles(:), invvol   !matrix for positions in the frame

    integer :: i,j,k
    integer :: idxatomtype, idx
    integer :: type1, type2, ntype, ntypepairs, nsyspart
    real*8  :: rmax, kval, rval, ival
    real*8  :: sqval, sqval2, sqbinpre, sqbinpre2, pirSinPir

    nsyspart = sum(npartInType)
    ntype = size(npartInType)
    ntypepairs = ntype**2

    sq = 0
    rmax = dr_bins * nbins
    do i=1,ntype
        sq = sq+scaleParticles(i)*scaleParticles(i)*npartInType(i)/(nsyspart*1.d0)
        write(6,*), 'ntype : ',scaleParticles(i),npartInType(i),nsyspart
    enddo
    sq2 = sq
    sqlorch = sq
    !write(6,*) 'sqval : ',sq
    !$OMP PARALLEL &
    !$OMP   DEFAULT (FIRSTPRIVATE) &
    !$OMP   SHARED (sq,sq2,sqPart,sqPartDirect,sqlorch,scaleParticles,npartInType,grsum)
    !$OMP DO
    do k=1, MAX_Q_FORM
      kval = k*kunit
      do j=1, ntypepairs
        type1 = (j-1)/ntype +1
        type2 = mod((j-1),ntype)+1
        sqval = scaleParticles(type1) * scaleParticles(type2)/(nsyspart*1.d0)
        sqval2 = sqval * invvol * npartInType(type1) * npartInType(type2)
        if (k==1) then
            write(6,*) type1, type2, scaleParticles(type1),scaleParticles(type2)
        endif
        do i=1, nbins-1
          ival = i + 0.5d0
          rval = dr_bins * ival
          pirSinPir = dsin(pi * rval / rmax ) / (pi * rval / rmax)
          sqbinpre = 4.d0 * pi * rval * dsin( kval * rval) / kval
          sqbinpre2 = sqbinpre * pirSinPir
          sqPartDirect(k,j) = sqPartDirect(k,j) + sqval * grout(j,i) * dsin( kval * rval) / (kval * rval)
          sq(k) = sq(k) + sqval * grout(j,i) * dsin( kval * rval) / (kval * rval)
          sqPart(k,j) = sqPart(k,j) + dr_bins * sqval2 * sqbinpre * (grsum(j,i) - 1.d0)
          sq2(k) = sq2(k) + dr_bins * sqval2 * sqbinpre * (grsum(j,i) - 1.d0)
          sqlorch(k) = sqlorch(k) + dr_bins * sqval2 * sqbinpre2 * (grsum(j,i) - 1.d0 )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine calc_sq_from_grsum
    
! Compute and store the exp(-i q r) for each snapshot
  subroutine get_sq_frame_kvec(xposframe, scaleParticles, kvec, sqsum)
    implicit none
    real*8, intent(in) :: xposframe(:,:),scaleParticles(:)   !matrix for positions in the frame
    real*8, intent(in) :: kvec(3)
    real*8, intent(out) :: sqsum

    integer :: iperm,idx,nperm
    logical :: bPermute(3), bRun(8)
    real*8  :: qr, cossum(8), sinsum(8)
    real*8,dimension(3)    :: xpos,kvecPermutation   ! temporary position difference vector for dx

    cossum = 0
    sinsum = 0
    nperm = 0
    bPermute(:) = (kvec(:) .gt. 0)
    bRun = .false.
    bRun(1:4) = (/ .true., bPermute(3), bPermute(2), (bPermute(2) .and. bPermute(3)) /)
    bRun(5:8) = bPermute(1) .and. bRun(1:4)
    do iperm=1,8
        kvecPermutation = kvec
        if (iperm > 4) kvecPermutation(1) = -kvec(1)
        if (mod(iperm-1,4) .ge. 2) kvecPermutation(2) = -kvec(2)
        if (mod(iperm,2) == 0) kvecPermutation(3) = -kvec(3)
        if (bRun(iperm)) then
            nperm = nperm + 1
            do idx=1, size(xposframe(1,:))
                if (scaleParticles(idx) .ne. 0) then
                    xpos = xposframe(:,idx)
                    qr = dot_product(kvecPermutation,xpos)
                    cossum(iperm) = cossum(iperm) + scaleParticles(idx)*dcos(qr)
                    sinsum(iperm) = sinsum(iperm) + scaleParticles(idx)*dsin(qr)
                endif
            enddo
        endif
    enddo
    sqsum = sum(cossum*cossum + sinsum*sinsum)/dble(nperm)
  end subroutine get_sq_frame_kvec
    
  subroutine get_sq_frame(xposframe, dk, scaleParticles, sqsum)
    use variables
    implicit none
    real*8, intent(in) :: xposframe(:,:)  !matrix for positions in the frame
    real*8, intent(in) :: dk(3)           ! kvector unit from box dimension
    real*8, intent(in) :: scaleParticles(:)   !matrix for atomic constribution (charge for charge-charge structure factor and form_factor for neutron or xray structure factor)
    real*8, intent(inout) :: sqsum(:)

    integer :: i,j,k,idx,nqtot
    real*8  :: sumval,sumval2,sumval3
    real*8,dimension(3)    :: kvec   ! temporary position difference vector for dx

    nqtot = nqgrid(1)*nqgrid(2)*nqgrid(3)
    do idx=1,nqtot
        i=idx/nqgrid(1)
        j=mod(idx,nqgrid(1))/nqgrid(2)
        k=mod(mod(idx,nqgrid(1)),nqgrid(2))
        kvec = (/ dk(1)*i, dk(2)*j, dk(3)*k /)
        call get_sq_frame_kvec(xposframe,scaleParticles,kvec,sumval)
        if ((i==j) .and. (j==k)) then
            sumval2 = sumval
            sumval3 = sumval
        else
            kvec = (/ dk(1)*k, dk(2)*i, dk(3)*j /)
            call get_sq_frame_kvec(xposframe,scaleParticles,kvec,sumval2)
            kvec = (/ dk(1)*j, dk(2)*k, dk(3)*i /)
            call get_sq_frame_kvec(xposframe,scaleParticles,kvec,sumval3)
        endif
        !$OMP CRITICAL
        sqsum(idx*3-2) = sqsum(idx*3-2) + sumval
        sqsum(idx*3-1) = sqsum(idx*3-2) + sumval2
        sqsum(idx*3) = sqsum(idx*3) + sumval3
        !$OMP END CRITICAL
    enddo
  end subroutine get_sq_frame
end module read_frame
