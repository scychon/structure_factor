!**************************************************
! these are global variables used by many subroutines throughout
! the program
!**************************************************

MODULE variables
   real*8, parameter :: pi = 3.141592654
   real*8, parameter :: DEBYE = 3.33564E-30, EPS0 = 8.85418781762E-12, BOLTZMANN_CONST = 1.380648813E-23
   real*8, parameter :: enm = 0.020819434, ELEMENTARY_CHARGE = 1.60217656535E-19
   real*8, parameter :: kunit = 0.1
   integer, parameter :: MAX_FN=100, MAX_ANAME=5, MAX_N_ATOM=90000 ,  MAX_N_ATOM_TYPE=20
   integer, parameter :: MAX_Q_FORM = 1000
   integer, dimension(3) :: nqgrid
   character*1, parameter :: creturn = achar(13)
   real*8  :: dq_form, dr_bins, sqRmax
   real*8 , dimension( max_q_form ) :: q_grid_form       ! array for list of k values
   real*8 , dimension( max_q_form ) :: form_mol          ! array for x_alp * f_alp(k)
   real*8 , dimension( max_q_form ) :: form_mol2         ! array for x_alp * f_alp(k)^2
   real*8 , dimension( max_q_form ) :: sq, sq2, sqlorch, sqlorch2,sqtest           ! array for S(k)
!   real*8 , dimension( max_q_form ) :: sinqr_over_qr
   real*8 , dimension(:,:), allocatable :: atomtype_form
   real*8 , dimension(:,:), allocatable :: sinqr_over_qr, sqPart, sqPartDirect

   logical            :: bCalcGr=.true., bCalcSq=.true., bGrCOM=.true.

   integer :: nskip,nskipgr     ! skip this number of trajectories between data points
   integer :: nbins     ! number of bins for g(r) and S(q) histograms
  integer :: nmolcat,nmolani,nxtc,nmolsys,nxtcfile,nmoltype
  integer,DIMENSION(:,:,:),ALLOCATABLE:: idxClustIons
  integer,DIMENSION(:,:),ALLOCATABLE:: idxPairIon

  real*8  :: temperature
  real*8,DIMENSION(:,:),ALLOCATABLE:: distPairIon
!  real*8,dimension(:),allocatable :: timestamp
  real*8,dimension(:),allocatable :: corr_IP,corr_IP_cat,corr_IP_ani
  real*8,dimension(:),allocatable :: corr_CP,corr_CP_cat,corr_CP_ani
!  real*8 :: dt,t0
  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strMSDFile, strIPFile, strCPFile,strConductFile,strRotACFFile,strRotACFP2File
  character(len=256)         :: strGrFile, strSqFile, strSqFile2, strAFFExt, strAFFDir
  character(len=256),allocatable,dimension(:) :: strXtcfiles
  integer,DIMENSION(:),ALLOCATABLE:: nXtcFrames

END MODULE variables
