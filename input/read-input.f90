!> This is an primitive fortran input reader 
!> right now it has flexibility with the ordering,
  !> and some values have defaults, so they need not be specified
  !> Comments are technically supported, any keywords that aren't recognized won't match anything
  ! and thus nothing will be done with them
  
  !> This module contains core parameters and should not need to be changed
  module params
    ! 64 and 32 bit integers
    integer, parameter :: i64 = selected_int_kind(15)
  !  integer, parameter :: i64 = selected_int_kind(6)
    integer, parameter :: i32 = selected_int_kind(6)
    ! explicit definition of double precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! quadruple precision, if ever needed
    integer, parameter :: qp = selected_real_kind(33, 4931)
    ! some complex numbers
    complex(dp), parameter :: c0 = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cr = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: c1 = (1.0_dp, 1.0_dp)
    ! pi to a nearly unnecessary amount of precision
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: thresh_time_zero = 1.0d-6
    real(dp), parameter :: thresh_parameter_zero = 1.0d-6    
    ! the mass of an electron in au
    !> never actually used
  !  real(dp), parameter :: mass = 1.0_dp
    ! atomic-time conversion to fs
    real(dp), parameter :: au2fs = 0.0241888_dp
    real(dp), parameter :: proelec_mratio = 1836.15267389_dp
    real(dp), parameter :: neutron_elec_mratio = 1838.683661589_dp
  end module params
  
  
  !> these are the input variables, this module declares all input variables needed to start MCEND
  !> it also contains subroutines that read in these variables from an input file
  module inputvars
    use params

  !  use hdf5
  
      implicit none
  
      ! these variables may at some point be merged into the standard input
      ! file, but for now they'll live here
  
      !> initial step-size (fs)
      real(dp), parameter :: dt0 = 1.0d-5
      !> adaptive step-size increase/decrease limits
      real(dp), parameter :: facdec = 0.333333_dp
      real(dp), parameter :: facinc = 6.0_dp
      !> further adaptive step-size parameter
      real(dp), parameter :: beta = 0.04_dp
      !> CAP parameters eta, naocap
      !> naocap: nr of first AO with CAP
      !> this is the electronic CAP
      real(dp), parameter :: etacap = 0.0_dp
      integer,  parameter :: naocap = 26
      !> density matrix regularization parameter
      real(dp), parameter :: epsreg0 = 1.0d-10
      ! tolerance for adaptive step-size control
      real(dp), parameter :: tol = 1.0d-10
      !> no excitations of core-electrons? <- part of non-existent code
  !    logical,  parameter :: nece = .false.
      integer :: nrfrorb
      !> energy shift of electronic one-particle Hamiltonian
      real(dp) :: shifte
      !> wgrid: internuclear distance at which the nuclear CAP starts
      real(dp) :: wgrid
      !> order of the CAP - usually 2 or 3
      integer  :: bn
      !> strength of the CAP - usually 2 or 3
      real(dp) :: etan
      !> sign of the CAP - has to be switched for imagt
      integer  :: wsign
      !> no of electronic basis functions, and no of electrons
      integer :: nrprime, nel
      
      !integer :: nrprime_alpha,nrprime_beta    
      !integer :: nrprime_spinorbital
      
      integer :: nel_alpha, nel_beta
      integer :: nel_spinorbital    
      
      !> no of nuclear grid points and single particle functions
      integer :: nrprimn, nrspf
      !> a list of nrprimn for each dof, N.B., to be merged
      !>
      integer, allocatable :: nrprimn_poly(:)      
      !> no of orbitals, nrorb is will be used for open-shell
      integer :: nrorb
      integer :: nrorb_alpha, nrorb_beta
      integer :: nrorb_spinorbital !nrorb_so,
      integer :: nrorb_init_alpha, nrorb_init_beta    
      
      !integer :: norb
      ! better to be unitied as norb -> nrorb
     ! integer :: norb_alpha, nrorb_beta
     ! integer :: norb_init_alpha, nrorb_init_beta    
     ! integer :: norb_spinorbital

      ! frozen_core and virtual 
      ! frozen virtual has not been tested !!
      integer :: nrorb_fc, nrorb_fc_spinorbital
      integer :: nrorb_fv, nrorb_fv_spinorbital

      !> number of explicit nuclear sampling points
      integer :: nrensp
      !> nr of independent A-vector elements
      ! cw: nrindep = nint(combinatorial(2*nrorb, nel)) 
      ! is this electronic part?
      ! Q:should nr of A becomes lA(nrindep*nrspf)?
      integer :: nrindep
     ! integer :: nrindep_alpha, nrindep_beta
      integer :: nrindep_spinorbital
      integer :: nrindep_frzorb, nrindep_frzorb_spinorbital
      !> starting point on nuclear grid
      integer :: inigrd
      ! read from h5 files 1 = yes, 0 = no
  !    integer :: readfromh5
      !> multiplicity -> not yet used
      integer :: mult
      ! use precomputed basis or moving basis for electrons
      !> moving basis implementation was not successful
      !> will only have first option until we wish to revisit it or remove it all together
      ! basis = 1 precomputed basis at nrensp positions
      ! basis = 2 atom-centered basis, nrensp=nrprimn
      integer :: basis
      ! read in initial state? (2 => restart, 0 => no restart,
      !                         1 => generate additional SPF)
      integer :: restart
      ! autocorrelation function? (0=no, 1=yes, 2=after lpw)
      integer :: aucofu
      integer :: max_nrindep
    !  integer :: max_nrindep_alpha, max_nrindep_beta
  !    integer :: max_nrindep_spinorbital
      integer :: max_nrindep_spinorbital
      
      integer :: max_nrindep_frzorb
      integer :: max_nrindep_frzorb_spinorbital
      
      
      integer :: max_nrindep_2
      integer :: max_nrindep_2_spinorbital
      integer :: max_nrindep_2_frzorb, max_nrindep_2_frzorb_spinorbital
      
      integer :: nrshf
     ! integer :: nrshf_alpha, nrshf_beta
     ! integer :: nrshf_spinorbital
      integer :: nrshf_spinorbital
      integer :: nrshf_frzorb,nrshf_frzorb_spinorbital   
      integer :: nrdhf, nrdhf_spinorbital
      integer :: nrdhf_frzorb,  nrdhf_frzorb_spinorbital    
      
      !> need to figure out what this stands for, this name is also used in MCTDH
      !> cw: seems nrshf = nrshdl 
      integer :: dgldim
      integer :: dgldim_spinorbital
      
      !> how many threads to use
      integer(i32) :: numthreads
      !> mean field hash table size
      integer :: hmfsize
      ! number of primitive nuclear basis functions, min and spacing of grid
      real(dp) :: dr, rmin
      ! I will merge them...
      ! it may be better to be named differently in the initial stage
      ! otherwise interrupt many existing code
      real(dp), allocatable :: dr_poly(:), rmin_poly(:)
      ! total propagation time (fs)
      real(dp) :: tfin0
      ! if this is a restart this can be used to pick up where the laser left off
      real(dp) :: t_initial
      ! output interval (fs)
      real(dp) :: tout0
      ! laser pulse parameters, sin**2*cos shape
      !> pulse width/duration (fs)
      real(dp), allocatable :: lpw0(:)
      ! pulse height (Eh/ea0)
      real(dp), allocatable :: lph(:)
      !> photon energy (hartree)
      real(dp), allocatable :: lfreq(:)
      !> atomic masses
      !> will be generalized into what kind of convention?
      !> m_poly(:)?
      real(dp) :: m2, m1
      real(dp), allocatable :: m_poly(:)      
      ! real(dp), allocatable :: amasses(:)
      !> nuclear charge
      real(dp) :: N2, N1
      real(dp), allocatable  :: N_poly(:)      
      ! real(dp), allocatable :: Ncharge(:)
      real(dp) :: A2N, A1N
      ! the need for python-mcend?
      real(dp), allocatable  ::  A_poly(:)         
      !> reduced nuclear mass
      real(dp) :: massn
      !> pulse polarization direction(s)
      real(dp), allocatable :: lpolar(:,:)
  
  !    type lpulse
  !      real(dp), allocatable :: x(:)
  !      real(dp), allocatable :: y(:)
  !      real(dp), allocatable :: z(:)
  !    end type
  !    type(lpulse) :: lpolar
      !> number of quantum degree of freedom
      integer :: ndof
      !> number of non-zero eri
      integer :: nint2e
      !> how to scale the hmfsize
      integer :: hmfscale
      ! name of the compound
      character(10) :: cmpdname
      ! folder with the integrals
      character(300) :: intdir
      !> file with the grid points
      character(300) :: rspfile
      !> initial guess file
      character(300) :: initguessf
      !> path to the hmat file with values generated from orthog routine
      character(300) :: scfv_path
      !> path to the two electron integrals
      character(300) :: n2int_path
      !> path to the mcend code top
      character(300) :: mcend_top
      character(300) :: path2ints
      ! print wave function every tout?
      logical :: logwavef
      ! diagonalize the CI matrix to get energy eigenstates
      logical :: logev
      ! imaginary time propagation? (imagt=(1,0) => normal propagation)
      !                             (imagt=(0,-1) => relaxation)
      complex(dp) :: imagt
      ! freeze nuclear orbitals
      logical :: freeze_nuc
      ! read integrals in the standard (old) way
      !> now the old way
      logical :: read_stand
      ! use the nuclear CAP
      logical :: use_wcap
      ! shift the energy towards zero automatically
      logical :: autoshifte
      ! control open-shell unrestricted or closed-shell restricted methods
      logical :: do_openshell
      ! control polyayomic module, debug at this stage
      logical :: do_polyatomics      
      ! move one particle operator, so far Te, outside meanfields
      logical :: mf_reduction      
      !> number of pulses, determined from pulse params
      integer :: nrpulses
      !> number of parallel electrons
      integer ::  nsz 
      !> similar to do_openshell, used in previous version; may be extended beyond 0/1 for configuration state function (CSF) approach
      integer :: flag_spinorbital
      !> polyatomic version or not, in progress
      !> may just use do_polyatomics
      integer :: flag_poly       
      !> control frozen core (works) and frozen virtual(not working)
      integer :: flag_fc, flag_fv
      !> control print level for intermediate wall time
      integer :: print_level
  !    open(10,'.mcendrc')
  !    read(10,*) mcend_top
  !    close(10)
  
    contains
  
    !> This subroutine reads the input file
    subroutine read_file(filename)
      implicit none
  
      character(len=*), intent(in) :: filename
      integer :: inputio, ios, ios2
      character(1) :: temp
      character(20) :: temp2
      character(300) :: iom
  !    character (len=300) :: line
      integer :: i, j, k
     ! integer :: lcounter, icounter
      integer :: input_len, nrows
      integer :: nrownum, nrowchar, nrownum2
      integer :: scfvals
      real(dp), allocatable :: scf_temp(:)
      character(1) :: line
  
      character(100), allocatable :: charvars(:)
      character(100), allocatable :: linevars(:) 
      character(6), allocatable :: lineoutvars(:)
      
      character(50), allocatable :: keywargs(:)
      character(50), allocatable :: keywargs2(:)
      real(dp), allocatable :: inpvars(:)
      integer :: in, natom
      
      allocate(nrprimn_poly(ndof))
      allocate(rmin_poly(ndof))
      allocate(dr_poly(ndof))
      
  !    call get_environment_variable("MCEND_TOP", mcend_top)
      call get_environment_variable("PWD", mcend_top)
  
      ! get the number of lines in a file w/o using a external bash call
      ios = 0
      ios2 = 0
      
      write (*,*) 'filename', filename
      
      open(newunit=inputio, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iom)
      if (ios /= 0) then
        write(*,*) 'Fatal error!!!', trim(iom)
        stop
      endif
  
      input_len = 0
      if (ios == 0) then
        do while (ios2 == 0)
          if (ios2 == 0) then
            read(inputio,*,iostat=ios2) temp
            input_len = input_len + 1
          end if
        enddo
      endif
      input_len = input_len - 1
      rewind(inputio)
  
      nrows = input_len
  
      allocate(keywargs(nrows))
      allocate(keywargs2(nrows))
  
      ! this is for implementation of comments, a bit too complex
      ! for the immediate purpose, but will revisit
      ReadComments: do
         read (inputio,'(a)') line
         if (line (1:1) /= "#") exit ReadComments
      end do ReadComments
      backspace (inputio)
  
      ! this is for getting the number of lines before the character declarations
      nrownum2 = 0
      do i=1, nrows
        read(inputio,*) keywargs2(i)
        if (trim(keywargs2(i)) == '@chars') then
          nrownum = i
  !        cycle
        else if (trim(keywargs2(i)) == '@pulseparams') then
          nrownum2 = i
  !      else if (trim(keywargs2(i)) == '')
        endif
      enddo
  
  !    nrowchar = nrows - nrownum - nrownum2
      nrowchar = nrows - nrownum! + nrownum2
  
      allocate(inpvars(nrownum))
      allocate(charvars(nrowchar))
      allocate(linevars(nrownum))
  
      rewind(inputio)
  
      ! load in integers and reals
      do j=1, nrownum - 1
  !    do while (trim(keywargs2(j)) /= '@chars')
  ! reading format, keywords, '=' for temp, and variables
        read(inputio, *) keywargs(j), temp, inpvars(j)
  !     write (*,*) 'aaa', j, inpvars(j)
      enddo
      
      rewind(inputio)
  
      ! linevars to contain more entries in N
      do j=1, nrownum - 1
  !    do while (trim(keywargs2(j)) /= '@chars')
  ! reading format, keywords, '=' for temp, and variables
        read(inputio,"(A80)") linevars(j) !, temp, inpvars(j)
   !     write (*,*) 'bbb', j, linevars(j)
      enddo
        
      
  
  !    rewind(inputio)
  !    do j=1, nrownum - 1
  !    do while (trim(keywargs2(j)) /= '@chars')
  ! reading format, keywords, '=' for temp, and variables
  !      write (inputio,*) 
  !      write(*,*) j, linevars(j)
  !    enddo  
  
  
  !    do while (trim(keywargs2(j+1+icounter)) /= '@end')
  !    if (trim(keywargs2(j+icounter)) == '@end') then
      ! read in that @chars line
      read(inputio,*) keywargs(nrownum)
  
      ! read in character variables
      k = 1
  !    do j = nrownum+1, nrows
      do j = nrownum+1, nrownum2-1
        read(inputio,*) keywargs(j), temp, charvars(k)
        k = k + 1
      enddo
  
      nrpulses = nrows - nrownum2 - 1
  
      allocate(lpolar(3,nrpulses))
      allocate(lph(nrpulses))
      allocate(lpw0(nrpulses))
      allocate(lfreq(nrpulses))
  
      lpolar(:,:) = 0.0_dp
  
  !    read(inputio,*) keywargs(nrownum2)
      ! discard first two lines as they are strings
      read(inputio,*) temp
      read(inputio,*) temp
      do k = 1, nrpulses
        read(inputio,*) lph(k), lpw0(k), lfreq(k), lpolar(1,k), lpolar(2,k), lpolar(3,k)
      enddo
  
      close(inputio)
  
      ! assign initial values of these variables
      nrfrorb = 0
      wgrid = 0.0_dp
      ! hard-coding this to always use static basis set until development resumes
      ! for moving basis set or we give up on that route
      basis = 1
      shifte = 0.0_dp
      !default singlet
      ndof = 1
      mult = 1
      !> us
      hmfscale = 3000
      t_initial = 0.0_dp
      A1N = 0.0_dp
      A2N = 0.0_dp
  
      nrorb_fc = 0
      nrorb_fc_spinorbital = 0
      nrorb_fv = 0      
      nrorb_fv_spinorbital = 0
      ! default parameters
      flag_spinorbital = 0
      flag_poly = 0
      nsz = 0
      print_level = 0
      ! check the key word arguments for the following,
      ! and assign reals, integers accordingly
      
      write (*,*) 'nrownum', nrownum 

      do i=1, nrownum - 1

        !write (*,*) 'trim',i, trim(keywargs(i)) 
        ! more often read from basis set information  
        if (trim(keywargs(i)) == 'nrprime') then
          nrprime = dint(inpvars(i))
          write (*,*) 'nrprime in read-input', nrprime
        elseif (trim(keywargs(i)) == 'nel') then
          nel = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'multiplicity') then
          mult= int(inpvars(i))   
          nsz = mult - 1
         ! write (*,*) 'read multiplicity', mult, nsz  
        elseif (trim(keywargs(i)) == 'natom') then
          natom = int(inpvars(i))         
        elseif (trim(keywargs(i)) == 'nrprimn') then
          nrprimn = dint(inpvars(i))
         
  !      elseif (trim(keywargs(i)) == 'nint2e') then
  !        nint2e = dint(inpvars(i))
  
        elseif (trim(keywargs(i)) == 'nrspf') then
          nrspf = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrorb') then
          nrorb = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrorb_fc') then
          nrorb_fc = int(inpvars(i))
          ! current convention
          nrorb_fc_spinorbital = nrorb_fc * 2
          ! nrorb_fc_spinorbital = 0
          !write (*,*) 'nrorb_fc', nrorb_fc  
          
       ! Oct-26 by cw: postpone this option for a while 
       ! elseif (trim(keywargs(i)) == 'nrorb_fv') then
       !   nrorb_fv = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrensp') then
          nrensp = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'dr') then
          dr = inpvars(i)
        elseif (trim(keywargs(i)) == 'rmin') then
          rmin = inpvars(i)
        elseif (trim(keywargs(i)) == 'shifte') then
          shifte = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'inigrd') then
          inigrd = dint(inpvars(i))
  
        ! WARNING: Nuclear CAP parameters
        elseif (trim(keywargs(i)) == 'wgrid') then
          wgrid = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'wstrength') then
          etan = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'worder') then
          bn = dint(inpvars(i))
  
        !> omitting moving basis set options until it can be revisited
  !      elseif (trim(keywargs(i)) == 'basis') then
  !        basis = int(inpvars(i))
  !        basis = 1
  
        !> omitting read hdf5 data until the dark day on which
        !> we must revist it
  !      elseif (trim(keywargs(i)) == 'readh5') then
  !        readfromh5 = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'restart') then
          restart = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'aucofu') then
          aucofu = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'numthreads') then
          numthreads = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'tfin0') then
          tfin0 = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 't_initial') then
          t_initial = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'tout0') then
          tout0 = dble(inpvars(i))
  
  !      elseif (trim(keywargs(i)) == 'lpw0') then
  !        lpw0(1) = inpvars(i)
  !      elseif (trim(keywargs(i)) == 'lph') then
  !        lph(1) = inpvars(i)
  !      elseif (trim(keywargs(i)) == 'lfreq') then
  !        lfreq(1) = inpvars(i)
        !> plan to be moved into input coordinate
        !> specify the proton mass
        !> to be done 
        elseif (trim(keywargs(i)) == 'N') then
          !N1 = inpvars(i)        
          !N_poly = inpvars(i) 
          !by default, only the first entry will be readed :(
          !write (*,*) 'reading N list', linevars(i)!, N_poly
          ! so far, only support space interval, not comma
          lineoutvars = string_to_integers(linevars(i), ' ') 
          !write (*,*) size(lineoutvars), lineoutvars
          allocate(N_poly(size(lineoutvars)-2))
          
          do j=3, size(lineoutvars)
              read( lineoutvars(j), * ) N_poly(j-2)  
          end do          
          write (*,*) 'N_poly', N_poly
        elseif (trim(keywargs(i)) == 'N1') then
          N1 = inpvars(i)
          write (*,*) 'reading N1 list', inpvars(i), 'pp', N1
        elseif (trim(keywargs(i)) == 'N2') then
          N2 = inpvars(i)
          
        elseif (trim(keywargs(i)) == 'A') then

          lineoutvars = string_to_integers(linevars(i), ' ') 
          allocate(A_poly(size(lineoutvars)-2))
          do j=3, size(lineoutvars)
              read( lineoutvars(j), * ) A_poly(j-2)  
          end do          
          write (*,*) 'A_poly', A_poly          
          
        !> specify the neutron mass
        elseif (trim(keywargs(i)) == 'A1N') then
          A1N = inpvars(i)
        elseif (trim(keywargs(i)) == 'A2N') then
          A2N = inpvars(i)
        !> [not yet fully working] Number of frozen orbitals
        elseif (trim(keywargs(i)) == 'nrfrorb') then
          nrfrorb = dint(inpvars(i))
          !> multiplicity 
        elseif (trim(keywargs(i)) == 'mult') then
          mult = dint(inpvars(i))
          !> scaling for hmf,
        elseif (trim(keywargs(i)) == 'hmfscale') then
          hmfscale = dint(inpvars(i))
          
        elseif (trim(keywargs(i)) == 'print_level') then
          !> 0: fewer print; 1: print time for intermediate steps
          print_level = int(inpvars(i))     
          write (*,*) 'print_level', print_level
        endif
        
        
      enddo
  
      !> crude fail safe for isotope assignment, so in case neutrons number is not
      !> given, assume num protons = number of neutrons, unless the proton count is
      !> 1 which is hydrogen and won't have a neutron in the most common isotope
      !> eventually should put this into an array form
      
      !if (A1N == 0.0_dp .and. N1 /= 1.0_dp) A1N = N1
      if ( dabs(A1N - 0.0_dp) < thresh_parameter_zero .and. dabs(N1 - 1.0_dp) > thresh_parameter_zero ) A1N = N1      
      !if (A2N == 0.0_dp .and. N2 /= 1.0_dp) A2N = N2
      if (dabs(A2N - 0.0_dp) < thresh_parameter_zero .and. dabs(N2 - 1.0_dp) > thresh_parameter_zero ) A2N = N2      
  
  
      !> assume nuclear orbitals are not frozen by default (for 1 SPF this equates
      !> to fixed nuclei)
      freeze_nuc = .false.
  !    read_stand = .true.
      !> now false by default, legacy feature for testing the old gamess basis
      !> sets
      read_stand = .false.
      use_wcap = .false.
      autoshifte = .false.
      do_openshell = .false.
      mf_reduction = .false. 
      ! check the keywords for the character vars, and assign based on input
      do i=1, nrowchar
        j = i + nrownum
        if (trim(keywargs(j)) == 'compound') then
          cmpdname = charvars(i)
        elseif (trim(keywargs(j)) == 'intpath') then
          intdir = charvars(i)
        elseif (trim(keywargs(j)) == 'autoshifte') then
          if (trim(charvars(i)) == 'yes') then
            autoshifte = .true.
          else
            autoshifte = .false.
          endif
        elseif (trim(keywargs(j)) == 'use_ncap' ) then
          if (trim(charvars(i)) == 'yes') then
            use_wcap = .true.
          else
            use_wcap = .false.
          endif
        elseif (trim(keywargs(j)) == 'itimeprop') then
          if (trim(charvars(i)) == 'yes') then
            imagt = (0.0_dp, -1.0_dp)
            wsign = -1
          else
            imagt = (1.0_dp, 0.0_dp)
            wsign = 1
          endif
  
        elseif (trim(keywargs(j)) == 'logev') then
          if (trim(charvars(i)) == 'yes') then
            logev = .true.
          else
            logev = .false.
          endif
        ! experimental
        elseif (trim(keywargs(j)) == 'freeze_nuc') then
          if (trim(charvars(i)) == 'yes') then
            freeze_nuc = .true.
          endif
        elseif (trim(keywargs(j)) == 'read_stand') then
          if (trim(charvars(i)) == 'no') then
            read_stand = .false.
          elseif (trim(charvars(i)) == 'yes') then
            read_stand = .true.
          endif
        elseif (trim(keywargs(j)) == 'logwavef') then
          if (trim(charvars(i)) == 'yes') then
            logwavef = .true.
          else
            logwavef = .false.
          endif
        elseif (trim(keywargs(j)) == 'do_openshell') then
          if (trim(charvars(i)) == 'yes') then
            !do_openshell = .true.
            !open-shell maybe done by spatial orbital spin-adapted vs spinorbital unrestricted
            !int type for the above two cases, though spin adapted spatial not implementated yet
            flag_spinorbital = 1     
          else
            !do_openshell = .false.
            flag_spinorbital = 0
          endif
        elseif (trim(keywargs(j)) == 'do_polyatomics') then
          if (trim(charvars(i)) == 'yes') then
             do_polyatomics = .true.     
          else
            !do_openshell = .false.
            do_polyatomics = .false.
          endif        
  !      elseif (trim(kewargs(j)) == 'basis_path') then
  !        basis_path_dir = charvars(i)
        elseif (trim(keywargs(j)) == 'mf_reduction') then
          if (trim(charvars(i)) == 'yes') then
            mf_reduction = .true.     
          else
            mf_reduction = .false.
          endif
  !      elseif (trim(kewargs(j)) == 'basis_path') then
  !        basis_path_dir = charvars(i)  
        endif
      enddo
  
      !> get current working directory
      intdir = 'basis_library/'//trim(intdir)
  
  !      xeriints = trim(cmpdname)//'-ints-01.dat'
      scfv_path = trim(intdir)//'/scf-guess-'//trim(cmpdname)//'.dat'
      n2int_path = trim(intdir)//'/info-'//trim(cmpdname)//'.dat'
      rspfile = trim(intdir)//'/rsp_'//trim(cmpdname)
  
      !> open info file to get these parameters
      !>cw: NEED TO EXTEND FOR POLYNAMICS
      !>cw: read as input
      !>cw: may-4, why call n2? nuclei n2?
      open(20,file=trim(n2int_path))
      read(20,*) temp2, temp, ndof
      if (temp2 == 'ndof') then
        ndof = ndof
        if (ndof > 1) then
          flag_poly = 1
        end if
        read(20,*) temp2, temp, nint2e 
        
      write (*,*) 'dof check', ndof  
        
      !adopt the previous version, only for diatomics
      ! if I evolve all format into polyatomic, the code could be more concise
      ! not sure if this is a good practice
      ! previous     read(20,*) temp2, temp, nint2e
      else if (temp2 == 'nint2e') then
        nint2e = ndof
      else
        write (*,*) 'not ndof/nint2e in the first line, not supported'
        stop
      end if
      
      
      ! so far, the precomputed basis for polyatomic is merged so to speak
      ! a single nrprime, for combinational positions
      read(20,*) temp2, temp, nrprime
      !again, to change or not to changel; do_polyatomics or flag_poly
      !it seems by ndof to control is better
      !if I change nrprimn into an array, the subsequent section needs to be changed
      !keep it for a while
      !then merge
      !it should, not rhf vs uhf
      if (flag_poly==0) then
        read(20,*) temp2, temp, nrprimn
      else
        read(20,*) temp2, temp, (nrprimn_poly(in), in=1,ndof)
      end if
      
      read(20,*) temp2, temp, nrensp
      
      if (flag_poly==0) then
        read(20,*) temp2, temp, rmin
        read(20,*) temp2, temp, dr
      else
        read(20,*) temp2, temp, (rmin_poly(in), in=1,ndof) 
        read(20,*) temp2, temp, (dr_poly(in), in=1,ndof)         
      end if
        
        
      close(20)
  
     ! necessary?
      nint2e = int(nint2e)
      nrprime = int(nrprime)
      nrprimn = int(nrprimn)
      nrensp = int(nrensp)
  
      ! not using this approach, to distingush spatial (norb/nrorb) and spinorbital (_spinorbital)
      ! it can be used to reduce passing arguments, but not sure if it is conceptual readable
      ! in open-shell, _spinorbital variables are used for unrestricted approach
      ! similar to psi4, rhf, uhf in different files
      !if (do_openshell) then
      !  norb = 2*nrorb
      !else
      !norb = nrorb
        
      ! use a series of indices, odd alpha even better, as the closed-shell convention
      ! then separate into alpha/beta entries
      nrorb_alpha = nrorb  
      nrorb_beta  = nrorb  
      nrorb_spinorbital = nrorb_alpha + nrorb_beta
        
        
      !norb_alpha  = nrorb_alpha
      !norb_beta   = nrorb_beta
      !norb_spinorbital = nrorb_alpha + nrorb_beta
      !endif
  
      !> Automatically determine a max energy value to shift to from the scfvals input
      if (autoshifte) then
        allocate(scf_temp(nrprimn))
        open(newunit=scfvals, file=trim(scfv_path), status='old')
        do i=1, nrprimn
          read(scfvals,*) scf_temp(i)
        enddo
  !      shifte = dble(ceiling(abs(minval(scf_temp))))
        shifte = dble(nint(abs(minval(scf_temp))))
  
        close(scfvals)
        deallocate(scf_temp)
      else
        shifte = shifte
      endif
  
  
      !> Assemble remaining global input constants from the input
  
      !> there was definitely an issue here...
      !> forget the mass of the neutrons! Awesome, well this is the right fix
      ! need to be extended to polyatomics
      m1 = N1*proelec_mratio + A1N*neutron_elec_mratio
      m2 = N2*proelec_mratio + A2N*neutron_elec_mratio
  
      ! use nint and not int, because any gfortran version and possibly older
      ! intel versions cannot get the precision right for the gamma values and int always rounds down.
      
  
      !cw: contemperory check for closed-shell
      !!nrorb_alpha = nrorb
      !nrorb_beta  = nrorb
      !nrorb_spinorbital    = nrorb_alpha + nrorb_beta    
  
      !nsz = 2
      nel_alpha = (nel + nsz)/2
      nel_beta  = (nel - nsz)/2     
      nel_spinorbital  = nel   
      
      write (*,*) 'check input', nel_alpha, nel_beta, nel, nsz, nel_spinorbital
      
      if ( nel_alpha < 0 ) then
        write (*,*) 'negative number of alpha-spin electrons, calculation stops'
        stop
      end if
      
      if ( nel_beta < 0 ) then
        write (*,*) 'negative number of beta-spin electrons, calculation stops'
        stop
      end if        
      
      if ( nel < 0 ) then
        write (*,*) 'negative number of electrons, calculation stops'
        stop
      end if    
      
      if ( flag_spinorbital == 0 .and.  nel > 2*nrorb ) then 
        write (*,*) 'number of electrons exceeds number of orbitals *2 in closed-shell calculations, calculation stop'
        stop
      end if
      
      
      if ( flag_spinorbital == 1 .and.  nel_alpha > nrorb_alpha ) then 
        write (*,*) 'number of alpha-spin electrons exceeds number of alpha-spin orbitals in open-shell & 
  & unrestricted calculations, calculation stop'
        stop
      end if    
      
      if ( flag_spinorbital == 1 .and.  nel_beta > nrorb_beta ) then 
        write (*,*) 'number of beta-spin electrons exceeds number of beta-spin orbitals in open-shell & 
& calculations, calculation stop'
        stop
      end if    
      
  
   
      if (flag_spinorbital == 1) then
        write (*,*) "Open-shell calculation"
        write (*,*) 'multiplicity', mult, 'nel', nel  
     
     !nsz is the number of parallel electrons
        if ( mod(nel,2)==0 .and. mod(nsz,2)==1 ) then
          write (*,*) 'even number of electrons and odd number of singly occupied parallel electrons, inconsistent'
          stop
        end if  
 
        if ( mod(nel,2)==1 .and. mod(nsz,2)==0 ) then
          write (*,*) 'odd number of electrons and even number of singly occupied parallel electrons, inconsistent'
          stop
        end if      
     
     
      else if (flag_spinorbital == 0) then
        write (*,*) 'Using closed-shell modules'
        nsz = 0
      else  
        write (*,*) 'not supported for flag_spinorbital not 0 nor 1'
        stop
      end if 
         
      
      !nel_spinorbital = nel
      !nevertheless, keep it for a while
      !the reason to adopt _spinorbital, is, unrestricted and restricted methods are different
      ! in psi4, rhf and uhf are in different files
      ! https://github.com/psi4/psi4/blob/master/psi4/src/psi4/libscf_solver/rhf.cc
      ! https://github.com/psi4/psi4/blob/master/psi4/src/psi4/libscf_solver/uhf.cc
      ! in MCEND, odd and even number represent alpha and beta electrons
      ! not sure if this is the optimal approach
      ! it is possible to either more merge or more separate
      ! leave with it for a while
      nrindep = nint(combinatorial(2*nrorb, nel))
      nrindep_spinorbital = nint(combinatorial(nrorb_spinorbital, nel_spinorbital))
      nrindep_frzorb = nint(combinatorial(2*nrorb - 2*nrorb_fc - 2*nrorb_fv, nel- 2*nrorb_fc))   
      nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital - nrorb_fv_spinorbital, &
& nel - nrorb_fc_spinorbital))   

! active frozen core
      nrindep = nrindep_frzorb
      nrindep_spinorbital = nrindep_frzorb_spinorbital
      
      
      write (*,*) 'nrindep_spinorbital', nrindep_spinorbital
      
      !> ISOMASS
      massn = (m2*m1)/(m2 + m1)
      !            A             phie               phin
      dgldim = nrindep*nrspf + nrorb*nrprime + nrspf*nrprimn
      dgldim_spinorbital = nrindep_spinorbital*nrspf + nrorb_spinorbital*nrprime + nrspf*nrprimn
      
  !    dgldim = nrindep*nrspf + nrorb*nrprime + nrspf*nrprimn
      hmfsize = hmfscale*nrindep
  
  !    max_nrindep = 14*nrindep
  
  
      !nrorb_spinorbital = nrorb * 2
      !nrorb_alpha + nrorb_beta
      
      write (*,*) 'check nr orb, alpha, beta, spinorbital', nrorb, nrorb_alpha, nrorb_beta !, nrorb_spinorbital
      write (*,*) 'check nel, alpha, beta', nel, nel_alpha, nel_beta
      
      
      max_nrindep = nint(combinatorial(2*nrorb, nel-1))
      ! fixme: need to take care h2+ case, when nel=1
      ! seems get 0
      max_nrindep_2 = nint(combinatorial(2*nrorb, nel-2))
      
      
    !  write (*,*) 'test section', nint(combinatorial(2, 0))
    !  write (*,*) 'test section', nint(combinatorial(2, -1))
      
    !  max_nrindep_alpha = nint(combinatorial(nrorb_alpha, nel_alpha-1))
    !  max_nrindep_beta = nint(combinatorial(nrorb_beta, nel_beta - 1))
    !  max_nrindep_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-1))   
      max_nrindep_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-1)) 
      max_nrindep_2_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-2)) 
      
      
      write (*,*) 'max_nrindep_2_spinorbital', max_nrindep_2_spinorbital
      
      if (nrorb_fc == 0) then
        max_nrindep_frzorb = nint(combinatorial(2*nrorb-2*nrorb_fv, nel-1)) 
        max_nrindep_2_frzorb = nint(combinatorial(2*nrorb-2*nrorb_fv, nel-2)) 
      else
      ! if 2e excited, no way to combine 1e back, the substraction part corresponds to the entries of this possibility
        max_nrindep_frzorb = nint(combinatorial(2*nrorb-2*nrorb_fv, nel-1) &
        &- combinatorial(2*nrorb - 2*nrorb_fc -  2*nrorb_fv, nel-1) )    
        max_nrindep_2_frzorb = nint(combinatorial(2*nrorb-2*nrorb_fv, nel-2)  &
        &- combinatorial(2*nrorb - 2*nrorb_fc -  2*nrorb_fv - 1, nel-2) )    
      end if
      
      
      if (nrorb_fc_spinorbital == 0) then
        max_nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital - nrorb_fv_spinorbital, nel-1))
        max_nrindep_2_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital - nrorb_fv_spinorbital, nel-2))
      else
        max_nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital - nrorb_fv_spinorbital, nel-1) &
        &- combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital - nrorb_fv_spinorbital, nel-1))
        
        max_nrindep_2_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital-nrorb_fv_spinorbital, nel-2) &
        &- combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital - nrorb_fv_spinorbital - 1, nel-2))        
      end if
      
      
      
      ! frzorb is for correlation, not orbital optimization
      ! active frozen core
      max_nrindep = max_nrindep_frzorb
      max_nrindep_spinorbital = max_nrindep_frzorb_spinorbital
  
      !> ok, this will give the exact value we need for the shdl array
      !> nrshf is set to the same thing, though I suppose it need not be so
      !> as it probably is bad practice to have two variables that represent
      !> the same thing.
      !
      ! cw: nrshf is the number of single-hole function (shf)
      ! shf is <i|psi>, 
      ! I think the nrshf is the number of orbitals. 
      ! In principle, nrshf != max_nrindep/ size of shdl array
      ! In practice, it is used as the size of shdl array
      
      nrshf = max_nrindep
      !nrshf_alpha = max_nrindep_alpha
      !nrshf_beta  = max_nrindep_beta
  !    nrshf_spinorbital  = max_nrindep_spinorbital
      nrshf_spinorbital  = max_nrindep_spinorbital
      
      nrshf_frzorb = max_nrindep_frzorb
      nrshf_frzorb_spinorbital = max_nrindep_frzorb_spinorbital
     ! nrshf_spinorbital = nrshf_frzorb_spinorbital
     ! nrshf = nrshf_frzorb
      
      nrdhf       = max_nrindep_2
      nrdhf_spinorbital  = max_nrindep_2_spinorbital  
      
      !write (*,*) 'check nrshf, alpha, beta', nrshf, nrshf_alpha, nrshf_beta
      
      
      write(*,'((a40,1x),(a7))') 'Compound: ', trim(cmpdname)
      write(*,'((a40,1x),(a60))') 'Integral Folder: ', trim(intdir)
      write(*,'((a40,1x),(a60))') 'Initial Hmat values read from: ', trim(scfv_path)
      write(*,'((a40,1x),(i7))') 'No. electronic bfn: ', nrprime
      write(*,'((a40,1x),(i7))') 'No. electrons: ', nel
      write(*,'((a40,1x),(i7))') 'No. molecular orbitals: ', 2*nrorb
      ! experimental -> not functioning
      if (do_openshell) then
        write(*,'((a40,1x),(i7))') 'No. alpha orbitals: ', nrorb_alpha
        write(*,'((a40,1x),(i7))') 'No. beta orbitals: ', nrorb_beta
      else
  !      write(*,'((a40,1x),(i7))') 'No. spatial orbitals: ', nrorb
        write(*,'((a40,1x),(i7))') 'No. spatial orbitals: ', nrorb
      endif
      write(*,'((a40,1x),(i7))') 'No. nuclear grid bfn.: ', nrprimn
      write(*,'((a40,1x),(i7))') 'Nuclear SPFs: ', nrspf
      write(*,'((a40,1x),(i7))') 'Nuclear sampling points: ', nrensp
      write(*,'((a40,1x),(i7))') 'Indepedent A elements: ', nrindep
      write(*,'((a40,1x),(i7))') 'Non-zero two-electron ints: ', nint2e
      write(*,'((a40,1x),(f15.7))')  'Initial grid point: ', rmin
      write(*,'((a40,1x),(f15.7))') 'Grid step size: ', dr
      write(*,'((a40,1x),(f15.3))') 'Energy shifted by: ', shifte
      write(*,'((a40,1x),(l7))') 'logev: ', logev
      write(*,'((a40,1x),(l7))') 'Frozen nuclei: ', freeze_nuc
      write(*,'((a40,1x),(f15.3))') 'Mass 1: ', m1
      write(*,'((a40,1x),(f15.3))') 'Mass 2: ', m2
  
      if (wsign < 0) then
        write(*,'((a40,1x),(f5.1),(f5.1)"i")') 'Imaginary time propagation, t: ', imagt
      else if (wsign > 0) then
        write(*,'((a40,1x),(f5.1),(f5.1)"i")') 'Real time propagation, t: ', imagt
        if (aucofu == 0) then
          write(*,*) 'No autocorrelation function used'
        else
          write(*,*) 'Autocorrelation function in use, parameter given: ', aucofu
        endif
  
        write(*,*) ''
        write(*,*) '                       Laser Parameters'
        write(*,*) '______________________________________________________________'
        write(*,*)
        write(*,'((a20,1x),(i7))') 'Number of pulses: ', nrpulses
        do j=1, nrpulses
          write(*,'((a30,1x),(i7))') ' Pulse ', j
          write(*,'((a30,1x),3(f4.2,1x))') 'Pulse direction (x,y,z): ', (lpolar(i,j),i=1,3)
          write(*,'((a30,1x),(f15.7))') 'Pulse height: ', lph(j)
          write(*,'((a30,1x),(f15.7))') 'Pulse width: ', lpw0(j)
          write(*,'((a30,1x),(f15.7))') 'Laser frequency: ', lfreq(j)
          write(*,*) ''
        enddo
        write(*,'((a40,1x),(f15.2))') 'Total Propa Time (fs) = ', tfin0
        write(*,'((a40,1x),(f15.2))') 'Time step (fs) = ', tout0
      endif
      if (use_wcap) then
        write(*,*) '****************** Warning CAP in use! ***********************'
        write(*,*) "The CAP implementation is not yet finalized"
        write(*,*) "Consistent results across Intel/gfortran versions are also totally guaranteed!!!!!"
        write(*,*) "Exercise use with caution!"
        write(*,*) ""
        write(*,*) 'CAP parameters are as follows:'
        write(*,'((a40,1x),(f11.5))') 'wgrid: ', wgrid
        write(*,'((a40,1x),(es23.16))') 'wstrength: ', etan
        write(*,'((a40,1x),(i3))') 'worder: ', bn
      endif
  
      deallocate(charvars,inpvars,keywargs,keywargs2)
      return
  
     end subroutine read_file
  
      ! function that computes the combinatorial, uses the gamma
      ! intrinsic instead of some self defined factorial
      ! as the factorial function messes up for higher values of n
      pure real(dp) function combinatorial(n, r)
        implicit none
        integer, intent(in) :: n, r
        real(dp) :: nr, rr
        real(dp) :: cnr
        nr = real(n)
        rr = real(r)
  
        if (0.0_dp <= rr .and. rr <= nr) then
          cnr = gamma(nr + 1.0_dp)/(gamma(rr + 1.0_dp)*gamma(nr - rr + 1.0_dp))
        else
          cnr = 0.0_dp
        end if
        combinatorial = cnr
  
      end function combinatorial
  
      ! convert a string into an integer, quite simple actually
      elemental subroutine str2int(str,intval)
        character(len=*), intent(in) :: str
        integer, intent(out)         :: intval
  
        read(str,*) intval
  
      end subroutine str2int
  
      !> new functions to get the integer string, should be better than int2str
      pure integer function str_ilen(i) result(sz)
        ! Returns the length of the string representation of 'i'
        !> gives the length of an integer string
        integer, intent(in) :: i
        integer, parameter :: str_max = 100
        character(str_max) :: s
        !> if you give a string that has 100 characters, it's not gonna happen
        write(s, '(i0)') i
        sz = len_trim(s)
      end function
  
      pure function stri(i) result(s)
  !    elemental function int2str(s, i)
        ! Converts integer "i" to string
        integer, intent(in) :: i
        character(len=str_ilen(i)) :: s
  
        write(s,'(i0)') i
  
      end function
  
  !    end subroutine
  
      ! [older] convert an integer into a string
      ! much more difficult to do in fortran, this routine
      !>
      elemental subroutine int2str(str,intval)
  
        implicit none
        integer, intent(in)           :: intval
        character(len=*), intent(out) :: str
        integer :: i, smin, smax
        character(len=7) :: int_format
        integer :: max_int_len
        character(len=1) :: int_str_list(9)
        int_str_list = ['1','2','3','4','5','6','7','8','9']
  
        ! go up to integers of size 1e8
        max_int_len = 9
  
        ! slightly more 'elegant way' to write this
        ! allows a more easy way to set the maximum value of iterations
        ! requires a manual integer string list however
        if (intval == 0) then
          write(str,'(i1)') intval
        else
          do i=1, max_int_len
            smin = int(10**(i-1))
            smax = int(10**(i))
            if (smin <= intval .and. intval < smax) then
              int_format = '(i'//int_str_list(i)//')'
              write(str,trim(int_format)) intval
              exit
            endif
          enddo
  
        endif
  
      end subroutine int2str
  
  !    subroutine read_input_file(filename, nrows, inputio)
  !      ! this opens a file reads the number of lines and then rewinds it
  !
  !      character (len=*), intent(in) :: filename
  !      integer, intent(inout) :: inputio
  !      integer :: ios, ios2, input_len
  !      character(len=1) temp
  !      character(len=20) temp2
  !      character(len=300) iom
  !
  !      ios = 0
  !      ios2 = 0
  !      open(newunit=inputio, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iom)
  !      if (ios /= 0) then
  !        write(*,*) 'Fatal error!!!', trim(iom)
  !        stop
  !      endif
  !
  !      input_len = 0
  !      if (ios == 0) then
  !        do while (ios2 == 0)
  !          if (ios2 == 0) then
  !            read(inputio,*, iostat=ios2) temp
  !            input_len = input_len + 1
  !          end if
  !        enddo
  !      endif
  !      input_len = input_len - 1
  !      rewind(inputio)
  !
  !      nrows = input_len
  !
  !    end subroutine
  

! from stackoverflow
! https://stackoverflow.com/questions/30006834/reading-a-file-of-lists-of-integers-in-fortran/30007516?noredirect=1#comment119141841_30007516
! to split string into integers
    function string_to_integers(str, sep) result(a)
      character(6), allocatable :: a(:)
      character(*) :: str
      character :: sep
      integer :: i, n_sep
    
      n_sep = 0
      
    !  write (*,*) 'aa', sep, 'bb'
    !  write (*,*) 'str', str, 'len', len(str)    
    
      do i = 1, len_trim(str)
        if (str(i:i)==sep) then
          n_sep = n_sep + 1
          str(i:i) = ','
         end if
      end do
      
    !  write (*,*) 'check', n_sep, str
      
      allocate(a(n_sep+1))
      read(str,*) a
    end function  
  
  
    
  
  
  end module inputvars
  
  ! unused
  !    pure integer function factorial(n)
  !      implicit none
  !      integer, intent(in) :: n
  !      integer :: fact, i
  !
  !      fact = 1
  !      do i = 1, n
  !        fact = fact*i
  !      enddo
  !      factorial = fact
  !    end function factorial
  !        cnr = factorial(n)/(factorial(r)*factorial(n - r))
  
  
