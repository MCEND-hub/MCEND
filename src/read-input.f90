!> This is an primitive fortran input reader, has flexibility with the ordering,
!! and some values have defaults, so they need not be specified. Comments are technically
!! supported, any keywords that aren't recognized won't match anything and thus nothing will be done with them.
!! This module contains core parameters and should not need to be changed.
!! @param i64, i32 64 and 32 bit integers
!! @param dp explicit definition of double precision
!! @param qp quadruple precision, if ever needed
!! @param c0, cr, ci, c1 some complex numbers
!! @param pi the number pi
!! @param au2fs atomic time conversion factor to fs
  module params
    integer, parameter :: i64 = selected_int_kind(15)
    integer, parameter :: i32 = selected_int_kind(6)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)
    complex(dp), parameter :: c0 = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cr = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: c1 = (1.0_dp, 1.0_dp)
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: thresh_time_zero = 1.0d-6
    real(dp), parameter :: thresh_parameter_zero = 1.0d-6
    real(dp), parameter :: au2fs = 0.0241888_dp
    real(dp), parameter :: proelec_mratio = 1836.15267389_dp
    real(dp), parameter :: neutron_elec_mratio = 1838.683661589_dp
  end module params


!> This module declares all input variables needed to start MCEND;
!! it also contains subroutines that read in these variables from
!! an input file.
!! @param dt0 initial step-size (fs)
!! @param facdec, facinc adaptive step-size increase/decrease limits
!! @param beta further adaptive step-size parameter
!! @param etacap Eta Cap parameter, electrons
!! @param naocap nr of first AO with CAP, electrons
!! @param epsreg0 density matrix regularization parameter
!! @param tol tolerance for adaptive step-size control
!! @param nrfrorb number of frozen core orbitals
!! @param shifte energy shift of electronic one-particle Hamiltonian
!! @param wgrid internuclear distance at which the nuclear CAP starts
!! @param bn order of the CAP - usually 2 or 3
!! @param etan strength of the CAP - usually 2 or 3
!! @param wsign sign of the CAP - has to be switched for imagt
!! @param nrprime no of electronic basis functions
!! @param nel no of electrons
!! @param nrprimn no of nuclear grid points
!! @param nrspf number of single particle functions
!! @param nrorb no of orbitals
!! @param nrorb_fc frozen core orbitals, same as nrfrorb
!! @param nrensp number of explicit nuclear sampling points
!! @param nrindep nr of independent A-vector elements
!! @param inigrid starting point on nuclear grid
!! @param mul multiplicity - not yet used
!! @param basis use precomputed basis 
!! basis = 1 precomputed basis at nrensp positions
!! @param restart read in initial state? (2 => restart, 0 => no restart, 1 => generate additional SPF)
!! @param aucofu autocorrelation function? (0=no, 1=yes, 2=after lpw, please add the additional options that are missing here)
!! @param numthreads how many openmp threads to use
!! @param hmfsize mean field hash table size
!! @param dr, rmin spacing of nuclear grid,  and smallest grid point
!! @param tfin0 total propagation time (fs)
!! @param t_initial if this is a restart this can be used to pick up where the laser left off
!! @param tout0 output interval (fs) - this significantly slows down the program when set too low
!! @param lpw0 laser pulse parameters, sin**2*cos shape; pulse width/duration (fs)
!! @param lph pulse height (Eh/ea0)
!! @param lpolar pulse polarization direction(s)
!! @param lfreq photon energy (hartree)
!! @param m1, m2 atomic masses
!! @param N1, N2 nuclear charge
!! @param massn reduced nuclear mass
!! @param ndof number of quantum degree of freedom
!! @param nint2e number of non-zero two-electron integrals
!! @param hmfscale scales hmfsize
!! @param cmpdname name of the compound
!! @param intdir folder with the electronic integrals
!! @param rspfile file with the nuclear grid points
!! @param initguessf initial guess file
!! @param scfv_path path to the hmat file with values generated from orthog routine
!! @param n2int_path path to the two electron integrals
!! @param mcend_top path to the current execution directory
!! @param integral_library_env variable that points to the path of the integral library
!! @param logwavef print wave function every tout?
!! @param logev diagonalize the CI matrix to get energy eigenstates
!! @param imagt imaginary time propagation (imagt=(1,0) => normal propagation) (imagt=(0,-1) => relaxation)
!! @param freeze_nuc freeze nuclear orbitals, may not be functional
!! @param read_stand read integrals in the standard (old) way
!! @param use_wcap use the nuclear CAP
!! @param autoshifte shift the energy towards zero automatically
!! @param do_openshell control open-shell unrestricted or closed-shell restricted methods
!! @param mf_reduction move one particle operator, so far Te, outside meanfields
!! @param nrpulses number of pulses, determined from pulse params
!! @param nsz number of parallel electrons
!! @param flag_spinorbital similar to do_openshell, used in previous version; may be extended beyond 0/1 for configuration state function (CSF) approach
!! @param flag_fc controls frozen core
!! @param print_level control print level for intermediate wall time
module inputvars
  use params
    implicit none

      real(dp), parameter :: dt0 = 1.0d-5
      real(dp), parameter :: facdec = 0.333333_dp
      real(dp), parameter :: facinc = 6.0_dp
      real(dp), parameter :: beta = 0.04_dp
      real(dp), parameter :: etacap = 0.0_dp
      integer,  parameter :: naocap = 26
      real(dp), parameter :: epsreg0 = 1.0d-10
      real(dp), parameter :: tol = 1.0d-10
      integer :: nrfrorb
      real(dp) :: shifte
      real(dp) :: wgrid
      integer  :: bn
      real(dp) :: etan
      integer  :: wsign
      integer :: nrprime, nel
      integer :: nel_alpha, nel_beta
      integer :: nel_spinorbital
      integer :: nrprimn, nrspf
      integer :: nrorb
      integer :: nrorb_alpha, nrorb_beta
      integer :: nrorb_spinorbital
      integer :: nrorb_init_alpha, nrorb_init_beta
      integer :: nrorb_fc, nrorb_fc_spinorbital
      integer :: nrensp
      integer :: nrindep
      integer :: nrindep_spinorbital
      integer :: nrindep_frzorb, nrindep_frzorb_spinorbital
      integer :: inigrd
      integer :: mult
      integer :: basis
      integer :: restart
      integer :: aucofu
      integer :: max_nrindep
      integer :: max_nrindep_spinorbital
      integer :: max_nrindep_frzorb
      integer :: max_nrindep_frzorb_spinorbital
      integer :: max_nrindep_2
      integer :: max_nrindep_2_spinorbital
      integer :: max_nrindep_2_frzorb, max_nrindep_2_frzorb_spinorbital
      integer :: nrshf
      integer :: nrshf_spinorbital
      integer :: nrshf_frzorb,nrshf_frzorb_spinorbital
      integer :: nrdhf, nrdhf_spinorbital
      integer :: nrdhf_frzorb,  nrdhf_frzorb_spinorbital
      integer :: dgldim
      integer :: dgldim_spinorbital
      integer(i32) :: numthreads
      integer :: hmfsize
      real(dp) :: dr, rmin
      real(dp) :: tfin0
      real(dp) :: t_initial
      real(dp) :: tout0
      real(dp), allocatable :: lpw0(:)
      real(dp), allocatable :: lph(:)
      real(dp), allocatable :: lfreq(:)
      real(dp), allocatable :: lpolar(:,:)
      real(dp) :: m2, m1
      real(dp) :: N2, N1
      real(dp) :: A2N, A1N
      real(dp) :: massn
      integer :: ndof
      integer :: nint2e
      integer :: hmfscale
      character(10) :: cmpdname
      character(300) :: intdir
      character(300) :: rspfile
      character(300) :: initguessf
      character(300) :: scfv_path
      character(300) :: n2int_path
      character(300) :: mcend_top
      character(300) :: integral_library_env
      character(300) :: path2ints
      logical :: logwavef
      logical :: logev
      complex(dp) :: imagt
      logical :: freeze_nuc
      logical :: read_stand
      logical :: use_wcap
      logical :: autoshifte
      logical :: do_openshell
      logical :: mf_reduction
      integer :: nrpulses
      integer ::  nsz
      integer :: flag_spinorbital
      integer :: flag_fc
      integer :: print_level

    contains

    !> This subroutine reads the input file
    subroutine read_file(filename)
      implicit none

      character(len=*), intent(in) :: filename
      integer :: inputio, ios, ios2
      character(1) :: temp
      character(20) :: temp2
      character(300) :: iom
      integer :: i, j, k
      integer :: input_len, nrows
      integer :: nrownum, nrowchar, nrownum2
      integer :: scfvals
      real(dp), allocatable :: scf_temp(:)
      character(1) :: line
      character(255) :: line_temp
      character(100), allocatable :: charvars(:)
      character(100), allocatable :: linevars(:)
      character(50), allocatable :: keywargs(:)
      character(50), allocatable :: keywargs2(:)
      real(dp), allocatable :: inpvars(:)
      integer :: natom

      call get_environment_variable("PWD", mcend_top)
      call get_environment_variable("MCEND_BASIS_LIBRARY", integral_library_env)

      ! get the number of lines in a file w/o using a external bash call
      ios = 0
      ios2 = 0

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
        endif
      enddo

      nrowchar = nrows - nrownum

      allocate(inpvars(nrownum))
      allocate(charvars(nrowchar))
      allocate(linevars(nrownum))

      rewind(inputio)

      ! load in integers and reals
      do j=1, nrownum - 1
  ! reading format, keywords, '=' for temp, and variables
        read(inputio, *) keywargs(j), temp, inpvars(j)
      enddo

      rewind(inputio)

      ! linevars to contain more entries in N
      do j=1, nrownum - 1
        read(inputio,"(A80)") linevars(j)
      enddo

      read(inputio,*) keywargs(nrownum)

      ! read in character variables
      k = 1
      do j = nrownum+1, nrownum2-1
        read(inputio,'(A)') line_temp !keywargs(j), temp, charvars(k) 
        call read_entries(line_temp, keywargs(j), temp, charvars(k)) 
        k = k + 1
      enddo

      nrpulses = nrows - nrownum2 - 1

      allocate(lpolar(3,nrpulses))
      allocate(lph(nrpulses))
      allocate(lpw0(nrpulses))
      allocate(lfreq(nrpulses))

      lpolar(:,:) = 0.0_dp

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
      ! default parameters
      flag_spinorbital = 0
      nsz = 0
      print_level = 0

     ! write (*,*) 'nrownum', nrownum

      do i=1, nrownum - 1
        ! more often read from basis set information
        if (trim(keywargs(i)) == 'nrprime') then
          nrprime = int(inpvars(i))
          write (*,*) 'nrprime in read-input', nrprime
        elseif (trim(keywargs(i)) == 'nel') then
          nel = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'multiplicity') then
          mult= int(inpvars(i))
          nsz = mult - 1
        elseif (trim(keywargs(i)) == 'natom') then
          natom = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrprimn') then
          nrprimn = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrspf') then
          nrspf = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrorb') then
          nrorb = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'nrorb_fc') then
          nrorb_fc = int(inpvars(i))
          ! current convention
          nrorb_fc_spinorbital = nrorb_fc * 2
        elseif (trim(keywargs(i)) == 'nrensp') then
          nrensp = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'dr') then
          dr = inpvars(i)
        elseif (trim(keywargs(i)) == 'rmin') then
          rmin = inpvars(i)
        elseif (trim(keywargs(i)) == 'shifte') then
          shifte = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'inigrd') then
          inigrd = int(inpvars(i))
        ! WARNING: Nuclear CAP parameters
        elseif (trim(keywargs(i)) == 'wgrid') then
          wgrid = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'wstrength') then
          etan = dble(inpvars(i))
        elseif (trim(keywargs(i)) == 'worder') then
          bn = int(inpvars(i))
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
        elseif (trim(keywargs(i)) == 'N1') then
          N1 = inpvars(i)
        elseif (trim(keywargs(i)) == 'N2') then
          N2 = inpvars(i)
        !> specify the neutron mass
        elseif (trim(keywargs(i)) == 'A1N') then
          A1N = inpvars(i)
        elseif (trim(keywargs(i)) == 'A2N') then
          A2N = inpvars(i)
        !> [not yet fully working] Number of frozen orbitals
        elseif (trim(keywargs(i)) == 'nrfrorb') then
          nrfrorb = int(inpvars(i))
          !> multiplicity
        elseif (trim(keywargs(i)) == 'mult') then
          mult = int(inpvars(i))
          !> scaling for hmf,
        elseif (trim(keywargs(i)) == 'hmfscale') then
          hmfscale = int(inpvars(i))
        elseif (trim(keywargs(i)) == 'print_level') then
          !> 0: fewer print; 1: print time for intermediate steps
          print_level = int(inpvars(i))
         ! write (*,*) 'print_level', print_level
        endif


      enddo

      !> crude fail safe for isotope assignment, so in case neutrons number is not
      !> given, assume num protons = number of neutrons, unless the proton count is
      !> 1 which is hydrogen and won't have a neutron in the most common isotope
      !> eventually should put this into an array form
      if ( dabs(A1N - 0.0_dp) < thresh_parameter_zero .and. dabs(N1 - 1.0_dp) > thresh_parameter_zero ) A1N = N1
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
  !      elseif (trim(kewargs(j)) == 'basis_path') then
  !        basis_path_dir = charvars(i)
        elseif (trim(keywargs(j)) == 'mf_reduction') then
          if (trim(charvars(i)) == 'yes') then
            mf_reduction = .true.
          else
            mf_reduction = .false.
          endif
        endif
      enddo


      !> get current working directory
      if ( intdir(1:1)  == '.' .or.  intdir(1:1) == '/' ) then  
        intdir = trim(intdir)
        write(*,*) "Relative or absolute path to integral library detected", intdir
      else if ( len(trim(integral_library_env)) >= 1 ) then
        write(*,*) "Environment variable for integral library detected", integral_library_env
        intdir = trim(integral_library_env)//'/'//trim(intdir)
      else
        intdir = 'basis_library/'//trim(intdir)
        write(*,*) "Using MCEND legacy basis_library folder that needs to be &
         & present in the current directory", intdir
      endif




      scfv_path = trim(intdir)//'/scf-guess-'//trim(cmpdname)//'.dat'
      n2int_path = trim(intdir)//'/info-'//trim(cmpdname)//'.dat'
      rspfile = trim(intdir)//'/rsp_'//trim(cmpdname)

      !> open info file to get these parameters
      open(20,file=trim(n2int_path))
      read(20,*) temp2, temp, ndof
      if (temp2 == 'ndof') then
        ndof = ndof
        read(20,*) temp2, temp, nint2e
      !adopt the previous version, only for diatomics
      ! if I evolve all format into polyatomic, the code could be more concise
      ! not sure if this is a good practice
      else if (temp2 == 'nint2e') then
        nint2e = ndof
        ndof   = 1
      else
        write (*,*) 'not ndof/nint2e in the first line, not supported'
        stop
      end if

    !  write (*,*) 'dof check', ndof


      ! so far, the precomputed basis for polyatomic is merged so to speak
      ! a single nrprime, for combinational positions
      read(20,*) temp2, temp, nrprime

      if (ndof > 1) then
        write (*,*) 'more than 1 dof, calculation stops'
        stop
      end if


      read(20,*) temp2, temp, nrprimn
      read(20,*) temp2, temp, nrensp

      read(20,*) temp2, temp, rmin
      read(20,*) temp2, temp, dr


      close(20)

      nint2e = int(nint2e)
      nrprime = int(nrprime)
      nrprimn = int(nrprimn)
      nrensp = int(nrensp)

      ! use a series of indices, odd alpha even better, as the closed-shell convention
      ! then separate into alpha/beta entries
      nrorb_alpha = nrorb
      nrorb_beta  = nrorb
      nrorb_spinorbital = nrorb_alpha + nrorb_beta

      !> Automatically determine a max energy value to shift to from the scfvals input
      if (autoshifte) then
        allocate(scf_temp(nrprimn))
        open(newunit=scfvals, file=trim(scfv_path), status='old')
        do i=1, nrprimn
          read(scfvals,*) scf_temp(i)
        enddo
        shifte = dble(nint(abs(minval(scf_temp))))

        close(scfvals)
        deallocate(scf_temp)
      else
        shifte = shifte
      endif

      !> Assemble remaining global input constants from the input
      m1 = N1*proelec_mratio + A1N*neutron_elec_mratio
      m2 = N2*proelec_mratio + A2N*neutron_elec_mratio


      nel_alpha = (nel + nsz)/2
      nel_beta  = (nel - nsz)/2
      nel_spinorbital  = nel

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

      ! in MCEND, odd and even number represent alpha and beta electrons
      nrindep = nint(combinatorial(2*nrorb, nel))
      nrindep_spinorbital = nint(combinatorial(nrorb_spinorbital, nel_spinorbital))
      nrindep_frzorb = nint(combinatorial(2*nrorb - 2*nrorb_fc, nel- 2*nrorb_fc))
      nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital, &
& nel - nrorb_fc_spinorbital))

! active frozen core
      nrindep = nrindep_frzorb
      nrindep_spinorbital = nrindep_frzorb_spinorbital

      !> ISOMASS
      massn = (m2*m1)/(m2 + m1)


      dgldim = nrindep*nrspf + nrorb*nrprime + nrspf*nrprimn
      dgldim_spinorbital = nrindep_spinorbital*nrspf + nrorb_spinorbital*nrprime + nrspf*nrprimn
      hmfsize = hmfscale*nrindep

      max_nrindep = nint(combinatorial(2*nrorb, nel-1))
      max_nrindep_2 = nint(combinatorial(2*nrorb, nel-2))
      max_nrindep_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-1))
      max_nrindep_2_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-2))

      if (nrorb_fc == 0) then
        max_nrindep_frzorb = nint(combinatorial(2*nrorb, nel-1))
        max_nrindep_2_frzorb = nint(combinatorial(2*nrorb, nel-2))
      else
      ! if 2e excited, no way to combine 1e back, the substraction part corresponds to the entries of this possibility
        max_nrindep_frzorb = nint(combinatorial(2*nrorb, nel-1) &
        &- combinatorial(2*nrorb - 2*nrorb_fc, nel-1) )
        max_nrindep_2_frzorb = nint(combinatorial(2*nrorb, nel-2)  &
        &- combinatorial(2*nrorb - 2*nrorb_fc - 1, nel-2) )
      end if

      if (nrorb_fc_spinorbital == 0) then
        max_nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-1))
        max_nrindep_2_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-2))
      else
        max_nrindep_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-1) &
        &- combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital, nel-1))

        max_nrindep_2_frzorb_spinorbital = nint(combinatorial(nrorb_spinorbital, nel-2) &
        &- combinatorial(nrorb_spinorbital - nrorb_fc_spinorbital - 1, nel-2))
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
      nrshf_spinorbital = max_nrindep_spinorbital
      nrshf_frzorb = max_nrindep_frzorb
      nrshf_frzorb_spinorbital = max_nrindep_frzorb_spinorbital
      nrdhf = max_nrindep_2
      nrdhf_spinorbital  = max_nrindep_2_spinorbital

      write(*,'((a40,1x),(a7))') 'Compound: ', trim(cmpdname)
      write(*,'((a40,1x),(a60))') 'Integral Folder: ', trim(intdir)
      write(*,'((a40,1x),(a60))') 'Initial Hmat values read from: ', trim(scfv_path)
      write(*,'((a40,1x),(i7))') 'No. electronic bfn: ', nrprime
      write(*,'((a40,1x),(i7))') 'No. electrons: ', nel
      write(*,'((a40,1x),(i7))') 'No. molecular orbitals: ', 2*nrorb
      if (do_openshell) then
        write(*,'((a40,1x),(i7))') 'No. alpha orbitals: ', nrorb_alpha
        write(*,'((a40,1x),(i7))') 'No. beta orbitals: ', nrorb_beta
      else
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

      !> function that computes the combinatorial, uses the gamma
      !! intrinsic instead of some self defined factorial
      !! as the factorial function messes up for higher values of n
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
        ! Converts integer "i" to string
        integer, intent(in) :: i
        character(len=str_ilen(i)) :: s

        write(s,'(i0)') i

      end function

      ! [older] convert an integer into a string
      ! much more difficult to do in fortran, this routine
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





      subroutine read_entries(input_line, var_1, var_2, var_3)
        implicit none
        character(len=*), intent(in)  :: input_line
        character(len=*), intent(out)  :: var_1, var_2, var_3
        character(len=255) :: temp_line
        integer :: i
    
        i = index(input_line, ' ') 
        if (i > 0) then
            var_1 = input_line(1: i - 1)
            temp_line = input_line(i + 1: )
        end if    
    
        i = index(temp_line, ' ') 
        if (i > 0) then
            var_2 = temp_line(1: i - 1)
            var_3 = temp_line(i + 1: )
        end if    
    
     !   write (*,*) var_1, var_2, var_3
    
      end subroutine read_entries
    
    
    
    







! from stackoverflow
! https://stackoverflow.com/questions/30006834/reading-a-file-of-lists-of-integers-in-fortran/30007516?noredirect=1#comment119141841_30007516
! to split string into integers
    function string_to_integers(str, sep) result(a)
      character(6), allocatable :: a(:)
      character(*) :: str
      character :: sep
      integer :: i, n_sep

      n_sep = 0

      do i = 1, len_trim(str)
        if (str(i:i)==sep) then
          n_sep = n_sep + 1
          str(i:i) = ','
         end if
      end do

      allocate(a(n_sep+1))
      read(str,*) a
    end function

  end module inputvars

  !> @file
  !> @brief contains the input reader
  !! and setting of global constants.
