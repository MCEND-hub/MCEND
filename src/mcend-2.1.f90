!> Main driver of the MCEND program.
 program MCElectronNuclearDynamics

   use params
   use globalvars
   use inputvars

   use omp_lib
   use moreadin
   use utils
   use propa
   use analyse

   implicit none

   character(32) :: arg
   ! get the date and time
   real(dp)              :: wtime
   character(8)          :: date
   character(10)         :: time
   character(5)          :: zone
   character(len=255)    :: myname
   integer :: values(8)


   call get_command_argument(1, arg)
   if (len_trim(arg) == 0) then
     write(*,*) 'No input specified!! Terminating Program!!!!'
   endif

   call date_and_time(date, time, zone, values)

   ! get hostname
   call hostnm(myname)
   wtime = omp_get_wtime()
   write(*,'(" "78("*")" ")')
   write(*,*) ""
   write(*,*) ""
   write(*,*) ""
   write(*,*) "               _____                                                          "
   write(*,*) "           ___|    _|__    ______    ______    ____   _    _____              "
   write(*,*) "          |    \  /  | |  |   ___|  |   ___|  |    \ | |  |     \             "
   write(*,*) "          |     \/   | |  |   |__   |   ___|  |     \| |  |      \            "
   write(*,*) "          |__/\__/|__|_|  |______|  |______|  |__/\____|  |______/            "
   write(*,*) "              |_____|                                                         "
   write(*,*) ""
   write(*,*) ""
   write(*,'(" "78("*")" ")')
   write(*,*) ""
   write(*,*) ""
   write(*,*) "-------------------MultiConfigurationElectronNuclearDynamics------------------"
   write(*,*) "---------------------------------MCEND v2.3a----------------------------------"
   write(*,*) ""
   write(*,'(" "78("*")" ")')
   write(*,*) ""
   write(*,*) "--------Mathias Nest, Inga S. Ulusoy, Lucas E. Aebersold, Cong Wang-----------"
   write(*,*) ""
   write(*,'(" "78("*")" ")')
   write(*,'(" "78("*")" ")')
   write(*,*) ""
   write(*,'("  Executed on: "a30)') trim(myname)
   write(*,'("  Date: "20(" "),i4"."i2.2"."i2.2" at "i2.2":"i2.2":"i2.2)') values(1:3), values(5:7)

   ! read input file
   write(*,*) "Reading input parameters from file: ", arg
   call read_file(arg)
   write(*,*) "Input parameters read successfully"
   write(*,*) "Using ", numthreads, " threads for parallelization."
   call omp_set_num_threads(numthreads)
  !uncomment the below line hwn using mkl=parallel
  !call mkl_set_num_threads(numthreads)

   write(*,'(a40,1x,a65)') 'MCEND location: ', mcend_top
   write(*,'(a40,1x,a32)') " Input file: ", arg
   write(*,'(a40,i8)') "The number of processors available: ", omp_get_num_procs()
   write(*,'(a40,i8)') "The number of threads called: ", omp_get_max_threads()
   write(*,*) ""

   call allocate_arrays()

   ! read in basis set
   if (read_stand) then
     call read_oldmethod()
   else
     call read_newmethod()
   endif

   call initstuff()
   !> initialize wavefunction
   call initwavefunc()
   !> propagation
   call propagation()
   call deallocate_arrays()
   !> end of mcend run

   call date_and_time(date, time, zone, values)
 !
 !  write(*,'(" "78("*")" ")')
 !  write(*,*) "TRON HALP. TRON PLS HALP! HALP PLS TRON"
 !  write(*,*) ""
 !  write(*,*) ""
   ! exit statement
   write(*,*) ""
   write(*,'(" "78("*")" ")')
   write(*,'(" "78("-")" ")')
   write(*,'("  Normal termination of MCEND on: "20(" "),i4"."i2.2"."i2.2" at "i2.2":"i2.2":"i2.2)') values(1:3), values(5:7)
   wtime = omp_get_wtime() - wtime
   write(*,'("  Wall time (s): "f15.2)') wtime
   write(*,*) ""
   write(*,*) ""
   write(*,*) " End transmission."
   write(*,'(" "78("-")" ")')
   write(*,'(" "78("*")" ")')

 end
 !> @file
 !> @brief contains the main driver
 !! for the input reading, allocation of arrays, and propagation.
