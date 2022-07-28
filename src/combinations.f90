!> Module containing the calculation of possible combinations in determinants for the A-vector.
module combinations
!  use iso_fortran_env
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use params
    use inputvars

    implicit none

    type comb_result
       integer, allocatable :: combs(:)
    end type comb_result

    type(comb_result), pointer :: rcom(:)
    type(comb_result), pointer :: rcom_alpha(:)
    type(comb_result), pointer :: rcom_beta(:)
    type(comb_result), pointer :: rcom_spinorbital(:)
    type(comb_result), pointer :: rcom_so(:)
    type(comb_result), pointer :: rcom2(:)
    type(comb_result), pointer :: rcom2_spinorbital(:)
    contains

!> Pick the other orbitals.
    integer function choose(n, k, err)
      integer, intent(in) :: n, k
      integer, optional, intent(out) :: err
      integer :: imax, i, imin, ie

      ie = 0
      if ((n < 0) .or. (k < 0)) then
        write (*,*) "negative in choose", n,k
        write(ERROR_UNIT,*) "negative in choose"
        choose = 0
        ie = 1
      else
        if (n < k) then
          choose = 0
        else if (n == k) then
          choose = 1
        else
          imax = max(k, n-k)
          imin = min(k, n-k)
          choose = 1
          do i = imax+1, n
            choose = choose*i
          enddo
          do i = 2, imin
            choose = choose/i
          enddo
        endif
      endif

      if (present(err)) err = ie

    end function choose

!> Find possible number of combinations without repetition.
    subroutine comb(n, k, co)
      integer, intent(in) :: n, k
      type(comb_result), pointer, intent(out) :: co(:)

      integer :: i,  s, ix, kx, hm, t !j,
      integer :: err

      hm = choose(n, k, err)

      if (err /= 0) then
        nullify(co)
        return
      end if

      allocate(co(0:hm-1))

      do i = 0, hm - 1
        allocate(co(i)%combs(0:k-1))
      enddo
      do i = 0, hm - 1
        ix = i; kx = k
        do s = 0, n - 1
          if (kx == 0) exit
          t = choose(n - (s+1), kx-1)
          if (ix < t) then
! hash-key mean number, electron
! orbital
            co(i)%combs(kx-1) = s + 1
! electron
            kx = kx - 1
          else
            ix = ix - t
          endif
        enddo
      enddo

    end subroutine comb

  end module

  !> @file
  !> @brief contains the calculation of possible combinations of electrons in
  !! orbitals for number of unique determinants and elements in the MCEND expansion.
