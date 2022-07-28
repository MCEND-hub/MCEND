!> Module containing the Akima interpolation for the electron-nuclear attraction integrals.
module akima_interpol
  use params

  implicit none

  contains

!> Smooth the curve for large distances between grid points.
    subroutine abs_smooth(x, delta_x, y)

      real(dp), intent(in) :: x, delta_x
      real(dp), intent(out) :: y

      if (x >= delta_x) then
        y = x
      elseif (x <= -delta_x) then
        y = -x
      else
        y = x**2/(2.0_dp*delta_x) + delta_x/2.0_dp
      end if

    end subroutine abs_smooth

!> Main Akima interpolation routine. This interpolates the electron-nuclear one-electron integrals
!! onto the finer nuclear grid that is used in the MCEND propagation.
!!@param npt number of points
    subroutine akima_interp(npt, n, x, xpt, ypt, y)

      real(dp), parameter   :: eps = 1.0d-30
      integer,  intent(in)  :: n
      integer,  intent(in)  :: npt
      real(dp), intent(in)  :: xpt(npt), ypt(npt)
      real(dp), intent(in)  :: x(n)
      real(dp), intent(out) :: y(n)
      integer  :: i, j
      real(dp) :: delta_x
      real(dp) :: p0(npt-1), p1(npt-1), p2(npt-1), p3(npt-1)
      real(dp) :: m(-1:npt+1)
      real(dp) :: t(npt)
      real(dp) :: w1, w2
      real(dp) :: t1, t2, dx

      delta_x = 0.1_dp

      do i = 1, npt-1
        m(i) = (ypt(i+1) - ypt(i))/(xpt(i+1) - xpt(i))
      enddo

      m(0)  = 2.0_dp*m(1) - m(2)
      m(-1) = 2.0_dp*m(0) - m(1)
      m(npt) = 2.0_dp*m(npt-1) - m(npt-2)
      m(npt+1) = 2.0_dp*m(npt) - m(npt-1)

      ! the slope at the points
      do i = 1, npt
        call abs_smooth(m(i+1) - m(i), delta_x, w1)
        call abs_smooth(m(i-1) - m(i-2), delta_x, w2)
        if (w1 < eps .and. w2 < eps) then
          ! special case so we don't divide by zero
          t(i) = 0.5_dp*(m(i-1) + m(i))
        else
          t(i) = (w1*m(i-1) + w2*m(i))/(w1 + w2)
        end if
      enddo

      ! polynomial coefficients
      do i = 1, npt-1
        dx = xpt(i+1) - xpt(i)
        t1 = t(i)
        t2 = t(i+1)
        p0(i) = ypt(i)
        p1(i) = t1
        p2(i) = (3.0_dp*m(i) - 2.0_dp*t1 - t2)/dx
        p3(i) = (t1 + t2 - 2.0_dp*m(i))/dx**2
      end do

      ! interpolate at each point
      do i = 1, n
        ! find location in array (use end segments if out of bounds)
        if (x(i) < xpt(1)) then
          j = 1
        else
          do j = npt-1, 1, -1
            if (x(i) >= xpt(j)) exit
          enddo
        endif
        dx = x(i) - xpt(j)
        y(i) = p0(j) + p1(j)*dx + p2(j)*dx**2 + p3(j)*dx**3

      enddo

    end subroutine akima_interp

end module akima_interpol

!> @file
!> @brief contains the Akima interpolation
!! for the one-electron integrals in AO basis.
