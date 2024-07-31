!> Module containing the analysis routine to write expectation values.
  module analyse

    use params
    use globalvars
    use inputvars
    use omp_lib
    use moreadin
    use utils
!   use hdf5

    contains

!> This routine handles the calculation of excpectation values.
!! @param val_input Contains the computed expectation values.
    subroutine analysis(Ao_input,phio_input, phino, phi_input, val_input, valreal_input, &
    & hel_input, rho_input, A_input, nrindep_input,nrorb_input, nrorb_spinorbital_input, detl_input,hel2_input, &
    & venmatmo_input)

      implicit none
      integer     :: nrindep_input, nrorb_input, nrorb_spinorbital_input
      complex(dp) :: Ao_input(nrindep_input*nrspf)
      complex(dp) :: phio_input(nrprime,nrorb_input)
      complex(dp) :: phino(nrprimn,nrspf)
      complex(dp) :: enmat(nrindep_input,nrindep_input,nrprimn)
      complex(dp) :: csum
      complex(dp) :: scme3(nrprimn)
      complex(dp) :: scme2
      complex(dp) :: val_input(13)
      real(dp)    :: valreal_input(13)
      integer     :: is, i1s, i1r, vorz
      integer     :: nrdiff, ms, ps, qs, ns, alpha, i, j, l, idet, jdet
      integer     :: l1r, l1s,  ix
      integer     :: detl_input(nel*nrindep_input)
      complex(dp) :: hel_input(nrorb_input,nrorb_input,5+nrensp), hel2_input(nrorb_input,nrorb_input,nrorb_input,nrorb_input),&
      & rho_input(nrorb_spinorbital_input,nrorb_spinorbital_input), A_input(nrindep_input*nrspf), venmatmo_input(nrorb_input, &
      & nrorb_input,nrprimn), phi_input(nrprime, nrorb_input)

      val_input(:) = c0
      valreal_input(:) = 0.0_dp

      do is=1, nrspf
        val_input(1) = val_input(1) + rhon(is,is)
      enddo

      csum = c0

      do is=1, nrorb_spinorbital_input
        csum = csum + rho_input(is,is)
      enddo

      valreal_input(1) = dreal(val_input(1))

      do i1s=1, nrorb_spinorbital_input
        if (flag_spinorbital == 0) then
          i1r = (i1s - 1)/2 + 1
        else
          i1r = i1s
        end if
        do l1s=1, nrorb_spinorbital_input
          if (flag_spinorbital == 0) then
            l1r = (l1s - 1)/2 + 1
          else
            l1r =  l1s
          end if
          val_input(2:4) = val_input(2:4) - rho_input(i1s,l1s)*hel_input(i1r,l1r,2:4)*deltaspin(i1s,l1s) 
        enddo
      enddo

      valreal_input(2:4) = dreal(val_input(2:4))/dreal(val_input(1))

      do ix=1,nrprimn
        do i=1,nrspf
          do j=1,nrspf
            val_input(5) = val_input(5) + r(ix)*dconjg(phin(ix,i))*phin(ix,j)*rhon(i,j)
          enddo
        enddo
      enddo

      valreal_input(5) = dreal(val_input(5))*dr/dreal(val_input(1))

      ! energies
      ! kinetic
      do i1s=1, nrorb_spinorbital_input
        if (flag_spinorbital == 0) then
          i1r = (i1s - 1)/2 + 1
        else
          i1r = i1s
        end if
        do l1s=1, nrorb_spinorbital_input
          if (flag_spinorbital == 0) then
            l1r = (l1s - 1)/2 + 1
          else
            l1r = l1s
          end if
          val_input(6) = val_input(6) + rho_input(i1s,l1s)*hel_input(i1r,l1r,1)*deltaspin(i1s,l1s)

        enddo
      enddo

      valreal_input(6) = nel*dreal(val_input(6))/dreal(val_input(1))

      do i=1, nrspf
        do j=1, nrspf
          val_input(7) = val_input(7) + heln(i,j,1)*rhon(j,i)
        enddo
      enddo
      valreal_input(7) = dreal(val_input(7))/dreal(val_input(1))

      ! repulsion terms
      do idet=1, nrindep_input
        do jdet=1, nrindep_input
          call maxcoinc(idet,jdet,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)
          call sc2(nrdiff,idet,jdet,ms,ns,ps,qs,vorz,scme2, detl_input, nrindep_input,&
  & hel2_input, nrorb_input)

          do j=1, nrspf

            val_input(8) = val_input(8) + dconjg(A_input(idet+(j-1)*nrindep_input)) &
  & *A_input(jdet+(j-1)*nrindep_input)*scme2

          enddo

        enddo
      enddo

      valreal_input(8) = dreal(val_input(8))/dreal(val_input(1))

      do i=1,nrspf
        do j=1,nrspf
          val_input(9) = val_input(9) + heln(i,j,2)*rhon(j,i)
        enddo
      enddo

      valreal_input(9) = dreal(val_input(9))/dreal(val_input(1))

      ! electron-nuclear attraction
      do idet=1, nrindep_input
        do jdet=1, nrindep_input
          call maxcoinc(idet,jdet,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)
          call sc1en(nrdiff,idet,ms,ps,vorz,scme3, detl_input, nrindep_input, venmatmo_input, &
          &nrorb_input)
          do alpha=1, nrprimn
            enmat(idet,jdet,alpha) = scme3(alpha)
          enddo
        enddo
      enddo

      do idet=1, nrindep_input
        do jdet=1, nrindep_input
          do j=1, nrspf
            do l=1, nrspf
              do alpha=1, nrprimn
                val_input(10) = val_input(10) + dconjg(A_input(idet+(j-1)*nrindep_input)) &
                          *A_input(jdet+(l-1)*nrindep_input)*enmat(idet,jdet,alpha) &
                          *dconjg(phin(alpha,j))*phin(alpha,l)*(r(2) - r(1))
              enddo
            enddo
          enddo
        enddo
      enddo

      valreal_input(10) = dreal(val_input(10))/dreal(val_input(1))
      valreal_input(11) = sum(valreal_input(6:10))
      val_input(12) = c0
      if (aucofu /= 0) then
        if (flag_spinorbital==0) then
          val_input(13) = overlap(Ao_input,phio_input, phino, nrorb_input, nrindep_input, detl_input, phi_input,&
          & A_input)
        end if
        if (flag_spinorbital==1) then
          val_input(13) = overlap(Ao_input,phio_input, phino, nrorb_input, nrindep_input, detl_input, phi_input,&
          & A_input)
        end if
      endif
      valreal_input(12) = dreal(val_input(13))
      valreal_input(13) = aimag(val_input(13))

      return

    end subroutine

!> Calculates the overlap of two wave functions.
    complex(dp) function overlap(Ao_input,phio_input, phino, nrorb_input, nrindep_input, detl_input, phi_input, &
    & A_input)
      implicit none
      complex(dp) :: bigovl_input(nrorb_input,nrorb_input), ovl(nel,nel)
      complex(dp) :: Ao_input(nrindep_input*nrspf)
      complex(dp) :: phio_input(nrprime,nrorb_input), deter, work(2*nel), overl
      complex(dp) :: w_input(nel), vl(nel,nel), vr(nel,nel)
      complex(dp) :: bigovln(nrspf,nrspf)
      complex(dp) :: phino(nrprimn,nrspf)
      real(dp)    :: rwork(2*nel)
      integer     :: lwork, info, ir, jr, ix, cis, cjs, ie, ie2, jn, ln
      integer     :: is, js
      character   :: jobvr, jobvl
      integer     :: nrorb_input, nrindep_input, detl_input(nel*nrindep_input)
      complex(dp) :: phi_input(nrprime,nrorb_input),  A_input(nrindep_input*nrspf)

      jobvr = "N"
      jobvl = "N"
      lwork = 2*nel
      info = 0
      do ir=1, nrorb_input
        do jr=1, nrorb_input
          bigovl_input(ir,jr) = c0
          do ix=1, nrprime
            bigovl_input(ir,jr) = bigovl_input(ir,jr) + dconjg(phio_input(ix,ir))*phi_input(ix,jr)
          enddo
        enddo
      enddo

      do jn=1, nrspf
        do ln=1, nrspf
          bigovln(jn,ln) = c0
          do ix=1, nrprimn
            bigovln(jn,ln) = bigovln(jn,ln) + dconjg(phino(ix,jn))*phin(ix,ln)
          enddo
          bigovln(jn,ln) = bigovln(jn,ln)*dr
        enddo
      enddo

      overl = c0

      do cis=1, nrindep_input
        do cjs=1, nrindep_input
          do ie=1, nel
            is = detl_input((cis-1)*nel+ie)
            if (flag_spinorbital == 0) then
              ir = (is - 1)/2 + 1
            else
              ir = is
            end if
            do ie2=1,nel
              js = detl_input((cjs-1)*nel+ie2)
              if (flag_spinorbital == 0) then
                jr = (js-1)/2 + 1
              else
                jr = js
              end if

              ovl(ie,ie2) = bigovl_input(ir,jr)*deltaspin(is,js)
            enddo
          enddo
          call zgeev(jobvl,jobvr,nel,ovl,nel,w_input,vl,nel,vr,nel,work,lwork,rwork,info)
          if (info /= 0) then
            write(*,*) "ERROR:"
            write(*,*) "info zgeev(overlap)= ", info
            stop
          endif
          deter = cr
          do ie=1, nel
            deter = deter*w_input(ie)
          enddo
          do jn=1,nrspf
            do ln=1,nrspf
              overl = overl + dconjg(Ao_input(cis+(jn-1)*nrindep_input)) &
                             *A_input(cjs+(ln-1)*nrindep_input)*deter*bigovln(jn,ln)
            enddo
          enddo
        enddo
      enddo

      overlap = overl

      return
    end function

  end module

!> @file
!> @brief contains the analysis routine
!! for the calculation of expectation values.
