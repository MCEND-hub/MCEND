!> This module contains utility functions such as the electric field generation.
module utils
  use params
  use globalvars
  use inputvars
  use omp_lib

  contains

    !> max coincidence for single hole functions (SHF)
    !! Phi_1 = |abcd> and Phi_2 = |crds>
    !! |crds> = -|crsd> = |srcd>
  subroutine maxcoincshf(cisshf,cjsshf,vorz,nrdiff,ms,ns,ps,qs, shdl_input, max_nrindep_input)

    implicit none
    integer :: cisshf, cjsshf, vorz, nrdiff, ms, ns, ps, qs
    integer :: ie, je, hh, iarr(nel-1), iarr2(nel-1), max_nrindep_input
    integer :: shdl_input( max_nrindep_input*(nel-1) )

    vorz = 1
    nrdiff = 0
    ms = 0
    ns = 0
    ps = 0
    qs = 0

    do ie=1, nel-1
      iarr(ie)  = shdl_input((cisshf-1)*(nel-1)+ie)
      iarr2(ie) = shdl_input((cjsshf-1)*(nel-1)+ie)
    enddo

    ! if the same number appears in different positions, then swap
    do ie=1, nel-1
      do je=1, nel-1
        if ((iarr(ie) == iarr2(je)) .and. (ie /= je)) then
          hh = iarr2(ie)
          iarr2(ie) = iarr2(je)
          iarr2(je) = hh
          vorz = -vorz
        endif
      enddo
    enddo

    do ie=1, nel-1
      if (iarr(ie) /= iarr2(ie)) then
        nrdiff = nrdiff + 1
        if (nrdiff == 1) then
          ms = iarr(ie)
          ps = iarr2(ie)
        else
          ns = iarr(ie)
          qs = iarr2(ie)
        endif
      endif
    enddo

    return
  end subroutine

!> Maximum number of coincidence
 subroutine maxcoinc(cis,cjs,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)

    implicit none
    integer :: cis, cjs, vorz, nrdiff, ms, ns, ps, qs
    integer :: ie, je, hh, iarr(nel), iarr2(nel)
    integer :: nrindep_input, detl_input(nel*nrindep_input)

    vorz = 1
    nrdiff = 0
    ms = 0
    ns = 0
    ps = 0
    qs = 0

    do ie=1,nel
      iarr(ie) = detl_input((cis-1)*nel+ie)
      iarr2(ie) = detl_input((cjs-1)*nel+ie)
    enddo

    ! if the same number appears in different positions, then swap
    do ie=1,nel
      do je=1,nel
        if ((iarr(ie) == iarr2(je)) .and. (ie /= je)) then
          hh = iarr2(ie)
          iarr2(ie) = iarr2(je)
          iarr2(je) = hh
          vorz = -vorz
        endif
      enddo
    enddo

    do ie=1,nel
      if (iarr(ie) /= iarr2(ie)) then
        nrdiff = nrdiff + 1
        if (nrdiff == 1) then
          ms = iarr(ie)
          ps = iarr2(ie)
        else
          ns = iarr(ie)
          qs = iarr2(ie)
        endif
      endif
    enddo

    return
  end subroutine

  !> calculates matrix elements of two particle operator (sum over
  !! e-e-repulsion), for nel electrons
  !!cw <shf_1 | Vee| shf_2>
  subroutine sc2(nrdiff,cis,cjs,ms,ns,ps,qs,vorz,scme2, detl_input, nrindep_input,&
  & hel2_input, nrorb_input)

  implicit none

  complex(dp) :: scme2
  integer :: cis, cjs
  integer :: is, ir, ie, je, ds, js, jr, ms, mr, ps, pr, qs, qr
  integer :: ns, nr, ds1, ds2, vorz, nrdiff
  integer :: nrindep_input, detl_input(nel*nrindep_input), nrorb_input, sum_sz_tmp, orb_temp, i
  complex(dp) :: hel2_input(nrorb_input,nrorb_input,nrorb_input,nrorb_input)

  scme2 = c0
  if (flag_spinorbital == 0) then
    sum_sz_tmp = 0
    do i = 1, nel
      orb_temp = detl( (cis - 1) *nel + i )
      if (mod(orb_temp,2) == 0) then
        sum_sz_tmp = sum_sz_tmp - 1
      else if (mod(orb_temp,2) == 1) then
        sum_sz_tmp = sum_sz_tmp + 1
      else
        write (*,*) 'error in mod 2'
        stop
      end if
    end do

    if (sum_sz_tmp .ne. nsz) then
      return
    end if

    sum_sz_tmp = 0
    do i = 1, nel
      orb_temp = detl( (cjs - 1) *nel + i )
      if (mod(orb_temp,2) == 0) then
        sum_sz_tmp = sum_sz_tmp - 1
      else if (mod(orb_temp,2) == 1) then
        sum_sz_tmp = sum_sz_tmp + 1
      else
        write (*,*) 'error in mod 2'
        stop
      end if
    end do
    if (sum_sz_tmp .ne. nsz) then
      return
    end if
  end if

  if (nrdiff == 0) then
    do ie=1, nel
      is = detl_input((cis-1)*nel+ie)
      if (flag_spinorbital == 0) then
        ir = (is - 1)/2 + 1
      else if (flag_spinorbital == 1) then
        ir = is
      else
        write (*,*) 'flag_spinorbital nor 0/1 not supported'
        stop
      end if
      do je=1, nel
        js = detl_input((cjs-1)*nel+je)
        ds = deltaspin(is,js)
        if (flag_spinorbital == 0) then
          jr = (js - 1)/2 + 1
        else if (flag_spinorbital == 1) then
          jr = js
        else
          write (*,*) 'flag_spinorbital nor 0/1 not supported'
          stop
        end if
        scme2 = scme2 + hel2_input(ir,jr,ir,jr) - hel2_input(ir,jr,jr,ir)*ds
      enddo
    enddo
    scme2 = 0.50_dp*scme2

  elseif (nrdiff == 1) then
    do ie=1, nel
      is = detl_input((cis-1)*nel+ie)
      if (flag_spinorbital == 0) then
        ir = (is - 1)/2 + 1
        mr = (ms - 1)/2 + 1
        pr = (ps - 1)/2 + 1
      else if (flag_spinorbital == 1) then
        ir = is
        mr = ms
        pr = ps
      else
        write (*,*) 'flag_spinorbital nor 0/1 not supported'
        stop
      end if
      ds1 = deltaspin(ms,ps)
      ds2 = deltaspin(ms,is)*deltaspin(ps,is)
      scme2 = scme2 + hel2_input(mr,ir,pr,ir)*ds1 - hel2_input(mr,ir,ir,pr)*ds2
    enddo
    scme2 = vorz*scme2
  elseif (nrdiff == 2) then
    if (flag_spinorbital == 0) then
      mr = (ms - 1)/2 + 1
      nr = (ns - 1)/2 + 1
      pr = (ps - 1)/2 + 1
      qr = (qs - 1)/2 + 1
    else if (flag_spinorbital == 1) then
      mr = ms
      nr = ns
      pr = ps
      qr = qs
    else
      write (*,*) 'flag_spinorbital nor 0/1 not supported'
      stop
    end if
    ds1 = deltaspin(ms,ps)*deltaspin(ns,qs)
    ds2 = deltaspin(ms,qs)*deltaspin(ns,ps)
    !> antisymmetrize
    scme2 = hel2_input(mr,nr,pr,qr)*ds1 - hel2_input(mr,nr,qr,pr)*ds2
    scme2 = vorz*scme2
  endif

  return

  end subroutine

  !> calculates matrix elements of single particle operator,
  !! for nel electrons, scme1: Slater-Condon matrix element
 subroutine sc1en(nrdiff,cis,ms,ps,vorz,scme3, detl_input, nrindep_input, venmatmo_input, nrorb_input)

    implicit none

    complex(dp) :: scme3(nrprimn)
    integer :: cis, ms, ps, nrdiff, vorz, is, ir
    integer :: ie, mr, pr, aa
    integer :: nrindep_input, detl_input(nel*nrindep_input), nrorb_input
    complex(dp) :: venmatmo_input(nrorb_input,nrorb_input,nrprimn)

    do aa=1,nrprimn
      scme3(aa) = c0
    enddo
    if (nrdiff == 0) then
      do ie=1, nel
        is = detl_input((cis-1)*nel+ie)
        if (flag_spinorbital==0) then
          ir = (is - 1)/2 + 1
        else if (flag_spinorbital==1) then
          ir = is
        else
          write (*,*) "flag_spinorbital nor 0/1 not suported"
          stop
        end if
        do aa=1,nrprimn
          scme3(aa) = scme3(aa) + venmatmo_input(ir,ir,aa)
        enddo
      enddo
    elseif (nrdiff == 1) then
      if (flag_spinorbital == 0) then
        mr = (ms - 1)/2 + 1
        pr = (ps - 1)/2 + 1
      else if (flag_spinorbital == 1) then
        mr = ms
        pr = ps
      else
        write (*,*) "flag_spinorbital nor 0/1 not suported"
        stop
      end if
      do aa=1, nrprimn
        scme3(aa) = vorz*venmatmo_input(mr,pr,aa)*deltaspin(ms,ps)
      enddo
    endif
    return
  end subroutine

  !> calculates matrix elements of single particle operator,
  !! for nel electrons, scme1 : Slater-Condon matrix element
  subroutine sc1t(nrdiff,idet,ms,ps,vorz,scme1, detl_input, nrindep_input, hel_input, nrorb_input)

    implicit none

    complex(dp) :: scme1
    integer :: idet, ms, ps, nrdiff, vorz, is, ir
    integer :: ie
    integer :: mr, pr
    integer :: nrindep_input, nrorb_input
    integer :: detl_input(nel*nrindep_input)
    complex(dp) :: hel_input(nrorb_input,nrorb_input,5+nrensp)

    scme1 = c0
    if (nrdiff == 0) then
      do ie=1, nel
        is = detl_input((idet-1)*nel+ie)
        if (flag_spinorbital == 0) then
          ir = (is - 1)/2 + 1
        else if (flag_spinorbital == 1) then
          ir = is
        else
          write (*,*) 'flag_spinorbital not 0 or 1 is unsupported'
          stop
        end if
        scme1 = scme1 + hel_input(ir,ir,1)
      enddo
    elseif (nrdiff == 1) then
      if (flag_spinorbital == 0) then
        mr = (ms - 1)/2 + 1
        pr = (ps - 1)/2 + 1
      else if (flag_spinorbital == 1) then
        mr = ms
        pr = ps
      else
        write (*,*) 'flag_spinorbital not 0 or 1 is unsupported'
        stop
      end if
      scme1 = vorz*hel_input(mr,pr,1)*deltaspin(ms,ps)
    endif
    return
  end subroutine

  !> calculates matrix elements of single particle operator,
  !! for nel electrons in the dipole moment part of the electric field,
  !! scme1 : Slater-Condon matrix element
  subroutine sc1mu(nrdiff,cis,ms,ps,vorz,scme1mu, detl_input, nrindep_input, hel_input, nrorb_input)

    implicit none

    complex(dp) :: scme1mu(3)
    integer :: cis, ms, ps, nrdiff, vorz, is, ir
    integer :: ie
    integer :: mr, pr
    integer :: nrindep_input, nrorb_input
    integer :: detl_input(nel*nrindep_input)
    complex(dp) :: hel_input(nrorb_input,nrorb_input,5+nrensp)

    scme1mu(1) = c0
    scme1mu(2) = c0
    scme1mu(3) = c0
    if (nrdiff == 0) then
      do ie=1,nel
        is = detl_input((cis - 1)*nel + ie)
        if (flag_spinorbital == 0) then
          ir = (is - 1)/2 + 1
        else if (flag_spinorbital == 1) then
          ir = is
        else
          write (*,*) 'flag_spinorbital not 0/1 is unsupported'
          stop
        end if
        scme1mu(1) = scme1mu(1) + hel_input(ir,ir,2)
        scme1mu(2) = scme1mu(2) + hel_input(ir,ir,3)
        scme1mu(3) = scme1mu(3) + hel_input(ir,ir,4)
      enddo
    elseif (nrdiff == 1) then
      if (flag_spinorbital == 0) then
        mr = (ms - 1)/2 + 1
        pr = (ps - 1)/2 + 1
      else if (flag_spinorbital == 1) then
        mr = ms
        pr = ps
      else
        write (*,*) 'flag_spinorbital not 0/1 is unsupported'
        stop
      end if
      scme1mu(1) = vorz*hel_input(mr,pr,2)*deltaspin(ms,ps)
      scme1mu(2) = vorz*hel_input(mr,pr,3)*deltaspin(ms,ps)
      scme1mu(3) = vorz*hel_input(mr,pr,4)*deltaspin(ms,ps)
    endif

    return
  end subroutine

!> Generate the electric field (laser pulse)
  subroutine getefield(t, ef)

    implicit none

    real(dp) ::  t, ef(3)
    real(dp) :: lpw(nrpulses)
    real(dp) :: t_now
    integer :: i

    lpw(:) = lpw0(:)/au2fs
    t_now = 0.0_dp
    ef(:) = 0.0_dp

    if (t < lpw(1)) then
      ef(:) = lpolar(:,1)*lph(1)*(dsin(pi*t/lpw(1))**2)*dcos(lfreq(1)*t)
      return
    elseif (t >= lpw(1) .and. t < sum(lpw(:))) then
      do i = 2, nrpulses
        t_now = sum(lpw(:i-1))
        if (t >= t_now .and. t < t_now + lpw(i)) then
          ef(:) = lpolar(:,i)*lph(i)*(dsin(pi*(t - t_now)/(lpw(i)))**2)*dcos(lfreq(i)*(t - t_now))
          return
        endif
      enddo
    elseif (t >= sum(lpw(:))) then
      ef(:) = 0.0_dp
      return
    endif

  end subroutine

!> Get the nuclear momentum for electronic translation factors
  subroutine get_nuc_momentum()

    implicit none

    complex(dp) :: phin3(nrprimn,nrspf)
    complex(dp) :: mymoment(nrspf)
    complex(dp) :: nuc_moment
    complex(dp) :: phibuff(nrprimn), phibuff2(nrprimn)
    integer :: i, j, ix
    integer :: fftw_forward=-1,fftw_backward=1
    integer(i64) :: planf, planb

    call dfftw_plan_dft_1d(planf,nrprimn,phibuff,phibuff2,fftw_forward,0)
    call dfftw_plan_dft_1d(planb,nrprimn,phibuff2,phibuff,fftw_backward,0)

    do i=1,nrspf
      do ix=1,nrprimn
        phibuff(ix) = phin(ix,i)
      enddo
      call dfftw_execute_dft(planf,phibuff,phibuff2)
      do ix=1,nrprimn
        phibuff2(ix) = kx(ix)*phibuff2(ix)
      enddo

      call dfftw_execute_dft(planb,phibuff2,phibuff)
      do ix=1,nrprimn
        phin3(ix,i) = phibuff(ix)/nrprimn
      enddo
    enddo

    mymoment(:) = c0
    nuc_moment = c0
    do i=1, nrspf
      do ix=1, nrprimn
        mymoment(i) = mymoment(i) + cdabs(phin3(ix,i))
      enddo
    enddo

    mymoment(:) = mymoment(:)*dr

    do i=1, nrspf
      do j=1, nrindep
        nuc_moment = nuc_moment + mymoment(i)*A(j+(i-1)*nrindep)
      enddo
    enddo

  end subroutine

    !> find out if the two SOs i and j have the same sz
    !! altering names of m1 and m2 as a precaution, since they are already assigned
    !! to the masses
    !! cw: assume +-+-...  extended for open-shell systems
    !! cw: e.g., +- ++
    integer function deltaspin(i,j)
      integer :: i, j, ijmax, ijmin

      ijmax = max(i,j)
      ijmin = min(i,j)
      if (mod(ijmax - ijmin, 2) == 0) then
        deltaspin = 1
      else
        deltaspin = 0
      endif
      return
    end function

    integer function detl_nsz(idet, nrindep_input, detl_input)
      integer :: idet, nrindep_input, detl_input(nel*nrindep_input)
      integer :: i, orb_temp, iarr(nel)

        detl_nsz = 0
        do i = 1, nel
          orb_temp = detl_input( (idet - 1) *nel + i )
          iarr(i) = orb_temp
        end do

        detl_nsz = checksz(iarr)

    end function

    !> maps A, phi -> psi for RK8 algorithm if flag=0
    !! maps psi -> A, phi otherwise
  subroutine mapformat(psi_input, flag, dgldim_input, nrindep_input, nrorb_input, A_input, phi_input)

  implicit none

  integer ::  dgldim_input, nrindep_input, nrorb_input
  complex(dp) :: A_input(nrindep_input*nrspf)
  complex(dp) :: phi_input(nrprime,nrorb_input)
  integer(8) :: flag
  complex(dp) ::  psi_input(dgldim_input)
  integer :: k, l, m

  k = nrindep_input*nrspf
  l = nrorb_input*nrprime
  m = nrspf*nrprimn
  if (flag == 0) then
    psi_input(:) = c0
    psi_input(1:k) = A_input(:)
    psi_input(k+1:k+l) = reshape(phi_input,[l])
    psi_input(k+l+1:k+l+m) = reshape(phin,[m])
  else
    A_input(:) = psi_input(1:k)
    phi_input(:,:) = reshape(psi_input(k+1:k+l), [nrprime,nrorb_input])
    phin(:,:) = reshape(psi_input(k+l+1:k+l+m), [nrprimn,nrspf])
  endif
  return

  end subroutine

  !> maps A, phi -> psi_spinorbital for RK8 algorithm if flag=0
  !! maps psi_spinorbital -> A, phi otherwise
  subroutine mapformat2(lA_input, psi_input, flag, phi2_input, nrindep_input, dgldim_input, nrorb_input )

    implicit none

    integer :: dgldim_input
    integer(8) :: flag
    complex(dp) ::  psi_input(dgldim_input)
    complex(dp) :: lA_input(nrindep_input*nrspf), phi2_input(nrprime,nrorb_input)
    integer :: i, j, k, l, m, nrindep_input, nrorb_input

    k = nrindep_input*nrspf
    l = nrorb_input*nrprime
    m = nrspf*nrprimn

    if (flag == 0) then
      psi_input(:) = c0
      psi_input(1:k) = lA_input(:)
      psi_input(k+1:k+l) = reshape(phi2_input,[l])
      do i=1, nrprime
        do j=1, nrorb_input
              if ( mod(j, 2) == 1 .and. j < nrorb_input ) then
                 if ( dabs( dreal( phi2_input(i,j) ) - dreal( phi2_input(i,j+1) ) ) > 1.e-9) then
                 end if
               end if
        end do
      end do
      psi_input(k+l+1:k+l+m) = reshape(phin2,[m])
    else
      lA_input(:) = psi_input(1:k)
      phi2_input(:,:) = reshape(psi_input(k+1:k+l), [nrprime,nrorb_input])
      phin2(:,:) = reshape(psi_input(k+l+1:k+l+m), [nrprimn,nrspf])
    endif

    return

  end subroutine

  

!> to verify matrix inverse of rho and overlap
  subroutine matrix_multiplication(a,b,c,dim_m)
    integer :: i,j,k,dim_m
    complex(dp)    :: a(dim_m,dim_m), b(dim_m,dim_m), c(dim_m,dim_m)

    do i=1,dim_m
      do j=1,dim_m
        c(i,j) = 0.0
      enddo
    enddo

    do i=1,dim_m
      do j=1,dim_m
        do k=1,dim_m
          c(i,j) = c(i,j) + a(i,k)*b(k,j)
        end do
      enddo
    enddo

    write (*,*) 'result of matrix_multiplication'

    do i=1,dim_m
      do j=1,dim_m
        write (*,*) i,j, c(i,j)
      enddo
    enddo

    return
  end subroutine

  subroutine print_phi(flag_print)
    integer        :: i,mu, flag_print
    write (*,*) 'print_phi'

    if (flag_print == 0) then
      do mu=1, nrprime
        do i=1, nrorb
          write (*,*) mu,i,dreal(phi(mu,i))
        enddo
      enddo

    else if (flag_print == 1) then
      do mu=1, nrprime
        do i=1, nrorb_spinorbital
          write (*,*) mu,i,dreal(phi_spinorbital(mu,i))
        enddo
      enddo

    else if (flag_print == 2) then
      do mu=1, nrprime
        do i=1, nrorb
          write (*,*) 'phi',mu,i,dreal(phi(mu,i))
        enddo
      enddo

      do mu=1, nrprime
        do i=1, nrorb_spinorbital
          write (*,*) 'phi_spinorbital', mu,i,dreal(phi_spinorbital(mu,i))
        enddo
      enddo

    !compare
    else if (flag_print == 3) then
      do mu=1, nrprime
        do i=1, nrorb
            if ( dabs( dreal(phi(mu,i)) -  dreal(phi_spinorbital(mu,2*i-1)) )>1.e-9 )  then
              write (*,*) 'warning in phi and phi_spinorbital', mu,i, dreal(phi(mu,i)), dreal(phi_spinorbital(mu,2*i-1))
            end if
            if ( dabs( dreal(phi(mu,i)) -  dreal(phi_spinorbital(mu,2*i)) )>1.e-9 )  then
              write (*,*) 'warning phi and phi_spinorbital', mu,i, dreal(phi(mu,i)), dreal(phi_spinorbital(mu,2*i))
            end if
        enddo
      enddo
    end if

  end subroutine



  subroutine print_phin()
    integer        :: i,ix
    write (*,*) 'print_phin'
    do i=1, nrspf
      do ix=1, nrprimn
        write (*,*) i, ix, dreal(phin(ix,i))
      enddo
    enddo
  end subroutine

  subroutine print_complex_array_2d(array, dim_1, dim_2)
    integer        :: i,j, dim_1, dim_2
    complex(dp)    :: array(dim_1,dim_2)
  ! 3 and 6 are used to reduce output
    write (*,*) 'print_array_2d'
    do i=1, dim_1
      do j=1, dim_2
        write (*,*)  i,j
        write (*,*) dreal(array(i,j))
      enddo
    enddo
  end subroutine

  subroutine print_mf1(flag_print, mf1_input, nrorb_input)
    integer :: i, j, mu, nu
    integer :: flag_print, unit_orbital, unit_spinorbital
    integer :: nrorb_input
    complex(dp) ::  mf1_input(nrorb_spinorbital,nrorb_spinorbital,nrprime,nrprime)


    if (flag_print == 0) then

      open(newunit=unit_orbital,file="mf1_output.txt")

      do i=1, nrorb_input
        do j=1, nrorb_input
          do mu=1, nrprime
            do nu=1, nrprime
              write (*,*) i, j, mu, nu, dreal(mf1_input(i,j,mu,nu))
              write (unit_orbital,*)  i, j, mu, nu, dreal(mf1_input(i,j,mu,nu))
            enddo
          enddo
        enddo
      enddo
      close(unit_orbital)

    end if

    if (flag_print == 1) then

      open(newunit=unit_spinorbital,file="mf1_spinorbital_output.txt")
      do i=1, nrorb_input
        do j=1, nrorb_input
          do mu=1, nrprime
            do nu=1, nrprime
              write (*,*) i, j, mu, nu, dreal(mf1_input(i,j,mu,nu))
              write (unit_spinorbital,*)  i, j, mu, nu, dreal(mf1_input(i,j,mu,nu))
            enddo
          enddo
        enddo
      enddo
      close(unit_spinorbital)

    end if

  end subroutine

  subroutine print_hel(flag_hel, flag_print)
    integer :: i, j, k, l, mu, nu
    integer :: flag_hel, flag_print, unit_orbital, unit_spinorbital

    if (flag_hel == 1) then
      if (flag_print == 0) then
        open(newunit=unit_orbital,file="hel_output.txt")
        do i=1, nrorb
          do j=1, nrorb
            write (*,*) i,j, hel(i,j,1)
            write (unit_orbital,*) i,j, hel(i,j,1)
          enddo
        enddo
        close(unit_orbital)
      end if

      if (flag_print == 1) then
        open(newunit=unit_spinorbital,file="hel_spinorbital_output.txt")
        do i=1, nrorb_spinorbital
          do j=1, nrorb_spinorbital
            write (unit_spinorbital,*) i,j, hel_spinorbital(i,j,1)
          enddo
        enddo
        close(unit_spinorbital)
      end if
    end if

    if (flag_hel == 2) then
      if (flag_print == 0) then
        open(newunit=unit_orbital,file="hel2_output.txt")
        do i=1, nrorb
          do j=1, nrorb
            do k=1, nrorb
              do l=1, nrorb
                write(*,*)  i, j, k, l, dreal(hel2(i,j,k,l))
                write(unit_orbital,*)  i, j, k, l, dreal(hel2(i,j,k,l))
              enddo
            enddo
          enddo
        enddo
        close(unit_orbital)
      end if

      if (flag_print == 1) then

        open(newunit=unit_spinorbital,file="hel2_spinorbital_output.txt")
        do i=1, nrorb_spinorbital
          do j=1, nrorb_spinorbital
            do k=1, nrorb_spinorbital
              do l=1, nrorb_spinorbital
                write(*,'(I3,I3,I3,I3,F15.9)')  i, j, k, l, dreal(hel2_spinorbital(i,j,k,l))
                write(unit_spinorbital,'(I3,I3,I3,I3,F15.9)')  i, j, k, l, dreal(hel2_spinorbital(i,j,k,l))
              enddo
            enddo
          enddo
        enddo
        close(unit_spinorbital)
      end if
    end if

    if (flag_hel == 3) then

      if (flag_print == 0) then
        open(newunit=unit_orbital,file="hel3_output.txt")
        do mu=1, nrprime
          do i=1, nrorb
            do nu=1, nrprime
              do j=1, nrorb
                write (*,*) mu,i,nu,j, hel3(mu,i,nu,j)
                write (unit_orbital,*) mu,i,nu,j, hel3(mu,i,nu,j)
              enddo
            enddo
          enddo
        enddo
        close(unit_orbital)
      end if

      if (flag_print == 1) then
        open(newunit=unit_spinorbital,file="hel3_spinorbital_output.txt")
        do mu=1, nrprime
          do i=1, nrorb_spinorbital
            do nu=1, nrprime
              do j=1, nrorb_spinorbital
                write (*,*) mu,i,nu,j, hel3_spinorbital(mu,i,nu,j)
                write (unit_spinorbital,*) mu,i,nu,j, hel3_spinorbital(mu,i,nu,j)
              enddo
            enddo
          enddo
        enddo
        close(unit_spinorbital)
      end if
    end if

  end subroutine

  subroutine compare_scalar_real(A,B)
    real(dp) ::  A,B
    real(dp) ::  threshold_compare

    threshold_compare  = 1.0d-9
    if ( dabs(A-B) > threshold_compare ) then
       write (*,*) 'warning,  A != B ', A,B
    end if

  end subroutine

  subroutine compare_vector_integer(A,B,dim_A,dim_B)
    integer  ::  dim_A, dim_B, i
    integer  ::  A(dim_A),B(dim_B)
    integer  ::  flag_consistent
    real(dp) ::  threshold_compare
    threshold_compare  = 1.0d-9
    flag_consistent = 0

    dim_A = size(A)
    dim_B = size(B)

    if ( dim_A /= dim_B ) then
      flag_consistent = 1
      write (*,*) 'inconsistent shape'
    else
      do i = 1, dim_A
        if  ( A(i) /= B(i)  ) then
            flag_consistent = 1
            write (*,*) 'warning, A(i) != B(i)', i, A(i), B(i)
        end if
       end do
    end if

    if (flag_consistent == 0) then
      write (*,*) 'consistent in compare_vector_integer'
    end if

  end subroutine

  subroutine compare_matrix_real(A,B,dim_A,dim_B)
    integer  ::  dim_A(2), dim_B(2), i,j
    real(dp) ::  A(dim_A(1),dim_A(2)),B(dim_B(1),dim_B(2))
    real(dp) ::  threshold_compare

    threshold_compare  = 1.0d-9

    dim_A = shape(A)
    dim_B = shape(B)

    if ( dim_A(1) /= dim_B(1) .or.  dim_A(2) /= dim_B(2)  ) then
      write (*,*) 'inconsistent shape'
    else
      write (*,*) 'to be filled'
      do i = 1, dim_A(1)
        do j = 1, dim_A(2)
          if ( dabs(A(i,j)-B(i,j)) > threshold_compare ) then
            write (*,*) 'warning, A(i,j) != B(i,j)', i, j, A(i,j), B(i,j)
          end if
        end do
      end do
    end if

  end subroutine


  subroutine check_detl_frozen(iarr,res)
    integer  ::  iarr(nel)
    integer  ::  i, counter_fc
    integer  ::  res

  ! flag/res = 1 means violation of frozen orbital condition
    flag_fc = 0
    counter_fc = 0

    do i = 1, nel
       if ( iarr(i) <= nrorb_fc*2) then
         counter_fc = counter_fc +1
       end if

    end do

    if (counter_fc .ne. nrorb_fc*2 )  then
       flag_fc = 1
    end if

    if (flag_fc == 0) then
      res = 0
    else
      res = 1
    end if

  end subroutine

  subroutine check_shdl_frozen(iarr,res)
    integer  ::  iarr(nel-1)
    integer  ::  i, counter_fc
    integer  ::  res

  ! flag/res = 1 means violation of frozen orbital condition
    flag_fc = 0
    counter_fc = 0

    do i = 1, nel-1
       if ( iarr(i) <= nrorb_fc*2) then
         counter_fc = counter_fc +1
       end if


    end do

    ! =nnumber of fc, =number fc-1 both ok
    ! not possible for counter_fc > nrorb_fc
    if (counter_fc < nrorb_fc*2 -1 )  then
       flag_fc = 1
    end if

    if (flag_fc == 0) then
      res = 0
    else
     ! write (*,*) 'warning detl fc inconsistent'
      res = 1
    end if

  end subroutine





  subroutine backup_finalpsi()
    integer :: n
    character(len=30) :: target_filename, backup_filename
    character(len=100) :: cp_command
  
    logical :: exists_file
  
    if (flag_spinorbital == 0) then
      target_filename = 'finalpsi'
    end if
    if (flag_spinorbital == 1) then
      target_filename = 'finalpsi_spinorbital'
    end if
  
    do n = 1, 255
      if (flag_spinorbital == 0) then
        write(backup_filename, '(a,i1)') 'finalpsi_', n
      end if
  
      if (flag_spinorbital == 1) then
        write(backup_filename, '(a,i1)') 'finalpsi_spinorbital_', n
      end if
                                      
      inquire(file = backup_filename, exist = exists_file)
    
      if (.not. exists_file) then
        write(cp_command, '(a,a,a,a)') 'cp ', target_filename, ' ', backup_filename
        call system(cp_command)
        exit
      end if  
    end do  
  
  end subroutine
  






  integer function checksz(iarr)
    implicit none
    integer :: iarr(nel), ie, cs

    cs = 0
    !cw> calculate total Sz by accounting odd/even orb number
    do ie=1, nel
      if (mod(iarr(ie), 2) == 0) then
        cs = cs - 1
      else
        cs = cs + 1
      end if
    enddo

    checksz = cs
  end function checksz

end module

!> @file
!> @brief contains utility functions.
