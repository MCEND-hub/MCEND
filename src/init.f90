!> Module containing the initialization of the MCEND wave function.
  module moreadin
    use params
    use globalvars
    use inputvars
    use combinations
    use omp_lib
    use utils

    implicit none

    integer :: numel
    character(len=300) :: iom
    character(len=255) :: workdir
    character(len=2) ::  grdnum
    integer(i64) :: total_dim
    integer :: oneet, p, q
    integer :: rhoss, hmflen
    integer :: hmflen_alpha, hmflen_beta
    integer :: hmflen_spinorbital
    integer :: hmflen_frzorb, hmflen_spinorbital_frzorb
    integer :: rhoss_alpha, rhoss_beta
    integer :: rhoss_spinorbital
    integer :: rhoss2, rhoss2_spinorbital
    integer :: flag_dhdl
    integer :: rhoss_frzorb,rhoss_frzorb_spinorbital

    contains

!> Read in the electronic integrals from the basis library.
      subroutine read_oldmethod()

        implicit none
        real(dp) ::  smoh(nrprime,nrprime) ! transformation matrix
        integer :: i, j, k, l
        character(len=3) ::  basisfn
        character(len=50) :: inpformat
        integer :: ios, iensp

        numel = nrprime*(nrprime+1)/2
        total_dim = 6_i64*(nrprime**2) + (nrprime**4)
        inpformat='(e23.16)'

        call get_environment_variable("PWD", workdir)
        path2ints = trim(workdir)//'/'//trim(intdir)//'/'//trim(cmpdname)//'-data2-'

        if (2*nrorb < nel) then
          write(*,*) "Not enough orbitals for HF! Aborting run..."
          stop
        endif

        ! find out which electronic basis to use
        write(grdnum,'(i2.2)') inigrd
        write(basisfn,'(i3.3)') nrprime

        select case(basis)
          case(1)
            write(*,*) "Using precomputed basis at ",nrensp," nuclear positions"
          case(2)
            write(*,*) "Using atom-centered basis updated for each nuclear position"
        end select


          open(24, status='old', file=trim(path2ints)//grdnum//'.dat', iostat=ios, iomsg=iom)
            if (ios /= 0) then
              write(*,*) 'Fatal error! could not open integral files'
              write(*,*) 'Iomessage: ', trim(iom)
              stop
            endif
          do i=1,nrprime
            do j=1,nrprime
               !> hmat(i,j,#) at whatever point we select.
              read(24,*) hmat(i,j,inigrd)
            enddo
          enddo

          ! electronic kinetic energy
          do i=1,nrprime
            do j=1,nrprime
              read(24,*) tmat(i,j)
              if (basis == 2) taux(i,j,inigrd) = tmat(i,j)
            enddo
          enddo

          ! electronic position
          do i=1,nrprime
            do j=1,nrprime
              read(24,*) x(1,i,j)
            enddo
          enddo

          do i=1,nrprime
            do j=1,nrprime
              read(24,*) x(2,i,j)
            enddo
          enddo

          do i=1,nrprime
            do j=1,nrprime
              read(24,*) x(3,i,j)
            enddo
          enddo

          ! Overlap matrix
          do i=1,nrprime
            do j=1,nrprime
              read(24,*) smoh(i,j)
            enddo
          enddo

          ! two-electron repulsion
          do i=1,nrprime
            do j=1,nrprime
              do k=1,nrprime
                do l=1,nrprime
                  read(24,*) veeao(i,j,k,l)
                enddo
              enddo
            enddo
          enddo
          close(24)

        do iensp=1, nrensp
          if (iensp /= inigrd) then
            write(grdnum,'(i2.2)') iensp
              open(24, status='old', file=trim(path2ints)//grdnum//'.dat', iostat=ios, iomsg=iom)
              do i=1, nrprime
                do j=1, nrprime
                  read(24,*) hmat(i,j,iensp)
                enddo
              enddo

              if (basis == 2) then
                do i=1, nrprime
                  do j=1, nrprime
                    read(24,*) taux(i,j,iensp)
                  enddo
                enddo
              endif
          endif
          close(24)
        enddo

        do i = 1, nrprime
          hmat(i,i,:) = hmat(i,i,:) + shifte/dble(nel)
        enddo

      end subroutine

      !> Reads in integrals from Psi4 instead of gamess.
      subroutine read_newmethod()
        implicit none
        real(dp), allocatable   :: h(:), t(:), s(:,:), rij(:)
        integer,  allocatable   :: rijindex(:,:)
        integer :: aoints, hmatf
        integer :: ios, ix
        integer :: i1, i2, i3, i4, hh
        integer :: i, j, k
        character(256) :: iom
        character(256) :: xerints
        character(256) :: h_vmatf
        character(356) :: path2ints1
        character(356) :: path2ints2

        call get_environment_variable("PWD", workdir)

        select case(basis)
          case(1)
            write(*,*) "Using precomputed basis at ",nrensp," nuclear positions"
          case(2)
            write(*,*) "Using atom-centered basis updated for each nuclear position"
        end select

        numel = nrprime*(nrprime + 1)/2

        write(grdnum,'(i2.2)') nrensp
        h_vmatf = trim(cmpdname)//'-hmat-r01-'//grdnum//'.dat'
        xerints = trim(cmpdname)//'-ints-01.dat'
        path2ints1 = trim(intdir)//'/'//trim(h_vmatf)
        path2ints2 = trim(intdir)//'/'//trim(xerints)
      !  write (*,*) "path2ints1", path2ints1
        open(newunit=hmatf, file=trim(path2ints1), status='old', iostat=ios, iomsg=iom)
        if (ios /= 0) then
          write(*,*) 'Fatal error, could not open '//trim(path2ints1)//', iomsg:'
          write(*,*) trim(iom)
          stop
        endif

        if (.not. allocated(h)) allocate(h(nrensp*numel))
        if (.not. allocated(t)) allocate(t(nrensp*numel))


  ! h under ao basis with grid from symmetric orthogonalized basis
        do i = 1, nrensp*numel
          read(hmatf,*) h(i), t(i)
        enddo

        k = 1
        !> number of explicit nuclear sampling points
        do ix=1, nrensp
        !> nrprime: number of electronic basis functions, and no of electrons
          do i=1, nrprime
            do j=1, i
              hmat(i,j,ix) = h(k)
              hmat(j,i,ix) = h(k)
              taux(i,j,ix) = t(k)
              taux(j,i,ix) = t(k)

              k = k + 1
            enddo
          enddo
        enddo

        !> shift energy closer to zero to increase numerical stability
        !> make sure loop only goes over only the diagonal
        !write (*,*) 'shifte', shifte 
        do i=1, nrprime
          hmat(i,i,:) = hmat(i,i,:) + shifte/dble(nel)
        enddo

        ! inigrd: starting point on nuclear grid
        tmat(:,:) = taux(:,:,inigrd)
        close(hmatf)

        !write (*,*) 'read step 2 path2ints2', path2ints2
        open(newunit=aoints, file=trim(path2ints2), status='old', iostat=ios, iomsg=iom)
        if (ios /= 0) then
          write(*,*) 'Fatal error, could not open '//trim(path2ints2)//', iomsg:'
          write(*,*) trim(iom)
          stop
        endif

        if (.not. allocated(s)) allocate(s(3,numel))
        if (.not. allocated(rijindex)) allocate(rijindex(4,nint2e))
        if (.not. allocated(rij)) allocate(rij(nint2e))


        ! electronic coordinate
        do i = 1, numel
          read(aoints,*) s(1,i)
        enddo
        do i = 1, numel
          read(aoints,*) s(2,i)
        enddo
        do i = 1, numel
          read(aoints,*) s(3,i)
        enddo

        do i = 1, nint2e
          read(aoints,*) rijindex(1,i), rijindex(2,i), &
                         rijindex(3,i), rijindex(4,i), rij(i)
        enddo

        close(aoints)

        ! rebuild
        k = 1
        do i = 1, nrprime
          do j = 1, i
            x(1,i,j) = s(1,k)
            x(1,j,i) = s(1,k)

            x(2,i,j) = s(2,k)
            x(2,j,i) = s(2,k)

            x(3,i,j) = s(3,k)
            x(3,j,i) = s(3,k)

            k = k + 1
          enddo
        enddo

        veeao(:,:,:,:) = 0.0_dp
        ! change to the physicists' notation
        do i=1, nint2e
          hh = rijindex(2,i)
          rijindex(2,i) = rijindex(3,i)
          rijindex(3,i) = hh
        enddo

        do i=1, nint2e

          i1 = rijindex(1,i)
          i2 = rijindex(2,i)
          i3 = rijindex(3,i)
          i4 = rijindex(4,i)
          veeao(i1,i2,i3,i4) = rij(i)
          veeao(i1,i4,i3,i2) = rij(i)
          veeao(i3,i2,i1,i4) = rij(i)
          veeao(i3,i4,i1,i2) = rij(i)
          veeao(i2,i1,i4,i3) = rij(i)
          veeao(i2,i3,i4,i1) = rij(i)
          veeao(i4,i1,i2,i3) = rij(i)
          veeao(i4,i3,i2,i1) = rij(i)

        enddo

        deallocate(h,t,s,rijindex,rij)

      end subroutine

!> This subroutine initializes the wave functions, and constructs the hash-tables
!! for the evaluation of Slater-Condon rules.
      subroutine initstuff()

      implicit none
        real(dp) ::  dk
        integer :: iarr(nel)
        integer :: iarr_spinorbital(nel), iarr3_spinorbital(nel)
        integer :: iarr_shdl(nel-1)
        integer :: i, j, k
        integer :: counter
        integer :: counter_spinorbital
        integer :: counter_frzorb,counter_frzorb_spinorbital
        integer :: counter2
        integer :: counter2_spinorbital
        integer :: ls
        integer :: iarr_temp(nel)
        integer :: nexcit
        integer :: nexcit_spinorbital
        integer :: nexcit2
        integer :: nexcit2_spinorbital
        integer :: vorz, ind, lr, ie
        integer :: vorz2, ind2, cis, js ! bookkeeping vectors
        integer ::  cid
        integer :: iarr3(nel)
        integer :: ms, ns, ps, qs
        integer :: nrdiff, jr, ds
        integer :: jshf, lshf, ind1, vorz1
        integer :: jshf_spinorbital
        integer :: rsp_f
        integer :: allowlist(2*nrorb-(nel-1)), nrallow, iensp, jn
        integer :: nrallow_spinorbital
        integer :: allowlist_spinorbital(nrorb_spinorbital-(nel-1))
        integer :: hmfcounter, hmfcounter2
        integer :: hmfcounter_spinorbital, hmfcounter2_spinorbital
        integer :: flag_frzorb
        character(len=10) :: ncol
        character(len=1)  :: ncol2


        do i=1, nrprimn
          r(i) = rmin + (i - 1)*dr
        enddo

        ! electron nuclear sampling points -> all
        select case(basis)
          case(1)
            open(newunit=rsp_f, file=rspfile)
            do iensp=1, nrensp
              read(rsp_f,*) rsp(iensp)
            enddo
            close(rsp_f)
          case(2)
            do i=1, nrensp
              rsp(i) = r(i)
            enddo
        end select

        ! nrprimn # of nuclear grids
        ! nuclear position
        dk = 2.0_dp*pi/((nrprimn - 1)*dr)
        do i=1, nrprimn/2 + 1
          kx(i) = dble(i - 1)*dk
        enddo

        do i=nrprimn/2 + 2, nrprimn
          kx(i) = dble(i - 1 - nrprimn)*dk
        enddo


        if (basis == 1) call akimpotintvenao()


        ! complex absorbing potential
        vcap(1:nrprime,1:nrprime) = 0.0_dp
        !> right now this does nothing since etacap is hardcoded to be zero
        if (etacap > 0.0_dp) then
          do i=naocap, nrprime
            do j=naocap, nrprime
              vcap(i,j) = etacap
            enddo
          enddo

        endif


        if (use_wcap) call getnuclcap()
        ncol = stri(nel)
        nrorb_alpha = nrorb
        nrorb_beta  = nrorb



        if (flag_spinorbital == 0) then
          call comb(2*nrorb, nel, rcom)

          counter = 0
          do i=0, choose(2*nrorb, nel) - 1
             do j=1, nel
               counter = counter + 1
               k = nel - j
             ! k went through the number of electrons
             enddo
             deallocate(rcom(i)%combs)
          enddo

          call comb(2*nrorb, nel, rcom)

          counter_frzorb = 0
          do i=0, choose(2*nrorb, nel) - 1
             do j=1, nel
               k = nel - j
             ! k went through the number of electrons
               iarr_temp(j) = rcom(i)%combs(k)
             enddo

             call check_detl_frozen(iarr_temp,flag_frzorb)

             if (flag_frzorb == 0) then
               do j=1, nel
                 counter_frzorb = counter_frzorb + 1
                 k = nel - j
             ! k went through the number of electrons
                 detl(counter_frzorb) = rcom(i)%combs(k)
               enddo
             end if

             deallocate(rcom(i)%combs)
          enddo

          deallocate(rcom)

          do i=1, nrindep
            do j=1, nel
              iarr(j) = detl((i-1)*nel+j)
            enddo

            hval3(i) = checksum(iarr,2*nrorb)

          enddo

          open (23,file="detl.ij")
          write(23,'(a6,1x,a'//stri(3*nel)//',2x,a8)') "Number", "Index", "Checksum"

          do i=1, nrindep
            write(23,'(i6,1x,'//stri(nel)//'(i3),2x,i9)') i, (detl((i-1)*nel+j),j=1,nel), hval3(i)
          enddo
          close (23)
        end if

        if (flag_spinorbital == 1) then
          ! write (*,*) 'comb spinorbital nrorb_spinorbital', nrorb_spinorbital
          call comb(nrorb_spinorbital, nel, rcom_spinorbital)
          counter_spinorbital = 0
          counter_frzorb_spinorbital = 0
          do i=0, choose(nrorb_spinorbital, nel) - 1
            do j=1, nel
              counter_spinorbital = counter_spinorbital + 1
              k = nel - j
              iarr_temp(j) = rcom_spinorbital(i)%combs(k)
              call check_detl_frozen(iarr_temp,flag_frzorb)
            enddo
            if ( flag_frzorb == 0 ) then
              do j=1, nel
                counter_frzorb_spinorbital = counter_frzorb_spinorbital + 1
                k = nel - j
                detl_spinorbital(counter_frzorb_spinorbital) = rcom_spinorbital(i)%combs(k)
              end do
            end if

            deallocate(rcom_spinorbital(i)%combs)
          enddo
          deallocate(rcom_spinorbital)
          ! write (*,*) 'nrindep_frzorb_spinorbital', nrindep_frzorb_spinorbital

          do i=1, nrindep_spinorbital
            do j=1, nel_spinorbital
              iarr_spinorbital(j) = detl_spinorbital((i-1)*nel_spinorbital+j)
            enddo
            hval3_spinorbital(i) = checksum(iarr_spinorbital, nrorb_spinorbital)
          enddo
          open (223,file="detl_spinorbital.ij")
          write(223,'(a6,1x,a'//stri(3*nel_spinorbital)//',2x,a8)') "Number", "Index", "Checksum"
          do i=1, nrindep_spinorbital
            write(223,'(i6,1x,'//stri(nel_spinorbital)//'(i3),2x,i9)') i, &
  & (detl_spinorbital((i-1)*nel_spinorbital+j),j=1,nel_spinorbital), hval3_spinorbital(i)
          enddo
          close (223)
        end if

        if (flag_spinorbital == 0) then
          ncol = stri(nel-1)
          nexcit = nel - 1
          call comb(2*nrorb, nexcit, rcom)
          counter = 0
          do i=0, choose(2*nrorb, nexcit) - 1
            do j=1, nexcit
              counter = counter + 1
              k = nexcit - j
            enddo
          enddo

          counter_frzorb = 0
          do i=0, choose(2*nrorb, nexcit) - 1
            do j=1, nexcit
              k = nexcit - j
              iarr_shdl(j) = rcom(i)%combs(k)
            enddo
            call check_shdl_frozen(iarr_shdl, flag_frzorb)
            if (flag_frzorb == 0) then
              do j=1, nexcit
                counter_frzorb = counter_frzorb + 1
                k = nexcit - j
                shdl(counter_frzorb) = rcom(i)%combs(k)
              end do
            end if
            deallocate(rcom(i)%combs)
          enddo
          deallocate(rcom)
          counter2 = 0
          nexcit2  = nel - 2
          if (nexcit2>=0) then
            call comb(2*nrorb, nexcit2, rcom2)
            do i=0, choose(2*nrorb, nexcit2) - 1
              do j=1, nexcit2
                counter2 = counter2 + 1
                k = nexcit2 - j
                dhdl(counter2) = rcom2(i)%combs(k)
              enddo
              deallocate(rcom2(i)%combs)
            enddo
            deallocate(rcom2)
          end if
        end if


        if (flag_spinorbital == 1) then
          nexcit_spinorbital = nel - 1
          call comb(nrorb_spinorbital, nexcit_spinorbital, rcom_spinorbital)
          counter2_spinorbital = 0
          nexcit2_spinorbital  = nel - 2
          counter = 0
          do i=0, choose(nrorb_spinorbital, nexcit_spinorbital) - 1
            do j=1, nexcit_spinorbital
              counter = counter + 1
              k = nexcit_spinorbital - j
            enddo
          enddo

          counter_frzorb = 0
          do i=0, choose(nrorb_spinorbital, nexcit_spinorbital) - 1
            do j=1, nexcit_spinorbital
              k = nexcit_spinorbital - j
              iarr_shdl(j) = rcom_spinorbital(i)%combs(k)
            enddo

            call check_shdl_frozen(iarr_shdl, flag_frzorb)
            if (flag_frzorb == 0) then
              do j=1, nexcit_spinorbital
                counter_frzorb = counter_frzorb + 1
                k = nexcit_spinorbital - j
                shdl_spinorbital(counter_frzorb) = rcom_spinorbital(i)%combs(k)
              end do
            end if

            deallocate(rcom_spinorbital(i)%combs)
          enddo
          deallocate(rcom_spinorbital)

          if (nexcit2_spinorbital>=0) then
            call comb(nrorb_spinorbital, nexcit2_spinorbital, rcom2_spinorbital)
            do i=0, choose(nrorb_spinorbital, nexcit2_spinorbital) - 1
              do j=1, nexcit2_spinorbital
                counter2_spinorbital = counter2_spinorbital + 1
                k = nexcit2_spinorbital - j
                dhdl_spinorbital(counter2_spinorbital) = rcom2_spinorbital(i)%combs(k)
              enddo
              deallocate(rcom2_spinorbital(i)%combs)
            enddo
            deallocate(rcom2_spinorbital)
          end if
        end if

        if (flag_spinorbital == 0) then

          write(*,*) "shdl-size: ", nrshf, " / ", max_nrindep*(nel-1)
         ! write(*,*) "shdl determinants", size(shdl) , counter
          write(*,*) "max_nrindep", max_nrindep
          if (nel > 1) then
            open(23,file="shdl.ij")
            do i=1, max_nrindep
                write(23,'(i5,'//stri(nel-1)//'(i3))') i, (shdl((i-1)*(nel-1)+j),j=1,nel-1)
            enddo
            close (23)
          end if

          if (nel > 2) then
            open(23,file="dhdl.ij")
            do i=1, max_nrindep_2
              write(23,'(i5,'//stri(nel-2)//'(i3))') i, (dhdl((i-1)*(nel-2)+j),j=1,nel-2)
            enddo
            close (23)
          end if
        end if

        if (flag_spinorbital == 1) then
          if ( nel > 1 ) then
            open(23,file="shdl_spinorbital.ij")
            do i=1, max_nrindep_spinorbital
              write(23,'(i5,'//stri(nel-1)//'(i3))') i, (shdl_spinorbital((i-1)*(nel-1)+j),j=1,nel-1 )
            enddo
            close (23)
          end if

          if (nel > 2) then
            open(23,file="dhdl_spinorbital.ij")
            do i=1, max_nrindep_2_spinorbital
              write(23,'(i5,'//stri(nel-2)//'(i3))') i, &
& (dhdl_spinorbital((i-1)*(nel-2)+j),j=1,nel-2)
            enddo
            close (23)
          end if
        end if

        if (flag_spinorbital == 0) then
          n_detl = 0
          open(23,file="allow.list")
          do jshf=1, nrshf
            call getallowed1st(jshf,nrallow,allowlist, nrorb*2, iarr, max_nrindep*(nel-1), shdl, nrindep, hval3)
            allow1(jshf,0) = nrallow
            call int2str(ncol2,nrallow)


            if (nrallow /= 0) then
              write(23,'(i5,'//ncol2//'(i3,2x))') jshf, (allowlist(i),i=1,nrallow)
            endif
            do i=1, nrallow
              allow1(jshf,3*(i-1)+1) = allowlist(i)
            enddo
            do i=1, nrallow
              iarr(1) = allow1(jshf,3*(i-1)+1)
              do ie=1, nel-1
                iarr(1+ie) = shdl((jshf-1)*(nel-1)+ie)
              enddo
              call getindex(iarr,ind,vorz, nrorb*2, hval3, nrindep )
              allow1(jshf,3*(i-1)+2) = ind
              allow1(jshf,3*(i-1)+3) = vorz
            enddo
          enddo
          close(23)
         end if

        if (flag_spinorbital == 1) then
          open(230,file="allow_spinorbital.list")
          do jshf_spinorbital=1, nrshf_spinorbital
            call getallowed1st(jshf_spinorbital,nrallow_spinorbital,allowlist_spinorbital, nrorb_spinorbital,&
          &iarr_spinorbital, max_nrindep_spinorbital*(nel-1), shdl_spinorbital, nrindep_spinorbital, &
          &hval3_spinorbital)
            allow1_spinorbital(jshf_spinorbital,0) = nrallow_spinorbital
            call int2str(ncol2,nrallow_spinorbital)
            if (nrallow_spinorbital /= 0) then
              write(230,'(i5,'//ncol2//'(i3,2x))') jshf_spinorbital, (allowlist_spinorbital(i),i=1,nrallow_spinorbital)
            else
            endif
            do i=1, nrallow_spinorbital
              allow1_spinorbital(jshf_spinorbital,3*(i-1)+1) = allowlist_spinorbital(i)
            enddo
            do i=1, nrallow_spinorbital
              iarr_spinorbital(1) = allow1_spinorbital(jshf_spinorbital,3*(i-1)+1)
              do ie=1, nel-1
                iarr_spinorbital(1+ie) = shdl_spinorbital((jshf_spinorbital-1)*(nel-1)+ie)
              enddo
              call getindex(iarr_spinorbital,ind,vorz, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
              allow1_spinorbital(jshf_spinorbital,3*(i-1)+2) = ind
              allow1_spinorbital(jshf_spinorbital,3*(i-1)+3) = vorz
            enddo
          enddo
          close(230)
        end if

        flag_dhdl = 1
        if (flag_spinorbital == 0) then
          counter = 1
          do js=1, 2*nrorb
            do ls=1, 2*nrorb
              do cis=1, nrshf
                iarr(1)  = js
                iarr3(1) = ls
                do i=2, nel
                  iarr(i) = shdl((cis-1)*(nel-1)+i-1)
                  iarr3(i) = iarr(i)
                enddo
                call getindex2(iarr,ind,vorz, nrorb*2, hval3, nrindep)
                call getindex2(iarr3,ind2,vorz2, nrorb*2, hval3, nrindep)
                if ((vorz /= 0) .and. (vorz2 /= 0)) then
                ! yes, these two determinants are valid and have one single SO difference
                ! between them
                  hvalrho(counter) = vorz*vorz2
                  hvalrho(counter+1) = ind
                  hvalrho(counter+2) = ind2
                  hvalrho(counter+3) = js
                  hvalrho(counter+4) = ls
                  counter = counter + 5
                endif
              enddo
            enddo
          enddo

        ! actual number of matrix elements between determinants - one SO difference
          write(*,*) "hvalrho-size:", counter - 1, " / ", 20*nrorb*nrorb*nrindep
          rhoss = (counter - 1)/5
          write(*,*) "rhoss:", rhoss
          flag_dhdl = 1

  ! rdm2
  ! seems got error in frozen core, block it for a while
          if (flag_dhdl == 1) then
            counter2 = 1
            do js=1, 2*nrorb
              do ls=1, 2*nrorb
                do ms=1, 2*nrorb
                  do ns=1, 2*nrorb
                    do cid=1, nrdhf
                      iarr(1)  = js
                      iarr(2)  = ls
                      iarr3(1) = ms
                      iarr3(2) = ns
                      do i=3, nel
                        iarr(i) = dhdl((cid-1)*(nel-2)+i-2)
                        iarr3(i) = iarr(i)
                      enddo
                      call getindex2(iarr,ind,vorz, nrorb*2, hval3, nrindep )
                      call getindex2(iarr3,ind2,vorz2, nrorb*2, hval3, nrindep )
                      if ((vorz /= 0) .and. (vorz2 /= 0)   ) then
                ! yes, these two determinants are valid and have one single SO difference
                ! between them
                        hvalrdm2(counter2)   = vorz*vorz2
                        hvalrdm2(counter2+1) = ind
                        hvalrdm2(counter2+2) = ind2
                        hvalrdm2(counter2+3) = js
                        hvalrdm2(counter2+4) = ls
                        hvalrdm2(counter2+5) = ms
                        hvalrdm2(counter2+6) = ns
                        counter2 = counter2 + 7
                      endif
                    end do
                  end do
                end do
              end do
            end do
            write(*,*) "hvalrdm2-size:", counter2 - 1, " / ", 20*nrorb*nrorb*nrorb*nrorb*nrindep
            rhoss2 = (counter2 - 1)/7
            write(*,*) "rhoss2:", rhoss2
          end if
       end if

       if (flag_spinorbital == 1) then
         if (flag_dhdl == 1) then
           counter2_spinorbital = 1
           do js=1, nrorb_spinorbital
             do ls=1, nrorb_spinorbital
               do ms=1, nrorb_spinorbital
                 do ns=1, nrorb_spinorbital
                   do cid=1, nrdhf_spinorbital
                     iarr_spinorbital(1)  = js
                     iarr_spinorbital(2)  = ls
                     iarr3_spinorbital(1) = ms
                     iarr3_spinorbital(2) = ns
                     do i=3, nel
                       iarr_spinorbital(i) = dhdl_spinorbital((cid-1)*(nel-2)+i-2)
                       iarr3_spinorbital(i) = iarr_spinorbital(i)
                     enddo
                     call getindex2(iarr_spinorbital,ind,vorz, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital )
                     call getindex2(iarr3_spinorbital,ind2,vorz2, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital )
                     if ((vorz /= 0) .and. (vorz2 /= 0)   ) then
                ! yes, these two determinants are valid and have one single SO difference
                ! between them
                       hvalrdm2_spinorbital(counter2_spinorbital)   = vorz*vorz2
                       hvalrdm2_spinorbital(counter2_spinorbital+1) = ind
                       hvalrdm2_spinorbital(counter2_spinorbital+2) = ind2
                       hvalrdm2_spinorbital(counter2_spinorbital+3) = js
                       hvalrdm2_spinorbital(counter2_spinorbital+4) = ls
                       hvalrdm2_spinorbital(counter2_spinorbital+5) = ms
                       hvalrdm2_spinorbital(counter2_spinorbital+6) = ns
                       counter2_spinorbital = counter2_spinorbital + 7
                     endif
                   end do
                 end do
               end do
             end do
           end do

           write(*,*) "hvalrdm2_spinorbital-size:", counter2_spinorbital - 1, " / ", &
           & 20*nrorb_spinorbital*nrorb_spinorbital*nrorb_spinorbital*nrorb_spinorbital*nrindep_spinorbital

           rhoss2_spinorbital = (counter2_spinorbital - 1)/7
           write(*,*) "rhoss2_spinorbital:", rhoss2_spinorbital

         end if

         counter_spinorbital = 1
         do js=1, nrorb_spinorbital
           do ls=1, nrorb_spinorbital
             do cis=1, nrshf_spinorbital
               iarr_spinorbital(1)  = js
               iarr3_spinorbital(1) = ls
               do i=2,nel_spinorbital
                 iarr_spinorbital(i)  = shdl_spinorbital((cis-1)*(nel_spinorbital-1)+i-1)
                 iarr3_spinorbital(i) = iarr_spinorbital(i)
               end do

               call getindex2(iarr_spinorbital,ind,vorz, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
               call getindex2(iarr3_spinorbital,ind2,vorz2, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)

               if ((vorz /= 0) .and. (vorz2 /= 0)) then
                 hvalrho_spinorbital(counter_spinorbital) = vorz*vorz2
                 hvalrho_spinorbital(counter_spinorbital+1) = ind
                 hvalrho_spinorbital(counter_spinorbital+2) = ind2
                 hvalrho_spinorbital(counter_spinorbital+3) = js
                 hvalrho_spinorbital(counter_spinorbital+4) = ls
                 counter_spinorbital = counter_spinorbital + 5
               endif

             end do
           end do
         end do

         write(*,*) "hvalrho_spinorbital-size:", counter - 1, " / ", 20*nrorb_spinorbital*nrorb_spinorbital*nrindep_spinorbital
         rhoss_spinorbital = (counter_spinorbital - 1)/5
         ! write(*,*) "rhoss_spinorbital:", rhoss_spinorbital
         ! write (*,*) 'compare hvalrho and hval_spinorbital'
        !write (*,*)  hvalrho
         ! write (*,*) 'so'

       end if

        if (flag_spinorbital == 0) then
          hmfcounter = 0
          hmfcounter2 = 0
          do jn=1, nrspf
            do jr=1, nrorb
              do lr=1, nrorb
                do jshf=1, nrshf
                  do lshf=1, nrshf
                  iarr(1) = 2*jr
                    do ie=2, nel
                      iarr(ie) = shdl((jshf-1)*(nel-1)+(ie-1))
                    enddo
                    call getindex2(iarr,ind1,vorz1, nrorb*2, hval3, nrindep)
                    iarr(1) = 2*lr
                    do ie=2, nel
                      iarr(ie) = shdl((lshf-1)*(nel-1)+(ie-1))
                    enddo
                  call getindex2(iarr,ind2,vorz2, nrorb*2, hval3, nrindep)
                   call maxcoincshf(jshf,lshf,vorz,nrdiff,ms,ns,ps,qs, shdl, max_nrindep)
                    ds = deltaspin(ms,ps)
                    if ((nrdiff == 0) .and. (vorz1*vorz2 /= 0)) then
                      hmfcounter = hmfcounter + 7
                    endif

                    if ((nrdiff == 1) .and. (vorz1*vorz2*ds /= 0)) then
                      hmfcounter2 = hmfcounter2 + 8
                    endif

                  enddo
                enddo
              enddo
            enddo
          enddo
        end if


        if (flag_spinorbital == 1) then

        hmfcounter_spinorbital = 0
        hmfcounter2_spinorbital = 0

        write (*,*) 'shdl_spinorbital'
        do i = 1, size(shdl_spinorbital)
        end do

        do jn=1, nrspf
          do jr=1, nrorb_spinorbital
            do lr=1, nrorb_spinorbital
              do jshf=1, nrshf_spinorbital
                do lshf=1, nrshf_spinorbital
                  iarr(1) = 1*jr
                  do ie=2, nel_spinorbital
                    iarr(ie) = shdl_spinorbital((jshf-1)*(nel_spinorbital-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind1,vorz1, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
                  iarr(1) = 1*lr
                  do ie=2, nel_spinorbital
                    iarr(ie) = shdl_spinorbital((lshf-1)*(nel_spinorbital-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind2,vorz2, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
                  call maxcoincshf(jshf,lshf,vorz,nrdiff,ms,ns,ps,qs, shdl_spinorbital, max_nrindep_spinorbital)
                  ds = deltaspin(ms,ps)
                  if ((nrdiff == 0) .and. (vorz1*vorz2 /= 0)) then
                    hmfcounter_spinorbital = hmfcounter_spinorbital + 7
                  endif
                  if ((nrdiff == 1) .and. (vorz1*vorz2*ds /= 0)) then
                    hmfcounter2_spinorbital = hmfcounter2_spinorbital + 8
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        endif

        if (flag_spinorbital == 0) then
        hmflen = max(hmfcounter, hmfcounter2)
        write(*,*) 'Estimated HMF matrix size: ', hmflen
        allocate(hmf(2,0:hmflen))
        hmf(:,:) = 0
        do jn=1, nrspf
          do jr=1, nrorb
            do lr=1, nrorb
              do jshf=1, nrshf
                do lshf=1, nrshf
                  iarr(1) = 2*jr
                  do ie=2, nel
                    iarr(ie) = shdl((jshf-1)*(nel-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind1,vorz1, nrorb*2, hval3, nrindep)
                  iarr(1) = 2*lr
                  do ie=2, nel
                    iarr(ie) = shdl((lshf-1)*(nel-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind2,vorz2, nrorb*2, hval3, nrindep)
                  !> now find out the maximum coincidence between those two determinants-1:
                  !> maximum coincidence between the single-hole functions
                  call maxcoincshf(jshf,lshf,vorz,nrdiff,ms,ns,ps,qs, shdl, max_nrindep)
                  !> jshf and lshf are the number of single-hole function
                  !> shdl contains all the single-hole functions
                  !> vorz tells if the sign changed due to a permutation of an SO
                  !> nrdiff tells how many SOs are different between the SHF: 0 up to nel-1
                  !> ms and ps are the indices of the first two different SOs
                  !> ns and qs are the indices of the second two different SOs
                  !> so the excitation goes: ms-->ps and ns-->qs
                  !> check that there is no spin flip: ds
                  ds = deltaspin(ms,ps)
                  !> everything is saved IF the determinants are allowed and
                  !> the single-hole part is equal = single excitation of the determinant
                  !> 6 different entires, spin, sign, mo index for 1 and 2
                  if ((nrdiff == 0) .and. (vorz1*vorz2 /= 0)) then
                    hmf(1,0) = hmf(1,0) + 1
                    if (6*hmf(1,0) > hmflen) then
                      write(*,*) "hmfsize1 too small!"
                      stop
                    endif
                    ind1 = (jn - 1)*nrindep + ind1
                    ind2 = (jn - 1)*nrindep + ind2
                    ! SHF is identical so only one index needs to be saved
                    hmf(1,6*(hmf(1,0)-1)+1) = jshf
                    ! MOs are saved in jr and lr
                    hmf(1,6*(hmf(1,0)-1)+2) = jr
                    hmf(1,6*(hmf(1,0)-1)+3) = lr
                    ! keep track of the sign - of all the permutations
                    ! would only need to be vorz1*vorz2 since for nrdiff=0 --> vorz=1
                    hmf(1,6*(hmf(1,0)-1)+4) = vorz*vorz1*vorz2
                    ! determinant number is saved in in1 and ind2
                    hmf(1,6*(hmf(1,0)-1)+5) = ind1
                    hmf(1,6*(hmf(1,0)-1)+6) = ind2
                  endif

                  ! now save if the determinants are allowed and there is one
                  ! difference in the SHFs -> double excitation of the determinant
                  ! also make sure that there is no spin flip --> ds
                  if ((nrdiff == 1) .and. (vorz1*vorz2*ds /= 0 )) then
                    hmf(2,0) = hmf(2,0) + 1
                    if (7*hmf(2,0) > hmflen) then
                      write(*,*) "hmfsize2 too small!"
                      stop
                    endif
                    ! push these indices up if there is more than one single-particle
                    ! function for the nuclei
                    ind1 = (jn - 1)*nrindep + ind1
                    ind2 = (jn - 1)*nrindep + ind2
                    ! SO number in which SHFs differ is saved in ms and ps
                    hmf(2,7*(hmf(2,0)-1)+1) = ms
                    hmf(2,7*(hmf(2,0)-1)+2) = ps
                    ! keep track of the sign - of all the permutations
                    hmf(2,7*(hmf(2,0)-1)+3) = vorz*vorz1*vorz2
                    ! determinant number is saved in ind1 and ind2
                    hmf(2,7*(hmf(2,0)-1)+4) = ind1
                    hmf(2,7*(hmf(2,0)-1)+5) = ind2
                    ! MOs are saved in jr and lr
                    hmf(2,7*(hmf(2,0)-1)+6) = jr
                    hmf(2,7*(hmf(2,0)-1)+7) = lr
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo


        end if

        if (flag_spinorbital == 1) then

        hmflen_spinorbital = max(hmfcounter_spinorbital, hmfcounter2_spinorbital)
        write(*,*) 'Estimated HMF_spinorbital matrix size: ', hmflen_spinorbital
        allocate(hmf_spinorbital(2,0:hmflen_spinorbital))

        hmf_spinorbital(:,:) = 0
        do jn=1, nrspf
          do jr=1, nrorb_spinorbital
            do lr=1, nrorb_spinorbital
              do jshf=1, nrshf_spinorbital
                do lshf=1, nrshf_spinorbital
                  iarr(1) = 1*jr
                  do ie=2, nel_spinorbital
                    iarr(ie) = shdl_spinorbital((jshf-1)*(nel_spinorbital-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind1,vorz1, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
                  iarr(1) = 1*lr
                  do ie=2, nel_spinorbital
                    iarr(ie) = shdl_spinorbital((lshf-1)*(nel_spinorbital-1)+(ie-1))
                  enddo
                  call getindex2(iarr,ind2,vorz2, nrorb_spinorbital, hval3_spinorbital, nrindep_spinorbital)
                  call maxcoincshf(jshf,lshf,vorz,nrdiff,ms,ns,ps,qs, shdl_spinorbital, max_nrindep_spinorbital)
                  ds = deltaspin(ms,ps)
                  if ((nrdiff == 0) .and. (vorz1*vorz2 /= 0)) then
                    hmf_spinorbital(1,0) = hmf_spinorbital(1,0) + 1
                    if (6*hmf_spinorbital(1,0) > hmflen_spinorbital) then
                      write(*,*) "hmfsize1 too small!"
                      stop
                    endif
                    ! push these indices up if there is more than one single-particle
                    ! function for the nuclei
                    ind1 = (jn - 1)*nrindep_spinorbital + ind1
                    ind2 = (jn - 1)*nrindep_spinorbital + ind2
                    ! SHF is identical so only one index needs to be saved
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+1) = jshf
                    ! MOs are saved in jr and lr
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+2) = jr
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+3) = lr
                    ! keep track of the sign - of all the permutations
                    ! would only need to be vorz1*vorz2 since for nrdiff=0 --> vorz=1
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+4) = vorz*vorz1*vorz2
                    ! determinant number is saved in in1 and ind2
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+5) = ind1
                    hmf_spinorbital(1,6*(hmf_spinorbital(1,0)-1)+6) = ind2
                  endif
                  ! now save if the determinants are allowed and there is one
                  ! difference in the SHFs -> double excitation of the determinant
                  ! also make sure that there is no spin flip --> ds
                  if ((nrdiff == 1) .and. (vorz1*vorz2*ds /= 0 )) then
                    hmf_spinorbital(2,0) = hmf_spinorbital(2,0) + 1
                    if (7*hmf_spinorbital(2,0) > hmflen_spinorbital) then
                      write(*,*) "hmfsize2 too small!"
                      stop
                    endif
                    ind1 = (jn - 1)*nrindep_spinorbital + ind1
                    ind2 = (jn - 1)*nrindep_spinorbital + ind2
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+1) = ms
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+2) = ps
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+3) = vorz*vorz1*vorz2
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+4) = ind1
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+5) = ind2
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+6) = jr
                    hmf_spinorbital(2,7*(hmf_spinorbital(2,0)-1)+7) = lr
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo


        ! open(newunit=fhmf_spinorbital1,file="hmf_spinorbital_1.txt")
        ! write (fhmf_spinorbital1,*) hmf_spinorbital(1,:)
        ! close(fhmf_spinorbital1)
        ! open(newunit=fhmf_spinorbital2,file="hmf_spinorbital_2.txt")
        ! write (fhmf_spinorbital2,*) hmf_spinorbital(2,:)
        ! close(fhmf_spinorbital2)
        end if

        return
      end subroutine

!> Initialize electronic and nuclear wave functions. Different options for restart:
!! no restart; read in startpsi; read in startpsi and generate additional SPF.
      subroutine initwavefunc()

        implicit none
        complex(dp) :: olap
        complex(dp) :: auxphin(nrprimn), auxphin2(nrprimn), auxphin3(nrprimn)
        complex(dp), allocatable :: projP(:,:)
        real(dp)     :: work(3*nrprime), w(nrprime)
        real(dp)     :: workn(3*nrprimn), wn(nrprimn)
        real(dp)     :: norm
        real(dp)     :: norm_spinorbital
        real(dp)     :: h2(nrprime,nrprime)
        real(dp)     :: hnuc(nrprimn,nrprimn)
        real(dp)     :: rep(nrorb), imp(nrorb), rea, ima
        real(dp)     :: rep_spinorbital(nrorb_spinorbital),imp_spinorbital(nrorb_spinorbital), rea_spinorbital, ima_spinorbital
        integer(i64) :: planf, planb
        integer      :: mu, in, ix, jx
        integer      :: i, j
        integer      :: scfvals
        integer      :: fftw_forward=-1,fftw_backward=1
        integer      :: info, lwork
        integer      :: spsi
        integer      :: init_coef_A, sum_sz_tmp, detl_spinorbital_count, nel_count
        character(3) :: jobz, uplo, ncol
        character(3) :: nnn

        if (restart == 2 .or. restart == 5) then
          nnn = stri(2*nrorb)
          if (flag_spinorbital==0) then
          write(*,*) "Restart requested - reading from startpsi"
            open(newunit=spsi,file="startpsi")
            do i=1, nrindep*nrspf
              read(spsi,*) rea, ima
              A(i) = dcmplx(rea,ima)
            end do


            do mu=1, nrprime
              read(spsi,'('//stri(2*nrorb)//'(e23.16,1x))') (rep(in),in=1,nrorb), (imp(in),in=1,nrorb)
              do in=1, nrorb
                phi(mu,in) = dcmplx(rep(in),imp(in))
              enddo
            enddo

            write(ncol,'(i1)') 2*nrspf

            do ix=1, nrprimn
              read(spsi,'('//stri(2*nrspf)//'(e23.16,1x))') (rep(in),in=1,nrspf), (imp(in),in=1,nrspf)
              do in=1, nrspf
                phin(ix,in) = dcmplx(rep(in),imp(in))
              enddo
            enddo
            close(spsi)

          else if (flag_spinorbital==1) then
          write(*,*) "Restart requested - reading from startpsi_spinorbital"
            open(newunit=spsi,file="startpsi_spinorbital")
            do i=1, nrindep_spinorbital*nrspf
              read(spsi,*) rea_spinorbital, ima_spinorbital
              A_spinorbital(i) = dcmplx(rea_spinorbital,ima_spinorbital)
            end do
            do mu=1, nrprime
              read(spsi,'('//stri(nrorb_spinorbital)//'(e23.16,1x))') (rep_spinorbital(in),in=1,nrorb_spinorbital)
              read(spsi,'('//stri(nrorb_spinorbital)//'(e23.16,1x))') (imp_spinorbital(in),in=1,nrorb_spinorbital)
              do in=1, nrorb_spinorbital
                phi_spinorbital(mu,in) = dcmplx(rep_spinorbital(in),imp_spinorbital(in))
              enddo
            enddo

            write(ncol,'(i1)') 2*nrspf

            do ix=1, nrprimn
              read(spsi,'('//stri(2*nrspf)//'(e23.16,1x))') (rep_spinorbital(in),in=1,nrspf), (imp_spinorbital(in),in=1,nrspf)
              do in=1, nrspf
                phin(ix,in) = dcmplx(rep_spinorbital(in),imp_spinorbital(in))
              enddo
            enddo

            close(spsi)

          else
            write (*,*) 'flag_spinorbital not 0/1, supported'
            stop
          end if

        !> If this isn't a restart then initialize the wavefunction from the scf guess vals
        else if (restart == 0) then
          A(:) = c0
          A(1) = cr
          if (nsz == 0) then
            A(1) = cr
          end if

          A_spinorbital(:) = c0
          sum_sz_tmp = 0
          detl_spinorbital_count = 1
          nel_count = 0
          if (flag_spinorbital == 1) then
            do i = 1, nrindep_spinorbital*nel
               nel_count = nel_count +1
               if (mod(detl_spinorbital(i),2) == 0) then
                 sum_sz_tmp = sum_sz_tmp - 1
               else if (mod(detl_spinorbital(i),2) == 1) then
                 sum_sz_tmp = sum_sz_tmp + 1
               else
                 write (*,*) 'error in mod 2'
                 stop
               end if

               if (sum_sz_tmp == nsz .and. nel_count == nel)  then
                 init_coef_A = detl_spinorbital_count
                 exit
               end if

               if (mod(i,nel) == 0) then
                 sum_sz_tmp = 0
                 nel_count = 0
                 detl_spinorbital_count = detl_spinorbital_count + 1
               end if

            end do

            if (nsz == 0) then
              A_spinorbital(1) = cr
            else
               A_spinorbital(init_coef_A) = cr
            end if
          end if

          ! generalized option, read in SCF energies from run on each grid point
          open(newunit=scfvals, file=trim(scfv_path), status='old')
          ! in init, we setup electronic Hamiltonian, now we set up nuclear part
          do ix=1, nrprimn
            do jx=1, nrprimn
              hnuc(ix,jx) = 0.0_dp
            enddo
            read(scfvals,*) hnuc(ix,ix)
          enddo
          close(scfvals)


          
          if (mod(nrprimn, 2) == 0) then
                 
            !write (*,*) 'enter even branch'
            ! so far uniform nrprimn, convienent to define array
            do ix = 1, nrprimn
              do jx = 1, nrprimn
                if (ix == jx) then
            ! Tannor, quantum mechanics time-dependent perspective, p 307, 11.172
                  hnuc(ix, ix) = hnuc(ix, ix) + pi ** 2 * (nrprimn ** 2 + 2) & 
                  & / (6.0_dp * massn * nrprimn ** 2 * dr ** 2)
                  
                else
                  hnuc(ix, jx) = ((- 1.0_dp) ** (ix - jx)) * pi ** 2 & 
                  & / ((nrprimn * dr * dsin(dble(ix - jx) * pi / dble(nrprimn))) ** 2 * massn )
                end if
              end do
            end do 
            
          else                
           ! write (*,*) 'enter odd branch'
            do ix = 1, nrprimn
              do jx = 1, nrprimn
                if (ix == jx) then
            ! Tannor, quantum mechanics time-dependent perspective, p 307, 11.172
                  hnuc(ix, ix) = hnuc(ix, ix) + pi ** 2 * (nrprimn ** 2 + 1) &
                  & / (6.0_dp * massn * nrprimn ** 2 * dr ** 2)
                else
                  hnuc(ix, jx) = ((- 1.0_dp) ** (ix - jx)) * pi ** 2 * dcos( dble(ix - jx) * pi / dble(nrprimn) ) &
& / (nrprimn * dr)**2 /  massn
                end if
              end do
            end do   


            do ix = 1, nrprimn
              do jx = 1, nrprimn
                if (ix == jx) then
                   continue
                else
                   hnuc(ix, jx) = hnuc(ix, jx)  / ( dsin(dble(ix - jx) * pi /dble(nrprimn)) )**2
                end if
              end do
            end do




          end if

          jobz = "V"
          uplo = "U"
          lwork = 3*nrprimn
          info = 0
          call dsyev(jobz,uplo,nrprimn,hnuc,nrprimn,wn,workn,lwork,info)
          phin(:,:) = c0
          do ix=1, nrprimn
            do i=1, nrspf
              phin(ix,i) = -hnuc(ix,i)
              phin_spinorbital(ix,i) = -hnuc(ix,i)
            enddo
          enddo

          do i=1, nrspf
            norm = 0.0_dp
            norm_spinorbital = 0.0_dp
            do ix=1, nrprimn
              norm             = norm + dconjg(phin(ix,i))*phin(ix,i)
              norm_spinorbital = norm_spinorbital + dconjg(phin_spinorbital(ix,i))*phin_spinorbital(ix,i)
            enddo
            norm = dsqrt(norm*dr)
            norm_spinorbital = dsqrt(norm_spinorbital*dr)
            do ix=1, nrprimn
              phin(ix,i) = phin(ix,i)/norm
              phin_spinorbital(ix,i) = phin_spinorbital(ix,i)/norm_spinorbital
            enddo
          enddo
          call compare_scalar_real(norm,norm_spinorbital)

          lwork = 3*nrprime
          h2(:,:) = hmat(:,:,inigrd)

          call dsyev(jobz,uplo,nrprime,h2,nrprime,w,work,lwork,info)

          do mu=1, nrprime
            do i=1, nrorb
              !> I suppose if we wanted to keep virtuals we could compute a full nrprime x nrprime matrix
              !> or in the case we want to swap orbitals around, which might need to happen
              phi(mu,i) = h2(mu,i)
            enddo
          enddo

          nrorb_init_alpha = nrorb_alpha
          nrorb_init_beta  = nrorb_beta

          do mu=1, nrprime
            do i=1, nrorb_init_alpha
              phi_alpha(mu,i) = h2(mu,i)
            enddo
          enddo
          do mu=1, nrprime
            do i=1, nrorb_init_beta
              phi_beta(mu,i) = h2(mu,i)
            enddo
          enddo

          do mu=1, nrprime
            do i=1, nrorb_init_beta
               phi_spinorbital(mu,(2*i-1) )   = h2(mu,i)
               phi_spinorbital(mu,2*i)        = h2(mu,i)
            enddo
          enddo

    ! normalize
          do i=1,nrorb
            norm = 0.0_dp
            do mu=1,nrprime
              norm = norm + dconjg(phi(mu,i))*phi(mu,i)
            enddo
            norm = dsqrt(norm)
          end do

          do i=1,nrorb_spinorbital
            norm_spinorbital = 0.0_dp
            do mu=1,nrprime
              norm_spinorbital = norm_spinorbital + dconjg(phi_spinorbital(mu,i))*phi_spinorbital(mu,i)
            enddo
            norm_spinorbital = dsqrt(norm_spinorbital)
          end do
          call compare_scalar_real(norm,norm_spinorbital)

        else if (restart == 1) then

          call dfftw_plan_dft_1d(planf,nrprimn,auxphin,auxphin2,fftw_forward,0)
          call dfftw_plan_dft_1d(planb,nrprimn,auxphin2,auxphin3,fftw_backward,0)

          auxphin(:) = c0
          auxphin2(:) = c0
          auxphin3(:) = c0
          A(:) = c0
          A_spinorbital(:) = c0

          write(nnn,'(i2)') 2*nrorb
          open(newunit=spsi,file="startpsi")
          !> read in A elements from converged single SPF
          do i=1, nrindep*(nrspf-1)
            read(spsi,*) rea, ima
            A(i) = dcmplx(rea,ima)
            ! Mar-7: need to correct
            A_spinorbital(i) = dcmplx(rea,ima)
          enddo

          !> set the A elements for the second SPF to 0 (complex)
          do mu=1, nrprime
            read(spsi,'('//stri(2*nrorb)//'(e23.16,1x))') (rep(in),in=1,nrorb), (imp(in),in=1,nrorb)
            do in=1, nrorb
              phi(mu,in) = dcmplx(rep(in),imp(in))
            enddo
          enddo
          write(ncol,'(i1)') 2*(nrspf-1)
          !> fill the n-1 SPF(s) with values from previous calculation
          do ix=1, nrprimn
            read(spsi,'('//stri(2*(nrspf-1))//'(e23.16,1x))') &
                  (rep(in),in=1,nrspf-1), (imp(in),in=1,nrspf-1)
            do in=1, nrspf-1
              phin(ix,in) = dcmplx(rep(in),imp(in))
            enddo
          enddo
          close(spsi)

          !> prepare the projection operator
          if (.not. allocated(projP)) allocate(projP(nrprimn,nrprimn))

          projP(:,:) = c0
          !> create projection operator from phin of converged part
          do in=1, nrspf-1
            do ix=1, nrprimn
              do jx=1, nrprimn
                projP(ix,jx) = projP(ix,jx) + phin(ix,in)*dconjg(phin(jx,in))
              enddo
            enddo
          enddo
          do i=1, nrprimn
            auxphin(i) = phin(i,nrspf-1)
          enddo
          !> do a FFT on auxphin result is stored in auxphin2
          call dfftw_execute_dft(planf,auxphin,auxphin2)
          !> multiply it by the k-space factor
          auxphin2(:) = kx(:)*auxphin2(:)

          !> now do an inverst FFT of auxphin2 to get auxphin3
          call dfftw_execute_dft(planb,auxphin2,auxphin3)

          !> take the norm as the sum of \phi_n^* \cdot \phi_n
          norm = sum(dreal(dconjg(auxphin3(:))*auxphin3(:)))
          norm = dsqrt(norm*dr)
          phin(:,nrspf) = auxphin3(:)/norm

          !> ok, now do some stuff with the projection operator
          do i=1, nrprimn
            auxphin(i) = phin(i,nrspf)
            do j=1, nrprimn
              auxphin(i) = auxphin(i) - projP(i,j)*phin(j,nrspf)*dr
            enddo
          enddo

          deallocate(projP)

          !> do the same thing again with the normalization with whatever auxphin has become
          norm = sum(dreal(dconjg(auxphin(:))*auxphin(:)))
          norm = dsqrt(norm*dr)
          !> this is our projected initial phin
          phin(:,nrspf) = auxphin(:)/norm

          !> destroy all evidence
          call dfftw_destroy_plan(planf)
          call dfftw_destroy_plan(planb)

        endif

        if (flag_spinorbital == 0) then 
          open(20,file="olap.ij")
          write(20,*) 'Electronic phi overlap'
          write(20,'(a5,a5,a15)') 'MO', 'MO', 'Overlap'
          do i=1, nrorb
            do j=1, nrorb
              olap = c0
              do mu=1, nrprime
                olap = olap + dconjg(phi(mu,i))*phi(mu,j)
              enddo
              write(20,'(i5,i5,f15.7)') i, j, cdabs(olap)
            enddo
          enddo
        
         write(20,*) 'Nuclear phi overlap'
         write(20,'(a5,a5,a15)') 'SPF', 'SPF', 'Overlap'

         do i=1, nrspf
           do j=1, nrspf
             olap = c0
             do ix=1, nrprimn
               olap = olap + dconjg(phin(ix,i))*phin(ix,j)
             enddo
             write(20,'(i5,i5,f15.7)') i, j, cdabs(olap)*dr
           enddo
         enddo
         close(20)

         end if 

        return
      end subroutine

    !> Helping functions and subroutines called within the main wavefunction init routines
    !cw> map detl into a number
    integer function checksum(iarr, nrorb_input)
      implicit none
      integer :: iarr(nel), base, ie, cs, nrorb_input

      cs = 0
      ! base gives the number of virtual orbitals plus one
      base = nrorb_input - nel + 1
      do ie=1, nel
        cs = cs + (iarr(nel-ie+1) - nel + ie - 1)*(base**ie)
      enddo

      checksum = cs

    end function checksum



    

!> Get index 1.
    subroutine getindex(iarr,ind,vorz, nrorb_input, hval3_input, nrindep_input)

      implicit none
      integer :: iarr(nel), ind, vorz, perm, h, ie, ie2, cs
      integer :: nrorb_input, nrindep_input
      integer :: hval3_input(nrindep_input)

      perm = 0

      do ie=1, nel-1
        do ie2=2, nel
          if (iarr(ie2-1) > iarr(ie2)) then
            h = iarr(ie2)
            iarr(ie2) = iarr(ie2-1)
            iarr(ie2-1) = h
            perm = perm + 1
          endif
        enddo
      enddo

      cs = checksum(iarr, nrorb_input)
      ind = 1
      do while ((hval3_input(ind) /= cs) .and. (ind <= nrindep_input))
        ind = ind + 1
      enddo
      if (ind > nrindep_input) then
        write(*,*) "Index not found!"
        write(*,*) (iarr(ie),ie=1,nel)
        stop
      endif

      if (mod(perm,2) == 0) then
        vorz = 1
      else
        vorz = -1
      endif

      return
    end subroutine

  ! cw> input: iarr, e.g., 1234 presents a determianent
  ! cw> output: nth determinant in list
  !> Get index 2.
    subroutine getindex2(iarr,ind,vorz, nrorb_input, hval3_input, nrindep_input)
      integer :: iarr(nel), ind, vorz, perm, h, ie, ie2, cs
      integer :: sum_sz_tmp, orb_temp, i

      integer :: nrorb_input, nrindep_input
      integer :: hval3_input(nrindep_input)

      perm = 0
      ! sort so that smallest SO stands left and largest right
      ! keep track of the required permutations
      do ie=1, nel-1
        do ie2=2, nel
          if (iarr(ie2-1) > iarr(ie2)) then
            h = iarr(ie2)
            iarr(ie2) = iarr(ie2-1)
            iarr(ie2-1) = h
            perm = perm + 1
          endif
        enddo
      enddo

      ! two of the same SOs - should not happen
      if (iszero(iarr)) then
        ind = 1
        vorz = 0
        return
      else
        cs = checksum(iarr, nrorb_input)
        ind = 1
        ! find out which number determinant this single-hole function+SO corresponds to
        do while (hval3_input(ind) /= cs)
          ind = ind + 1
          if (ind > nrindep_input) then
                    vorz = 0
                    return
          endif
        enddo
      endif

      if (mod(perm,2) == 0) then
        vorz = 1
      else
        vorz = -1
      endif

     ! two of the same SOs - should not happen
      if (iszero(iarr)) then
        ind = 1
        vorz = 0
        return
      else
        cs = checksum(iarr, nrorb_input)
        ind = 1
       ! find out which number determinant this single-hole function+SO corresponds to
        do while (hval3_input(ind) /= cs)
            ind = ind + 1
            if (ind > nrindep_input) then
              vorz = 0
              return
            end if

        enddo
      endif

      if (mod(perm,2) == 0) then
        vorz = 1
      else
        vorz = -1
      endif

      sum_sz_tmp = 0
      do i = 1, nel
        orb_temp = iarr(  i )

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
        vorz = 0
      end if

      return
    end subroutine


  !cw>: input: iarr, e.g., 1234
  !cw>: output: the number of iarr in detl list
  !> Find the indices of determinants.
    subroutine getindexnece(iarr,ind,vorz,nrorb_input,hval3_input,nrindep_input)

      integer :: iarr(nel), ind, vorz, perm, h, ie, ie2, cs
      integer :: nrorb_input,nrindep_input
      integer :: hval3_input(nrindep_input)

      ! find out if permutation of electrons leads to same determinant
      perm = 0
      do ie=1, nel-1
        do ie2=2, nel
          if (iarr(ie2-1) > iarr(ie2)) then
            h = iarr(ie2)
            iarr(ie2) = iarr(ie2-1)
            iarr(ie2-1) = h
            perm = perm + 1
          endif
        enddo
      enddo

      if (iszero(iarr)) then
        ind = 1
        vorz = 0
        return

      else
        cs = checksum(iarr, nrorb_input)
        ind = 1
        do while (hval3_input(ind) /= cs)
          ind = ind + 1
          if (ind > nrindep_input) then
            ind = 1
            vorz = 0
            return
          endif
        enddo
      endif

      if (mod(perm,2) == 0) then
        vorz = 1
      else
        vorz = -1
      endif

      return
    end subroutine

!> Find out if difference between all other electrons is zero.
    logical function iszero(iarr)

      integer :: iarr(nel)
      integer :: ie

      iszero = .false.

      do ie=1, nel-1
        if (iarr(ie) == iarr(ie+1)) iszero = .true.
      enddo

      return
    end function

!> Subroutine to get the nuclear CAP. Not completed at this point.
    subroutine getnuclcap()

      implicit none
      integer :: i, rc

      rc = nrprimn
      ! find out grid point where CAP starts
      i = 1
      do while(r(i) < wgrid)
        i = i + 1
      enddo
      rc = i

      write(*,*) "Found CAP starting point ", rc, " for position ",wgrid

      wcap(:) = c0
      ! only right of onset
      do i=rc, nrprimn
        wcap(i)  = -ci*wsign*etan*(r(i) - r(rc))**bn
      enddo

      write(*,*) "grid            Wcap"
      do i=rc, nrprimn
        write(*,'("wcap"(i6),(es18.11,1x),(es18.11)"i")') i, wcap(i)
      enddo

    end subroutine

!> Get list of allowed excitations between determinants.
    subroutine getallowed1st(cisshf_input, nrallow_input, allowlist_input, nrorb_input, iarr_input, &
    & dim_shdl_input, shdl_input, nrindep_input, hval3_input)

      implicit none
      integer :: nrorb_input, dim_shdl_input, nrindep_input
      integer :: iarr_input(nel), allowlist_input(nrorb_input-(nel-1)), shdl_input(dim_shdl_input), hval3_input(nrindep_input)
      integer :: nrallow_input, i, is, cisshf_input, ind, vorz
      logical :: found
      integer :: sum_nsz

      nrallow_input = 0

      do is=1, nrorb_input
        found = .false.
        do i=1, nel-1
          if (is == shdl_input((cisshf_input-1)*(nel-1)+i)) found = .true.
        enddo

        ! build a determinant using this SO and the SOs from the
        ! single-hole determinant
        iarr_input(1) = is
        do i=1, nel-1
          iarr_input(i+1) = shdl_input((cisshf_input-1)*(nel-1)+i)
        enddo

        call getindexnece(iarr_input,ind,vorz,nrorb_input,hval3_input,nrindep_input)

        sum_nsz = 0
        do i=1, nel
            if (mod(iarr_input(i),2) == 0) then
                sum_nsz = sum_nsz - 1
            else
                sum_nsz = sum_nsz + 1
            endif
        enddo

        if (flag_spinorbital == 0) then

          if ((.not. found) .and. (mod(is,2) == 0) .and. (vorz /= 0) .and. (sum_nsz == 0) ) then
            nrallow_input = nrallow_input + 1
            allowlist_input(nrallow_input) = is
          endif
        else

          if ((.not. found) .and. (vorz /= 0) .and. (sum_nsz == nsz)) then

            nrallow_input = nrallow_input + 1
            allowlist_input(nrallow_input) = is
          endif
        end if
      enddo

      return
    end subroutine

!> Interpolates venmat = hmat - tmat on a grid for the nuclear coord.
!! using akima interpolation
    subroutine akimpotintvenao()

      use akima_interpol

      implicit none
      real(dp) :: ypt(nrensp)
      real(dp) :: tempvals(nrprimn)

      integer  :: j, ix, mu, nu
      character(40) :: printform
      printform = '(2(i5,1x),(f13.7,1x),(f18.12))'

      ! do akima interpolations
      do mu = 1, nrprime
        do nu = 1, nrprime
          ! avoid array temporaries
          do ix = 1, nrensp
            tempvals(ix) = 0.0_dp
            ypt(ix) = hmat(mu,nu,ix) - tmat(mu,nu)
          enddo
          call akima_interp(nrensp, nrprimn, r, rsp, ypt, tempvals)
          venmat(mu,nu,:) = tempvals(:)
        enddo
      enddo

      write (*,*) 'matrix element of interpolated <lambda mu |Ven|mu \lambda>'
      write (*,*) 'mu lambda R_lambda <lambda mu |Ven|mu \lambda>'
      do mu = 1, nrprime
        do j = 1, nrprimn
          write(*,printform) mu, j, r(j), venmat(mu,mu,j)
        enddo
        write(*,*) ""
      enddo

    end subroutine

!> Interpolates venmat = hmat - tmat on a grid for the nuclear coord.
!! using cubic spline.
    subroutine potintvenao()

      implicit none
      real(dp) ::  denom(nrensp), nom
      integer  :: i, j, ix, mu, nu
      character(50) ::  printform
      printform = '(2(i5,1x),(f13.7,1x),(f18.12))'

      do i=1,nrensp
        denom(i) = 1.0_dp
      enddo

      do i=1,nrensp
        do j=1,nrensp
          if (i /= j) denom(i) = denom(i)*(rsp(i) - rsp(j))
        enddo
      enddo

      do mu=1,nrprime
        do nu=1,nrprime
          do ix=1,nrprimn
            venmat(mu,nu,ix) = 0.0_dp
            do i=1, nrensp
              nom = 1.0_dp
              do j=1,nrensp
                if (i /= j) nom = nom*(r(ix) - rsp(j))
              enddo
              venmat(mu,nu,ix) = venmat(mu,nu,ix) + nom*(hmat(mu,nu,i) - tmat(mu,nu))/denom(i)
            enddo
          enddo
        enddo
      enddo

      write(*,*) "# Printing all V_en (electron-nuclear) matrix elements"
      do i=1,nrensp
        do mu=1,nrprime
          write(*,printform) i, mu, rsp(i), hmat(mu,mu,i) - tmat(mu,mu)
        enddo
      enddo
      do mu=1,nrprime
        do ix=1,nrprimn
          write(*,printform) mu, ix, r(ix), venmat(mu,mu,ix)
        enddo
        write(*,*) ""
      enddo

      return

    end subroutine

    

  end module

  !> @file
  !> @brief contains the initialization of the wave function and
  !! Hamiltonian elements.
