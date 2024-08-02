!> Contains all things related to real and imaginary time propagation.
module propa
    use params
    use moreadin
    use globalvars
    use inputvars
    use utils
    use omp_lib
    use analyse

    complex(dp), allocatable :: irho(:,:)
    complex(dp), allocatable :: irho_alpha(:,:), irho_beta(:,:)
    complex(dp), allocatable :: irho_spinorbital(:,:)
    complex(dp), allocatable :: irho_input(:,:)
    complex(dp), allocatable :: mf1(:,:,:,:)
    complex(dp), allocatable :: mf1_spinorbital(:,:,:,:)
    complex(dp), allocatable :: mf1_input(:,:,:,:)
    real(dp), allocatable :: np(:)
    real(dp), allocatable :: np_alpha(:),np_beta(:)
    real(dp), allocatable :: np_spinorbital(:)
    real(dp), parameter ::                          &
     c2  =    0.526001519587677318785587544488d-01, &
     c3  =    0.789002279381515978178381316732d-01, &
     c4  =    0.118350341907227396726757197510d+00, &
     c5  =    0.281649658092772603273242802490d+00, &
     c6  =    0.333333333333333333333333333333d+00, &
     c7  =    0.25d+00,                             &
     c8  =    0.307692307692307692307692307692d+00, &
     c9  =    0.651282051282051282051282051282d+00, &
     c10 =    0.6d+00,                              &
     c11 =    0.857142857142857142857142857142d+00, &
     b1  =    5.42937341165687622380535766363d-02,  &
     b6  =    4.45031289275240888144113950566d+00,  &
     b7  =    1.89151789931450038304281599044d+00,  &
     b8  =   -5.80120396001058478146721142270d+00,  &
     b9  =    3.11164366957819894408916062370d-01,  &
     b10 =   -1.52160949662516078556178806805d-01,  &
     b11 =    2.01365400804030348374776537501d-01,  &
     b12 =    4.47106157277725905176885569043d-02,  &
     bh1 =    0.244094488188976377952755905512d+00, &
     bh2 =    0.733846688281611857341361741547d+00, &
     bh3 =    0.220588235294117647058823529412d-01
    real(dp), parameter ::                          &
     e1  =    0.1312004499419488073250102996d-01,   &
     e6  =   -0.1225156446376204440720569753d+01,   &
     e7  =   -0.4957589496572501915214079952d+00,   &
     e8  =    0.1664377182454986536961530415d+01,   &
     e9  =   -0.3503288487499736816886487290d+00,   &
     e10 =    0.3341791187130174790297318841d+00,   &
     e11 =    0.8192320648511571246570742613d-01,   &
     e12 =   -0.2235530786388629525884427845d-01
    real(dp), parameter ::                          &
     a21 =    5.26001519587677318785587544488d-02,  &
     a31 =    1.97250569845378994544595329183d-02,  &
     a32 =    5.91751709536136983633785987549d-02,  &
     a41 =    2.95875854768068491816892993775d-02,  &
     a43 =    8.87627564304205475450678981324d-02,  &
     a51 =    2.41365134159266685502369798665d-01,  &
     a53 =   -8.84549479328286085344864962717d-01,  &
     a54 =    9.24834003261792003115737966543d-01,  &
     a61 =    3.70370370370370370370370370370d-02,  &
     a64 =    1.70828608729473871279604482173d-01,  &
     a65 =    1.25467687566822425016691814123d-01,  &
     a71 =    3.7109375d-02,                        &
     a74 =    1.70252211019544039314978060272d-01,  &
     a75 =    6.02165389804559606850219397283d-02,  &
     a76 =   -1.7578125d-02
    real(dp), parameter ::                          &
     a81 =    3.70920001185047927108779319836d-02,  &
     a84 =    1.70383925712239993810214054705d-01,  &
     a85 =    1.07262030446373284651809199168d-01,  &
     a86 =   -1.53194377486244017527936158236d-02,  &
     a87 =    8.27378916381402288758473766002d-03,  &
     a91 =    6.24110958716075717114429577812d-01,  &
     a94 =   -3.36089262944694129406857109825d+00,  &
     a95 =   -8.68219346841726006818189891453d-01,  &
     a96 =    2.75920996994467083049415600797d+01,  &
     a97 =    2.01540675504778934086186788979d+01,  &
     a98 =   -4.34898841810699588477366255144d+01,  &
     a101 =   4.77662536438264365890433908527d-01,  &
     a104 =  -2.48811461997166764192642586468d+00,  &
     a105 =  -5.90290826836842996371446475743d-01,  &
     a106 =   2.12300514481811942347288949897d+01,  &
     a107 =   1.52792336328824235832596922938d+01,  &
     a108 =  -3.32882109689848629194453265587d+01,  &
     a109 =  -2.03312017085086261358222928593d-02
    real(dp), parameter ::                          &
     a111 =  -9.37142430085987325717040216580d-01,  &
     a114 =   5.18637242884406370830023853209d+00,  &
     a115 =   1.09143734899672957818500254654d+00,  &
     a116 =  -8.14978701074692612513997267357d+00,  &
     a117 =  -1.85200656599969598641566180701d+01,  &
     a118 =   2.27394870993505042818970056734d+01,  &
     a119 =   2.49360555267965238987089396762d+00,  &
     a1110 = -3.04676447189821950038236690220d+00,  &
     a121 =   2.27331014751653820792359768449d+00,  &
     a124 =  -1.05344954667372501984066689879d+01,  &
     a125 =  -2.00087205822486249909675718444d+00,  &
     a126 =  -1.79589318631187989172765950534d+01,  &
     a127 =   2.79488845294199600508499808837d+01,  &
     a128 =  -2.85899827713502369474065508674d+00,  &
     a129 =  -8.87285693353062954433549289258d+00,  &
     a1210 =  1.23605671757943030647266201528d+01,  &
     a1211 =  6.43392746015763530355970484046d-01

    contains

!> Main propagation driver
    subroutine propagation()
      !> .t extensions mean -> time-dependent files
      ! F.No: File Name
      !  30 : expec.t
      !  31 : steps.t
      !  33 : density.xt <- never used
      !  35 : npop.t
      !  36 : wavef.chk/finalpsi/finalrho -> wavef_f, fpsi
      !  36 : finalpsi/finalrho
      !  37 : acfphi.t
      !  38 : rhoe.chk
      !  39 : rhon.chk
      !  40 : rhon0.Rt
      !  41 : rhondiff.Rt
      !  42 : trrho2.t
      !  43 : eigvals.t
      !  45 : efield.t -> efield_f
      !  46 : wavee.chk -> wavee_f
      !  47 : waven.chk -> waven_f
      !  350: npop_spinorbital.t
      !-------------------------------------
      ! possible conventions
      ! rhon + .dat -> rhon_d
      ! rhon + file -> rhon_f
      ! rhon + .log -> rhon_l
      ! rhon + .t

      implicit none
      complex(dp) :: phino(nrprimn,nrspf)
      complex(dp) :: psi(dgldim), val(13), acfspf(nrorb)
      complex(dp) :: psi_spinorbital(dgldim_spinorbital), val_spinorbital(13),acfspf_spinorbital(nrorb_spinorbital)
      complex(dp) :: Ao(nrindep*nrspf)
      complex(dp) :: Ao_spinorbital(nrindep_spinorbital*nrspf)
      complex(dp) :: phio(nrprime,nrorb)
      complex(dp) :: phio_spinorbital(nrprime,nrorb_spinorbital)
      complex(dp) :: trrho2
      complex(dp) :: trrho2_spinorbital
      complex(dp) :: irhon(nrspf,nrspf)
      complex(dp) :: rhonRt(nrprimn), rhonRt0(nrprimn)
      complex(dp) :: rhonRtbuff
      real(dp)    :: valreal(13)
      real(dp)    :: valreal_spinorbital(13)
      real(dp)    :: dt, tbeg, tend, facold, t_init
      real(dp)    :: tfin, tout, ef(3)
      real(dp)    :: npn(nrspf)
      real(dp)    :: npn_spinorbital(nrspf)
      integer     :: ix, in, i, it, nt
      integer     :: is, js, j, mu
      integer     :: newgrd, oldgrd, indepA, indepA_spinorbital
      integer     :: waven_f, wavee_f, wavee_f_spinorbital
      integer     :: fpsi, frho, fexpec, fentro
      integer     ::  fpsi_spinorbital,  fexpec_spinorbital, fentro_spinorbital
      character     :: ncol
      character(10)  :: npf, nnn,   ccol
      character(50) :: formatv1, formatv2, Formtt
      character(35) :: form_wavee, form_wavec
      character(35) :: form_waven, form_waven2
      integer :: efieldf
      integer :: efieldf_spinorbital
      integer :: trrho2f
      integer :: trrho2f_spinorbital, flag_print, spsi_0, fpsi_spinorbital_0
      real(dp)     :: rep(nrorb), imp(nrorb), rea, ima, wtime_1
      real(dp)     :: rep_spinorbital(nrorb_spinorbital),imp_spinorbital(nrorb_spinorbital), rea_spinorbital, ima_spinorbital

      formatv1 = "(1(f8.3,1x),4(f12.7,1x),10(f16.10,1x))"
      formatv2 = "(2(f16.10,1x))"
      formtt = "(1(a8,1x),(a9,1x),3(a12,1x),10(a16,1x))"

      npf = stri(2*nrorb + 1)
      ncol = stri(2*nrspf)
      write(nnn,'(i2)') 2*nrorb
      write(ccol,'(i2)') nrorb
      g_counter = 0
      hel_counter = 0
      hel2_counter = 0
      hel2_transformation_counter = 0
      hel2a_counter = 0
      hel2a2_counter = 0
      hel3_counter = 0
      hamelem_counter = 0

      runrk8_counter = 0
      runrk8_spinorbital_counter = 0
      runrk8_general_counter = 0

      t_counter = 0
      t_spinorbital_counter = 0
      t_general_counter = 0
      heln_spinorbital_counter = 0
      heln_counter = 0
      flag_print = 3

      tend = 0.0_dp

      form_wavee = '('//nnn//'(e23.16,1x))'
      form_wavec = '('//nnn//'(e23.16,1x))'
      form_waven = '((i5,1x),(F15.9,1x),'//ncol//'(e23.16,1x))'
      form_waven2 = '('//ncol//'(e23.16,1x))'

      if (flag_spinorbital == 0) then
        open(newunit=efieldf,file="efield.t")
        open(newunit=fexpec,file="expec.t")
        open(newunit=fentro,file="entropy.t")
        write(fexpec,'(A6,I3,A6,I3,A6,I3,A6,I3,A6,I3,1x,A120)')  "nel", nel, "nsz", nsz, "nrorb" , &
        & nrorb,"nrfrz" ,nrorb_fc, "nrspf", nrspf, intdir
        write(fexpec,formtt) "time","norm","x","y","z","R","Te","Tn","Vee","Vnn","Ven","Htot","Re(Acf)","Im(Acf)"

        if (restart == 5 .and. aucofu == 2) then
          write (*,*) 'read spsi_0, restart acf may require psi0 from the very begining-1'
          open(newunit=spsi_0,file="startpsi_0")
          write (*,*) 'read spsi_0, restart acf may require psi0 from the very begining-2'

          do i=1, nrindep*nrspf
            read(spsi_0,*) rea, ima
            Ao(i) = dcmplx(rea,ima)
          end do

          do mu=1, nrprime
            read(spsi_0,'('//stri(2*nrorb)//'(e23.16,1x))') (rep(in),in=1,nrorb), (imp(in),in=1,nrorb)
            do in=1, nrorb
              phio(mu,in) = dcmplx(rep(in),imp(in))
            enddo
          enddo

          write(ncol,'(i1)') 2*nrspf

          do ix=1, nrprimn
            read(spsi_0,'('//stri(2*nrspf)//'(e23.16,1x))') (rep(in),in=1,nrspf), (imp(in),in=1,nrspf)
            do in=1, nrspf
              phino(ix,in) = dcmplx(rep(in),imp(in))
            enddo
          enddo

          close(spsi_0)
        else
          phio(:,:) = phi(:,:)
          phino(:,:) = phin(:,:)
          Ao(:) = A(:)
        end if

      else if (flag_spinorbital == 1) then
        open(newunit=efieldf_spinorbital,file="efield_spinorbital.t")
        open(newunit=fexpec_spinorbital,file="expec_spinorbital.t")
        open(newunit=fentro_spinorbital,file="entropy_spinorbital.t")
        write(fexpec_spinorbital,'(A6,I3,A6,I3,A6,I3,A6,I3,A6,I3,1x,A50)')  "nel", nel, "nsz", nsz, "nrorb" , &
&nrorb,"nrfrz" ,nrorb_fc, "nrspf", nrspf, intdir
        write(fexpec_spinorbital,formtt) "time","norm","x","y","z","R","Te","Tn","Vee","Vnn","Ven","Htot","Re(Acf)","Im(Acf)"
        ! restart acf calculation and read initial psi for acf
        if ( restart == 5 .and. aucofu == 2) then
          open(newunit=spsi_0,file="startpsi_spinorbital_0")
          do i=1, nrindep_spinorbital*nrspf
            read(spsi_0,*) rea_spinorbital, ima_spinorbital
            Ao_spinorbital(i) = dcmplx(rea_spinorbital,ima_spinorbital)
          end do
          write (*,*) 'read spsi_0, restart acf may require psi0 from the very begining'
          do mu=1, nrprime
            read(spsi_0,'('//stri(nrorb_spinorbital)//'(e23.16,1x))') (rep_spinorbital(in),in=1,nrorb_spinorbital), &
            & (imp_spinorbital(in),in=1,nrorb_spinorbital)
            do in=1, nrorb_spinorbital
              phio_spinorbital(mu,in) = dcmplx(rep_spinorbital(in),imp_spinorbital(in))
            enddo
          enddo
          write(ncol,'(i1)') 2*nrspf
          do ix=1, nrprimn
            read(spsi_0,'('//stri(2*nrspf)//'(e23.16,1x))') (rep_spinorbital(in),in=1,nrspf), (imp_spinorbital(in),in=1,nrspf)
            do in=1, nrspf
              phino(ix,in) = dcmplx(rep_spinorbital(in),imp_spinorbital(in))
            enddo
          enddo
          close(spsi_0)
        else
          phio_spinorbital(:,:) = phi_spinorbital(:,:)
          phino(:,:) = phin(:,:)
          Ao_spinorbital(:) = A_spinorbital(:)
        end if

      else
         write (*,*) 'flag_spinorbital 0 /1 not supprted'
         stop
      end if

! IU: this should be done in the allocation subroutine
      if (.not. allocated(irho)) allocate(irho(2*nrorb, 2*nrorb))
      if (.not. allocated(irho_alpha)) allocate(irho_alpha(nrorb_alpha, nrorb_alpha))
      if (.not. allocated(irho_beta)) allocate(irho_beta(nrorb_beta, nrorb_beta))
      if (.not. allocated(irho_spinorbital)) allocate(irho_spinorbital(nrorb_spinorbital, nrorb_spinorbital))
      if (.not. allocated(np)) allocate(np(2*nrorb))
      if (.not. allocated(np_spinorbital)) allocate(np_spinorbital(nrorb_spinorbital))

      if (flag_spinorbital == 0) then
        call hamelem(phi, nrorb, hel, hel2, hel3, venmatmo)
        call calc_rhoe(rho, nrorb, nrorb*2, nrindep, hvalrho, rhoss, A  )
        wtime_1 = omp_get_wtime()
        call calc_rhoe2()
      end if

      if (flag_spinorbital == 1) then
        call hamelem(phi_spinorbital, nrorb_spinorbital, hel_spinorbital, hel2_spinorbital, hel3_spinorbital, &
        & venmatmo_spinorbital)
        call calc_rhoe(rho_spinorbital, nrorb_spinorbital, nrorb_spinorbital, nrindep_spinorbital, &
        & hvalrho_spinorbital, rhoss_spinorbital, A_spinorbital  )
        wtime_1 = omp_get_wtime()
        call calc_rhoe2()
      end if

      if (flag_spinorbital == 0) then
        call invert_rhoe(irho, np, rho,  nrorb*2 )
      end if

      if (flag_spinorbital == 1) then
        call invert_rhoe(irho_spinorbital, np_spinorbital, rho_spinorbital,  nrorb_spinorbital )
      end if

      if (flag_spinorbital == 0) then
         call calc_rhon(A,nrindep)
      end if

      if (flag_spinorbital == 1) then
        call calc_rhon(A_spinorbital,nrindep_spinorbital)
      end if

      call invert_rhon(irhon,npn)

      if (flag_spinorbital == 0) then
        call analysis(Ao,phio,phino, phi, val,valreal, hel, rho, A, nrindep, nrorb, nrorb*2, &
        & detl, hel2, venmatmo)
      end if

      if (flag_spinorbital == 1) then
         call analysis(Ao_spinorbital,phio_spinorbital,phino, phi_spinorbital, val_spinorbital, valreal_spinorbital, &
        & hel_spinorbital, rho_spinorbital, A_spinorbital, nrindep_spinorbital, nrorb_spinorbital, nrorb_spinorbital, &
        & detl_spinorbital, hel2_spinorbital, venmatmo_spinorbital)

      end if

      ! write out all expectation values.
      if (flag_spinorbital == 0) then
       write(fexpec,formatv1) 0.0_dp, (valreal(i), i=1,13)
        call flush(fexpec)
        write(fentro,'(A6,I3,A6,I3,A6,I3,A6,I3,1x,A50)')  "nel", nel, "nsz", nsz, "nrorb" ,nrorb, "nrspf", nrspf, intdir
        call rdm1analysis(0.0_dp, fentro,  nrorb*2, rho)
      end if

      if (flag_spinorbital == 1) then
        write(fexpec_spinorbital,formatv1) t_initial, (valreal_spinorbital(i), i=1,13)
        call flush(fexpec_spinorbital)
        write(fentro_spinorbital,'(A6,I3,A6,I3,A6,I3,A6,I3,1x,A50)')  "nel", nel, "nsz", nsz, "nrorb" ,nrorb, "nrspf", nrspf, intdir
        call rdm1analysis(0.0_dp, fentro_spinorbital,  nrorb_spinorbital, rho_spinorbital)


      end if

      if (logwavef) then
        if (flag_spinorbital == 0) then 
          open(newunit=indepA, file='indepA_elem.chk')
          open(newunit=wavee_f,file="wavee.chk")
        end if

        if (flag_spinorbital == 1) then 
          open(newunit=indepA_spinorbital, file='indepA_elem_spinorbital.chk')
          open(newunit=wavee_f_spinorbital,file="wavee_spinorbital.chk")
        end if


        open(newunit=waven_f,file="waven.chk")

        if (flag_spinorbital == 0) then
          do i=1, nrindep*nrspf
            write(indepA,'(2(e23.16,1x))') dreal(A(i)), dimag(A(i))
          enddo
        end if

        if (flag_spinorbital == 1) then
          do i=1, nrindep_spinorbital*nrspf
            write(indepA_spinorbital,'(2(e23.16,1x))') dreal(A_spinorbital(i)), dimag(A_spinorbital(i))
          enddo
        end if

        ! need so version
        if (flag_spinorbital == 0) then
          do mu=1, nrprime
            write(wavee_f,form_wavee) (dreal(phi(mu,in)),in=1,nrorb), (dimag(phi(mu,in)),in=1,nrorb)
          enddo
          do ix=1, nrprimn
            write(waven_f,form_waven) 0, r(ix), (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
          enddo

          call flush(wavee_f)
          call flush(waven_f)
        end if


        if (flag_spinorbital == 1) then
          do mu=1, nrprime
            write(wavee_f_spinorbital,form_wavee) (dreal(phi_spinorbital(mu,in)),in=1,nrorb_spinorbital),&
            & (dimag(phi_spinorbital(mu,in)),in=1,nrorb_spinorbital)
          enddo
          do ix=1, nrprimn
            write(waven_f,form_waven) 0, r(ix), (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
          enddo

          call flush(wavee_f_spinorbital)
          call flush(waven_f)
        end if

        if (flag_spinorbital == 0) then
          open(38,file="rhoe.chk")
          do i=1, 2*nrorb
            do j=1, 2*nrorb
              write(38,formatv2) dreal(rho(i,j)), dimag(rho(i,j))
            enddo
          enddo

          call flush(38)

          open(39,file="rhon.chk")
          do i=1, nrspf
            do j=1, nrspf
              write(39,formatv2) dreal(rhon(i,j)), dimag(rhon(i,j))
            enddo
          enddo
          call flush(39)

        else
          open(38,file="rhoe_spinorbital.chk")
          do i=1,nrorb_spinorbital
            do j=1, nrorb_spinorbital
              write(38,formatv2) dreal(rho_spinorbital(i,j)), dimag(rho_spinorbital(i,j))
            enddo
          enddo

          call flush(38)

          open(39,file="rhon_spinorbital.chk")
          do i=1, nrspf
            do j=1, nrspf
              write(39,formatv2) dreal(rhon(i,j)), dimag(rhon(i,j))
            enddo
          enddo
          call flush(39)

        end if

        open(40,file="rhon.Rt")
        do ix=1, nrprimn
          rhonRtbuff = c0
          do i=1, nrspf
            do j=1, nrspf
              rhonRtbuff = rhonRtbuff + phin(ix,i)*rhon(i,j)*dconjg(phin(ix,j))
            enddo
          enddo
          rhonRt0(ix) = rhonRtbuff
        enddo

        do ix=1, nrprimn
          write(40,'(3(f16.10,1x))') r(ix), dreal(rhonRt0(ix)), dimag(rhonRt0(ix))
        enddo

        rhonRtbuff = c0
        do ix=1, nrprimn
          rhonRtbuff = rhonRtbuff + rhonRt0(ix)
        enddo
        call flush(40)

      endif

      if (flag_spinorbital == 0) then
        open(37, file="acfphi.t")
        call acfphi(phio, acfspf, nrorb, phi)
        write(37,'('//trim(npf)//'(f18.12))') 0.0, (dreal(acfspf(in)), dimag(acfspf(in)),in=1,nrorb)
        call flush(37)
      end if

      if (flag_spinorbital == 1) then
        open(370, file="acfphi_spinorbital.t")
        call acfphi(phio_spinorbital,acfspf_spinorbital, nrorb_spinorbital, phi_spinorbital)
        write(370,'('//trim(npf)//'(f18.12))') 0.0, (dreal(acfspf_spinorbital(in)), &
& dimag(acfspf_spinorbital(in)),in=1,nrorb_spinorbital)
        call flush(370)
      end if

      if (logev) open(43,file="eigvals.t")

      open(31, file="steps.t")

      if (flag_spinorbital == 0) then
        open(35, file="npop.t")
        open(newunit=trrho2f, file="trrho2.t")
      end if

      if (flag_spinorbital == 1) then
        open(350, file="npop_spinorbital.t")
        open(newunit=trrho2f_spinorbital, file="trrho2_spinorbital.t")
      end if





      !> ok yeah, so this does the \Tr \rho**2
      if (flag_spinorbital == 0) then
      !>  take trace of rho squared -> tr[rho]^2
        trrho2 = c0
        do is=1, 2*nrorb
          do js=1, 2*nrorb
            trrho2 = trrho2 + rho(is,js)*rho(js,is)
          enddo
        enddo
  
        write(trrho2f,'(1(f8.3,1x),2(F15.10,2x))') 0.0, dreal(trrho2), dimag(trrho2)
        call flush(trrho2f)

        write(35,'('//trim(npf)//'(F15.10,1x))') tend*au2fs, (nel*np(i),i=1,2*nrorb)
        call flush(35)
      end if

      if (flag_spinorbital == 1) then
        trrho2_spinorbital = c0
        do is=1, nrorb_spinorbital
          do js=1, nrorb_spinorbital
            trrho2_spinorbital = trrho2_spinorbital + rho_spinorbital(is,js)*rho_spinorbital(js,is)
          enddo
        enddo

        write(trrho2f_spinorbital,'(1(f8.3,1x),2(F15.10,2x))') 0.0, dreal(trrho2_spinorbital), dimag(trrho2_spinorbital)
        call flush(trrho2f_spinorbital)

        write(350,'('//trim(npf)//'(F15.10,1x))') tend*au2fs, (nel*np_spinorbital(i),i=1,nrorb_spinorbital)
        call flush(350)        
      end if


      tfin = tfin0/au2fs
      tout = tout0/au2fs
      dt = dt0/au2fs
      nt = dint(tfin/tout)
      facold = 1.0d-4
      t_init = t_initial/au2fs
      oldgrd = inigrd
      newgrd = oldgrd

      !************************!
      !>> PROPAGATION BEGINS <<!
      !************************!

      !> beginning of time loop
      do it=1, nt
        !tout is the output time step
        !nt is the number of time step
        tbeg = (it - 1)*tout
        tend = it*tout
        if (restart == 2 .or. restart == 5 ) then
          tbeg = tbeg + t_init
          tend = tend + t_init
        endif

        if (flag_spinorbital == 0) then
          call mapformat(psi, 0_i64, dgldim, nrindep, nrorb, A, phi)
        end if

        if (flag_spinorbital == 1) then
          call mapformat(psi_spinorbital, 0_i64, dgldim_spinorbital, nrindep_spinorbital, &
          & nrorb_spinorbital, A_spinorbital, phi_spinorbital)
        end if

        if ( flag_spinorbital == 0 ) then
           call runrk8(psi, tbeg, tend, facold, dt, np, npn, &
    &  nrindep, dgldim,&
    & nrorb, nrorb*2, max_nrindep, hvalrho, rhoss, &
    & detl,  allow1, phi, phi2, A, rho,  hel, hel2, hel3, &
    & venmatmo)
        end if

        if ( flag_spinorbital == 1 ) then
           call runrk8(psi_spinorbital, tbeg, tend, facold, dt, np_spinorbital, npn_spinorbital, &
    &  nrindep_spinorbital, dgldim_spinorbital,&
    & nrorb_spinorbital, nrorb_spinorbital, max_nrindep_spinorbital, hvalrho_spinorbital, rhoss_spinorbital, &
    & detl_spinorbital, allow1_spinorbital, phi_spinorbital, phi2_spinorbital, A_spinorbital, &
    & rho_spinorbital,  hel_spinorbital, hel2_spinorbital, hel3_spinorbital, venmatmo_spinorbital)
        end if

        if (flag_spinorbital == 0) then
          call mapformat(psi, 1_i64, dgldim, nrindep, nrorb, A, phi)
        end if

        if (flag_spinorbital == 1) then
          call mapformat(psi_spinorbital, 1_i64, dgldim_spinorbital, nrindep_spinorbital, &
          & nrorb_spinorbital, A_spinorbital, phi_spinorbital)

        end if

  !     this line will set Ao up to the final time after lpw0
        if ((aucofu == 2) .and. (tend < sum(lpw0(:))/au2fs)) then
          phio(:,:) = phi(:,:)
          phio_spinorbital(:,:) = phi_spinorbital(:,:)
          phino(:,:) = phin(:,:)
          Ao(:) = A(:)
          Ao_spinorbital(:) = A_spinorbital(:)
        endif

         if ( flag_spinorbital == 0 ) then
           call hamelem(phi, nrorb, hel, hel2, hel3, venmatmo)
           call analysis(Ao,phio,phino, phi, val, valreal, hel, rho, A, nrindep, nrorb, nrorb*2, &
        & detl, hel2, venmatmo)
        end if

        if ( flag_spinorbital == 1 ) then
          call hamelem(phi_spinorbital, nrorb_spinorbital, hel_spinorbital, hel2_spinorbital, hel3_spinorbital, &
        & venmatmo_spinorbital)
         call analysis(Ao_spinorbital,phio_spinorbital,phino, phi_spinorbital, val_spinorbital, valreal_spinorbital, &
        & hel_spinorbital, rho_spinorbital, A_spinorbital, nrindep_spinorbital, nrorb_spinorbital, nrorb_spinorbital, &
        & detl_spinorbital, hel2_spinorbital, venmatmo_spinorbital)

        end if

        call getefield(tend,ef)

        if ( flag_spinorbital == 0 ) then
          write(efieldf,'((f9.2,1x),3(es23.15,1x))') tend*au2fs, ef(1:3)
          call flush(efieldf)
        end if

        if ( flag_spinorbital == 1 ) then
          write(efieldf_spinorbital,'((f9.2,1x),3(es23.15,1x))') tend*au2fs, ef(1:3)
          call flush(efieldf_spinorbital)
        end if

        if ( flag_spinorbital == 0 ) then
          write(fexpec,formatv1) tend*au2fs, (valreal(i),i=1,13)
          call flush(fexpec)
          call rdm1analysis(tend*au2fs, fentro,  nrorb*2, rho)
          flush(fentro)
        end if

        if ( flag_spinorbital == 1 ) then
          write(fexpec_spinorbital,formatv1) tend*au2fs, (valreal_spinorbital(i),i=1,13)
          call flush(fexpec_spinorbital)
          call rdm1analysis(tend*au2fs, fentro_spinorbital,  nrorb_spinorbital, rho_spinorbital)
          flush(fentro_spinorbital)
        end if

        call flush(31)

        if ( flag_spinorbital == 0 ) then
          write(35,'('//trim(npf)//'(F15.10,1x))') tend*au2fs, (nel*np(i),i=1,2*nrorb)
          call flush(35)
        end if

        if ( flag_spinorbital == 1 ) then
          write(350,'('//trim(npf)//'(F15.10,1x))') tend*au2fs, (nel*np_spinorbital(i),i=1,nrorb_spinorbital)
          call flush(350)
        end if

        !> write out info for each step in finalpsi format so it checkpoints every round
        !> this section will be modified to only checkpoint after n number of time steps
        !> as well as the option to forgo this completely
        if (flag_spinorbital == 0) then
          open(newunit=fpsi,file="finalpsi")
          do i=1, nrindep*nrspf
            write(fpsi,'(2(e23.16,1x))') dreal(A(i)), dimag(A(i))
          enddo
          do ix=1, nrprime
            write(fpsi,form_wavee) (dreal(phi(ix,in)),in=1,nrorb), (dimag(phi(ix,in)),in=1,nrorb)
          enddo
          do ix=1, nrprimn
            write(fpsi,form_waven2) (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
          enddo
          call flush(fpsi)
          close(fpsi)
! IU: this throws an error as fpsi_0 seems not to have been openend
!          do ix=1, nrprime
!            write(fpsi_0,form_wavee) (dreal(phio(ix,in)),in=1,nrorb), (dimag(phio(ix,in)),in=1,nrorb)
!          enddo
!          do ix=1, nrprimn
!            write(fpsi_0,form_waven2) (dreal(phino(ix,in)),in=1,nrspf), (dimag(phino(ix,in)),in=1,nrspf)
!          enddo
!          call flush(fpsi_0)
!          close(fpsi_0)

        else if (flag_spinorbital == 1) then
          open(newunit=fpsi_spinorbital,file="finalpsi_spinorbital")
          do i=1, nrindep_spinorbital*nrspf
            write(fpsi_spinorbital,'(2(e23.16,1x))') dreal(A_spinorbital(i)), dimag(A_spinorbital(i))
          enddo

          do ix=1, nrprime
            write(fpsi_spinorbital,form_wavee) (dreal(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital), &
  &  (dimag(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital)
          enddo

          do ix=1, nrprimn
            write(fpsi_spinorbital,form_waven2) (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
          enddo

          call flush(fpsi_spinorbital)
          close(fpsi_spinorbital)

          open(newunit=fpsi_spinorbital_0,file="finalpsi_spinorbital_0")

          do i=1, nrindep_spinorbital*nrspf
            write(fpsi_spinorbital_0,'(2(e23.16,1x))') dreal(Ao_spinorbital(i)), dimag(Ao_spinorbital(i))
          enddo

          do ix=1, nrprime
            write(fpsi_spinorbital_0,form_wavee) (dreal(phio_spinorbital(ix,in)),in=1,nrorb_spinorbital), &
            & (dimag(phio_spinorbital(ix,in)),in=1,nrorb_spinorbital)
          enddo

          do ix=1, nrprimn
            write(fpsi_spinorbital_0,form_waven2) (dreal(phino(ix,in)),in=1,nrspf), &
            & (dimag(phino(ix,in)),in=1,nrspf)
          enddo

          call flush(fpsi_spinorbital_0)
          close(fpsi_spinorbital_0)

        else
          write (*,*) 'flag_spinorbital not 0/1, supported'
        end if

        if (logwavef) then

          if (flag_spinorbital == 0) then
            do i=1, nrindep*nrspf
              write(indepA,'(2(e23.16,1x))') dreal(A(i)), dimag(A(i))
            enddo
          end if

          if (flag_spinorbital == 1) then
            do i=1, nrindep_spinorbital*nrspf
              write(indepA_spinorbital,'(2(e23.16,1x))') dreal(A_spinorbital(i)), dimag(A_spinorbital(i))
            enddo
          end if

          if (flag_spinorbital == 0) then
            do ix=1, nrprime
              write(wavee_f,form_wavee) (dreal(phi(ix,in)),in=1,nrorb), (dimag(phi(ix,in)),in=1,nrorb)
            enddo
            call flush(wavee_f)
          end if

          if (flag_spinorbital == 1) then
            do ix=1, nrprime
              write(wavee_f_spinorbital,form_wavee) (dreal(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital), &
              & (dimag(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital)
            enddo
            call flush(wavee_f_spinorbital)
          end if

          do ix=1, nrprimn
            write(waven_f,form_waven) it, r(ix), (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
          enddo

          call flush(waven_f)

          if (flag_spinorbital == 0) then
            do i=1, 2*nrorb
              do j=1, 2*nrorb
                write(38,formatv2) dreal(rho(i,j)), dimag(rho(i,j))
              end do
            enddo
            call flush(38)

          else
            do i=1, nrorb_spinorbital
              do j=1, nrorb_spinorbital
                write(38,formatv2) dreal(rho_spinorbital(i,j)), dimag(rho_spinorbital(i,j))
              end do
            enddo
            call flush(38)

          end if

          if (flag_spinorbital == 0) then
            do i=1, nrspf
              do j=1, nrspf
                write(39,formatv2) dreal(rhon(i,j)), dimag(rhon(i,j))
              enddo
            enddo
            call flush(39)

          else
            do i=1, nrspf
              do j=1, nrspf
                write(39,formatv2) dreal(rhon(i,j)), dimag(rhon(i,j))
              enddo
            enddo
            call flush(39)

          end if

          open(41,file="rhondiff.Rt")
          do ix=1, nrprimn
            rhonRtbuff = c0
            do i=1, nrspf
              do j=1, nrspf
                rhonRtbuff = rhonRtbuff + phin(ix,i)*rhon(i,j)*dconjg(phin(ix,j))
              enddo
            enddo
            rhonRt(ix) = rhonRtbuff
            rhonRt(ix) = rhonRt(ix) - rhonRt0(ix)
          enddo

          do ix=1,nrprimn
            write(41,'(2(f10.5,1x),2(F15.10,1x))') tend*au2fs, r(ix), &
                                      dreal(rhonRt(ix)), dimag(rhonRt(ix))
          enddo

          write(41,*) "   "

          rhonRtbuff = c0
          do ix=1, nrprimn
            rhonRtbuff = rhonRtbuff + rhonRt(ix)
          enddo
          call flush(41)
        endif



        if (flag_spinorbital == 0) then
          trrho2 = c0
          do is=1, 2*nrorb
            do js=1, 2*nrorb
              trrho2 = trrho2 + rho(is,js)*rho(js,is)
            enddo
          enddo

          write(trrho2f,'(1(f8.3,1x),2(F15.10,2x))') tend*au2fs, dreal(trrho2), dimag(trrho2)
          call flush(trrho2f)
          
          call acfphi(phio, acfspf, nrorb, phi)
          write(37,'('//trim(npf)//'(f18.12))') tend*au2fs, (dreal(acfspf(in)), dimag(acfspf(in)),in=1,nrorb)
          call flush(37)
        end if

        if (flag_spinorbital == 1) then

          call acfphi(phio_spinorbital,acfspf_spinorbital, nrorb_spinorbital, phi_spinorbital)
          write(370,'('//trim(npf)//'(f18.12))') tend*au2fs, (dreal(acfspf_spinorbital(in)), &
& dimag(acfspf_spinorbital(in)),in=1,nrorb_spinorbital)
          call flush(370)
        end if

      enddo
      !!!> END OF PROPAGATION LOOP

      if (logwavef) then
        if (flag_spinorbital == 0) then
          close(indepA)
        end if

        if (flag_spinorbital == 1) then
          close(indepA_spinorbital)
        end if

! IU: These need to be fixed to correct units
        close(38)
        close(39)
        close(40)
        close(41)

        if (flag_spinorbital == 0) then
          close(wavee_f)
        end if

        if (flag_spinorbital == 1) then
          close(wavee_f_spinorbital)
        end if
 

        close(waven_f)
      endif
      if (logev) close(43)


      if (flag_spinorbital == 0) then
          close(fexpec)
          close(trrho2f)
      end if

      if (flag_spinorbital == 1) then
          close(fexpec_spinorbital)
          close(trrho2f_spinorbital)
      end if    

      close(31)

      if (flag_spinorbital == 0) then
        close(35)
        open(newunit=fpsi,file="finalpsi")
        do i=1, nrindep*nrspf
          write(fpsi,'(2(e23.16,1x))') dreal(A(i)), dimag(A(i))
        enddo
        do ix=1, nrprime
          write(fpsi,form_wavee) (dreal(phi(ix,in)),in=1,nrorb), (dimag(phi(ix,in)),in=1,nrorb)
        enddo
        do ix=1, nrprimn
          write(fpsi,form_waven2) (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
        enddo
        call flush(fpsi)
        close(fpsi)


        ! write out final density
        open(newunit=frho,file="finalrho")
        do i=1, 2*nrorb
          do j=1, 2*nrorb
            write(frho,*) i, j, dreal(rho(i,j)), dimag(rho(i,j))
          enddo
        enddo
        call flush(frho)
        close(frho)

        if (allocated(irho)) deallocate(irho)
        if (allocated(np)) deallocate(np)
      end if


      if (flag_spinorbital == 1) then

        close(350)
        open(newunit=fpsi_spinorbital,file="finalpsi_spinorbital")
        do i=1, nrindep_spinorbital*nrspf
          write(fpsi_spinorbital,'(2(e23.16,1x))') dreal(A_spinorbital(i)), dimag(A_spinorbital(i))
        enddo
        do ix=1, nrprime
          write(fpsi_spinorbital,form_wavee) (dreal(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital), &
          & (dimag(phi_spinorbital(ix,in)),in=1,nrorb_spinorbital)
        enddo
        do ix=1, nrprimn
          write(fpsi_spinorbital,form_waven2) (dreal(phin(ix,in)),in=1,nrspf), (dimag(phin(ix,in)),in=1,nrspf)
        enddo
        call flush(fpsi_spinorbital)
        close(fpsi_spinorbital)


        if (allocated(irho_spinorbital)) deallocate(irho_spinorbital)
        if (allocated(np_spinorbital)) deallocate(np_spinorbital)


      end if

      return
    end subroutine

!> Runge-Kutta solver of 8th order with adaptive step size. Taken from Fortran ODE library.
    subroutine runrk8(psi_input, tbeg, tend, facold, h, np_input, npn_input, &
    &  nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input, rhoss_input, &
    & detl_input, allow1_input, phi_input, phi2_input, A_input, rho_input,  hel_input, hel2_input, hel3_input, &
    & venmatmo_input)
      !> Runge-Kutta 8th order numerical integration of ODE method
      !> Challenge: figure out where all these coefficients came from
      !> Solived it!!! It comes straight from the MCDTDH code! Almost verbatim
      !> mctdh/source/lib/ode/odesim/oded85.f
      !> which itself came from some other people, and stated that code was copyrighted, so...
      !> thankfully it was reformulated without any goto statements
      !> in the MCTDH code they give reference to Hairer, Norsett, and Wanner
      !> Solving Ordinary Differential Equations I, 2 ed. Springer-Verlag, 1993
      !> the authors were apparenetly obsessed with cats eg. pgs 100, 104, 122

      implicit none

      integer :: dgldim_input
      complex(dp) :: psi_input(dgldim_input)
      complex(dp) :: k1(dgldim_input), k2(dgldim_input), k3(dgldim_input), k4(dgldim_input)
      complex(dp) :: k5(dgldim_input), k6(dgldim_input), k7(dgldim_input), k8(dgldim_input)
      complex(dp) :: k9(dgldim_input), k10(dgldim_input), kerr(dgldim_input), kerr2(dgldim_input)
      complex(dp) :: ztmp, y(dgldim_input), y1(dgldim_input)
      real(dp), allocatable, intent(inout) :: np_input(:)
      real(dp) :: npn_input(nrspf)
      real(dp) :: tend, tbeg, h, facold, expo1
      real(dp) ::  fac, fac1, err, err2, deno, wt(dgldim_input), hnew, t
      logical :: lreject
      integer     :: nrindep_input, nrorb_input, nrorb_spinorbital_input, &
      & max_nrindep_input, allow1_input(max_nrindep_input*(nel-1),0:6*nrorb_input), detl_input(nel*nrindep_input), &
      & hvalrho_input(20*nrorb_input*nrorb_input*nrindep_input), rhoss_input
      complex(dp) :: A_input(nrindep_input*nrspf), phi_input(nrprime, nrorb_input), hel_input(nrorb_input,nrorb_input,5+nrensp), &
      & hel2_input(nrorb_input,nrorb_input,nrorb_input,nrorb_input), hel3_input(nrprime, nrorb_input, nrprime, nrorb_input), &
      & rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input), phi2_input(nrprime, nrorb_input), &
      & venmatmo_input(nrorb_input, nrorb_input, nrprimn)

      t = tbeg
      expo1 = 0.125_dp - beta*0.75_dp
      lreject = .false.
      y = psi_input
      k5 = y

      runrk8_general_counter = runrk8_general_counter + 1

      call derivs(y,t,np_input,npn_input, k4, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      do while (t < tend)
        t_general_counter = t_general_counter + 1

        if (lreject .eqv. .false.) then
          y = k5
          k1 = k4
        endif

        ! adjust h, if integrator is close to the end
        if ((t + 1.010_dp*h) > tend) then
          h = tend - t
        endif

        y1 = y + dcmplx(a21*h)*k1

        call derivs(y1,t+c2*h,np_input,npn_input,k2, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

        y1 = dcmplx(a31*h)*k1
        ztmp = dcmplx(a32*h)
        y1 = y1 + dcmplx(a32*h)*k2
        y1 = y1 + cr*y

        call derivs(y1,t+c3*h,np_input,npn_input,k3, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

        !>  step 3
        ztmp = dcmplx(a41*h)
        y1 = ztmp*k1
        ztmp = dcmplx(a43*h)
        y1 = y1 + ztmp*k3
        y1 = y1 + cr*y

        call derivs(y1,t+c4*h,np_input,npn_input,k4, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !>  step 4
      ztmp = dcmplx(a51*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a53*h)
      y1 = y1 + ztmp*k3
      ztmp = dcmplx(a54*h)
      y1 = y1 + ztmp*k4
      y1 = y1 + cr*y

      call derivs(y1,t+c5*h,np_input,npn_input,k5, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !>  step 5
      ztmp = dcmplx(a61*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a64*h)
      y1 =  y1 + k4*ztmp
      ztmp = dcmplx(a65*h)
      y1 = y1 + ztmp*k5
      y1 = y1 + cr*y

      call derivs(y1,t+c6*h,np_input,npn_input,k6, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      ztmp = dcmplx(a71*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a74*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a75*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a76*h)
      y1 = y1 + ztmp*k6
      y1 = y1 + cr*y

      call derivs(y1,t+c7*h,np_input,npn_input,k7, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !>  step 7
      ztmp = dcmplx(a81*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a84*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a85*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a86*h)
      y1 = y1 + ztmp*k6
      ztmp = dcmplx(a87*h)
      y1 = y1 + ztmp*k7
      y1 = y1 + cr*y

      call derivs(y1,t+c8*h,np_input,npn_input,k8, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !>  step 8
      ztmp = dcmplx(a91*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a94*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a95*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a96*h)
      y1 = y1 + ztmp*k6
      ztmp = dcmplx(a97*h)
      y1 = y1 + ztmp*k7
      ztmp = dcmplx(a98*h)
      y1 = y1 + ztmp*k8
      y1 = y1 + cr*y

      call derivs(y1,t+c9*h,np_input,npn_input,k9, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !> step 9
      ztmp = dcmplx(a101*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a104*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a105*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a106*h)
      y1 = y1 + ztmp*k6
      ztmp = dcmplx(a107*h)
      y1 = y1 + ztmp*k7
      ztmp = dcmplx(a108*h)
      y1 = y1 + ztmp*k8
      ztmp = dcmplx(a109*h)
      y1 = y1 + ztmp*k9
      y1 = y1 + cr*y

      call derivs(y1,t+c10*h,np_input,npn_input,k10, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !> step 10
      ztmp = dcmplx(a111*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a114*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a115*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a116*h)
      y1 = y1 + ztmp*k6
      ztmp = dcmplx(a117*h)
      y1 = y1 + ztmp*k7
      ztmp = dcmplx(a118*h)
      y1 = y1 + ztmp*k8
      ztmp = dcmplx(a119*h)
      y1 = y1 + ztmp*k9
      ztmp = dcmplx(a1110*h)
      y1 = y1 + ztmp*k10
      y1 = y1 + cr*y

      call derivs(y1,t+c11*h,np_input,npn_input,k2, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !> step 11
      ztmp = dcmplx(a121*h)
      y1 = ztmp*k1
      ztmp = dcmplx(a124*h)
      y1 = y1 + ztmp*k4
      ztmp = dcmplx(a125*h)
      y1 = y1 + ztmp*k5
      ztmp = dcmplx(a126*h)
      y1 = y1 + ztmp*k6
      ztmp = dcmplx(a127*h)
      y1 = y1 + ztmp*k7
      ztmp = dcmplx(a128*h)
      y1 = y1 + ztmp*k8
      ztmp = dcmplx(a129*h)
      y1 = y1 + ztmp*k9
      ztmp = dcmplx(a1210*h)
      y1 = y1 + ztmp*k10
      ztmp = dcmplx(a1211*h)
      y1 = y1 + ztmp*k2
      y1 = y1 + cr*y

      call derivs(y1,t+h,np_input,npn_input,k3, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)

      !> step 12
      ztmp = dcmplx(b1)
      k4 = ztmp*k1
      ztmp = dcmplx(b6)
      k4 = k4 + ztmp*k6
      ztmp = dcmplx(b7)
      k4 = k4 + ztmp*k7
      ztmp = dcmplx(b8)
      k4 = k4 + ztmp*k8
      ztmp = dcmplx(b9)
      k4 = k4 + ztmp*k9
      ztmp = dcmplx(b10)
      k4 = k4 + ztmp*k10
      ztmp = dcmplx(b11)
      k4 = k4 + ztmp*k2
      ztmp = dcmplx(b12)
      k4 = k4 + ztmp*k3
      k5 = y
      ztmp = dcmplx(h)
      k5 = k5 + ztmp*k4

      ! set error-weight-vector according to desired accuracy
      call setewv(y,k5,wt, dgldim_input)

      !> error step 1
      ztmp = dcmplx(e1)
      kerr = ztmp*k1
      ztmp = dcmplx(e6)
      kerr = kerr + ztmp*k6
      ztmp = dcmplx(e7)
      kerr = kerr + ztmp*k7
      ztmp = dcmplx(e8)
      kerr = kerr + ztmp*k8
      ztmp = dcmplx(e9)
      kerr = kerr + ztmp*k9
      ztmp = dcmplx(e10)
      kerr = kerr + ztmp*k10
      ztmp = dcmplx(e11)
      kerr = kerr + ztmp*k2
      ztmp = dcmplx(e12)
      kerr = kerr + ztmp*k3

      call wgtnrm(kerr,wt,err, dgldim_input)

      !> error step 2
      ztmp = dcmplx(-bh1)
      kerr2 = k4 + ztmp*k1
      ztmp = dcmplx(-bh2)
      kerr2 = kerr2 + ztmp*k9
      ztmp = dcmplx(-bh3)
      kerr2 = kerr2 + ztmp*k3
      call wgtnrm(kerr2,wt,err2, dgldim_input)
      err = err**2
      deno = err + 0.010_dp*(err2**2)
      err  = h*err/dsqrt(deno)
      fac1 = err**expo1
      fac  = fac1/facold**beta
      fac  = max(1.0_dp/facinc, min(1.0_dp/facdec, fac/0.90_dp))
      hnew = h/fac
      if (err < 1.0_dp) then

        write(31,*) t*au2fs, h*au2fs
        t = t + h
        h = hnew

        lreject = .false.
        if ( dabs(dreal(imagt) - 1.0_dp) > thresh_time_zero ) then
          call mapformat(k5, 1_i64, dgldim_input, nrindep_input, &
          & nrorb_input, A_input, phi_input)
          call renorm(rho_input, nrorb_input, nrorb_spinorbital_input, A_input, nrindep_input,&
          & phi_input)
          call mapformat(k5, 0_i64, dgldim_input, nrindep_input, &
          & nrorb_input, A_input, phi_input)
        endif
        call derivs(k5,t,np_input,npn_input,k4, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input,&
    & rhoss_input, detl_input, allow1_input, phi_input,&
    & phi2_input, A_input, rho_input, hel_input, hel2_input, hel3_input, &
    & venmatmo_input)
        facold = max(err,1.0d-4)

      else
        hnew   = h/min(1.0_dp/facdec,fac1/0.90_dp)
        h = hnew

        lreject = .true.
        if (h < 1.0d-10) then
          write(*,*) "Time-step too small!!! tstep =  ", h
          write(*,*) "Halting MCEND run due to integrator instabilites"
          stop
        endif
      endif
    enddo

    psi_input = k5
    h = hnew

    return
    end subroutine

!> Main driver for the calculation of derivatives.
    subroutine derivs(psi_input, t, np_input, npn_input, dtpsi_input, nrindep_input, dgldim_input,&
    & nrorb_input, nrorb_spinorbital_input, max_nrindep_input, hvalrho_input, rhoss_input, &
    & detl_input, allow1_input, phi_input, phi2_input, A_input, rho_input,  hel_input, hel2_input, hel3_input, &
    & venmatmo_input)
      implicit none

      integer     :: dgldim_input, allow1_input_dim_1
      complex(dp) :: lA_input(nrindep_input*nrspf)
      complex(dp) :: psi_input(dgldim_input), dtpsi_input(dgldim_input)
      complex(dp) :: irhon(nrspf,nrspf)
      complex(dp) :: mfn(nrprimn,nrspf,nrspf)
      real(dp), allocatable, intent(inout) :: np_input(:)
      real(dp)    :: npn_input(nrspf)
      real(dp)    :: t, wtime_1
      integer     :: nrindep_input, nrorb_input, nrorb_spinorbital_input, &
      & max_nrindep_input, allow1_input(max_nrindep_input*(nel-1),0:6*nrorb_input), detl_input(nel*nrindep_input), &
      & hvalrho_input(20*nrorb_input*nrorb_input*nrindep_input), rhoss_input
      complex(dp) :: A_input(nrindep_input*nrspf), phi_input(nrprime, nrorb_input), hel_input(nrorb_input,nrorb_input,5+nrensp), &
      & hel2_input(nrorb_input,nrorb_input,nrorb_input,nrorb_input), hel3_input(nrprime, nrorb_input, nrprime, nrorb_input), &
      & rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input), phi2_input(nrprime, nrorb_input), &
      & venmatmo_input(nrorb_input, nrorb_input, nrprimn)

      call mapformat(psi_input, 1_i64, dgldim_input, nrindep_input, &
          & nrorb_input, A_input, phi_input)

      call hamelem(phi_input, nrorb_input, hel_input, hel2_input, hel3_input, &
        & venmatmo_input)

     call calc_rhoe(rho_input, nrorb_input, nrorb_spinorbital_input, nrindep_input, &
     & hvalrho_input, rhoss_input, A_input  )

        wtime_1 = omp_get_wtime()
        call calc_rhoe2()

      if (.not. allocated(irho_input)) allocate(irho_input(nrorb_spinorbital_input,nrorb_spinorbital_input))
      if (.not. allocated(mf1_input)) allocate(mf1_input(nrorb_input,nrorb_input,nrprime,nrprime))
      call invert_rhoe(irho_input, np_input, rho_input,  nrorb_spinorbital_input )

      if (nel > 1) then
        allow1_input_dim_1 = max_nrindep_input*(nel-1)
      else
        allow1_input_dim_1 = 1
      end if

      call calc_meanfe(mf1_input, rho_input, nrorb_input, nrorb_spinorbital_input, &
    & A_input, nrindep_input, allow1_input, allow1_input_dim_1)

      call calc_rhon(A_input,nrindep_input)
      call invert_rhon(irhon,npn_input)
      call calc_meanfn(mfn, nrindep_input, venmatmo_input, nrorb_input, &
      & A_input, allow1_input, allow1_input_dim_1)

      call dadt(t, lA_input, nrindep_input, detl_input, hel_input, nrorb_input, &
      & hel2_input, venmatmo_input, A_input)
      call dphiedt(t,irho_input, nrorb_input,  phi_input, phi2_input, &
      & mf1_input)
      call dphindt(t,mfn,irhon)
      call mapformat2(lA_input,dtpsi_input,0_i64, phi2_input, nrindep_input, &
      & dgldim_input, nrorb_input )

      return
    end subroutine

    !> EOM for the A elements.
    subroutine dadt(t, lA_input, nrindep_input, detl_input, hel_input, nrorb_input, hel2_input, venmatmo_input, A_input)
      !cw> passing argument lA_input to support both spatial and spinorbital
      !cw> and not using too many flag_spinorbital/if control
      !cw> not sure if there is more elegant solution
      implicit none
      complex(dp) :: summ(nrspf,nrspf)
      complex(dp) :: lA_input(nrindep_input*nrspf)
      complex(dp) :: scme1, scme2
      complex(dp) :: enmat(nrindep_input,nrindep_input,nrprimn)
      complex(dp) :: scme1mu(3), scme3(nrprimn)
      complex(dp) :: neham(nrindep_input,nrindep_input)
      real(dp) ::  ef(3), t
      integer :: idet, jdet, vorz, nrdiff, ms, ps, ns, qs
      integer :: aa, alpha, l, j,  sum_sz_tmp
      integer :: nrindep_input, nrorb_input
      integer :: detl_input(nel*nrindep_input)
      complex(dp) :: hel_input(nrorb_input,nrorb_input,5+nrensp), hel2_input(nrorb_input,nrorb_input,nrorb_input,nrorb_input), &
      & venmatmo_input(nrorb_input, nrorb_input, nrprimn), A_input(nrindep_input*nrspf)
      real(dp) :: wtime_1, wtime_2

      do idet=1,nrindep_input*nrspf
        lA_input(idet) = c0
      enddo

      call getefield(t,ef)
      ! electronic Hamiltonian
      g_counter = g_counter + 1
      wtime_1 = omp_get_wtime()
      do idet=1, nrindep_input
        sum_sz_tmp = detl_nsz(idet, nrindep_input, detl_input)
        if (sum_sz_tmp .ne. nsz) then
          cycle
        end if
        do jdet=1, nrindep_input
          sum_sz_tmp = detl_nsz(jdet, nrindep_input, detl_input)
          if (sum_sz_tmp .ne. nsz) then
            cycle
          end if
          wtime_2 = omp_get_wtime()
          call maxcoinc(idet,jdet,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)
          scme1   = c0
          scme1mu = c0
          scme2   = c0
          if (nrdiff>=3) then
            neham(idet,jdet) = scme1 + sum(scme1mu(:)*ef(:)) + scme2
            cycle
          end if
          if (print_level==1) then
            write (*,*)  'time for maxcoinc in dadt', omp_get_wtime() - wtime_2
          end if
          wtime_2 = omp_get_wtime()
          if (nrdiff<=1) then
            call sc1t(nrdiff,idet,ms,ps,vorz,scme1, detl_input, nrindep_input, hel_input,&
            & nrorb_input)
            if (print_level==1) then
              write (*,*)  'time for sc1t in dadt', omp_get_wtime() - wtime_2
            end if
            wtime_2 = omp_get_wtime()
            call sc1mu(nrdiff,idet,ms,ps,vorz,scme1mu, detl_input, nrindep_input, hel_input,&
            &nrorb_input)
            if (print_level==1) then
              write (*,*)  'time for sc1mu in dadt', omp_get_wtime() - wtime_2
            end if
            wtime_2 = omp_get_wtime()
          end if
  ! Mar-9: NB: Q: will nrdiff=0 lead to linear dependence problem?
          call sc2(nrdiff,idet,jdet,ms,ns,ps,qs,vorz,scme2, detl_input, nrindep_input,&
  & hel2_input, nrorb_input)
          if (print_level==1) then
            write (*,*)  'time for sc2 in dadt', omp_get_wtime() - wtime_2
          end if
          wtime_2 = omp_get_wtime()
          neham(idet,jdet) = scme1 + sum(scme1mu(:)*ef(:)) + scme2
        enddo
      enddo

      if (print_level==1) then
        write (*,*)  'time for part1 in dadt', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()
      do idet=1, nrindep_input
        sum_sz_tmp = detl_nsz(idet, nrindep_input, detl_input)
        if (sum_sz_tmp .ne. nsz) then
          cycle
        end if
        do jdet=1, nrindep_input
          sum_sz_tmp = detl_nsz(jdet, nrindep_input, detl_input)
          if (sum_sz_tmp .ne. nsz) then
            cycle
          end if
          do j=1, nrspf
            lA_input(idet+(j-1)*nrindep_input) = lA_input(idet+(j-1)*nrindep_input)                 &
                                     + neham(idet,jdet)*A_input(jdet+(j-1)*nrindep_input)
          enddo
        enddo
      enddo

      if (print_level==1) then
        write (*,*)  'time for part2 in dadt', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()

      ! nuclear Hamiltonian
      !>>> heln(:,:,5) is often only multiplied by ef(3), the z component of the efield,
      !>>> since orientation of the axis is always set to z
      do j=1, nrspf
        do l=1, nrspf
          do idet=1, nrindep_input
            sum_sz_tmp = detl_nsz(idet, nrindep_input, detl_input)
            if (sum_sz_tmp .ne. nsz) then
              cycle
            end if
            if (use_wcap) then
              lA_input(idet+(j-1)*nrindep_input) = lA_input(idet+(j-1)*nrindep_input) &
                 + (heln(j,l,4) + ef(3)*heln(j,l,5) + heln(j,l,6))*A_input(idet+(l-1)*nrindep_input)
            else
              lA_input(idet+(j-1)*nrindep_input) = lA_input(idet+(j-1)*nrindep_input) &
                 + (heln(j,l,4) + ef(3)*heln(j,l,5))*A_input(idet+(l-1)*nrindep_input)
            endif
          enddo
        enddo
      enddo

      ! V_en
      do idet=1, nrindep_input
        sum_sz_tmp = detl_nsz(idet, nrindep_input, detl_input)
        if (sum_sz_tmp .ne. nsz) then
          cycle
        end if
        do jdet=1, nrindep_input
          sum_sz_tmp = detl_nsz(jdet, nrindep_input, detl_input)
          if (sum_sz_tmp .ne. nsz) then
            cycle
          end if
          call maxcoinc(idet,jdet,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)
          if (nrdiff >= 2) then
             scme3 = c0
             enmat(idet,jdet,:) =  c0
             cycle
          end if
          call sc1en(nrdiff,idet,ms,ps,vorz,scme3, detl_input, nrindep_input, venmatmo_input,&
          & nrorb_input)
          do alpha=1, nrprimn
            enmat(idet,jdet,alpha) = scme3(alpha)
          enddo
        enddo
      enddo

      do idet=1, nrindep_input
        sum_sz_tmp = detl_nsz(idet, nrindep_input, detl_input)
        if (sum_sz_tmp .ne. nsz) then
          cycle
        end if
        do jdet=1, nrindep_input
          sum_sz_tmp = detl_nsz(jdet, nrindep_input, detl_input)
          if (sum_sz_tmp .ne. nsz) then
            cycle
          end if
          do j=1, nrspf
            do l=1, nrspf
              summ(j,l) = c0
            enddo
          enddo
          call maxcoinc(idet,jdet,vorz,nrdiff,ms,ns,ps,qs, detl_input, nrindep_input)
          if (nrdiff >= 2) then
             cycle
          end if
          do aa=1, nrprimn
            do j=1, nrspf
              do l=1, nrspf
                summ(j,l) = summ(j,l) + dconjg(phin(aa,j))*phin(aa,l)*enmat(idet,jdet,aa)
              enddo
            enddo
          enddo
          do j=1, nrspf
            do l=1, nrspf
              lA_input(idet+(j-1)*nrindep_input) = lA_input(idet+(j-1)*nrindep_input) &
                                        + summ(j,l)*A_input(jdet+(l-1)*nrindep_input)*dr
            enddo
          enddo
        enddo
      enddo

      do idet=1, nrindep_input
        sum_sz_tmp =  detl_nsz(idet, nrindep_input, detl_input)
        if (sum_sz_tmp .ne. nsz) then
          cycle
        end if
        do j=1, nrspf
          lA_input(idet+(j-1)*nrindep_input) = -ci*imagt*lA_input(idet+(j-1)*nrindep_input)
        enddo
      enddo

      if (print_level==1) then
        write (*,*)  'time for part3 in dadt', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()
      return
    end subroutine

    !> derivative for the electronic phi (EOM)
    subroutine dphiedt(t,irho_input, nrorb_input,  phi_input, phi2_input, mf1_input)

      implicit none

      complex(dp), allocatable, intent(inout) :: irho_input(:,:)
      complex(dp) :: lphi2_input(nrprime,nrorb_input)
      real(dp) :: t
      real(dp) :: ef(3)
      integer :: ir, mu, la, jr, i, j, nu, kr
      integer :: nrorb_input
      complex(dp) :: phi_input(nrprime, nrorb_input), phi2_input(nrprime, nrorb_input), &
      &  mf1_input(nrorb_input, nrorb_input, nrprime, nrprime), id_matrix_temp(nrorb_spinorbital,nrorb_spinorbital), &
      &  phi2_temp(nrprime,nrorb_input), lphi2_temp(nrprime,nrorb_input), lphi2_temp2(nrprime,nrorb_input)

      phi2_input  = c0
      lphi2_input = c0
      phi2_temp   = c0
      lphi2_temp  = c0
      lphi2_temp2  = c0
      id_matrix_temp = c0

      if (flag_spinorbital==1) then
        do i=1, nrprime
          do j=1, nrorb_input
            if ( mod(j, 2) == 1 .and. j < nrorb_input ) then
              if ( dabs( dreal( phi2_input(i,j) ) - dreal( phi2_input(i,j+1) ) ) > 1.e-9) then
                write (*,'(A,I3,I3)')  'detect alpha/beta wf difference in phi2_spinorbital in dphiedt', i, j
                write (*,'(F15.9,F15.9)')  dreal( phi2_input(i,j) ), dreal( phi2_input(i,j+1) )
              end if
            end if
          end do
        end do
      end if

      do la=1,nrprime
        do mu=1,nrprime
          do jr=1,nrorb_input
            do ir=1,nrorb_input
              phi2_input(mu,ir) = phi2_input(mu,ir) + mf1_input(ir,jr,mu,la) * phi_input(la,jr)
            enddo
          enddo
        enddo
      enddo

      if (flag_spinorbital == 0) then
        do ir=1,nrorb_input
          do jr=1,nrorb_input
            do mu=1,nrprime
              lphi2_input(mu,jr)  = lphi2_input(mu,jr) + irho_input(2*jr,2*ir)*phi2_input(mu,ir)
            enddo
          enddo
        enddo

        if (mf_reduction) then
          do ir=1,nrorb_input
            do nu=1,nrprime
              do mu=1,nrprime
                lphi2_temp2(mu,ir) = lphi2_temp2(mu,ir) + tmat(mu,nu)*phi_input(nu,ir)
              enddo
            enddo
          enddo

          lphi2_input = lphi2_input + lphi2_temp2

          if (print_level == 1) then
            write (*,*) 'rho matrix in dphiedt'
            do ir=1,nrorb*2
              do jr=1,nrorb*2
                write (*,*) ir, jr, dreal(rho(ir,jr))
              enddo
            enddo

            do ir=1,nrorb*2
              do jr=1,nrorb*2
                do kr=1, nrorb*2
                  id_matrix_temp(ir,jr) = id_matrix_temp(ir,jr) + irho_input(ir,kr) * rho(kr,jr)
                enddo
              enddo
            enddo

            write (*,*) 'id matrix'
            do ir=1,nrorb*2
              do jr=1,nrorb*2
                write (*,*) ir, jr, dreal(id_matrix_temp(ir,jr))
              enddo
            enddo
          end if
        end if
      end if

      if (flag_spinorbital == 1) then
        do ir=1,nrorb_input
          do jr=1,nrorb_input
            do mu=1,nrprime
              lphi2_input(mu,jr) = lphi2_input(mu,jr) + irho_input(jr,ir)*phi2_input(mu,ir)
            enddo
          enddo
        enddo

        if (mf_reduction) then
          do ir=1,nrorb_input
            do nu=1,nrprime
              do mu=1,nrprime
                lphi2_temp2(mu,ir) = lphi2_temp2(mu,ir) + tmat(mu,nu)*phi_input(nu,ir)
              enddo
            enddo
          enddo
          lphi2_input = lphi2_input + lphi2_temp2
        end if
      end if

      call getefield(t,ef)
      phi2_input(:,:) = c0

      do jr=nrfrorb+1, nrorb_input
        do mu=1, nrprime
          phi2_input(mu,jr) = c0
          do la=1, nrprime
            phi2_input(mu,jr) = phi2_input(mu,jr) + (sum(ef(:)*x(:,mu,la)) + ci*vcap(mu,la))*phi_input(la,jr)
          enddo
        enddo
      enddo

      do ir=nrfrorb+1, nrorb_input
        do mu=1, nrprime
          phi2_input(mu,ir) = phi2_input(mu,ir) + lphi2_input(mu,ir)
        enddo
      enddo

      call projectione(phi_input, phi2_input, nrorb_input)

      do ir=nrfrorb+1, nrorb_input
        do mu=1, nrprime
          phi2_input(mu,ir) = -ci*imagt*phi2_input(mu,ir)
        enddo
      enddo

      return
    end subroutine

    !> derivative of nuclear phi, (EOM) ->  d/dt phi_n
    subroutine dphindt(t,mfn,irhon)

      implicit none

      complex(dp) :: lphin2(nrprimn,nrspf)
      complex(dp) :: mfn(nrprimn,nrspf,nrspf), irhon(nrspf,nrspf)
      complex(dp) :: phibuff(nrprimn), phibuff2(nrprimn)
      real(dp) :: t
      real(dp) :: ef(3)
      integer :: i, j, ix, ir, jr !, mu
      integer :: fftw_forward=-1, fftw_backward=1
      integer(i64) :: planf, planb

      call dfftw_plan_dft_1d(planf,nrprimn,phibuff,phibuff2,fftw_forward,0)
      call dfftw_plan_dft_1d(planb,nrprimn,phibuff2,phibuff,fftw_backward,0)

      phin2(:,:) = c0
      lphin2(:,:) = c0

      do i=1,nrspf
        do j=1,nrspf
          do ix=1,nrprimn
            phin2(ix,i) = phin2(ix,i) + mfn(ix,i,j)*phin(ix,j)
          enddo
        enddo
      enddo

      do ir=1,nrspf
        do jr=1,nrspf
          do ix=1,nrprimn
          ! Eq. (14) in MCEND 2012 paper
            lphin2(ix,jr) = lphin2(ix,jr)+irhon(jr,ir)*phin2(ix,ir)
          enddo
        enddo
      enddo

      phin2(:,:) = c0

      do i=1,nrspf
        do ix=1,nrprimn
          phibuff(ix) = phin(ix,i)
        enddo
        call dfftw_execute_dft(planf,phibuff,phibuff2)
        do ix=1,nrprimn
          phibuff2(ix) = kx(ix)**2*phibuff2(ix)
        enddo
        call dfftw_execute_dft(planb,phibuff2,phibuff)
        do ix=1, nrprimn
          phin2(ix,i) = phibuff(ix)/(2.0_dp*nrprimn*massn)
        enddo
      enddo

      call getefield(t,ef)

      !> calculate nuclear dipole moment
      do j=1, nrspf
        do ix=1, nrprimn
          ! SDI: nuclear dipole moment - nuclear charge
          !> ISOMASS change
          if (use_wcap) then
            phin2(ix,j) = phin2(ix,j) + N1*N2*phin(ix,j)/r(ix) &
                        - r(ix)*(N2*m1 - N1*m2)*ef(3)*phin(ix,j)/(m1 + m2) + wcap(ix)*phin(ix,j)
          else
            phin2(ix,j) = phin2(ix,j) + N1*N2*phin(ix,j)/r(ix) &
                        - r(ix)*(N2*m1 - N1*m2)*ef(3)*phin(ix,j)/(m1 + m2)
          endif
        enddo
      enddo

      phin2(:,:) = phin2(:,:) + lphin2(:,:)

      call projectionn()

      do i=1, nrspf
        do ix=1, nrprimn
          ! no motion of nuclear SPF
          ! don't allow nuclear orbitals to change with time
          if (freeze_nuc) then
            phin2(ix,i) = c0
          else
            phin2(ix,i) = -ci*imagt*phin2(ix,i)*dr
          endif
        enddo
      enddo

      call dfftw_destroy_plan(planf)
      call dfftw_destroy_plan(planb)

      return
    end subroutine

    !> electronic projection operator
    subroutine projectione(phi_input, phi2_input, nrorb_input)

      implicit none

      complex(dp) :: phi3_input(nrprime,nrorb_input), phi_input(nrprime,nrorb_input), &
  & phi2_input(nrprime,nrorb_input), ovl_alpha(nrorb_alpha,nrorb_alpha), ovl_beta(nrorb_beta, nrorb_beta), &
& ovl(nrorb_input,nrorb_input)
      complex(dp) :: h1(nrorb_input,nrorb_input), h2(nrorb_input,nrorb_input)
      integer :: i, j, k, mu
      integer :: nrorb_input

      if (flag_spinorbital == 0) then
        do i=1, nrorb_input
          do j=1, nrorb_input
            ovl(i,j) = c0
            do mu=1, nrprime
              ovl(i,j) = ovl(i,j) + dconjg(phi_input(mu,i))*phi_input(mu,j)
            enddo
          enddo
        enddo
      end if

      if (flag_spinorbital == 1) then
        do i=1, nrorb_beta
          do j=1, nrorb_beta
            ovl_alpha(i,j) = c0
            do mu=1, nrprime
              ovl_alpha(i,j) = ovl_alpha(i,j) + dconjg(phi_input(mu,2*i-1))*phi_input(mu,2*j-1)
            enddo
          enddo
        enddo

        do i=1, nrorb_beta
          do j=1, nrorb_beta
            ovl_beta(i,j) = c0
            do mu=1, nrprime
              ovl_beta(i,j) = ovl_beta(i,j) + dconjg(phi_input(mu,2*i))*phi_input(mu,2*j)
            enddo
          enddo
        enddo

      end if

      if (flag_spinorbital == 0) then
        call invovle(ovl,nrorb_input)
      end if

      if (flag_spinorbital == 1) then
        call invovle(ovl_alpha,nrorb_alpha)
        call invovle(ovl_beta,nrorb_beta)
      end if

      if (flag_spinorbital == 0) then
        do i=1, nrorb_input
          do j=1, nrorb_input
            h1(i,j) = c0
            do mu=1, nrprime
              h1(i,j) = h1(i,j) + dconjg(phi_input(mu,i))*phi2_input(mu,j)
            enddo
          enddo
        enddo
      end if

      if (flag_spinorbital == 1) then
        do i=1, nrorb_input
          do j=1, nrorb_input
            h1(i,j) = c0
            do mu=1, nrprime
              h1(i,j) = h1(i,j) + dconjg(phi_input(mu,i))*phi2_input(mu,j) * deltaspin(i,j)
            enddo
          enddo
        enddo
      end if

      if (flag_spinorbital == 0) then
        do i=1, nrorb_input
          do j=1, nrorb_input
            h2(i,j) = c0
            do k=1, nrorb_input
              h2(i,j) = h2(i,j) + ovl(i,k)*h1(k,j)
            enddo
          enddo
        enddo
      end if

  ! ref. 2.151 - 2.152 in Eur. Phys. J. Special Topics 223, 177–336 (2014)
      if (flag_spinorbital == 1) then
        do i=1, nrorb_input
          do j=1, nrorb_input
            h2(i,j) = c0

            do k=1, nrorb_input
              if (mod(k, 2) == 1  .and.  mod(i, 2) == 1  .and. mod(j, 2) == 1   ) then
                h2(i,j) = h2(i,j) + ovl_alpha((i+1)/2,(k+1)/2)*h1(k,j) * deltaspin(i,j)  * deltaspin(k,j)
              end if
              if (mod(k, 2) == 0 .and.  mod(i, 2) == 0 .and.  mod(j, 2) == 0 ) then
                h2(i,j) = h2(i,j) + ovl_beta((i)/2,(k)/2)*h1(k,j)      * deltaspin(i,j)   * deltaspin(k,j)
              end if
            enddo
          enddo
        enddo
      end if

      if (flag_spinorbital == 0) then
        do i=1,nrorb_input
          do mu=1,nrprime
          phi3_input(mu,i) = c0
            do j=1,nrorb_input
              phi3_input(mu,i) = phi3_input(mu,i) + phi_input(mu,j)*h2(j,i)
            enddo
          enddo
        enddo
      end if

     if (flag_spinorbital == 1) then
       do i=1,nrorb_input
         do mu=1,nrprime
           phi3_input(mu,i) = c0
           do j=1,nrorb_input
             phi3_input(mu,i) = phi3_input(mu,i) + phi_input(mu,j)*h2(j,i) * deltaspin(i,j)
           enddo
         enddo
       enddo
     end if

      phi2_input(:,nrfrorb+1:nrorb_input) = phi2_input(:,nrfrorb+1:nrorb_input) &
  & - phi3_input(:,nrfrorb+1:nrorb_input)

      return
    end subroutine

    !> nuclear projection operator
    subroutine projectionn()
      implicit none

      complex(dp) :: phin3(nrprimn,nrspf), ovl(nrspf,nrspf)
      complex(dp) :: h1(nrspf,nrspf), h2(nrspf,nrspf)
      integer :: i, j, k, mu

      do i=1,nrspf
        do j=1,nrspf
          ovl(i,j) = c0
          do mu=1,nrprimn
            ovl(i,j) = ovl(i,j) + dconjg(phin(mu,i))*phin(mu,j)
          enddo
        enddo
      enddo

      call invovln(ovl)

      do i=1,nrspf
        do j=1,nrspf
          h1(i,j) = c0
          do mu=1,nrprimn
            h1(i,j) = h1(i,j) + dconjg(phin(mu,i))*phin2(mu,j)
          enddo
        enddo
      enddo

      do i=1, nrspf
        do j=1, nrspf
          h2(i,j) = c0
          do k=1, nrspf
            h2(i,j) = h2(i,j) + ovl(i,k)*h1(k,j)
          enddo
        enddo
      enddo

      do i=1, nrspf
        do mu=1, nrprimn
          phin3(mu,i) = c0
          do j=1 ,nrspf
            phin3(mu,i) = phin3(mu,i) + phin(mu,j)*h2(j,i)
          enddo
        enddo
      enddo

      phin2(:,:) = phin2(:,:) - phin3(:,:)

      return

    end subroutine

    !> invert overlap for electrons
    subroutine invovle(ovl,dim_ovl)
      implicit none

      integer :: lwork, info, i, j, k, dim_ovl
      complex(dp) :: ovl(dim_ovl,dim_ovl), work(2*dim_ovl), h(dim_ovl,dim_ovl)
      real(dp) ::  w(dim_ovl), rwork(3*dim_ovl)
      character ::  jobz, uplo

      jobz = "V"
      uplo = "U"
      lwork = 3*dim_ovl
      info = 0
      call zheev(jobz,uplo,dim_ovl,ovl,dim_ovl,w,work,lwork,rwork,info)

      if (info /= 0) then
        write(*,*) 'Error inverting electronic overlap'
        write(*,*) "info zheev(involve)= ", info
        stop
      endif

      do i=1, dim_ovl
        do j=1, dim_ovl
          h(i,j) = c0
          do k=1, dim_ovl
            h(i,j) = h(i,j) + dconjg(ovl(j,k))*ovl(i,k)/w(k)
          enddo
        enddo
      enddo

      do i=1, dim_ovl
        do j=1, dim_ovl
          ovl(i,j) = h(i,j)
        enddo
      enddo

      return
    end subroutine

!> Invert overlap for nuclei
    subroutine invovln(ovl)

      implicit none
      complex(dp) :: ovl(nrspf,nrspf), work(2*nrspf), h(nrspf,nrspf)
      real(dp) ::  w(nrspf), rwork(3*nrspf)
      integer :: lwork, info, i, j, k
      character ::  jobz, uplo

      jobz = "V"
      uplo = "U"
      lwork = 3*nrspf
      info=0
      call zheev(jobz,uplo,nrspf,ovl,nrspf,w,work,lwork,rwork,info)

      if (info /= 0) then
        write(*,*) 'Error inverting nuclear overlap'
        write(*,*) "info zheev(involvn)= ", info
        stop
      endif

      do i=1,nrspf
        do j=1,nrspf
          h(i,j) = c0
          do k=1,nrspf
            h(i,j) = h(i,j) + dconjg(ovl(j,k))*ovl(i,k)/w(k)
          enddo
        enddo
      enddo

      ovl(:,:) = h(:,:)

      return

    end subroutine

    !> calculates matrix elements of the Hamiltonian
    !> matrix elements with spatial molecular orbitals
    !> hel(i,j,1)  = <phi_i | tmat | \phi_j>
    !> hel(i,j,2-4) = <phi_i | x,y,z | phi_j>
    !> hel(i,j,5) = <phi_i | vcap | phi_j >
    !> this here changes, and becomes hel(i,j,6) <- nothing ever happens with these
    !> w/ or w/o the efield
    !> hel(i,j,6,...) = <phi_i | hmat - tmat | phi_j> (WITHOUT efield!)
    !> nuclear degrees of freedom:
    !> heln(i,j,1) = kinetic energy
    !> heln(i,j,2) = n-n repulsion
    !> heln(i,j,3) = momentum
    !> heln(i,j,4) = kinetic energy + n-n repulsion
    !> heln(i,j,5) = dipole
    !> heln(i,j,6) = wcap
    subroutine hamelem(phi_input, nrorb_input, hel_input, hel2_input, hel3_input, venmatmo_input)

      implicit none

      complex(dp), allocatable :: hel2a_input(:,:,:,:)
      complex(dp) :: phin3(nrprimn,nrspf)
      complex(dp) ::  dens_input
      complex(dp) :: phibuff(nrprimn), phibuff2(nrprimn)
      integer :: i, j, k, mu, nu, si, l, ix, iensp, nrorb_input
      integer :: fftw_forward=-1, fftw_backward=1
      integer(i64) :: planf, planb
      complex(dp) :: hel_input(nrorb_input,nrorb_input,5+nrensp), &
      & venmatmo_input(nrorb_input, nrorb_input, nrprimn), phi_input(nrprime,nrorb_input), &
      & hel2_input(nrorb_input, nrorb_input, nrorb_input, nrorb_input), hel3_input(nrprime, nrorb_input, nrprime, nrorb_input)
      real(dp) :: wtime_1, wtime_2, wtime_6
      real(dp) :: alpha, beta, beta_dgemm
      integer :: INCX, INCY
      real(dp), allocatable :: vee_vector_temp(:), re_hel_vector_temp(:), re_phi_vector(:), &
      & im_hel_vector_temp(:), re_phi_input(:,:), im_phi_input(:,:), re_vector_temp(:), im_vector_temp(:), &
      & re_hel2a(:,:,:,:), im_hel2a(:,:,:,:), re_hel2a_2(:,:,:,:), im_hel2a_2(:,:,:,:)
      complex(dp), allocatable :: vector_temp(:), hel_vector_temp(:) !, combine_hel2a(:,:,:,:), combine_hel2a_2(:,:,:,:)
      integer :: nthr, nstart, nend

      save dens_input
      !$omp threadprivate(dens_input)
      wtime_1 = omp_get_wtime()
      alpha = 1.0_dp
      beta  = 1.0_dp

      INCX = 1
      INCY = 1

      if ( .not. allocated(vee_vector_temp))  allocate(vee_vector_temp(nrprime))
      if ( .not. allocated(hel_vector_temp))  allocate(hel_vector_temp(nrprime))
      if ( .not. allocated(re_hel_vector_temp))  allocate(re_hel_vector_temp(nrprime))
      if ( .not. allocated(im_hel_vector_temp))  allocate(im_hel_vector_temp(nrprime))
      if ( .not. allocated(re_phi_input))  allocate(re_phi_input(nrprime,nrorb_input))
      if ( .not. allocated(im_phi_input))  allocate(im_phi_input(nrprime,nrorb_input))
      if ( .not. allocated(vector_temp))  allocate(vector_temp(nrprime))
      if ( .not. allocated(re_vector_temp))  allocate(re_vector_temp(nrprime))
      if ( .not. allocated(im_vector_temp))  allocate(im_vector_temp(nrprime))
      if ( .not. allocated(re_phi_vector)) allocate(re_phi_vector(nrprime))
      if ( .not. allocated(re_hel2a)) allocate(re_hel2a(nrprime,nrprime,nrprime,nrorb_input))
      if ( .not. allocated(im_hel2a)) allocate(im_hel2a(nrprime,nrprime,nrprime,nrorb_input))
      if ( .not. allocated(re_hel2a_2)) allocate(re_hel2a_2(nrprime,nrprime,nrprime,nrorb_input))
      if ( .not. allocated(im_hel2a_2)) allocate(im_hel2a_2(nrprime,nrprime,nrprime,nrorb_input))

      re_hel2a   = 0.0_dp
      im_hel2a   = 0.0_dp
      re_hel2a_2 = 0.0_dp
      im_hel2a_2 = 0.0_dp
      wtime_6 =  omp_get_wtime()
      re_phi_input = dreal(phi_input)
      im_phi_input = dimag(phi_input)
      re_hel_vector_temp = 0.0_dp
      im_hel_vector_temp = 0.0_dp
      beta_dgemm = 0.0_dp
      vector_temp = 0.0_dp

      call dfftw_plan_dft_1d(planf,nrprimn,phibuff,phibuff2,fftw_forward,0)
      call dfftw_plan_dft_1d(planb,nrprimn,phibuff2,phibuff,fftw_backward,0)

      wtime_2 =  omp_get_wtime() - wtime_1

      if (print_level == 1) then
        write (*,*) 'time for dfftw_plan_dft_1d', wtime_2
      end if

      wtime_1 = omp_get_wtime()

      if (.not. allocated(hel2a_input)) allocate(hel2a_input(nrprime,nrprime,nrprime,nrorb_input))

      !$omp parallel
      !$omp workshare
      hel_input(:,:,1:5+nrensp) = c0
      venmatmo_input(:,:,:) = c0
      !$omp end workshare



      if (flag_spinorbital == 0) then

      !$omp do private(i,j,nu,mu,ix,iensp) schedule(dynamic) !! reduction(+:hel_input,venmatmo_input)
        do i=1, nrorb_input
          do j=i, nrorb_input
            do nu=1, nrprime
              do mu=1, nrprime
                dens_input = dconjg(phi_input(mu,i))*phi_input(nu,j)
                hel_input(i,j,1)   = hel_input(i,j,1)   + dens_input*tmat(mu,nu)
                hel_input(i,j,2:4) = hel_input(i,j,2:4) + dens_input*x(1:3,mu,nu)
                hel_input(i,j,5)   = hel_input(i,j,5)   + dens_input*vcap(mu,nu)
                do iensp=1, nrensp
                  hel_input(i,j,5+iensp) = hel_input(i,j,5+iensp) + dens_input &
  & *(hmat(mu,nu,iensp) - tmat(mu,nu))
                enddo
                do ix=1, nrprimn
                  venmatmo_input(i,j,ix) = venmatmo_input(i,j,ix) + dens_input*venmat(mu,nu,ix)
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do
        
      end if

      if (flag_spinorbital == 1) then
       !$omp do private(i,j,nu,mu,ix,iensp) schedule(dynamic) !!reduction(+:hel_input,venmatmo_input)
        do i=1, nrorb_input
          do j=i, nrorb_input
            if ( mod( i-j, 2) /= 0 ) then
              cycle
            end if
            do nu=1, nrprime
              do mu=1, nrprime
                dens_input = dconjg(phi_input(mu,i))*phi_input(nu,j)
                hel_input(i,j,1)   = hel_input(i,j,1)   + dens_input*tmat(mu,nu)
                hel_input(i,j,2:4) = hel_input(i,j,2:4) + dens_input*x(1:3,mu,nu)
                hel_input(i,j,5)   = hel_input(i,j,5)   + dens_input*vcap(mu,nu)
                do iensp=1, nrensp
                  hel_input(i,j,5+iensp) = hel_input(i,j,5+iensp) + dens_input &
  & *(hmat(mu,nu,iensp) - tmat(mu,nu))
                enddo
                do ix=1, nrprimn
                  venmatmo_input(i,j,ix) = venmatmo_input(i,j,ix) + dens_input*venmat(mu,nu,ix)
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      end if

      !$omp do schedule(dynamic) private(k,i,j)
      do k=1, 5+nrensp
        do i=2, nrorb_input
          do j=1, i-1
            hel_input(i,j,k) = dconjg(hel_input(j,i,k))
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel

      if (print_level == 1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for hel step1', wtime_2
        wtime_1 = omp_get_wtime()
      end if


      !$omp parallel
      !$omp do schedule(dynamic) private(i,j,k)
      do k=1, nrprimn
        do i=2, nrorb_input
          do j=1, i-1
            venmatmo_input(i,j,k) = dconjg(venmatmo_input(j,i,k))
          enddo
        enddo
      enddo
      !$omp end do
      

      !$omp workshare
      hel2a_input(:,:,:,:)   = c0
      !$omp end workshare
      !$omp end parallel

      wtime_1 = omp_get_wtime()
      !$omp parallel private(nthr, nstart, nend)
      nthr = omp_get_thread_num()
      nstart = nthr * nrprime / omp_get_num_threads()
      nend = (nthr + 1) * nrprime / omp_get_num_threads()
        call DGEMM( 'N', 'N', nrprime*nrprime*(nend - nstart), nrorb_input, nrprime, alpha, veeao(1 , 1, 1+ nstart, 1), &
      & Size(veeao, Dim = 1)**3,  &
      & re_phi_input, Size( re_phi_input, Dim = 1 ), beta_dgemm, re_hel2a_2(1, 1, 1 + nstart, 1),  Size(re_hel2a_2, Dim = 1)**3 )
        call DGEMM( 'N', 'N', nrprime*nrprime*(nend - nstart), nrorb_input, nrprime, alpha, veeao(1 , 1, 1+ nstart, 1 ), &
      & Size(veeao, Dim = 1)**3, &
      & im_phi_input, Size( im_phi_input, Dim = 1 ), beta_dgemm, im_hel2a_2(1 , 1, 1+ nstart, 1), Size(im_hel2a_2, Dim = 1)**3 )
     !$omp end parallel
      
      hel2a_input = dcmplx(re_hel2a_2, im_hel2a_2)
      if (print_level == 1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for hel2a step2 parallel method 2', wtime_2
        wtime_1 = omp_get_wtime()
      end if

     !$omp parallel
      if (flag_spinorbital == 0) then
       !$omp do private(l, j, si, nu, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=1, l
              do mu=1, si
                hel3_input(mu,j,si,l) = c0
                do nu=1, nrprime
                  hel3_input(mu,j,si,l) = hel3_input(mu,j,si,l) + dconjg(phi_input(nu,j)) &
  & *hel2a_input(mu,nu,si,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=1, nrorb_input
              do mu=si+1, nrprime
                  hel3_input(mu,j,si,l) = hel3_input(si,j,mu,l)
              enddo
            enddo
          enddo
        enddo
    !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=l+1, nrorb_input
              do mu=1, si
                hel3_input(mu,j,si,l) = dconjg(hel3_input(mu,l,si,j))
              enddo
            enddo
          enddo
        enddo
      !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=l+1, nrorb_input
              do mu=si+1, nrprime
                hel3_input(mu,j,si,l) = dconjg(hel3_input(si,l,mu,j))
              enddo
            enddo
          enddo
        enddo
      !$omp end do

      end if

      if (flag_spinorbital == 1) then
     !$omp do private(l, j, si, nu, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=1, l
              if ( mod( j-l, 2) /= 0 ) then
                cycle
              end if
              do mu=1, si
                hel3_input(mu,j,si,l) = c0
                do nu=1, nrprime
                  hel3_input(mu,j,si,l) = hel3_input(mu,j,si,l) + dconjg(phi_input(nu,j)) &
  & *hel2a_input(mu,nu,si,l)
                enddo

              enddo
            enddo
          enddo
        enddo
      !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=1, nrorb_input
              do mu=si+1, nrprime
                hel3_input(mu,j,si,l) = hel3_input(si,j,mu,l)
              enddo
            enddo
          enddo
        enddo
      !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=l+1, nrorb_input
              do mu=1, si
                hel3_input(mu,j,si,l) = dconjg(hel3_input(mu,l,si,j))
              enddo
            enddo
          enddo
        enddo
      !$omp end do

     !$omp do private(l, j, si, mu) schedule(dynamic)
        do l=1, nrorb_input
          do si=1, nrprime
            do j=l+1, nrorb_input
              do mu=si+1, nrprime
                hel3_input(mu,j,si,l) = dconjg(hel3_input(si,l,mu,j))
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      end if
     !$omp end parallel


      
      if (print_level == 1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for hel3 step3', wtime_2
        wtime_1 = omp_get_wtime()
      end if




      !$omp parallel
      !$omp do private(l, k, j, mu, si)
      do l=1,nrorb_input
        do k=1,nrorb_input
          do j=1,nrorb_input
            do mu=1,nrprime
              hel2a_input(mu,j,k,l) = c0
              do si=1,nrprime
                hel2a_input(mu,j,k,l) = hel2a_input(mu,j,k,l) + phi_input(si,k) &
  & *hel3_input(mu,j,si,l)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel
     

      if (print_level == 1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for hel2a step4', wtime_2
        wtime_1 = omp_get_wtime()
      end if


      !$omp parallel     
      if (flag_spinorbital == 0) then
      !$omp do private(l,k,j,i,mu)
        do l=1, nrorb_input
          do k=1, nrorb_input
            do j=1, nrorb_input
              do i=1, nrorb_input
                hel2_input(i,j,k,l) = c0
                do mu=1, nrprime
                   hel2_input(i,j,k,l) = hel2_input(i,j,k,l) &
  & + dconjg(phi_input(mu,i))*hel2a_input(mu,j,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo
       !$omp end do
      end if

      if (flag_spinorbital == 1) then
      !$omp do private(l,k,j,i,mu)
        do l=1, nrorb_input
          do k=1, nrorb_input
            do j=1, nrorb_input
              do i=1, nrorb_input
                if ( mod(i-k,2) /= 0 ) then
                  cycle
                end if
                hel2_input(i,j,k,l) = c0
                do mu=1, nrprime
                   hel2_input(i,j,k,l) = hel2_input(i,j,k,l) &
  & + dconjg(phi_input(mu,i))*hel2a_input(mu,j,k,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      end if
      !$omp end parallel

      if (print_level==1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for hel2 step5', wtime_2
        wtime_1 = omp_get_wtime()
      end if

      deallocate(hel2a_input)

      do i=1, nrspf
        do ix=1, nrprimn
          phibuff(ix) = phin(ix,i)
        enddo
        call dfftw_execute_dft(planf,phibuff,phibuff2)
        do ix=1, nrprimn
          phibuff2(ix) = kx(ix)*phibuff2(ix)
        enddo
        call dfftw_execute_dft(planb,phibuff2,phibuff)
        do ix=1, nrprimn
          phin3(ix,i) = phibuff(ix)/nrprimn
        enddo
      enddo

      !$omp parallel do private(i,j,ix) schedule(dynamic)
      do i=1, nrspf
        do j=1, nrspf
          heln(i,j,1) = c0
          heln(i,j,2) = c0
          heln(i,j,3) = c0
          heln(i,j,5) = c0
          do ix=1,nrprimn
            heln(i,j,1) = heln(i,j,1) + dconjg(phin3(ix,i))*phin3(ix,j)/(2.0_dp*massn)
            heln(i,j,2) = heln(i,j,2) + N1*N2*dconjg(phin(ix,i))*phin(ix,j)/r(ix)
            heln(i,j,3) = heln(i,j,3) + dconjg(phin(ix,i))*phin3(ix,j)
            heln(i,j,5) = heln(i,j,5) - r(ix)*(N2*m1 - N1*m2)*dconjg(phin(ix,i))*phin(ix,j)/(m2 + m1)
          enddo
        enddo
      enddo
      !$omp end parallel do

      if (use_wcap) then
        heln(:,:,6) = c0
        do i=1, nrspf
          do j=1, nrspf
            do ix=1, nrprimn
              heln(i,j,6) = heln(i,j,6) + dconjg(phin(ix,i))*wcap(ix)*phin(ix,j)
            enddo
          enddo
        enddo
        heln(:,:,6) = heln(:,:,6)*dr
      endif

      heln(:,:,1:3) = heln(:,:,1:3)*dr
      heln(:,:,4) = heln(:,:,1) + heln(:,:,2)
      heln(:,:,5) = heln(:,:,5)*dr

      call dfftw_destroy_plan(planf)
      call dfftw_destroy_plan(planb)

      if (print_level==1) then
        wtime_2 =  omp_get_wtime() - wtime_1
        write (*,*) 'time for misc', wtime_2
      end if

      deallocate(vee_vector_temp,hel_vector_temp, re_hel_vector_temp, im_hel_vector_temp, vector_temp)
      deallocate(re_vector_temp,im_vector_temp, re_phi_vector)
      deallocate(re_phi_input, im_phi_input)
      deallocate(re_hel2a, im_hel2a)
      deallocate(re_hel2a_2, im_hel2a_2)

      return

    end subroutine


    !> calculating rdm2 (electronic reduced density matrix, second order)
    ! not sure if i should rename it
    ! not sure if i should put nuclear-electron density from calc_meanfe to here
    subroutine calc_rhoe2()

      implicit none

      complex(dp) :: summa
      integer :: js, ls, cis, ind1, ind2, vorz1
      integer :: ii, j
      integer :: ms, ns

      rdm2 = c0
      rdm2_spinorbital = c0

      do j=1, nrspf
        ii = 1
        do cis=1, rhoss2
          vorz1 = hvalrdm2(ii)
          ind1 = hvalrdm2(ii+1)
          ind2 = hvalrdm2(ii+2)
          js = hvalrdm2(ii+3)
          ls = hvalrdm2(ii+4)
          ms = hvalrdm2(ii+5)
          ns = hvalrdm2(ii+6)
          ii = ii + 7
          summa = vorz1*dconjg(A(ind1+(j-1)*nrindep))*A(ind2+(j-1)*nrindep)
          rdm2(js,ls,ms,ns) = rdm2(js,ls,ms,ns) + summa
        enddo
      enddo

      do j=1, nrspf
        ii = 1
        do cis=1, rhoss2_spinorbital
          vorz1 = hvalrdm2_spinorbital(ii)
          ind1 = hvalrdm2_spinorbital(ii+1)
          ind2 = hvalrdm2_spinorbital(ii+2)
          js = hvalrdm2_spinorbital(ii+3)
          ls = hvalrdm2_spinorbital(ii+4)
          ms = hvalrdm2_spinorbital(ii+5)
          ns = hvalrdm2_spinorbital(ii+6)
          ii = ii + 7
          summa = vorz1*dconjg(A_spinorbital(ind1+(j-1)*nrindep_spinorbital))*A_spinorbital(ind2+(j-1)*nrindep_spinorbital)
  ! 1 2 1 2 order
          rdm2_spinorbital(js,ls,ms,ns) = rdm2_spinorbital(js,ls,ms,ns) + summa
        enddo
      enddo

    end subroutine

!> Electronic reduced density matrix, first order
    subroutine calc_rhoe(rho_input, nrorb_input, nrorb_spinorbital_input, nrindep_input, hvalrho_input, &
    & rhoss_input, A_input  )

      implicit none

      complex(dp) :: summa
      integer :: js, ls, cis, ind1, ind2, vorz1 !i,
      integer :: ii, j
      integer ::  nrorb_input,  nrorb_spinorbital_input, nrindep_input, rhoss_input
      complex(dp) :: rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input), A_input(nrindep_input*nrspf)
      integer :: hvalrho_input(20*nrorb_input*nrorb_input*nrindep_input)

      do js=1, nrorb_spinorbital_input
        do ls=1, nrorb_spinorbital_input
          rho_input(js,ls) = c0
        enddo
      enddo

      do j=1, nrspf
        ii = 1
        do cis=1, rhoss_input
          vorz1 = hvalrho_input(ii)
          ind1 = hvalrho_input(ii+1)
          ind2 = hvalrho_input(ii+2)
          js = hvalrho_input(ii+3)
          ls = hvalrho_input(ii+4)
          ii = ii + 5
          summa = vorz1*dconjg(A_input(ind1+(j-1)*nrindep_input))*A_input(ind2+(j-1)*nrindep_input)
          rho_input(js,ls) = rho_input(js,ls) + summa
        enddo
      enddo

      do js=1, nrorb_spinorbital_input
        do ls=1, nrorb_spinorbital_input
          rho_input(js,ls) = rho_input(js,ls)/dble(nel)
        enddo
      enddo

      if (print_level == 1) then
        do js=1, nrorb_spinorbital_input
          do ls=1, nrorb_spinorbital_input
            write (*,*) 'rho_input', js, ls, rho_input(js,ls)
          enddo
        enddo
      end if

      return

    end subroutine

    !> Nuclear density matrix
   subroutine calc_rhon(A_input,nrindep_input)
      !cw> integrating out electronic dof from the coefficient tensor A

      implicit none

      integer :: j, l, idet
      integer :: nrindep_input
      complex(dp) :: A_input(nrindep_input*nrspf)

      rhon(:,:) = c0

      do j=1, nrspf
        do l=1, nrspf
          do idet=1, nrindep_input
            rhon(j,l) = rhon(j,l) + dconjg(A_input(idet+(j-1)*nrindep_input)) &
  & *A_input(idet+(l-1)*nrindep_input)
          enddo
        enddo
      enddo

    end subroutine

!> Invert electronic density matrix.
    subroutine invert_rhoe(irho_input, np_input, rho_input,  nrorb_spinorbital_input )

      implicit none
      real(dp), allocatable, intent(inout) :: np_input(:)
      complex(dp), allocatable, intent(inout) :: irho_input(:,:)
      complex(dp) :: rhod_input(nrorb_spinorbital_input)
      complex(dp), allocatable :: rhoh_input(:,:)
      real(dp) ::  epsreg_input, trace_input
      integer :: is, js, ls, ir, jr, kr
      integer ::  nrorb_spinorbital_input
      complex(dp) ::  rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input)
      ! lapack vars:
      complex(dp) :: work_input(2*nrorb_spinorbital_input)
      real(dp) ::  w_input(nrorb_spinorbital_input), rwork_input(3*nrorb_spinorbital_input)
      integer ::  info
      integer :: lwork_input
      character :: jobz, uplo
      complex(dp) :: id_matrix_temp(nrorb_spinorbital, nrorb_spinorbital)

      id_matrix_temp = c0
      jobz = "V"
      uplo = "U"
      lwork_input = 3*nrorb_spinorbital_input
      w_input(:) = 0.0_dp
      rwork_input(:) = 0.0_dp
      work_input(:) = c0
      info = 0
      if (.not. allocated(rhoh_input)) allocate(rhoh_input(nrorb_spinorbital_input, nrorb_spinorbital_input))
      rhoh_input(:,:) = - rho_input(:,:)

      call zheev(jobz, uplo, nrorb_spinorbital_input, rhoh_input, nrorb_spinorbital_input, &
  & w_input, work_input, lwork_input, rwork_input, info)

      if (info /= 0) then
        write(*,*) "Info zheev(invert_rhoe): ", info
        stop
      endif

      do is=1, nrorb_spinorbital_input
        rhod_input(is) = c0
        do ls=1, nrorb_spinorbital_input
          rhod_input(is) = rhod_input(is) + dconjg(rhoh_input(is,ls))*rhoh_input(is,ls)*w_input(ls)
        enddo
      enddo

      np_input(:) = -dreal(rhod_input(:))

      trace_input = 0.0_dp
      do is=1, nrorb_spinorbital_input
        trace_input = trace_input + dreal(rho_input(is,is))
      enddo

      epsreg_input = epsreg0*trace_input
      ! regularization
      do is=1, nrorb_spinorbital_input
        w_input(is) = -w_input(is)
        if (w_input(is) < 64.0_dp*epsreg_input) then
          w_input(is) = w_input(is) + epsreg_input*dexp(-w_input(is)/epsreg_input)
        end if
      enddo

      w_input(:) = 1.0_dp/w_input(:)
      do is=1, nrorb_spinorbital_input
        do js=1, nrorb_spinorbital_input
          irho_input(is,js) = c0
          do ls=1, nrorb_spinorbital_input
            irho_input(is,js) = irho_input(is,js) + dconjg(rhoh_input(js,ls)) &
  & *rhoh_input(is,ls)*w_input(ls)
          enddo
        enddo
      enddo

      if (print_level == 2) then
        if (flag_spinorbital == 0) then
          do ir=1,nrorb*2
            do jr=1,nrorb*2
              do kr=1, nrorb*2
                id_matrix_temp(ir,jr) = id_matrix_temp(ir,jr) + irho_input(ir,kr) * rho(kr,jr)
              enddo
            enddo
          enddo

          do ir=1,nrorb*2
            do jr=1,nrorb*2
              write (*,*) 'inverse check', ir, jr,  id_matrix_temp(ir,jr)
            enddo
          enddo
        end if
      end if

      return
    end subroutine

!> Invert nuclear density matrix
    subroutine invert_rhon(irhon,npn)

      implicit none

      complex(dp) :: rhoh(nrspf,nrspf)
      complex(dp) :: irhon(nrspf,nrspf), rhod(nrspf)
      real(dp) ::  epsreg, trace, npn(nrspf)
      integer :: is, js, ls
      ! lapack vars:
      complex(dp) :: work(2*nrspf)
      real(dp) ::  w(nrspf), rwork(3*nrspf)
      integer :: lwork, info
      character ::  jobz, uplo

      jobz = "V"
      uplo = "U"
      lwork = 3*nrspf
      rwork(:) = 0.0_dp
      w(:) = 0.0_dp
      work(:) = c0
      info = 0
      rhoh(:,:) = -rhon(:,:)

      call zheev(jobz,uplo,nrspf,rhoh,nrspf,w,work,lwork,rwork,info)

      if (info /= 0) then
        write(*,*) "Info zheev(invert_rhon): ", info
        stop
      endif

      do is=1, nrspf
        rhod(is) = c0
        do ls=1, nrspf
          rhod(is) = rhod(is) + dconjg(rhoh(is,ls))*rhoh(is,ls)*w(ls)
        enddo
      enddo

      npn(:) = -dreal(rhod(:))

      trace = 0.0_dp
      do is=1,nrspf
        trace = trace + dreal(rhon(is,is))
      enddo

      epsreg = epsreg0*trace

      ! regularization
      do is=1, nrspf
        w(is) = -w(is)
        if (w(is) < 64.0_dp*epsreg) w(is) = w(is) + epsreg*dexp(-w(is)/epsreg)
      enddo

      do is=1, nrspf
        w(is) = 1.0_dp/w(is)
      enddo

      do is=1, nrspf
        do js=1, nrspf
          irhon(is,js) = c0
          do ls=1, nrspf
           irhon(is,js) = irhon(is,js) + dconjg(rhoh(js,ls))*rhoh(is,ls)*w(ls)
          enddo
        enddo
      enddo

      return

    end subroutine

    !> Calculate the mean field for the electrons
    !> this is the main bottle neck, openmp has been able to tame him a bit
    subroutine calc_meanfe(mf1_input, rho_input, nrorb_input, nrorb_spinorbital_input, &
    & A_input, nrindep_input, allow1_input,  allow1_input_dim_1)

      implicit none

      complex(dp) ::  venmat_enmo(nrspf,nrspf,nrprime,nrprime), density_enmo(nrorb_input,nrorb_input,nrspf,nrspf)
      integer :: mu, la, jr, lr, ind1, ind2
      integer :: vorz, nu, jn, ln, j, mr2, pr2, mr3, pr3
      integer :: jr0, lr0, mr0, pr0
      integer :: jshf, i, iall, lall, alpha
      real(dp) :: wtime_1, wtime_2

      integer :: nrorb_input, nrorb_spinorbital_input, nrindep_input, &
      & allow1_input_dim_1
      integer :: allow1_input(allow1_input_dim_1,0:6*nrorb_input)
      complex(dp) ::  A_input(nrindep_input*nrspf), &
      & mf1_tmp_input(nrorb_input, nrorb_input, nrprime, nrprime), temp_rdm, &
      & rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input), &
      & mf1_t_tmp_input(nrorb_input, nrorb_input, nrprime, nrprime)
      complex(dp), allocatable :: mf1_input(:,:,:,:), mf1_2e_tmp_input(:,:,:,:)
      save vorz, ind1, ind2, jr, lr, mr2, pr2, mr3, pr3

      !$omp threadprivate(vorz, ind1, ind2, jr, lr, mr2, pr2, mr3, pr3)
      wtime_1 = omp_get_wtime()
      if ( .not. allocated(mf1_input) )    allocate( mf1_input(nrorb_input, nrorb_input, nrprime, nrprime) )
      if ( .not. allocated(mf1_2e_tmp_input) )    allocate( mf1_2e_tmp_input(nrorb_input, nrorb_input, nrprime, nrprime) )
      wtime_1 = omp_get_wtime()
      !$omp parallel

      !$omp workshare
      mf1_input = c0
      mf1_2e_tmp_input = c0
      mf1_t_tmp_input = c0
      mf1_tmp_input = c0
      venmat_enmo   = c0
      density_enmo  = c0
      !$omp end workshare
      !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime init', omp_get_wtime() - wtime_1
        wtime_1 = omp_get_wtime()
      end if


      wtime_1 = omp_get_wtime()
      !$omp parallel
   if (flag_spinorbital == 0) then
    !$omp do private(jr0, lr0, mu, la, mr0, pr0) schedule(dynamic) reduction(+: temp_rdm)
     do la = 1, nrprime
       do mu = 1, la
         do lr0 = 1, nrorb
           do jr0 = 1, lr0
            do pr0 = 1, nrorb
              do mr0 = 1, nrorb
                  mr2 = mr0*2-1
                  pr2 = pr0*2-1
                  mr3 = mr0*2
                  pr3 = pr0*2
                  temp_rdm = rdm2(jr0*2-1,mr2,lr0*2-1,pr2)  + rdm2(jr0*2-1,mr2,lr0*2-1,pr3) &
                  & + rdm2(jr0*2-1,mr3,lr0*2-1,pr2)  + rdm2(jr0*2-1,mr3,lr0*2-1,pr3)
                  mf1_2e_tmp_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,mu,la) + temp_rdm * &
                  & hel3(mu,mr0,la,pr0)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
     do la = 1, nrprime
       do mu = la+1, nrprime
         do lr0 = 1, nrorb
           do jr0 = 1, lr0
              mf1_2e_tmp_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,la,mu)
          enddo
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
     do la = 1, nrprime
       do mu = 1, la
         do lr0 = 1, nrorb
           do jr0 = lr0+1, nrorb
              mf1_2e_tmp_input(jr0,lr0,mu,la) = DConjg( mf1_2e_tmp_input(lr0,jr0,mu,la) )
          enddo
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
     do la = 1, nrprime
       do mu = la+1, nrprime
         do lr0 = 1, nrorb
           do jr0 = lr0+1, nrorb
              mf1_2e_tmp_input(jr0,lr0,mu,la) = DConjg(  mf1_2e_tmp_input(lr0,jr0,la,mu) )
          enddo
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
     do la = 1, nrprime
       do mu = 1, nrprime
         do lr0 = 1, nrorb
           do jr0 = 1, nrorb
             mf1_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,mu,la) /dble(nel)
          enddo
        enddo
      enddo
    enddo
    !$omp end do

    end if

    if (flag_spinorbital == 1) then
  ! overall 6 loops, 2 nrprime 4 nrorb_input <- for larger system, one may do some density-fitting/hyper tensor decomposition
  ! <t> simple product 4 loops: 2 nrprime 2 nrorb_input, can be migrated to dphiedt-> 3 loops
  ! <ven>  6 loops, 2 nrprime, 2 nrorb_input, 2 nrspf <- ==1 in most cases -> can be migrated to dphiedt-> 3 loops
  ! maybe dgemm will help, but the current contraction is not adject index
    !$omp do private(jr0, lr0, mu, la, mr0, pr0) schedule(dynamic) !!reduction(+: mf1_2e_tmp_input)
      do la = 1, nrprime
        do mu = 1, la
          do lr0 = 1, nrorb_spinorbital
            do jr0 = 1, lr0
              if ( mod(jr0,2) .ne. mod(lr0,2) ) then
                cycle
              end if
              do pr0 = 1, nrorb_spinorbital
                if ( lr0 == pr0 ) then
                  cycle
                end if
                do mr0 = 1, nrorb_spinorbital
                  if ( jr0 == mr0 ) then
                    cycle
                  end if
                  if ( mod(mr0,2) .ne. mod(pr0,2) ) then
                    cycle
                  end if
                   mf1_2e_tmp_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,mu,la) + &
                    & rdm2_spinorbital(jr0,mr0,lr0,pr0) * hel3_spinorbital(mu,mr0,la,pr0)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
     !$omp end do

      !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
      do la = 1, nrprime
        do mu = la+1, nrprime
          do lr0 = 1, nrorb_spinorbital
            do jr0 = 1, lr0
                mf1_2e_tmp_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,la,mu)
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
      do la = 1, nrprime
        do mu = 1, la
          do lr0 = 1, nrorb_spinorbital
            do jr0 = lr0+1, nrorb_spinorbital
                mf1_2e_tmp_input(jr0,lr0,mu,la) = dconjg( mf1_2e_tmp_input(lr0,jr0,mu,la) )
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
      do la = 1, nrprime
        do mu = la+1, nrprime
          do lr0 = 1, nrorb_spinorbital
            do jr0 = lr0+1, nrorb_spinorbital
                mf1_2e_tmp_input(jr0,lr0,mu,la) = dconjg( mf1_2e_tmp_input(lr0,jr0,la,mu)  )
            enddo
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do private(jr0, lr0, mu, la) schedule(dynamic)
      do la = 1, nrprime
        do mu = 1, nrprime
          do lr0 = 1, nrorb_spinorbital
            do jr0 = 1, nrorb_spinorbital
                mf1_input(jr0,lr0,mu,la) = mf1_2e_tmp_input(jr0,lr0,mu,la) /dble(nel)
            enddo
          enddo
        enddo
      enddo
      !$omp end do
    end if
    !$omp end parallel

    if (print_level == 1) then
      wtime_2 =  omp_get_wtime() - wtime_1
      write(*,*) 'time for mean-field rdm2 calc vee', wtime_2
      wtime_1 = omp_get_wtime()
    end if


    !$omp parallel
    ! mean field due to V_en + T
    if (flag_spinorbital == 0) then
      if (.not. mf_reduction) then
      !$omp do private(jr0, lr0, mu, nu) schedule(dynamic) !!reduction(+: mf1_input)
        do nu=1,nrprime
          do mu=1,nrprime
            do lr0=1,nrorb_input
              do jr0=1,nrorb_input
                mf1_input(jr0,lr0,mu,nu) = mf1_input(jr0,lr0,mu,nu) +  rho_input(jr0*2,lr0*2) * tmat(mu,nu)
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      end if
    end if

    if (flag_spinorbital == 1) then
      if (.not. mf_reduction ) then
        !$omp do private(jr0, lr0, mu, nu) schedule(dynamic) !! reduction(+: mf1_input)
        do jr0=1,nrorb_input
          do lr0=1,nrorb_input
            do mu=1,nrprime
              do nu=1,nrprime
                mf1_input(jr0,lr0,mu,nu) = mf1_input(jr0,lr0,mu,nu) + rho_input(jr0,lr0) * tmat(mu,nu)
              enddo
            enddo
          enddo
        enddo
        !$omp end do
      end if
    end if
    !$omp end parallel


    if (print_level == 1) then
      write (*,*) 'wtime t step-1', omp_get_wtime() - wtime_1
      wtime_1 = omp_get_wtime()
    end if


      if (flag_spinorbital==0) then
      !$omp parallel
      !$omp do private(jn, ln, mu, nu, alpha) schedule(dynamic) !!reduction(+:venmat_enmo)
      do nu=1, nrprime
        do mu=1, nrprime
          do ln=1, nrspf
            do jn=1, nrspf
              do alpha=1, nrprimn
                venmat_enmo(jn,ln,mu,nu) =  venmat_enmo(jn,ln,mu,nu) + dconjg(phin(alpha,jn)) &
                & *phin(alpha,ln)*venmat(mu,nu,alpha)
              end do
            end do
          end do
        end do
      end do
      !$omp end do
    !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime ven step-1', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()

      !$omp parallel
      !$omp do private(jshf, iall, lall, jn, ln) schedule(dynamic) reduction(+:density_enmo)
      do jshf=1, nrshf
        do iall=1, allow1_input(jshf,0)
          jr = (allow1_input(jshf,3*(iall-1)+1) - 1)/2 + 1
          ind1 = allow1_input(jshf,3*(iall-1)+2)
          do lall=1, allow1_input(jshf,0)
            lr =  (allow1_input(jshf,3*(lall-1)+1) - 1)/2 + 1
            ind2 = allow1_input(jshf,3*(lall-1)+2)
            vorz = allow1_input(jshf,3*(iall-1)+3)*allow1_input(jshf,3*(lall-1)+3)
            do jn=1, nrspf
              if ( ind1+(jn-1)*nrindep_input == 0  ) then
                cycle
              end if
              do ln=1, nrspf
                 if ( ind2+(ln-1)*nrindep_input == 0 ) then
                   cycle
                 end if
                density_enmo(jr,lr,jn,ln) =  density_enmo(jr,lr,jn,ln)  +  vorz*dconjg(A_input(ind1+(jn-1)*nrindep_input)) &
                                 &     *A_input(ind2+(ln-1)*nrindep_input)
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime ven step-2', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()


      !$omp parallel
      !$omp do private(i, j, mu, nu, jn, ln) schedule(dynamic) !!reduction(+:mf1_tmp_input)
      do nu=1, nrprime
        do mu=1, nrprime
          do j=1, nrorb_input
            do i=1,nrorb_input
              do ln=1, nrspf
                do jn=1, nrspf
                     mf1_tmp_input(i,j,mu,nu) =  mf1_tmp_input(i,j,mu,nu) + density_enmo(i,j,jn,ln) * venmat_enmo(jn,ln,mu,nu)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime ven step-3', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()

      !$omp parallel
      !$omp do private(i, j, mu, nu)  schedule(dynamic) !!reduction(+:mf1_input)
      do nu=1, nrprime
        do mu=1, nrprime
          do j=1, nrorb_input
            do i=1,nrorb_input
                 mf1_input(i,j,mu,nu) = mf1_input(i,j,mu,nu) + mf1_tmp_input(i,j,mu,nu)*dr/dble(nel)
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      !$omp end parallel

       if (print_level == 1) then
         write (*,*) 'wtime ven step-4', omp_get_wtime() - wtime_1
       end if
       wtime_1 = omp_get_wtime()
       end if

      !!! open-shell new section !!!
      if (flag_spinorbital == 1) then
      !$omp parallel
      !$omp do private(jn, ln, mu, nu, alpha) schedule(dynamic) !!reduction(+:venmat_enmo)
      do nu=1, nrprime
        do mu=1, nrprime
          do ln=1, nrspf
            do jn=1, nrspf
              do alpha=1, nrprimn
                venmat_enmo(jn,ln,mu,nu) =  venmat_enmo(jn,ln,mu,nu) + dconjg(phin(alpha,jn)) &
                & *phin(alpha,ln)*venmat(mu,nu,alpha)
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime ven_spinorbital step-1', omp_get_wtime() - wtime_1
      end if

      wtime_1 = omp_get_wtime()

      !$omp parallel
      !$omp do private(jshf, iall, lall, jn, ln) schedule(dynamic) reduction(+:density_enmo)
      do jshf=1, nrshf
        do iall=1, allow1_input(jshf,0)
          jr = allow1_input(jshf,3*(iall-1) + 1 )
          ind1 = allow1_input(jshf,3*(iall-1)+2)
          do lall=1, allow1_input(jshf,0)
            lr = allow1_input(jshf,3*(lall-1)+1)
            ind2 = allow1_input(jshf,3*(lall-1)+2)
            vorz = allow1_input(jshf,3*(iall-1)+3)*allow1_input(jshf,3*(lall-1)+3)
            if ( mod(jr-lr,2) /=0 ) then
              cycle
            end if
            do jn=1, nrspf
              do ln=1, nrspf
                density_enmo(jr,lr,jn,ln) =  density_enmo(jr,lr,jn,ln)  +  vorz*dconjg(A_input(ind1+(jn-1)*nrindep_input)) &
                                 &     *A_input(ind2+(ln-1)*nrindep_input)
              end do
            end do
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

      if (print_level == 1) then
        write (*,*) 'wtime ven_spinorbital step-2', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()

      !$omp parallel
      !$omp do private(i, j, mu, nu, jn, ln) schedule(dynamic) !!reduction(+:mf1_tmp_input)
      do nu=1, nrprime
        do mu=1, nrprime
          do j=1, nrorb_input
            do i=1,nrorb_input
              do ln=1, nrspf
                do jn=1, nrspf
                     mf1_tmp_input(i,j,mu,nu) =  mf1_tmp_input(i,j,mu,nu) &
                     & + density_enmo(i,j,jn,ln) * venmat_enmo(jn,ln,mu,nu)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      !$omp end do
      !$omp end parallel


      if (print_level == 1) then
        write (*,*) 'wtime ven_spinorbital step-3', omp_get_wtime() - wtime_1
      end if
      wtime_1 = omp_get_wtime()
      !$omp parallel
      !$omp do private(i, j, mu, nu)  schedule(dynamic) !!reduction(+:mf1_input)
         do nu=1, nrprime
           do mu=1, nrprime
             do j=1, nrorb_input
               do i=1,nrorb_input
                 mf1_input(i,j,mu,nu) = mf1_input(i,j,mu,nu) + mf1_tmp_input(i,j,mu,nu)*dr/dble(nel)
               enddo
             enddo
           enddo
         enddo
      !$omp end do
      !$omp end parallel

         if (print_level == 1) then
           write (*,*) 'wtime ven_spinorbital step-4', omp_get_wtime() - wtime_1
         end if
         wtime_1 = omp_get_wtime()

       end if

      deallocate(mf1_2e_tmp_input)

      return
    end subroutine

    !> calculate mean field for the nuclei
    subroutine calc_meanfn(mfn, nrindep_input, venmatmo_input, nrorb_input, A_input, allow1_input,&
    & allow1_input_dim_1)

      implicit none
      complex(dp) :: mfn(nrprimn,nrspf,nrspf)
      integer :: alpha
      integer :: vorz
      integer :: nrindep_input, nrorb_input
      complex(dp) :: venmatmo_input(nrorb_input,nrorb_input,nrprimn), A_input(nrindep_input*nrspf)
      complex(dp) ::  density_nemo(nrspf,nrspf,nrorb_input,nrorb_input)
      integer :: allow1_input_dim_1
      integer :: allow1_input(allow1_input_dim_1,0:6*nrorb_input)
      integer :: jshf, iall, lall, jn, ln, ind1, ind2, jr, lr, jr0, lr0
      complex(dp) :: mfn_temp(nrprimn,nrspf,nrspf)
      real(dp) :: wtime_1

      mfn(:,:,:) = c0
      if (.not. freeze_nuc) then
        wtime_1 = omp_get_wtime()
        density_nemo(:,:,:,:) = c0
        mfn_temp(:,:,:)       = c0
        if ( flag_spinorbital == 0 ) then
          do jshf=1, nrshf
            do iall=1, allow1_input(jshf,0)
              jr = (allow1_input(jshf,3*(iall-1)+1) - 1)/2 + 1
              ind1 = allow1_input(jshf,3*(iall-1)+2)
              do lall=1, allow1_input(jshf,0)
                lr =  (allow1_input(jshf,3*(lall-1)+1) - 1)/2 + 1
                ind2 = allow1_input(jshf,3*(lall-1)+2)
                vorz = allow1_input(jshf,3*(iall-1)+3)*allow1_input(jshf,3*(lall-1)+3)
                do jn=1, nrspf
                  if ( ind1+(jn-1)*nrindep_input == 0  ) then
                    cycle
                  end if
                  do ln=1, nrspf
                    if ( ind2+(ln-1)*nrindep_input == 0 ) then
                      cycle
                    end if
                    density_nemo(jn,ln,jr,lr) =  density_nemo(jn,ln,jr,lr)  +  vorz*dconjg(A_input(ind1+(jn-1)*nrindep_input)) &
                                 &     *A_input(ind2+(ln-1)*nrindep_input)
                  end do
                end do
              end do
            end do
          end do

          do alpha=1, nrprimn
            do jr=1, nrorb_input*2
              do lr=1, nrorb_input*2
                do ln=1, nrspf
                  do jn=1, nrspf
                    jr0 = (jr+1)/2
                    lr0 = (lr+1)/2
                    if (mod(jr-lr,2) /= 0 ) then
                      cycle
                    end if
                    mfn_temp(alpha,jn,ln) = mfn_temp(alpha,jn,ln) + density_nemo(jn,ln,jr0,lr0) * venmatmo_input(jr0,lr0,alpha)
                  enddo
                enddo
              enddo
            enddo
          enddo
        end if

        if ( flag_spinorbital == 1 ) then
          do jshf=1, nrshf
            do iall=1, allow1_input(jshf,0)
              jr = allow1_input(jshf,3*(iall-1)+1)
              ind1 = allow1_input(jshf,3*(iall-1)+2)
              do lall=1, allow1_input(jshf,0)
                lr =  allow1_input(jshf,3*(lall-1)+1)
                ind2 = allow1_input(jshf,3*(lall-1)+2)
                vorz = allow1_input(jshf,3*(iall-1)+3)*allow1_input(jshf,3*(lall-1)+3)
                do jn=1, nrspf
                  if ( ind1+(jn-1)*nrindep_input == 0  ) then
                    cycle
                  end if
                  do ln=1, nrspf
                    if ( ind2+(ln-1)*nrindep_input == 0 ) then
                      cycle
                    end if
                    density_nemo(jn,ln,jr,lr) =  density_nemo(jn,ln,jr,lr)  +  vorz*dconjg(A_input(ind1+(jn-1)*nrindep_input)) &
                                 &     *A_input(ind2+(ln-1)*nrindep_input)
                  end do
                end do
              end do
            end do
          end do

          do alpha=1, nrprimn
            do jr=1, nrorb_input
              do lr=1, nrorb_input
                do ln=1, nrspf
                  do jn=1, nrspf
                    jr0 = jr
                    lr0 = lr
                    if (mod(jr-lr,2) /= 0 ) then
                      cycle
                    end if
                    mfn_temp(alpha,jn,ln) = mfn_temp(alpha,jn,ln) + density_nemo(jn,ln,jr0,lr0) * venmatmo_input(jr0,lr0,alpha)
                  enddo
                enddo
              enddo
            enddo
          enddo
        end if
        mfn = mfn_temp
      endif

    end subroutine

!> Analysis of the reduced density matrix: Compute Shannon entropy.
    subroutine rdm1analysis(time, file_unit,  nrorb_spinorbital_input, rho_input)
      ! proto routine adapted from TDCI
      ! though not quite optimal as it requires the eigenvalues and vectors to be computed
      ! which they are in other routines so it shouldn't need to be used, just
      ! need to identify what variables those are
      !> ok so there is something else that does this, earlier, so this would just
      !> compute the shannon entro

      implicit none
      complex(dp), allocatable :: rho2_input(:,:)
      complex(dp), allocatable :: work_input(:)
      real(dp), allocatable :: w_input(:), rwork_input(:)
      real(dp) :: norm_input, entro_input, cor_input
      real(dp) :: time
      integer :: lwork_input, info_input, i, nso_input
      character :: jobz_input, uplo_input
      integer :: file_unit,  nrorb_spinorbital_input
      complex(dp) :: rho_input(nrorb_spinorbital_input, nrorb_spinorbital_input)
      allocate(work_input(12*nrorb_spinorbital_input))
      allocate(w_input(nrorb_spinorbital_input), rwork_input(18*nrorb_spinorbital_input))
      jobz_input = "V"
      uplo_input = "U"
      lwork_input = 18*nrorb_spinorbital_input
      allocate(rho2_input(nrorb_spinorbital_input,nrorb_spinorbital_input))
      nso_input = nrorb_spinorbital_input
      rho2_input(:,:) = rho_input(:,:)

      call zheev(jobz_input,uplo_input,nso_input,rho2_input, &
  & nso_input, w_input, work_input, lwork_input,rwork_input,info_input)
      norm_input = 0.0_dp
      do i=1, nso_input
        if (dabs(w_input(i)) < 1.0e-10_dp) then
          w_input(i) = 1.0e-10_dp
        endif
        norm_input = norm_input + w_input(i)
      enddo
      do i=1, nso_input
        w_input(i) = w_input(i)/(norm_input)
      enddo
      entro_input= 0.0_dp
      do i=1, nso_input
        if ( w_input(i) < 1.e-10_dp ) then
          cycle
        endif
        entro_input = entro_input - dlog(w_input(i))*w_input(i)
      enddo

      cor_input = 1.0_dp
      do i=1, nso_input
        cor_input = cor_input - w_input(i)**2/dble(nel)
      enddo
      write(file_unit,*) time, entro_input, cor_input
      deallocate(rho2_input,work_input,w_input,rwork_input)

    end subroutine

    subroutine wgtnrm(kerr, wt, err, dgldim_input)
      implicit none

      integer  :: dgldim_input
      complex(dp) :: kerr(dgldim_input)
      real(dp) :: tmp
      real(dp) :: err
      real(dp) :: wt(dgldim_input)
      integer :: i

      tmp = 0.0_dp
      do i=1, dgldim_input
        tmp = tmp + (cdabs(kerr(i))/wt(i))**2
      enddo
      err = dsqrt(tmp/dgldim_input)
      return
    end subroutine

    subroutine renorm(rho_input, nrorb_input, nrorb_spinorbital_input, A_input, nrindep_input, phi_input)
      implicit none

      complex(dp) :: cnorm
      real(dp) ::  norm
      integer :: i, js, j
      integer :: nrorb_input, nrorb_spinorbital_input, nrindep_input
      complex(dp) :: rho_input(nrorb_spinorbital_input,nrorb_spinorbital_input), A_input(nrindep_input*nrspf), &
      & phi_input(nrprime,nrorb_input)

      cnorm = c0
      do js=1,nrorb_spinorbital_input
        cnorm = cnorm + rho_input(js,js)
      enddo

      A_input(:) = A_input(:)/dsqrt(dreal(cnorm))

      !> norm phi_e
      do i=1, nrorb_input
        norm = 0.0_dp
        do j=1, nrprime
          norm = norm + dconjg(phi_input(j,i))*phi_input(j,i)
        enddo
        norm = dsqrt(norm)

        do j=1, nrprime
          phi_input(j,i) = phi_input(j,i)/norm
        enddo

      enddo

      !> norm phi_n
      do i=1, nrspf
        norm = 0.0_dp
        do j=1, nrprimn
          norm = norm + dconjg(phin(j,i))*phin(j,i)
        enddo
        norm = dsqrt(norm*dr)

        do j=1, nrprimn
          phin(j,i) = phin(j,i)/norm
        enddo
      enddo

      return
    end subroutine

    subroutine setewv(amat, bmat, cmat, dgldim_input)
      implicit none

      integer :: dgldim_input
      complex(dp) :: amat(dgldim_input), bmat(dgldim_input)
      real(dp) ::  cmat(dgldim_input)
      integer :: i

      do i=1, dgldim_input
        cmat(i) = tol + tol*max(cdabs(amat(i)), cdabs(bmat(i)), 4.44d-16)
      enddo

      return
    end subroutine

    subroutine acfphi(phi0, acfspf, nrorb_input, phi_input)

      implicit none

      complex(dp) :: phi0(nrprime,nrorb_input), acfspf(nrorb_input), phi_input(nrprime,nrorb_input)
      integer :: ix, ir, nrorb_input

      do ir=1, nrorb_input
        acfspf(ir) = c0
        do ix=1, nrprime
          acfspf(ir) = acfspf(ir) + dconjg(phi0(ix,ir))*phi_input(ix,ir)
        enddo
      enddo

      return
    end subroutine

      subroutine zcopy_wr(xin,xcopy)
        !> BLAS zcopy wrapper
        implicit none
        complex(dp) :: xin(dgldim), xcopy(dgldim)

        call zcopy(dgldim,xin,1_i64,xcopy,1_i64)

        return
      end

      subroutine zscal_wr(xscaler,xscaled)
        !> BLAS zscal wrapper
        implicit none
        complex(dp) :: xscaler, xscaled(dgldim)

        call zscal(dgldim,xscaler,xscaled,1_i64)

        return
      end

      subroutine zaxpy_wr(xscaler,xscaled,xout)
        !> BLAS zaxpy wrapper
        implicit none
        complex(dp) :: xscaler, xscaled(dgldim), xout(dgldim)

        call zaxpy(dgldim,xscaler,xscaled,1_i64,xout,1_i64)

        return
      end

  end module
  !> @file
  !> @brief contains the main propagation
  !! and evaluation of Hamiltonian elements/mean fields.

