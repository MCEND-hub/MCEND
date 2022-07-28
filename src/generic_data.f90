!> Module containing the definition of the data types.
!! @param tmat Electronic kinetic energy in AO basis.
!! @param hmat One-electron integrals in AO basis.
!! @param x Electronic position integrals in AO basis.
!! @param veeao Two-electron integrals in AO basis.
!! @param venmat Electron-nuclear attraction integrals in AO basis for the given nuclear positions.
!! @param vcap Complex absorbing potential on the electronic grid.
!! @param wcap Complex absorbing potential on the nuclear grid.
!! @param venmatmo Electron-nuclear attraction in MO basis.
!! @param A The A-vector that determines the weight of each configuration in the MCEND wave function.
!! @param phi The electronic wave function.
!! @param phin The nuclear wave function.
   module globalvars
   use params
   use omp_lib
   use inputvars

     real(dp), allocatable :: tmat(:,:)
     real(dp), allocatable :: hmat(:,:,:), x(:,:,:)
     real(dp), allocatable :: taux(:,:,:)
     real(dp), allocatable :: veeao(:,:,:,:)
     real(dp), allocatable :: venmat(:,:,:)
     real(dp), allocatable :: venmat_spinorbital(:,:,:)
     real(dp), allocatable :: vcap(:,:)
     complex(dp), allocatable :: wcap(:)
     complex(dp), allocatable :: venmatmo(:,:,:)
     complex(dp), allocatable :: venmatmo_spinorbital(:,:,:)
     complex(dp), allocatable :: A(:)
     complex(dp), allocatable :: A_spinorbital(:)
     complex(dp), allocatable :: phi(:,:)
     complex(dp), allocatable :: phi_alpha(:,:), phi_beta(:,:)
     complex(dp), allocatable :: phi_spinorbital(:,:)
     complex(dp), allocatable :: phi2(:,:)
     complex(dp), allocatable :: phi2_spinorbital(:,:)
     complex(dp), allocatable :: dens_mo(:,:,:,:)
     complex(dp), allocatable :: dens_mo_alpha(:,:,:,:),dens_mo_beta(:,:,:,:)
     complex(dp), allocatable :: dens_mo_spinorbital(:,:,:,:)
     complex(dp), allocatable :: mf1_t_temp(:,:,:,:)
     complex(dp), allocatable :: mf1_t_spinorbital_temp(:,:,:,:)
     complex(dp), allocatable :: phin(:,:)
     complex(dp), allocatable :: phin2(:,:)
     complex(dp), allocatable :: phin_spinorbital(:,:)
     complex(dp), allocatable :: phin2_spinorbital(:,:)
     complex(dp), allocatable :: hel(:,:,:)
     complex(dp), allocatable :: hel_spinorbital(:,:,:)
     complex(dp), allocatable :: hel2(:,:,:,:)
     complex(dp), allocatable :: hel2_spinorbital(:,:,:,:)
     complex(dp), allocatable :: hel3(:,:,:,:)
     complex(dp), allocatable :: hel3_spinorbital(:,:,:,:)
     complex(dp), allocatable :: heln(:,:,:)
     complex(dp), allocatable :: rho(:,:)
     complex(dp), allocatable :: rho_spinorbital(:,:)
     complex(dp), allocatable :: rho_temp(:,:)
     complex(dp), allocatable :: rho_spinorbital_temp(:,:)
     complex(dp), allocatable :: rdm2(:,:,:,:)
     complex(dp), allocatable :: rdm2_spinorbital(:,:,:,:)
     complex(dp), allocatable :: rhon(:,:)
     complex(dp), allocatable :: id_matrix(:,:)
     complex(dp), allocatable :: id_matrix_spinorbital(:,:)
     real(dp),    allocatable  :: r(:), rsp(:)
     real(dp),    allocatable  :: kx(:)
     real(dp),    allocatable  :: ssqmat(:,:)
     integer,     allocatable  :: detl(:)
     integer,     allocatable  :: detl_spinorbital(:)
     integer,     allocatable  :: hval3(:)
     integer,     allocatable  :: hval3_spinorbital(:)
     integer,     allocatable  :: shdl(:)
     integer,     allocatable  :: shdl_spinorbital(:)
     integer,     allocatable  :: shdl_so(:)
     integer,     allocatable  :: dhdl(:)
     integer,     allocatable  :: dhdl_spinorbital(:)
     integer,     allocatable  :: hvalrho(:)
     integer,     allocatable  :: hvalrho_alpha(:), hvalrho_beta(:)
     integer,     allocatable  :: hvalrho_spinorbital(:)
     integer,     allocatable  :: hvalrdm2(:)
     integer,     allocatable  :: hvalrdm2_spinorbital(:)
     integer,     allocatable  :: allow1(:,:)
     integer,     allocatable  :: allow1_spinorbital(:,:)
     integer,     allocatable  :: hmf(:,:)
     integer,     allocatable  :: hmf_spinorbital(:,:)
     integer,     allocatable  :: spin_sign_list(:)
     integer,     allocatable  :: list_detl(:), list_detl_spinorbital(:)
     integer                   :: g_counter, hamelem_counter, hel_counter, hel2_counter, &
 & hel2_transformation_counter, hel2a_counter,  hel2a2_counter, hel3_counter, t_counter, &
 & t_spinorbital_counter, t_general_counter, runrk8_counter, runrk8_spinorbital_counter, n_detl, n_detl_spinorbital,  &
 & heln_spinorbital_counter, heln_counter, runrk8_general_counter

     contains

!> Routine to allocate the arrays.
       subroutine allocate_arrays()

         implicit none
         if (.not. allocated(hmat))     allocate(hmat(nrprime,nrprime,nrensp))
         if (.not. allocated(tmat))     allocate(tmat(nrprime,nrprime))
         if (.not. allocated(x))        allocate(x(3,nrprime,nrprime))
         if (.not. allocated(veeao))    allocate(veeao(nrprime,nrprime,nrprime,nrprime))
         if (.not. allocated(venmat))   allocate(venmat(nrprime,nrprime,nrprimn))
         if (.not. allocated(venmatmo)) allocate(venmatmo(nrorb,nrorb,nrprimn))
         if (.not. allocated(venmatmo_spinorbital)) allocate(venmatmo_spinorbital(nrorb_spinorbital,nrorb_spinorbital,nrprimn))
         if (.not. allocated(venmat_spinorbital))   allocate(venmat_spinorbital(nrprime,nrprime,nrprimn))
         if (.not. allocated(hel))      allocate(hel(nrorb,nrorb,5+nrensp))
         if (.not. allocated(hel_spinorbital))       allocate(hel_spinorbital(nrorb_spinorbital,nrorb_spinorbital,5+nrensp))
         if (.not. allocated(hel2))     allocate(hel2(nrorb,nrorb,nrorb,nrorb))
         if (.not. allocated(hel2_spinorbital))      allocate(hel2_spinorbital(nrorb_spinorbital, &
         & nrorb_spinorbital,nrorb_spinorbital,nrorb_spinorbital))
         if (.not. allocated(hel3))     allocate(hel3(nrprime,nrorb,nrprime,nrorb))
         if (.not. allocated(hel3_spinorbital))      allocate(hel3_spinorbital(nrprime,nrorb_spinorbital,nrprime,nrorb_spinorbital))
         if (.not. allocated(taux))     allocate(taux(nrprime,nrprime,nrensp))
         if (.not. allocated(phi))      allocate(phi(nrprime,nrorb))
         if (.not. allocated(phi_alpha))      allocate(phi_alpha(nrprime,nrorb_spinorbital))
         if (.not. allocated(phi_beta))      allocate(phi_beta(nrprime,nrorb_spinorbital))
         if (.not. allocated(phi_spinorbital))      allocate(phi_spinorbital(nrprime,nrorb_spinorbital))
         if (.not. allocated(dens_mo))     allocate(dens_mo(nrorb,nrorb,nrprime,nrprime))
         if (.not. allocated(dens_mo_spinorbital))    &
 & allocate(dens_mo_spinorbital(nrorb_spinorbital,nrorb_spinorbital,nrprime,nrprime))
         if (.not. allocated(mf1_t_temp))  &
 &   allocate(mf1_t_temp(nrorb,nrorb,nrprime,nrprime))
         if (.not. allocated(mf1_t_spinorbital_temp))  &
 &   allocate(mf1_t_spinorbital_temp(nrorb_spinorbital,nrorb_spinorbital,nrprime,nrprime))
         if (.not. allocated(phi2))     allocate(phi2(nrprime,nrorb))
         if (.not. allocated(phi2_spinorbital))     allocate(phi2_spinorbital(nrprime,nrorb_spinorbital))
         if (.not. allocated(phin))     allocate(phin(nrprimn,nrspf))
         if (.not. allocated(phin2))    allocate(phin2(nrprimn,nrspf))
         if (.not. allocated(phin_spinorbital))     allocate(phin_spinorbital(nrprimn,nrspf))
         if (.not. allocated(phin2_spinorbital))    allocate(phin2_spinorbital(nrprimn,nrspf))
         if (.not. allocated(vcap))     allocate(vcap(nrprime,nrprime))
         if (.not. allocated(wcap))     allocate(wcap(nrprimn))
         if (.not. allocated(ssqmat))   allocate(ssqmat(nrindep,nrindep))
         if (.not. allocated(rsp))      allocate(rsp(nrensp))
         if (.not. allocated(r))        allocate(r(nrprimn))
         if (.not. allocated(kx))       allocate(kx(nrprimn))
         if (.not. allocated(hval3))    allocate(hval3(nrindep))
         if (.not. allocated(hval3_spinorbital))    allocate(hval3_spinorbital(nrindep_spinorbital))
         if (.not. allocated(detl))     allocate(detl(nel*nrindep))
         if (.not. allocated(detl_spinorbital))  allocate(detl_spinorbital(nel_spinorbital*nrindep_spinorbital))
         if (.not. allocated(shdl))     allocate(shdl(max_nrindep*(nel-1)))
         if (.not. allocated(dhdl))     allocate(dhdl(max_nrindep_2*(nel-2)))
         if (.not. allocated(shdl_spinorbital))     allocate(shdl_spinorbital(max_nrindep_spinorbital*(nel-1)))
         if (.not. allocated(dhdl_spinorbital))     allocate(dhdl_spinorbital(max_nrindep_2_spinorbital*(nel-1)))
         if (.not. allocated(hvalrho))  allocate(hvalrho(20*nrorb*nrorb*nrindep))
         if (.not. allocated(hvalrdm2))  allocate(hvalrdm2(60*nrorb*nrorb*nrorb*nrorb*nrindep))
         if (.not. allocated(hvalrho_spinorbital)) &
 & allocate(hvalrho_spinorbital(20*nrorb_spinorbital*nrorb_spinorbital*nrindep_spinorbital))
         if (.not. allocated(hvalrdm2_spinorbital)) &
 & allocate(hvalrdm2_spinorbital(20*nrorb_spinorbital*nrorb_spinorbital*nrorb_spinorbital*nrorb_spinorbital*nrindep_spinorbital))
         if (nel>1) then
           if (.not. allocated(allow1))   allocate(allow1(max_nrindep*(nel-1),0:6*nrorb))
           if (.not. allocated(allow1_spinorbital))  &
 & allocate(allow1_spinorbital(max_nrindep_spinorbital*(nel-1),0:6*nrorb_spinorbital))
         else if (nel == 1) then
           allocate(allow1(1,0:6*nrorb))
           allocate(allow1_spinorbital(1,0:6*nrorb))
         else
           write (*,*) 'number of electron less than 1, not supported'
           stop
         end if
         if (.not. allocated(A))        allocate(A(nrindep*nrspf))
         if (.not. allocated(A_spinorbital))        allocate(A_spinorbital(nrindep_spinorbital*nrspf))
         if (.not. allocated(rho))      allocate(rho(2*nrorb,2*nrorb))
         if (.not. allocated(rho_spinorbital))   allocate(rho_spinorbital(nrorb_spinorbital,nrorb_spinorbital))
         if (.not. allocated(rho_temp))      allocate(rho_temp(2*nrorb,2*nrorb))
         if (.not. allocated(rho_spinorbital_temp))   allocate(rho_spinorbital_temp(nrorb_spinorbital,nrorb_spinorbital))
         if (.not. allocated(id_matrix))      allocate(id_matrix(2*nrorb,2*nrorb))
         if (.not. allocated(id_matrix_spinorbital))      allocate(id_matrix_spinorbital(nrorb_spinorbital,nrorb_spinorbital))
         if (.not. allocated(rdm2))      allocate(rdm2(2*nrorb,2*nrorb,2*nrorb,2*nrorb))
         if (.not. allocated(rdm2_spinorbital))    &
         & allocate(rdm2_spinorbital(nrorb_spinorbital,nrorb_spinorbital,nrorb_spinorbital,nrorb_spinorbital))
         if (.not. allocated(rhon))     allocate(rhon(nrspf,nrspf))
         if (use_wcap) then
           if (.not. allocated(heln))   allocate(heln(nrspf,nrspf,6))
         else
           if (.not. allocated(heln))   allocate(heln(nrspf,nrspf,5))
         endif
         if (.not. allocated(spin_sign_list))      allocate(spin_sign_list(nrorb_spinorbital))
         if (.not. allocated(list_detl))              allocate(list_detl(nel*nrindep))
         if (.not. allocated(list_detl_spinorbital))      allocate(list_detl_spinorbital(nel*nrindep_spinorbital))
       end subroutine

!> Routine to deallocate the arrays.
       subroutine deallocate_arrays()
         deallocate(x, hmat, tmat, veeao, venmat, venmatmo, taux)
         deallocate(venmatmo_spinorbital)
         deallocate(venmat_spinorbital)
         deallocate(phi, phi2, phin, phin2)
         deallocate(phin_spinorbital)
         deallocate(phin2_spinorbital)
         deallocate(phi_alpha,phi_beta)
         deallocate(phi_spinorbital)
         deallocate(phi2_spinorbital)
         deallocate(dens_mo)
         deallocate(dens_mo_spinorbital)
         deallocate(mf1_t_temp,mf1_t_spinorbital_temp)
         deallocate(hel, hel2, hel3, detl)
         deallocate(hel_spinorbital)
         deallocate(hel2_spinorbital)
         deallocate(hel3_spinorbital)
         deallocate(detl_spinorbital)
         deallocate(vcap, ssqmat, r, rsp, kx, hval3)
         deallocate(hval3_spinorbital)
         deallocate(shdl, hvalrho, allow1)
         deallocate(hvalrdm2)
         deallocate(dhdl,dhdl_spinorbital)
         if (flag_spinorbital == 0) then
           deallocate(hmf)
         end if
         if (flag_spinorbital == 1) then
           deallocate(hmf_spinorbital)
         end if
         deallocate(allow1_spinorbital)
         deallocate(shdl_spinorbital)
         deallocate(hvalrho_spinorbital)
         deallocate(hvalrdm2_spinorbital)
         deallocate(wcap)
         deallocate(rho, rhon)
         deallocate(rho_spinorbital)
         deallocate(rho_temp)
         deallocate(rho_spinorbital_temp)
         deallocate(rdm2)
         deallocate(rdm2_spinorbital)
         deallocate(id_matrix)
         deallocate(id_matrix_spinorbital)
         deallocate(heln)
         deallocate(A)
         deallocate(A_spinorbital)
         deallocate(spin_sign_list)

         deallocate(list_detl)
         deallocate(list_detl_spinorbital)


       end subroutine

 end module globalvars
 !> @file
 !> @brief contains the global variables that are set by the program
 !! and not read from input.
