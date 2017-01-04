program polar_matrix_PIbb
  !!
  !!   NAME
  !!     
  !!   SYNOPSIS
  !!     polar_matrix
  !!   FUNCTION
  !!     Program to calculate the polarisation matrix FT[\PI(t-t')] for the harmonic bath,
  !!     described by a dynamical matrix D_bb' in the bath index b and b'
  !!     See notes GLE pages X-XX
  !!
  !!   INPUTS
  !!     (main routine)
  !!    
  !!   OUTPUT
  !!     (main routine)
  !!    
  !!    
  !!    
  !!    
  !!     
  !!   
  !!   NOTES
  !!     Do not run on Tuesdays!
  !!   BUGS
  !!     None known (ha ha)
  !!   SEE ALSO
  !!     
  !!***
  !!  You can use this space for remarks that should not be included
  !!  in the documentation.
  !! /
  !____________________________________________________________________________
  !
  !Program to calculate the polarisation matrix FT[\PI(t-t')] for the harmonic bath,
  !described by a dynamical matrix D_bb' in the bath index b and b'
  !See notes GLE pages X-XX
  !
  ! Date started: 08 Oct 2013
  !
  !____________________________________________________________________________

  !use matrix_operations
  use parameters_and_co

  implicit none

  integer :: i, ii, j, k, kk, ip, jp, n, m, alpha, beta, tNdimE
  integer :: ilat, jlat, indx, meig, extra, mmeig
  integer :: info, lwork, dimwrk
  
  real(kind=8) :: tmp1, tmp2, dw, w, dEe, tmp1b, tmp2b, mini
  real(kind=8) :: X, der1, der2, Etot, nrm
  real(kind=8) :: sumxx, sumxy, sumxz, sumyy, sumyz, sumzz
  real(kind=8) :: conv_factor, conv_factor_old, Delta_conv, gradientAk, gradientauk, cte_steepgrad, percent
  real(kind=8) :: units_conversion

  real(kind=8), allocatable ::  VLLp(:,:), lapackDbb(:,:), eigenval(:)
  real(kind=8), allocatable ::  xnp1(:), xn(:), xnm1(:), vectmp1(:), vectmp2(:)
  real(kind=8), allocatable ::  Arec(:), Brec(:)
  real(kind=8), allocatable ::  work(:)

  real(kind=8), allocatable ::  tmpvec(:), vecL(:,:), Mm(:)
  real(kind=8), allocatable ::  projDOSp(:), projDOSm(:), fnc1PI(:), fnc2PI(:)
  real(kind=8), allocatable ::  derfnc(:), wgrid(:), res(:)
  real(kind=8), allocatable ::  eigfreqk(:), Ak(:), tauk(:), Ak1(:), tauk0(:)
  real(kind=8), allocatable ::  rdiff(:), rdiff_plus(:) 

  real(kind=8), allocatable ::  AkMat(:,:), taukMat(:,:), modCMat(:,:)

  integer     , allocatable ::  histo(:), histo_seg(:), eigfreq_ind(:)

  complex(kind=8) :: zw

  character :: TypeAt(2)
  character ::  strNip(80), strNjp(80), strNipNjp(80)

  logical :: flagprt
  integer :: iip, jjp, nn, inddex, Nbath, NbatReduced, countplus, countminus
  integer     , allocatable ::  typeNbath(:), typeNbathRed(:), NbathRed(:), indexfit(:)
  real(kind=8), allocatable ::  mass(:)
  real(kind=8) :: xx, yy, zz, twoPI
  logical     , allocatable ::  posign(:,:,:)
  integer     , allocatable ::  truesign(:,:,:)
  real(kind=8), allocatable ::  pmsign(:,:), cbk_xyz(:)

  !********** external functions ********************************************
  real(kind=8) :: ddot






  !**************************************************************************
  !  READ input parameters / def matrix,vector size
  !**************************************************************************

	call read_input_parameters

  	open(unit=25,file="dynmat_coeff.dat")
  	read(25,*) Ndim
  	close(25)

  	tNdim = 3 * Ndim
  	dimwrk = tNdim * int(dexp(0.8d0*dlog(1.d0*tNdim)))
  	write(*,*) "Nb lattice sites Ndim = ",Ndim
  	write(*,*) "tot Ndim = ",tNdim
  	write(*,*) " "


  !**************************************************************************
  !  Allocate arrays
  !**************************************************************************

  allocate( tmpvec(3) )

  allocate( VLLp(tNdim,tNdim), lapackDbb(tNdim,tNdim), eigenval(tNdim) )
  allocate( xnp1(tNdim), xn(tNdim), xnm1(tNdim), vectmp1(tNdim), vectmp2(tNdim) )
  allocate( vecL(Ndim,3), Mm(Ndim) )
  allocate( work(dimwrk) )

  allocate( Arec(Nbcoef), Brec(Nbcoef) )
  allocate( projDOSp(nbp), projDOSm(nbp), fnc1PI(nbp), fnc2PI(nbp), derfnc(nbp), histo(nbp), histo_seg(nbp) )
  allocate( wgrid(nbp), res(nbp) )
  allocate( eigfreqk(nbp), Ak(nbp), tauk(nbp), Ak1(nbp), tauk0(nbp), eigfreq_ind(nbp) )
  allocate( rdiff(nbp), rdiff_plus(nbp) )

  allocate( AkMat(tNdim,nbp), taukMat(tNdim,nbp), modCMat(tNdim,nbp) )

  allocate( typeNbath(Ndim), typeNbathRed(Ndim), NbathRed(Ndim), indexfit(tNdim), mass(tNdim) )

  allocate( posign(tNdim,tNdim,nbp), truesign(tNdim,tNdim,nbp), pmsign(tNdim,nbp), cbk_xyz(3*nbp) )

  !**************************************************************************
  !**************************************************************************
  !		START
  !**************************************************************************
  !**************************************************************************

  	write(*,*) " "
  	write(*,*) " "
  
  !**************************************************************************
  !		READ Dynamical Matrix from file
  !**************************************************************************
  	open(unit=25,file="dynmat_coeff.dat")
  	read(25,*) i

  	do alpha = 1, tNdim
     		do beta = 1, tNdim
        		read(25,*) VLLp(alpha,beta)
     		enddo
  	enddo
  	close(25)

  !**************************************************************************
  !		Units conversion factor
  !**************************************************************************
! [23.05.16]
  	if (LAMMPSunits.eq.1) then
!
! multiply Dyn Mat calculated by LAMMPS in units of [g/mol]^{-1} [eV][\AA]^{-2}
! by number ftm2v = 9648.4484d0
! to get eigenvalues w_k^2 in [ps]^{-2}
!
		units_conversion = 9648.4484d0
  	else
		if (LAMMPSunits.eq.2) then
!
! multiply Dyn Mat calculated by LAMMPS in units of [g/mol]^{-1} [ Kcal/mole][\AA]^{-2}
! by number ftm2v = 0.0004184
! to get eigenvalues w_k^2 in [fs]^{-2}
!
			units_conversion = 0.0004184d0
		else
			write(*,*) "LAMMPS units not know!"
			stop
		endif
  	endif

  !**************************************************************************
  !		make sure Dynamical Matrix is symmetric
  !**************************************************************************
  	do i = 1, tNdim
     		do j = 1, i
!
! eigen val srqt(w^2) = w in [cm^-1]
! tmp1 = ( VLLp(i,j) + VLLp(j,i) ) * 0.5d0 * 5.30887255D-12 * 5.30887255D-12
!
!	cte 5.30887255D-12 = [c *100 * 2\pi]^{-1} transform [s^-1] in [cm^-1] 

!
! eigen val srqt(w^2) = w in [eV]
! tmp1 = ( VLLp(i,j) + VLLp(j,i) ) * 0.5d0 * 6.58212d-16 * 6.58212d-16
!
!	cte 6.58212d-16 = [\hbar]^{-1}

! eigen val srqt(w^2) = w in [meV]
! tmp1 = ( VLLp(i,j) + VLLp(j,i) ) * 0.5d0 * 6.58212d-13 * 6.58212d-13
!
!
!			tmp1 = ( VLLp(i,j) + VLLp(j,i) ) * 0.5d0 * 5.30887255D-12 * 5.30887255D-12


! [08.07.14]
! multiply Dyn Mat calculated by LAMMPS in units of [g/mol]^{-1}[eV][\AA]^{-2}
! by number ftm2v = 9648.4484d0
! to get eigenvalues w_k^2 in [ps]^{-2}
!
			tmp1 = ( VLLp(i,j) + VLLp(j,i) ) * 0.5d0 * units_conversion

			VLLp(i,j) = tmp1
			VLLp(j,i) = tmp1
     		enddo 
  	enddo

!  open(unit=25,file="Dbb_matrix.dat")
!  do alpha = 1, tNdim
!     do beta = 1, t Ndim
!        write(25,*) "Dbb(",alpha,",",beta,") = ", VLLp(alpha,beta)
!     enddo
!  enddo
!  close(25)
!  write(*,*) " "
  !**************************************************************************


  !**************************************************************************
  !		READ Index for REDUCED BATH ATOMS
  !		do the fit for PIbb only on these bath atoms
  !**************************************************************************
  	open(unit=25,file="coord_bath_at_sort.dat")
	read(25,*) 
	read(25,*) 
	read(25,*) 
	read(25,*) Nbath
	read(25,*)
	read(25,*)
	read(25,*)
	read(25,*)
	read(25,*)

  	if (Nbath.ne.Ndim) then
		write(*,*) "WARNING: Nb of atoms in bath region NOT equal to total nb of atoms"
		write(*,*) "WARNING: Please check that Ndim-Nbath=",Ndim-Nbath
		write(*,*) "WARNING: is equal to the nb of atoms in the system region"
		write(*,*) "WARNING: Otherwise the mapping of PIbb' IS MEANINGLESS"
		write(*,*) " "
  	endif

  	do i = 1, Nbath
		read(25,*) typeNbath(i)
  	enddo
  	close(25)

  	open(unit=25,file="coord_bathreduced_at_sort.dat")
	read(25,*) 
	read(25,*) 
	read(25,*) 
	read(25,*) NbatReduced
	read(25,*)
	read(25,*)
	read(25,*)
	read(25,*)
	read(25,*)

  	if (NbatReduced.gt.Nbath) then
		write(*,*) "Incompatible sets for baths atoms and reduced bath atoms"
		write(*,*) "Nb atoms in bath_reduced region LARGER than in bath region"
		stop
  	endif

  	do i = 1, NbatReduced
		read(25,*) typeNbathRed(i), xx, yy, zz, mass(i)
  	enddo
  	close(25)


  !
  !		Map index of reduced bath atoms versus full bath atoms
  !		and Dynamical Matrix
  !
	write(*,*) "Nb atom in         bath = ",Ndim
	write(*,*) "Nb atom in reduced bath = ",NbatReduced
	write(*,*) " "
	write(*,*) "Mapping between bath atoms and reduced bath atoms from LAMMPS"
	write(*,*) " "
	do i = 1, Ndim
		flagprt = .FALSE.
		do j = 1, NbatReduced
			if ( typeNbathRed(j).eq.typeNbath(i) ) then
				NbathRed(j) = i
				write(*,*) "Bath atom i=", i, " , Lammps ID=", typeNbath(i), " , Reduced bath ip=", NbathRed(j)
				flagprt = .TRUE.
			endif
		enddo
                if (.not.flagprt) write(*,*) "Bath atom i=", i, " , Lammps ID=", typeNbath(i)
	enddo
	write(*,*) " "

  !
  !		Calculate proper index for DynMat / PIbb elements
  !		from reduced bath atoms indexes
  !
  	inddex = 0
  	do i = 1, NbatReduced
		do j = 1, 3
			inddex = inddex + 1
			ip = 3*( NbathRed(i) - 1 ) + j
			indexfit(inddex) = ip
		enddo
  	enddo

!	do i = 1, NbatReduced
!		write(*,*) "i=",i," read reduced bath at: ",NbathRed(i)
!	enddo

	write(*,*) "        Bath atoms: ", Ndim," * 3 = ", tNdim, "indexes" 
	write(*,*) "Reduced Bath atoms: ", NbatReduced," * 3 = ", 3*NbatReduced, "indexes" 
	do i = 1, 3*NbatReduced
		write(*,*) "ip=",i," calc reduced bath at index: ",indexfit(i)
	enddo
	write(*,*) " "
	
  !
  !		Array indexfit[1->NbatReduced*3]
  !		contains the correct indexes for fit of PIbb
  !		versus indexes in DynMat from LAMMPS	
  !**************************************************************************






  !**************************************************************************
  !		START Solving Eigenvalues/vectors problem with LAPACK
  !**************************************************************************
  lapackDbb = VLLp
  call dsyev("V","U",tNdim,lapackDbb,tNdim,eigenval,work,-1,info)
  lwork=work(1)
  write(*,*) "LWORK for DSYEV is          :", lwork
  write(*,*) "Declared size for work() is :", dimwrk
  if (lwork.gt.dimwrk) stop

  lapackDbb = VLLp
  call dsyev("V","U",tNdim,lapackDbb,tNdim,eigenval,work,lwork,info)
  write(*,*) "info call dsyev is: ", info
  write(*,*) " "
  write(*,*) " "

  if (info.ne.0) stop

  do n = 1, tNdim
     write(*,*) "Eigenval(",n,") = ",eigenval(n)  
  enddo
  write(*,*) " "

  !**************************************************************************
  !		Print local and off-diag DOS-like fnc
  !**************************************************************************
!  tmp1 = eigenval(2) - eigenval(1)
!  do n = 2, Ndim-1
!     dE = eigenval(n+1) - eigenval(n)
!     if ( (dE.le.tmp1) .and. (dE.gt.1.d-07) ) tmp1 = dE
!  enddo
!  write(*,*) "Minimum dE between eigenval : ", tmp1
!
!  dE = tmp1 / 2.d0

  	tmp1 = 0.d0
  	i = 0
  	do n = 1, tNdim-1
     		dEe = eigenval(n+1) - eigenval(n)
     		tmp1 = tmp1 + dEe
     		if (dEe.gt.1.d-07) i = i + 1
  	enddo
  	dEe = tmp1 / (2.d0*i)

  	if (dE.lt.1.d-8) dE = dEe



  !**************************************************************************
  !		print out eigenval(n)
  !**************************************************************************
  	open(unit=25,file="Eigenval.dat")
  	open(unit=26,file="posEigenfreq.dat")
  	do n = 1, tNdim
     		write(25,*) eigenval(n), 0.0
     		if (eigenval(n).gt.0.d0) write(26,*) dsqrt(eigenval(n)), 0.0
  	enddo
  	close(25)
  	close(26)
  
  	dw = eigenval(tNdim)*1.02d0/(nbp-1)
  
  	write(*,*) "Estimate for Lorentzian broadening : ", dEe
  	write(*,*) "Value    for Lorentzian broadening : ", dE
  	write(*,*) "Value    for dw grid : ", dsqrt(dw)
  	write(*,*) " "
  	write(*,*) "total Etot = ", Etot
  	write(*,*) " "


  !**************************************************************************
  !		DO loop for calc local DIAG PIbb matrix elements
  !		extract peak position w_k
  !**************************************************************************

  	histo = 0

  	do iip = 1, 3*NbatReduced

		ip = indexfit(iip)

		write(*,*) "DIAG Site/coord (ilat/xyz) indices ip = ", ip

  !**************************************************************************
  !		check sum rule on dyn mat
  !**************************************************************************

  !**************************************************************************
  !		print out local DOS
  !             with lorentzians around eigenval(n) & width 1/2 min Delta E
  !**************************************************************************

  !**************************************************************************
  !		LOOP on starting vec for Lanczos tri-diagonalisation
  !**************************************************************************

		open(unit=25,file="imDbbrecur_vs_eig.dat")
		open(unit=26,file="imDbbrecur_vs_w.dat")
		open(unit=27,file="imDbbrecur_ovw_vs_w.dat")
!		open(unit=22,file="abcoef_diag.dat")

		write(*,*) "Diagonal - Lanczos recursion"

		n = ( ip - 1 ) /3
		write(*,*) "Initial vector - index ", ip ," / ilat = ", n+1 , " / (xyz)=(123): ",ip-3*n 
     		write(*,*) " "
     
		xn = 0
        	xn(ip) = 1.d0
 
		call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)

!		do n = 1, Nbcoef
!        		write(*,*) "Coef a(",n,") = ", Arec(n), " #    b(",n,") = ", Brec(n)  
!			write(22,*) n, Arec(n), Brec(n)
!     		enddo
!     		write(*,*) " "

     		do j = 1, nbp-1
        		w = dw * j
                	tmp2 = sqrt(w)
        		call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!        		call cont_frac2(Nbcoef,Arec,Brec,eigenval(tNdim),w,Typecontfrac,dE,zw) 

!
!	see Notes NoteBook page 109 / 19mar2014
!
			fnc1PI(j) = - 2.d0 * dimag(zw) / tmp2
!
!
!			fnc1PI(j) = - dimag(zw) / pi 
!        		write(25,*) w, - dimag(zw) / pi
!        		write(26,*) tmp2, - dimag(zw) / pi
!        		write(27,*) tmp2, fnc1PI(j)
			tmp1 = tmp1 + fnc1PI(j)
     		enddo
!		close(22)

  !**************************************************************************
  !		START fitting procude of the PI matrix elements : find peak positions
  !**************************************************************************
     		do j = 2, nbp-2
        		derfnc(j) = ( fnc1PI(j+1) - fnc1PI(j-1) ) / 2.d0
			derfnc(j) = derfnc(j) * (nbp - 2) / tmp1 
     		enddo
		do j = nbp-2, 3, -1
			tmp2 = derfnc(j) * derfnc(j-1)
			if (tmp2.lt.0.d0) then
				if ( derfnc(j) .lt. 0.d0 ) then
					if ( dabs( derfnc(j) ) .gt. dabs( derfnc(j-1) ) ) then
						histo(j-1) = histo(j-1) + 1
					else
						histo(j) = histo(j) + 1
					endif
				endif
			endif
		enddo

  !**************************************************************************
  !		END fitting procude of the PI matrix elements : find peak positions
  !**************************************************************************
  !		END loop for calc diag elements PIbb
  !**************************************************************************
  	enddo

	close(25)
	close(26)
	close(27)

  !**************************************************************************
  !		Fitting procude of the PI matrix elements : treat all w_k peaks
  !**************************************************************************
	open(unit=25,file="histoEigenFreq.dat")
  	do j = 2, nbp-2
		w = dw * j
                tmp2 = dsqrt(w)
        	write(25,*) tmp2, histo(j)
	enddo
	close(25)

        n = (nbp - 2) / segment
	histo_seg = 0

	write(*,*) "Loop segmentation / integration"
	write(*,*) "nbp = ",nbp
	do i = 1, n
		write(*,*) "i= ",i
		alpha = 0
		beta = 0
		do j = 1, segment
			k = (i-1)*segment
			write(*,*) "k+j= ",k+j
			alpha = alpha + histo(k+j)
			beta = beta + (k+j) * histo(k+j)
		enddo
		if (alpha.gt.0) then
			beta = beta / alpha
			histo_seg(beta) = alpha
		endif
	enddo

	open(unit=25,file="histosegEigenFreq.dat")
	meig = 0
  	do j = 2, nbp-2
		w = dw * j
                tmp2 = sqrt(w)
        	write(25,*) tmp2, histo_seg(j)
                if ( histo_seg(j) .gt. 0 ) meig = meig + 1 
	enddo
	close(25)

	write(*,*) "Nb of calculated peaks  = ", meig
  	n = 1
	do j = 2, nbp-2
		w = dw * j
		if ( histo_seg(j) .gt. 0 ) then
			eigfreqk(n) =  dsqrt(w)
			eigfreq_ind(n) =  j
			n = n + 1
		endif
	enddo

	open(unit=25,file="calcEigenFreq.dat")
  	do n = 1, meig
        	write(25,*) eigfreqk(n), 0.0
	enddo
	close(25)
  !**************************************************************************
  !		Fitting procude of the PI matrix elements : ADD extra w_k peaks
  !		MAKE A SUB ROUTINE of this part
  !		not only linear space, but quadratic around w = 0, or other
  !**************************************************************************

  	if (type_extra.gt.0) then
		call extra_peaks(meig,eigfreqk,type_extra,extra_coef,extra)
	else
		extra = 0
	endif

	meig = meig + extra

	open(unit=25,file="chosen_wk.dat")
  	do n = 1, meig
        	write(25,*) eigfreqk(n), 0.0
	enddo
	close(25)

	write(*,*) "Nb of chosen w_k peaks = ", meig

	do n = 1, meig
!		tauk0(n) = 20.d0/dsqrt(dE)
		tauk0(n) = scaling_tauks / dsqrt(dE)
	enddo
	
	
  !**************************************************************************
  !		START least-square fitting of PI matrix elements : height & width
  !		diagonal elements 
  !**************************************************************************

	do iip = 1, 3*NbatReduced

		ip = indexfit(iip)

  !**************************************************************************
  !		Recalculate PI matrix elements from Lanczos 
  !**************************************************************************
		n = ( ip - 1 ) /3
		write(*,*) "Initial vector - index ", ip ," / ilat = ", n+1 , " / (xyz)=(123): ",ip-3*n 
     		write(*,*) " "
     
		call intstr(ip,strNip,4)

		xn = 0
        	xn(ip) = 1.d0

		call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)
!		call write_real_array1(Nbcoef,histo,Arec,"recur_coef_A",strNip,4)
!		call write_real_array1(Nbcoef,histo,Brec,"recur_coef_A",strNip,4)

     		do j = 1, nbp-1
        		w = dw * j
                	wgrid(j) = dsqrt(w)
        		call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!
!	see Notes NoteBook page 109 / 19mar2014
!
			fnc1PI(j) = - 2.d0 * dimag(zw) / wgrid(j)
!
!
    		enddo

		call write_real_array2(nbp-1,wgrid,fnc1PI,"imDbbrecur_ovw",strNip,4)

		
  !**************************************************************************
  !		Choose initial peak amplitude & width for fit 
  !**************************************************************************
		mmeig = meig - extra

		do n = 1, mmeig
			Ak1(n) = scaling_amplitude * fnc1PI(eigfreq_ind(n))
!			Ak1(n) = 0.5d0 * fnc1PI(eigfreq_ind(n))
!			Ak1(n) = 0.2d0 * fnc1PI(eigfreq_ind(n)) / tauk0(1)
			histo(n) = n
		enddo

        	do n = mmeig+1, meig-1
			w = eigfreqk(n)
			histo(n) = n
	
			ii = 1
			do j = 2, nbp-1
				if ( (w.gt.wgrid(j-1)) .and. (w.lt.wgrid(j)) ) then
					ii = j
				endif
			enddo
			Ak1(n) = scaling_amplitude * fnc1PI(ii)
!			Ak1(n) = 0.5d0 * fnc1PI(ii)
!                	Ak1(n) = 0.2d0 * fnc1PI(ii) / tauk0(1)

!               write(*,*) "new peak w_k(", n+meig,")=", eigfreqk(n+meig)
!		write(*,*) "between ",wgrid(ii-1)," and ",wgrid(ii)
!		write(*,*) "Extra Ak(", n+meig,") =", Ak(n+meig), " for ii=", ii
!		write(*,*) " "

        	enddo
		histo(meig)  = meig

!	stop

  !**************************************************************************
  !		Calculate initial fnc_value and residual		/ ip index
  !**************************************************************************
		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk0,res)
		call write_real_array2(nbp-1,wgrid,res,"initfit_imDbbrecur_ovw",strNip,4)

		rdiff = fnc1PI - res
		conv_factor = dsqrt( ddot(nbp,rdiff,1,rdiff,1) )
		nrm = dsqrt( ddot(nbp,fnc1PI,1,fnc1PI,1) )
        	conv_factor = conv_factor / nrm 

		write(*,*) "FOR FIT: initial convergence factor = ", conv_factor

!	do n = 1, meig
!	write(*,*) "FOR FIT: initial Ak(",n,")   = ", Ak(n), "tauk(",n,") = ", tauk(n)
!	enddo

!	open(unit=22,file="initfit_coef_Ak.dat")
!	open(unit=23,file="initfit_coef_tauk.dat")
!	do n = 1, meig
!        	write(22,*) n, Ak(n)
!        	write(23,*) n, tauk(n)
!	enddo
!	close(22)
!	close(23)
	
  !**************************************************************************
  !		DO the fit 					/ ip index
  !**************************************************************************
		tauk = tauk0
		call minimization1(nbp,wgrid,fnc1PI,conv_factor,nrm,meig,eigfreqk,Ak1,tauk)
		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk,res)

  !**************************************************************************
  !		Print out results 				/ ip index
  !**************************************************************************
		call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw",strNip,4)
!		call write_real_array1(meig,histo,Ak1,"fit_coef_Ak",strNip,4)
!		call write_real_array1(meig,histo,tauk1,"fit_coef_tauk",strNip,4)

		do n = 1, meig
			AkMat(ip,n) = Ak1(n)
			taukMat(ip,n) = tauk(n)
			modCMat(ip,n) = Ak1(n) / tauk(n)
			write(*,*) "index ip=",ip," tauk(n=",n,")= ",tauk(n)
		enddo
  !**************************************************************************
  !		END LOOP on diagonal elements 
  !**************************************************************************
	enddo
 
  !**************************************************************************
  !		POST TREATMENT tau_k should be independent of index ip,jp
  !**************************************************************************
	do n = 1, meig
		mini = taukMat(indexfit(1),n)
		tauk(n) = mini

		do iip = 2, 3*NbatReduced
			i = indexfit(iip)
			tmp1 = taukMat(i,n)
			if ( tmp1.lt.mini ) then
				mini = tmp1
				tauk(n) = mini
			endif
		enddo

		write(*,*) "min tauk(n=",n,")= ",tauk(n)
	enddo

	call intstr(0,strNip,1)
!	call write_real_array1(meig,histo,tauk,"fit_coef_tauk",strNip,1)

  !**************************************************************************
  !		POST TREATMENT for Amplitudes
  !**************************************************************************
	write(*,*) "Post treatment: averaged tauk"
	write(*,*) "tauk(n)=",tauk(1:meig)

	do iip = 1, 3*NbatReduced
		i = indexfit(iip)
		do n = 1, meig
			AkMat(i,n) = modCMat(i,n) * tauk(n)
			write(*,*) "modCMat(ip=", i, ",n=", n, ")=", modCMat(i,n)
			Ak1(n) = AkMat(i,n) 
		enddo
		
!		write(*,*) "ip= ",i
!		write(*,*) "Ak(n)=",Ak1(1:meig)
		write(*,*)

		call intstr(i,strNip,4)
		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk,res)
		call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_modCavtauk",strNip,4)

	enddo



  !**************************************************************************
  !		END least-square fitting of chose diag PI matrix elements 
  !**************************************************************************

!	write(*,*) "READY for NEXT STEP?"

  !**************************************************************************
  !		DO loops for calc off-diga PIbb' matrix elements
  !**************************************************************************

	do iip = 1, 3*NbatReduced

		ip = indexfit(iip)

		call intstr(ip,strNip,4)
		strNipNjp(1:4) = strNip(1:4)

		do jjp = iip+1, 3*NbatReduced

			jp = indexfit(jjp)

			call intstr(jp,strNjp,4)
			strNipNjp(5:8) = strNjp(1:4)

  !**************************************************************************
  !		Get target function from Lanczos tri-diag
  !**************************************************************************
			write(*,*) "Off Diagonal - Lanczos recursion"
			n = ( ip - 1 ) /3
			write(*,*) "Initial vector - index ", ip ," / ilat = ", n+1 , " / (xyz)=(123): ",ip-3*n 
			n = ( jp - 1 ) /3
			write(*,*) "Initial vector - index ", jp ," / ilat = ", n+1 , " / (xyz)=(123): ",jp-3*n 
     			write(*,*) " "
     
			xn = 0
        		xn(ip) = 1.d0 / sqrt(2.0)
        		xn(jp) = 1.d0 / sqrt(2.0)
 
			call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)

			write(*,*) "Symmetric superposition ip + jp"
			do j = 1, nbp-1
        			w = dw * j
        			call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!
!	see Notes NoteBook page 109 / 19mar2014
!
				projDOSp(j) = - 2.d0 * dimag(zw) / wgrid(j)
!
     			enddo

			xn = 0
        		xn(ip) =  1.d0 / sqrt(2.0)
        		xn(jp) = -1.d0 / sqrt(2.0)
 
			call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)

			write(*,*) "Antisymmetric superposition ip - jp"
	     		do j = 1, nbp-1
        			w = dw * j
        			call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!
!	see Notes NoteBook page 109 / 19mar2014
!
				projDOSm(j) = - 2.d0 * dimag(zw) / wgrid(j)
!
     			enddo

     			do j = 1, nbp-1
				fnc1PI(j) = ( projDOSp(j) - projDOSm(j) ) / 2.d0
     			enddo
			nrm = dsqrt( ddot(nbp,fnc1PI,1,fnc1PI,1) ) 	
!			call write_real_array2(nbp-1,wgrid,fnc1PI,"imDbbrecur_ovw_IJ",strNipNjp,8)

  !**************************************************************************
  !		Get crossed amplitudes
  !**************************************************************************
			do n = 1, meig
				Ak(n) = dsqrt( modCMat(ip,n) * modCMat(jp,n) ) * tauk(n)
			enddo

  !**************************************************************************
  !		Get Initial sign for the amplitudes
  !**************************************************************************
			do n = 1, meig-extra
				if ( fnc1PI(eigfreq_ind(n)) .lt. 0 ) then
					Ak(n) = - Ak(n)
!				write(*,*) "Negative amplitude: Ak(",n,") =",Ak(n)
				endif
			enddo

      			do n = meig-extra+1, meig 
				w = eigfreqk(n)
				ii = 2
				do j = 2, nbp-1
					if ( w.gt.wgrid(j-1) )  then
						if ( w.lt.wgrid(j) )  ii = j
					endif
				enddo
				if ( fnc1PI(ii-1) .lt. 0.d0 ) then
					Ak(n) = - Ak(n)
!				write(*,*) "Negative amplitude: Ak =",Ak(n)
				endif
!			write(*,*) " "
        		enddo

			call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
!			call write_real_array2(nbp-1,wgrid,res,"initfit_imDbbrecur_ovw_IJ_avtauk",strNipNjp,8)

			conv_factor = 100.d0 
			call minimization_IJ_2(nbp,wgrid,fnc1PI,conv_factor,nrm,meig,eigfreqk,Ak,tauk)
			call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
!			call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_IJ_avtauk_min2",strNipNjp,8)
!			call write_real_array1(meig,histo,Ak,"fit_coef_Ak_after_min2",strNipNjp,8)

  !**************************************************************************
  !		POST TREATMENT for off-diag elements
  !		Get the proper sign of sqrt ( ModCMat(ip/jp,n) )
  !**************************************************************************
			do n = 1, meig
				if (Ak(n).ge.0.d0) then
					posign(ip,jp,n)=.TRUE.
				else
					posign(ip,jp,n)=.FALSE.
				endif
			enddo
			posign(jp,ip,n) = posign(ip,jp,n)
			
  !**************************************************************************
  !		END loops for calc PIbb'
  !**************************************************************************
		enddo		! LOOP on jp
	enddo			! LOOP on ip


  !**************************************************************************
  !		POST TREATMENT for off-diag elements ip <> jp
  !**************************************************************************
  !		FIND THE CORRECT SIGN OF EACH ELEMENT Cbk
  !		from the knowledge of the sign of the cross product Cbk * Cb'k
  !**************************************************************************

	open(unit=26,file="SignMixedCoef.dat")
	do n = 1, meig
		write(26,*) "Peak No: n=",n
		write(26,*) " "
		do iip = 1, 3*NbatReduced
			ip = indexfit(iip)
			do jjp = iip+1, 3*NbatReduced
				jp = indexfit(jjp)
				write(26,*) "Coef(ip=",ip,",jp=",jp,") == ",posign(ip,jp,n)
			enddo
		write(26,*) " "
		enddo
		write(26,*) " "
	enddo
	close(26)

	open(unit=26,file="SignEachCoef.dat")
	do n = 1, meig

		write(26,*) "Peak No: n=",n
		write(26,*) " "

		iip = 1
		ip = indexfit(iip)
		write(26,*) "Sign Coef(ip=",ip,",jp) == ", posign(ip,indexfit(1:3*NbatReduced),n)
		write(26,*) " "

		truesign(ip,ip,n) = +1 		
		pmsign(ip,n) = +1.d0
		do ii = 2, 3*NbatReduced
			if ( posign(ip,indexfit(ii),n) ) then
				truesign(ip,indexfit(ii),n) = +1
			else
				truesign(ip,indexfit(ii),n) = -1
			endif
		enddo
		
		do jjp = 2, 3*NbatReduced
			jp = indexfit(jjp)

			if ( posign(ip,jp,n) ) then
!			    	write(26,*) "Sign Coef(ip=",jp,",jp) == ", posign(jp,indexfit(jjp+1:3*NbatReduced),n)
			    	write(26,*) "Sign Coef(ip=",jp,",jp) == ", posign(jp,indexfit(1:3*NbatReduced),n)
				do ii = 1, 3*NbatReduced
					if ( posign(jp,indexfit(ii),n) ) then
						truesign(jp,indexfit(ii),n) = +1
					else
						truesign(jp,indexfit(ii),n) = -1
					endif
				enddo
			else
			    	write(26,*) "Sign Coef(ip=",jp,",jp) == ", .not.posign(jp,indexfit(1:3*NbatReduced),n)
				do ii = 1, 3*NbatReduced
					if ( posign(jp,indexfit(ii),n) ) then
						truesign(jp,indexfit(ii),n) = -1
					else
						truesign(jp,indexfit(ii),n) = +1
					endif
				enddo
			endif

!			if ( posign(1,ip,n) ) then
!				do jjp = 1, iip
!					jp = indexfit(jjp)
!			    		write(26,*) "Sign Coef(ip=",ip,",jp=)", jp," == ", posign(ip,jp,n)
!				enddo
!			else
!				do jjp = 1, iip
!					jp = indexfit(jjp)
!			    		write(26,*) "Sign Coef(ip=",ip,",jp=)", jp," == ", .not.posign(ip,jp,n)
!				enddo
!			endif

		write(26,*) " "
		enddo
		write(26,*) " "
	enddo

	do n = 1, meig

		write(26,*) "Peak No: n=",n
		write(26,*) " "

		do jjp = 2, 3*NbatReduced
			jp = indexfit(jjp)
			countplus = 0
			countminus = 0

			do iip = 1, jjp-1
				ip = indexfit(iip)

				if ( truesign(ip,jp,n).eq.1 ) then
					countplus  = countplus  + 1
				else
					countminus = countminus + 1
				endif
			enddo
			write(26,*) "Sign coef jp = ", jp, "   (+):", countplus, " (-):",countminus, "   sign=",countplus-countminus

			if ( (countplus-countminus).gt.0) pmsign(jp,n) = +1.d0
			if ( (countplus-countminus).lt.0) pmsign(jp,n) = -1.d0
			if ( (countplus-countminus).eq.0) then
				if ( posign(indexfit(1),jp,n) ) then
					pmsign(jp,n) = +1.d0
				else
					pmsign(jp,n) = -1.d0
				endif
			endif
			
		enddo
		write(26,*) " "
	enddo
	close(26)

	open(unit=27,file="SignPMONE_EachCoef.dat")
	do n = 1, meig

		write(27,*) "Peak No: n=",n
		write(27,*) " "

		do jjp = 1, 3*NbatReduced
			jp = indexfit(jjp)
			write(27,*) "Sign Coef(ip=",jp,",jp) == ", truesign(jp,indexfit(1:3*NbatReduced),n)
			write(27,*) " "
		enddo
		write(27,*) " "
	enddo
	close(27)



  !**************************************************************************
  !		RECALCULATED OFF-DIAG ELEMENS PIbb'
  !		with SIGN attribution for each Cbk element
  !**************************************************************************
	do iip = 1, 3*NbatReduced

		ip = indexfit(iip)

		call intstr(ip,strNip,4)
		strNipNjp(1:4) = strNip(1:4)

		do jjp = iip+1, 3*NbatReduced

			jp = indexfit(jjp)

			call intstr(jp,strNjp,4)
			strNipNjp(5:8) = strNjp(1:4)

  !**************************************************************************
  !		Get new crossed amplitudes
  !**************************************************************************
			do n = 1, meig
				Ak(n) = pmsign(ip,n) * dsqrt( modCMat(ip,n) * modCMat(jp,n) ) * pmsign(jp,n) * tauk(n)
			enddo

			call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
!			call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_IJ_avtauk_withsign",strNipNjp,8)
		enddo
	enddo




  !**************************************************************************
  !		GENERATE OUTPUT FILE FOR LAMMPS GLE DYNAMICS, need:
  !		ONE FILE FOR eigfreqk(n) AND tauk(n)
  !     	and
  !		ONE FILE FOR cbk = pmsign(ip,n) * dsqrt( modCMat(ip,n) )
  !		
  !		according to Chris format for LAMMPS GLE code
  !		[21.05.14]
  !**************************************************************************
!
!	then print out resutls w_k, t_k and pmsign(ip,n) * dsqrt( modCMat(ip,n) ) according to Chris format for Lammps code
!
!
  !**************************************************************************
  !		Include the 2PI factors in w_k and cbk coef
  !**************************************************************************
	twoPI = 8.d0 * datan(1.d0)
	write(*,*) "2PI =", twoPI
	write(*,*) " "
  !**************************************************************************

	open(unit=27,file="fit_coef_taukomegak_lmp.dat")
!	write(27,*) NbatReduced
	do n = 1, meig
		write(27,*) (n-1), tauk(n), twoPI * eigfreqk(n)
	enddo
	close(27)

	open(unit=27,file="fit_coef_cbk_lmp.dat")
!	write(27,*) NbatReduced
	write(27,*) (meig+1)*3
	do i = 1, NbatReduced
!
!		FIND the CORRECT IP NUMBER !
!
		ip = 3*( NbathRed(i) - 1 )

		nn = 1
		do n = 1, meig
			cbk_xyz(nn+0) =  pmsign(ip+1,n) * dsqrt( modCMat(ip+1,n) )
			cbk_xyz(nn+1) =  pmsign(ip+2,n) * dsqrt( modCMat(ip+2,n) )
			cbk_xyz(nn+2) =  pmsign(ip+3,n) * dsqrt( modCMat(ip+3,n) )
			nn = nn + 3
		enddo 
		write(27,*) (i-1), typeNbathRed(i), mass(i), cbk_xyz(1:3*meig) / twoPI
	enddo
	close(27)

111	format(i5," ",i5,"    ",f6.3," ",5000(f18.15," ") )



  !*************************************************************************
  !     DEALLOCATE ARRAYS

  deallocate( VLLp, xnp1, xn, xnm1, vectmp1, vectmp2, Mm )
  deallocate( Arec, Brec, lapackDbb, work, eigenval )
  deallocate( projDOSp, projDOSm )

  deallocate( fnc1PI, fnc2PI, derfnc, histo, histo_seg, wgrid, res )
  deallocate( eigfreqk, Ak, tauk, Ak1, tauk0, eigfreq_ind, rdiff, rdiff_plus )

  deallocate( AkMat, taukMat, modCMat )

  deallocate( NbathRed, indexfit, posign, truesign)



end program polar_matrix_PIbb




