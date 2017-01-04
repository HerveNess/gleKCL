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

  integer :: i, ii, j, k, kk, ip, jp, n, m, alpha, beta
  integer :: ilat, jlat, indx, meig, extra
  integer :: info, lwork, dimwrk
  
  real(kind=8) :: tmp1, tmp2, dw, w, dEe, tmp1b, tmp2b
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
  real(kind=8), allocatable ::  eigfreqk(:), Ak(:), tauk(:), Ak1(:), tauk1(:), Ak2(:), tauk2(:)
  real(kind=8), allocatable ::  modC1(:), modC2(:) 
  real(kind=8), allocatable ::  rdiff(:), rdiff_plus(:) 

  integer     , allocatable ::  histo(:), histo_seg(:), eigfreq_ind(:)

  complex(kind=8) :: zw

  character :: TypeAt(2)
  character ::  strNip(80), strNjp(80), strNipNjp(80)

  !********** external functions ********************************************
  real(kind=8) :: ddot






  !  READ input parameters / def matrix,vector size

  call read_input_parameters


  open(unit=25,file="dynmat_coeff.dat")
  read(25,*) Ndim
  close(25)

  tNdim = 3 * Ndim
  dimwrk = tNdim * int(dexp(0.8d0*dlog(1.d0*tNdim)))
  write(*,*) "Nb lattice sites Ndim = ",Ndim
  write(*,*) "tot Ndim = ",tNdim
  write(*,*) " "


  !  Allocate arrays

  allocate( tmpvec(3) )

  allocate( VLLp(tNdim,tNdim), lapackDbb(tNdim,tNdim), eigenval(tNdim) )
  allocate( xnp1(tNdim), xn(tNdim), xnm1(tNdim), vectmp1(tNdim), vectmp2(tNdim) )
  allocate( vecL(Ndim,3), Mm(Ndim) )
  allocate( work(dimwrk) )

  allocate( Arec(Nbcoef), Brec(Nbcoef) )
  allocate( projDOSp(nbp), projDOSm(nbp), fnc1PI(nbp), fnc2PI(nbp), derfnc(nbp), histo(nbp), histo_seg(nbp) )
  allocate( wgrid(nbp), res(nbp) )
  allocate( eigfreqk(nbp), Ak(nbp), tauk(nbp), Ak1(nbp), tauk1(nbp), Ak2(nbp), tauk2(nbp), eigfreq_ind(nbp) )
  allocate( rdiff(nbp), rdiff_plus(nbp), modC1(nbp), modC2(nbp) )

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
  !		DO loop for calc and print out local DIAG DOS_bb / PIbb
  !             matrix elements
  !**************************************************************************

  histo = 0

  do ip = 1, tNdim

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

!	open(unit=25,file="imDbbrecur_vs_eig.dat")
!	open(unit=26,file="imDbbrecur_vs_w.dat")
!	open(unit=27,file="imDbbrecur_ovw_vs_w.dat")

     	
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
			fnc1PI(j) = - dimag(zw) / ( pi * tmp2 )
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

  !**************************************************************************
  !		END loop for calc diag elements PIbb
  !**************************************************************************
  enddo

!	close(25)
!	close(26)
!	close(27)

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
!			tauk1(n) = 20.d0/dsqrt(dE)
!			tauk2(n) = 20.d0/dsqrt(dE)
			tauk1(n) = scaling_tauks / dsqrt(dE)
			tauk2(n) = scaling_tauks / dsqrt(dE)
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
		write(*,*) "Calculating extra peaks!"
		call extra_peaks(meig,eigfreqk,type_extra,extra_coef,extra)
	else
		extra = 0
	endif

	open(unit=25,file="chosen_wk.dat")
  	do n = 1, meig+extra
        	write(25,*) eigfreqk(n), 0.0
	enddo
	close(25)
	write(*,*) "Nb of chosen w_k peaks = ", meig+extra
	
	
  !**************************************************************************
  !		START least-square fitting of PI matrix elements : height & width 
  !**************************************************************************
  ip = 666
  jp = 666

!  do while ( (ip.ne.0) )

	write(*,*) "FOR FIT: Site/coord (ilat/xyz) indices ip, jp = ?"
	read(5,*) ip, jp
	
  !**************************************************************************
  !		Recalculate PI matrix elements from Lanczos 
  !**************************************************************************
	n = ( ip - 1 ) /3
	write(*,*) "Initial vector - index ", ip ," / ilat = ", n+1 , " / (xyz)=(123): ",ip-3*n 
     	write(*,*) " "
     
	xn = 0
        xn(ip) = 1.d0
	call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)

     	do j = 1, nbp-1
        	w = dw * j
                wgrid(j) = dsqrt(w)
        	call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
		fnc1PI(j) = - dimag(zw) / ( pi * wgrid(j) )
!		fnc1PI(j) = - dimag(zw) / pi
     	enddo

	call intstr(ip,strNip,4)
	call write_real_array2(nbp-1,wgrid,fnc1PI,"imDbbrecur_ovw",strNip,4)

	
	n = ( jp - 1 ) /3
	write(*,*) "Initial vector - index ", jp ," / jlat = ", n+1 , " / (xyz)=(123): ",jp-3*n 
     	write(*,*) " "
     
	xn = 0
        xn(jp) = 1.d0
	call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)
     	do j = 1, nbp-1
        	w = dw * j
        	call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
		fnc2PI(j) = - dimag(zw) / ( pi * wgrid(j) )
!		fnc2PI(j) = - dimag(zw) / pi
     	enddo

	call intstr(jp,strNjp,4)
	call write_real_array2(nbp-1,wgrid,fnc2PI,"imDbbrecur_ovw",strNjp,4)

	
  !**************************************************************************
  !		print out local DOS
  !             with lorentzians around eigenval(n) & width 1/2 min Delta E
  !**************************************************************************
	open(unit=27,file="DOSoverw_vs_w.dat")
     	do i = 1, nbp-1
        	w = dw * i
        	tmp1  = 0.d0
        	tmp1b = 0.d0
        	do n = 1, tNdim
           		tmp2 = dE / ( ( w - eigenval(n) )**2 + dE**2 )
           		tmp2  = tmp2 * ( lapackDbb(ip,n) * lapackDbb(jp,n) ) / PI
           		tmp1  = tmp1  + tmp2  
        	enddo
		tmp2 = sqrt(w)
	        res(i) = tmp1/tmp2
     	enddo
	call intstr(ip,strNip,4)
	call write_real_array2(nbp-1,wgrid,res,"DOSoverw_vs_w",strNip,4)


  !**************************************************************************
  !		Choose initial peak amplitude & width for fit 
  !**************************************************************************
	do n = 1, meig
!		Ak1(n) = 0.2d0 * fnc1PI(eigfreq_ind(n)) / tauk1(n)
!		Ak2(n) = 0.2d0 * fnc2PI(eigfreq_ind(n)) / tauk2(n)
		Ak1(n) = scaling_amplitude * fnc1PI(eigfreq_ind(n))
		Ak2(n) = scaling_amplitude * fnc2PI(eigfreq_ind(n))
		histo(n) = n
	enddo

        do n = 1, extra-1
		tauk1(n+meig) = tauk1(1)
		tauk2(n+meig) = tauk2(1)

		w = eigfreqk(n+meig)

		histo(n+meig) = n+meig

		ii = 1
		do j = 2, nbp-1
			if ( (w.gt.wgrid(j-1)) .and. (w.lt.wgrid(j)) ) then
				ii = j
			endif
		enddo
!                Ak1(n+meig) = 0.2d0 * fnc1PI(ii) / tauk1(1)
!                Ak2(n+meig) = 0.2d0 * fnc2PI(ii) / tauk2(1)
                Ak1(n+meig) = scaling_amplitude * fnc1PI(ii)
                Ak2(n+meig) = scaling_amplitude * fnc2PI(ii)

!                write(*,*) "new peak w_k(", n+meig,")=", eigfreqk(n+meig)
!		write(*,*) "between ",wgrid(ii-1)," and ",wgrid(ii)
!		write(*,*) "Extra Ak(", n+meig,") =", Ak(n+meig), " for ii=", ii
!		write(*,*) " "
        enddo

!	stop


	meig = meig + extra
	tauk1(meig) = tauk1(1)
	tauk2(meig) = tauk2(1)

!	cte_steepgrad = 0.00001d0

  !**************************************************************************
  !		Calculate initial fnc_value and residual		/ ip index
  !**************************************************************************
	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk1,res)

	call write_real_array2(nbp-1,wgrid,res,"initfit_imDbbrecur_ovw",strNip,4)



	rdiff = fnc1PI - res
	conv_factor = dsqrt( ddot(nbp,rdiff,1,rdiff,1) )
	nrm = dsqrt( ddot(nbp,fnc1PI,1,fnc1PI,1) )
        conv_factor = conv_factor / nrm 

!        conv_factor = conv_factor / 1000.0

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
	call minimization1(nbp,wgrid,fnc1PI,conv_factor,nrm,meig,eigfreqk,Ak1,tauk1)

	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk1,res)

  !**************************************************************************
  !		Print out results 				/ ip index
  !**************************************************************************
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw",strNip,4)

	call write_real_array1(meig,histo,Ak1,"fit_coef_Ak",strNip,4)
	call write_real_array1(meig,histo,tauk1,"fit_coef_tauk",strNip,4)



  !**************************************************************************
  !		Calculate initial fnc_value and residual 	/ jp index
  !**************************************************************************
	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak2,tauk2,res)
	call write_real_array2(nbp-1,wgrid,res,"initfit_imDbbrecur_ovw",strNjp,4)

	rdiff = fnc2PI - res
	conv_factor = dsqrt( ddot(nbp,rdiff,1,rdiff,1) )
	nrm = dsqrt( ddot(nbp,fnc2PI,1,fnc2PI,1) )
        conv_factor = conv_factor / nrm 

	write(*,*) "FOR FIT: initial convergence factor = ", conv_factor
	
  !**************************************************************************
  !		DO the fit 					/ jp index
  !**************************************************************************
	call minimization1(nbp,wgrid,fnc2PI,conv_factor,nrm,meig,eigfreqk,Ak2,tauk2)

	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak2,tauk2,res)

  !**************************************************************************
  !		Print out results 				/ jp index
  !**************************************************************************
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw",strNjp,4)

	call write_real_array1(meig,histo,Ak2,"fit_coef_Ak",strNjp,4)
	call write_real_array1(meig,histo,tauk2,"fit_coef_tauk",strNjp,4)


 
  !**************************************************************************
  !		POST TREATMENT tau_k should be independent of index ip,jp
  !**************************************************************************
	do n = 1, meig
		if ( tauk1(n).lt.tauk2(n) ) then
			tauk(n) = tauk1(n)
		else
			tauk(n) = tauk2(n)
		endif
	enddo

!	tauk = 0.5d0 * (tauk1 + tauk2 )

	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak1,tauk,res)
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_avtauk",strNip,4)
	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak2,tauk,res)
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_avtauk",strNjp,4)


  !**************************************************************************
  !		ANOTHER POSSIBILITY OF 
  !		POST TREATMENT for Amplitudes and tau_k
  !		[18.01.14]
  !**************************************************************************
	do n = 1, meig
		modC1(n) = Ak1(n) / tauk1(n)
		modC2(n) = Ak2(n) / tauk2(n)
	enddo

	Ak = modC1 * tauk

	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_modCavtauk",strNip,4)

	Ak = modC2 * tauk

	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_modCavtauk",strNjp,4)

  !**************************************************************************
  !		POST TREATMENT for off-diag elements ip <> jp
  !**************************************************************************
	do n = 1, meig
!		Ak(n) = dsqrt( Ak1(n) * Ak2(n) )
		Ak(n) = dsqrt( modC1(n) * modC2(n) ) * tauk(n)
	enddo
	call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
	strNipNjp(1:4) = strNip(1:4)
	strNipNjp(5:8) = strNjp(1:4)
	call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_ovw_IJ_avtauk",strNipNjp,8)

  !**************************************************************************
  !		END least-square fitting of chose diag PI matrix elements 
  !**************************************************************************
!  enddo

!	write(*,*) "READY for NEXT STEP?"
!	read(5,*) ii

  !**************************************************************************
  !		while loop for calc and print out local DOS_bb' / PIbb'
  !             matrix elements
  !**************************************************************************

!  ip = 666
!  jp = 666
!  do while ( (ip.ne.0).and.(jp.ne.0) )
   	if (ip.ne.jp) then

!	write(*,*) "Site/coord (ilat/xyz) indices ip, jp = ?"
!	read(5,*) ip,jp

!		open(unit=22,file="abcoef_odgp.dat")
!		open(unit=23,file="abcoef_odgm.dat")

!		open(unit=30,file="projDOSm_vs_eig.dat")
!		open(unit=31,file="projDOSp_vs_eig.dat")

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
!		do n = 1, Nbcoef
!        		write(*,*) "Coef a(",n,") = ", Arec(n), " #    b(",n,") = ", Brec(n)  
!			write(22,*) n, Arec(n), Brec(n)
!     		enddo
!     		write(*,*) " "

     		do j = 1, nbp-1
        		w = dw * j
        		call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!			projDOSp(j) = - dimag(zw) / pi
			projDOSp(j) = - dimag(zw) / ( pi * wgrid(j) )
     		enddo

		xn = 0
        	xn(ip) = 1.d0 / sqrt(2.0)
        	xn(jp) = -1.d0 / sqrt(2.0)
 
		call Lanczos_rec(tNdim,VLLp,xn,Nbcoef,Arec,Brec)

		write(*,*) "Antisymmetric superposition ip - jp"
!		do n = 1, Nbcoef
!        		write(*,*) "Coef a(",n,") = ", Arec(n), " #    b(",n,") = ", Brec(n) 
!			write(23,*) n, Arec(n), Brec(n)
!     		enddo
!     		write(*,*) " "

     		do j = 1, nbp-1
        		w = dw * j
        		call cont_frac(Nbcoef,Arec,Brec,w,Typecontfrac,dE,zw) 
!			projDOSm(j) = - dimag(zw) / pi
			projDOSm(j) = - dimag(zw) / ( pi * wgrid(j) )
     		enddo

     		do j = 1, nbp-1
!        		w = dw * j
!                	tmp2 = sqrt(w)
			tmp1 = ( projDOSp(j) - projDOSm(j) ) / 2.d0
!			fnc1PI(j) = tmp1 / wgrid(j)
			fnc1PI(j) = tmp1
!        		write(25,*) w, tmp1
!        		write(26,*) tmp2, tmp1
!        		write(27,*) tmp2, tmp1 / tmp2
!			write(30,*) w, projDOSp(j)
!			write(31,*) w, projDOSm(j)
     		enddo
 	
		call write_real_array2(nbp-1,wgrid,fnc1PI,"imDbbrecur_ovw_IJ",strNipNjp,8)

!		close(22)
!		close(23)

!		close(30)
!		close(31)
  !**************************************************************************
  !		DO the fit 					/ jp index
  !**************************************************************************
!		rdiff = fnc1PI - res
!		conv_factor = dsqrt( ddot(nbp,rdiff,1,rdiff,1) )
!		nrm = dsqrt( ddot(nbp,fnc1PI,1,fnc1PI,1) )
!        	conv_factor = conv_factor / nrm 
!		do n = 1, meig, 2
!			Ak(n) = - Ak(n)
!		enddo

 		Ak1 = Ak

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
 !               	write(*,*) "peak w_k(", n,")=", eigfreqk(n)
!			write(*,*) "between ",wgrid(ii-1)," and ",wgrid(ii)
!			write(*,*) "PI(ii-1)= ",fnc1PI(ii-1)," and PI(ii)=",fnc1PI(ii)
			if ( fnc1PI(ii-1) .lt. 0.d0 ) then
				Ak(n) = - Ak(n)
!				write(*,*) "Negative amplitude: Ak =",Ak(n)
			endif
!			write(*,*) " "
        	enddo

!		Ak = Ak1
		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
		call write_real_array2(nbp-1,wgrid,res,"initfit_imDbbrecur_ovw_IJ_avtauk",strNipNjp,8)

		conv_factor = 100.d0 
		call minimization2(nbp,wgrid,fnc1PI,conv_factor,nrm,meig,eigfreqk,Ak,tauk)

		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
		call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_IJ_ovw_min2",strNipNjp,8)

		Ak = Ak1
		conv_factor = 100.d0 
		call minimization3(nbp,wgrid,fnc1PI,conv_factor,nrm,meig,extra,eigfreqk,Ak,tauk)

		call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk,res)
		call write_real_array2(nbp-1,wgrid,res,"fit_imDbbrecur_IJ_ovw_min3",strNipNjp,8)



  !**************************************************************************
  !		END while loop for calc and print out local DOS_bb' / PIbb'
  !**************************************************************************
 	endif 
!  enddo





  !*************************************************************************
  !     DEALLOCATE ARRAYS

  deallocate( VLLp, xnp1, xn, xnm1, vectmp1, vectmp2, Mm )
  deallocate( Arec, Brec, lapackDbb, work, eigenval )
  deallocate( projDOSp, projDOSm )


end program polar_matrix_PIbb




