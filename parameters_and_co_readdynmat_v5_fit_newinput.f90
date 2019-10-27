module parameters_and_co
  !!****h* polar_matrix/parameters
  !!   NAME
  !!     parameters_and_co
  !!   SYNOPSIS
  !!     Module containing parameter definitions and subroutines Mat/Vec multiplication
  !!   FUNCTION
  !!     Module to read input parameters and set up global values
  !!
  !!   INPUTS
  !!     none
  !!    
  !!   OUTPUT
  !!     see separate subroutines
  !!    
  !!    
  !!    
  !!    
  !!     
  !!   
  !!   NOTES
  !!     
  !!   BUGS
  !!     
  !!   SEE ALSO
  !!     
  !!*** 
  implicit none
  
  integer :: Nread, Ndim, tNdim, Nbcoef, Typecontfrac, nbp, type_extra
  integer :: LAMMPSunits, segment

  real(kind=8) :: dE 
  real(kind=8) :: tmp, extra_coef, scaling_amplitude, scaling_tauks

  real(kind=8), parameter :: pi = 3.14159265358979323846
    
contains
  


  subroutine read_input_parameters
  !!****f* sssm-main/parameters/read_input_parameters
  !!   NAME
  !!     read_input_parameters
  !!   SYNOPSIS
  !!     Read the input parameters
  !!   FUNCTION
  !!     Subroutine to  read the input parameters from a user-generated file, 
  !!     default name "input.in".
  !!
  !!   INPUTS
  !!     External file "input.in"
  !!     Variables are:
  !!     *   Ndim : dimension of matrix NxN and vectors N
  !!     *   Nbcoef : Nb of coef for recursion method
  !!
  !!   OUTPUT
  !!     No local variables passed
  !!    
  !!    
  !!    
  !!    
  !!     
  !!   
  !!   NOTES
  !! 
  !!   BUGS
  !!     
  !!   SEE ALSO
  !!     
  !!*** 
    !**************************************************************************
    !     READ INPUTS
    !
    !***************************************************************************

    implicit none

    ! read input parameters
    write (*,*) "Reading input variables"
   
    open(unit=25,file="input_new_rdynmat.in")
    
    read(25,*) LAMMPSunits
    
    read(25,*) Nbcoef
    read(25,*) Typecontfrac
    read(25,*) dE
    read(25,*) nbp
    read(25,*) segment
    read(25,*) type_extra, extra_coef
    read(25,*) scaling_amplitude
    read(25,*) scaling_tauks

    scaling_amplitude = scaling_amplitude / 100.0

    close(25)

    ! echo variables
    write (*,*) "  Units in LAMMPS = ",LAMMPSunits
    write (*,*) "  Nb recur coef = ",Nbcoef
    write (*,*) "  Nb grid pts   = ",nbp
    write (*,*) "  Length segment   = ",segment
    write (*,*) "  extra peaks = ",type_extra
    write (*,*) "  extra coef = ",extra_coef
    write (*,*) "  Scaling amplitude for init_fit = ",scaling_amplitude
    write (*,*) "  Scaling tau_k for init_fit = ",scaling_tauks
    write (*,*) "  "
  

  end subroutine read_input_parameters





!********************************************************************************
!*	Matrix / vector multiplication  
!********************************************************************************

    subroutine multMvec(Ndim,A,x,b)
  
    implicit none

    !Input variables
    integer, intent(in) :: Ndim
    real(kind=8), intent(in) :: A(Ndim,Ndim), x(Ndim)
    real(kind=8), intent(out) :: b(Ndim)
    
    ! Local variables
    integer :: i, j
    real(kind=8) :: tmp

    do i = 1, Ndim
       tmp = 0.d0
       do j = 1, Ndim
          tmp = tmp + A(i,j) * x(j)
       enddo
       b(i) = tmp
    end do
    
    end subroutine multMvec





!********************************************************************************
!********************************************************************************
    subroutine multvecvec(Ndim,u,v,tmp)
  
    implicit none

    !Input variables
    integer, intent(in) :: Ndim
    real(kind=8), intent(in) :: u(Ndim), v(Ndim)
    real(kind=8), intent(out) :: tmp
    
    ! Local variables
    integer :: i

    tmp = 0.d0
    do i = 1, Ndim
       tmp = tmp + u(i) * v(i)
    end do
    
    end subroutine multvecvec





!********************************************************************************
!********************************************************************************
    subroutine Lanczos_rec(Ndim,M,x0,Nbcoef,a,b)
  
    implicit none

    !Input variables
    integer, intent(in) :: Ndim, Nbcoef
    real(kind=8), intent(in) :: M(Ndim,Ndim), x0(Ndim)
    real(kind=8), intent(out) :: a(Ndim), b(Ndim)
    
    ! External functions
    real(kind=8) :: ddot

    ! Local variables
    integer :: n
    real(kind=8) :: tmp1
    real(kind=8) :: xn(Ndim), xnp1(Ndim), xnm1(Ndim), vectmp1(Ndim)

  !**************************************************************************
  !		Recursion loop
  !**************************************************************************
	xn = x0

  !********* CALC M * x_1 ***************************************************
	call dgemv("N",Ndim,Ndim,1.d0,M,Ndim,xn,1,0.d0,vectmp1,1)

  !********* CALC Mat a_1 ***************************************************
	a(1) = ddot(Ndim,xn,1,vectmp1,1)

  !********* CALC x_2 *****************************************************
	xnp1 = vectmp1 - a(1) * xn
     
  !********* CALC Mat b_1 ***************************************************
	tmp1 = ddot(Ndim,xnp1,1,xnp1,1)
	b(1) = dsqrt(tmp1)
        
  !********* Normalise x_2 **************************************************
  	xnp1 = xnp1 / b(1)

  	do n = 2, Nbcoef
  !**************************************************************************
  !		Swap vectors for next step
  !**************************************************************************
		xnm1 = xn
     		xn = xnp1
  !**************************************************************************

     !********* CALC M * x_n ***************************************************
     		call dgemv("N",Ndim,Ndim,1.d0,M,Ndim,xn,1,0.d0,vectmp1,1)

     !********* CALC Mat a_n ***************************************************
     		a(n) = ddot(Ndim,xn,1,vectmp1,1)

     !********* CALC x_n+1 *****************************************************
     		xnp1 = vectmp1 - a(n) * xn - b(n-1) * xnm1

     !********* CALC Mat b_n ***************************************************
     		tmp1 = ddot(Ndim,xnp1,1,xnp1,1)
     		b(n) = dsqrt(tmp1)

     !********* Normalise x_n+1 ************************************************
     		xnp1 = xnp1 / b(n)

  	enddo
  !**************************************************************************
  !		END Lanczos tri-diagonalisation of Dbb
  !**************************************************************************    
    end subroutine Lanczos_rec





!********************************************************************************
!********************************************************************************
    subroutine der_ppot(Ntype,R,res,resder1,resder2,Ar,rho,Br,cutoff)
  
    implicit none

    !Input variables
    integer, intent(in) :: Ntype
    real(kind=8), intent(in) :: R, Ar, rho, Br, cutoff
    real(kind=8), intent(out) :: res, resder1, resder2
    
    ! Local variables
    real(kind=8) :: s6, x7, x14, u, power6

!********************************************************************************
! 	use conventional notation
!	Ntype = 1 : Lennard-Jones as	 4 e [ (s/X)**12 - (s/X)**6 ]
!	Ntype = 2 : Buckingham as	 A e^{-X/rho} - C/X**6 ]
!********************************************************************************

	s6 = Br**6
   	x7 = R ** 7
        x14 = x7**2
    	u = R / rho
   	power6 = 6.d0 * Br / x7

	if ((Ntype.eq.1).or.(Ntype.eq.2)) then
		res = 4.d0 * Ar * ( s6/R**6 - 1.d0 ) * s6 / R**6
		resder1 = 24.d0 * Ar * s6 * ( 1/x7 - 2.d0 * s6 * R / x14 )
		resder2 = 24.d0 * Ar * s6 * ( 26.d0 * s6 / x14 - 7.d0 / (R * x7) )

		if ((Ntype.eq.2).and.(R.gt.cutoff)) then
			res = 0.d0
			resder1 = 0.d0
			resder2 = 0.d0
		endif
    	else
		res = Ar * dexp(-u) - Br * R / x7   
		resder1 = -u * Ar * dexp(-u) + power6
		resder2 =  u * u * Ar * dexp(-u) - 7.d0 * power6 / R
    	endif
    
    end subroutine der_ppot





!********************************************************************************
!********************************************************************************
      	subroutine cont_frac(localNcoef,Arec,Brec,w,Ntype,dE,zw) 
      
    	implicit none

    !Input variables
    	integer, intent(in) :: localNcoef, Ntype
    	real(kind=8), intent(in) :: Arec(Nbcoef), Brec(Nbcoef), dE, w
    	complex(kind=8), intent(out) :: zw
    
    ! Local variables
      	integer i,j, n
    	real(kind=8) :: einf, esup, x, y, eim, binfty, delta
    	complex(kind=8) :: rr


	Binfty = Brec(localNcoef) * 1.02d0

      	einf = Arec(localNcoef) - 2.d0 * Binfty 
      	esup = Arec(localNcoef) + 2.d0 * Binfty 

        delta = ( w - Arec(localNcoef) )**2 - 4.d0 * Binfty**2

      	if (Ntype.eq.0) then
         	x = ( w - Arec(localNcoef) )
         	y = dE
        	eim = dE
      	endif

      	if (Ntype.eq.1) then
         	x = ( w - Arec(localNcoef) ) / 2.d0
!         	eim = 0.d0
         	eim = dE

         	if (w.gt.esup) then
            		x = x + dsqrt( delta ) / 2.d0
            		y = 0.d0
         	else
            		if (w.lt.einf) then
               			x = x - dsqrt( delta ) / 2.d0
               			y = 0.d0
            		else
               			y = dsqrt( - delta ) / 2.d0
            		endif
         	endif
      	endif

      	rr = dcmplx(x,y)
      	do i = localNcoef-1, 1, -1
         	rr = dcmplx( w - Arec(i) , 0.d0 ) - dcmplx( Brec(i) * Brec(i) , 0.d0 ) / rr
         	rr = rr + dcmplx(0.0,eim)
      	enddo
      	zw = 1.d0 / rr
      

      	end subroutine cont_frac





!********************************************************************************
!********************************************************************************
      	subroutine cont_frac_inv(localNcoef,Arec,Brec,w,Ntype,dE,zw) 
      
    	implicit none

    !Input variables
    	integer, intent(in) :: localNcoef, Ntype
    	real(kind=8), intent(in) :: Arec(Nbcoef), Brec(Nbcoef), dE, w
    	complex(kind=8), intent(out) :: zw
    
    ! Local variables
      	integer i,j, n
    	real(kind=8) :: einf, esup, x, y, eim, binfty, delta
    	complex(kind=8) :: rr


	Binfty = Brec(localNcoef) * 1.02d0

      	einf = Arec(localNcoef) - 2.d0 * Binfty 
      	esup = Arec(localNcoef) + 2.d0 * Binfty 

        delta = ( w - Arec(localNcoef) )**2 - 4.d0 * Binfty**2

      	if (Ntype.eq.0) then
         	x = ( w - Arec(localNcoef) )
         	y = dE
        	eim = dE
      	endif

      	if (Ntype.eq.1) then
         	x = ( w - Arec(localNcoef) ) / 2.d0
!         	eim = 0.d0
         	eim = dE

         	if (w.gt.esup) then
            		x = x + dsqrt( delta ) / 2.d0
            		y = 0.d0
         	else
            		if (w.lt.einf) then
               			x = x - dsqrt( delta ) / 2.d0
               			y = 0.d0
            		else
               			y = dsqrt( - delta ) / 2.d0
            		endif
         	endif
      	endif

      	rr = dcmplx(x,y)
      	do i = localNcoef-1, 1, -1
         	rr = dcmplx( w - Arec(i) , 0.d0 ) - dcmplx( Brec(i) * Brec(i) , 0.d0 ) / rr
         	rr = rr + dcmplx(0.0,eim)
      	enddo
      	zw = rr
      

      	end subroutine cont_frac_inv





!********************************************************************************
!********************************************************************************
      	subroutine cont_frac2(localNcoef,Arec,Brec,Maxw2,w,Ntype,dE,zw) 
      
    	implicit none

    !Input variables
    	integer, intent(in) :: localNcoef, Ntype
    	real(kind=8), intent(in) :: Arec(Nbcoef), Brec(Nbcoef), dE, w, Maxw2
    	complex(kind=8), intent(out) :: zw
    
    ! Local variables
      	integer i,j, n
    	real(kind=8) :: einf, esup, x, y, eim, binfty, delta, emid
    	complex(kind=8) :: rr


	Binfty = Maxw2 * 1.02d0

      	einf = 0.d0 
      	esup = Maxw2
        emid = 0.d0 

        delta = ( w - emid )**2 - 4.d0 * Binfty**2

      	if (Ntype.eq.0) then
         	x = ( w - Arec(localNcoef) )
         	y = dE
         	eim = dE
      	endif

      	if (Ntype.eq.1) then
         	x = ( w - emid ) / 2.d0
         	eim = dE

         	if (w.gt.esup) then
            		x = x + dsqrt( delta ) / 2.d0
            		y = 0.d0
         	else
            		if (w.lt.einf) then
               			x = x - dsqrt( delta ) / 2.d0
               			y = 0.d0
            		else
               			y = dsqrt( - delta ) / 2.d0
            		endif
         	endif
      	endif

      	rr = dcmplx(x,y)
      	do i = localNcoef-1, 1, -1
         	rr = dcmplx( w - Arec(i) , 0.d0 ) - dcmplx( Brec(i) * Brec(i) , 0.d0 ) / rr
         	rr = rr + dcmplx(0.0,eim)
      	enddo
      	zw = 1.d0 / rr
      

      	end subroutine cont_frac2


!********************************************************************************
!********************************************************************************
      	subroutine fnc_fit(Ngrid,xi,Nvar,wk,Xamp,Xtau,out) 
      
    	implicit none

    !Input variables
    	integer, intent(in) :: Ngrid, Nvar
    	real(kind=8), intent(in) :: xi(Ngrid), wk(Nvar), Xamp(Nvar), Xtau(Nvar)
    	real(kind=8), intent(out) :: out(Ngrid)
    
    ! Local variables
      	integer i, n
    	real(kind=8) :: tmp1, tmp2, tmp2b


      	do i = 1, Ngrid-1
        	tmp1  = 0.d0
        	do n = 1, Nvar
           		tmp2b = 1.d0 + ( xi(i) - wk(n) )**2 * Xtau(n)**2 
                        tmp2  = Xamp(n) / tmp2b

!		REMEMBER THAT IS FOR THE FITTING PROCEDURE
!		CHECK REAL EXPRESSION FOR \PI_bb in terms of Lorentzians
!
!                        tmp2  = Xamp(n) * Xtau(n) / tmp2b

           		tmp1  = tmp1  + tmp2  
        	enddo
		out(i) = tmp1
     	enddo

     

      	end subroutine fnc_fit
      
!  ===================================================================



!********************************************************************************
!********************************************************************************
      	subroutine extra_peaks(meig,eigfreqk,type_extra,extra_coef,extra)
      
    	implicit none

    !Input variables
    	integer, intent(in) :: meig, type_extra
	real(kind=8), intent(in) :: extra_coef
    	integer, intent(out) :: extra
    
    ! Local variables
      	integer n
    	real(kind=8) :: tmp1, tmp2, eigfreqk(nbp)


	extra = meig * extra_coef
	tmp2 = eigfreqk(1)*0.5d0

        if (type_extra.eq.1) then
	
		write(*,*) "Extra peaks: linear spacing"
	        tmp1 = ( 1.01d0*eigfreqk(meig) - tmp2 ) / (extra - 1)
		do n = 1, extra
           		eigfreqk(n+meig) = tmp1 * n + tmp2
		enddo

	else
	
		write(*,*) "Extra peaks: quadratic space / more points close to 0"
	        tmp1 = ( 1.01d0*eigfreqk(meig) - tmp2 ) / (extra - 1)**2
		do n = 1, extra
           		eigfreqk(n+meig) = tmp1 * (n-1)**2 + tmp2
		enddo


	endif
	
      	end subroutine extra_peaks
      
!  ===================================================================




!********************************************************************************
!********************************************************************************
	subroutine minimization1(nbp,wgrid,fncPI,conv_factor,nrm,meig,eigfreqk,Ak,tauk)

    	implicit none

    !Input variables
    	integer, intent(in) :: nbp, meig
    	real(kind=8), intent(in) :: nrm
    	real(kind=8), intent(in) :: wgrid(nbp), fncPI(nbp), eigfreqk(nbp)
	
    ! Local variables
      	integer k, n, kk
    	real(kind=8) :: tmp1, tmp2, percent, Delta_conv, conv_factor_old, conv_factor
	real(kind=8) :: Ak_plus(nbp), tauk_plus(nbp), res(nbp), rdiff_plus(nbp) 

    	real(kind=8) :: Ak(nbp), tauk(nbp) 


  !********** external functions ********************************************
  	real(kind=8) :: ddot



  !**************************************************************************
  !		START OUTER LOOP FOR THE FIT 
  !**************************************************************************
	Delta_conv = 10.d0

        percent = 2.d0

        kk = 0
	
        Ak_plus = Ak
	tauk_plus = tauk
	
	do while (kk.le.15)
		write(*,*) "Outerloop Iteration kk = ", kk

  !**************************************************************************
  !		START INNER LOOP : fit with changing Ak parameters 
  !**************************************************************************
	k = 0
!        percent = 2.d0
	do while (Delta_conv.gt.1.d-8)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on Ak : k = ", k
		write(*,*) "percent = ", percent
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

			do n = 1, meig
				write(*,*) " n= ", n

				Ak_plus(n)   = Ak(n) * (1.d2 + percent)/1.d2
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
!				write(*,*) "Min along Ak path +1% / conv = ", tmp1

				if (tmp1.lt.conv_factor) then 
!					write(*,*) "Min OK"
					Ak(n) = Ak_plus(n)
					conv_factor = tmp1
				else
					Ak_plus(n)   = Ak(n) * (1.d2 - percent)/1.d2
					call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
					rdiff_plus = fncPI - res
					tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
					tmp1 = tmp1 / nrm
!					write(*,*) "Min along Ak path -1% / conv = ", tmp1

					if (tmp1.lt.conv_factor) then 
!						write(*,*) "Min OK"
						Ak(n) = Ak_plus(n)
						conv_factor = tmp1
					endif
				endif

			enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
		Delta_conv = conv_factor_old - conv_factor

!		percent = percent / 2.d0

	enddo
  !**************************************************************************
  !		END INNER LOOP on Ak parameters 
  !**************************************************************************
!		do j = 1, nbp-1
!        			write(27,*) wgrid(j), res(j)
!     		enddo
!		write(27,*) " "


  !**************************************************************************
  !		START INNER LOOP : fit with changing Ak & tauk parameters 
  !**************************************************************************
	k = 0
!	percent = 2.d0
	Delta_conv = 1.d0
	do while (Delta_conv.gt.1.d-8)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on tauk : k = ", k
		write(*,*) "percent = ", percent
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

		do n = 1, meig
			write(*,*) " n= ", n

			Ak_plus(n)   = Ak(n) * (1.d2 + percent)/1.d2
			tauk_plus(n) = tauk(n) * (1.d2 + percent)/1.d2
			call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk_plus,res)
			rdiff_plus = fncPI - res
			tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
			tmp1 = tmp1 / nrm
!			write(*,*) "Min along path Ak + 1% & tauk +1% / conv = ", tmp1

			if (tmp1.lt.conv_factor) then 
!				write(*,*) "Min OK"
				tauk(n) = tauk_plus(n)
				Ak(n)   = Ak_plus(n)
				conv_factor = tmp1
			else
				Ak_plus(n)   = Ak(n) * (1.d2 + percent)/1.d2
				tauk_plus(n) = tauk(n) * (1.d2 - percent)/1.d2
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk_plus,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
!				write(*,*) "Min along path Ak + 1% & tauk -1% / conv = ", tmp1

				if (tmp1.lt.conv_factor) then 
!					write(*,*) "Min OK"
					tauk(n) = tauk_plus(n)
					Ak(n)   = Ak_plus(n)
					conv_factor = tmp1
				else
					Ak_plus(n)   = Ak(n) * (1.d2 - percent)/1.d2
					tauk_plus(n) = tauk(n) * (1.d2 + percent)/1.d2
					call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk_plus,res)
					rdiff_plus = fncPI - res
					tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
					tmp1 = tmp1 / nrm
!					write(*,*) "Min along path Ak - 1% & tauk +1% / conv = ", tmp1

					if (tmp1.lt.conv_factor) then 
!						write(*,*) "Min OK"
						tauk(n) = tauk_plus(n)
						Ak(n)   = Ak_plus(n)
						conv_factor = tmp1
					else
						Ak_plus(n)   = Ak(n) * (1.d2 - percent)/1.d2
						tauk_plus(n) = tauk(n) * (1.d2 - percent)/1.d2
						call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk_plus,res)
						rdiff_plus = fncPI - res
						tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
						tmp1 = tmp1 / nrm
!						write(*,*) "Min along path Ak - 1% & tauk -1% / conv = ", tmp1

						if (tmp1.lt.conv_factor) then 
!							write(*,*) "Min OK"
							tauk(n) = tauk_plus(n)
							Ak(n)   = Ak_plus(n)
							conv_factor = tmp1
						endif

					endif

				endif
			endif
		enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
!	percent = percent / 2.d0
		Delta_conv = conv_factor_old - conv_factor

	enddo
  !**************************************************************************
  !		END INNER LOOP on Ak & tauk parameters 
  !**************************************************************************

!		do j = 1, nbp-1
!        		write(27,*) wgrid(j), res(j)
!     		enddo
!		write(27,*) " "
  !**************************************************************************
  !		START INNER LOOP : fit with changing tauk parameters 
  !**************************************************************************
	k = 0
!	percent = 2.d0
	Delta_conv = 1.d0
	do while (Delta_conv.gt.1.d-8)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on tauk : k = ", k
		write(*,*) "percent = ", percent
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

		do n = 1, meig
			write(*,*) " n= ", n

			tauk_plus(n) = tauk(n) * (1.d2 + percent)/1.d2
			call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk_plus,res)
			rdiff_plus = fncPI - res
			tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
			tmp1 = tmp1 / nrm
!			write(*,*) "Min along tauk path +1% / conv = ", tmp1

			if (tmp1.lt.conv_factor) then 
!				write(*,*) "Min OK"
				tauk(n) = tauk_plus(n)
				conv_factor = tmp1
			else

				tauk_plus(n) = tauk(n) * (1.d2 - percent)/1.d2
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak,tauk_plus,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
!				write(*,*) "Min along tauk path -1% / conv = ", tmp1

				if (tmp1.lt.conv_factor) then 
!					write(*,*) "Min OK"
					tauk(n) = tauk_plus(n)
					conv_factor = tmp1
				endif
			endif
		enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
		Delta_conv = conv_factor_old - conv_factor
!	percent = percent / 2.d0
	enddo
  !**************************************************************************
  !		END INNER LOOP on tauk parameters 
  !**************************************************************************

!		do j = 1, nbp-1
!        		write(27,*) wgrid(j), res(j)
!     		enddo
!		write(27,*) " "




	kk = kk + 1 
	percent = percent / 2.d0
	enddo
!	close(27)
  !**************************************************************************
  !		END OUTER LOOP for fit 
  !**************************************************************************
	
      	end subroutine minimization1
      
!  ===================================================================


!********************************************************************************
!********************************************************************************
	subroutine minimization2(nbp,wgrid,fncPI,conv_factor,nrm,meig,eigfreqk,Ak,tauk)

    	implicit none

    !Input variables
    	integer, intent(in) :: nbp, meig
    	real(kind=8), intent(in) :: nrm
    	real(kind=8), intent(in) :: wgrid(nbp), fncPI(nbp), eigfreqk(nbp), tauk(nbp)
	
    ! Local variables
      	integer k, n, kk
    	real(kind=8) :: tmp1, tmp2, percent, Delta_conv, conv_factor_old, conv_factor
	real(kind=8) :: Ak_plus(nbp), res(nbp), rdiff_plus(nbp) 

    	real(kind=8) :: Ak(nbp)


  !********** external functions ********************************************
  	real(kind=8) :: ddot



  !**************************************************************************
  !		START OUTER LOOP FOR THE FIT 
  !**************************************************************************
        kk = 0
	
        Ak_plus = Ak
	
	do while (kk.le.0)
		write(*,*) "Outerloop Iteration kk = ", kk

  !**************************************************************************
  !		START INNER LOOP : fit with changing sing of Ak parameters 
  !**************************************************************************
	k = 0

	Delta_conv = 1.d0

	do while (Delta_conv.gt.1.d-6)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on +- Ak : k = ", k
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

			do n = meig, 1, -1 
				write(*,*) " n= ", n

				Ak_plus(n)   = - Ak(n)
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
				write(*,*) "Min along +- Ak path = ", tmp1

				if (tmp1.lt.conv_factor) then 
					write(*,*) "Min OK"
					Ak(n) = Ak_plus(n)
					conv_factor = tmp1
				endif

			enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
		Delta_conv = conv_factor_old - conv_factor


	enddo
  !**************************************************************************
  !		END INNER LOOP on Ak parameters 
  !**************************************************************************
!		do j = 1, nbp-1
!        			write(27,*) wgrid(j), res(j)
!     		enddo
!		write(27,*) " "

  !**************************************************************************
  !		START INNER LOOP : fit with changing sing of 3-Ak parameters 
  !**************************************************************************
	k = 0

	Delta_conv = 1.d0

	do while (Delta_conv.gt.1.d-6)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on +- Ak : k = ", k
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

			do n = meig, 2, -1 
				write(*,*) " n= ", n

				Ak_plus(n)   = - Ak(n)
				Ak_plus(n-1)   = - Ak(n-1)
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
				write(*,*) "Min along +- 2 Ak path = ", tmp1

				if (tmp1.lt.conv_factor) then 
					write(*,*) "Min OK"
					Ak(n) = Ak_plus(n)
					Ak(n-1) = Ak_plus(n-1)
					conv_factor = tmp1
				endif

			enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
		Delta_conv = conv_factor_old - conv_factor


	enddo
  !**************************************************************************
  !		END INNER LOOP on Ak parameters 
  !**************************************************************************



	kk = kk + 1 

	enddo
  !**************************************************************************
  !		END OUTER LOOP for fit 
  !**************************************************************************
	
      	end subroutine minimization2
      
!  ===================================================================



!********************************************************************************
!********************************************************************************
	subroutine minimization3(nbp,wgrid,fncPI,conv_factor,nrm,meig,extra,eigfreqk,Ak,tauk)

    	implicit none

    !Input variables
    	integer, intent(in) :: nbp, meig, extra
    	real(kind=8), intent(in) :: nrm
    	real(kind=8), intent(in) :: wgrid(nbp), fncPI(nbp), eigfreqk(nbp), tauk(nbp)
	
    ! Local variables
      	integer k, n, nn, m, n_interv_start, n_interv, nb_interv
    	real(kind=8) :: tmp1, tmp2, percent, Delta_conv, conv_factor_old, conv_factor
	real(kind=8) :: Ak_plus(nbp), res(nbp), rdiff_plus(nbp) 

    	real(kind=8) :: Ak(nbp)


  !********** external functions ********************************************
  	real(kind=8) :: ddot



  !**************************************************************************
  !		START on the original peaks
  !**************************************************************************
	
        Ak_plus = Ak

	Delta_conv = 1.d0

	do while (Delta_conv.gt.1.d-6)
		conv_factor_old = conv_factor
		write(*,*) "Iteration on +- Ak (original peaks) : k = ", k
		write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

			do n = meig-extra, 1, -1 
				write(*,*) " n= ", n

				Ak_plus(n)   = - Ak(n)
				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
				write(*,*) "Min along +- Ak path = ", tmp1

				if (tmp1.lt.conv_factor) then 
					write(*,*) "Min OK"
					Ak(n) = Ak_plus(n)
					conv_factor = tmp1
				endif

			enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

		k = k + 1
		Delta_conv = conv_factor_old - conv_factor


	enddo


	n_interv_start = 1

	do n_interv = n_interv_start, 1, -1

		write(*,*) "Outerloop Iteration n_interv = ", n_interv

  !**************************************************************************
  !		START INNER LOOP : fit with changing sing of Ak parameters 
  !**************************************************************************
		k = 0

		Delta_conv = 1.d0

		do while (Delta_conv.gt.1.d-6)
			conv_factor_old = conv_factor
			write(*,*) "Iteration on +- Ak (extra peaks) : k = ", k
			write(*,*) "FOR FIT: convergence factor = ", conv_factor_old

			if (n_interv.ge.2) then
				nb_interv = extra / (n_interv-1)
			else
				nb_interv = extra
			endif 

			do n = 1, nb_interv
 
				write(*,*) " n= ", n

				do nn = 1, n_interv
					m = (n-1)*nb_interv + nn
					Ak_plus(m+(meig-extra))   = - Ak(m+(meig-extra))
				enddo

				call fnc_fit(nbp,wgrid,meig,eigfreqk,Ak_plus,tauk,res)
				rdiff_plus = fncPI - res
				tmp1 = dsqrt( ddot(nbp,rdiff_plus,1,rdiff_plus,1) )
				tmp1 = tmp1 / nrm
				write(*,*) "Min along +- grp of Ak path = ", tmp1

				if (tmp1.lt.conv_factor) then 
					write(*,*) "Min OK"

					do nn = 1, n_interv
						m = (n-1)*nb_interv + nn
						Ak(m+(meig-extra))   = Ak_plus(m+(meig-extra))
					enddo

					conv_factor = tmp1
				endif

			enddo

!		if (conv_factor.ge.conv_factor_old) then
!			write(*,*) "No minimisation along gradient path!!!"
!			stop
!		endif

			k = k + 1
			Delta_conv = conv_factor_old - conv_factor


		enddo
  !**************************************************************************
  !		END INNER LOOP on Ak parameters 
  !**************************************************************************
!		do j = 1, nbp-1
!        			write(27,*) wgrid(j), res(j)
!     		enddo
!		write(27,*) " "

	enddo
  !**************************************************************************
  !		END OUTER LOOP for fit 
  !**************************************************************************
	
      	end subroutine minimization3
      
!  ===================================================================



!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

	subroutine write_real_array1(nbp, w, array, filename, str, nlength)

  	implicit none
  
  	integer, intent(in) :: w(nbp)
  	real(kind=8), intent(in) :: array(nbp)
  	integer :: nbp, i, j, long, nlength
  	real :: tmp
  	character(len=*) ::  filename
  	character(len=80) ::  fnom
  	character :: str(80)


  	long=LEN(trim(filename))
  	fnom(1:long)=filename(1:long)
  	fnom(long+1:long+1)="."
  	do i=1,nlength
     		j=long+1+i
     		fnom(j:j)=str(i)
  	enddo
  	fnom(long+nlength+2:long+nlength+5)=".dat"

!  write(*,*) "fnom     = ", fnom(1:long+nlength+5)

     	open(unit=22,file=fnom(1:long+nlength+5))
     	do i = 1, nbp
        tmp = array(i)
        	write(22,*) w(i), tmp
     	enddo
	close(22)

	end subroutine write_real_array1

!--------------------------------------------------------------------------------

	subroutine write_real_array2(nbp, w, array, filename, str, nlength)

  	implicit none
  
  	real(kind=8), intent(in) :: array(nbp), w(nbp)
  	integer :: nbp, i, j, long, nlength
  	real :: tmp
  	character(len=*) ::  filename
  	character(len=80) ::  fnom
  	character :: str(80)


  	long=LEN(trim(filename))
  	fnom(1:long)=filename(1:long)
  	fnom(long+1:long+1)="."
  	do i=1,nlength
     		j=long+1+i
     		fnom(j:j)=str(i)
  	enddo
  	fnom(long+nlength+2:long+nlength+5)=".dat"

!  write(*,*) "fnom     = ", fnom(1:long+nlength+5)

     	open(unit=22,file=fnom(1:long+nlength+5))
     	do i = 1, nbp
        tmp = array(i)
        	write(22,*) w(i), tmp
     	enddo
	close(22)

	end subroutine write_real_array2
!********************************************************************************
!*  
!!     translate a integer value into string
!!     see http://www.ccl.net/cca/software/SOURCES/FORTRAN/string-library/libstr/string.shtml
!!
!!     Intstr - Translate a integer value into string
!!     Synopsis:
!!
!!	subroutine intstr(int,char,nchar)
!!	character *N char
!!	integer nchar,int
!!
!!      Description:
!!      This routine translates the int integer into char string. Nchar is
!!      the output string length to obtain a fixed size of output (example:
!!      001,002,....,023,etc). If nchar is zero, the routine generates a 
!!      string with the exact size of the int number 
!!      (example:1,2,....,3423,etc).
!!

      subroutine intstr(num,str,lun)

      integer num,num_tmp,j,n,lun
      character ::  str(80)
      character *1 cifra(10)
      logical segno

      data cifra /'0','1','2','3','4','5','6','7','8','9'/

!!     translate the integer num
      num_tmp=num
      do j=1,lun
         n=num_tmp/10**(lun-j)
         num_tmp=num_tmp-(n*10**(lun-j))
         str(j:j)=cifra(n+1)
      end do


      return
      end subroutine intstr
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

 


end module parameters_and_co

