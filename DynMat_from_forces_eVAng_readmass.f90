program DynMat_from_forces

	implicit none

  	integer :: i, j, l, n, ip, jp, alpha, Natom, Ndim
	integer :: id, im, jm
	real :: conversionfactor, convLP, delta, tmp1, fatom(3)
	real, allocatable ::  Force_plus(:,:), Force_minus(:,:), DynMat(:,:), mass(:)

	real(kind=8) :: tmp_dynmat


  !**************************************************************************
  !**************************************************************************
  !		START
  !**************************************************************************
  !**************************************************************************

	delta = 0.01

	conversionfactor = 3.757383e25

	convLP = 16.02 * 2.14417766e25

	open(unit=25,file='forces.dat')

  !  READ un-used 3 lines in Lammps forces.dat file output
	read(25,*)
	read(25,*)
	read(25,*)

  !  READ Nb atoms in Lammps forces.dat file output
	read(25,*) Natom
  	Ndim = 3 * Natom
	
  	allocate( Force_plus(Ndim,3), Force_minus(Ndim,3), DynMat(Ndim,Ndim) )
	allocate( mass(Natom) )

  !  READ un-used 5 + Natom lines in Lammps forces.dat file output
	do n = 1, 5+Natom
		read(25,*)
	end do

	do i = 1, Ndim
  !  READ un-used 9 lines in Lammps forces.dat file output
		do n = 1, 9
			read(25,*)
		end do
		do j = 1, Natom
			read(25,*) id, fatom(1), fatom(2), fatom(3), mass(j)
			Force_plus(j,:) = fatom
!			write(*,*) "deplac + : atom i=",j, " force = ", Force_plus(j,:)
		end do
  !  READ un-used 9 lines in Lammps forces.dat file output
		do n = 1, 9
			read(25,*)
		end do
		do j = 1, Natom
			read(25,*) id, fatom(1), fatom(2), fatom(3), mass(j)
			Force_minus(j,:) = fatom
!			write(*,*) "deplac - : atom i=",j, " force = ", Force_minus(j,:)
		end do

!		write(*,*) " "

		do l = 1, Natom
			do alpha = 1, 3
				j = 3*(l-1) + alpha
				tmp1 = ( Force_minus(l,alpha) - Force_plus(l,alpha) ) / (2.d0*delta)

!!				DynMat(i,j) = tmp1 * conversionfactor / mass

!!				DynMat(i,j) = tmp1 / sqrt( mass(i)*mass(j) )

				DynMat(i,j) = tmp1

!				write(*,*) "DynMat(",i,",",j,") = ", tmp1 * convLP
!				DynMat(i,j) = tmp1 * convLP
			end do
		end do

	end do	

	close(25)

  !**************************************************************************
  !**************************************************************************
  !		Post treatment
  !**************************************************************************
  !**************************************************************************
  !		Include masses
  !**************************************************************************
	do i = 1, Ndim
		im = (i-1)/3
		do j = 1, Ndim
			jm = (j-1)/3
			tmp1 = DynMat(i,j) / sqrt( mass(im+1) * mass(jm+1) )
			DynMat(i,j) = tmp1
		end do
	end do

  !**************************************************************************
	open(unit=35,file='dynmat_coeff_unsym.dat')
	write(35,*) Natom
	do i = 1, Ndim
		do j = 1, Ndim
			write(35,*) DynMat(i,j)
		end do
	end do
	close(35)


  !************** Symmetrize Dyn Mat ************************************
	do i = 1, Ndim
		do j = 1,i
			tmp1 = ( DynMat(i,j) + DynMat(j,i) ) * 0.5d0 
			DynMat(i,j) = tmp1
			DynMat(j,i) = tmp1
		end do
	end do

	open(unit=35,file='dynmat_coeff.dat')
	write(35,*) Natom
	do i = 1, Ndim
		do j = 1, Ndim
			write(35,*) DynMat(i,j)
		end do
	end do
	close(35)


  !************** CHECK Acoustic Sum Rule *******************************
	write(*,*) "ASR per xy-comp = "
	do i = 1, Natom
		ip = 3*(i-1) + 1
		tmp1 = 0.0
		do j = 1, Natom
			jp = 3*(j-1) + 2
			tmp1 = tmp1 + DynMat(ip,jp)
		end do
		write(*,*) "ASR line : atom i = ", i, " : xy-comp = ",tmp1
	end do


  !************** CHECK Acoustic Sum Rule *******************************
	write(*,*) "ASR full line-column = "
	do i = 1, Natom
		tmp1 = 0.0
		do j = 1, Natom
			tmp1 = tmp1 + DynMat(i,j)
		end do
		write(*,*) "ASR line i = ", i, " : = ",tmp1
	end do


  !*************************************************************************
  !     DEALLOCATE ARRAYS

  deallocate( Force_plus, Force_minus, DynMat )


end program DynMat_from_forces




