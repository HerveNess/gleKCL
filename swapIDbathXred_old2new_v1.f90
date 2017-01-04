program swapID_old2new
  !
  !
  !____________________________________________________________________________
  !
  !Program to calculate new ID for bath(L,R)reduced form old -> new files 
  !
  ! Date started: 30 Mar 2016
  !
  !____________________________________________________________________________

  implicit none

  integer :: nbtotatom
  integer :: i, j, Natom_old, Natom_new
  
  logical :: found

  real :: xshift, yshift, zshift, tmpx, tmpy, tmpz, tmpdist

  integer     , allocatable ::  old_id(:), new_id(:)
  real        , allocatable ::  oldx(:), oldy(:), oldz(:)
  real        , allocatable ::  newx(:), newy(:), newz(:)

  !**************************************************************************
  !**************************************************************************
  nbtotatom = 200
  !**************************************************************************
  !**************************************************************************
  allocate( old_id(nbtotatom), new_id(nbtotatom) )
  allocate( oldx(nbtotatom), oldy(nbtotatom), oldz(nbtotatom) )
  allocate( newx(nbtotatom), newy(nbtotatom), newz(nbtotatom) )
  !**************************************************************************
  !**************************************************************************


  !**************************************************************************
  !  READ input parameters
  !**************************************************************************

	write(*,*) "read OLD file covelocXforcord_bathLred_at_sort.dat"
        open(unit=27,file="old_coord_bathXred_at_sort.dat")
	read(27,*) 
	read(27,*) 
	read(27,*) 
	read(27,*) Natom_old
	read(27,*)
	read(27,*)
	read(27,*)
	read(27,*)
	read(27,*)
	do i = 1, Natom_old
		read(27,*) old_id(i), oldx(i), oldy(i), oldz(i)
	enddo
	close(27)

	write(*,*) "read NEW file coord_bathLred_at_sort.dat"
        open(unit=27,file="new_coord_bathXred_at_sort.dat")
	read(27,*) 
	read(27,*) 
	read(27,*) 
	read(27,*) Natom_new
	read(27,*)
	read(27,*)
	read(27,*)
	read(27,*)
	read(27,*)
	do i = 1, Natom_new
		read(27,*) new_id(i), newx(i), newy(i), newz(i)
	enddo
	close(27)

	if (Natom_old.ne.Natom_new) then
		write(*,*) "NOT the same Nb of atoms - STOP"
		stop
	endif

	write(*,*) "Enter shift in x,y,z-dir"
	write(*,*) "x-shift = "
	read(5,*) xshift
	write(*,*) "y-shift = "
	read(5,*) yshift
	write(*,*) "z-shift = "
	read(5,*) zshift

	open(unit=28,file="ID_old2new_coord_bathXred_at_sort.dat")
	write(28,*) Natom_new

	do i = 1, Natom_old
		found = .FALSE.
		do j = 1, Natom_new
			tmpx = oldx(i) + xshift
			tmpy = oldy(i) + yshift
			tmpz = oldz(i) + zshift
			
			tmpdist = ( tmpx - newx(j) )**2 + ( tmpy - newy(j) )**2 + ( tmpz - newz(j) )**2 
			tmpdist = sqrt( tmpdist )
			if ( tmpdist.lt.1.e-6 ) then
				found = .TRUE.
				write(28,*) old_id(i), oldx(i), oldy(i), oldz(i),"  new ID = ", new_id(j), tmpx, tmpy, tmpz 		
			endif
		enddo
		if (found.eqv..FALSE.) then
			write(*,*) "Atom old ID =",old_id(i)," coord ", oldx(i), oldy(i), oldz(i)
			write(*,*) "NOT found !"
			stop
		endif
	enddo
	close(28)
  !*************************************************************************
  !     DEALLOCATE ARRAYS

  deallocate( old_id, new_id, oldx, oldy, oldz, newx, newy, newz )


end program swapID_old2new





