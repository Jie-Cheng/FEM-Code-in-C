program solidsolver
	
	! **********************************************************************
	! Main program to conduct the analysis
	! Currently 4 modes are available: static, dynamic, residual and debug
	! Mixed formulation is in use
	! **********************************************************************
	
	use read_file
	
	implicit none
	
	character(80) :: filepath
	integer :: ct, ct_rate, ct_max, ct1
	real(8) :: time_elapsed  
	call timestamp()
	filepath = '/Users/jiecheng/Documents/SolidResults/'
	call system_clock(ct,ct_rate,ct_max)
	call read_input(10, 'input.txt', mode, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, &
					materialprops, gravity, isbinary, penalty)
	call read_mesh(nsd, ned, nn, nel, nen, coords, connect, bc1, bc2, share)
	if (mode == 0) then 
		call statics(filepath)
	else if (mode == 1) then
		call dynamics(filepath)
	end if
	call system_clock(ct1)
	call timestamp()
	time_elapsed = dble(ct1-ct)/dble(ct_rate)
	if (time_elapsed < 60) then
		write(*,'("time elapsed:", f12.2, 3x, "seconds")'), time_elapsed
	else if (time_elapsed < 3600) then
		time_elapsed = time_elapsed/60
		write(*,'("time elapsed:", f12.2, 3x, "minutes")'), time_elapsed
	else 
		time_elapsed = time_elapsed/3600
		write(*,'("time elapsed:", f12.2, 3x, "hours")'), time_elapsed
	end if
	
end program solidsolver
	
subroutine statics(filepath)
		
	use read_file
	use integration
	use face
	use shapefunction
	use mixed_material
	use mixed_externalforce
	use symmetric_solver
	use mixed_output
	use mixed_internalforce
	use mixed_tangentstiffness
	
	implicit none
	
	integer :: i, nit, row, col, j, k
	real(8), allocatable, dimension(:) :: Fext, Fint, R, w, w1, dw
	real(8), allocatable, dimension(:,:) ::  A
	real(8) :: loadfactor, increment, err1, err2
	real(8), dimension(size(bc1, 2)) :: constraint
	real(8), dimension(size(bc1, 2), nn*ned) :: der_constraint
	character(80), intent(in) :: filepath
	
	allocate(Fext(nn*ned+nel))
	allocate(Fint(nn*ned+nel))
	allocate(R(nn*ned+nel))
	allocate(w(nn*ned+nel))
	allocate(w1(nn*ned+nel))
	allocate(dw(nn*ned+nel))
	allocate(A(nn*ned+nel,nn*ned+nel))
	
	! initialize w
	w = 0.
	constraint = 0.
	der_constraint = 0.
	
	nprint=1
	nsteps=1
	dt=1.
	loadfactor = 0.
	increment = firststep
	step = 0
	call write_results(filepath,w)
	
	! If the external load is traction, then the external force doesn't change
	if (size(bc2, 1) /= 3) then
		call force_traction(Fext)
	end if
	
	do while (loadfactor < 1.)
		step = step + 1
		if (loadfactor + increment > 1.) then
			increment = 1. - loadfactor
		end if
		w1 = w
		loadfactor  = loadfactor + increment
		err1 = 1.
		err2 = 1.
		nit = 0
		write(*,'(A, A, i5, 5x, A, e12.4, A)') repeat('=', 30), 'Step', step, 'Load', loadfactor, repeat('=', 36)
		do while (((err1 > tol) .or. (err2 > tol)) .and. (nit < maxit))
			nit = nit + 1
			Fint = force_internal(w)
			! If the external load is pressure, update the external force
			if (size(bc2, 1) == 3) then
				Fext = force_pressure(w)
			end if
			A = tangent_internal(w)
			R = Fint - loadfactor*Fext
			! penalty
			do i = 1, size(bc1,2)
				row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
				constraint(i) = w(row) - bc1(3,i)
				der_constraint(i,row) = 1. 
			end do
			do i = 1, size(bc1,2)
				row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
				A(row,row) = A(row,row) + penalty*der_constraint(i,row)*der_constraint(i,row)
				R(row) = R(row) + penalty*der_constraint(i,row)*constraint(i)
			end do
			! Solve
			!call ma41ds(nn*ned+nel,A,-R,dw)
			call ma57ds(nn*ned+nel,A,-R,dw)
			w = w + dw
			! check convergence
			err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
			err2 = sqrt(dot_product(R,R))/(ned*nn+nel)
			write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
		end do
		
		if (nit == maxit) then
			w = w1
			loadfactor = loadfactor - increment
			increment = increment / 2
		else if (nit > 6) then
			increment = increment*(1-adjust)
		else if (nit < 6) then
			increment = increment*(1+adjust)
		end if
	end do
	step = nprint
	call write_results(filepath,w)
	write(*,*) repeat("=", 95) 
	
	deallocate(Fext)
	deallocate(Fint)
	deallocate(R)
	deallocate(w)
	deallocate(w1)
	deallocate(dw)
	deallocate(A)
	
end subroutine statics

subroutine dynamics(filepath)
		
	use read_file
	use integration
	use face
	use shapefunction
	use mixed_material
	use mixed_externalforce
	use symmetric_solver
	use mixed_output
	use mixed_internalforce
	use mixed_tangentstiffness
	use mass
	
	implicit none
	
	integer :: i, nit, row, col, j, k
	real(8), allocatable, dimension(:) :: Fext, Fint, F, R, un, un1, vn, vn1, an, an1, w, dw
	real(8), allocatable, dimension(:,:) ::  A, M
	real(8) :: loadfactor, increment, err1, err2, gamma, beta
	real(8), dimension(size(bc1, 2)) :: constraint
	real(8), dimension(size(bc1, 2), nn*ned) :: der_constraint
	character(80), intent(in) :: filepath
	
	allocate(Fext(nn*ned+nel))
	allocate(Fint(nn*ned+nel))
	allocate(F(nn*ned+nel))
	allocate(R(nn*ned+nel))
	allocate(w(nn*ned+nel))
	allocate(dw(nn*ned+nel))
	allocate(un(nn*ned+nel))
	allocate(un1(nn*ned+nel))
	allocate(vn(nn*ned+nel))
	allocate(vn1(nn*ned+nel))
	allocate(an(nn*ned+nel))
	allocate(an1(nn*ned+nel))
	allocate(A(nn*ned+nel,nn*ned+nel))
	allocate(M(nn*ned+nel,nn*ned+nel))
	
	! initialization
	w = 0.
	un = 0.
	un1 = 0.
	vn = 0.
	vn1 = 0.
	an = 0.
	an1 = 0.
	constraint = 0.
	der_constraint = 0.
	gamma = 0.5
	beta = 0.25
	nprint = 20
	nsteps = 100
	dt = 0.1
	loadfactor = 0.
	increment = firststep
	step = 0
	call write_results(filepath,w)
	call mass_matrix(M)
	
	! If the external load is traction, then the external force doesn't change
	if (size(bc2, 1) /= 3) then
		call force_traction(Fext)
	end if
	
	F = Fext
	do i = 1, size(bc1,2)
		row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
		do col = 1, nn*ned
			M(row, col) = 0.
			M(col, row) = 0.
		end do
		M(row, row) = 1.
		F(row) = 0.
	end do
	call ma57ds(nn*ned+nel,M,F,an)
	
	do step = 1, nsteps
		w = un
		err1 = 1.
		err2 = 1.
		nit = 0
		write(*,'(A, A, i5, 5x, A, e12.4, A)') repeat('=', 30), 'Step', step, 'Time', step*dt, repeat('=', 36)
		do while (((err1 > tol) .or. (err2 > tol)) .and. (nit < maxit))
			nit = nit + 1
			an1 = (w - (un + dt*vn + 0.5*dt**2*(1 - 2*beta)*an))/(beta*dt**2)
			vn1 = vn + (1 - gamma)*dt*an + gamma*dt*an1
			Fint = force_internal(w)
			! If the external load is pressure, update the external force
			if (size(bc2, 1) == 3) then
				Fext = force_pressure(w)
			end if
			F = Fext - Fint - damp*vn1
			A = tangent_internal(w)
			R = matmul(M, an1) - F
			! penalty
			do i = 1, size(bc1,2)
				row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
				constraint(i) = w(row) - bc1(3,i)
				der_constraint(i,row) = 1. 
			end do
			do i = 1, size(bc1,2)
				row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
				A(row,row) = A(row,row) + penalty*der_constraint(i,row)*der_constraint(i,row)
				R(row) = R(row) + penalty*der_constraint(i,row)*constraint(i)
			end do
			A = M/beta*dt**2 + A
			do i = 1, nn*ned
				A(i,i) = A(i,i) + damp*gamma*dt
			end do
			! Solve
			!call ma41ds(nn*ned+nel,A,-R,dw)
			call ma57ds(nn*ned+nel,A,-R,dw)
			w = w + dw
			! check convergence
			err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
			err2 = sqrt(dot_product(R,R))/(ned*nn+nel)
		end do
		write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
		vn = vn1
		un = w
		an = an1
		
		if (MOD(step, nprint) == 0) then
			call write_results(filepath, w)
		end if
	end do
	
	write(*,*) repeat("=", 95) 
	
	deallocate(Fext)
	deallocate(Fint)
	deallocate(F)
	deallocate(R)
	deallocate(w)
	deallocate(dw)
	deallocate(un)
	deallocate(un1)
	deallocate(vn)
	deallocate(vn1)
	deallocate(an)
	deallocate(an1)
	deallocate(A)
	deallocate(M)
	
end subroutine dynamics

subroutine timestamp ( )

! *****************************************************************************
!
!  TIMESTAMP prints the current YMDHMS date as a time stamp.
! 
!  Author:
!
!    John Burkardt
! *****************************************************************************

  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp