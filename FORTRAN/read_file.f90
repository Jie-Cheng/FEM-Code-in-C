module read_file
	
	implicit none
	
	integer :: mode, maxit, nsteps, nprint, step, isbinary, materialtype
	real(8) :: firststep, adjust, tol, dt, damp, penalty
	real(8) :: materialprops(5), gravity(3)
	integer :: nsd, nen, nn, nel, bc_size, load_size, load_type
	integer, allocatable :: connect(:, :), bc_num(:, :), load_num(:, :)
	real(8), allocatable :: coords(:, :), bc_val(:), load_val(:, :)
	integer, allocatable :: share(:)
	
	integer :: no_nonzeros
	integer, allocatable :: col_ind(:), row_ptr(:), row_ind(:)
	real(8), allocatable :: nonzeros(:)
	
	save
	
contains
	subroutine read_input(mode, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, &
		materialtype, materialprops, gravity, isbinary, penalty)
						  
		implicit none
	
		character(100) :: text
		character(1) :: flag = ':'
		integer :: i, j, k, l, ios
		real(8) :: temp(20)
	
		integer, intent(out) :: mode, maxit, nsteps, nprint, isbinary, materialtype
		real(8), intent(out) :: firststep, adjust, tol, dt, damp, materialprops(5), gravity(3), penalty
	
		open(10, file='input.txt')
		i = 0
		ios = 0
		do while (ios == 0)
			read(10, '(a)', iostat = ios) text
			j = index(text, flag)
			l = 0
			if (j /= 0) then
				i = i + 1
				do k = j + 1, len_trim(text)
					if (text(k:k) /= ' ') then
						l = l + 1
					end if
				end do
				read(text(j+1 : j+1+l), *) temp(i)
			end if	
		end do
	
		mode = int(temp(1))
		isbinary = int(temp(2))
		tol = temp(3)
		penalty = temp(4)
		maxit = int(temp(5))
		gravity(:) = temp(6:8)
		
		firststep = temp(9)
		adjust = temp(10)
		
		nsteps = int(temp(11))
		dt = temp(12)
		nprint = int(temp(13))
		damp = temp(14)
		
		materialtype = int(temp(15))
		materialprops(:) = temp(16:20)
	
		close(10)
	end subroutine read_input
	
	subroutine read_mesh(nsd, nn, nel, nen, coords, connect, bc_size, bc_num, bc_val, &
		load_size, load_type, load_num, load_val, share)
		
		implicit none
		
		integer, intent(out) :: nsd, nen, nn, nel, bc_size, load_size, load_type
		integer, allocatable, intent(out) :: connect(:, :), bc_num(:, :), load_num(:, :)
		real(8), allocatable, intent(out) :: coords(:,:), bc_val(:), load_val(:, :)
		integer, allocatable :: share(:)
		integer :: i, j, temp
	
		open(10, file = 'coords.txt')
		read(10, *) nsd, nn
		allocate(coords(nsd, nn)) 
		do i=1, nn
			read(10,*) coords(:, i)
		end do
		close(10)
	
		open(10, file = 'connect.txt')
		read(10, *) nel, nen
		allocate(connect(nen, nel))
		do i=1, nel
			read(10, *) connect(:, i)
		end do
		close(10)
	
		open(10, file = 'bc.txt')
		read(10, *) bc_size
		allocate(bc_num(2, bc_size))
		allocate(bc_val(bc_size))
		do i = 1, bc_size
			read(10, *) bc_num(:, i), bc_val(i)
		end do
		close(10)
	
		open(10, file = 'load.txt')
		read(10, *) load_size, load_type
		load_type = load_type - 2
		allocate(load_num(2, load_size))
		allocate(load_val(load_type, load_size))
		do i = 1, load_size
			read(10, *) load_num(:, i), load_val(:, i)
		end do
		close(10)
		
		allocate(share(nn))
		share = 0
		do i = 1, nel
			do j = 1, nen
				share(connect(j,i)) = share(connect(j,i)) + 1
			end do
		end do
		
	    close(10)
	
	end subroutine read_mesh
	
	subroutine read_CRS(no_nonzeros, col_ind, row_ind, nonzeros, row_ptr)
		! Assuming the input CRS matrix is a full matrix
		implicit none
		integer, intent(inout) :: no_nonzeros
		integer, allocatable, intent(inout) :: col_ind(:), row_ind(:), row_ptr(:)
		real(8), allocatable, intent(inout) :: nonzeros(:)
		
		integer :: no_rows, i, j, k, full_no_nonzeros, row, col
		integer, allocatable :: full_col_ind(:), full_row_ptr(:), full_row_ind(:)
		
		open(unit=42, form='unformatted', access='stream', file='CRS.bin')
		read(42) full_no_nonzeros
	
		allocate(full_col_ind(full_no_nonzeros))
		read(42) full_col_ind
		read(42) no_rows
		allocate(full_row_ptr(no_rows))
		read(42) full_row_ptr
		close(42)
		allocate(full_row_ind(full_no_nonzeros))
		
		if (no_rows /= nn*nsd + nel + 1) then
			write(*, *) 'Wrong length of row_ptr in CRS input data!'
			stop
		end if
		
		k = 0
		do i = 2, no_rows
			! The (i-1)th row has full_row_ptr(i) - full_row_ptr(i-1) entries
			! every entry has a full_row_ind of i-1
			do j = 1, full_row_ptr(i) - full_row_ptr(i-1)
				k = k + 1
				full_row_ind(k) = i - 1
			end do
		end do
		
		if (k /= full_no_nonzeros) then
			write(*, *) 'Wrong length of full_row_ind!'
			stop 
		end if
		
		! Truncate the matrix, only keep the upper triangle
		no_nonzeros = 0
		do k = 1, full_no_nonzeros
			col = full_col_ind(k)
			row = full_row_ind(k)
			if (col >= row) then
				no_nonzeros = no_nonzeros + 1
			end if
		end do
		
		if ((full_no_nonzeros + no_rows - 1)/2 /= no_nonzeros) then
			write(*, *) 'Wrong no_nonzeros!'
			stop
		end if
		
		allocate(col_ind(no_nonzeros))
		allocate(row_ind(no_nonzeros))
		allocate(nonzeros(no_nonzeros))
		allocate(row_ptr(no_rows))
		row_ptr = 0
		
		i = 0
		do k = 1, full_no_nonzeros
			col = full_col_ind(k)
			row = full_row_ind(k)
			if (col >= row) then
				i = i + 1
				col_ind(i) = col
				row_ind(i) = row
			end if
		end do
		
		if (i /= no_nonzeros) then
			write(*, *) 'Wrong length of col_ind and row_ind!'
			stop
		end if
		
		i = 1
		row_ptr(1) = full_row_ptr(1)
		do k = 2, no_nonzeros
			if (row_ind(k) /= row_ind(k-1)) then
				i = i + 1
				row_ptr(i) = k
			end if
		end do
		row_ptr(i+1) = no_nonzeros + 1
		
		if (i + 1 /= no_rows) then
			write(*, *) 'Wrong length of row_ptr!'
			stop
		end if
		
		do i = 1, no_rows - 1
			if (row_ind(row_ptr(i)) /= i .or. row_ind(row_ptr(i)) /= i ) then
				write(*, *) 'Missing diagonal entry!'
				stop
			end if
		end do
		
		deallocate(full_col_ind)
		deallocate(full_row_ind)
		deallocate(full_row_ptr)
	end subroutine read_CRS
	
	! Add val at (row, col) into a SYMMETRIC matrix
	subroutine addValueSymmetric(nonzeros, row, col, val)
		implicit none
		integer, intent(in) :: row, col
		real(8), dimension(:) :: nonzeros
		real(8), intent(in) :: val
		integer :: i, pos
		
		if (col >= row) then
			pos = 0
			do i = row_ptr(row), row_ptr(row+1)
				if (col == col_ind(i)) then
					pos = i
					exit
				end if
			end do
			nonzeros(pos) = nonzeros(pos) + val
		end if
	end subroutine addValueSymmetric
	
	! Set val at (row, col) into a SYMMETRIC matrix
	subroutine setValueSymmetric(nonzeros, row, col, val)
		implicit none
		integer, intent(in) :: row, col
		real(8), dimension(:) :: nonzeros
		real(8), intent(in) :: val
		integer :: i, pos
		
		if (col >= row) then
			pos = 0
			do i = row_ptr(row), row_ptr(row+1)
				if (col == col_ind(i)) then
					pos = i
					exit
				end if
			end do
			nonzeros(pos) = val
		end if
	end subroutine setValueSymmetric
	
	! Get val at (row, col)/(col, row) from a SYMMETRIC matrix
	function getValueSymmetric(nonzeros, row, col)
		implicit none
		integer, intent(in) :: row, col
		real(8), dimension(:) :: nonzeros
		real(8) :: val, getValueSymmetric
		integer :: i, pos, rw, cl
		
		if (col >= row) then
			rw = row
			cl = col
		else
			rw = col
			cl = row
		end if
		
		pos = 0
		do i = row_ptr(row), row_ptr(row+1)
			if (col == col_ind(i)) then
				pos = i
				exit
			end if
		end do
		getValueSymmetric = nonzeros(pos)
	end function getValueSymmetric
	
end module read_file