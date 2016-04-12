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
	integer, allocatable :: irn(:), jcn(:), map(:)
	real(8), allocatable :: nonzeros(:)
	
	save
	
contains
	subroutine read_input(unitnum, filename, mode, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, &
		materialtype, materialprops, gravity, isbinary, penalty)
						  
		implicit none
		
		integer, intent(in) :: unitnum
		character(len = *), intent(in) :: filename
	
		character(100) :: text
		character(1) :: flag = ':'
		integer :: i, j, k, l, ios
		real(8) :: temp(20)
	
		integer, intent(out) :: mode, maxit, nsteps, nprint, isbinary, materialtype
		real(8), intent(out) :: firststep, adjust, tol, dt, damp, materialprops(5), gravity(3), penalty
	
		open(unit = unitnum, file = filename)
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
	
	function position(row, col, ndofs)
		implicit none
		integer, intent(in) :: row, col, ndofs
		integer :: position
		position = (2*ndofs - row + 2)*(row - 1)/2 + col - row + 1
	end function position
	
	subroutine initial_compress(ndofs, no_nonzeros, irn, jcn, map, nonzeros)
		implicit none
		integer, intent(in) :: ndofs
		integer, intent(inout) :: no_nonzeros
		integer, allocatable, intent(inout) :: irn(:), jcn(:), map(:)
		real(8), allocatable, intent(inout) :: nonzeros(:)
		
		integer, allocatable :: flag(:)
		integer :: i, j, a, b, row, col, ele, pos, n, count, maxsize
		
		no_nonzeros = 0
		count = 0		
		maxsize = (ndofs**2 + ndofs)/2
		
		allocate(map( maxsize ))
		allocate(flag( maxsize ))
		map = 0
		flag = 0
		
		! Count the no_nonzeros in the upper triangle
		do ele = 1, nel
			do a = 1, nen
				do i = 1, nsd
					row = nsd*(connect(a, ele) - 1) + i
					do b = 1, nen
						do j = 1, nsd
							col = nsd*(connect(b, ele) - 1) + j
							if (col >= row) then
								pos = position(row, col, ndofs)
								map(pos) = 1
							end if
						end do
					end do
					col = nsd*nn + ele
					if (col >= row) then
						pos = position(row, col, ndofs)
						map(pos) = 1
					end if
				end do
			end do
			row = nsd*nn + ele
			col = nsd*nn + ele
			pos = position(row, col, ndofs)
			map(pos) = 1
		end do
		
		do i = 1, maxsize
			if (map(i) == 1) then
				no_nonzeros = no_nonzeros + 1
			end if
		end do
		
		allocate(nonzeros(no_nonzeros))
		allocate(irn(no_nonzeros))
		allocate(jcn(no_nonzeros))
		
		! Modify map so that map(pos) gives the number in the nonzeros
		! flag is used to indicate whether a pos has been taken before
		do ele = 1, nel
			do a = 1, nen
				do i = 1, nsd
					row = nsd*(connect(a, ele) - 1) + i
					do b = 1, nen
						do j = 1, nsd
							col = nsd*(connect(b, ele) - 1) + j
							if (col >= row) then
								pos = position(row, col, ndofs)
								if (flag(pos) /= 1) then
									count = count + 1
									map(pos) = count
									irn(count) = row
									jcn(count) = col
									flag(pos) = 1
								end if
							end if
						end do
					end do
					col = nsd*nn + ele
					if (col >= row) then
						pos = position(row, col, ndofs)
						if (flag(pos) /= 1) then
							count = count + 1
							map(pos) = count
							irn(count) = row
							jcn(count) = col
							flag(pos) = 1
						end if
					end if
				end do
			end do
			row = nsd*nn + ele
			col = nsd*nn + ele
			pos = position(row, col, ndofs)
			if (flag(pos) /= 1) then
				count = count + 1
				map(pos) = count
				irn(count) = row
				jcn(count) = col
				flag(pos) = 1
			end if
		end do
		
		nonzeros = 0.0
		
		deallocate(flag)
		
	end subroutine initial_compress
	
end module read_file