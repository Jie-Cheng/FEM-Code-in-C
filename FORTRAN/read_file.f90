module read_file
	
	implicit none
	
	integer :: mode, maxit, nsteps, nprint, step, isbinary, materialtype
	real(8) :: firststep, adjust, tol, dt, damp, penalty, stretch
	real(8) :: materialprops(5), gravity(3)
	integer :: nsd, nen, nn, nel, bc_size, load_size, load_type
	integer, allocatable :: connect(:, :), bc_num(:, :), load_num(:, :)
	real(8), allocatable :: coords(:, :), bc_val(:), load_val(:, :)
	integer, allocatable :: share(:)
	
	integer :: NE
	integer, allocatable :: IRN(:), JCN(:)
	real(8), allocatable :: NONZEROS(:)
	
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
		real(8) :: temp(21)
	
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
	
	subroutine analyze_pattern(nsd, nn, nel, nen, NE, IRN, JCN, NONZEROS)
		implicit none
		integer, intent(in) :: nsd, nn, nel, nen
		integer, intent(inout) :: NE
		integer, allocatable, intent(inout) :: IRN(:), JCN(:)
		real(8), allocatable, intent(inout) :: NONZEROS(:)
		integer :: i, j, a, b, row, col, ele
		
		NE = 0
		
		do ele = 1, nel
			do a = 1, nen
				do i = 1, nsd
					row = nsd*(connect(a, ele) - 1) + i
					do b = 1, nen
						do j = 1, nsd
							col = nsd*(connect(b, ele) - 1) + j
							NE = NE + 1
						end do
					end do
					col = nsd*nn + ele
					NE = NE + 1
				end do
			end do
			row = nsd*nn + ele
			col = nsd*nn + ele
			NE = NE + 1
		end do
		
		allocate(IRN(NE))
		allocate(JCN(NE))
		allocate(NONZEROS(NE))
		
		NE = 0
		do ele = 1, nel
			do a = 1, nen
				do i = 1, nsd
					row = nsd*(connect(a, ele) - 1) + i
					do b = 1, nen
						do j = 1, nsd
							col = nsd*(connect(b, ele) - 1) + j
							NE = NE + 1
							IRN(NE) = row
							JCN(NE) = col
						end do
					end do
					col = nsd*nn + ele
					NE = NE + 1
					IRN(NE) = row
					JCN(NE) = col
				end do
			end do
			row = nsd*nn + ele
			col = nsd*nn + ele
			NE = NE + 1
			IRN(NE) = row
			JCN(NE) = col
		end do
							
		write(*,'(e12.4)') NE*1.0/((nsd*nn + nel)**2)
	end subroutine analyze_pattern
	
end module read_file