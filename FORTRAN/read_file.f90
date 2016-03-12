module read_file
	
	implicit none
	
	integer :: mode, maxit, nsteps, nprint, step, isbinary
	real(8) :: firststep, adjust, tol, dt, damp, penalty
	real(8) :: materialprops(5), gravity(3)
	integer :: nsd, nen, nn, nel, ned
	integer, allocatable :: connect(:,:)
	real(8), allocatable :: coords(:,:), bc2(:,:), bc1(:,:)
	integer, allocatable :: share(:)
	
	save
	
contains
	subroutine read_input(unitnum, filename, mode, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, materialprops, &
						  gravity, isbinary, penalty)
		
		! ********************************************************************************
		! Read in the user-specified input parameters
		! unitnum and filename are the unit number and file name of the input file
		! mode: the type of analysis. 0 -> static, 1 -> dynamic, 2 -> residual, 3 -> debug
		! maxit: the maximum number of iterations in one timestep
		! firststep: the step size of the first step
		! adjust: the adjust factor of automatic stepping
		! nsteps: number of steps (used only in dynamics analysis)
		! nprint: interval of printing
		! tol: tolerance of residual
		! dt: timestep
		! damp: damping
		! materialprops: type, mu1, mu2, K1, rho
		! gravity: always a 3d vector
		! isbinary: the format of output file. 1 -> binary, 0 -> ASCII
		! penalty: penalty parameter
		! ********************************************************************************
						  
		implicit none
		
		integer, intent(in) :: unitnum
		character(len = *), intent(in) :: filename
	
		character(100) :: text
		character(1) :: flag = ':'
		integer :: i, j, k, l, ios
		real(8) :: temp(19)
	
		integer, intent(out) :: mode, maxit, nsteps, nprint, isbinary
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
	
		mode=int(temp(1))
		tol=temp(2)
		maxit=int(temp(3))
		firststep=temp(4)
		adjust=temp(5)
		nsteps=int(temp(6))
		dt=temp(7)
		nprint=int(temp(8))
		damp=temp(9)
		materialprops(:)=temp(10:14)
		gravity(:)=temp(15:17)
	    isbinary = temp(18)
		penalty = temp(19)
		
		close(10)
	end subroutine read_input
	
	subroutine read_mesh(nsd, ned, nn, nel, nen, coords, connect, bc1, bc2, share)
		
		! ************************************************************************************************************************
		! nsd: number of dimensions
		! ned: number of nodal degree of freedom. ned = nsd currently.
		! nn: number of nodes
		! nel: number of elements
		! nen: number of nodes in one element
		! coords(nsd, nn): coordinates of nodes
		! connect(nen, nel): connectivity
		! bc1(3, *): boundary conditions. bc1(1, *) -> node number, bc1(2, *) -> dof number, bc1(3, *) -> value
		! bc2(*, *): loads. pressure load -> bc2(3, *), traction load -> bc2(2 + nsd, *). element number, face number, value
		! share(nn): number of elements that share a specific node
		! ************************************************************************************************************************
		implicit none
		
		integer, intent(out) :: nsd, ned, nen, nn, nel
		integer, allocatable, intent(out) :: connect(:,:)
		real(8), allocatable, intent(out) :: coords(:,:), bc2(:,:), bc1(:,:)
		integer, allocatable :: share(:)
		integer :: no_bc1, no_bc2, i, j, temp
	
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
	
		open(10, file='bc.txt')
		read(10, *) no_bc1
		if (no_bc1 /= 0) then
			allocate(bc1(3, no_bc1))
			do i=1, no_bc1
				read(10, *) bc1(:, i)
			end do
		else
			allocate(bc1(3, 1)) ! no boundary condition
			bc1 = -1.0
		end if	
		close(10)
	
		open(10, file = 'load.txt')
		read(10, *) no_bc2, temp
		if (no_bc2 /= 0) then
			allocate(bc2(temp, no_bc2)) ! temp = 3 -> pressure load; temp = 4 or 5 -> traction load
			do i=1, no_bc2
				read(10, *) bc2(:, i)
			end do
		else
			allocate(bc2(2+nsd, 1)) ! no load
			bc2 = -1.0
		end if
		close(10)
		
		ned = nsd ! in the current program, ned is always equval to nsd
		
		allocate(share(nn))
		share = 0
		do i = 1, nel
			do j = 1, nen
				share(connect(j,i)) = share(connect(j,i)) + 1
			end do
		end do
		
	    close(10)
	
	end subroutine read_mesh
	
end module read_file
	
	