module symmetric_solver
	implicit none
	
contains
	subroutine ma57ds(nonzeros, ndofs, R, dw)
		use read_file, only: row_ind, col_ind, no_nonzeros
		implicit none
		integer, intent(in) :: ndofs
		real(8), dimension(no_nonzeros), intent(in) :: nonzeros
		real(8), dimension(ndofs), intent(in) :: R
		real(8), dimension(ndofs), intent(out) :: dw
		! MA57 parameters
		integer :: n, ne, job, lrhs, lwork, nrhs, lkeep, lifact, lfact
		integer, allocatable :: keep(:), iwork(:), ifact(:)
		real(8), allocatable :: fact(:), work(:)
		real(8), dimension(5) :: cntl
		integer, dimension(20) :: icntl
		integer, dimension(40) :: info
		real(8), dimension(20) :: rinfo
		! trivial parameters
		integer :: i, j, k 
		
		external MA57ID, MA57AD, MA57BD, MA57CD
	
		! Set default values for control parameters
		call MA57ID(cntl,icntl)
		n = ndofs
		ne = no_nonzeros
		
		! Set printing option
		icntl(5) = 1
		! Analyse sparsity pattern
		lkeep = 7*n + ne + max(n, ne) + 42
		allocate(keep(lkeep))
		allocate(iwork(5*n))
		call MA57AD(n, ne, row_ind, col_ind, lkeep, keep, iwork, icntl, info, rinfo)
		! Factorize matrix
		lfact = 2*info(9)
		allocate(fact(lfact))
		lifact = 2*info(10)
		allocate(ifact(lifact))
		call MA57BD(n, ne, nonzeros, fact, lfact, ifact, lifact, lkeep, keep, iwork, icntl, cntl, info, rinfo)
		! Solve the equations
		job = 1
		lrhs = n
		nrhs = 1
		lwork = n*nrhs
		allocate(work(lwork))
		call MA57CD(job, n, fact, lfact, ifact, lifact, 1, R, lrhs, work, lwork, iwork, icntl, info)
		dw = R
		deallocate(keep)
		deallocate(iwork)
		deallocate(ifact)
		deallocate(fact)
		deallocate(work)
	end subroutine ma57ds
end module symmetric_solver