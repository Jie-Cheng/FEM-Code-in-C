module symmetric_solver
	implicit none
	
contains
	subroutine ma57ds(ndofs,Ain,R,dw)
		implicit none
		integer, intent(in) :: ndofs
		real(8), dimension(ndofs,ndofs), intent(in) :: Ain
		real(8), dimension(ndofs), intent(in) :: R
		real(8), dimension(ndofs), intent(out) :: dw
		! MA57 parameters
		integer :: n, ne, job, lrhs, lwork, nrhs, lkeep, lifact, lfact
		integer, allocatable :: irn(:), jcn(:), keep(:), iwork(:), ifact(:)
		real(8), allocatable :: rhs(:), fact(:), work(:), a(:)
		real(8), dimension(5) :: cntl
		integer, dimension(20) :: icntl
		integer, dimension(40) :: info
		real(8), dimension(20) :: rinfo
		! trivial parameters
		integer :: i, j, k 
		
		external MA57ID, MA57AD, MA57BD, MA57CD
	
		! Set default values for control parameters
		call MA57ID(cntl,icntl)
		
		! Count the nonzeros
		n = ndofs
		ne = 0
		k = 0
		do j=1,n
			do i=j,n
				if (Ain(i,j) /= 0.) then
					ne = ne + 1
				end if
			end do
		end do
		allocate(irn(ne))
		allocate(jcn(ne))
		allocate(a(ne))
		allocate(rhs(n))
		rhs = R
		do j=1,n
			do i=j,n
				if (Ain(i,j) /= 0.) then
					k = k + 1
					irn(k) = i
					jcn(k) = j
					a(k) = Ain(i,j)
				end if
			end do
		end do
		
		! Set printing option
		icntl(5) = 1
		
		! Analyse sparsity pattern
		lkeep = 7*n+ne+max(n,ne)+42
		allocate(keep(lkeep))
		allocate(iwork(5*n))
		call MA57AD(n,ne,irn,jcn,lkeep,keep,iwork,icntl,info,rinfo)
	
		! Factorize matrix
		lfact = 100*info(9)
		allocate(fact(lfact))
		lifact = 100*info(10)
		allocate(ifact(lifact))
		call MA57BD(n,ne,a,fact,lfact,ifact,lifact,lkeep,keep,iwork,icntl,cntl,info,rinfo)
	
		! Solve the equations
		job = 1
		lrhs = n
		nrhs = 1
		lwork = n*nrhs
		allocate(work(lwork))
		call MA57CD(job,n,fact,lfact,ifact,lifact,1,rhs,lrhs,work,lwork,iwork,icntl,info)
	
		dw = rhs
	
		deallocate(irn)
		deallocate(jcn)
		deallocate(keep)
		deallocate(iwork)
		deallocate(ifact)
		deallocate(rhs)
		deallocate(fact)
		deallocate(work)
		deallocate(a)
	end subroutine ma57ds
end module symmetric_solver