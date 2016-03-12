module mass
	
	implicit none
	
contains
	subroutine mass_matrix(M)
		use read_file, only:  nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		
		implicit none
		
		real(8), dimension(nn*ned+nel,nn*ned+nel), intent(out) :: M
		real(8), dimension(nen*ned,nen*ned) :: mele
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nsd,nsd) :: dxdxi
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		real(8), dimension(nen) :: N
		integer :: ele,a,b,i,j,npt,k,intpt,row,col
		real(8) :: det, rho
		real(8), dimension(nsd) :: xi
		real(8), dimension(nen,nsd) :: dNdxi 
		
		rho = materialprops(5)
		
		! allocate
		npt = int_number(nsd,nen,0)
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
			
		! initialize
		M = 0.
		
		! loop over elements
		do ele=1,nel
			! extract coords of nodes, and dofs for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
			end do
			! fully integration
			! set up integration points and weights
			xilist = int_points(nsd,nen,npt)
			weights = int_weights(nsd,nen,npt)
			! initialize
			mele = 0.
			! loop over integration points
			do intpt=1, npt
				xi = xilist(:,intpt)
				dNdxi = sfder(nen,nsd,xi)
				N = sf(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi = matmul(elecoord,dNdxi)
				! determinant
				if (nsd == 2) then
					det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				else if (nsd == 3) then
					det = dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) - dxdxi(1,1)*dxdxi(3,2)*dxdxi(2,3) &
						  - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3) + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1) &
						  + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2) - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
				end if
				! compute the element mass matrix
				do a = 1, nen
					do b = 1, nen
						do i = 1, ned
							row = ned*(a-1) + i
							col = ned*(b-1) + i
							mele(row,col) = mele(row,col) + N(a)*N(b)*rho*det*weights(intpt)
						end do
					end do
				end do
			end do
			! scatter
			do a=1,nen
				do i=1,ned
					do b=1,nen
						do k=1,ned
							row = ned*(connect(a,ele)-1)+i
							col = ned*(connect(b,ele)-1)+k
							M(row,col) = M(row,col) + mele(ned*(a-1)+i,ned*(b-1)+k)
						end do
					end do
				end do
			end do
		end do
		
		deallocate(xilist)
		deallocate(weights)
		
	end subroutine mass_matrix
end module