module mixed_externalforce
	implicit none
	
contains
	function cross(a, b)
		real(8), dimension(3) :: cross
		real(8), dimension(3), intent(in) :: a, b
		
		cross(1) = a(2)*b(3) - a(3)*b(2)
		cross(2) = a(3)*b(1) - a(1)*b(3)
		cross(3) = a(1)*b(2) - a(2)*b(1)
	end function cross
	
	function force_pressure(dofs)
		! returns the external force due to external pressure
		use read_file, only: nsd,ned,nn,nel,nen,coords,connect,bc1,bc2, gravity, materialprops
		use face
		use shapefunction
		use integration
		
		implicit none
		
		real(8), dimension(nn*ned+nel) :: force_pressure
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		integer :: j,i, no_bc2, ele, faceid, nfacenodes, a, npt, intpt, row
		integer, allocatable, dimension(:) :: nodelist
		real(8), allocatable, dimension(:,:) :: coord, xilist, dNdxi, eledof
		real(8), allocatable, dimension(:) :: fele, weights, N
		real(8), dimension(nsd) :: traction
		real(8), dimension(nsd - 1) :: xi
		real(8), dimension(nsd,nsd-1) :: dydxi
		real(8) :: external_pressure
		
		no_bc2 = size(bc2,2)
		nfacenodes = face_nodes_no(nsd,nen)
		npt = int_number(nsd-1, nfacenodes, 0)
		
		! allocate
		allocate(nodelist(nfacenodes))
		allocate(coord(nsd,nfacenodes))
		allocate(fele(ned*nfacenodes))
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		allocate(dNdxi(nfacenodes,nsd-1))
		allocate(N(nfacenodes))
		allocate(eledof(nsd,nfacenodes))
		
		! initialize Fext
		force_pressure = 0.
		
		! loop over faces with prescribed pressure
		do j=1,no_bc2
			! extract the coords and traction of that face
			ele = int(bc2(1,j))
			faceid = int(bc2(2,j))
			nfacenodes = face_nodes_no(nsd,nen)
			nodelist(:) = face_nodes(nsd,nen,nfacenodes,faceid)
			do a=1,nfacenodes
				coord(:,a) = coords(:,connect(nodelist(a),ele))
				do i = 1, ned
					eledof(i,a) = dofs(ned*(connect(nodelist(a),ele)-1)+i)
				end do
			end do
			external_pressure = bc2(3,j)
			! compute the force on the face
			fele = 0.
			! set up integration points and weights
			npt = int_number(nsd-1, nfacenodes, 0)
			xilist = int_points(nsd-1,nfacenodes,npt)
			weights = int_weights(nsd-1,nfacenodes,npt)
			! loop over the integration points
			do intpt=1,npt
				xi(:) = xilist(:,intpt)
				dNdxi(:,:) = sfder(nfacenodes,nsd-1,xi)
				N(:) = sf(nfacenodes,nsd-1,xi)
				! set up the jacobian matrix
				dydxi(:,:) = matmul(coord+eledof, dNdxi)
				traction = cross(dydxi(:,1), dydxi(:,2))
				traction = traction*external_pressure
				! compute the force
				do a=1,nfacenodes
					do i=1,nsd
						row = ned*(a-1)+i
						fele(row) = fele(row) + traction(i)*N(a)*weights(intpt)
					end do
				end do
			end do
			! scatter the force into the global external force
			do a=1,nfacenodes
				do i=1,ned
					row = (connect(nodelist(a),ele)-1)*ned + i
					force_pressure(row) = force_pressure(row) + fele((a-1)*ned+i)
				end do
			end do
		end do
		
		deallocate(nodelist)
		deallocate(coord)
		deallocate(fele)
		deallocate(xilist)
		deallocate(weights)
		deallocate(dNdxi)
		deallocate(N)
		deallocate(eledof)
			
	end function force_pressure
	
	subroutine force_traction(Fext)
		! Returns the external force due to traction
		use read_file, only: nsd,ned,nn,nel,nen,coords,connect,bc1,bc2, gravity, materialprops
		use face
		use shapefunction
		use integration
		
		implicit none
		
		real(8), dimension(nn*ned+nel), intent(inout) :: Fext
		
		integer :: j,i, no_bc2, ele, faceid, nfacenodes, a, npt, intpt, row
		integer, allocatable, dimension(:) :: nodelist
		real(8), allocatable, dimension(:,:) :: coord, xilist, dNdxi
		real(8), allocatable, dimension(:) :: f, weights, N
		real(8), dimension(nsd) :: traction, xi
		real(8), dimension(nsd,nsd-1) :: dxdxi
		real(8) :: det, rho
		
		no_bc2 = size(bc2,2)
		nfacenodes = face_nodes_no(nsd,nen)
		npt = int_number(nsd-1, nfacenodes, 0)
		
		! allocate
		allocate(nodelist(nfacenodes))
		allocate(coord(nsd,nfacenodes))
		allocate(f(ned*nfacenodes))
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		allocate(dNdxi(nfacenodes,nsd-1))
		allocate(N(nfacenodes))
		
		! initialize Fext
		Fext = 0.
		
		! loop over faces with prescribed traction
		do j = 1, no_bc2
			! extract the coords and traction of that face
			ele = int(bc2(1,j))
			faceid = int(bc2(2,j))
			nfacenodes = face_nodes_no(nsd,nen)
			nodelist(:) = face_nodes(nsd,nen,nfacenodes,faceid)
			do a=1,nfacenodes
				coord(:,a) = coords(:,connect(nodelist(a),ele))
			end do
			do i=1,nsd
				traction(i) = bc2(i+2,j)
			end do
			! compute the force on the face
			f = 0.
			! set up integration points and weights
			xilist = int_points(nsd-1,nfacenodes,npt)
			weights = int_weights(nsd-1,nfacenodes,npt)
			! loop over the integration points
			do intpt=1,npt
				xi(:) = xilist(:,intpt)
				dNdxi(:,:) = sfder(nfacenodes,nsd-1,xi)
				N(:) = sf(nfacenodes,nsd-1,xi)
				! set up the jacobian matrix
				dxdxi(:,:) = matmul(coord,dNdxi)
				if (nsd == 2) then
					det = sqrt((dxdxi(1,1))**2 + (dxdxi(2,1)**2))
				else if (nsd == 3) then
					det = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))**2 &
                		  + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))**2 &
               		      + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))**2)
				end if
				! compute the force
				do a = 1, nfacenodes
					do i = 1, nsd
						row = ned*(a-1)+i
						f(row) = f(row) + traction(i)*N(a)*weights(intpt)*det
					end do
				end do
			end do
			! scatter the force into the global external force
			do a = 1, nfacenodes
				do i = 1, ned
					row = (connect(nodelist(a),ele)-1)*ned + i
					Fext(row) = Fext(row) + f((a-1)*ned+i)
				end do
			end do
		end do
		
		deallocate(nodelist)
		deallocate(coord)
		deallocate(f)
		deallocate(xilist)
		deallocate(weights)
		deallocate(dNdxi)
		deallocate(N)
			
	end subroutine force_traction
	
	
	subroutine force_body(Fext)
		! returns external force due to gravity
		use read_file, only: nsd,ned,nn,nel,nen,coords,connect,bc1,bc2, gravity, materialprops
		use face
		use shapefunction
		use integration
		
		implicit none
		
		real(8), dimension(nn*ned+nel), intent(inout) :: Fext
		integer :: ele, a, npt, i, intpt, row
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned*nen) :: f
		real(8), dimension(nsd) :: xi
		real(8), dimension(nen,nsd) :: dNdxi
		real(8), dimension(nen) :: N
		real(8), dimension(nsd,nsd) :: dxdxi
		real(8) :: det, rho
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		
		rho = materialprops(5)
		npt = int_number(nsd, nen, 0)	
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		
		! initialize Fext
		Fext = 0.
		
		! loop over elements
		do ele=1,nel
			! extract coords of nodes
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
			end do
			! compute the force on this element
			do i=1,ned*nen
				f(i) = 0.
			end do
			! set up integration points and weights
			xilist = int_points(nsd ,nen,npt)
			weights = int_weights(nsd,nen,npt)
			! loop over integration points
			do intpt=1,npt
				xi = xilist(:,intpt)
				dNdxi = sfder(nen,nsd,xi)
				N = sf(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi(:,:) = matmul(elecoord,dNdxi)
				if (nsd == 2) then
					det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				else if (nsd == 3) then
					det = dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) - dxdxi(1,1)*dxdxi(3,2)*dxdxi(2,3) &
						  - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3) + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1) &
						  + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2) - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
				end if
				! compute the force
				do a = 1, nen
					do i = 1, ned
						row = ned*(a-1) + i
						f(row) = f(row) + gravity(i)*N(a)*weights(intpt)*rho*det
					end do
				end do
			end do
			! scatter
			do a = 1, nen
				do i = 1, ned
					row = ned*(connect(a,ele)-1)+i;
					Fext(row) = Fext(row) + f(ned*(a-1)+i);
				end do
			end do
		end do
		
		deallocate(xilist)
		deallocate(weights)
		
	end subroutine force_body
			
end module mixed_externalforce
					  
				
			