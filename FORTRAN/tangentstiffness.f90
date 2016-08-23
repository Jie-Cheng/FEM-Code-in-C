module tangentstiffness
	implicit none
contains
	function cross(a, b)
		real(8), dimension(3) :: cross
		real(8), dimension(3), intent(in) :: a, b
		
		cross(1) = a(2)*b(3) - a(3)*b(2)
		cross(2) = a(3)*b(1) - a(1)*b(3)
		cross(3) = a(1)*b(2) - a(2)*b(1)
	end function cross
	
	subroutine tangent_internal(dofs)
		use read_file
		use shapefunction
		use integration
		use material
		implicit none
		real(8), dimension(nn*nsd+nel), intent(in) :: dofs
		real(8), dimension(nsd,nsd,nsd,nsd) :: C
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(nsd,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(nsd,nsd) :: stress
		real(8), dimension(nen*nsd+1,nen*nsd+1) :: kint
		real(8), dimension(nsd) :: xi, intcoord ! intcoord is the coordinates of the integration points, necessary for anisotropic models
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, Finv, B, eye
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		integer :: ele,a,i,npt,j,row,intpt,l,d,k,col,pos,begin,current
		real(8) :: det, Ja, pressure
		real(8), dimension(nsd) :: work ! for lapack inverse
		integer, dimension(nsd) :: ipiv ! for lapack inverse
		integer :: info, n1 ! for lapack inverse
		
		external DGETRF
		external DGETRI
		n1 = nsd
		
		! square matrix
		do i=1,nsd
			do j=1,nsd
				if(i==j) then
					eye(i,j) = 1.
				else
					eye(i,j) = 0.
				end if
			end do
		end do
		
		! initialize 
		nonzeros = 0.0
				
		! allocate
		npt = int_number(nsd,nen,0)
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		xilist = int_points(nsd,nen,npt)
		weights = int_weights(nsd,nen,npt)
		
		! loop over elements
		do ele=1,nel
			! extract coords of nodes, and dofs for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
				do i=1,nsd
					eledof(i,a) = dofs(nsd*(connect(a,ele)-1)+i)
				end do
			end do
			! extract the pressure
			pressure = dofs(nsd*nn+ele)
			! initialize
			kint = 0.
			! loop over integration points
			do intpt=1, npt
				xi = xilist(:, intpt)
				intcoord = matmul(elecoord, sf(nen,nsd,xi))
				dNdxi = sfder(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi = matmul(elecoord, dNdxi)
				! inverse matrix and determinant
				dxidx = dxdxi
				call DGETRF(n1,n1,dxidx,n1,ipiv,info)
				call DGETRI(n1,dxidx,n1,ipiv,work,n1,info)
				if (nsd == 2) then
					det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				else if (nsd == 3) then
					det = dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) - dxdxi(1,1)*dxdxi(3,2)*dxdxi(2,3) &
						  - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3) + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1) &
						  + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2) - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
				end if
				! compute dNdx
				dNdx = matmul(dNdxi, dxidx)
				! deformation gradient, F_ij = delta_ij + dU_i/dx_j
				F = eye + matmul(eledof ,dNdx)
				if (nsd == 2) then
					Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
				else if (nsd == 3) then
					Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
						  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
						  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
				end if
				! compute dNdy, in which y is the coord. after deformation
				! inverse of F
				Finv = F
				call DGETRF(n1,n1,Finv,n1,ipiv,info)
				call DGETRI(n1,Finv,n1,ipiv,work,n1,info)
				dNdy = matmul(dNdx, Finv)
				! compute the Kirchhoff stress
				call Kirchhoffstress(nsd, intcoord, F, pressure, materialtype, materialprops, stress)
				! compute the material stiffness C_ijkl
				call materialstiffness(nsd, intcoord, F, pressure, materialtype, materialprops, C)
				! compute the element internal force
				do a = 1, nen
					do i = 1, nsd
						do d = 1, nen
							do k = 1, nsd
								row = (a-1)*nsd + i
								col = (d-1)*nsd + k
								do j = 1,nsd
									do l = 1,nsd
										kint(row,col) = kint(row,col) + &
											            (eye(i,k)*stress(j,l)+C(i,j,k,l))*dNdy(a,j)*dNdy(d,l)*weights(intpt)*det;			
									end do
								end do
							end do
						end do
						! displacement + pressure term
						row = (a-1)*nsd + i
						col = nsd*nen + 1
						kint(row,col) = kint(row,col) + dNdy(a,i)*Ja*weights(intpt)*det
						kint(col,row) = kint(row,col)
					end do
				end do
				! pressure + pressure term
				kint(nsd*nen+1, nsd*nen+1) = kint(nsd*nen+1, nsd*nen+1) - 1/materialprops(2)*weights(intpt)*det
			end do
			! scatter the element kint into the global Kglo
			do a = 1, nen
				do i = 1, nsd
					row = nsd*(connect(a, ele) - 1) + i
					do d = 1, nen
						do k = 1, nsd
							col = nsd*(connect(d, ele) - 1) + k
							if (col >= row) then
								begin = row_ptr(row)
								current = begin + col - col_ind(begin)
								nonzeros(current) = nonzeros(current) + kint(nsd*(a-1)+i,nsd*(d-1)+k)
							end if
						end do
					end do
					col = nsd*nn + ele
					if (col >= row) then
						begin = row_ptr(row)
						current = begin + col - col_ind(begin)
						nonzeros(current) = nonzeros(current) + kint(nsd*(a-1)+i,nsd*nen+1)
					end if
				end do
			end do
			row = nn*nsd + ele
			col = nn*nsd + ele
			begin = row_ptr(row)
			current = begin + col - col_ind(begin)
			nonzeros(current) = nonzeros(current) + kint(nsd*nen+1,nsd*nen+1)
		end do
		deallocate(xilist)
		deallocate(weights)
	end subroutine tangent_internal
	
	subroutine tangent_external(dofs, Kglo)
		use read_file, only: nsd,nn,nel,nen,coords,connect, load_size, load_num, load_val
		use shapefunction
		use integration
		use face
		implicit none
		! input argument
		real(8), dimension(nn*nsd+nel), intent(in) :: dofs
		! output
		real(8), dimension(nn*nsd+nel,nn*nsd+nel) :: Kglo
		! Local variables
		real(8), allocatable, dimension(:,:) :: kext ! Element stiffness 
		! nodelist
		! face_coord: Coord. of the nodes in a face
		! face_xi: Local coord. of the integration points in a face
		! dNdxi: Shape function derivatives of a 2d face
		! face_dof: The dofs of the nodes in a face
		real(8), allocatable, dimension(:) :: weights, N
		real(8), allocatable, dimension(:,:) :: face_coord, face_xi, dNdxi, face_dof
		real(8), dimension(nsd - 1) :: xi
		real(8), dimension(nsd,nsd-1) :: dydxi
		real(8) :: external_pressure
		integer, allocatable, dimension(:) :: nodelist
		integer :: i, j, ele, faceid, nfacenodes, a, npt, intpt, b, k, row, col, l
		real(8), dimension(nsd,nsd,nsd) :: epsilon
		real(8), dimension(nsd) :: normal
		
		! initialize
		Kglo = 0.
		nfacenodes = face_nodes_no(nsd,nen)
		npt = int_number(nsd-1, nfacenodes, 0)
		epsilon = 0.
		epsilon(1,2,3) = 1.
		epsilon(3,1,2) = 1. 
		epsilon(2,3,1) = 1.
		epsilon(2,1,3) = -1. 
		epsilon(3,2,1) = -1.
		epsilon(1,3,2) = -1.
		
		allocate(nodelist(nfacenodes))
		allocate(face_coord(nsd,nfacenodes))
		allocate(face_dof(nsd,nfacenodes))
		allocate(face_xi(nsd,npt))
		allocate(weights(npt))
		allocate(dNdxi(nfacenodes,nsd-1))
		allocate(N(nfacenodes))
		allocate(kext(nfacenodes*nsd,nfacenodes*nsd))
		
		do l = 1, load_size
			ele = load_num(1, l)
			faceid = load_num(2, l)
			nodelist(:) = face_nodes(nsd,nen,nfacenodes,faceid)
			do a = 1, nfacenodes
				face_coord(:,a) = coords(:,connect(nodelist(a),ele))
				do i = 1, nsd
					face_dof(i,a) = dofs(nsd*(connect(nodelist(a),ele)-1)+i)
				end do 
			end do
			external_pressure = load_val(1, l)
			! element external tangent stiffness
			kext = 0.
			! set up integration points and weights
			npt = int_number(nsd-1, nfacenodes, 0)
			face_xi = int_points(nsd-1,nfacenodes,npt)
			weights = int_weights(nsd-1,nfacenodes,npt)
			! loop over the integration points
			do intpt=1,npt
				xi(:) = face_xi(:,intpt)
				dNdxi(:,:) = sfder(nfacenodes,nsd-1,xi)
				N(:) = sf(nfacenodes,nsd-1,xi)
				! set up the jacobian matrix
				dydxi(:,:) = matmul(face_coord+face_dof, dNdxi)
				normal = cross(dydxi(:,1), dydxi(:,2))
				do a = 1, nfacenodes
					do i = 1, nsd
						row = nsd*(a-1) + i
						do b = 1, nfacenodes
							do k = 1, nsd
								col = nsd*(b-1) + k
								kext(row, col) = kext(row, col) - (normal(i)*dydxi(k,1)*N(b)*dNdxi(a,1) &
								- normal(k)*dydxi(i,1)*N(b)*dNdxi(a,1) + normal(i)*dydxi(k,2)*N(b)*dNdxi(a,2) - normal(k)*dydxi(i,2)*N(b)*dNdxi(a,2)) &
								*weights(intpt)*external_pressure
							end do
						end do
					end do
				end do
			end do
			! scatter
			do a = 1, nfacenodes
				do i = 1, nsd
					row = (connect(nodelist(a),ele)-1)*nsd + i
					do b = 1, nfacenodes
						do k = 1, nsd
							col = (connect(nodelist(b),ele)-1)*nsd + k
							Kglo(row, col) = Kglo(row,col) + kext(nsd*(a-1)+i, nsd*(b-1)+k)
						end do
					end do
				end do
			end do
			
		end do
		
		deallocate(nodelist)
		deallocate(face_coord)
		deallocate(face_dof)
		deallocate(face_xi)
		deallocate(weights)
		deallocate(dNdxi)
		deallocate(N)
		deallocate(kext)
	end subroutine tangent_external
end module tangentstiffness