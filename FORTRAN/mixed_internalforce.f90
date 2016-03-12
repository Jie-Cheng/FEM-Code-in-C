module mixed_internalforce	
	implicit none
contains
	function force_internal(dofs)
		use read_file, only: nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		use mixed_material
		implicit none
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		real(8), dimension(nn*ned+nel) :: force_internal
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(ned,nsd) :: stress
		real(8), dimension(nen*ned+1) :: fele
		real(8), dimension(nsd) :: xi, intcoord ! intcoord is the coordinates of the integration points, necessary for anisotropic models
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, Finv, B, eye
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		integer :: ele,a,i,npt,j,row,intpt
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
		
		! initialize force_internal
		force_internal = 0.
		
		! allocate
		npt = int_number(nsd,nen,0)
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		
		! loop over elements
		do ele=1,nel
			! extract coords of nodes, and dofs for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
				do i=1,ned
					eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i)
				end do
			end do
			! extract the element pressure
			pressure = dofs(ned*nn+ele)
			! compute the internal force
			! initialize
			fele = 0.
			! fully integration
			! set up integration points and weights
			npt = int_number(nsd,nen,0)
			xilist = int_points(nsd,nen,npt)
			weights = int_weights(nsd,nen,npt)
			! loop over integration points
			do intpt=1,npt
				xi = xilist(:,intpt)
				intcoord = matmul(elecoord,sf(nen,nsd,xi))
				dNdxi = sfder(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi = matmul(elecoord,dNdxi)
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
				dNdx = matmul(dNdxi,dxidx)
				! deformation gradient, F_ij = delta_ij + dU_i/dx_j
				F = eye + matmul(eledof,dNdx)
				! left Cauchy-Green tensor, B = FF^T and Ja = det(F)
				B = matmul(F,transpose(F))
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
				dNdy = matmul(dNdx,Finv)
				! compute the Kirchhoff stress
				stress = Kirchhoffstress(nsd,ned,intcoord,F,pressure,materialprops)
				! compute the element internal force
				do a=1,nen
					do i=1,nsd
						row=(a-1)*ned+i
						do j=1,nsd
							fele(row) = fele(row) + stress(i,j)*dNdy(a,j)*weights(intpt)*det
						end do
					end do
				end do
				fele(nen*ned+1) = fele(nen*ned+1) + ((Ja-1)-pressure/materialprops(4))*weights(intpt)*det
			end do
			! scatter the element internal force into the global internal force
			do a=1,nen
				do i=1,ned
					row = ned*(connect(a,ele)-1) + i;
					force_internal(row) = force_internal(row) + fele(ned*(a-1)+i);
				end do
			end do
			! scatter the element pressure
			row = nn*ned + ele
			force_internal(row) = force_internal(row) + fele(ned*nen+1)
		end do
		
		deallocate(xilist)
		deallocate(weights)
		
	end function force_internal	
end module mixed_internalforce