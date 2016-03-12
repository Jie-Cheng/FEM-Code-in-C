module mixed_output
	implicit none
	
contains
	subroutine write_results(filepath,dofs)
		use read_file, only: mode, step, nn, ned, nel, isbinary
		implicit none
		character(80) :: filepath
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		if (step == 0) then
			call write_case(filepath)
		end if
		call write_geometry(filepath,dofs)
		call write_displacement(filepath,dofs)
		call write_stress(filepath,dofs)
		
	end subroutine write_results
	subroutine write_case(filepath)
		use read_file, only: nsteps, nprint, dt, ned, nn
		implicit none
		
		character(80) :: filepath, filename
		real(8), dimension(nsteps/nprint+1) :: time
		integer :: i
		
		filename=trim(filepath)//'solid.case'
		open(unit=10,file=trim(filename))
		write(10,'("FORMAT",/)') 
		write(10,'("type:",12x,"ensight gold",/)') 
		write(10,'("GEOMETRY",/)')
		write(10,'("model:",12x,"solid.geo******",12x,"change_coords_only",/)')
		write(10,'("VARIABLE",/)') 
		write(10,'("vector per node:",12x,"displacement",12x,"solid.dis******")')
		write(10,'("tensor symm per node:",12x,"stress",12x,"solid.sig******",/)')
		write(10,'("TIME",/)')
		write(10,'("time set:",12x,i10)') 1
		write(10,'("number of steps:",12x,i10)') nsteps/nprint + 1
		write(10,'("filename start number:",12x,i10)') 0
		write(10,'("filename increment:",12x,i10)') 1
		write(10,'("time values:")')
		do i=1,(nsteps/nprint+1)
			time(i) = (i-1)*nprint*dt
			write(10,'(f12.3)')  time(i)
		end do
		close(10)
	end subroutine write_case
	
	subroutine write_geometry(filepath,dofs)
		use read_file, only: step, nn, nsd,ned, nen, nel, connect, coords, nprint, isbinary
		implicit none
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		character(80) :: filepath, filename, buffer
		integer :: i,j
		integer, dimension(nn) :: nodeid
		integer, dimension(nel) :: eleid
		real(8), dimension(3,nn) :: coords1
		character(6) :: temp
		
		write(temp,'(i6.6)') step/nprint
		filename = trim(filepath)//'solid.geo'//trim(temp)
		
		
		do i=1,nsd
			do j=1,nn
				coords1(i,j) = coords(i,j) + dofs((j-1)*ned+i)
			end do
		end do
		if (nsd==2) then
			do i=1,nn
				coords1(3,i) = 0.
			end do
		end if
		
		if (isbinary == 1) then
			open(unit=11,file=trim(filename),form='unformatted')
			buffer = 'Fortran Binary'
			write(11) buffer
			buffer = 'This is a geometry file'
			write(11) buffer
			buffer = 'Useless discription line'
			write(11) buffer
			buffer = 'node id given'
			write(11) buffer
			buffer = 'element id given'
			write(11) buffer
			buffer = 'part'
			write(11) buffer
			i = 1
			write(11) i
			buffer = 'solid'
			write(11) buffer
			buffer = 'coordinates'
			write(11) buffer
			write(11) nn
			do i = 1, nn
				nodeid(i) = i
			end do
			write(11) nodeid
			write(11) sngl(coords1(1,:))
			write(11) sngl(coords1(2,:))
			write(11) sngl(coords1(3,:))
			if (nsd==2) then
			    if (nen==3) then
					buffer = 'tria3'
			    elseif (nen==4) then
					buffer = 'quad4'
			    end if
			elseif (nsd==3) then
			    if (nen==4) then
					buffer = 'tetra4'
			    elseif (nen==8) then
					buffer = 'hexa8'
			    end if
			end if
			write(11) buffer
			write(11) nel
			do i = 1, nel
				eleid(i) = i
			end do
			write(11) eleid
			write(11) connect
			close(11)
		else
			open(unit=10,file=trim(filename),form='FORMATTED')
			write(10,'(A)') 'Ensight Model Geometry File'
			write(10,'(A)') 'Ensight Model Geometry File'
			write(10,'(A)') 'node id given'
			write(10,'(A)') 'element id given'
			write(10,'(A)') 'part'
			write(10,'(i10)') 1
			write(10,'(A)') 'solid'
			write(10,'(A)') 'coordinates'
			write(10,'(i10)') nn
			do i=1,nn
				write(10,'(i10)') i
			end do
			do i=1,3
				do j=1,nn
					write(10,'(e12.5)') coords1(i,j)
				end do
			end do
			if (nsd==2) then
			    if (nen==3) then
					write(10,'(A)') 'tria3'
			    elseif (nen==4) then
					write(10,'(A)') 'quad4'
			    end if
			elseif (nsd==3) then
			    if (nen==4) then
					write(10,'(A)') 'tetra4'
			    elseif (nen==8) then
					write(10,'(A)') 'hexa8'
			    end if
			end if
			write(10,'(i10)') nel
			do i=1,nel
				write(10,'(i10)') i
			end do
			do i=1,nel
				write(10,'(*(i10))') connect(:,i)
			end do
			close(10)
		end if	
	end subroutine write_geometry
	
	
	subroutine write_displacement(filepath,dofs)
		use read_file, only: step, nsd, ned, nn, coords, nel, nen, connect, nprint, isbinary
		implicit none
		
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		character(80) :: filepath, filename, buffer
		integer :: i,j,row
		character(6) :: temp
		real(8), dimension(3,nn) :: displacement
		
		do i=1,nsd
			do j=1,nn
				row=(j-1)*nsd+i
				displacement(i,j)=dofs(row)
			end do
		end do
		if (nsd==2) then
			do i=1,nn
				displacement(3,i) = 0.
			end do
		end if
		
		write(temp,'(i6.6)') step/nprint
		filename = trim(filepath)//'solid.dis'//trim(temp)
		
		if (isbinary == 1) then
			open(unit=11,file=trim(filename),form='unformatted')
			buffer = 'This is a vector per node file for displacement'
			write(11) buffer
			buffer = 'part'
			write(11) buffer
			i = 1
			write(11) i
			buffer = 'coordinates'
			write(11) buffer
			do i = 1, 3
				write(11) sngl(displacement(i,:))
			end do
			close(11)
		else
			open(unit=10,file=trim(filename),form='FORMATTED')
			write(10,'(A)') 'This is a vector per node file for displacement'
			write(10,'(A)') 'part'
			write(10,'(i10)') 1
			write(10,'(A)') 'coordinates'
			do i=1,3
				do j=1,nn
					write(10,'(e12.5)') displacement(i,j)
				end do
			end do
			close(10)
		end if
	end subroutine write_displacement
		
	
	subroutine write_stress(filepath,dofs)
		use read_file, only: step, nsd, ned, nn, coords, nel, nen, connect, materialprops, share, nprint, isbinary
		use shapefunction
		use integration
		use mixed_material
		
		implicit none
		
		real(8), dimension(nn*ned+nel), intent(in) :: dofs
		character(80) :: filepath, filename, buffer
		character(6) :: temp
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(ned,nsd) :: stress
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, B, eye
		real(8), allocatable, dimension(:,:) :: xilist
		integer :: ele,a,i,npt,j,intpt
		real(8) :: Ja, pressure, avg_vcr
		real(8), dimension(nel) :: vcr, abs_vcr ! volume change rate
		real(8), dimension(nsd) :: xi, intcoord ! intcoord is the coordinates of the integration points, necessary for anisotropic models
		real(8), dimension(nsd) :: work ! for lapack inverse
		integer, dimension(nsd) :: ipiv ! for lapack inverse
		integer :: info, n1 ! for lapack inverse
		
		real(8), dimension(6) :: temp_sigma, sigma
		real(8), dimension(6,nn) :: sum_sigma ! the sigma of a particular node added by elements

		write(temp,'(i6.6)') step/nprint
		filename = trim(filepath)//'solid.sig'//trim(temp)
		
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
		do i=1,6
			do j=1,nn
				sum_sigma(i,j) = 0.
			end do
		end do
		vcr = 0.
		! allocate
		npt = int_number(nsd,nen,0)
		allocate(xilist(nsd,npt))
		! loop over elements
		do ele=1,nel
			! extract coords of nodes, and dofs for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
				do i=1,ned
					eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i)
				end do
			end do
			! extract pressure
			pressure = dofs(nn*ned+ele)
			! fully integration
			! set up integration points and weights
			npt = int_number(nsd,nen,0)
			xilist = int_points(nsd,nen,npt)
			! initialize
			do i=1,6
				temp_sigma(i) = 0.
			end do
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
				! compute dNdx
				dNdx = matmul(dNdxi,dxidx)
				! deformation gradient, F_ij = delta_ij + dU_i/dx_j
				F = eye + matmul(eledof,dNdx)
				! left Cauchy-Green tensor, B = F^TF and Ja = det(F)
				B = matmul(F,transpose(F))
				!B = matmul(F,transpose(F))
				if (nsd == 2) then
					Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
				else if (nsd == 3) then
					Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
						  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
						  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
				end if
				vcr(ele) = vcr(ele) + Ja
				! compute the Kirchhoff stress
				stress = Kirchhoffstress(nsd,ned,intcoord,F,pressure,materialprops)
				! Cauchy stress
				stress = stress/Ja
				! vectorize
				if (nsd==2) then
					sigma = [stress(1,1),stress(2,2),dble(0.),stress(1,2),dble(0.),dble(0.)]
				else
					sigma = [stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3)]
				end if
				! average the stress
				temp_sigma = temp_sigma + sigma
			end do
			vcr(ele) = vcr(ele)/npt
			sigma = temp_sigma/npt ! average sigma in this element
			do a=1,nen
				sum_sigma(:,connect(a,ele)) = sum_sigma(:,connect(a,ele)) + sigma
			end do
		end do
		vcr = vcr - 1.
		abs_vcr = abs(vcr)
		avg_vcr = 0.
		i = maxloc(abs_vcr,1)
		j = minloc(abs_vcr,1)
		do a = 1, nel
			avg_vcr = avg_vcr + vcr(a)
		end do
		avg_vcr = avg_vcr/nel
		write(*,'("The maximum volume change is",e12.5,1x,"in element",i10)'),vcr(i),i
		write(*,'("The minimum volume change is",e12.5,1x,"in element",i10)'),vcr(j),j
		write(*,'("The average volume change is",e12.5)'),avg_vcr
		! sigma per node
		do i=1,nn
			sum_sigma(:,i) = sum_sigma(:,i)/dble(share(i))
		end do
		! write to file
		if (isbinary == 1) then
			open(unit=11,file=trim(filename),form='unformatted')
			buffer = 'This is a symm tensor per node file for stress'
			write(11) buffer
			buffer = 'part'
			write(11) buffer
			i = 1
			write(11) i
			buffer = 'coordinates'
			write(11) buffer
			do i = 1, 6
				write(11) sngl(sum_sigma(i,:))
			end do
			close(11)
		else
			open(unit=10,file=trim(filename),form='FORMATTED')
			write(10,'(A)') 'This is a symm tensor per node file for stress'
			write(10,'(A)') 'part'
			write(10,'(i10)') 1
			write(10,'(A)') 'coordinates'
			do i=1,6
				do j=1,nn
					write(10,'(e12.5)') sum_sigma(i,j)
				end do
			end do
			close(10)	
		end if
		deallocate(xilist)
	end subroutine write_stress
end module mixed_output	