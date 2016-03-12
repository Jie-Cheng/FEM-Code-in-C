module integration
	
	implicit none
	
contains
	
	function int_number(nsd, nen, reduced)
		!  Return the number of integration points
		implicit none
		
		integer, intent(in) :: nsd, nen, reduced
		integer :: int_number
		
		if (reduced == 1) then
			if (nsd == 1) then
				if (nen == 2) then
					int_number = 2
				end if
			else if (nsd == 2) then
				if (nen == 3) then
					int_number = 1
				else if (nen == 4) then
					int_number = 1
				end if
			else if (nsd == 3) then
				if (nen == 4) then
					int_number = 1
				else if (nen == 8) then
					int_number = 1
				end if
			end if
		else if (reduced == 0) then
			if (nsd == 1) then
				if (nen == 2) then
					int_number = 2
				end if
			else if (nsd == 2) then
				if (nen == 3) then
					int_number = 1
				else if (nen == 4) then
					int_number = 4
				end if
			else if (nsd == 3) then
				if (nen == 4) then
					int_number = 1
				else if (nen == 8) then
					int_number = 8
				end if
			end if
		end if
		
    end function int_number
	
	function int_weights(nsd, nen, npt)
		! Return the weights of integration points
		implicit none
		
		integer, intent(in) :: nsd, nen, npt
		real(8), dimension(npt) :: int_weights
		
		if (nsd == 1) then
			if (npt == 2) then
				int_weights = [1., 1.]
			end if
		else if (nsd == 2) then
			if (nen == 3) then
				if (npt == 1) then
					int_weights = 0.5
				end if
			else if (nen == 4) then
				if (npt == 1) then
					int_weights = 4.
				else if (npt == 4) then
					int_weights = [1.,1.,1.,1.]
				end if
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				if (npt == 1) then
					int_weights = 1/6.
				end if
			else if (nen == 8) then
				if (npt == 1) then
					int_weights = 8.
				else if (npt == 8) then
					int_weights = [1.,1.,1.,1.,1.,1.,1.,1.]
				end if
			end if
		end if
		
	end function int_weights
	
	function int_points(nsd, nen, npt)
		! Return the coordinates of integration points
		implicit none
		
		integer, intent(in) :: nsd, nen, npt
		real(8), dimension(nsd,npt) :: int_points
		real(8) :: temp = 1/sqrt(3.)
		
		if (nsd == 1) then
			if (npt == 2) then
				int_points = reshape([-temp,temp],shape(int_points))
			end if
		else if (nsd == 2) then
			if (nen == 3) then
				if (npt == 1) then
					int_points = reshape([1/3., 1/3.],shape(int_points))
				end if 
			else if (nen == 4) then
				if (npt == 1) then
					int_points = reshape([0., 0.],shape(int_points))
				else if (npt == 4) then
					int_points = reshape([-temp,-temp,temp,-temp,-temp,temp,temp,temp],shape(int_points))
				end if
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				if (npt == 1) then
					int_points = reshape([0.25,0.25,0.25],shape(int_points))
				end if
			else if (nen == 8) then
				if (npt == 1) then
					int_points = reshape([0.,0.,0.],shape(int_points))
				else if (npt == 8) then
					int_points = reshape([-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,-1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1.],shape(int_points))
					int_points = temp*int_points
				end if
			end if
		end if
	end function int_points
			
end module integration
					
					
	