module shapefunction
	
	implicit none
	! this module contains shape functions and the derivatives of shape functions
	
contains
	
	function sf(nen, nsd, xi)
		implicit none
		
		integer, intent(in) :: nen, nsd
		real(8), dimension(nsd), intent(in) :: xi
		real(8), dimension(nen) :: sf
		
		if (nsd == 1) then
			if (nen == 2) then
				sf(1) = 0.5*(1.+xi(1));
				sf(2) = 0.5*(1.-xi(1));
			end if
		else if (nsd == 2) then
			if (nen == 3) then
				sf(1) = xi(1);
				sf(2) = xi(2);
				sf(3) = 1.-xi(1)-xi(2); 
			else if (nen == 4) then
				sf(1) = 0.25*(1.-xi(1))*(1.-xi(2));
				sf(2) = 0.25*(1.+xi(1))*(1.-xi(2));
	            sf(3) = 0.25*(1.+xi(1))*(1.+xi(2));
		        sf(4) = 0.25*(1.-xi(1))*(1.+xi(2));
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				sf(1) = xi(1);
				sf(2) = xi(3);
				sf(3) = xi(2);
			    sf(4) = 1.-xi(1)-xi(2)-xi(3);
			else if (nen == 8) then
				sf(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
				sf(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
		        sf(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
		        sf(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
		        sf(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
		        sf(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
		        sf(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
		        sf(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
			end if
		end if
	end function sf
	
	function sfder(nen, nsd, xi)
		implicit none
		
		integer, intent(in) :: nen, nsd
		real(8), dimension(nsd), intent(in) :: xi
		real(8), dimension(nen,nsd) :: sfder
		
		if (nsd == 1) then
			if (nen == 2) then
				sfder = reshape([0.5, -0.5],shape(sfder))
			end if
		else if (nsd == 2) then
			if (nen == 3) then
				sfder = reshape([1.,0.,-1.,0.,1.,-1.],shape(sfder))
			else if (nen == 4) then
				sfder(1,1) = -0.25*(1.-xi(2));
				sfder(1,2) = -0.25*(1.-xi(1));
		        sfder(2,1) = 0.25*(1.-xi(2));
		        sfder(2,2) = -0.25*(1.+xi(1));
		        sfder(3,1) = 0.25*(1.+xi(2));
		        sfder(3,2) = 0.25*(1.+xi(1));
		        sfder(4,1) = -0.25*(1.+xi(2));
		        sfder(4,2) = 0.25*(1.-xi(1));
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				sfder(1,1) = 1.;
				sfder(2,3) = 1.;
		        sfder(3,2) = 1.;
		        sfder(4,1) = -1.;
		        sfder(4,2) = -1.;
		        sfder(4,3) = -1.;
			else if (nen == 8) then
				sfder(1,1) = -(1.-xi(2))*(1.-xi(3))/8.;
		        sfder(1,2) = -(1.-xi(1))*(1.-xi(3))/8.;
		        sfder(1,3) = -(1.-xi(1))*(1.-xi(2))/8.;
		        sfder(2,1) = (1.-xi(2))*(1.-xi(3))/8.;
		        sfder(2,2) = -(1.+xi(1))*(1.-xi(3))/8.;
		        sfder(2,3) = -(1.+xi(1))*(1.-xi(2))/8.;
		        sfder(3,1) = (1.+xi(2))*(1.-xi(3))/8.;
		        sfder(3,2) = (1.+xi(1))*(1.-xi(3))/8.;
		        sfder(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
		        sfder(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
		        sfder(4,2) = (1.-xi(1))*(1.-xi(3))/8.;
		        sfder(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
		        sfder(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
		        sfder(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
		        sfder(5,3) = (1.-xi(1))*(1.-xi(2))/8.;
		        sfder(6,1) = (1.-xi(2))*(1.+xi(3))/8.;
		        sfder(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
		        sfder(6,3) = (1.+xi(1))*(1.-xi(2))/8.;
		        sfder(7,1) = (1.+xi(2))*(1.+xi(3))/8.;
		        sfder(7,2) = (1.+xi(1))*(1.+xi(3))/8.;
		        sfder(7,3) = (1.+xi(1))*(1.+xi(2))/8.;
		        sfder(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
		        sfder(8,2) = (1.-xi(1))*(1.+xi(3))/8.;
		        sfder(8,3) = (1.-xi(1))*(1.+xi(2))/8.;
			end if
		end if
	end function sfder
	
end module shapefunction
		
		
		