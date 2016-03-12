module mixed_material
	
	implicit none

contains
	function materialstiffness(nsd,ned,intcoord,F,pressure,materialprops)
		! Compute the material stiffness tensor
		implicit none
	
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: pressure
		real(8), dimension(nsd,nsd), intent(in) :: F
		real(8), dimension(nsd), intent(in) :: intcoord
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd,ned,nsd):: materialstiffness
	
		integer :: i, j, k, l, p
		real(8) :: mu1, mu2, K1, I1, I2, c1, c2, c3, Ja, I4, I6, lambda, kk1, kk2, der14, der24, der16, der26, R, beta
		real(8), dimension(nsd,nsd) :: B, Bbar, BB, eye, C
		real(8), dimension(nsd) :: a0, g0, a, g
	
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
		materialstiffness = 0.
		
		! Change B to Bbar and compute Bbar squared
		if (nsd == 2) then
			Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
		else if (nsd == 3) then
			Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
				  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
				  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
		end if
		B = matmul(F,transpose(F))
		Bbar = B/Ja**(2/3.)
		BB = matmul(Bbar,Bbar)
		! I1 I2 denote I1bar I2bar actually
		I1 = 0.
		do i=1,nsd
			I1 = I1 + Bbar(i,i)
		end do
	    I2 = I1**2
		do i=1,nsd
			do j=1,nsd
				I2 = I2 - Bbar(i,j)**2
			end do 
		end do
		if (nsd == 2) then
			I1 = I1 + 1.
			I2 = I2 - 1.
		end if
		I2 = I2/2.
		
		if (dAbs(materialprops(1)-4) < 1d-4) then ! Yeoh
			c1 = materialprops(2)
			c2 = materialprops(3)
			c3 = materialprops(5)
			do i = 1, ned
				do j = 1, nsd
					do k = 1, ned
						do l = 1, nsd
							materialstiffness(i,j,k,l) = (8*c2 + 24*c3*(I1 - 3))*Bbar(i,j)*Bbar(k,l) &
													   - (4/3.*c1 + (16/3.*I1 - 8)*c2 + 12*(I1-1)*(I1-3)*c3)*(eye(i,j)*Bbar(k,l) + eye(k,l)*Bbar(i,j)) &
													   + (4/9.*I1*c1 + (16/9.*I1**2 - 8/3.*I1)*c2 + 4*I1*(I1-1)*(I1-3)*c3)*eye(i,j)*eye(k,l) &
													   + (2/3.*I1*c1 + 4/3.*I1*(I1-3)*c2 + (2*I1*(I1-3)**2)*c3)*(eye(i,k)*eye(j,l) + eye(i,l)*eye(j,k)) &
													   + (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
						end do
					end do
				end do
			end do
		else if (dAbs(materialprops(1)-5) < 1d-4) then ! HGO
			! Material parameters are hard-coded
			beta = 50*3.14159/180
			R = sqrt((1+(tan(beta))**2)*((intcoord(1))**2 + (intcoord(2))**2))
			!a0 = [0.,1.,0.]
			!g0 = [0.,1.,0.]
			a0 = [-intcoord(2)/R,intcoord(1)/R,tan(beta)/sqrt(1+(tan(beta))**2)]
			g0 = [-intcoord(2)/R,intcoord(1)/R,-tan(beta)/sqrt(1+(tan(beta))**2)]
			kk1 = 0.5158d6
			kk2 = 0.9

			C = matmul(transpose(F),F)
			I4 = dot_product(a0,matmul(C,a0))
			lambda = sqrt(I4)
			a = matmul(F,a0)/lambda
			I4 = I4/Ja**(2/3.)
			I6 = dot_product(g0,matmul(C,g0))
			lambda = sqrt(I6)
			g = matmul(F,g0)/lambda
			I6 = I6/Ja**(2/3.)
			! 1st order derivative of psi w.r.t. I4/I6
			der14 = kk1*(I4-1.)*exp(kk2*(I4-1.)**2)
			der16 = kk1*(I6-1.)*exp(kk2*(I6-1.)**2)
			! 2nd order derivative of psi w.r.t. I4/I6
			der24 = kk1*(1+2*kk2*(I4-1.)**2)*exp(kk2*(I4-1.)**2)
			der26 = kk1*(1+2*kk2*(I6-1.)**2)*exp(kk2*(I6-1.)**2)
			mu1 = materialprops(2)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					do k=1,ned
						do l=1,nsd
							materialstiffness(i,j,k,l) = - 2/3.*(mu1)*(Bbar(k,l)*eye(i,j)+Bbar(i,j)*eye(k,l)) &
										 				 + 2/9.*(mu1*I1)*eye(i,j)*eye(k,l) &
										 				 + 1/3.*(mu1*I1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)) &
										 				 + (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
							! anisotropic part 1
							materialstiffness(i,j,k,l) = materialstiffness(i,j,k,l) + 4*I4**2*der24*a(i)*a(j)*a(k)*a(l) &
													   - 4/3.*(I4*der24+der14)*I4*(a(i)*a(j)*eye(k,l) + eye(i,j)*a(k)*a(l)) &
													   + (4/9.*I4**2*der24 + 4/9.*I4*der14)*eye(i,j)*eye(k,l) &
													   + 2/3.*I4*der14*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))
							! anisotropic part 2
							materialstiffness(i,j,k,l) = materialstiffness(i,j,k,l) + 4*I6**2*der26*g(i)*g(j)*g(k)*g(l) &
													   - 4/3.*(I6*der26+der16)*I6*(g(i)*g(j)*eye(k,l) + eye(i,j)*g(k)*g(l)) &
													   + (4/9.*I6**2*der26 + 4/9.*I6*der16)*eye(i,j)*eye(k,l) &
													   + 2/3.*I6*der16*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))
						end do
					end do
				end do
			end do
		else !neo-Hookean and Mooney-Rivlin
			mu1 = materialprops(2)
			mu2 = materialprops(3)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					do k=1,ned
						do l=1,nsd
							materialstiffness(i,j,k,l) = mu2*(2*Bbar(i,j)*Bbar(k,l)-(Bbar(i,k)*Bbar(j,l)+Bbar(i,l)*Bbar(j,k))) &
										 - 2/3.*(mu1+2*mu2*I1)*(Bbar(k,l)*eye(i,j)+Bbar(i,j)*eye(k,l)) &
										 + 4/3.*mu2*(BB(i,j)*eye(k,l)+BB(k,l)*eye(i,j)) &
										 + 2/9.*(mu1*I1+4*mu2*I2)*eye(i,j)*eye(k,l) &
										 + 1/3.*(mu1*I1+2*mu2*I2)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)) &
										 !+ (Ja*(2*Ja-1)*eye(i,j)*eye(k,l) - Ja*(Ja-1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*K1
										 + (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
						end do
					end do
				end do
			end do
		end if
	
		
	end function materialstiffness
	
	function Kirchhoffstress(nsd,ned,intcoord,F,pressure,materialprops)
		! Compute the Kirchhoff stress
		implicit none
		
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: pressure
		real(8), dimension(nsd,nsd), intent(in) :: F
		real(8), dimension(nsd), intent(in) :: intcoord
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd) :: Kirchhoffstress
		real(8), dimension(nsd,nsd) :: B, Bbar, BB, eye, C
		integer :: i, j, k
		real(8) :: mu1, mu2, K1, I1, I2, c1, c2, c3, Ja, I4, I6, lambda, kk1, kk2, R, beta
		real(8), dimension(nsd) :: a0, g0, a, g
		
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
		! Change B to Bbar and compute Bbar squared
		if (nsd == 2) then
			Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
		else if (nsd == 3) then
			Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
				  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
				  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
		end if
		B = matmul(F,transpose(F))
		Bbar = B/Ja**(2/3.)
		BB = matmul(Bbar,Bbar)
		! I1 I2 denote I1bar I2bar actually
		I1 = 0.
		do i=1,nsd
			I1 = I1 + Bbar(i,i)
		end do
	    I2 = I1**2
		do i=1,nsd
			do j=1,nsd
				I2 = I2 - Bbar(i,j)**2
			end do 
		end do
		if (nsd == 2) then
			I1 = I1 + 1.
			I2 = I2 - 1.
		end if
		I2 = I2/2.
		Kirchhoffstress = 0.
		
		if (dAbs(materialprops(1)-4) < 1d-4) then ! Yeoh
			c1 = materialprops(2)
			c2 = materialprops(3)
			c3 = materialprops(5)
			do i = 1, ned
				do j = 1, nsd
					Kirchhoffstress(i,j) = (2*c1 + 4*c2*(I1 - 3) + 6*c3*(I1-3)**2)*(Bbar(i,j) - 1/3.*I1*eye(i,j)) + Ja*pressure*eye(i,j)
				end do
			end do
		else if (dAbs(materialprops(1)-5) < 1d-4) then ! HGO
			! Material parameters are hard-coded
			beta = 50*3.14159/180
			R = sqrt((1+(tan(beta))**2)*((intcoord(1))**2 + (intcoord(2))**2))
			!a0 = [0.,1.,0.]
			!g0 = [0.,1.,0.]
			a0 = [-intcoord(2)/R,intcoord(1)/R,tan(beta)/sqrt(1+(tan(beta))**2)]
			g0 = [-intcoord(2)/R,intcoord(1)/R,-tan(beta)/sqrt(1+(tan(beta))**2)]
			kk1 = 0.5158d6
			kk2 = 0.9
			
			C = matmul(transpose(F),F)
			I4 = dot_product(a0,matmul(C,a0))
			lambda = sqrt(I4)
			a = matmul(F,a0)/lambda
			I4 = I4/Ja**(2/3.)
			I6 = dot_product(g0,matmul(C,g0))
			lambda = sqrt(I6)
			g = matmul(F,g0)/lambda
			I6 = I6/Ja**(2/3.)
			
			mu1 = materialprops(2)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					Kirchhoffstress(i,j) = -1/3.*mu1*I1*eye(i,j) + mu1*Bbar(i,j) + Ja*pressure*eye(i,j);
					! anisotropic part
					Kirchhoffstress(i,j) = Kirchhoffstress(i,j) + 2*kk1*( (I4-1)*exp(kk2*(I4-1)**2)*(a(i)*a(j) - 1/3.*eye(i,j))*I4 ) &
					                     + 2*kk1*( (I6-1)*exp(kk2*(I6-1)**2)*(g(i)*g(j) - 1/3.*eye(i,j))*I6 )
				end do
			end do
		else !neo-Hookean and Mooney-Rivlin
			mu1 = materialprops(2)
			mu2 = materialprops(3)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					Kirchhoffstress(i,j) = -1/3.*(mu1*I1+2*mu2*I2)*eye(i,j) - mu2*BB(i,j) + (mu1+mu2*I1)*Bbar(i,j) + Ja*pressure*eye(i,j);
				end do
			end do
		end if
	end function Kirchhoffstress	
	
end module mixed_material