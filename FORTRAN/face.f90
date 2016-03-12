module face
	implicit none
	
contains
	
	function face_number(nsd, nen)
		! Return the number of faces in an element
		implicit none
		
		integer, intent(in) :: nsd, nen
		integer :: face_number
		
		if (nsd == 2) then
			if (nen == 3) then
				face_number = 3
			else if (nen == 4) then
				face_number = 4
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				face_number = 4
			else if (nen == 8) then
				face_number = 6
			end if
		end if
	end function face_number
	
	function face_nodes_no(nsd, nen)
		! Return the number of nodes on a face
		implicit none
		
		integer, intent(in) :: nsd, nen
		integer :: face_nodes_no
		
		if (nsd == 2) then
			face_nodes_no = 2
		else if (nsd == 3) then
			if (nen == 4) then
				face_nodes_no = 3
			else if (nen == 8) then
				face_nodes_no = 4
			end if
		end if
	end function face_nodes_no
	
	function face_nodes(nsd, nen, nodes, face)
		! Return the ids of the nodes on a face
		implicit none
		
		integer, intent(in) :: nsd,nen,nodes,face
		integer, dimension(nodes) :: face_nodes
		integer, dimension(3) :: i3
		integer, dimension(4) :: i4
		
		i3 = [2,3,1]
		i4 = [2,3,4,1]
		if (nsd == 2) then
			if (nen == 3) then
				face_nodes = [face, i3(face)]
			else if (nen == 4) then
				face_nodes = [face, i4(face)]
			end if
		else if (nsd == 3) then
			if (nen == 4) then
				if (face == 1) then
					face_nodes = [1,2,3]
				else if (face == 2) then
					face_nodes = [1,4,2]
				else if (face == 3) then
					face_nodes = [2,4,3]	
				else if (face == 4) then
					face_nodes = [3,4,1]
				end if
			else if (nen == 8) then
				if (face == 1) then
					face_nodes = [1,2,3,4]
				else if (face == 2) then
					face_nodes = [5,8,7,6]
				else if (face == 3) then
					face_nodes = [1,5,6,2]
				else if (face == 4) then
					face_nodes = [2,6,7,3]
				else if (face == 5) then
					face_nodes = [3,7,8,4]
				else if (face == 6) then
					face_nodes = [4,8,5,1]
				end if
			end if
		end if
	end function face_nodes
	
end module face
		
		