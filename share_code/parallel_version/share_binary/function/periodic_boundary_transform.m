function periodic_trans_coords = periodic_boundary_transform(relative_coordinate,lattice_constant,total_atom_num)
supercell_dublicate_size =(total_atom_num/4)^(1/3);% 4 is the number of atoms in a fcc cell 
for num_col=1:L_col
    for num_row = 1:L_row
        if relative_coordinate(num_row,num_col)>lattice_constant*supercell_dublicate_size/2
            relative_coordinate(num_row,num_col)= relative_coordinate(num_row,num_col)-lattice_constant*supercell_dublicate_size;
            
%         elseif relative_coordinate(num_row,num_col)<=-lattice_constant*supercell_dublicate_size/2
        elseif relative_coordinate(num_row,num_col)<-lattice_constant*supercell_dublicate_size/2
            relative_coordinate(num_row,num_col)= relative_coordinate(num_row,num_col)+lattice_constant*supercell_dublicate_size;
        elseif abs(relative_coordinate(num_row,num_col))==lattice_constant*supercell_dublicate_size/2
             if relative_coordinate(num_row,num_col)<0
                 disp([num2str(relative_coordinate(num_row,num_col)),'=-',num2str(lattice_constant*supercell_dublicate_size/2')])
             else
                 disp([num2str(relative_coordinate(num_row,num_col)),'=+',num2str(lattice_constant*supercell_dublicate_size/2')])
             end
        else
        end
    end    
end
periodic_trans_coords = relative_coordinate;