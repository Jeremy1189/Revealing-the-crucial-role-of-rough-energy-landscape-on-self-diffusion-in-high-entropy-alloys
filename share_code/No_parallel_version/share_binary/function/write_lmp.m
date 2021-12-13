function write_lmp(filename,new_data,other_file_infor)
fid = fopen(filename,'w+');
    if fid==-1
        error(['Error opening ' filename]); 
    end
    L_other_file_infor = length(other_file_infor);
  for i = 1:L_other_file_infor
      fprintf(fid,[other_file_infor{i} '\n']);
  end
  fprintf(fid, '\n');
  fprintf(fid, '%19.0f %2.0f %f %f %f\n', new_data');
  fclose(fid);