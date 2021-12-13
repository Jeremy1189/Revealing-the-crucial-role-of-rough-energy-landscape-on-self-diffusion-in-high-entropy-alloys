function y= dec2bin_alloy2(x)
[r,c]= size(x);
temp_y= zeros(r,c);
for j= 1:c
    
        index=find(x(:,j)==1);
        temp_y(index,j)=0;
                
        index=find(x(:,j)==2);
        temp_y(index,j)=1;
      
end
y=temp_y;