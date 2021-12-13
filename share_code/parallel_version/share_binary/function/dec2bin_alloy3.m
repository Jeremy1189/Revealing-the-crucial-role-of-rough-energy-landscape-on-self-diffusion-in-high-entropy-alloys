function y= dec2bin_alloy3(x)
[r,c]= size(x);
temp_y= zeros(r,2*c);
for j= 1:c
    
        index=find(x(:,j)==1);
        temp_y(index,2*j-1)=0;
         temp_y(index,2*j)=0;
        
        index=find(x(:,j)==2);
        temp_y(index,2*j-1)=0;
        temp_y(index,2*j)=1;
        
        index=find(x(:,j)==3);
         temp_y(index,2*j-1)=1;
        temp_y(index,2*j)=0;
    
end
y=temp_y;