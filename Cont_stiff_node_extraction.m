n_i=2965;
n_col = 19;
n_row =51;
node_mat_stiff=zeros(51,19);
for i=1:n_row
    for j=1:n_col
       node_mat_stiff(i,j)=2965+(j-1)*51+(i-1)*1;
    end
end