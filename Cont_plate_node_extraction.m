n_i=275;
n_col = 19;
n_row =51;
node_mat=zeros(51,19);
for i=1:n_row
    for j=1:n_col
       node_mat(i,j)=275+(j-1)*1+(i-1)*147;
    end
end