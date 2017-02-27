function Y = simulate_cell(theta,X,rows,cols,num_obs,sim_iters)
% simulate_cell simulates the full set of cell presence matrices
% based on the inputted parameters for cell motility and cell proliferation
%
% INPUT:
% theta - the parameters of the model (cell motility and cell proliferation)
% X - the signed 32 bit integer version of the initial cell presence matrix
% rows - the signed 32 bit integer version for the row locations where a cell is present
% cols - the signed 32 bit integer version for the column locations where a cell is present
% num_obs - the number of time periods after the initial time
% sim_iters - 
%
% OUTPUT:
% Y - Cell presence matrices for all remaining time periods

S = X;
rowsS = rows; %Rows which have non-zero values
colsS = cols; %Columns which have non-zero values

Y = simulate_mex(S,rowsS,colsS,theta(1),theta(2),sim_iters,num_obs,randi(1e6));

end