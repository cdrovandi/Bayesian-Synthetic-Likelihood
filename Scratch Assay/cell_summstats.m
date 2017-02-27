function summ_stat = cell_summstats(Y,Yinit)
% cell_summstats is used for calculating the summary statistics of data
% The summary statistics are the Hamming distances and total cells at the final time period.
%
% INPUT:
% Y - the binary matrices for cell presence at all observed times after time zero.
% Yinit - the initial binary matrix for cell presence
%
% OUTPUT:
% summ_stat - the summary statistics for the given data

num_obs = size(Y,3);
summ_stat = zeros(num_obs+1,1); 

% Hamming distances between cell locations across time
summ_stat(1) = sum(sum(abs(Yinit-Y(:,:,1)))); 
for i = 2:num_obs
    summ_stat(i) = sum(sum(abs(Y(:,:,i-1)-Y(:,:,i))));
end

% Total number of cells in the final time period
summ_stat(end) = sum(sum(Y(:,:,end)));

end

