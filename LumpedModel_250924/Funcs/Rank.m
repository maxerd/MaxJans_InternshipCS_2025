function rnk = Rank(A)
%% rnk = Rank(A)
% Determine the rank of the matrix A, without thresholding
%       Matlab's rank() might reduce the rank due to too much difference in
%       the number's order. I.e. when the smallest numbers become ~1e10
%       times smaller than the largest numbers, the smallest numbers will
%       be considered 0, this causes the Observability matrix to lose rank
%
% See also rank

    valueTreshold = 0;
    rnk = sum((svd(A)>valueTreshold));

end