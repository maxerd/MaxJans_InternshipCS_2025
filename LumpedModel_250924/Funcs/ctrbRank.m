function rnk = ctrbRank(A,B)
%% rnk = ctrbRank(A,B)
% Determine the rank of the observability matrix using A and B
%       Matlab's rank() might reduce the rank due to too much difference in
%       the number's order. I.e. when the smallest numbers become ~1e10
%       times smaller than the largest numbers, the smallest numbers will
%       be considered 0, this causes the Observability matrix to lose rank
%
% See also ctrb, rank

    C = ctrb(A,B);

    valueTreshold = 0;
    rnk = sum((svd(C)>valueTreshold));

end












