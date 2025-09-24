function rnk = obsvRank(A,C)
%% rnk = obsvRank(A,C)
% Determine the rank of the observability matrix using A and C
%       Matlab's rank() might reduce the rank due to too much difference in
%       the number's order. I.e. when the smallest numbers become ~1e10
%       times smaller than the largest numbers, the smallest numbers will
%       be considered 0, this causes the Observability matrix to lose rank
%
% See also obsv, rank

    % states = size(A,1);
    % 
    % O = [C];
    % for i=1:states
    %     O = [O ; O(end-states+1:end,:)*A];
    % end

    O = obsv(A,C);

    valueTreshold = 0;
    rnk = sum((svd(O)>valueTreshold));

end












