function Nb = biotNumber(h,k,Lc)
%% Calculate the Biot number
% ratio between conductive and convection resistances
%   Nb = Rcond/Rconv
%      = h*Lc/k
%   if <0.1, all material in solid can be considered the same temperature
%   if >0.1, different temperatures in solid
%
% Where Lc = characteristic length
%          volume/surfaceArea
%          Lc = thickness for a plate
%          Lc = thickness/2 for a fin
%          Lc = thickness/4 for a long cylinder
%          Lc = thickness/6 for a sphere
%       h = convection constant
%       k = conduction constant

Nb = h*Lc/k;
end