function [ L_bf ] = cal_Lbf( f_MHz, d_km )
%   Free-space basic transmission loss
%   radio frequency f is in megahertz
%   distance d is in kilometers
%   f_MHz can be an array

L_bf = 32.45 + 20*log10(f_MHz)+20*log10(d_km);


end

