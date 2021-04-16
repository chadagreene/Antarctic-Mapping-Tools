function [Pe] = peclet(d,thck,sfz,varargin)
% peclet calculates the Peclet number along a glacier flowline using the 
% formulation by Felikson et al., 2017.
% 
%% Syntax 
% 
%  Pe = peclet(d,thck,sfz)
%  Pe = peclet(...,'m',m)
%  Pe = peclet(...,'CouplingLength',Nthck) 
%  Pe = peclet(...,'endpoints','fill') 
%  
%% Description
% 
% Pe = peclet(d,thck,sfz) returns the Peclet number along a glacier flowline
% where d is distance along the flowline in meters, and thck and sfz are the 
% corresponding thickness and surface elevation in meters. The dimensions of 
% d, thck, and sfz must all match. Tip: You can use pathdistps or pathdistpsn 
% to get the distance along the flowline d in meters. 
% 
% Pe = peclet(...,'m',m) specifies a positive exponent m that relates to basal
% sliding (see Supp. Eq. 4 of Felikson et al., 2017.). By default, m = 1.
% 
% Pe = peclet(...,'CouplingLength',Nthck) specifies a longitudinal 
% coupling length as a multiple of ice thickesses. This is equivalent to l/H
% in Kamb & Echelmeyer's  paper cited below. Important: Nthick is not the same 
% thing as the total window width. The Kamb & Echelmeyer paper describes it 
% in detail, but the "averaging length" is the full width of a boxcar window 
% and is equal to 4*l. In this function, the default value of Nthck is 2.5, 
% which is equivalent to a moving average window width of 10 ice thicknesses. 
% 
%   For guidance on choosing a value of Nthck, Kamb & Echelmeyer state that
%   "l/H ranges from about 1.5 to 10...for temperate valley glaciers, with f 
%   near 0.5 and with longitudinal strain-rates typically of order 0.01-0.05 /yr, 
%   l/H should be in the range from about 1 to 3, whereas for ice sheets ...
%   the expected l/H is in the range from about 4 to 10, distinctly higher
%   than for valley glaciers."
% 
% Pe = peclet(...,'endpoints','fill') sets endpoints to NaN when performing  
% the moving window average on the thickness and surface profiles, in a manner
% equivalent to the 'fill' option in the movmean function. By default, the 
% moving window shrinks at each end, to provide continuous measurements at the
% edges. This option will result in missing data at each end of the profile, as
% well as near each NaN datapoint. 
% 
%% Examples 
% 
% For examples, type 
% 
%  showdemo peclet_documentation 
% 
%% Citing this function 
% The formulas in this function are taken directly from Felikson et al., 
% 2017, so if you use this function, please cite Denis' paper! And at least 
% for accountability's sake, it's probably prudent to cite my Antarctic Mapping 
% Tools paper too. Here are the citations: 
% 
% Felikson, Denis, Timothy C. Bartholomaus, Ginny A. Catania, Niels J. Korsgaard,
% Kurt H. Kjær, Mathieu Morlighem, Brice Noël et al. "Inland thinning on the Greenland 
% ice sheet controlled by outlet glacier geometry." Nature Geoscience 10, no. 5 
% (2017): 366-369. https://doi.org/10.1038/ngeo2934
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 

%% Error checks: 

narginchk(3,Inf)
assert(isequal(size(d),size(thck),size(sfz)),'Dimensions of d, thck, and sfz must all agree.') 
assert(issorted(d,'ascend'),'Input vector d must be monotonic and increasing.') 
assert(exist('ice_profile_smoother.m','file')==2,'Cannot find the ice_profile_smoother.m function. Make sure it is in your filepath.') 

%% Input parsing: 

% Set defaults: 
Nthck = 2.5; % number of ice thicknesses for smoothing window 
endpoints = 'shrink'; % option for moving average window 
m = 1; 

tmp = strcmpi(varargin,'m'); 
if any(tmp)
   m = varargin{find(tmp)+1}; 
   assert(isscalar(m),'Exponent m must be a positive scalar.') 
end

tmp = strncmpi(varargin,'CouplingLength',3); 
if any(tmp)
   Nthck = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'endpoints',3); 
if any(tmp)
   endpoints = varargin{find(tmp)+1}; 
end

%% Moving averages 

if isequal(Nthck,0)
   thck0 = thck; 
   sfz0 = sfz; 
else
   thck0 = ice_profile_smoother(d,thck,thck,'CouplingLength',Nthck,'EndPoints',endpoints); 
   sfz0 = ice_profile_smoother(d,sfz,thck,'CouplingLength',Nthck,'EndPoints',endpoints); 
end

% Along-flow slope: 
d_d = gradient(d); 
alpha0 = gradient(sfz0)./d_d; 

%% Calculate Peclet 

C0 = (m+1) .* (thck0.^m) .* alpha0.^m; % Supplemental Eq. 14 (this is actually C0/Kb)
D0 = m .* thck0.^(m+1) .* alpha0.^(m-1); % Supp Eq. 15 (this is actually D0/Kb)

Pe = d.*(C0 - gradient(D0)./d_d)./D0; % Eq. 4 of main text. (Assumes di=0 @ terminus. Here my di is equivalent to Felikson's l.)

end