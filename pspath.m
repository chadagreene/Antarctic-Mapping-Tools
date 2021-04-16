function [out1,out2] = pspath(lat_or_x,lon_or_y,spacing,varargin)
% pspath returns coordinates a path with equal spacing in polar stereographic 
% coordinates. This function might be used to find even spacing for common 
% interpolation points along a satellite ground track. 
% 
%% Syntax
% 
% [lati,loni] = pspath(lat,lon,spacing) 
% [xi,yi] = pspath(x,y,spacing) 
% [...] = pspath(...,'method',InterpolationMethod) 
% 
%% Description 
% 
% [lati,loni] = pspath(lat,lon,spacing) connects the geographic points lat,lon 
% by a path whose points lati,loni are separated by spacing meters. If input
% coordinates are geo coordinates, output coodinates are also geo coordinates. 
% 
% [xi,yi] = pspath(x,y,spacing) connects the polar stereographic points x,y 
% by a path whose points xi,yi are separated by spacing meters. If input
% coordinates are polar stereographic coordinates, output coodinates are also 
% polar stereographic coordinates. 
% 
% [...] = pspath(...,'method',InterpolationMethod) specifies an interpolation 
% method for path creation. Default is 'linear'. 
% 
%% Example
% AMT comes with some sample data.  Consider this sample flight line: 
% 
% D = load('samplegrid.mat'); 
% lat = D.lat(4000:10:5000); 
% lon = D.lon(4000:10:5000); 
% z = D.z(4000:10:5000); 
% plotps(lat,lon,'ko-')
% 
% That flight line is not spaced very evenly.  Perhaps you want more equal spacing--
% let's get a path with a point every 500 m: 
% 
% [lati,loni] = pspath(lat,lon,500); 
% plotps(lati,loni,'rx-')
% 
% Maybe you prefer not to use linear interpolation: 
% [lati_spline,loni_spline] = pspath(lat,lon,500,'method','spline'); 
% plotps(lati_spline,loni_spline,'b^-')
% axis([-1499557   -1495441    -468396    -465326])
% 
%% Citing Antarctic Mapping Tools
% This function was developed for Antarctic Mapping Tools for Matlab (AMT). If AMT is useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
% @article{amt,
%   title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%   author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%   journal={Computers \& Geosciences},
%   year={2017},
%   volume={104},
%   pages={151--157},
%   publisher={Elsevier}, 
%   doi={10.1016/j.cageo.2016.08.003}, 
%   url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
% }
%   
%% Author Info
% This function was created by Chad A. Greene of the University of Texas at Austin
% Institute for Geophysics (UTIG), April 2016. 
% 
% See also: psgrid and pathdistps.  

%% Initial error checks: 

narginchk(3,inf) 
nargoutchk(2,2) 
assert(isvector(lat_or_x)==1,'Input error: input coordinates must be vectors of matching dimensions.') 
assert(isequal(size(lat_or_x),size(lon_or_y))==1,'Input error: dimensions of input coordinates must match.') 
assert(isscalar(spacing)==1,'Input error: spacing must be a scalar.') 
assert(exist('islatlon.m','file')==2,'Error: Cannot find Antarctic Mapping Tools (AMT). Make sure you have an up-to-date version of AMT which can be found on the Mathworks File Exchange site and make sure Matlab can find the functions in AMT.') 

%% Input parsing: 

% Determine input coordinates:
if islatlon(lat_or_x,lon_or_y) 
   geoin = true; 
   [x,y] = ll2ps(lat_or_x,lon_or_y); 
else
   geoin = false; 
   x = lat_or_x; 
   y = lon_or_y; 
end

% Did user specify an interpolation method? 
tmp = strncmpi(varargin,'method',4); 
if any(tmp) 
   method = varargin{find(tmp)+1}; 
else
   method = 'linear'; 
end   

%% Mathematics: 

% Calculate distance along the path given by input coordinates:  
d = pathdistps(x,y); 

% Interpolate xi and yi values individually to common spacing along the path: 
xi = interp1(d,x,0:spacing:d(end),method); 
yi = interp1(d,y,0:spacing:d(end),method); 

%% Package outputs: 

% Columnate xi,yi for consistent behavior: 
xi = xi(:); 
yi = yi(:); 

% Transpose outputs if inputs were row vectors: 
if isrow(lat_or_x) 
   xi = xi'; 
   yi = yi'; 
end

% Convert to geo coordinates if inputs were geo coordinates: 
if geoin
   [out1,out2] = ps2ll(xi,yi); 
else
   out1 = xi; 
   out2 = yi;
end


end