function [vx,vy] = uv2vxvy(lat_or_x,lon_or_y,u,v) 
% uv2vxvy transforms georeferenced (zonal and meridional) vector
% components to cartesian (polar stereographic) coordinate components. 
% 
%% Syntax
% 
%  [vx,vy] = uv2vxvy(lat,lon,U,V)
%  [vx,vy] = uv2vxvy(x,y,U,V)
% 
%% Description 
% 
% [vx,vy] = uv2vxvy(lat,lon,U,V) transforms zonal U and meridional V
% components of a vector field to cartesian horizontal vx and vertical
% vy components. Inputs lat and lon define locations of each point in
% U and V. 
% 
% [vx,vy] = uv2vxvy(x,y,U,V)  transforms zonal U and meridional V
% components of a vector field to cartesian horizontal vx and vertical
% vy components. Inputs lat and lon define locations of each point in
% vx and vy. Polar stereographic coordinates are automatically determined 
% if any value in the first two inputs of uv2vxvy exceed normal geographic 
% coordinate values. 
% 
%% Example 
% Let's create a wind-like vector field with a zonal maximum at 65 deg S
% and convert to cartesian coordinates: 
% 
% [lat,lon] = psgrid('south pole',8000,50); 
% 
% U = 10*ones(size(lat)).*cosd((lat+65)*3).^3; 
% V = 3*sind(lon); 
% speed = hypot(U,V);
% 
% antmap('northernlimit',-55)
% pcolorm(lat,lon,speed) 
% bedmap2 'patchshelves'
% bedmap2 'patchgl' 
% quivermc(lat,lon,U,V,'density',15) 
% 
% [vx,vy] = uv2vxvy(lat,lon,U,V); 
% 
% [x,y] = ll2ps(lat,lon); 
% 
% figure
% quiver(x,y,vx,vy) 
% 
%% Citing Antarctic Mapping Tools
% If AMT is useful for you, please cite the following paper: 
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
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG). September 2015.  
% http://www.chadagreene.com
% 
% See also vxvy2uv, ll2ps, ps2ll, and psgrid. 

%% Input Checks: 

narginchk(4,4) 
nargoutchk(2,2) 
assert(isnumeric(lat_or_x)==1,'All inputs for uv2vxvy must be numeric.') 
assert(isnumeric(lon_or_y)==1,'All inputs for uv2vxvy must be numeric.') 
assert(isnumeric(u)==1,'All inputs for uv2vxvy must be numeric.') 
assert(isnumeric(v)==1,'All inputs for uv2vxvy must be numeric.') 
assert(isequal(size(lat_or_x),size(lon_or_y),size(u),size(v))==1,'All inputs to uv2vxvy must be of equal dimensions.') 

%% Parse inputs: 

% Determine whether inputs are geo coordinates or polar stereographic meters 
if islatlon(lat_or_x,lon_or_y)
    lon = lon_or_y; % lat is really just a placeholder to make the function a little more intuitive to use. It is not necessary for calculation. 
else
    [~,lon] = ps2ll(lat_or_x,lon_or_y); 
end

%% Perform calculation: 

vx = u.*cosd(lon) + v.*sind(lon); 
vy = - u.*sind(lon) + v.*cosd(lon); 


end
