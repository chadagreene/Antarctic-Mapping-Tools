function [u,v] = vxvy2uv(lat_or_x,lon_or_y,vx,vy)
% vxvy2uv transforms polar stereographic vector components to
% georeferenced (zonal and meridional) vector components. 
% 
%% Syntax
% 
%  [U,V] = vxvy2uv(lat,lon,vx,vy)
%  [U,V] = vxvy2uv(x,y,vx,vy)
% 
%% Description 
% 
% [U,V] = vxvy2uv(lat,lon,vx,vy) transforms polar stereographic vector
% components vx, vy referenced to the geographic locations in lat and lon to
% geographic zonal and meridional components. 
% 
% [U,V] = vxvy2uv(x,y,vx,vy) transforms polar stereographic vector
% components vx, vy referenced to the polar stereographic locations in x and y to
% geographic zonal and meridional components. Polar stereographic
% coordinates are automatically determined if any value in the first two
% inputs of vxvy2uv exceed normal geographic coordinate values. 
% 
%% Example: Continental scale motion
% Consider a vector field depicting motion from left to right on a standard polar
% stereographic map projection.  Make its x component 1.5 magnitude everywhere, and we'll
% add a slight downward-dipping y component too, just to make things
% interesting.  To set up the grid, us psgrid to create a 6000 km wide grid at 
% 250 km resolution, centered on the South Pole: 
% 
%   [x,y] = psgrid('south pole',6000,250,'xy'); 
%   vx = 1.5*ones(size(x)); 
%   vy = -0.3*ones(size(x)); 
%   [u,v] = vxvy2uv(x,y,vx,vy); 
% 
% You can now plot quiver(x,y,vx,vy) in cartesian coordinates or 
% do [lat,lon] = ps2ll(x,y) then plot quivermc(lat,lon,u,v) on a map. 
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
% See also uv2vxvy, ll2ps, ps2ll, and psgrid. 

%% Input Checks: 

narginchk(4,4) 
nargoutchk(2,2) 
assert(isnumeric(lat_or_x)==1,'All inputs for vxvy2uv must be numeric.') 
assert(isnumeric(lon_or_y)==1,'All inputs for vxvy2uv must be numeric.') 
assert(isnumeric(vx)==1,'All inputs for vxvy2uv must be numeric.') 
assert(isnumeric(vy)==1,'All inputs for vxvy2uv must be numeric.') 
assert(isequal(size(lat_or_x),size(lon_or_y),size(vx),size(vy))==1,'All inputs to vxvy2uv must be of equal dimensions.') 

%% Parse inputs: 

% Determine whether inputs are geo coordinates or polar stereographic meters 
if islatlon(lat_or_x,lon_or_y)
    lon = lon_or_y; % lat is really just a placeholder to make the function a little more intuitive to use. It is not necessary for calculation. 
else
    [~,lon] = ps2ll(lat_or_x,lon_or_y); 
end

%% Perform coordinate transformations 

u = vx .* cosd(lon) - vy.* sind(lon);
v = vy .* cosd(lon) + vx.* sind(lon);


end