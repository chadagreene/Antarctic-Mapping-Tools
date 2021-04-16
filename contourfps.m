function [C,h] = contourfps(lat,lon,Z,varargin)
% contourfps works just like Matlab's contourf function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71°S).
% For example, 
% 
%  contourfps(lat,lon,Z) 
% 
% is equivalent to 
% 
%  [x,y] = ll2ps(lat,lon); 
%  contourf(x,y,Z) 
% 
%% Syntax
% 
%  contourfps(lat,lon,Z)
%  contourfps(lat,lon,Z,n)
%  contourfps(lat,lon,Z,v)
%  contourfps(...,LineSpec)
%  contourfps(...,'km')
%  contourfps(...,'meridian',meridian)
%  [C,h] = contourfps(...)
% 
%% Description 
% 
% contourfps(lat,lon,Z) draws contours of Z at gridded locations lat,
% lon. 
% 
% contourfps(lat,lon,Z,n) specifies a number of contour levels n if n is a scalar. 
% 
% contourfps(lat,lon,Z,v) draws a contour plot of matrix Z with contour lines at the 
% data values specified in the monotonically increasing vector v. The number of contour 
% levels is equal to length(v). To draw a single contour of level i, use v = [i i]; 
% 
% contourfps(...,LineSpec) draws the contours using the line type and color specified 
% by LineSpec. contour ignores marker symbols.
% 
% contourfps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
% 
% contourfps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% [C,h] = contourfps(...) returns contourfps matrix C and handle h of the contour object created.
% 
%% Examples
% 
% % Load some test data: 
% [lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); 
% 
% % Simplest case: 
% contourfps(lat,lon,sfz) 
% 
% % In kilometer coordinates: 
% contourfps(lat,lon,sfz,'km') 
% 
% % Specify levels: 
% contourfps(lat,lon,sfz,0:200:5000,'km') 
% cb = colorbar; 
% ylabel(cb,'surface elevation (m)')
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
% This function was written by Chad A. Greene of the University of Texas
% Institute for Geophysics (UTIG), July 2015, for inclusion in the
% Antarctic Mapping Tools package. 
%
% See also: scatter, scatterm, ll2ps, pcolorps, and plotps. 

%% Input checks: 

assert(nargin>2,'The contourfps function requires at least three input: lat and lon and Z.') 
assert(isnumeric(lat)==1,'contourfps requires numeric inputs first.') 
assert(isnumeric(lon)==1,'contourfps requires numeric inputs first.') 
assert(isnumeric(Z)==1,'contourfps requires numeric Z input.') 
assert(isvector(lat)==0,'lat must be 2D grid.') 
assert(isvector(lon)==0,'lon must be 2D grid.') 
assert(isvector(Z)==0,'Z must be 2D grid.') 
assert(max(abs(lat(:)))<=90,'I suspect you have entered silly data into contourfps because some of your latitudes have absolute values exceeding 90 degrees.') 

%% Parse inputs

plotkm = false; % by default, plot in meters 
meridian = 0;   % top of the map is Fimbul Ice Shelf

% Has user requested plotting in kilometers? 
if nargin > 2
   tmp = strcmpi(varargin,'km'); 
   if any(tmp)
      plotkm = true; 
      varargin = varargin(~tmp); 
   end
   
   tmp = strcmpi(varargin,'meridian'); 
   if any(tmp)
      meridian = varargin{find(tmp)+1}; 
      assert(isscalar(meridian)==1,'Error: meridian must be a scalar longitude.') 
      tmp(find(tmp)+1) = true; 
      varargin = varargin(~tmp); 
   end
end

%% Get initial conditions 

da = daspect; 
da = [1 1 da(3)]; 
hld = ishold; 
hold on

%% Convert units and plot: 

[x,y] = ll2ps(lat,lon,'meridian',meridian); 

% Convert to kilometers if user requested:
if plotkm
    x = x/1000; 
    y = y/1000; 
end

[C,h] = contourf(x,y,Z,varargin{:}); 
hold on; 

%% Put things back the way we found them: 

daspect(da)

if ~hld
   hold off
end

%% Clean up: 

if nargout==0
    clear C h
end

end

