function h = scatterps(lat,lon,varargin)
% scatterps works just like Matlab's scatter function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71 S).
% 
%  scatterps(lat,lon,S,C) 
% 
% is equivalent to 
% 
%  [x,y] = ll2ps(lat,lon); 
%  scatter(x,y,S,C)
% 
%% Syntax
% 
%  scatterps(lat,lon,S,C)
%  scatterps(lat,lon)
%  scatterps(lat,lon,S)
%  scatterps(...,markertype)
%  scatterps(...,'filled')
%  scatterps(...,'PropertyName',PropertyValue)
%  scatterps(...,'km')
%  scatterps(...,'meridian',meridian) 
%  h = scatterps(...)
% 
%% Description 
% 
% scatterps(lat,lon,S,C) displays colored circles at the locations specified by 
% the vectors lat and lon (which must be the same size), plotted in Antarctic 
% polar stereographic eastings and northings. S can be a vector the same length 
% as lat and lon or a scalar. If S is a scalar, MATLAB draws all the markers the 
% same size. If S is empty, the default size is used.
% 
% scatterps(lat,lon) draws the markers in the default size and color.
% 
% scatterps(lat,lon,S) draws the markers at the specified sizes (S) with a single 
% color. This type of graph is also known as a bubble plot.
% 
% scatterps(...,markertype) uses the marker type specified instead of 'o' (see 
% LineSpec for a list of marker specifiers).
%
% scatterps(...,'filled') fills the markers. 
% 
% scatterps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
% 
% scatterps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% h = scatterps(...) returns the handle of the scattergroup object created.
% 
%% Example
% 
% gl = load('asaid_gl.mat'); 
% plotps(gl.lat,gl.lon,'red','km')
% 
% % some random data: 
% lat = -80+randn(15,1); 
% lon = -120+8*randn(15,1); 
% z = rand(15,1); 
% 
% scatterps(lat,lon,50,z,'filled','km')
% 
% xlabel('eastings (km)') 
% ylabel('northings (km)')
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

assert(nargin>1,'The scatterps function requires at least two input: lat and lon.') 
assert(isnumeric(lat)==1,'scatterps requires numeric inputs first.') 
assert(isnumeric(lon)==1,'scatterps requires numeric inputs first.') 
assert(max(abs(lat(:)))<=90,'I suspect you have entered silly data into scatterps because some of your latitudes have absolute values exceeding 90 degrees.') 

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

h = scatter(x,y,varargin{:}); 

%% Put things back the way we found them: 

daspect(da)

if ~hld
   hold off
end

%% Clean up: 

if nargout==0
    clear h
end

end

