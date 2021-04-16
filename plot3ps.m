function h = plot3ps(lat,lon,z,varargin)
% plot3ps works just like Matlab's plot3 function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71°S).
% 
%  plot3ps(lat,lon,z) 
% 
% is equivalent to: 
% 
%  [x,y] = ll2ps(lat,lon); 
%  plot3(x,y,z) 
% 
%% Syntax
% 
%  plot3ps(lat,lon,z)
%  plot3ps(...,LineSpec)
%  plot3ps(...,'PropertyName',PropertyValue,...)
%  plot3ps(...,'km')
%  plot3ps(...,'meridian',meridian)
%  h = plot3ps(...)
% 
%% Description 
% 
% plot3ps(lat,lon,z) plots georeferenced data in Antarctic polar stereographic
% eastings and northings. 
% 
% plot3ps(...,LineSpec) specifies line or marker style. 
% 
% plot3ps(...,'PropertyName',PropertyValue,...) specifies any number of
% line or marker properties. 
% 
% plot3ps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
% 
% plot3ps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% h = plot3ps(...) returns a column vector of handles to lineseries objects, 
% one handle per line.
% 
%% Example
% Plot gridded bed elevation data as blue dots: 
% 
% [centerlat,centerlon] = scarloc('mount vinson'); 
% [lat,lon,bed] = bedmap2_data('bed',centerlat,centerlon,50); 
% 
% plot3ps(lat,lon,bed,'b.')
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
% Institute for Geophysics (UTIG), November 2015, for inclusion in the
% Antarctic Mapping Tools package. 
% 
% See also: plot3, plot3m, and surfps. 

%% Input checks: 

assert(nargin>2,'The plot3ps function requires at least three inputs: lat, lon, and z.') 
assert(islatlon(lat,lon)==1,'I suspect you have entered silly data into plot3ps because some of the lats or lons fall outside the normal range of geo coordinates.') 

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

%% Convert units and plot3: 

[x,y] = ll2ps(lat,lon,'meridian',meridian); 

% Convert to kilometers if user requested:
if plotkm
    x = x/1000; 
    y = y/1000; 
end

h = plot3(x,y,z,varargin{:}); 
axis tight
set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) diff(get(gca, 'ZLim'))])
hold on; 

%% Clean up: 

if nargout==0
    clear h
end

end

