function h = textps(lat,lon,string,varargin)
% textps works just like Matlab's text or textm functions, but places georeferenced
% text labels in Antarctic polar stereographic cartesian coordinates (true latitude 71°S).
% 
%  textps(lat,lon,'My text') 
% 
% is equivalent to: 
% 
%  [x,y] = ll2ps(lat,lon); 
%  plot(x,y,'My text') 
% 
%% Syntax
% 
%  textps(lat,lon,'string')
%  textps(...,'PropertyName',PropertyValue,...)
%  textps(...,'km')
%  textps(...,'meridian',meridian)
%  h = textps(...)
% 
%% Description 
% 
% textps(lat,lon,'string') adds the string in quotes to the location specified 
% by the point (lat,lon) polar stereographic eastings and northings. 
% 
% textps(...,'PropertyName',PropertyValue,...) formats the string with any
% name-value pairs of text properties. 
% 
% textps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
% 
% textps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% h = textps(...) returns a column vector of handles to text objects, one handle 
% per object.
% 
%% Example
% Place text labels in East and West Antarctica: 
% 
% antbounds('gl') 
% textps(-80,-105,'West Antarctica') 
% textps(-75,90,'East Antarctica',...
%   'fontangle','italic','color','red')
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
%% Author Info
% This function was written by Chad A. Greene of the University of Texas
% Institute for Geophysics (UTIG), July 2015 for inclusion in the
% Antarctic Mapping Tools package. 
% May 2017: fixed an alignment bug. 
%
% See also: text, textm, ll2ps, and plotps. 

%% Input checks: 

assert(nargin>2,'The textps function requires at least three inputs: lat, lon, and a text string.') 
assert(isnumeric(lat)==1,'textps requires numeric inputs first.') 
assert(isnumeric(lon)==1,'textps requires numeric inputs first.') 
assert(max(abs(lat(:)))<=90,'I suspect you have entered silly data into textps because some of your latitudes have absolute values exceeding 90 degrees.') 

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

h = text(x,y,string,'horiz','center',varargin{:}); 
hold on; 

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

