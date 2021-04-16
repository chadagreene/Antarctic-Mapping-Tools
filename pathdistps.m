function d = pathdistps(lat_or_x,lon_or_y,varargin)
% pathdistps returns the cumulative distance along a path in polar stereographic 
% coordinates (true lat 71 S). 
% 
%% Syntax
% 
%  d = pathdistps(lat,lon)
%  d = pathdistps(x,y)
%  d = pathdistps(...,'km')
%  d = pathdistps(...,'ref',[reflat reflon])
%  d = pathdistps(...,'ref',[refx refy])
%  d = pathdistps(...,'FixDistortion',false) 
% 
%% Description 
% 
% d = pathdistps(lat,lon) returns the cumulative distance d in meters along the path 
% specified by geo coordinates lat,lon. Coordinates must be vectors of equal size. 
%
% d = pathdistps(x,y) returns the cumulative distance d in meters along the path 
% specified by polar stereographic coordinates x,y where x and y are vectors of equal
% size in ps71 meters.  
%
% d = pathdistps(...,'km') simply divides output by 1000 to give distance in kilometers. 
% 
% d = pathdistps(...,'ref',[reflat reflon]) references the output to the track coordinate
% nearest to the location given by a two-element vector [reflat reflon].  This might be 
% useful when analyzing distance along a satellite ground track relative to a point of 
% interest such as a grounding line. 
% 
% d = pathdistps(...,'ref',[refx refy]) references the output as above, but using polar
% stereogprahic (ps71) coordinates. 
% 
% d = pathdistps(...,'FixDistortion',false) calculates distances in pure polar stereographic 
% coordinates without accounting for distortion induced by the projection. In Nov 2018, 
% the default behavior changed from 'FixDistortion',false to 'FixDistortion',true, 
% meaning distances now reflect true distances rather than distances in the polar 
% stereographic projection. 
% 
%% Example
% Use reftrack from the ICESat reference tracks toolbox and clip to the eastern hemisphere: 
% 
% [lat,lon] = reftrack(1304); 
% lat = lat(lon>0); 
% lon = lon(lon>0); 
% 
% % Calculate the total distance in meters: 
% d = pathdistps(lat,lon); 
% 
% % Or calculate total distance in kilometers, referenced to (66.8575 S, 143.5678 E), which 
% % is a point near the ICESat track's intersection with the grounding line: 
% d = pathdistps(lat,lon,'km','ref',[-66.8575 143.5678]);
% 
% plot(d,lon)
% xlabel 'distance relative to the grounding line (km)'
% ylabel 'longitude (deg)' 
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
% at Austin Institute for Geophysics (UTIG), April 2016. 
% http://www.chadagreene.com
% 
% See also: pathdist, pspath, and cumsum. 

%% Initial error checks: 

narginchk(2,inf) 
assert(exist('islatlon.m','file')==2,'Error: Cannot find Antarctic Mapping Tools (AMT). Make sure you have an up-to-date version of AMT which can be found on the Mathworks File Exchange site and make sure Matlab can find the functions in AMT.') 
assert(isvector(lat_or_x)==1,'Input Error: input coordinates must be vectors.') 
assert(isequal(size(lat_or_x),size(lon_or_y))==1,'Input error: Dimensions of input coordinates must agree.') 

%% Set defaults: 

kmout = false; 
ref = false; 

%% Parse inputs: 

% Convert geo coordinates to ps71 if necessary: 
if islatlon(lat_or_x,lon_or_y)
   lat = lat_or_x; 
   [x,y] = ll2ps(lat_or_x,lon_or_y); 
else
   x = lat_or_x; 
   y = lon_or_y; 
   [lat,~] = ps2ll(x,y); % only need lat for distortion 
end

% Does the user want output in kilometers? 
if any(strcmpi(varargin,'km'))
   kmout = true; 
end

% Does the user want the output distances in reference to some location?  
tmp = strcmpi(varargin,'ref'); 
if any(tmp)
   ref = true; 
   refcoord = varargin{find(tmp)+1}; 
   assert(numel(refcoord)==2,'Input Error: Reference coordinate must be a two-element vector.') 
   if islatlon(refcoord(1),refcoord(2))
      [refx,refy] = ll2ps(refcoord(1),refcoord(2)); 
   else
      refx = refcoord(1); 
      refy = refcoord(2); 
   end
end
   

%% Perform mathematics: 

m = psdistortion(lat(2:end)); 

% Cumulative sum of distances: 
if isrow(x)
   d = [0,cumsum(hypot(diff(x)./m,diff(y)./m))]; 
else
   d = [0;cumsum(hypot(diff(x)./m,diff(y)./m))];
end

% Reference to a location: 
if ref
   dist2refpoint = hypot(x-refx,y-refy); 
   [~,mind] = min(dist2refpoint); 
   d = d - d(mind(1)); 
end

% Convert to kilometers if user wants it that way: 
if kmout
   d = d/1000; 
end


end