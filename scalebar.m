function [h] = scalebar(varargin) 
% SCALEBAR places a graphical reference scale on a map. This function was
% designed as a simpler alternative to the built-in scaleruler function. 
% 
% This function requires Matlab's Mapping Toolbox. 
% 
%% Syntax 
% 
%  scalebar
%  scalebar('length',ScaleLength)
%  scalebar('units',LengthUnit)
%  scalebar('location',LocationOnMap)
%  scalebar('orientation',ScalebarOrientation)
%  scalebar('TextProperty',TextValue)
%  scalebar('LineProperty',LineValue)
%  h = scalebar(...)
% 
%% Description 
% 
% scalebar places a graphical reference scale at the lower left-hand
% corner of a map. Length of the scale is determined automatically based on 
% current extents of the map. 
% 
% scalebar('length',ScaleLength) specifies the length of the scalebar. 
% Default ScaleLength is approximately one fifth of the width of the current map.
% 
% scalebar('units',LengthUnit) specifies a length unit. Most common length
% units are supported. Text of the scalebar label matches input length
% unit--i.e., If LengthUnit is 'mi', a scalebar may show 100 mi as its label. If 
% you enter 'miles', the scalebar will show 100 miles as its label. Default 
% LengthUnit is 'km'. 
% 
% scalebar('location',LocationOnMap) specifies location of the scalebar
% on the map. Location can be 
%           'southwest' or 'sw' (lower left) {default} 
%           'northwest' or 'nw' (upper left) 
%           'northeast' or 'ne' (upper right)
%           'southeast' or 'se' (lower right) 
%
% scalebar('orientation',ScalebarOrientation) specifies a 'vertical' or
% 'horizontal' scalebar. Default ScalebarOrientation is 'horizontal'. 
%
% scalebar('TextProperty',TextValue) specifies properties of text. 
%
% scalebar('LineProperty',LineValue) specifies properties of the reference
% scale line. 
%
% h = scalebar(...) returns a handle for the scalebar. 
% 
%% Examples 
% 
% EXAMPLE 1: 
% figure; usamap('texas')
% states = shaperead('usastatelo.shp','UseGeoCoords',true);
% geoshow(states, 'DisplayType', 'polygon')
% scalebar % <--SIMPLEST CASE
% scalebar('length',68,'units','miles','color','b','location','se','fontangle','italic') %<--FANCY CASE 
% 
% EXAMPLE 2: (Requires Antarctic Mapping Tools)
% load coast
% antmap
% patchm(lat,long,[.588 .976 .482])
% scalebar
% 
% EXAMPLE 3: (Requires Bedmap2 Toolbox)
% bedmap2 'patchgl'
% bedmap2('patchshelves','oceancolor',[0.0118 0.4431 0.6118])
% mapzoom 'Mertz Glacier Tongue'
% scarlabel('Mertz Glacier Tongue','fontangle','italic')
% scalebar
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
%% Author Info.  
% This function was created by Chad A. Greene of the University of Texas 
% Institute for Geophysics in 2013. This function was originally designed
% for the Bedmap2 Toolbox for Matlab, but has been slightly updated for 
% inclusion in the Antarctic Mapping Tools package.  Although this function 
% was designed for Antarctic maps, it should work for other maps as well.
% 
% Updated July 2015 to automatically select length and allow any user-defined 
% length unit. 
% 
% See also scalebarps and scaleruler. 

%% Error checks:

assert(license('test','map_toolbox')==1,'Sorry, the scalebar function requires Matlab''s Mapping Toolbox.') 
assert(ismap(gca)==1,'A map must be open to use the scalebar function.') 

%% Set defaults:

lngth = 'auto'; % scalebar length automatically set as 1/5 mapwidth 
location = 'southwest'; % default location
orientation = 'horizontal'; % default orientation
units = 'km'; 

%% Parse inputs: 

% Check for user-declared location: 
tmp = strncmpi(varargin,'loc',3); 
if any(tmp)
    location = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isnumeric(location)==0,'scalebar location must be a string.')
end

% Check for user-declared length (also accept "width" or "scale")
tmp = strncmpi(varargin,'len',3)|strncmpi(varargin,'wid',3)|strcmpi(varargin,'scale'); 
if any(tmp)
    lngth = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isscalar(lngth)==1,'Scalebar Length must be a scalar value in kilometers.')
end

 
% Check for user-declared orientation: 
tmp = strncmpi(varargin,'orient',6); 
if any(tmp)
    orientation = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isnumeric(orientation)==0,'Scalebar orientation can only be vertical or horizontal.')
end

tmp = strncmpi(varargin,'units',4); 
if any(tmp)
    units = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    try
        validateLengthUnit(units); 
    catch
        error('Invalid length unit.') 
    end
end

%% Get size scale: 
        
[lat1,lon1]=minvtran(0,0);
[lat2,lon2]=minvtran(.1,0);
dstpermapunit = 10*pathdist([lat1;lat2],[lon1;lon2],units);
dstpermapunit = dstpermapunit(2); 

xl = get(gca,'xlim');
yl = get(gca,'ylim'); 


% Set automatic scale size as one approximately fifth of the current map width:   
if strcmpi(lngth,'auto') 
    
    % left and right extents of the current map: 
    [leftlat,leftlon]=minvtran(xl(1),mean(yl));
    [rightlat,rightlon]=minvtran(xl(2),mean(yl));
    
    % width of the current map: 
    mapwidth = pathdist([leftlat;rightlat],[leftlon;rightlon],units);
    
    % Reasonable values of auto scalebar length:
    lengths = [0.001 0.01 0.1 0.25 0.5 1 2 5 10 20 25 50 100 200 250 500 1000 2000 2500 5000 10000 20000 25000 50000]; 
 
    % Set auto scalebar length as one fifth of the current map width: 
    lngth = interp1(lengths,lengths,mapwidth(2)/5,'nearest');
end

switch lower(orientation)
    case 'horizontal'
    switch lower(location)
        case {'southwest','sw'}
            x1 = .05*(xl(2)-xl(1))+xl(1); 
            x2 = x1+lngth/dstpermapunit; 
            y1 = .05*(yl(2)-yl(1))+yl(1); 
            y2 = y1; 

        case {'southeast','se'}
            x1 = .95*(xl(2)-xl(1))+xl(1); 
            x2 = x1-lngth/dstpermapunit; 
            y1 = .05*(yl(2)-yl(1))+yl(1); 
            y2 = y1; 

        case {'northwest','nw'}
            x1 = .05*(xl(2)-xl(1))+xl(1); 
            x2 = x1+lngth/dstpermapunit; 
            y1 = .94*(yl(2)-yl(1))+yl(1); 
            y2 = y1;         

        case {'northeast','ne'}
            x1 = .95*(xl(2)-xl(1))+xl(1); 
            x2 = x1-lngth/dstpermapunit; 
            y1 = .94*(yl(2)-yl(1))+yl(1); 
            y2 = y1;  
            
        otherwise
            error('Invalid location string for scalebar.')
    end
    
    case 'vertical'
        switch lower(location)
        case {'southwest','sw'}
            x1 = .05*(xl(2)-xl(1))+xl(1); 
            x2 = x1;
            y1 = .05*(yl(2)-yl(1))+yl(1); 
            y2 = y1 + lngth/dstpermapunit; 

        case {'southeast','se'}
            x1 = .95*(xl(2)-xl(1))+xl(1); 
            x2 = x1; 
            y1 = .05*(yl(2)-yl(1))+yl(1); 
            y2 = y1 + lngth/dstpermapunit; 

        case {'northwest','nw'}
            x1 = .05*(xl(2)-xl(1))+xl(1); 
            x2 = x1; 
            y1 = .95*(yl(2)-yl(1))+yl(1); 
            y2 = y1 - lngth/dstpermapunit;        

        case {'northeast','ne'}
            x1 = .95*(xl(2)-xl(1))+xl(1); 
            x2 = x1; 
            y1 = .95*(yl(2)-yl(1))+yl(1); 
            y2 = y1 - lngth/dstpermapunit; 
            
        otherwise
            error('Invalid location string for scalebar.')
        end
end

% In Jan 2015, I added a z value to the scale bar and the text
% because sometimes the scale bar was getting buried underneath
% transparent layers:
h(1)=line([x1 x2],[y1 y2],9999*[1 1],'color','k','linewidth',2);

h(2) = text(mean([x1 x2]),mean([y1 y2]),9999,[num2str(lngth),' ',units],...
    'horizontalalignment','center',...
    'verticalalignment','bottom');
    
% This is brute-force, but it says let's try to set everything that can
% be set, be them text properties or line properties:
for k = 1:2:length(varargin)
    try
        set(h(1),varargin{k},varargin{k+1})
    end
    try
        set(h(2),varargin{k},varargin{k+1})
    end
end

% Return the title handle only if it is desired: 
if nargout==0
    clear h; 
end


