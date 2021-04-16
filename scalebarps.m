function [h_scalebar,h_text] = scalebarps(varargin) 
% scalebarps places a graphical reference scale in polar stereographic units. 
% 
%% Syntax 
% 
%  scalebarps
%  scalebarps(...,'length',ScaleLength)
%  scalebarps(...,'location',LocationOnMap)
%  scalebarps(...,'units',LengthUnits)
%  scalebarps(...,'TextProperty',TextValue)
%  [h_scalebar,h_scalebartext] = scalebarps(...)
% 
%% Description 
% 
% scalebarps places a graphical reference scale at the lower left-hand
% corner of a map. Length of the scale is determined automatically based on 
% current extents of the map. 
% 
% scalebarps(...,'length',ScaleLength) specifies the length of the scalebar.
% ScaleLength is in the specified units.  If no units are specified, kilometers 
% are assumed. 
% 
% scalebarps(...,'location',LocationOnMap) specifies location of the scalebar
% on the map. Location can be 
%           'southwest' or 'sw' (lower left) {default} 
%           'northwest' or 'nw' (upper left) 
%           'northeast' or 'ne' (upper right)
%           'southeast' or 'se' (lower right) 
% 
% scalebarps(...,'units',LengthUnits) specifies length units as 
%           'kilometers' or 'km' (default) 
%           'miles' or 'mi' 
%           'nautical miles' or 'nmi' 
%           'meters' or 'm' 
%           'feet' or 'ft' 
%           'yards' or 'yd'
% 
% scalebarps(...,'TextProperty',TextValue) specifies properties of the scalebar's 
% text label. 
%
% [h_scalebar,h_text]  = scalebarps(...) returns a handle for the scalebar
% and a handle for the text label. 
% 
%% Examples: 
% Start by plotting the extents of the continent as patch objects with bedmap2: 
% 
% bedmap2('patch coast','xy')
% 
% % Zoom in on the Amundsen sea: 
% axis(1000*[-2300 -600 -1600 -100])
% 
% SIMPLEST-CASE EXAMPLE: 
% % Place a scalebar using automatic sizing: 
% scalebarps
% 
% MOST COMPLICATED-CASE EXAMPLE: 
% % Place a red scalebar with large bold italic text exactly 100 nautical
% % miles long in the lower right-hand corner of the map: 
% scalebarps('color','red',... % make it red
%     'fontsize',20,...        % embiggen the text 
%     'fontweight','bold',...  % put some meat on that text's bones 
%     'fontangle','italic',... % lean in like MSNBC
%     'length',100,...         % 100 units long. What units exactly? 
%     'units','nmi',...        % oh, nautical miles.  Those units. 
%     'location','se')         % place it in the lower right hand corner. 
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
% This function and supporting documentation were written by Chad A. Greene
% of the University of Texas at Austin's Institute for Geophysics (UTIG), 
% November 2015. 
% http://www.chadagreene.com
% 
% See also plotps and scarlabel.

%% Set defaults:

autolength = true; % scalebarps length automatically calculated as ~1/5 mapwidth 
location = 'southwest'; % default location
units = 'km';
color = 'k'; 

%% Error checks: 

if license('test','map_toolbox')==1
    assert(ismap(gca)==0,'You have attempted to call scalebarps for map axes.  I think you want the plain old scalebar function instead.') 
end

%% Parse inputs: 

% Check for user-declared location: 
tmp = strncmpi(varargin,'loc',3); 
if any(tmp)
    location = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isnumeric(location)==0,'scalebarps location must be a string.')
end

% Check for user-declared length: 
tmp = strncmpi(varargin,'length',3); 
if any(tmp)
    autolength = false; 
    lngth_newunits = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isnumeric(lngth_newunits)==1,'Length must be numeric.')
    assert(isscalar(lngth_newunits)==1,'Length must be a scalar.')
end

% Check for user-declared units: 
tmp = strncmpi(varargin,'units',4); 
if any(tmp)
    units = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(isnumeric(units)==0,'Input error: Units must be a string.')
end

% Check for user-declared color: 
tmp = strncmpi(varargin,'color',3); 
if any(tmp) 
    color = varargin{find(tmp)+1}; 
end

%% Get size scale: 
       
switch lower(units)
    case {'km','kilometers','kilometer'}
        scalefactor = 1000; % converts plotted meters to kilometers; 
        
    case {'nmi','nautical miles','nm'} 
        scalefactor = 1852; 
        
    case {'feet','foot','ft'} 
        scalefactor = 0.3048; 
        
    case {'mi','mile','miles'} 
        scalefactor = 1609.35; 
        
    case {'m','meter','meters'} 
        scalefactor = 1; 
        
    case {'yard','yards','yd'} 
        scalefactor = 0.9144; 
        
    otherwise 
        error(['Unrecognized unit type ',units,'.'])
end
        
xl = get(gca,'xlim');
yl = get(gca,'ylim'); 

% Set automatic scale size as one approximately fifth of the current map width:   
if autolength
    
    % width of the current map: 
    mapwidth = abs(diff(xl));
    
    % Reasonable values of auto scalebarps length:
    lengths = [0.001 0.01 0.1 0.25 0.5 1 2 5 10 20 25 50 100 200 250 500 1000 2000 2500 5000 10000 20000 25000 50000]; 
 
    % Set auto scalebarps length as one fifth of the current map width: 
    lngth_newunits = interp1(lengths,lengths,(mapwidth/5)/scalefactor,'nearest');
end

lngth_m = lngth_newunits*scalefactor; 
    
%% Get coordinates for placement: 

switch lower(location)
    case {'southwest','sw'}
        x1 = .05*(xl(2)-xl(1))+xl(1); 
        x2 = x1+lngth_m; 
        y1 = .05*(yl(2)-yl(1))+yl(1); 
        y2 = y1; 

    case {'southeast','se'}
        x1 = .95*(xl(2)-xl(1))+xl(1); 
        x2 = x1-lngth_m; 
        y1 = .05*(yl(2)-yl(1))+yl(1); 
        y2 = y1; 

    case {'northwest','nw'}
        x1 = .05*(xl(2)-xl(1))+xl(1); 
        x2 = x1+lngth_m; 
        y1 = .93*(yl(2)-yl(1))+yl(1); 
        y2 = y1;         

    case {'northeast','ne'}
        x1 = .95*(xl(2)-xl(1))+xl(1); 
        x2 = x1-lngth_m; 
        y1 = .93*(yl(2)-yl(1))+yl(1); 
        y2 = y1;  

    otherwise
        error('Invalid location string for scalebarps.')
end

%% Plot scalebar and text: 

h_scalebar=line([x1 x2],[y1 y2],'color',color,'linewidth',2);

h_text = text(mean([x1 x2]),mean([y1 y2]),[num2str(lngth_newunits),' ',units],...
    'horizontalalignment','center',...
    'verticalalignment','bottom','fontsize',10,varargin{:});
 
%% Clean up

% Return the title handle only if it is desired: 
if nargout==0
    clear h_scalebar; 
end

end



