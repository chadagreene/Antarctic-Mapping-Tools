function [] = mapzoom(varargin)
%MAPZOOM zooms a map to a specified location and scale.
% 
% This function is part of Chad Greene's Antarctic Mapping Tools and
% requires Matlab's Mapping Toolbox. The Antarctic Mapping Tools package
% can be found at http://www.mathworks.com/matlabcentral/fileexchange/47638.
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
%% Syntax
% 
% mapzoom(lat,lon) 
% mapzoom(x,y)
% mapzoom('SCAR location') 
% mapzoom(...,'mapwidth',widthkm)
% mapzoom(...,widthkm)
% mapzoom(...,'inset')
% mapzoom(...,'inset','location')
% mapzoom(...,'insetsize',sizefraction)
% mapzoom(...,'frame','off')
% mapzoom('default') 
% 
%% Description 
% 
% mapzoom(lat,lon) centers a 500 km by 500 km map about the georeferenced 
% location given by lat, lon. 
% 
% mapzoom(x,y) centers a 500 km by 500 km map about the polar stereographic 
% eastings and northings x and y.  Polar stereographic coordinates are
% automatically determined by the islatlon function. 
% 
% mapzoom('SCAR location') uses scarloc to find the lat,lon corresponding
% to the string 'SCAR location'. 
% 
% mapzoom(...,'mapwidth',widthkm) specifies width of the map in kilometers.
% Default map width is 500 km. 
%
% mapzoom(...,widthkm) is shorthand for 'mapwidth',widthkm.
% 
% mapzoom(...,'inset') places an inset map in the lower right hand corner
% of the map to show geographic context. 
%
% mapzoom(...,'inset','location') specifies location of the inset. Location
% can be 
%           'southeast' or 'se'  lower right corner (default)
%           'northwest' or 'nw'  upper left corner
%           'northeast' or 'ne'  upper right corner
%           'southwest' or 'sw'  lower left corner
%
% mapzoom(...,'insetsize',sizefraction) specifies size of the inset as a
% fraction of the width of the current map. Default sizefraction is 0.25. 
%
% mapzoom(...,'frame','off') removes frame from the inset. 
%
% mapzoom('default') returns a map to default extents.
% 
%% EXAMPLES:
% 
% antmap
% load coast
% patchm(lat,long,[.588 .976 .482])
% mapzoom('antarctic peninsula')
% 
% figure
% bedmap2('patchshelves')
% bedmap2('patchgl')
% mapzoom('Wilkins Ice shelf')
% scarlabel('Wilkins Ice Shelf','fontangle','italic')
% 
% figure
% bedmap2('velocity','caxis',[0 1500])
% mapzoom('recovery glacier','mapwidth',1000,...
%     'inset','southeast','insetsize',.4,...
%     'frame','off')
% 
% figure
% bedmap2('patchshelves')
% bedmap2('patchgl')
% mapzoom(-80,-120,2000)
% 
% figure
% bedmap2 bed
% mapzoom('thwaites glacier','inset','northwest',...
%   'mapwidth',2000)
% scalebar('length',500,'color','white')
% 
%
%% Author info and function history
% 
% September 2013: This function was originally created as bedmap2_zoom.
% 
% August 2014: The function was renamed mapzoom and made more efficient
% and robust. A major improvement to the inset feature has been made possible 
% by Kelly Kearney's plotboxpos function (http://www.mathworks.com/matlabcentral/fileexchange/9615)
% which is included as a subfunction in mapzoom.  
% 
% January 2015: A minor change: Inset map axes are now given a tag 'insetmap'
% so if mapzoom is called multiple times, previous insets are found and
% deleted before placing a new inset map. 
% 
% Created by Chad A. Greene
% Institute for Geophysics
% The University of Texas at Austin
% 
% See also scarloc and mapzoomps. 

%% Check state, initialize a map if needed

assert(license('test','map_toolbox')==1,'Cannot find an active license for the Mapping Toolbox. Matlab''s Mapping Toolbox is required for use of the mapzoom function.  But you can use mapzoomps instead!') 


try % Is a map already initialized? 
    mapinitialized = getm(gca,'MapProjection');
catch
    antmap; 
end

% Get viewing angle of map: 
[az,el] = view; 

%% Parse inputs: 

inset=false; 
centeredmethod=false; % refers to whether the map is centered on a location 

% Check for inset declaration: 
tmp = strcmpi(varargin,'inset');
if any(tmp) 
    inset = true; 
    insetsize = 0.25; 
    frameon = true; 

    try
        tmpi = find(tmp)+1; 
        if strcmpi(varargin{tmpi},'southwest')||strcmpi(varargin{tmpi},'northwest')||...
            strcmpi(varargin{tmpi},'southeast')||strcmpi(varargin{tmpi},'northeast')||...
            strcmpi(varargin{tmpi},'sw')||strcmpi(varargin{tmpi},'nw')||...
            strcmpi(varargin{tmpi},'se')||strcmpi(varargin{tmpi},'ne')
            location = varargin{tmpi};
            tmp(tmpi)=1; 
        end
    catch
        location = 'southeast';
    end
    varargin = varargin(~tmp); 
end


% check for frame declaration: 
tmp =  strcmpi(varargin,'box')|strcmpi(varargin,'frame');
if any(tmp)
    if strcmpi(varargin{find(tmp)+1},'off')||strcmpi(varargin{find(tmp)+1},'none');
        frameon = false; 
    end
    tmp(find(tmp)+1) = 1; 
    varargin = varargin(~tmp); 
end

nk = 1:length(varargin); 
% check for inset size declaration: 
for k = 1:length(varargin)
    if strcmpi(varargin{k},'size')||strcmpi(varargin{k},'insetsize')
        insetsize = varargin{k+1}; 
        nk(k:k+1)=[];
    end
end
varargin = varargin(nk);


mapwidthkm = 500; % sets default map width to 500 km
for k = 1:length(varargin)
    if strcmpi(varargin{k},'mapwidth')||strcmpi(varargin{k},'mapwidthkm')||...
            strcmpi(varargin{k},'width')
        mapwidthkm = varargin{k+1}; 
    end
end


if length(varargin)>2 
    mapwidthkm = varargin{end};
elseif ischar(varargin{1}) && ~ischar(varargin{end})
    mapwidthkm = varargin{end};
end


% Find and delete old inset if one exists:
try
    hinsetmap = findobj(gcf,'tag','insetmap'); 
    delete(hinsetmap); 
end


if ischar(varargin{1})

switch lower(varargin{1})  
    case {'unzoom','default','continent','antarctica'}
        set(gca,'xlim',[-0.4802 0.4802],'ylim',[-0.4802 0.4802]);
                
    otherwise
        [centerlat,centerlon] = scarloc(varargin{1}); 
        centeredmethod=true; 
        
end

else
    if islatlon(varargin{1},varargin{2})
        centerlat = varargin{1}; 
        centerlon = varargin{2}; 
    else
        [centerlat,centerlon] = ps2ll(varargin{1},varargin{2}); 
    end
    centeredmethod=true; 
end

if centeredmethod
    [lat1,lon1]=minvtran(0,0);
    [lat2,lon2]=minvtran(.1,0);
    kmpermapunit = 10*distance(lat1,lon1,lat2,lon2,6378.137); 

    [centerx,centery] = mfwdtran(centerlat,centerlon); 
    unitstoedge = .5*mapwidthkm/kmpermapunit;

    set(gca,'xlim',[centerx-unitstoedge centerx+unitstoedge],...
        'ylim',[centery-unitstoedge centery+unitstoedge])
end

if inset
    load AMTdata.mat % contains crude grounding line and ice shelf outlines, approximated and downsampled significantly from Bedmap2
    
    gcah = gca; 
    gp = plotboxpos(gca);  
    
    insetwidth = insetsize*gp(3); 
    insetheight = insetsize*gp(4); 
    
    switch lower(location)
        case {'southwest','sw'}
            insetx = gp(1);
            insety = gp(2);   
            
        case {'northeast','ne'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'northwest','nw'}
            insetx = gp(1); 
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'southeast','se'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2);            
            
        otherwise 
            error('Unrecognized inset location.')
    end
    

    axes('position',[insetx insety insetwidth insetheight],'tag','insetmap');
  
        antmap('oceancolor','white')
        for k = 1:30
            patchm(slat{k},slon{k},[.8 .8 .8])
        end

        for k = 1:30
            patchm(glat{k},glon{k},[.5 .5 .5])

        end

    if ~frameon
        antmap('frame','off')
    end
    
    plot([centerx-unitstoedge centerx-unitstoedge centerx+unitstoedge...
        centerx+unitstoedge centerx-unitstoedge],...
        [centery-unitstoedge centery+unitstoedge centery+unitstoedge...
        centery-unitstoedge centery-unitstoedge],'r-','linewidth',1);
    
    view([az,el]); % Sets inset rotation to match main map rotation
    axes(gcah); % returns to original axes
    uistack(gcah,'down');
end




%% Kelly Kearney's plotboxpos function: 

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);
        