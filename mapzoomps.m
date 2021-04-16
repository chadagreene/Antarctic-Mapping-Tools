function h = mapzoomps(varargin)
% mapzoomps zooms a south polar stereographic map to a specified location and extent.
% This is an adaptation of mapzoom, but does not require Matlab's Mapping Toolbox. 
% Syntax for mapzoomps is similar to mapzoom, but differs slightly for some options.
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
%  mapzoomps
%  mapzoomps(lat,lon) 
%  mapzoomps(x,y)
%  mapzoomps('SCAR location') 
%  mapzoomps(...,'size',mapsizekm)
%  mapzoomps(...,InsetLocation)
%  mapzoomps(...,'insetsize',sizefraction)
%  mapzoomps(...,'frame','off')
%  mapzoomps(...,'km')
%  h = mapzoomps(...) 
% 
%% Description 
% 
% mapzoomps(lat,lon) centers a 500 km wide map about the georeferenced 
% location given by lat, lon. 
% 
% mapzoomps(x,y) centers a 500 km wide map about the polar stereographic 
% eastings and northings x and y.  Polar stereographic coordinates are
% automatically determined by the islatlon function. 
% 
% mapzoomps('SCAR location') uses scarloc to find the coordinates corresponding
% to the string 'SCAR location'. 
% 
% mapzoomps(...,'size',mapsizekm) specifies size of the map in kilometers given 
% mapsizekm, which can be a scalar to create a square map or a two-element array
% for a rectangular map in the form [mapwidthkm mapheightkm], where mapwidthkm and
% mapheightkm are the dimensions of the map in kilometers. 
%
% mapzoomps(...,InsetLocation) creates an inset map at the location InsetLocation,
% which can be 
%           'southeast' or 'se'  lower right corner 
%           'northwest' or 'nw'  upper left corner
%           'northeast' or 'ne'  upper right corner
%           'southwest' or 'sw'  lower left corner
%
% mapzoomps(...,'insetsize',sizefraction) specifies size of the inset as a
% fraction of the width of the current map. Default sizefraction is 0.25. 
%
% mapzoomps(...,'frame','off') removes frame from the inset. 
%
% mapzoomps(...,'km') is for plots in polar stereographic kilometers rather than the default meters.
% 
% h = mapzoomps(...) returns a handle h of inset map axes. 
% 
%% Example 1 
% Initialize a 500 km by 500 km map centered on Pitman Fracture Zone, plot ibcso
% bathymetry, and place a graphical reference scale with scalebarps.  When calling mapzoomps
% specify 'ne' to place an inset map in the upper right-hand corner. 
% 
%   mapzoomps('pitman fracture zone','ne')
%   ibcso('image','xy')
%   scalebarps
% 
%% Example 2 
% Plot Bedmap2 bed elevation, apply relief shading with shadem, and center a 1600 km wide
% by 800 km tall map on the Lambert Glacier/Amery Ice Shelf system at (71°S,69°N).  Place 
% frameless inset map in the upper left hand corner: 
% 
%   bedmap2('bed','xy') 
%   shadem([225 60],3)
%   mapzoomps(-71,69,'size',[1600 800],'nw','frame','off')
%   scalebarps
%   scarlabel({'Lambert Glacier','Cape Child','Landing Bluff'},...
%      'fontangle','italic','fontsize',8)
%
%% Author Info 
% This function and supporting documentation were written by Chad A. Greene of the 
% University of Texas at Austin's Institute for Geophysics (UTIG), February 2016. 
% Feel free to contact me if you have any questions or comments. 
% http://www.chadagreene.com
% 
% See also scarloc, scarlabel, scarclick, mapzoom, and scalebarps. 

%% Set defaults: 

inset = false; 
insetsize = 0.25; 
frameon = true; 
location = 'northeast'; 
usekm = false; 
meridian = 0; 
if nargin==0 
   UseCurrentExtents = true; 
else
   UseCurrentExtents = false; 
   mapsize = [500 500]; % sets default map size to 500 km by 500 km
end

%% Parse inputs: 

% Inset location: 
tmp = strcmpi(varargin,'southwest')|strcmpi(varargin,'northwest')|...
      strcmpi(varargin,'southeast')|strcmpi(varargin,'northeast')|...
      strcmpi(varargin,'sw')|strcmpi(varargin,'nw')|...
      strcmpi(varargin,'se')|strcmpi(varargin,'ne'); 
if any(tmp)
   inset = true; 
   location = varargin{tmp}; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end

% Check for inset size declaration: 
tmp = strcmpi(varargin,'insetsize'); 
if any(tmp) 
   inset = true; 
   insetsize = varargin{find(tmp)+1}; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end
   
% Check for frame declaration: 
tmp = strcmpi(varargin,'frame');
if any(tmp)
   inset = true; 
   if strcmpi(varargin{find(tmp)+1},'off')||strcmpi(varargin{find(tmp)+1},'none');
      frameon = false; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
   end
end

% Map width: 
tmp = strcmpi(varargin,'size')|strcmpi(varargin,'mapsize')|strcmpi(varargin,'mapwidth')|strcmpi(varargin,'mapwidthkm')|strcmpi(varargin,'width'); 
if any(tmp)
   mapsize = varargin{find(tmp)+1}; 
   assert(isnumeric(mapsize)==1,'Map size must be numeric.'); 
   if isscalar(mapsize) 
      mapsize = [mapsize mapsize]; 
   end
   assert(numel(mapsize)==2,'Map size must be a one- or two-element numeric value.') 
end

% Polar stereographic kilometers or meters? 
tmp = strcmpi(varargin,'km'); 
if any(tmp) 
   usekm = true; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end
   
% Center location declaration: 
if ~UseCurrentExtents
   % Determine center location: 
   if ~isnumeric(varargin{1})

      % User has entered a location by string: 
      if strcmpi(varargin{1},'default') 
         set(gca,'xlim',3000*[-1 1],'ylim',3000*[-1 1]);
      else
         [xc,yc] = scarloc(varargin{1},'xy'); 
%          [xc,yc] = ll2ps(tmplat,tmplon,'meridian',meridian); % no, this won't work yet. 
      end
      else

      % User has entered location by coordinates: 
      if islatlon(varargin{1},varargin{2})
         [xc,yc] = ll2ps(varargin{1},varargin{2},'meridian',meridian);
         
      else
         xc = varargin{1}; 
         yc = varargin{2}; 
      end
   end
   
   % Check for inputs of lat,lon limits: 
   try
      if isequal(numel(varargin{1}),numel(varargin{2}),2) & islatlon(varargin{1},varargin{2}) 
         latlim = varargin{1}; 
         lonlim = varargin{2}; 
         meridian = mean(lonlim);

      end
   end
   
end

tmp = strcmpi(varargin,'meridian'); 
if any(tmp)
   meridian = varargin{find(tmp)+1}; 
   assert(isscalar(meridian)==1,'Error: meridian must be a scalar.') 
   warning('specifying meridian might not work yet.') 
end     

%% Set axes of map: 

% Find and delete old inset if one exists:
try
   hinsetmap = findobj(gcf,'tag','insetmap'); 
   delete(hinsetmap); 
end

gcah = gca; % handle of initial plot

if UseCurrentExtents
   ax = axis; 
else
   if meridian==0
      axis equal xy
      ax = [xc-mapsize(1)*500 xc+mapsize(1)*500 yc-mapsize(2)*500 yc+mapsize(2)*500]; 
      axis(ax); 
   end
end

% Define x,y coordinates of axis perimeter: 
if meridian==0
   axx = [ax(1) ax(2) ax(2) ax(1) ax(1)];
   axy = [ax(3) ax(3) ax(4) ax(4) ax(3)]; 
else
   lat = [linspace(latlim(1),latlim(1),200),linspace(latlim(2),latlim(2),200),latlim(1)]; 
   lon = [linspace(lonlim(1),lonlim(2),200),linspace(lonlim(2),lonlim(1),200),lonlim(1)];
   [axx,axy] = ll2ps(lat,lon,'meridian',meridian); 
   
   axis equal xy
   axis([min(axx) max(axx) min(axy) max(axy)])
end
   
%% Place an inset map: 

if inset
    load AMTdata.mat % contains crude grounding line and ice shelf outlines, approximated and downsampled significantly from Bedmap2
      
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
    

   % Create new set of axes for inset map: 
   h = axes('position',[insetx insety insetwidth insetheight],'tag','insetmap');
   hold on
   
      
   % Plot ice shelves as patch objects: 
   for k = 1:30
      [tmpx,tmpy] = ll2ps(slat{k},slon{k},'meridian',meridian); 
      patch(tmpx,tmpy,[.8 .8 .8],'linewidth',0.25)
   end

   % Plot grounding lines as patch objects: 
   for k = 1:30
      [tmpx,tmpy] = ll2ps(glat{k},glon{k},'meridian',meridian); 
      patch(tmpx,tmpy,[.5 .5 .5],'linewidth',0.25)
   end
  
   
   % Plot red box:
   if usekm
      plot(axx*1000,axy*1000,'r-','linewidth',1); 
   else
      plot(axx,axy,'r-','linewidth',1); 
   end
      
   axis equal tight
   
   
   % Set final dimensions after plotting the inset: 
   gpinset = plotboxpos(gca); 
   insetwidth = gpinset(3); 
   insetheight = gpinset(4); 
    
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
   
   
   % Format inset axes: 
   set(gca,'xtick',[],'ytick',[],'position',[insetx insety insetwidth insetheight])
   if frameon
      box on
   else
      axis off
   end
   
  
   % Make the original map axes current: 
   axes(gcah); 
   
   % Ensure inset map is on top of the stack: 
   uistack(gcah,'down');
   
   
  % Clean up: 
  if nargout==0 
     clear h
  end
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
        