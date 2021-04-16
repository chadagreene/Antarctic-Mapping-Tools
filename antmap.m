function [h] = antmap(varargin)
% ANTMAP initializes a polar stereographic map of the southern 
% hemisphere. This function is an adaptation of Andrew Roberts' 
% ncpolarm function. This function requires Matlab's Mapping Toolbox. 
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
%% Syntax
% 
% antmap
% antmap('northernlimit',latitude)
% antmap('grid')
% antmap('lats',Latitudes)
% antmap('lons',Longitudes)
% antmap(...,'LineProperty',LineValue) 
% antmap(...,'Frame','off') 
% antmap(...,'OceanColor',Color)
% antmap('lambert',latlim,lonlim)
% h = antmap(...)
% 
% 
%% Description 
% 
% antmap initializes a polar stereographic projection map of Antarctica.
% If a figure is current, the new map is placed in the current figure. If 
% no figure is current, a new figure is opened and a map is initialized. 
% 
% antmap('northernlimit',latitude) sets the northern limit as a scalar 
% latitude value between -90 and 0, northernmost limit of the map. Default 
% is -63. If a positive value is declared, it is automatically switched to
% a negative value to indicate southern hemisphere. 
%
% antmap('grid') places a grid or graticule on the map. 
% 
% antmap('lats',Latitudes) places grid lines a specified latitudes. 
% 
% antmap('lons',Longitudes) places grid lines a specified longitudes. 
% 
% antmap(...,'LineProperty',LineValue) sets line properties of a grid. 
% 
% antmap(...,'Frame','off') removes frame from map perimeter. Default
% value is 'Frame', 'off'. Note that turning the frame off also removes
% ocean color.  
% 
% antmap(...,'OceanColor',Color) specifies background color by Matlab's
% color names (e.g., 'red','blue',etc), abbreviations (e.g., 'r','b'), or
% RGB triple. 
% 
% antmap('lambert',latlim,lonlim) initializes a Lambert conic projection
% map of the region given by latlim, lonlim.  latlim and lonlim must be in
% the form latlim = [minlat maxlat] and lonlim = [minlon maxlon]. Note that
% this usage is not fully supported yet.  It might work, or it might not.
% Matlab has trouble updating maps to a different projection, so you'll
% have to initialize a lambert projection map and stick with it. 
% 
% h = antmap(...) returns a handle h of new axes.
% 
% 
%% Example 1: Plot an ugly and inaccurate rendition of Antarctica's coast line: 
%
% load coast 
% antmap
% plotm(lat,long) 
% 
% 
%% Example 2: Plot a graticule under and over that ugly continent, and specify 
% locations and LineSpec for another grid: 
%
% figure('position',[100 100 1000 400])
% 
% subplot(1,3,1)
% antmap('grid')              % initializes left map and plots grid
% patchm(lat,long,'y')        % plots yellow continent atop  grid
% 
% subplot(1,3,2)  
% antmap                      % initializes center map
% patchm(lat,long,'b')        % plots blue continent
% antmap('grid','color','m')  % overlays magenta grid
% 
% subplot(1,3,3) 
% antmap                      % initializes right map
% patchm(lat,long,[.5 .5 .5]) % plot gray continent
% antmap('lats',-80:10:-50,...% plots lines at 50, 60, 70, & 80 S
%     'lons',0:45:180,...     % and 0, 45, 95, etc longitude
%     'color','red',...       % in red
%     'linewidth',2,...       % kinda thick
%     'linestyle',':')        % and dotted. 
% 
% 
%% Example 3: Lambert projection (not fully supported)
%
% figure('position',[100 100 900 400])
% antmap('lambert',[-73 -68],[-15 15])
% bedmap2 patchshelves
% bedmap2 patchgl
% outlineashelf('fimbul ice shelf','linewidth',3) 
% scarlabel('Fimbul Ice Shelf','fontangle','italic',...
%     'color','b','fontweight','bold')
% 
%% Author Info
% This function has been adapted by Chad Greene from Andrew 
% Roberts' ncpolarm, which can be found here:
% http://www.mathworks.com/matlabcentral/fileexchange/30414
% 
% See also WORLDMAP and AXESM. 

%% Parse inputs: 

assert(license('test','map_toolbox') ==1,'It looks like you do not have the Mapping Toolbox.  I''ve heard that some folks have been able to modify this toolset to work with m_Map, thereby side-stepping the purchase of the Mapping Toolbox.  Maybe that''ll work for you?') 

% Set defaults: 
northernlimit = -63;
plotlats = false; 
plotlons = false; 
frameon = true; 
ChangeOceanColor = false; 
currentMapProjection = [];
makeLambert = false; 
lons = -170:10:180; 
lats = -85:5:0; 

% Reset default values if declared: 
if nargin>0
    
    tmp = strcmpi(varargin,'lambert'); 
    if any(tmp)
        tmpind = find(tmp); 
        latlim = varargin{tmpind+1}; 
        lonlim = varargin{tmpind+2};
        assert(isnumeric(latlim)==1&&numel(latlim)==2,'latlim must be a numeric array containing two values.')
        assert(isnumeric(lonlim)==1&&numel(lonlim)==2,'lonlim must be a numeric array containing two values.')
        latlim = -abs(latlim); % ensures southern hemisphere
        makeLambert = true; 
        tmp(tmpind+1:tmpind+2) = 1; 
        varargin = varargin(~tmp); 
        frameon = false;         
    end
    
    % Determine if northernlimit is declared: 
    tmp = strncmpi(varargin,'nor',3)|strcmpi(varargin,'limit');
    if any(tmp)
        northernlimit = varargin{find(tmp)+1}; 
        assert(isscalar(northernlimit)==1,'Northern Limit of antmap must be a scalar value.')
        northernlimit = -abs(northernlimit); % assumes southern hemisphere value
        tmp(find(tmp)+1)=1; 
        varargin = varargin(~tmp); 
    end
    
    % Determine if grid or graticule is declared: 
    tmp = strcmpi(varargin,'grid')|strncmpi(varargin,'grat',4); 
    if any(tmp)
        plotlats = true; 
        plotlons = true; 
        varargin = varargin(~tmp); 
    end
    
    % Determine if user declared lats:
    tmp = strncmpi(varargin,'lat',3);
    if any(tmp) 
        lats = varargin{find(tmp)+1};
        plotlats = true; 
        tmp(find(tmp)+1)=1; 
        varargin = varargin(~tmp); 
        
        % Reshape if necessary: 
        if length(lats)==1
            lats = [lats lats];
        end
        if iscolumn(lats)
            lats=lats';
        end
        lats = -abs(lats); % ensures southern hemisphere
    end
        
    % Determine if user declared lons: 
    tmp = strncmpi(varargin,'lon',3);
    if any(tmp) 
        lons = varargin{find(tmp)+1};
        plotlons = true; 
        tmp(find(tmp)+1)=1; 
        varargin = varargin(~tmp); 
        
        % reshape if necessary: 
        if length(lons)==1
            lons = [lons lons];
        end

        if iscolumn(lons)
            lons = lons';
        end
        lons = [lons lons-180]; 
    end
    
    % Determine if frame or box is declared: 
    tmp = strcmpi(varargin,'frame')|strcmpi(varargin,'box'); 
    if any(tmp)
        if strcmpi(varargin{find(tmp)+1},'off')||...
                strcmpi(varargin{find(tmp)+1},'none')
            frameon = false; 
        else
            frameon = true; 
        end 
        tmp(find(tmp)+1)=1; 
        varargin = varargin(~tmp); 
    end
    
    % Check if ocean color is declared: 
    tmp = strncmpi(varargin,'ocean',5);
    if any(tmp)
        oceancolor = varargin{find(tmp)+1}; 
        ChangeOceanColor = true; 
        tmp(find(tmp)+1)=1; 
        varargin = varargin(~tmp); 
    end
end


%% Initialize Map


try % Is a map already initialized? 
    currentMapProjection = getm(gca,'MapProjection');
end

if ~isempty(currentMapProjection) && ~strcmpi(get(gca,'visible'),'on') && ~frameon
    frameon = false; 
end
      

% Make map: 
if makeLambert && ~strcmpi(currentMapProjection,'lambert')
    parallelstep = (max(latlim)-min(latlim))/4; 
    mapparallels = [min(latlim)+parallelstep min(latlim)+3*parallelstep];
    h=axesm('MapProjection','lambert',...
        'mapparallels',mapparallels,...
        'maplatlimit',latlim,...
        'maplonlimit',lonlim,...
        'flatlimit',latlim,...
        'flonlimit',lonlim);
    
elseif isempty(currentMapProjection) 
    h=axesm('MapProjection','stereo',...
      'AngleUnits','degrees',...
      'Aspect','normal',...
      'FalseNorthing',0,...
      'FalseEasting',0,...
      'MapLatLimit',[-90 0],...
      'Geoid',[1 0],...
      'Origin',[-90 0 0],...
      'Scalefactor',1,...
      'Frame','off',...
      'FFill',2000,...
      'FLatLimit',[-Inf -90],...
      'FEdgeColor',.4*[1 1 1],...
      'FFaceColor','white',...
      'FLineWidth',0.25);

    tightmap 
    xlim=get(gca,'xlim');
    ylim=get(gca,'ylim');
    [~,leftlon]=minvtran(xlim(1),ylim(1)+0.5*diff(ylim));
    [~,rightlon]=minvtran(xlim(2),ylim(1)+0.5*diff(ylim));
    [~,bottomlon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(1));
    [~,toplon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(2));
    [x1] = mfwdtran(northernlimit,leftlon);
    [x2] = mfwdtran(northernlimit,rightlon);
    [~,y1] = mfwdtran(northernlimit,bottomlon);
    [~,y2] = mfwdtran(northernlimit,toplon);
    set(gca,'Xlim',[x1 x2]); set(gca,'Ylim',[y1 y2]);
end


%% Plot graticule; 

if plotlats
    lath=plotm(repmat(lats,361,1),repmat((-180:180)',1,length(lats)),...
        'color',.8*[1 1 1],'linewidth',.5);
    if ~isempty(varargin)
        set(lath,varargin{:});
    end
end

if plotlons
    [longrid,latgrid] = meshgrid(lons,lats); 
        lonh=plotm(latgrid,longrid,...
       'color',.8*[1 1 1],'linewidth',.5);
    if ~isempty(varargin)
        set(lonh,varargin{:});
    end
end

%% Set frame 

if ChangeOceanColor
    set(gca,'color',oceancolor); % Sets ocean color
end

if ~frameon
    set(gca,'visible','off'); % Removes the box frame if desired. 
    set(findall(gca, 'type', 'text'), 'visible', 'on');
else set(gca,'visible','on');
end

%% Clean up: 

if nargout==0 ; 
    clear h; 
end
 

