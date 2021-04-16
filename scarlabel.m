function [h,featurelat,featurelon] = scarlabel(featurename,varargin)
% SCARLABEL labels Antarctic features on a map. Feature names and 
% locations correspond to 25,601 locations identified by the 
% Scientific Committee on Antarctic Research (SCAR).
% http://www.scar.org/cga
% 
%% Syntax
% 
% scarlabel('FeatureName')
% scarlabel('FeatureName','TextProperty',PropertyValue)
% scarlabel('FeatureName','Marker','MarkerStyle') 
% scarlabel(...,'km') 
% h = scarlabel(...)
% [h,featurelat,featurelon] = scarlabel(...)
% 
%% Description 
% 
% scarlabel('FeatureName') labels Antarctic features core location if a map
% is current. Multiple FeatureNames may be entered as a cell array (e.g.,
% {'McMurdo Station','Palmer Station','Casey Station'}.If the current axes are not map axes, polar stereographic cartesian
% coordinates are assumed.  
% 
% scarlabel('FeatureName','TextProperty',PropertyValue) formats the label with
% text properties as name-value pairs. 
% 
% scarlabel('FeatureName','Marker','MarkerStyle') places a marker at the
% center location. 
% 
% scarlabel(...,'km') places labels on cartesian coordinates in polar
% stereographic kilometers.  
% 
% h = scarlabel(...) returns the handle h of the label. 
% 
% [h,featurelat,featurelon] = scarlabel(...) returns geographic coordinates of
% the feature. 
%
%% Example 1
% % Before labelling we'll initialize a map with Matlab's ugly built-in coast line: 
% 
% antmap
% load coast 
% patchm(lat,long,[.5 .5 .5])
% 
% % Now let's start labelling.
% 
% scarlabel('Totten Glacier')
% 
% scarlabel('Recovery Glacier',...
%    'fontsize',10,...
%    'HorizontalAlignment','left',...
%    'marker','x',...
%    'color','blue')
% 
% scarlabel('Ross Ice Shelf',...
%     'edgecolor','magenta',...
%     'fontangle','italic',...
%     'backgroundcolor',[.4 .6 .8],...
%     'fontsize','20')
% 
% scarlabel('Marie Byrd Land','rotation',45)
% 
% scarlabel('South Pole',...
%     'fontweight','bold',...
%     'color','red',...
%     'marker','rp',...
%     'markersize',50)
% 
%% Example 2: Multiple locations at once:
% This example uses the modismoa function available on the Mathworks File Exchange site
% 
% modismoa('antarctic peninsula',1000)
% myplaces = {'Renaud Island','Longing Cape','Palmer Station'}; 
% scarlabel(myplaces,'color','b','backgroundcolor','w')
% scalebar('length',200,'color','w')
% 
%% Example 3: Get location
% % You can use |scarlabel| to return geographic locations of features found
% % by |scarloc|: 
% 
% [~,featurelat,featurelon] = scarlabel('elephant island');
% 
%% Example 4: Polar stereographic cartesian kilometers
% If current axes are not map coordinates, labeling is performed in polar
% stereographic cartesian coordinates (meters).  If, however, your data are plotted
% in polar stereographic kilometers, ensure correct labeling with the 'km'
% tag: 
% 
% load AMTdata
% plotps(glat{1},glon{1},'km')
% xlabel 'eastings (km)'
% ylabel 'northings (km)' 
% scarlabel('Recovery Glacier','km') 
% 
%% Known Issues
% Unfortunately, scarloc and scarlabel are plagued by problems with special
% characters. This includes letters with accents and umlauts. My sincerest
% apologies to our French, Norwegian, and just about everything aside from
% American friends.  
% 
%% Cite this dataset as
% 
% Secretariat SCAR (1992, updated 2015). Composite Gazetteer of Antarctica,
% Scientific Committee on Antarctic Research. GCMD Metadata 
% <http://gcmd.nasa.gov/records/SCAR_Gazetteer.html>
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
% Created by Chad A. Greene, July 2013. 
% The University of Texas Institute for Geophysics.
% Updated February 2015 to allow labeling non-georeferenced maps. 
% Updated July 2015 to allow polar stereographic kilometers. 
% 
% See also scarloc, scarclick, textm, corelabel, antmap, and coreloc.

%% Input checks: 

assert(nargin>0,'The scarlabel requires at least one input.  What are you trying to label?')
assert(isnumeric(featurename)~=1,'Scarlabel''s input featurename must be a string or a cell array of strings') 

%% Is this on a georeferenced map, or a polar stereographic map?  

mapit = false; 
try
    if ismap
        mapit = true; 
    end
end

plotkm = false; % if ps units, use meters by default
tmp = strcmpi(varargin,'km'); 
if any(tmp) 
    plotkm = true; 
    varargin = varargin(~tmp); 
end

%% Get location(s): 

if strcmpi(featurename,'glaciers')
    load scarnames.mat lat lon names featureType;  
    ind = strcmpi(featureType,'glacier'); 
    featurelat = lat(ind); 
    featurelon = lon(ind); 
    featurename = names(ind); 
elseif strcmpi(featurename,'ice shelves')
    load scarnames.mat lat lon names featureType;  
    ind = strcmpi(featureType,'ice shelf'); 
    featurelat = lat(ind); 
    featurelon = lon(ind); 
    featurename = names(ind);    
else
   [featurelat,featurelon] = scarloc(featurename); 
end

if ~mapit
    [fx,fy] = ll2ps(featurelat,featurelon); 
    if plotkm
        fx = fx/1000; 
        fy = fy/1000; 
    end
end

%% Place text label: 

% Apply the text label: 
if mapit
    h = textm(featurelat,featurelon,ones(size(featurelat)),featurename,'horiz','center','vert','middle');
else
    h = text(fx,fy,ones(size(featurelat)),featurename,'horiz','center','vert','middle','Clipping','on');
end

%% Place a marker if requested:  

tmp = strcmpi(varargin,'marker'); 
if any(tmp)
    if mapit
        h(2) = plotm(featurelat,featurelon,varargin{find(tmp)+1});
    else
        h(2) = plot(fx,fy,varargin{find(tmp)+1});
    end

    varargin{find(tmp)+1}=1; 
    varargin = varargin(~tmp); 
    set(h(1),'VerticalAlignment','bottom'); 
end

%% Format text and marker

% This is brute force, but it works:
for k = 1:length(h) 
    for k2 = 1:length(varargin)
        try
             set(h(k),varargin{k2},varargin{k2+1}); 
        end
    end
end

%% Clean up: 

% Delete output arguments if none are requested:
if nargout==0
    clear h featurelat featurelon; 
end
