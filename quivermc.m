function [hvectors,cb] = quivermc(lat,lon,u,v,varargin)
% QUIVERMC is an adapted version of Andrew Roberts' ncquiverref. 
% Dr. Roberts' function and this function fix a couple of problems with Matlab's quiverm function.
% The two primary issues with quiverm are as follows: 
% 
% 1. Matlab's quiverm confoundingly mixes up u and v.  By convention in the
% Earth sciences, u is the zonal component where eastward is positive and
% v is the meridional component where northward is positive.  Matlab gets this wrong, but the
% quivermc function described here gets it right. 
% 
% 2. For reasons related to ship travel and some old legacy code from Navy
% guys decades ago, Matlab's quiverm scales vectors in a strange way that
% depends on latitude.  If you're plotting some absolute field like wind
% vectors, there is no physical reason that you would want to scale vectors
% in such a way that their zonal components shrink to zero at the poles.  
% 
% In addition to fixing the problems described above, quivermc also
% allows a few extra options including color settings, arrow density, and 
% options for converging or diverging flow.  
% 
%% Syntax 
% 
%  quivermc(lat,lon,u,v)
%  quivermc(...,'units',unitString)
%  quivermc(...,'color',arrowcolor)
%  quivermc(...,'colormap',colorMap)
%  quivermc(...,'colorbar','ColorbarLocation')
%  quivermc(...,'density',densityVal)
%  quivermc(...,'arrowstyle',arrowStyle)
%  quivermc(...,'linewidth',lineWidth)
%  quivermc(...,'reference',referenceScale)
%  h = quivermc(...)
%  [h,cb] = quivermc(...)
% 
%% Description 
% 
% quivermc(lat,lon,u,v) plots vectors of zonal and meridional components
% u and v at locations given by lat and lon. 
% 
% quivermc(...,'units',unitString) prints any user-specified units
% alongside a reference vector. 
%
% quivermc(...,'color',arrowcolor) sets all arrows to the color given by
% arrowcolor, which can be a string (e.g. 'blue', or 'r') or RGB value.
% Default color is black. 
%
% quivermc(...,'colormap',colorMap) colors vectors scaled relative to
% their magnitude using any colormap such as jet or autumn(256). 
% 
% |quivermc(...,'colorbar','ColorbarLocation')| places a colorbar at a
% specified location. The argument |'colorbar','on'| may be used to place a
% colorbar to the outside right of the plot, or a location may be set as 
% 
%   * 'EastOutside', 'vertical', or 'on' Outside right
%   * 'SouthOutside' or 'horizontal'  Outside bottom
%   * 'North Inside' plot box near top
%   * 'South Inside' bottom
%   * 'East Inside' right
%   * 'West Inside' left
%   * 'NorthOutside' Outside plot box near top
%   * 'SouthOutside' Outside bottom
%   * 'WestOutside' Outside left
%
% quivermc(...,'density',densityVal) allows user-declared downsampling of
% input data. By default, if the input grid is larger than about 50x50,
% quivermc will attempt to downsample your data to some dimensions close
% to 50x50.  This is because large datasets take time to plot, and the
% arrows become so small they're hard to see.  If your dataset is 400x400
% and you would like to plot about 100 arrows by 100 arrows instead of the 
% default 50x50, use the name-value pair 'density',25 to specify 25
% percent of the data are to be plotted in each dimension.  
%
% quivermc(...,'arrowstyle',arrowStyle) specifies whether the arrow's
% tail or tip is located at its respective data point. By default, arrows
% are centered about their data points.  To pin arrow tails at their data
% points, declare 'arrowstyle','tail' or 'arrowstyle','divergent'. To
% pin arrow heads at their data points choose 'arrowstyle','tip' or
% 'arrowstyle','convergent'.  
%
% quivermc(...,'linewidth',lineWidth) sets the linewidth of plotted
% arrows in points.  
%
% quivermc(...,'reference',referenceScale) declares the scale by which arrows are 
% plotted. Can be 'median', which scales vectors relative to the median magnitude, 
% 'max', which scales vectors relative to the maximum magnitude, or
% 'equal' to make all vectors equal in size.  The referenceScale may
% also be a scalar value of your choosing. Default is 'max'. 
%
% h = quivermc(...) returns the vector handles of plotted arrow objects. 
% 
% |[h,cb] = quivermc(...)| returns the handle |cb| of the colorbar. 
%
%
%% Examples: 
% 
% % Make up some data: 
% load wind u v 
% u = squeeze(u(:,:,1)); 
% v = squeeze(v(:,:,1)); 
% u(1:30,:) = u(1:30,:)+5;     
% u(31:end,:) = u(31:end,:)-3;    
% lat = repmat((20:-1:-14)',1,41);
% lon = repmat(-160:-120,35,1);
% 
% % Initialize a map: 
% worldmap([min(lat(:))-1 max(lat(:))+1],[min(lon(:))-1 max(lon(:))+1]);cla; 
%  
% % Use the quivermc function:
% quivermc(lat,lon,u,v)
% quivermc(lat,lon,u,v,'units','miles per hour') % specifies units
% quivermc(lat,lon,u,v,'color','r') % sets all arrows to 'red' 
% quivermc(lat,lon,u,v,'colormap',hot(256),'units','m/s') % colors arrows with a hot colormap
% quivermc(lat,lon,u,v,'density',33.3) % downsamples data by a factor of 3 
% quivermc(lat,lon,u,v,'density',33.3,'arrowstyle','divergent','color',[.4 .6 .2]) % pins arrows at their tails 
% quivermc(lat,lon,u,v,'density',33.3,'arrowstyle','divergent','color',[.4 .6 .2],'reftype','median') % sizes arrows relative to median magnitude 
% 
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
% Original function ncquiverref written by Andrew Roberts
% Naval Postgraduate School
% Tested using MATLAB Version 7.11.0.584 (R2010b)
%
% $Id: ncquiverref.m 254 2011-02-09 05:26:08Z aroberts $
% 
% Adapted into quivermc by Chad A. Greene
% Institute for Geophysics
% The University of Texas at Austin
% July 2014
% Tested on Matlab 2012b with Mac OSX 10.8.5
% 
% Updated August 2014 to include John Barber's calcticks function 
% Calcticks can be found here: http://www.mathworks.com/matlabcentral/fileexchange/30671
% 
% See also quiver, quiverm, and quiver3m.

%% Input checking: 

% make sure the mapping toolbox is present
h=ver('map') ; if isempty(h) ; error('Mapping toolbox not installed') ; end

% error checking of inputs
assert(nargin>=4,'Not enough inputs.'); 
assert(numel(lat)==numel(lon)&numel(lat)==numel(u)&numel(lat)==numel(v),'Dimensions of lat, lon, u, and v must agree.')
assert(isscalar(lat)==0,'Input lat, lon, u, and v must be a grid.') 
assert(isvector(lat)==0,'Input lat, lon, u, and v must be a grid.') 


%% Set defaults and change them depending on user preferences: 

refvec = false; 
colorbarOn = false; 

units = ''; 
tmp = strncmpi(varargin,'unit',4); 
if any(tmp) 
    units = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
    refvec = true; 
end

reftype = 'max'; 
tmp = strncmpi(varargin,'ref',3)|strncmpi(varargin,'scale',5); 
if any(tmp) 
    reftype = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
end

tmp = strcmpi(varargin,'colormap'); 
if any(tmp) 
    cmap = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
end

tmp = strcmpi(varargin,'colorbar'); 
if any(tmp) 
    colorbarOn = true; 
    cbLocation = varargin{find(tmp)+1}; 
    cbLocation = strrep(cbLocation,'vertical','WestOutside'); 
    cbLocation = strrep(cbLocation,'on','EastOutside'); 
    cbLocation = strrep(cbLocation,'horizontal','SouthOutside'); 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
end

if nargout>1 && colorbarOn==false
    error('Too many output arguments.')
end

veccol='k'; 
tmp = strcmpi(varargin,'color'); 
if any(tmp) 
    veccol = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
end

arrowStyle = 'centered';
tmp = strcmpi(varargin,'arrowstyle'); 
if any(tmp) 
    tmp2 = strncmpi(varargin,'div',3)|strcmpi(varargin,'tail');
    if any(tmp2)
        arrowStyle = 'divergent'; 
    end
    tmp2 = strncmpi(varargin,'con',3)|strcmpi(varargin,'tip')|strcmpi(varargin,'head');
    if any(tmp2)
        arrowStyle = 'convergent'; 
    end
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
end

changelinewidth = false; 
tmp = strncmpi(varargin,'linewi',6)|strncmpi(varargin,'wid',3)|...
     strcmpi(varargin,'arrowwidth')|strcmpi(varargin,'arrow width'); 
if any(tmp)
    lineWidth = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp);
    assert(isnumeric(lineWidth)==1,'Line width must be a scalar.')
    changelinewidth = true; 
end

skipstep = 1;  % If the data set is small, do not skip any vectors 
defaultWidth = 50; % but if the data set is large, plot about 50x50 arrows for a large square dataset
if sqrt(numel(lat))/defaultWidth>1
    skipstep = round(sqrt(numel(lat))/defaultWidth); 
end
    
tmp = strncmpi(varargin,'dens',4)|strncmpi(varargin,'arrowdens',8)|strncmpi(varargin,'downsa',6);
if any(tmp)
    declaredDensityVal = varargin{find(tmp)+1};
    assert(isnumeric(declaredDensityVal)==1,'Arrow density value must be numeric.') 
    assert(declaredDensityVal>0&&declaredDensityVal<=100,'Arrow density must be greater than 0 and no more than 100.') 
    skipstep = round(100/declaredDensityVal); 
end

lat = lat(1:skipstep:end,1:skipstep:end); 
lon = lon(1:skipstep:end,1:skipstep:end); 
u = u(1:skipstep:end,1:skipstep:end); 
v = v(1:skipstep:end,1:skipstep:end); 
 
% get current axis 
%h=get(gcf,'CurrentAxes');
h = gca; 
assert(ismap(h)==1,'Current axes must be map axes.') 

% If plotting on a matlab map, determine if the axes are map or cartesian
% coordinates, and if the former calculate mapping to plot axis, and 
% then do vector field otherwise just plot the vector field.

%% Begin heavy lifting: 

% get x and y location on the map
sz=size(lat);
mstruct=gcm;
[x,y] = mfwdtran(mstruct,lat,lon,h,'none');
xz=size(x);
if sz~=xz
error('Change in size of x using mfwdtran. Try changing surface to none in the code')
end

% get angle on the map, but do not distort the length according to the projection
% so that all vectors can use the same reference vector.  DO NOT project
% the length of the vector to be different in x and y directions.
[th,z] = cart2pol(u,v);
[thproj,~] = vfwdtran(mstruct,lat,lon,90*ones(size(lat)));
[u,v] = pol2cart(th+deg2rad(thproj),z);

% remove masked grid points from the input by filling coordinates with NaN;
x(isnan(u))=NaN;
y(isnan(u))=NaN;
magnitude = hypot(u,v); 


% Scale the vectors according to the reference arrow vector length based on
% the mean distance between grid points. This is a good measure, as it remains 
% constant for multiple plots using the same grid with different values.
x1=abs(diff(x')); x2=abs(diff(x)); 
y1=abs(diff(y')); y2=abs(diff(y));
[~,z1] = cart2pol(x1,y1); [~,z2] = cart2pol(x2,y2);
scalelength=min(mean(z1(~isnan(z1))),mean(z2(~isnan(z2))));

% Calculate reference vector length based on rounded median
% or maximum value of plot.  The default is median based.
if isnumeric(reftype)
	refval=reftype;
elseif strncmpi(reftype,'median',3)
        z(z==0)=NaN;
	refval=median(z(~isnan(z)));
elseif strcmpi(reftype,'max') 
	refval=max(z(~isnan(z)));
elseif strcmpi(reftype,'equal')
    magnitude = hypot(u,v); 
    z(z==0)=NaN;
	refval=median(z(~isnan(z)));
    u = u./magnitude; 
    v = v./magnitude; 
else
end

% Remove NaN values that will not be plotted
% and turn points into a row of coordinates
u=u(~isnan(x))';
v=v(~isnan(x))';
y=y(~isnan(x))';
x=x(~isnan(x))';

% Set arrow size (1= full length of vector)
arrow=0.40;

% set scale value based on refval and scale length
roundp=floor(log10(refval));
refval=floor(refval/(10^roundp))*(10^roundp);
scale=scalelength/refval;


% adjust whether arrows are centered, 
switch arrowStyle
    case 'centered' % (centered is the default) 
        % Center vectors over grid points
        xstart=x-0.5*scale*u;
        xend=x+0.5*scale*u;
        ystart=y-0.5*scale*v;
        yend=y+0.5*scale*v;
        
    case 'divergent'
        % diverging arrows: 
        xstart=x;
        xend=x+scale*u;
        ystart=y;
        yend=y+scale*v;
        
    case 'convergent'
        % arrows converge upon a point
        xstart=x-scale*u;
        xend=x;
        ystart=y-scale*v;
        yend=y;
end
        

% Get x coordinates of each vector plotted
lx = [xstart; x; ...
  xstart+(1-arrow/3)*(xend-xstart); ...
  xend-arrow*(scale*u+arrow*(scale*v)); ...
  xend; ...
  xend-arrow*(scale*u-arrow*(scale*v)); ...
  xstart+(1-arrow/3)*(xend-xstart); ...
  NaN(size(x))];

% Get y coordinates of each vector plotted
ly = [ystart; y; ...
  ystart+(1-arrow/3)*(yend-ystart); ...
  yend-arrow*(scale*v-arrow*(scale*u)); ...
  yend; ...
  yend-arrow*(scale*v+arrow*(scale*u)); ...
  ystart+(1-arrow/3)*(yend-ystart); ...
  NaN(size(y))];

% Plot the vectors
hvectors = line(lx,ly,'Color',veccol);

%% Color vectors with cmap if requested by user: 

if exist('cmap','var') 

    magnitude = hypot(u,v); 
    minmag = min(magnitude(:))
    maxmag = max(magnitude(:))
    magind = 1+round((length(cmap(:,1))-1)*(magnitude(:)-minmag)./(maxmag - minmag));
    if colorbarOn
        colormap(gca,cmap); 
        cb = colorbar('location',cbLocation); 
        [newticks,newticklabels] = calcticks([minmag maxmag]);
        newticks = (length(cmap(:,1))-1)*(newticks-minmag)/(maxmag-minmag);
        switch lower(cbLocation)
            case {'east','west','eastoutside','westoutside'}
                set(cb,'ytick',newticks,'yticklabel',newticklabels);
                ylabel(cb,units); 
                
            case {'north','south','northoutside','southoutside'}
                set(cb,'xtick',newticks,'xticklabel',newticklabels);
                xlabel(cb,units); 
                
            otherwise
                error('Invalid colorbar location.')
        end
    end
    
    for k = 1:length(magnitude(:)) 
        set(hvectors(k),'Color',cmap(magind(k),:))
    end
end
 
 
%% Draw the reference vector key at altitude 2 above the map and grid
if refvec

 % Get the reference text string, formatted to powers of ten if required
  reftext=[num2str(refval),' ',units,' '];
 
 % Get the current axis limits
 xlim=get(gca,'xlim'); xp1=xlim(1); xp2=xlim(2);
 ylim=get(gca,'ylim'); yp1=ylim(1); yp2=ylim(2);

 % set padding around the reference vector
 padx=diff(xlim)/100; 
 pady=diff(ylim)/100;

 % Set x position of reference vector
 xend=xp2-padx;
 xstart=xend-scalelength;

 % Plot reference text in lower right hand corner
 ht=text(xstart,yp1+pady,reftext,'Visible','off','Parent',gca,'FontSize',8.5,...
        'VerticalAlignment','Bottom','HorizontalAlignment','Right');
 textextent=get(ht,'Extent');

 % Draw patch over area of vector key 
 xl=textextent(1)-padx;
 xr=xp2;
 yb=yp1;
 yt=textextent(2)+textextent(4)+pady;
 hp=patch([xl; xl; xr; xr],[yb; yt; yt; yb],[2; 2; 2; 2],'w',...
          'LineWidth',0.5,'Parent',gca);
 uistack(hp,'top');

 % Redraw reference text on top of patch
 ht=text(xstart,(yb+yt)/2,2.1,reftext,'Parent',gca,'FontSize',8.5,...
         'VerticalAlignment','Middle','HorizontalAlignment','Right');

 % Set y position of reference vector
 yend=textextent(2)+textextent(4)/2;
 ystart=yend;

 % Get x coordinates of reference vector plotted
 lx = [xstart; ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      xend-arrow*scalelength; ...
      xend; ...
      xend-arrow*scalelength; ...
      xstart+(1-arrow/3)*(xend-xstart); ...
      NaN];

 % Get y coordinates of reference vector plotted
 ly = [ystart; ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      yend+arrow*(arrow*scalelength); ...
      yend; ...
      yend-arrow*(arrow*scalelength); ...
      ystart+(1-arrow/3)*(yend-ystart); ...
      NaN];

 % Get z coordinates of reference vector
 lz = 2*ones(size(ly));

 % Plot the reference vector
 hrefvec = line(lx,ly,lz,'Color',veccol);

end

%% Change arrow width if requested by user: 

if changelinewidth
    set(hvectors,'linewidth',lineWidth)
    if refvec
        set(hrefvec,'linewidth',lineWidth)
    end
end
    

%% Clean up: 

if nargout==0
    clear hvectors
end
end % function



































%% CALCTICKS FUNCTION BY JOHN BARBER: 

function [ticks,tickLabels,scaleStr,minorTicks,overhang] = ...
         calcticks(limits,orientation,varargin)
% Calculate ticks and ticklabels for specified limits and text size
%
% SYNTAX
%
%   TICKS = CALCTICKS
%   TICKS = CALCTICKS(LIMITS)
%   TICKS = CALCTICKS(LIMITS,ORIENTATION)
%   TICKS = CALCTICKS(...,TEXTSIZE)
%   TICKS = CALCTICKS(...,SCALE)
%   TICKS = CALCTICKS(...,SEPARATEEXPONENT)
%   TICKS = CALCTICKS(...,EXPONENTFONTSIZE)
%   TICKS = CALCTICKS(...,MAXCHARS)
%   TICKS = CALCTICKS(AXHANDLE,...)
%   [TICKS,TICKLABELS,SCALESTR] = CALCTICKS(...)
%   [...,MINORTICKS] = CALCTICKS(...)
%   [...,OVERHANG] = CALCTICKS(...)
%
% DESCRIPTION
%
% TICKS = CALCTICKS Calculate ticks for the y-axis of the current axes, 
% using the axes' limits and text properties.
%
% TICKS = CALCTICKS(LIMITS) Calculate ticks for the y-axis of the current
% axes, using the specified limits instead of the axes limits.
%
% TICKS = CALCTICKS(LIMITS,ORIENTATION) Calculate ticks for the x or y axis
% of the current axes, using the specified limits. ORIENTATION can be any
% of 'x','h','horizontal' to get ticks for the x-axis, and any of 'y','v',
% or 'vertical' for the y-axis.
%
% TICKS = CALCTICKS(...,TEXTSIZE) Calculate ticks using the specified text
% size.  TEXTSIZE should be the size (height or width, depending on the
% selected orientation) of the string '2', in data units, using the desired
% font properties and axes size.  TEXTSIZE is used to determine the maximum 
% number of ticks that will fit in the specified data limits.  See the
% REMARKS section for more information about determining the correct value
% for TEXTSIZE.  If TEXTSIZE is not specified, CALCTICKS will calculate its
% value using the specified limits and text and position properties of the
% specified axes, or the current axes.
%
% TICKS = CALCTICKS(...,SCALE) Calculate ticks using the specified axis
% scaling.  Valid inputs are 'linear' and 'log'.  If SCALE is not 
% specified, CALCTICKS will use the value of the 'XScale' or 'YScale'
% property of the specified axes or the current axes. 
%
% TICKS = CALCTICKS(...,SEPARATEEXPONENT) If TRUE, calculate ticks and
% ticklabels, returning a separate string containing the data scale when
% the ticklabels use exponential notation.  The scale string is of the form
% 'x 10^NN' where NN is the scale of the maximum absolute value of the
% limits, and the tick labels will be of the form '-1.2345', normalized to
% 10^NN.  
%
% CALCTICKS determines automatically when to use exponential notation,
% and this setting will have no effect if the algorithm selects standard
% notation.  If SEPARATEEXPONENT is FALSE and the determination is to use
% exponential notation, the ticklabels will be of the form '1.234e+011'.
%
% TICKS = CALCTICKS(...,EXPONENTFONTSIZE) If TRUE, include TEX markup in
% the ticklabel strings to set the font size of exponents to 7 (the
% default). If set to a number, use that value for the font size of
% exponents.  If FALSE or 0, CALCTICKS will not include TEX markup to
% change the font size of exponents.
%
% TICKS = CALCTICKS(...,MAXCHARS) Set the maximum length (in characters) of
% ticklabels. This value determines the precision displayed in the
% ticklabels, and for horizontal (x) orientation, affects the tick spacing.
% The default maximum is 9 characters.  Setting the value of MAXCHARS too
% low can result in invalid outputs.  The actual label lengths are
% determined by the size of the tick interval relative to the data scale.
%
% [...,TICKLABELS] = CALCTICKS(...) Return a cell array of tick
% labels
%
% [...,SCALESTR] = CALCTICKS(...) In addition to ticks and ticklabels, 
% return a separate string containing the scale for ticklabels displayed 
% using exponential notation.  If the ticklabels are not displayed using
% exponential notation (as determined internally by CALCTICKS), SCALESTR
% will be the empty string.
%
% [...,MINORTICKS] = CALCTICKS(...) For 'log' scale, return a vector of
% minor ticks spaced at the [2:9] points in each decade.  If the major
% ticks are spaced at intervals greater than 3 decades, the minor ticks
% will be placed at the 'missing' decades.  For 'linear' scale, or if the 
% limits span less than a decade, MINORTICKS will be empty.
%
% [...,OVERHANG] = CALCTICKS(...) Return a 1x2 vector OVERHANG containing
% the distances (in data units) that the outermost tick labels extend from
% the lower and upper axes limits.  If the chosen tick interval results in
% the outermost ticks being inset from the data limits by at least one half
% of the label size, OVERHANG will be zero.
%
% 
% REMARKS
%
% Note that CALCTICKS does not draw the calculated ticks and ticklabels,
% but simply returns their values.  See below for a usage example.
%
% If the actual values of LIMITS are chosen by CALCTICKS as the outermost
% ticks, those ticks will be exactly the values of LIMITS.  Interior tick
% values can vary slightly from exact intervals due to floating point 
% precision limitations.  Ticks within 10*eps(min(abs(LIMITS))) of zero
% are rounded to zero.
%
% If the precision needed to display the tick values is greater than the
% number of characters specified by MAXCHARS, the values returned in the
% TICKLABEL strings will be truncated to MAXCHARS.
%
% To determine text size in data units, the following method can be used:
%
% First, ensure that the x or y limits of the axes are set to the desired
% value, and that the axes and parent figure are the intended size. As an
% alternative, perform the following with limits of [0 1], then multiply
% the result by the difference of the desired limits.
% 
% hTest = text(1,1,'2','units','data');
% ext = get(hTest,'Extent');
% delete(hTest)
%
% % For horizontal (x) orientation:
% textSize = ext(3);
%
% % For vertical (y) orientation:
% textSize = ext(4);
%
% IMPORTANT: The text size in data units is (by definition) relative to
% the limits of the data.  If the axes limits change after getting the 
% text size, the value for textSize will be incorrect, and should be 
% renormalized to the new limits.  
%     If the axes position changes (e.g. due to a figure resize), the
% value for textSize will be incorrect because the axes size (in data 
% units) remains the same while both the axes size (in absolute units)
% and the text size in points (absolute units) do not.  Be aware that
% events such as changing the axes limits, ticks, x or y label, title, or
% other axes properties often causes the axes to resize automatically.
%
%
% EXAMPLE
%
% Results will vary depending on monitor resolution.  On a monitor running 
% at 1280x1024 pixels at 96 dpi, the values shown here result in the axes'
% default xticklabels overlapping, and the yticklabels are not displayed
% with sufficient precision.
% 
% % Create a figure, set limits and plot a curve
% figure('Position',[360   502   480   360])
% hAx = axes;
% set(hAx,'FontSize',12)
% s = get(hAx);
% 
% xlimits = [1200 1200.003];
% ylimits = [1.2e-6 1.20003e-6];
% 
% x = linspace(xlimits(1),xlimits(2),101);
% y = sin(1e4*x)*0.4*diff(ylimits) + mean(ylimits);
% plot(x,y)
% 
% set(hAx,'XLim',xlimits,'YLim',ylimits)
% 
% % Get ticks for the x axis using calcticks
% [xTicks,xTickLabels] = calcticks(xlimits,'x');
% 
% % Plot the calculated ticks and labels
% tw = diff(ylimits)*.02;
% dy = ylimits(1)+[tw;2*tw];
% hXTicks = line(repmat(xTicks,2,1),repmat(dy,1,length(xTicks)),'Color','b');
% 
% hXTickLabels = text(xTicks',repmat(dy(2),length(xTicks),1),xTickLabels,...
%     'Color','b','HorizontalAlignment','center','verticalAlignment',...
%     'bottom','FontAngle',s.FontAngle,'FontName',s.FontName,'FontSize',...
%     s.FontSize,'FontWeight',s.FontWeight);
% 
% % Get y ticks using calcticks
% [yTicks,yTickLabels,scaleStr] = calcticks;
% 
% % Plot y ticks and labels
% tw = diff(xlimits)*.02;
% dx = xlimits(1)+[tw;2*tw];
% hYTicks = line(repmat(dx,1,length(yTicks)),repmat(yTicks,2,1),'color','r');
% 
% hYTickLabels = text(repmat(dx(2),length(yTicks),1),yTicks',yTickLabels,...
%     'Color','r','HorizontalAlignment','left','VerticalAlignment',...
%     'middle','FontAngle',s.FontAngle,'FontName',s.FontName,'FontSize',...
%     s.FontSize,'FontWeight',s.FontWeight);
% 
% % Now, set the new values as the axes' ticks and ticklabels and delete the
% % temporary text and labels.  
% set(hAx,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,...
%     'YTickLabel',yticklabels)
% delete([hXTicks; hXTickLabels; hYTicks; hYTickLabels])
% 
% % If the scale is not displayed for the y-axis, manually place the 
% % scale string
% text(xlimits(1),ylimits(2),scaleStr,'HorizontalAlignment','left',...
%     'VerticalAlignment','bottom','FontAngle',s.FontAngle,'FontName',...
%     s.FontName,'FontSize',s.FontSize,'FontWeight',s.FontWeight)
%
%
% See also AXES

% $$FileInfo
% $Filename: calcticks.m
% $Path: $toolboxroot/
% $Product Name: calcticks
% $Product Release: 1.1
% $Revision: 1.1.5
% $Toolbox Name: Custom Plots Toolbox
% $$
%
% Copyright (c) 2010-2011 John Barber.
%
% Release History:
% v 1.0 : 2011-Mar-08
%       - Initial release
% v 1.1 : 2011-Mar-29
%       - Fixed bug that caused log-scale ticklabels to be truncated
%       - Improved interval selection for log-scale ticks
%       - Moved into Custom Plots Toolbox 


%% Constants

% Default limit on label length (in characters)
defMaxChars = 9;

% Minimum value of the (upper) limit on label length (characters). Setting 
% this value too small will cause problems.
% Note: This value does not affect the minimum length of the ticklabels.
minChars = 6;

% Default font size for exponents (assumes that font units are 'points')
defExpFontSize = 7;

% Upper limit on number of ticks returned by CALCTICKS
initMaxTicks = 11;

% Multiplier for textSize for vertical orientation when scale is 'log', to
% account for labels using exponential notation.
vertExpScale = 1.3;

%% Parse inputs

nargs = nargin;

% Handle empty input
if nargin == 0
    limits = [];
    orientation = 'v';
    varargin = cell(0,0);
end

% Check for an axes handle as first argument
if isscalar(limits) && ishandle(limits) && strcmp(get(limits,'Type'),'axes')
    hAx = limits;
    if nargs == 2
        limits = orientation;
        orientation = [];
    elseif nargs > 2
        limits = orientation;
        orientation = varargin{1};
        varargin(1) = [];
    end
    nargs = nargs-1;
else
    % Leave hAx empty unless we absolutely need it
    hAx = [];
end

% Validate orientation first so we can get the right axes limits if needed
% Orientation: {'v'} or 'h', also 'x' or 'y'
if nargs < 2 || isempty(orientation) || ...
        ~any(strcmpi(orientation(1),{'h','x'}))
    orientation = 'v';
    axLim = 'YLim';
    axScale = 'YScale';
else
    orientation = 'h';
    axLim = 'XLim';
    axScale = 'XScale';
end

if nargs == 0 || isempty(limits)
    % Get limits from an axes handle passed in as first argument, or gca.
    if isempty(hAx)
        hAx = gca;
    end
    limits = get(hAx,axLim);
elseif ~isreal(limits) || ~all(size(limits) == [1 2]) || ...
       (limits(1) >= limits(2))
    eID = [mfilename ':InvalidLimits'];
    eStr = '''limits'' must be a 1x2 vector with limits(2) > limits(1).';
    error(eID,eStr)
end

% Text size
if nargs < 3 || isempty(varargin{1})
    if isempty(hAx)
        hAx = gca;
    end
    textSize = getTextSize(limits,orientation,hAx);
else
    textSize = varargin{1};
end

% Scale
if nargs < 4 || isempty(varargin{2})
    if isempty(hAx)
        hAx = gca;
    end
    scale = get(hAx,axScale);
else
    scale = varargin{2};
end

% Exponent string style
if nargs < 5 || isempty(varargin{3})
    separateExp = true;
else
    separateExp = varargin{3};
    if ischar(separateExp)
        % Accept 'y(es)', 't(rue)', 'o(n)', 's(eparate)' as true
        separateExp = any(strcmpi(separateExp(1),{'y','t','o','s'}));
    else
        separateExp = logical(separateExp(1));
    end
end

% Handle expFontSize, set a flag to use or not use this value
if nargs < 6 || isempty(varargin{4})
    smallExp = true;
    expFontSize = defExpFontSize;
else
    expFontSize = varargin{4};
    if islogical(expFontSize)
        smallExp = expFontSize;
        expFontSize = defExpFontSize;
    else
        smallExp = ~isnan(expFontSize) && expFontSize > 0;
    end
end

% Maximum number of characters in label string.  Determines numerical
% precision displayed by the labels, and also affects the tick spacing for
% horizontal orientation.
if nargs < 7 || isempty(varargin{5}) || ~(isnumeric(varargin{5}) && ...
        isscalar(varargin{5}))
    maxChars = defMaxChars;
else
    maxChars = varargin{5};
    if maxChars < minChars
        maxChars = minChars;
    end
end
        

%% Initial calculations

% Bypass to logticks calculation if scale is 'log'
if ~isempty(scale) && strcmpi(scale(1:2),'lo')
    [ticks,tickLabels,overhang,minorTicks] = logticks(limits,...
        textSize,orientation,smallExp,expFontSize,maxChars,vertExpScale);
    scaleStr = '';
    return
else
    % No minor ticks for linear scale
    minorTicks = [];
end

% Data range
range = diff(limits);

% Get eps values for rounding
lEps = eps(limits(1));
uEps = eps(limits(2));
minEps = min(lEps,uEps);

% Vector of allowed tick counts
testTickCounts = 2:initMaxTicks;

% Make a list of rough intervals as a starting point
roughInts = (range./(testTickCounts-1))';

% Vector of 'nice' intervals
niceVec = [1 2 5 10];

%% Find nice intervals

% Normalize rough intervals by their scale
decRoughInts = floor(log10(roughInts));
normRoughInts = roughInts./10.^decRoughInts;

% Get the distances to nice intervals, pick the shortest
deltas = abs(repmat(normRoughInts,1,length(niceVec)) - ...
             repmat(niceVec,length(normRoughInts),1));
[trash,idx] = min(deltas,[],2); %#ok<ASGLU>

% Get the nice intervals and scores
niceInts = niceVec(idx)'.*10.^decRoughInts;  

% Remove duplicates
niceInts = unique(niceInts);

% Get upper and lower limits, fixed by the list of nice intervals.  Round
% out to make sure we get ticks at the original limits.
lLims = floor(limits(1)./niceInts).*niceInts;
uLims = ceil(limits(2)./niceInts).*niceInts;

% Get tick counts using the list of nice intervals and limits
nTicks = floor(1 + (uLims - lLims + 10*minEps)./niceInts);

% Shrink nice limits that are outside of original limits
idx = lLims < limits(1) - 10*eps(limits(1));
nTicks(idx) = nTicks(idx)-1;
lLims(idx) = lLims(idx) + niceInts(idx);
idx = uLims > limits(2) + 10*eps(limits(1));
nTicks(idx) = nTicks(idx)-1;
uLims(idx) = uLims(idx) - niceInts(idx);

% Set values that are almost exactly the original limits to be the original
% limit value.
idx = abs(lLims - limits(1)) < 10*eps(limits(1));
lLims(idx) = limits(1);
idx = abs(uLims - limits(2)) < 10*eps(limits(2));
uLims(idx) = limits(2);

% Discard values where the limits are reversed or equal
idx = (lLims >= uLims);
lLims(idx)=[];
uLims(idx)=[];
nTicks(idx)=[];
niceInts(idx)=[];

%% Determine label size for each interval

% Get the decade span of the limits and the decade of the intervals
maxAbs = max(abs([lLims uLims]),[],2);
decMax = floor(log10(maxAbs));% - nDec;
decInts = floor(log10(niceInts));% - nDec;

% Get the number of characters needed for tick labels for normal notation
labelChars = max(decMax+1,1) + (decInts<0).*(1-decInts);
labelChars(labelChars > maxChars - 1) = maxChars - 1;

% Handle exponential notation

% Determine whether or not to use exponential notation
if separateExp
    % Large numbers:
    isExp = (decMax > 6) | (decMax == 6 & decInts > 0);
    % Small numbers:
    isExp = isExp | (decMax < -3) | (decMax == -3 & decInts < -5);
else
    % Large numbers:
    isExp = decMax > 6 | (decMax == 6 & decInts > 3); 
    % Small numbers:
    isExp = isExp | (decMax < -3) | (decMax == -3 & decInts < -5);   
end

% Get length of exponential labels depending on style
if separateExp
    expChars = 2 + max(0,min(maxChars-3,decMax-decInts));
%     scaleSignChar = decMax < 0;
%     scaleChars = 5 + scaleSignChar + max(0,floor(log10(abs(decMax))));
else
    expChars = 2 + max(0,min(max(0,maxChars-7),decMax-decInts)) + 4;
end

% Select between normal and exponential label lengths
labelChars(isExp) = expChars(isExp);

% For consistency, always include space for a negative sign, regardless of
% the sign of the actual limits
labelChars = labelChars + 1; 
% % Uncomment to not include the negative sign space if it isn't needed
% labelChars = labelChars + (lLims < 0);

% Get label size based on textSize, orientation and length of label string
if strcmp(orientation,'h')
    labelSize = textSize * labelChars;
else
    labelSize = textSize * ones(size(labelChars));
end

%% Choose the best interval

% Maximum number of ticks without overlapping labels
nMax = floor((uLims-lLims+10*minEps)./labelSize) + 1;

% Modify this value based on initMaxTicks
nMaxScore = min(nMax,initMaxTicks + 0.25*(nMax-initMaxTicks));

% Calculate a score based on the number of ticks relative to the maximum.
nTickScore = 1-(nTicks./nMaxScore - 0.7).^2;

% Severe penalty for more than nMax ticks
penalty = 4*nTicks./nMax;
idx = nTicks./nMax <= 1;
penalty(idx) = 0;

% Penalty for more than 0.75*nMax but less than nMax
penaltyScale = nTicks./nMax;
idx = penaltyScale > 0.75 & idx;
penalty(idx) = 1*nTickScore(idx).*penaltyScale(idx);

% Test for intervals that divide the limits exactly
rangeTest = range./niceInts;
isInt = abs(rangeTest-round(rangeTest)) < 1e-6;

% Test for intervals that land exactly on the endpoints
hitsEnds = (abs(lLims-limits(1))<100*lEps) & ...
           (abs(uLims-limits(2))<100*uEps);

% Compute a score using the above tests and the tick score and penalty
scores = isInt + 0.75*hitsEnds + nTickScore - penalty;

% Penalize for too few ticks
idx = (nMax > 5) & (nTicks < 4);
scores(idx) = scores(idx) - 0.5*scores(idx);

idx = (nMax > 4) & (nTicks == 2);
scores(idx) = scores(idx) - 0.5*scores(idx);

% Find the highest score
[maxScore,bestIdx] = max(scores);

if maxScore > -1
    % Use the limits and interval with the best score
    lLims = lLims(bestIdx);
    uLims = uLims(bestIdx);
    interval = niceInts(bestIdx);
else
    % If the best score is too low, just return two ticks at the limits
    lLims = limits(1);
    uLims = limits(2);
    interval = range;
end

% Create the vector of ticks, making sure to hit lLims and uLims exactly
ticks = [lLims:interval:(uLims-interval/2) uLims];

% Handle 0 as a special case
zeroIdx = abs(ticks)<10*minEps;
ticks(zeroIdx) = 0;

% Calculate how far the label overhangs (in data units) from both ends
% (0 <= overhang <= labelSize/2)
labelSize = labelSize(bestIdx);
overhang = max(0,limits(1)+labelSize/2-lLims);
overhang(2) = max(0,uLims-limits(2)+labelSize/2);

%% Create tick labels

% Get the decade span of the limits and the decade of the intervals
maxAbs = max(abs([ticks(1) ticks(end)]));
decMax = floor(log10(maxAbs));
decInt = floor(log10(interval));

if isExp(bestIdx) && separateExp
    % Exponential using 'nice' notation w/ separate exponent string:
    % Label: '1.23'  scaleStr: 'x 10^34'
    
    % Determine maximum number of decimals and set formatting string
    n = max(1,min(max(0,maxChars-3),decMax-decInt));
    fStr = ['%.' num2str(n) 'f'];
    
    % Normalize the ticks to the correct scale and create the labels
    normTicks = ticks/10^(decMax);
    tickLabels = strtrim(cellstr(num2str(normTicks',fStr)))';
    
    % Create the scale string
    % Use smaller font for exponent if requested
    if smallExp
        fs = ['\fontsize{' num2str(expFontSize) '}'];
    else
        fs = '';
    end
    
    scaleStr = ['x 10^{' fs num2str(decMax) '}'];

elseif isExp(bestIdx) 
    % Exponential using 'ugly' notation: '-2.34e+301'
    
    % Determine maximum number of decimals and set formatting string
    n = max(0,min(max(0,maxChars-7),decMax-decInt));
    fStr = ['%.' num2str(n) 'e'];
    
    % Create tick labels
    tickLabels = strtrim(cellstr(num2str(ticks',fStr)))';
    
    % Output 'e' for the scale string as a flag that we are using
    % exponential notation.
    scaleStr = 'e';
    
else
    % Normal (fixed point) notation
    
    % Determine maximum number of decimals and set formatting string
    n = max(0,min(-decInt,max(0,maxChars-3)-decMax));
    fStr = ['%.' num2str(n) 'f'];
    
    % Create tick labels
    tickLabels = strtrim(cellstr(num2str(ticks',fStr)))';
    scaleStr = '';
   
    tickLabels(zeroIdx) = {'0'};
    
end

end % End of calcticks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ticks,tickLabels,overhang,minorTicks] = logticks(limits,...
    textSize,orientation,smallExp,expFontSize,maxChars,vertExpScale)
% Calculate ticks, etc. using a log scale.  Also returns minor ticks at
% [2:9] points in each decade.  Limits must be non-zero and positive.

%%  Input check 
% Error if limits are <= 0.  We already ensured that they are increasing
% before calling this subfunction.
if any(limits <= 0)
    eID = [mfilename ':InvalidLimits'];
    eStr = '''limits'' must be non-zero and positive for log scale.';
    error(eID,eStr)
end

%% Calculate ticks

% Range in linear space
range = limits(2) - limits(1);

% Range in log space
limits10 = log10(limits);
normRange = limits10(2) - limits10(1);

% Number of decades spanned by the limits
decMax = floor(limits10(2));
decMin = ceil(limits10(1));
decRange = max(0,decMax - decMin);

% Force decRange to be decRange+1 if either of the limits is on a decade
% boundary.
if any(floor(limits10)==ceil(limits10))
    decRange = decRange + 1;
end

% Normalize textSize to the linearized data range.  
normTextSize = textSize*normRange/range;

% Get number of characters needed for label
if normRange > 0.5
    nDigits = 0;
    expChars = max(1,max(floor(log10(abs(decMin))),...
        floor(log10(abs(decMax)))));
    labelChars = 2 + expChars + (decMax < 0 | decMin < 0);
else
    nDigits = max(1,min(maxChars-4,-floor(log10(normRange))+1));
    expChars = 1 + max(1,floor(log10(abs(decMax)))) + nDigits;
    labelChars = 2 + expChars + (decMax < 0 | decMin < 0);
end

% Determine label size
if strcmp(orientation,'h')
    labelSize = normTextSize*labelChars;
else
    % Scale text size by vertExpScale to make room for exponents
    labelSize = vertExpScale*normTextSize;
end

% Get maximum number of ticks assuming decade intervals.  (Will recalculate
% this value using the actual limits for case (4)).
maxTicks = round(decRange/(2*labelSize));

% Initial values
decadeTicks = true;
lLim = decMin;
uLim = decMax;

% Cases:
% (1) decRange > 2 and maxTicks >= 2
% (2) decRange > 2 and maxTicks < 2
% (3) decRange > 1 or floor(log10(limits(1))) ~= decMax
% (4) decRange < 1 and no decade in the interval

if decRange >= 2 && maxTicks >= 2
    % Case (1): Span is multiple decades.  Get ticks at each decade, or at
    % multiple-decade intervals.
    if decRange <= maxTicks
        % If they all fit, just return a vector of ticks at each decade in
        % the interval
        logTicks = decMin:decMax;       
    else
        % Select a nice interval and return a vector of ticks spaced by
        % that interval (in decades)
        niceVec = [1 2 5 10 20 50 100 200 500];
        roughInt = decRange/(maxTicks-1);
        deltas = abs(niceVec - roughInt);
        interval = niceVec(deltas == min(deltas));
        interval = interval(1);
        lLim = ceil(decMin/(interval))*interval;
        uLim = floor(decMax/(interval))*interval;
        logTicks = lLim:interval:uLim;
    end
elseif decRange >= 2 && maxTicks < 2
    % Case (2): Span is multiple decades but there is only room for 1 tick.
    % Return 1 tick at the center
    logTicks = round(mean([decMin decMax]));
    
elseif floor(limits10(1)) ~= decMax
    % Case (3): Span is < 2 decades, but crosses at least 1 decade 
    % boundary.  Return tick(s) at the decade boundary(s).
    if decMin > floor(limits10(1))
        logTicks = unique([decMin decMax]);
    else
        logTicks = decMax;
    end
else
    % Case (4): Span is < 1 decade and does not cross a decade boundary.
    % Return ticks at non-decade intervals, based on maxTicks and the data
    % range.
    decadeTicks = false;
    
    nDigits = max(1,min(maxChars-4,-floor(log10(normRange))));
    expChars = 1 + max(1,floor(log10(abs(decMax)))) + nDigits;
    labelChars = 2 + expChars + (decMax < 0 | decMin < 0);
 
    if strcmp(orientation,'h')
        labelSize = normTextSize*labelChars;
    end  
    
    maxTicks = round(normRange/(1.5*labelSize));
    if maxTicks < 3
        % Return 2 ticks at the limits
        logTicks = limits10;
    else
        % Get 'nice' intervals for the exponent
        roughInt = normRange/(maxTicks-1);
        scale = floor(log10(roughInt));
        niceVec = [1 2 4 5 10];
        niceInts = niceVec*10^scale;
        nTicks = floor(limits10(2)./niceInts + 1e6*eps(limits10(2))) - ...
                 ceil(limits10(1)./niceInts - 1e6*eps(limits10(1))) + 1;
        tooMany = nTicks > maxTicks;
        [trash,idx] = sort(nTicks,'descend');  %#ok<ASGLU>
        niceInts = niceInts(idx);
        tooMany = tooMany(idx);
        bestInt = niceInts(find(~tooMany,1,'first'));
        if isempty(bestInt)
            logTicks = limits10;
        else
            lLim = ceil(limits10(1)/bestInt - 1e6*eps(limits10(1)))*bestInt;
            uLim = floor(limits10(2)/bestInt + 1e6*eps(limits10(2)))*bestInt;
            logTicks = lLim:bestInt:uLim;
        end
    end
end

% Convert to linear space
ticks = 10.^logTicks;

% Calculate how far the label overhangs (in data units) from both ends
% (0 <= overhang <= labelSize/2)
overhang = max(0,limits10(1)+labelSize/2-logTicks(1));
overhang(2) = max(0,logTicks(end)-limits10(2)+labelSize/2);

%% Tick labels

% Formatting string
if decadeTicks
    fStr = '%.0f';
else
    fStr = ['%.' num2str(nDigits) 'f'];
end

% Use smaller font for exponents if requested
if smallExp
    fs = ['\fontsize{' num2str(expFontSize) '}'];
else
    fs = '';
end

baseStr = '10^{';
endStr = '}';

tickLabels = strtrim(cellstr(num2str(logTicks',fStr)));
tickLabels = strcat(baseStr,fs,tickLabels,endStr);


%% Minor ticks
% Get minor ticks at the [2:9] values in each decade.  If major ticks skip
% decades, get minor ticks at the decades that were skipped.  If the major
% ticks are not at decade boundaries, do not return any minor ticks.

if ~decadeTicks
    minorTicks = [];
else
    if length(logTicks) > 1
        skip = logTicks(2) - logTicks(1);
    else
        skip = 0;
    end
    
    if skip < 3
        % Do minor ticks at [2:9] in each decade
        mVec = floor(log10(limits(1))):ceil(log10(limits(2)));
        minorTicks = 10.^(sort(repmat(mVec,1,8)) + repmat(log10(2:9),1,length(mVec)));
        minorTicks(minorTicks < limits(1)) = [];
        minorTicks(minorTicks > limits(2)) = [];
        
        if skip == 2
            minorTicks = [minorTicks 10.^((lLim+1):2:uLim)];
            minorTicks = sort(minorTicks);       
        elseif skip == 3
            minorTicks = [minorTicks 10.^((lLim+1):3:uLim)];
            minorTicks = [minorTicks 10.^((lLim+2):3:uLim)];            
            minorTicks = sort(minorTicks);                   
        end
    else
        % Do minor ticks at the decades not included in logTicks
        minorLogTicks = decMin:decMax;
        for k = 1:length(logTicks)
            minorLogTicks(minorLogTicks == logTicks(k)) = [];
        end
        minorTicks = 10.^minorLogTicks;
    end
end


end % End of calcticks/logticks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function textSize = getTextSize(limits,orientation,hAx)
% Determine the size of text in data units.  This value is dependent on
% font size, axes data limits and the size of the axes on screen. 

% Notes:
% - There doesn't appear to be a problem with the value returned for
%   'Extent' if the text object is located outside of the axes limits.  
% - 'Extent' does not appear to change when the axes scale is set to 'log'.
% - Requires a valid axes handle.
% - 'Extent' is not valid for 3D views.

% Get axes properties
s = get(hAx);

% Get text size in data units
hTest = text(1,1,'2','Units','data','FontUnits',s.FontUnits,...
    'FontAngle',s.FontAngle,'FontName',s.FontName,'FontSize',s.FontSize,...
    'FontWeight',s.FontWeight,'Parent',hAx);
textExt = get(hTest,'Extent');
delete(hTest)
textHeight = textExt(4);
textWidth = textExt(3);

% If using a proportional font, shrink text width by a fudge factor to
% account for kerning.
if ~strcmpi(s.FontName,'FixedWidth')
    textWidth = textWidth*0.8;
end

% Restore axes limits and set output
if strcmp(orientation,'h')
    textSize = textWidth*(limits(2)-limits(1))/(s.XLim(2)-s.XLim(1));
else
    textSize = textHeight*(limits(2)-limits(1))/(s.YLim(2)-s.YLim(1));
end    

end % End of calcticks/getTextSize

