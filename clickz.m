function varargout = clickz(varargin)
% clickz temporarily prints z values corresponding to clicked points
% on a surface or image.  If multiple surfaces or images exist in the
% same axes, clickz first looks for surfaces and determines the one on the 
% top of the graphical stack is the relevant one. If no surfaces exist, 
% but an image or multiple images exist, clickz probes the image on the
% top of the graphical stack. 
% 
% To find "z" values of an existing image or pcolor plot, simply type
% clickz and start clicking on areas of interest. Instead of clicking, 
% you may also hit any key on the keyboard except the keys listed below
% which perform the following functions:
% 
%   Carriage return   Terminates data entry.
%            + or z   Zooms in, centered on current cursor location. 
%            - or x   Zooms out, centered on current cursor location.  
% 
% Above I put "z" in quotation marks for a reason. This function does 
% not actually return plotted z values, but it probes the color data of surfaces
% and images and assumes color corresponds to z values. This is often
% true, but not always. For example, you may have some x,y,z data plotted
% as a surface and you're letting cdata represent some other variable. Accordingly, 
% this function is probably best used with 2D pcolor or image plots, although it 
% can sometimes be used with 3D surface plots. 
% 
% Occasionally, linear interpolation between plotted data points fails,
% such as when xdata and ydata are not perfectly monotonic and plaid. I've
% only run into this problem when plotting polar stereographic data
% referenced one latitude, onto a polar stereographic map referenced to
% another latitude.  If for some reason linear interpolation fails, clickz
% will attempt to find the nearest x,y data points by Euclidean distance,
% and will print z data at that nearby location.  In such a case, any output x,y,z
% data will reflect the nearby point instead of the clicked point. 
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
% clickz
% clickz(N)
% clickz(...,ax)
% clickz(...,'keep') 
% clickz(,...'TextProperty',TextValue)
% z = clickz(...)
% [x,y,z] = clickz(...)
% [x,y,z,h] = clickz(...)
%
%% Description
% 
% clickz temporarily prints a z value on a plot at each click location.
% Previous printed values are deleted with each new click. clickz continues to
% run until user hits Return on the keyboard. 
%
% clickz(N) performs clickz N times or until the user hits Return on the
% keyboard. If N is specified with 'keep' and/or text formatting, N must 
% be the first input argument. 
% 
% clickz(...,ax) specifies an axis handle ax on which to use clickz. 
%
% clickz(...,'keep') does not delete printed points. 
%
% clickz(,...'TextProperty',TextValue) formats printed text. Multiple text
% properties and values may be specified, including fontsize, color,
% background, etc. 
%
% z = clickz(...) returns an array of clicked z values. 
%
% [x,y,z] = clickz(...) returns clicked x, y, and z values. 
% 
% [x,y,z,h] = clickz(...) also returns text object handle h when the 'keep'
% command is used and all four outputs are requested by the user. 
% 
%% Examples
% Any of the inputs described above can be combined or used individually. For
% example, after plotting some data like this: 
% 
%   pcolor(peaks);
%   colormap(jet(256))
%   colorbar
%   shading interp
% 
% you may say, "Wow, what's the value at that interesting light blue spot near
% (x=14,y=26)?" To answer this, you might look at the colorbar and your eyes
% will dart back and forth between the image and the colorbar, trying to
% narrow in on where exactly the light blue of (x=14,y=26) lies on the
% colorbar.  
% 
% Then, if you're like me, you'll start second-guessing your ability to 
% distinguish between colors on that smooth gradient. So you might try
% clicking the Data Cursor tool in the Figure Window, which will in turn
% tell you that X:14, Y:26, Z:0. That's because Matlab hasn't plotted this
% colormap in vertical z space--the whole surface lies in the z=0 plane. 
% 
% The next line trick I sometimes try is to start setting caxis
% limits. By inspection I can see that the light blue lies somewhere
% between -4 and -1, so I might try caxis([-4 -1]) and then I see that the 
% interesting point is probably between -3.5 and -2.5, and I can keep
% narrowing the caxes from there to hone in on the actual value by
% inspecting the colors.  Or, I could simply type
% 
%    clickz
% 
% and click on the interesting spot. 
% 
% If you know you want 5 spots before you start clicking, type this: 
% 
%    clickz(5)
% 
% If you want the text to be blue, this does the trick: 
% 
%    clickz('color','blue')
% 
% If you don't want the printed values to disappear, use the 'keep'
% command:
% 
%    clickz(3,'color','blue','fontweight','bold','keep')
% 
% The command above is how I created the example image for this function.
% To return arrays of clicked x,y,z values, tell clickz that you want data
% like this: 
% 
%    [x,y,z] = clickz;
% 
%
%% Author Info
% 
% This function is a modification of Matlab version 2012b's ginput function.
% clickz was written by Chad A. Greene of the University of 
% Texas at Austin's Institute for Geophysics (UTIG) in December 2014.  
% Come see me sometime at http://www.chadagreene.com. 
% 
% Where xdata or ydata include NaNs, John D'Errico's inpaint_nans function
% is employed. 
% 
% See also GINPUT, GTEXT, WAITFORBUTTONPRESS.


% Set defaults: 
N = []; 
keepPoints = false; 
fig = gcf;
figure(fig);
axes_handle = gca;


% Set axes if user declares anything other than gca: 
if nargin>0
    tmp = isax(varargin{:}); 
    if any(tmp) 
        axes_handle = varargin{find(tmp)}; 
        varargin = varargin(~tmp); 
    end
end

% Set number of clicks if user wants to limit N: 
if ~isempty(varargin)
    if isscalar(varargin{1})
        N = varargin{1}; 
        varargin(1)=[]; 
    end
end


% determine whether labels should be kept on the figure
tmp = strcmpi(varargin,'keep');
if any(tmp) 
    keepPoints = true; 
    varargin = varargin(~tmp); 
end

x = []; y = []; 
c = computer;
if ~strcmp(c(1:2),'PC')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
    if ~isempty(N)
        x = trmginput(N);
    else
        x = trmginput;
    end
else
    

InitialAxes = axis; 
imagecase = false; % Only look for images if no surfaces are found.

axch = findobj(axes_handle,'type','surface'); 
if isempty(axch) 
    axch = findobj(axes_handle,'type','image'); 
    imagecase = true; 
end
assert(isempty(axch)==0,'Cannot find any surface or image objects on the current axes.')
xdata = get(axch(1),'xdata'); 
ydata = get(axch(1),'ydata'); 
zdata = get(axch(1),'cdata'); % <- Notice we get cdata and call it zdata. 


% Now, if xdata or ydata have any nans in them, interpolation will not
% work. So we'll use John D'Errico's wonderful inpaint_nans function if
% necessary.  Unfortunately, the output of inpaint nans may differ from
% strictly monotonically-increasing, equally-spaced values by a few eps. 
% A crude way of getting around that is to conver to single, then back
% to double precision: 
if sum(isnan(xdata(:)))>0
    xdata = double(single(inpaint_nans(xdata))); 
end

if sum(isnan(ydata(:)))>0
    ydata = double(single(inpaint_nans(ydata))); 
end

% If the best layer is an image, convert the limits of xdata and ydata to arrays of pixel values: 
if imagecase && numel(xdata)==2
    xdata = xdata(1):xdata(end); 
    ydata = ydata(1):ydata(end); 
end

if isempty(N)
    how_many = -1; % an oddly-named variable from ginput, how_many means if the user specified a number of points to enter, how many of those points still need to be entered? how_many counts down from N.  
else
    how_many = N;
    if  ischar(how_many) ...
            || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
            || ~(fix(how_many) == how_many) ...
            || how_many < 0
        error(message('MATLAB:clickz:NeedPositiveInt'))
    end
    if how_many == 0
        % If input argument is equal to zero points,
        % give a warning and return empty for the outputs.

        warning (message('MATLAB:clickz:InputArgumentZero'));
    end
end

% Setup the figure to disable interactive modes and activate pointers. 
initialState = setupFcn(fig);

% onCleanup object to restore everything to original state in event of
% completion, closing of figure errors or ctrl+c. 
c = onCleanup(@() restoreFcn(initialState));


% We need to pump the event queue on unix
% before calling WAITFORBUTTONPRESS
drawnow
char = 0;

while how_many ~= 0
    % Use no-side effect WAITFORBUTTONPRESS
    waserr = 0;
    try
        keydown = wfbp;
    catch %#ok<CTCH>
        waserr = 1;
    end
    if(waserr == 1)
        if(ishghandle(fig))
            cleanup(c);
            error(message('MATLAB:clickz:Interrupted'));
        else
            cleanup(c);
            error(message('MATLAB:clickz:FigureDeletionPause'));
        end
    end
    % g467403 - clickz failed to discern clicks/keypresses on the figure it was
    % registered to operate on and any other open figures whose handle
    % visibility were set to off
    figchildren = allchild(0);
    if ~isempty(figchildren)
        ptr_fig = figchildren(1);
    else
        error(message('MATLAB:clickz:FigureUnavailable'));
    end
    %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
    %         clicked figure has handlevisibility set to callback
    if(ptr_fig == fig)
        if keydown
            char = get(fig, 'CurrentCharacter');
        end
        if ismember(char,[122 90 43 61]) % 122 and 90 are z and Z, 43 and 61 are + and =.  
            lim = axis; % get current axis limits
            pti=get(axes_handle, 'CurrentPoint');  % get current cursor point
            axis([pti(1,1)+diff(lim(1:2))/2*[-1 1] pti(1,2)+diff(lim(3:4))/2*[-1 1]]); % reset axis limits
            zoom(2) % zoom in

        elseif ismember(char,[120 88 45]) % lowercase x is 120, uppercase is 88 minus is 45
            lim = axis; % get current axis limits
            pti=get(axes_handle, 'CurrentPoint');  % get current cursor point
            axis([pti(1,1)+diff(lim(1:2))/2*[-1 1] pti(1,2)+diff(lim(3:4))/2*[-1 1]]); % reset axis limits
            zoom(.5)

        else
        drawnow;


        pt = get(axes_handle, 'CurrentPoint');


        if(char == 13) % & how_many ~= 0)
            % if the return key was pressed, char will == 13,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.
            % If this was an early breakout, don't include
            % the <Return> key info in the return arrays.
            % We will no longer count it if it's the last input.
            break;
        end


        % Get current x/y click points: 
        x = [x;pt(1,1)]; %#ok<AGROW>
        y = [y;pt(1,2)]; %#ok<AGROW>
        
        
        if x(end)<InitialAxes(1) || x(end)>InitialAxes(2) || y(end)<InitialAxes(3) || y(end)>InitialAxes(4)
            % do nothing if clicks are outside axis range. 
        else
            how_many = how_many-1; 

            try
                % Linearly interpolate z at those click points:     
                z(length(x)) = interp2(xdata,ydata,zdata,x(end),y(end)); 
            catch
                % If interpolation didn't work, return the nearest data point:
                nearInd = near2(xdata,ydata,x(end),y(end));
                x(end) = xdata(nearInd); 
                y(end) = ydata(nearInd); 
                z(length(x)) = zdata(nearInd); 
                
            end

            % Place text 
            if keepPoints
                hp(length(z)) = text(x(end),y(end),num2str(z(end)),'vert','middle','horiz','center',varargin{:}); 
            else
                try; delete(hp); end
                hp = text(x(end),y(end),num2str(z(end)),'vert','middle','horiz','center',varargin{:}); 
            end
            set(fig,'Pointer','crosshair')
        end

        end
        char = []; 
    end
end

% Cleanup and Restore 
cleanup(c);


%% Outputs 

% columnate: 
x = x(:); 
y = y(:); 
z = z(:); 

if ~keepPoints
    delete(hp)
end

axis(InitialAxes)

switch nargout
    case 1
        varargout{1} = z; 
    case 3
        varargout{1} = x; 
        varargout{2} = y; 
        varargout{3} = z; 
    case 4
        varargout{1} = x; 
        varargout{2} = y; 
        varargout{3} = z; 
        varargout{4} = hp; 
        
    case 0 
        varargout{1} = z; 
    otherwise
        error('Wrong number of output arguments. You can request only z, or x,y,z, or x,y,z,h--no other combination.') 
end
        
  
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok<NASGU>

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error(message('MATLAB:clickz:Interrupted'));
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function initialState = setupFcn(fig)

% Store Figure Handle. 
initialState.figureHandle = fig; 

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

% Setup FullCrosshair Pointer without warning. 
oldwarnstate = warning('off', 'MATLAB:hg:Figure:Pointer');

set(fig,'Pointer','crosshair');

warning(oldwarnstate);

% Adding this to enable automatic updating of currentpoint on the figure 
set(fig,'WindowButtonMotionFcn',@(o,e) dummy());

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
end

function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    set(initialState.figureHandle,'WindowButtonMotionFcn','');
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
end

function dummy()
% do nothing, this is there to update the clickz WindowButtonMotionFcn. 
end

function cleanup(c)
if isvalid(c)
    delete(c);
end
end








function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
  method = 0;
elseif ~ismember(method,0:5)
  error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
 case 0
  % The same as method == 1, except only work on those
  % elements which are NaN, or at least touch a NaN.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really a 1-d case
    work_list = nan_list(:,1);
    work_list = unique([work_list;work_list - 1;work_list + 1]);
    work_list(work_list <= 1) = [];
    work_list(work_list >= nm) = [];
    nw = numel(work_list);
    
    u = (1:nw)';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
      repmat([1 -2 1],nw,1),nw,nm);
  else
    % a 2-d case
    
    % horizontal and vertical neighbors only
    talks_to = [-1 0;0 -1;1 0;0 1];
    neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
    
    % list of all nodes we have identified
    all_list=[nan_list;neighbors_list];
    
    % generate sparse array with second partials on row
    % variable for each element in either list, but only
    % for those nodes which have a row index > 1 or < n
    L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    else
      fda=spalloc(n*m,n*m,size(all_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    end
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 1
  % least squares approach with del^2. Build system
  % for every array element as an unknown, and then
  % eliminate those which are knowns.

  % Build sparse matrix approximating del^2 for
  % every element in A.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % a 1-d case
    u = (1:(nm-2))';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
      repmat([1 -2 1],nm-2,1),nm-2,nm);
  else
    % a 2-d case
    
    % Compute finite difference for second partials
    % on row variable first
    [i,j]=ndgrid(2:(n-1),1:m);
    ind=i(:)+(j(:)-1)*n;
    np=(n-2)*m;
    fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      repmat([1 -2 1],np,1),n*m,n*m);
    
    % now second partials on column variable
    [i,j]=ndgrid(1:n,2:(m-1));
    ind=i(:)+(j(:)-1)*n;
    np=n*(m-2);
    fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
      repmat([1 -2 1],np,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 2
  % Direct solve for del^2 BVP across holes

  % generate sparse array with second partials on row
  % variable for each nan element, only for those nodes
  % which have a row index > 1 or < n
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really just a 1-d case
    error('Method 2 has problems for vector input. Please use another method.')
    
  else
    % a 2-d case
    L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    else
      fda=spalloc(n*m,n*m,size(nan_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    end
    
    % fix boundary conditions at extreme corners
    % of the array in case there were nans there
    if ismember(1,nan_list(:,1))
      fda(1,[1 2 n+1])=[-2 1 1];
    end
    if ismember(n,nan_list(:,1))
      fda(n,[n, n-1,n+n])=[-2 1 1];
    end
    if ismember(nm-n+1,nan_list(:,1))
      fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
    end
    if ismember(nm,nan_list(:,1))
      fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
    end
    
    % eliminate knowns
    rhs=-fda(:,known_list)*A(known_list);
    
    % and solve...
    B=A;
    k=nan_list(:,1);
    B(k)=fda(k,k)\rhs(k);
    
  end
  
 case 3
  % The same as method == 0, except uses del^4 as the
  % interpolating operator.
  
  % del^4 template of neighbors
  talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
      0 1;0 2;1 -1;1 0;1 1;2 0];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with del^4, but only
  % for those nodes which have a row & column index
  % >= 3 or <= n-2
  L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
  nl=length(L);
  if nl>0
    % do the entire template at once
    fda=sparse(repmat(all_list(L,1),1,13), ...
        repmat(all_list(L,1),1,13) + ...
        repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
        repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % on the boundaries, reduce the order around the edges
  L = find((((all_list(:,2) == 2) | ...
             (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
           (((all_list(:,3) == 2) | ...
             (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,5), ...
      repmat(all_list(L,1),1,5) + ...
        repmat([-n,-1,0,+1,n],nl,1), ...
      repmat([1 1 -4 1 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,2) == 1) | ...
             (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-n,0,n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,3) == 1) | ...
             (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-1,0,1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 4
  % Spring analogy
  % interpolating operator.
  
  % list of all springs between a node and a horizontal
  % or vertical neighbor
  hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for i=1:4
    hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
  
 case 5
  % Average of 8 nearest neighbors
  
  % generate sparse array to average 8 nearest neighbors
  % for each nan element, be careful around edges
  fda=spalloc(n*m,n*m,size(nan_list,1)*9);
  
  % -1,-1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,-1
  L = find(nan_list(:,3) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,-1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,0
  L = find(nan_list(:,2) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,0
  L = find(nan_list(:,2) < n);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,+1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,+1
  L = find(nan_list(:,3) < m);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,+1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);
end

% ====================================================
%      end of main function
% ====================================================
% ====================================================
%      begin subfunctions
% ====================================================
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for i=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(i,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end

end

function [rowi,coli] = near2(X,Y,xi,yi)
%near2 returns indices of values in X and Y that are close to some x,y point. It's similar to 
%find, for nearest neighbors, on a 2D grid. 
%
%% Syntax
% 
%  ind = near2(X,Y,xi,yi)
%  [rowi,coli] = near2(X,Y,xi,yi)
% 
%% Description
% 
% ind = near2(X,Y,xi,yi) returns a linear index corresponding to values
% in X and Y that are closest to the point given by ( xi, yi ). X
% and Y must be 2D grids of equal size and xi and yi must be scalar. 
%
% [rowi,coli] = near2(X,Y,xi,yi) returns the row and column indices 

%% Examples
% % First create a grid of X and Y data:
% 
% [X,Y] = meshgrid(1:5,10:-1:5)
% 
% % Now we can look for X, Y indices corresponding to the point closest
% to (2.1,9.724): 
% 
% [rowi,coli] = near2(X,Y,2.1,9.724)
%
% % Or if you'd prefer linear indices, 
% 
% ind = near2(X,Y,2.1,9.724)
% 
% % If xi,yi describes some point equidistant between multiple grid points
% % in X and Y, all valid minimum-distance indices are returned: 
% 
% [rowi,coli] = near2(X,Y,2.5,9.1)
% 
%% Author Info
% 
% Chad A. Greene of the University of Texas at Austin's Institute for Geophysics
% wrote this in October of 2014.  You can visit Chad over at his internet
% website, www.chadagreene.com.  
% 
% See also find, interp2, ind2sub, sub2ind, and hypot. 

%% Some input checks:

assert(nargin==4,'near2 requires exactly four inputs: gridded X and Y, and scalar values for xi and yi.')
assert(numel(X)==numel(Y),'X and Y must be the same size.') 
assert(size(X,1)==size(Y,1),'X and Y must be the same size.') 
assert(size(X,2)>1,'X and Y must be 2D grids of equal size.')
assert(isscalar(xi)==1,'xi must be a scalar.') 
assert(isscalar(yi)==1,'yi must be a scalar.') 

% extrapolation warnings: 
if xi<min(X(:)) || xi>max(X(:))
    warning('xi lies outside the range of X.')
end
if yi<min(Y(:)) || yi>max(Y(:))
    warning('yi lies outside the range of Y.')
end

%% Calculate distance between (xi,yi) and all points in (X,Y) 

dst = hypot(X-xi,Y-yi);

%% Return index of shortest distance: 

if nargout<2
    % if only one output, return linear index: 
    rowi = find(dst==min(dst(:))); 
else
    [rowi,coli] = find(dst==min(dst(:))); 
end

end


function [tf] = isax(varargin)
% isax determines whether inputs are axis handles.
% 
%% Syntax & Description 
% 
% isax(H) returns an array whose elements are 1 where the elements of
%         H are valid axis handles, and 0 where they are not.
% 
%% Author Info
% This function was written by Chad A. Greene of the University of 
% Texas at Austin's Institute for Geophysics (UTIG). December 2014. 
% Come say hi at http://www.chadagreene.com. 
% 
% See also gca, axes, axis, ishandle, ishghandle, and isa.

fax = findobj('type','axes');
if nargin>1
    tf = zeros(size(varargin)); 
    for k = 1:length(varargin) 
        try
            if ismember(varargin{k},fax)
                tf(k)=1; 
            end
        end
    end
    else
tf = ismember(varargin{:},fax); 
end

end


    

