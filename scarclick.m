function [htext] = scarclick(varargin)
% scarclick labels features on a map of Antarctica using mouse clicks. 
% This function searches the Scienfitic Committee on Antarctic Research
% (SCAR) Composite Gazetteer database for feature names. Information 
% for citing this dataset is provided below. 
% 
%%  * * * SCARCLICK Operators * * *
% 
%       Mouse click   Labels nearest point. 
%   Carriage return   Terminates data entry.
%         Backspace   Deletes previous label. 
%                 +   Zooms in, centered on current cursor location. 
%                 -   Zooms out, centered on current cursor location.  
% 
% 
%% Syntax
% 
%  scarclick
%  scarclick('FeatureType') 
%  scarclick('keep') 
%  scarclick('TextProperty',TextValue)
%  h = scarclick(...)
%
%% Description 
% 
% scarclick labels whichever features are closest to where you click on a
% map. Each click deletes the previous label and places a new one. When 
% you're done clicking, hit return to exit scarclick. 
%  
% scarclick('FeatureType') limits labeled objects to a specific feature
% type.  Note that not all features have a type associated with them,
% so if you select 'ice stream' as a feature type and click directly
% on Kamb Ice Stream, a label will appear on the closest identified ice
% stream, Nimrod Ice Stream. Using scarlabel('kamb ice stream') works
% fine, but only three features are officially identified as ice streams,
% despite many more ice streams having "ice stream" in their names. 
% 
% The following are the most common FeatureTypes: 
%     'Anchorage'
%     'Archipelago'
%     'Bank'
%     'Basin'
%     'Bay'
%     'Beach'
%     'Bight'
%     'Bluff'
%     'Building'
%     'Butte'
%     'Buttress'
%     'Camp'
%     'Canyon'
%     'Cape'
%     'Channel'
%     'Cirque'
%     'Cliff'
%     'Coast'
%     'Cone'
%     'Corner'
%     'Cove'
%     'Crag'
%     'Crater'
%     'Crevasse'
%     'Dome'
%     'Escarpment'
%     'Fjord'
%     'Glacier'
%     'Gulf'
%     'Gully'
%     'Harbour'
%     'Head'
%     'Heights'
%     'Hill'
%     'Hillock'
%     'Ice field'
%     'Ice front'
%     'Ice rise'
%     'Ice shelf'
%     'Icefall'
%     'Inlet'
%     'Island'
%     'Knob'
%     'Knoll'
%     'Land'
%     'Massif'
%     'Mesa'
%     'Moraine'
%     'Mountain'
%     'Neve'
%     'Nunatak'
%     'Pass'
%     'Passage'
%     'Peak'
%     'Peninsula'
%     'Piedmont'
%     'Plain'
%     'Plateau'
%     'Platform'
%     'Point'
%     'Promontory'
%     'Range'
%     'Reef'
%     'Ridge'
%     'Rock'
%     'Saddle'
%     'Sea'
%     'Shoal'
%     'Slope'
%     'Snowfield'
%     'Sound'
%     'Spit'
%     'Spur'
%     'Stack'
%     'Station'
%     'Strait'
%     'Terrace'
%     'Tongue'
%     'Valley'
%     'Wall'
%     'Water body'
% 
% scarclick('keep') does not automatically delete labels. 
% 
% scarclick('TextProperty',TextValue) formats text labels with
% property-value pairs. (e.g., 'fontsize',30,'color','b')
% 
% h = scarclick(...) returns a handle h of text labels. 
%
%% Examples: 
% 
% First, initialize a map of the grounding line. If you don't have 
% the Bedmap2 Toolbox, you'll have to find some other way to make 
% a map of Antarctica:
% 
%    bedmap2 gl 
% 
% Then type clickz and start clicking away: 
%   
%    clickz   % (Hit return when you get bored with clicking.)
% 
% Don't automatically delete labels: 
%
%    clickz('keep') 
% 
% Combine commands if you'd like. Here we only label glaciers, using
% boldface red font, and don't automatically delete them. If you 
% accidentally click on something you do not want labeled, hit backspace
% to remove the last label: 
% 
%    clickz('glacier','keep','fontweight','bold','color','red')
% 
%% Known Issues: 
% 
% 1. This function is awfully slow when you start clicking around on a
% large dataset. If you've plotted a modismoa image and overlain measures
% ice speeds, scarclick will probably respond slowly. 
% 
% 2. When specifying a feature type, it may be tempting to write the plural
% form of the feature. That would certainly be more intuitive, but for now,
% you'll have to specify 'glacier', not 'glaciers'.  If you accidentally
% type a plural form, or the singular form of some feature that is not officially
% categorized by SCAR (i.e., 'Pizza Hut'), there will be no helpful error
% message. Sorry about that. 
% 
% 3. Not all features are associated with a FeatureType. 
% 
% 4. In some cases, the SCAR Composite Gazetteer contains multiple entries
% for the same feature, and sometimes the coordinates of the multiple
% entries do not match exactly. I have dealt with this issue by arbitrarily 
% picking one entry for each unique feature name and deleting duplicates. I did the
% same kind of arbitrary culling for a previous release of Antarctic Mapping 
% Tools, but I did not necessarily choose to keep the same entries  in both 
% cullings.  As a result, if you have used a previous version of scarloc
% or scarlabel, it is possible that a few features may appear to have moved in 
% this release. 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG), January 2015. Drop me a line
% if you'd like. If you do a minimal amount of digging, you can find my contact  
% info via http://www.chadagreene.com. 
% 
%% References
% 
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
% The SCAR database can be found here: http://www.scar.org/cga
% SCAR has kindly requested that Gazetteer data should be cited as:
% 
% Secretariat SCAR (1992, updated 2015). Composite Gazetteer of Antarctica,
% Scientific Committee on Antarctic Research. GCMD Metadata 
% http://gcmd.nasa.gov/records/SCAR_Gazetteer.html
% 
% See also inputm, ginput, gdistm, clickz, scarlabel, and scarloc. 

%% Error checks: 

assert(license('test','map_toolbox') ==1,'It looks like you do not have the Mapping Toolbox. As such, this will not work.') 
assert(ismap==1,'A map must be open to use scarclick.') 

%% Load SCAR data 

sn = load('scarnames.mat'); 

%% Set defaults: 


keepPoints = false; 
htext = [];

if nargin>0
    
    % determine whether labels should be kept on the figure
    tmp = strcmpi(varargin,'keep');
    if any(tmp) 
        keepPoints = true; 
        varargin = varargin(~tmp); 
    end
    
    % Is a feature type requested? 
    if ~isempty(varargin)
        typeDeclared = strcmpi(varargin{1},sn.featureType);
        if any(typeDeclared)
            sn.x = sn.x(typeDeclared); 
            sn.y = sn.y(typeDeclared); 
            sn.names = sn.names(typeDeclared); 
            sn.lat = sn.lat(typeDeclared); 
            sn.lon = sn.lon(typeDeclared); 
            varargin(1) = []; 
        end
    end
end

%% Begin work 


fig = gcf; 
axes_handle = gca;
how_many = -1; % This is a counter that is leftover from ginput. It never changes in the scarclick function. 
InitialAxes = axis; 


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
            
        elseif char==8
            try
                delete(htext(end));
                htext(end)=[];
            end

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
        [latclick,lonclick] = minvtran(pt(1,1),pt(1,2));
        [xclick,yclick] = ll2ps(latclick,lonclick); 

        % Get indices of nearest feature: 
        [~,nearind] = min(hypot(sn.x -xclick,sn.y-yclick));

        % Place text 
        if keepPoints
            htext(length(htext)+1) = textm(sn.lat(nearind),sn.lon(nearind),1,sn.names(nearind),'vert','middle','horiz','center'); 

        else
            try; delete(htext); end
            htext = textm(sn.lat(nearind),sn.lon(nearind),1,sn.names(nearind),'vert','middle','horiz','center'); 
        end
        
        if ~isempty(varargin)
            try
               set(htext,varargin{:})
            end
        end
        set(fig,'Pointer','crosshair')

        end
        char = []; 
    end
end

% Cleanup and Restore 
cleanup(c);


%% Clean up: 

% Reset initial axes if user zoomed in or out: 
axis(InitialAxes); 

% delete labels unless 'keep' was requested by user: 
if ~keepPoints
    delete(htext)
    clear htext
end

% delete text handles unless requested by user: 
if nargout==0
    clear htext
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

