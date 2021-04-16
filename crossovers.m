function varargout = crossovers(varargin)
% crossovers efficiently finds crossover locations of georeferenced or polar
% stereographic coordinates. This function works for self-intesecting flight 
% lines or multiple flight lines concatenated into a single array. 
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
%  [xi,yi] = crossovers(x,y) 
%  [lati,loni] = crossovers(lat,lon)
%  [...,zi,ti] = crossovers(...,z,t)
%  [...] = crossovers(...,'maxdt',max_time_diff)
%  [...] = crossovers(...,'mindt',min_time_diff)
%  [...] = crossovers(...,'maxdist',max_distance)
%  [...] = crossovers(...,'inrange',latlim_or_xlim,lonlim_or_ylim)
%  [...] = crossovers(...,'tile',TilePreference)
%  [...] = crossovers(...,'interp',InterpolationMethod)
% 
%% Description
% 
% [xi,yi] = crossovers(x,y) returns location(s) of crossover points in the arrays x,y, 
% where x and y are polar stereographic (ref -71) coordinates. x and y must be column 
% vectors.  
%
% [lati,loni] = crossovers(lat,lon) returns georeferenced location(s) of crossovers
% if inputs are geo coordinates. Geo coordinates are assumed if no
% element in the first argument has an absolute value greater than 90 and no element  
% in the second argument has an absolute value greater than 360. lat and lon must 
% be column vectors. 
% 
% [...,zi,ti] = crossovers(...,z,t) returns elevations and times of crossovers. zi and ti 
% are both two-column vectors whose first column corresponds to the elevation and time of
% the first measurement in time, and column 2 corresponds to the second pass. Elevation 
% change at crossovers dzi/dti may be calculated as diff(zi,1,2)./diff(ti,1,2).  
%
% [...] = crossovers(...,'maxdt',max_time_diff) limits results to crossovers separated in time 
% by less than or equal to max_time_diff. Using a max_time_diff value of 150 is one way
% of considering only crossovers from the same field campaign year.  Units of max_time_diff
% are days. 
% 
% [...] = crossovers(...,'mindt',min_time_diff) limits results to crossovers separated in time
% by at least min_time_diff. Using min_time_diff = 1/96 (15 minutes, that is) is one way
% to prevent the intersection-finding function from getting confused. Units of min_time_diff
% are days. 
% 
% [...] = crossovers(...,'maxdist',max_distance) limits results to crossovers that are within 
% max_distance meters of the nearest measurement in both passes. Using max_distance = 100 
% ensures that every value in (xi,yi) or (lati,loni) lies within 100 meters of a measurement 
% from pass 1 and pass 2.  
% 
% [...] = crossovers(...,'inrange',latlim_or_xlim,lonlim_or_ylim) limits analysis to a subset
% geoquad or polar stereographic bounding box. If input locations are geo coordinates, bounding
% limits are assumed to be geo coordinates.  If inputs coordinates are polar stereographic meters, 
% xlim and ylim are assumed.  For example, crossovers(lat,lon,'inrange',[-70 -65],[110 115]) 
% subsets lat,lon to the geoquad given by 70S to 65S and 110E to 115E.   
% 
% [...] = crossovers(...,'tile',TilePreference) specifies a solver setting as 'on', 'off', or 
% 'auto'. For inputs of more than 15,000 points, calculating path intersections all at once with
% tiling off crashes my computer. This function is designed to intelligently solve sections at
% a time by breaking the domain into multiple tiles.  Default TilePreference is 'auto'.   
% 
% [...] = crossovers(...,'interp',InterpolationMethod) specifies method of interpolation for  
% zi and ti. Default method is 'linear', but you may also choose 'cubic'. Cubic interpolation
% requires more points in each direction and it's very slightly slower than linear.  
% 
%% Example: Plot data
% % Consider this self-intersecting flight line of ~9000 data points:
% 
% load samplegrid 
% modismoa(mean(lat),mean(lon),50,'contrast','lc')
% plotm(lat,lon,'.','color',rgb('light blue'))
% scalebar('length',10)
% 
% % To find all crossover locations the syntax is easy: 
% 
% [latc,lonc] = crossovers(lat,lon); 
% plotm(latc,lonc,'o','color',rgb('dark blue'))
% 
% % Sometimes you may want to specify conditions for crossovers.  Often times
% % we want to filter out processes that take place on short or long time
% % scales.  Here let's say we only want crossovers where the different legs 
% % of the flight were more than one hour apart and less than 2.5 hours
% % apart.  Further, we want to make sure that each crossover location lies
% % within 45 meters of a measurement on _each_ leg of the crossover: 
% 
% [latc,lonc,zc] = crossovers(lat,lon,z,t,...
%     'mindt',1/24,...   % crossovers separated by > 1 hr
%     'maxdt',2.5/24,... % crossovers separated by < 2.5 hrs
%     'maxdist',45,...   % within 45 m of a measurement   
%     'interp','cubic'); % use cubic interpolation
% 
% plotm(latc,lonc,'gx')
% 
%% Author Info: 
% Chad A. Greene, The University of Texas at Austin's Institute for 
% Geophysics, June 2015. This function includes Douglas M. Schwarz's 
% intersections function as a subfunction.  
% 
% June 26, 2015: Changed this input parse: 
% tmp = strncmpi(varargin,'inrange',2);  to  tmp = strncmpi(varargin,'inrange',3); 
% because 'interp' as an input was being misunderstood as 'inrange'.  
% 
% See also pathcrossingsps71. 

%% Error checks: 

narginchk(2,inf) 
nargoutchk(2,4)
assert(isnumeric(varargin{1})==1,'crossovers input 1 must be numeric lat or x.') 
assert(isnumeric(varargin{2})==1,'crossovers input 2 must be numeric lon or y.') 
assert(exist('ll2ps','file')==2,'Cannot find the necessary function ll2ps, which is a function in Antarctic Mapping Tools.  You can find AMT here: http://www.mathworks.com/matlabcentral/fileexchange/47638.') 
assert(license('test','map_toolbox')==1,'Crossovers function requires Matlab''s Mapping Toolbox.') 
assert(numel(varargin{1})==numel(varargin{2}),'Dimensions of first two inputs must match.') 

%% Set defaults: 

maxdt = inf; 
mindt = 0; 
maxdist = inf; 
TilePref = 'auto'; 
interptype = 'linear'; 

%% Parse inputs: 

if islatlon(varargin{1},varargin{2})
    [x,y] = ll2ps(varargin{1},varargin{2}); 
    outputgeo = true; 
else
    x = varargin{1};
    y = varargin{2}; 
    outputgeo = false; 
end

% Columnate if necessary: 
if isrow(x) 
    x = x'; 
    y = y'; 
end

% Define dummy z and t in case user does not: 
t = zeros(size(x)); 
z = zeros(size(x)); 


if nargin>2
    if isnumeric(varargin{3})
        assert(nargin>3 & isnumeric(varargin{4})==1,'Having trouble parsing inputs to the crossovers function. The third input is numeric, (presumably z), meaning you''ll also have to enter an input for time.') 
        assert(numel(varargin{3})==numel(varargin{1}),'Dimension of z vector must match location vector.') 
        assert(numel(varargin{4})==numel(varargin{1}),'Dimension of time vector must match location vector.') 
        z = varargin{3}; 
        t = varargin{4}; 
    end 
    
    % Define maximum time separation between crossovers: 
    tmp = strcmpi(varargin,'maxdt'); 
    if any(tmp)
        maxdt = varargin{find(tmp)+1}; 
    end
    
    % Define minimum time separation between crossovers:  
    tmp = strcmpi(varargin,'mindt'); 
    if any(tmp)
        mindt = varargin{find(tmp)+1}; 
    end
    
    % Define maximum distance a crossover can be from data in each direction:  
    tmp = strcmpi(varargin,'maxdist'); 
    if any(tmp)
        maxdist = varargin{find(tmp)+1}; 
    end    
    
    % Cull data to a user-specified geoquad or xy range:  
    tmp = strncmpi(varargin,'inrange',3); 
    if any(tmp) 
        lim1 = varargin{find(tmp)+1}; 
        lim2 = varargin{find(tmp)+2}; 
        assert(numel(lim1)==2,'Input latlim or xlim must be a two element vector in the form [min max].')
        assert(numel(lim2)==2,'Input lonlim or ylim must be a two-element vector in the form [min max].')
        ind = varargin{1}>=lim1(1) & varargin{1} <= lim1(2) & varargin{2}>=lim2(1) & varargin{2} <= lim2(2); 
        x = x(ind); 
        y = y(ind); 
        z = z(ind); 
        t = t(ind); 
    end
    
    % Solve by several tiled bits or solve everything in one go? 
    tmp = strncmpi(varargin,'tile',3); 
    if any(tmp) 
        TilePref = varargin{find(tmp)+1}; 
    end
        
    % Interpolation method: 
    tmp = strncmpi(varargin,'interp',6);
    if any(tmp)
        interptype = lower(varargin{find(tmp)+1}); 
    end
end


%% Status message for large datasets: 

if numel(x)>100e3
    fprintf(['Solving for ',num2str(numel(x)),' points. I find that 400,000 points takes about a minute to process \n',...
        'on my Mac Mini from 2012. Processing time scales somewhat linearly with number of points,\n',...
        'but can vary depending on how clustered the data points are. The best way to improve processing\n',...
        'time is to limit geographic region using the ''inrange'' option described in the crossover \n',...
        'function header. \n'])
end

%% Determine whether to solve by tiling or all in one go: 

switch lower(TilePref) 
    case 'auto' 
        if numel(x)>5e3
            tilesoln = true; 
        else
            tilesoln = false; 
        end
        
    case 'on'
        tilesoln = true; 
        
    case 'off'
        tilesoln = false; 
        
    otherwise
        error('Unrecognized tiling preference. Check inputs to crossovers function.') 
end


%% Find all possible crossover points: 

if tilesoln
    
    % For large number of points, tile the solver to prevent memory problems:  
    xgrid = linspace(min(x),max(x),ceil(numel(x)/5e3)); 
    ygrid = linspace(min(y),max(y),ceil(numel(x)/5e3)); 

    xi = []; 
    yi = []; 
    I = []; 
    J = []; 
    for kx = 1:length(xgrid)-1
        for ky = 1:length(ygrid)-1
            
            % Get all input points within tile plus a 1 km buffer:  
            ind = x>= xgrid(kx)-1e3 & x<= xgrid(kx+1)+1e3 & y>=ygrid(ky)-1e3 & y<=ygrid(ky+1)+1e3;
            
            % Only solve for intersections of data if data exist in this tile: 
            if sum(ind)>3
                
                % Temporarily log indices so we can properly reference I,J to x,y later:   
                tmpind = find(ind); 
                
                [xk,yk,Ik,Jk] = intersections(x(ind),y(ind)); 
                
                % Sometimes lines run away and never intersect. I think keeping just the finite ones'll fix it:   
                isf = isfinite(Ik); 
                xk = xk(isf); 
                yk = yk(isf); 
                Ik = Ik(isf); 
                Jk = Jk(isf); 
                
                xi = [xi;xk]; 
                yi = [yi;yk];
                
                I = [I;tmpind(floor(Ik)) + rem(Ik,1)]; 
                J = [J;tmpind(floor(Jk)) + rem(Jk,1)]; 
            end
        end
    end
    
    % Delete duplicate data: 
    C = unique([xi yi I J],'rows'); 
    xi = C(:,1); 
    yi = C(:,2); 
    I = C(:,3); 
    J = C(:,4); 
    
else
    % Do it all at once if not tiled: 
    [xi,yi,I,J] = intersections(x,y); 
    
end

%% Cull crossover points to only include points within user-specified distance of measurements and temporal separation:   

% Spatial distance from (xi,yi) to nearest data point in first crossover:   
di = hypot(xi-x(round(I)),yi-y(round(I)));

% Spatial distance from (xi,yi) to nearest data point in second crossover:  
dj = hypot(xi-x(round(J)),yi-y(round(J)));

% Temporal separation between crossovers:  
dt = abs(t(round(I)) - t(round(J))); 

% Get measurement times on each side of ti for each pass (if it's hours between measurements, that's probably different flights):   
dti = abs(t(ceil(I))-t(floor(I))); 
dtj = abs(t(ceil(J))-t(floor(J))); 
sameflightdt = 30*median(diff(t(isfinite(t)))); % assumes any long temporal separation between data points are different flights, and thererfore not really crossovers.    

% Indices of crossover points meeting spatial and temporal criteria:  
ind = dt>=mindt & dt<=maxdt & di<=maxdist & dj<=maxdist & dti<=sameflightdt & dtj<=sameflightdt;


% Delete all data that do not meet spatio-temporal criteria:
xi = xi(ind); 
yi = yi(ind); 

%% Calculate crossover times and elevations: 

if nargout>2 
    I = I(ind); 
    J = J(ind); 

    switch interptype
        case {'linear','lin'}
            % Column 1 corresponds to pass I; column 2 is pass J: 
            zI = z(floor(I)) + rem(I,1).*(z(ceil(I))-z(floor(I))); 
            zJ = z(floor(J)) + rem(J,1).*(z(ceil(J))-z(floor(J))); 
            tI = t(floor(I)) + rem(I,1).*(t(ceil(I))-t(floor(I))); 
            tJ = t(floor(J)) + rem(J,1).*(t(ceil(J))-t(floor(J))); 
            
        case {'cubic','cub'}
            zI = NaN(size(I)); 
            zJ = zI;     
            tI = NaN(size(I)); 
            tJ = tI; 

            for k = 1:length(I)
                if abs(t(round(I(k)-2)) - t(round(I(k)+2)))< mindt/2 && abs(t(round(J(k)-2)) - t(round(J(k)+2)))< mindt/2
                try
                zI(k) = interp1(round(I(k))+(-2:2),z(round(I(k))+(-2:2)),I(k),'cubic'); 
                zJ(k) = interp1(round(J(k))+(-2:2),z(round(J(k))+(-2:2)),J(k),'cubic'); 
                tI(k) = interp1(round(I(k))+(-2:2),t(round(I(k))+(-2:2)),I(k),'cubic'); 
                tJ(k) = interp1(round(J(k))+(-2:2),t(round(J(k))+(-2:2)),J(k),'cubic'); 
                catch 
                    k
                end
                end
            end
            
        otherwise
            error('Unrecognized interpolation type.') 
    end
    
    
    zi = [zI zJ]; 
    varargout{3} = zi; 
    

    if nargout==4
        ti = [tI tJ]; 
        varargout{4} = ti; 
    end
end


%% Package results into proper outputs: 

if outputgeo
    [varargout{1},varargout{2}] = ps2ll(xi,yi); 
else
    varargout{1} = xi; 
    varargout{2} = yi; 
end


end






function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Adjustments when fewer than five arguments are supplied.
switch nargin
	case 2
		robust = true;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 3
		robust = x2;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 4
		robust = true;
		self_intersect = false;
	case 5
		self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
		length(x1) ~= length(y1)
	error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
		length(x2) ~= length(y2)
	error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
	repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
	repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
	repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
	repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
	repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Force i and j to be column vectors, even when their length is zero, i.e.,
% we want them to be 0-by-1 instead of 0-by-0.
i = reshape(i,[],1);
j = reshape(j,[],1);

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

if robust
	overlap = false(n,1);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		for k = 1:n
			T(:,k) = AA(:,:,k)\B(:,k);
			[~,last_warn] = lastwarn;
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
	catch err
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	% For overlapping segment pairs the algorithm will return an
	% intersection point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	[xy0,index] = unique(xy0,'rows');
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
		iout = i(sel) + T(1,sel).';
		jout = j(sel) + T(2,sel).';
	end
else % non-robust option
	for k = 1:n
		[L,U] = lu(AA(:,:,k));
		T(:,k) = U\(L\B(:,k));
	end
	
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		iout = i(in_range) + T(1,in_range).';
		jout = j(in_range) + T(2,in_range).';
	end
end

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');
end



function tf = islatlon(lat,lon)
% islatlon determines whether lat,lon is likely to represent geographical
% coordinates. 
% 
%% Syntax
% 
% tf = islatlon(lat,lon) returns true if all values in lat are numeric
% between -90 and 90 inclusive, and all values in lon are numeric between 
% -180 and 360 inclusive. 
% 
%% Example 1: A single location
% 
% islatlon(110,30)
%    = 0
% 
% because 110 is outside the bounds of latitude values. 
% 
%% Example 2: A grid
% 
% [lon,lat] = meshgrid(-180:180,90:-1:-90); 
% 
% islatlon(lat,lon)
%    = 1 
% 
% because all values in lat are between -90 and 90, and all values in lon
% are between -180 and 360.  What if it's really, really close? What if
% just one value violates these rules? 
% 
% lon(1) = -180.002; 
% 
% islatlon(lat,lon)
%    = 0
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG). http://www.chadagreene.com. 
% March 30, 2015. 
% 
% See also wrapTo180, wrapTo360, projfwd, and projinv.  

% Make sure there are two inputs: 
narginchk(2,2)

% Set default output: 
tf = true; 

%% If *any* inputs don't look like lat,lon, assume none of them are lat,lon. 

if ~isnumeric(lat)
    tf = false; 
    return
end

if ~isnumeric(lon)
    tf = false; 
    return
end
if any(abs(lat(:))>90)
    tf = false; 
    return
end

if any(lon(:)>360)
    tf = false; 
    return
end    

if any(lon(:)<-180)
    tf = false; 
end

end


