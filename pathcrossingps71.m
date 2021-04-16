function [lati,loni] = pathcrossingps71(lat1,lon1,lat2,lon2,ClipOption)
%PATHCROSSINGPS71 returns crossing point(s) of two paths in or around
%Antarctica. Calculation is performed by converting input lat/lon arrays to
%polarstereo (-71) x,y coordinates, finding intersections, and then
%converting intersection point(s) back to georeferenced coordinates. 
%
%% Syntax 
% 
%  [lati,loni] = pathcrossingps71(lat1,lon1,lat2,lon2)
%  [lati,loni] = pathcrossingps71(lat1,lon1,lat2,lon2,ClipOption)
% 
%% Description 
% 
% [lati,loni] = pathcrossingps71(lat1,lon1,lat2,lon2) returns lat/lon pair(s) 
% for intersection of paths given by lat1,lon1 and lat2,lon2.  
% 
% [lati,loni] = pathcrossingps71(lat1,lon1,lat2,lon2,ClipOption) allows
% turning off of automatic data clipping for large data sets. Use the
% string 'noclip' to turn off clipping.  Clipping is on by default
% because calculating the intersection(s) of a million-point flight line
% with another multi-million point flight line would be very slow. 
% 
%% Example: Intersection of snowmobile track with grounding line
% Here we'll grab (very approximate) grounding line data from Bedmap2 and
% see where it intersects a snowmobile track given by four waypoints: 
% 
% [glat,glon] = bedmap2_data('gl'); 
% tracklat = [-75.4491  -75.4657  -75.3771  -75.1145];
% tracklon = [-94.4371  -95.4879  -97.8617 -101.0862]; 
% 
% %% 
% % For context, let's plot these data: 
% 
% modismoa('pine island glacier',350,'inset','northeast')
% bedmap2('gl','color','green','linewidth',2)
% plotm(tracklat,tracklon,'ro-','linewidth',2)
% 
% %% 
% % Now find the intersection and plot it as a big blue star: 
% 
% [intersectlat,intersectlon] = pathcrossingps71(glat,glon,tracklat,tracklon); 
% plotm(intersectlat,intersectlon,'bp','markersize',20,'markerfacecolor','blue')
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
% Chad A. Greene of the University of Texas at Austin's Intitute for
% Geophysics (UTIG).  August 2014.  This function includes the InterX function 
% described here: http://www.mathworks.com/matlabcentral/fileexchange/22441


%% Input checking: 

assert(isvector(lat1)==1,'Input lat1 must be a vector.')
assert(numel(lat1)==numel(lon1),'Input lat1 and lon1 must be the same size.') 
assert(isvector(lat2)==1,'Input lat2 must be a vector.')
assert(numel(lat2)==numel(lon2),'Input lat2 and lon2 must be the same size.') 

ClipData = true; 
if exist('ClipOption','var')==1
    if strncmpi(ClipOption,'no',2)||strcmpi(ClipOption,'off')
        ClipData = false;
    end
end

%% Begin work: 

% Columnate: 
lat1 = lat1(:); 
lon1 = lon1(:); 
lat2 = lat2(:); 
lon2 = lon2(:); 

% Transform to polar stereo coordinates with standard parallel at 71 S: 
[x1,y1] = ll2ps(lat1,lon1); 
[x2,y2] = ll2ps(lat2,lon2); 


%% Delete faraway points before performing InterX function for large data sets. 

if ClipData
    if numel(x1)*numel(x2)>1e6
        % First find standard deviations of each array so we'll know which outliers
        % we can delete: 
        for k = 1:2; 
        stdx1 = std(x1(2:end)-x1(1:end-1)); 
        stdy1 = std(y1(2:end)-y1(1:end-1)); 
        stdx2 = std(x2(2:end)-x2(1:end-1)); 
        stdy2 = std(y2(2:end)-y2(1:end-1)); 


        % Clip outliers on the low x side: 
        if min(x1)>min(x2)
            y2(x2<min(x1)-stdx2) = []; 
            x2(x2<min(x1)-stdx2) = []; 
        end
        if min(x1)<min(x2)
            y1(x1<min(x2)-stdx1) = []; 
            x1(x1<min(x2)-stdx1) = []; 
        end

        % Clip outliers on the high x side: 
        if max(x1)>max(x2)
            y1(x1>max(x2)+stdx1) = []; 
            x1(x1>max(x2)+stdx1) = []; 
        end
        if max(x2)>max(x1)
            y2(x2>max(x1)+stdx2) = []; 
            x2(x2>max(x1)+stdx2) = []; 
        end

        % Clip outliers on the low y side: 
        if min(y1)>min(y2)
            x2(y2<min(y1)-stdy2) = []; 
            y2(y2<min(y1)-stdy2) = []; 
        end
        if min(y1)<min(y2)
            x1(y1<min(y2)-stdy1) = []; 
            y1(y1<min(y2)-stdy1) = []; 
        end

        % Clip outliers on the high y side: 
        if max(y1)>max(y2)
            x1(y1>max(y2)+stdy2) = []; 
            y1(y1>max(y2)+stdy2) = []; 
        end
        if max(y2)>max(y1)
            x2(y2>max(y1)+stdy1) = []; 
            y2(y2>max(y1)+stdy1) = []; 
        end
        end
    end
end

%% Find intersection x,y point(s): 

P = InterX([x1 y1]',[x2 y2]'); 

% Transform back to lat/lon space: 
[lati,loni] = ps2ll(P(1,:),P(2,:)); 


end 


%% Below is the InterX function written by "NS": 
% http://www.mathworks.com/matlabcentral/fileexchange/22441.m


function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    narginchk(1,2);
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end