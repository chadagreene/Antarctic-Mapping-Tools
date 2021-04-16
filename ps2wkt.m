function wkt = ps2wkt(lati_or_xi,loni_or_yi,filename)
% ps2wkt converts polar stereographic (or heck, geographic too) coordinates 
% to a well-known-text polygon. 
% Currently only works for a single bounding box area. 
% 
%% Syntax
% 
%  wkt = ps2wkt(lati,loni) 
%  wkt = ps2wkt(xi,yi)
%  ps2wkt(xi,yi,filename)
% 
%% Description
% 
% wkt = ps2wkt(lati,loni) converts the geographic coordinates lati,loni
% into wkt format, returned as a string. This will typically work best if 
% lati,loni draw a single polygon shape, and it's always best if they do 
% it in a counterclocwise fashion, but who knows, clockwise might work too.
%
% wkt = ps2wkt(xi,yi) the same as above, but for south polar stereographic
% input coordinates xi,yi in meters. 
%
% ps2wkt(xi,yi,filename) does not return a string, but writes the output 
% to a filename such as 'myfile.wkt'. 
% 
%% Example 
% Follow along with this example, one step at a time: 
% 
% 1. Initialize a map of Antarctica: 
% 
%  figure
%  bedmachine
% 
% 2. Now use coord to click around on the map to get coordinates: 
% 
%  [xi,yi] = coord
%
% 3. Let's say you clicked on these xi,yi coordinates around the 
% Pine Island Glacier. Convert them to a .wkt polygon: 
% 
%  xi = [-1634539    -1803223    -1620482    -1325287]; 
%  yi = [-43829     -212512     -577993     -381196]; 
%
%  ps2wkt(xi,yi)
%  ans =
%     'POLYGON((-91.5359774956602 -96.7213734589033 -109.630284519675 -106.046970607966 -75.0331055800544 -73.4012772678685 -74.26104695933 -77.3574435910396))'
% 
% 4. Alternatively, convert write that polygon to a .wkt text file: 
% 
%  ps2wkt(xi,yi,'myfile.wkt')
% 
%% Author Info
% Chad Greene, September 2020. 

%% Parse inputs: 
 
narginchk(2,3) 
assert(isequal(size(lati_or_xi),size(loni_or_yi)),'Error: Dimensions of input coordinates do not agree.')

if ~islatlon(lati_or_xi,loni_or_yi)
   [lati,loni] = ps2ll(lati_or_xi,loni_or_yi); 
else 
   lati = lati_or_xi;
   loni = loni_or_yi;    
end

%% Convert coordinates to text: 

coords = mat2str([loni lati]);
coords = strrep(coords, '[', '(');
coords = strrep(coords, ']', ')');
coords = strrep(coords, ';', ', ');
coords = regexprep(coords, '(, NaN NaN, )+', '), (');

wkt = ['POLYGON(',coords,')']; 

%% Write to .wkt file 

if nargin==3
   fid = fopen(filename,'w'); 
   fprintf(fid,wkt);
   fclose(fid); 
   
   % If the user wants to write a file and didn't specifically ask for a variable output, clear wkt.  
   if nargout==0
      clear wkt
   end
end


end