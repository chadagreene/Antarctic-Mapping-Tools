%% |contourfps| documentation
% |contourfps| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017).  <List_of_functions.html Click here
% for a complete list of functions in AMT.>
% 
% The |contourfps| works just like Matlab's |contourf| or |contourfm| functions, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71°S).
% For example, 
% 
%  contourfps(lat,lon,Z) 
% 
% is equivalent to 
% 
%  [x,y] = ll2ps(lat,lon); 
%  contour(x,y,Z) 
% 
%% Syntax
% 
%  contourfps(lat,lon,Z)
%  contourfps(lat,lon,Z,n)
%  contourfps(lat,lon,Z,v)
%  contourfps(...,LineSpec)
%  contourfps(...,'km')
%  contourfps(...,'meridian',meridian) 
%  [C,h] = contourfps(...)
% 
%% Description 
% 
% |contourfps(lat,lon,Z)| draws contours of |Z| at gridded locations |lat, lon|. 
% 
% |contourfps(lat,lon,Z,n)| specifies a number of contour levels |n| if |n| is a scalar. 
% 
% |contourfps(lat,lon,Z,v)| draws a contour plot of matrix |Z| with contour lines at the 
% data values specified in the monotonically increasing vector |v|. The number of contour 
% levels is equal to |length(v)|. To draw a single contour of level |i|, use |v = [i i]|; 
% 
% |contourfps(...,LineSpec)| draws the contours using the line type and color specified 
% by |LineSpec|. |contourfps| ignores marker symbols.
% 
% |contourfps(...,'km')| plots in polar stereographic kilometers instead of the
% default meters. 
% 
% |contourfps(...,'meridian',meridian)| specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default |meridian| is |0|. 
% 
% |[C,h] = contourfps(...)| returns contour matrix |C| and handle |h| of the contour object created.
% 
%% Example 1: Simple case 
% Contour some surface elevation data from <http://www.mathworks.com/matlabcentral/fileexchange/42353-bedmap2-toolbox-for-matlab/content/bedmap2_toolbox_v4.0/Bedmap2_documentation/html/bedmap2_data_documentation.html
% |bedmap2|>:

[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); 

contourfps(lat,lon,sfz) 
 
%% Example 2: Polar stereographic kilometers
% Using data from Example 1, plot in polar stereographic kilometers.

figure
[C,h] = contourfps(lat,lon,sfz,'km');

xlabel 'eastings (km)' 
ylabel 'northings (km)' 

clabel(C,h,'LabelSpacing',300,'fontsize',8)
colormap(jet)

%% Example 3: Specify levels 
% Specify every 200 meters surface elevation from sea level to 5 km:

figure
contourfps(lat,lon,sfz,0:200:5000,'km') 
cb = colorbar; 
ylabel(cb,'surface elevation (m)')
 
%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function was written by <http://www.chadagreene.com Chad Greene> of the University of Texas
% Institute for Geophysics (UTIG), May 2017, for inclusion in the Antarctic Mapping Tools package. 
