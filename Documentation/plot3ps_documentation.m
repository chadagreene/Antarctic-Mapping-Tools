%% |plot3ps| documentation
% |plot3ps| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017). Click <List_of_functions.html here>
% for a complete list of functions in AMT. 
% 
% No Mapping Toolbox?  No problem.  |plot3ps| works just like Matlab's |plot3| function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71°S).
% 
%% Syntax
% 
%  plot3ps(lat,lon,z)
%  plot3ps(...,LineSpec)
%  plot3ps(...,'PropertyName',PropertyValue,...)
%  plot3ps(...,'km')
%  plot3ps(...,'meridian',meridian) 
%  h = plot3ps(...)
% 
%% Description 
% 
% |plot3ps(lat,lon,z)| plots georeferenced data in Antarctic polar stereographic
% eastings and northings. 
% 
% |plot3ps(...,LineSpec)| specifies line or marker style. 
% 
% |plot3ps(...,'PropertyName',PropertyValue,...)| specifies any number of
% line or marker properties. 
% 
% |plot3ps(...,'km')| plots in polar stereographic kilometers instead of the
% default meters. 
% 
% |plot3ps(...,'meridian',meridian)| specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default |meridian| is |0|. 
% 
% |h = plot3ps(...)| returns a column vector of handles to lineseries objects, 
% one handle per line.
% 
%% Example
% Plot gridded bed elevation data as blue dots: 

[centerlat,centerlon] = scarloc('mount vinson'); 
[lat,lon,bed] = bedmap2_data('bed',centerlat,centerlon,50); 

plot3ps(lat,lon,bed,'b.')

%% 
% For more 3D plotting capabilities, check out <surfps_documentation.html
% |surfps|> and see how even you can create neat 3D renderings of Antarctic 
% topography like this one: 
% 
% <<plot3ps.png>>

%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University of Texas
% Institute for Geophysics (UTIG), November 2015, for inclusion in the Antarctic Mapping Tools package.