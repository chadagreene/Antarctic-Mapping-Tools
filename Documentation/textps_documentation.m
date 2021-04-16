%% |textps| documentation
% |textps| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017).  <List_of_functions.html Click here
% for a complete list of functions in AMT.>
% 
% The |textps| function works just like Matlab's |text| or |textm| functions, but places georeferenced
% text labels in Antarctic polar stereographic cartesian coordinates (true
% latitude 71°S). For example, 
% 
%  textps(lat,lon,'My text') 
% 
% is equivalent to: 
% 
%  [x,y] = ll2ps(lat,lon); 
%  plot(x,y,'My text') 
% 
%% Syntax
% 
%  textps(lat,lon,'string')
%  textps(...,'PropertyName',PropertyValue,...)
%  textps(...,'km')
%  textps(...,'meridian',meridian) 
%  h = textps(...)
% 
%% Description 
% 
% |textps(lat,lon,'string')| adds the string in quotes to the location specified 
% by the point (lat,lon) polar stereographic eastings and northings. 
% 
% |textps(...,'PropertyName',PropertyValue,...)| formats the string with any
% name-value pairs of text properties. 
% 
% |textps(...,'km')| plots in polar stereographic kilometers instead of the
% default meters. 
% 
% |textps(...,'meridian',meridian)| specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default |meridian| is |0|. 
% 
% |h = textps(...)| returns a column vector of handles to text objects, one handle 
% per object.
% 
%% Example 
% Place text labels in East and West Antarctica: 

antbounds('gl') 
textps(-80,-105,'West Antarctica') 
textps(-75,90,'East Antarctica',...
  'fontangle','italic','color','red')
xlabel 'eastings (m)' 
ylabel 'northings (m)' 

%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University of Texas
% Institute for Geophysics (UTIG), July 2015, for inclusion in the
% Antarctic Mapping Tools package.