function [] = amt(functionName) 
% amt displays a help window with documenation for Antarctic Mapping Tools functions. 
% 
%% Syntax 
% 
%  amt 
%  amt functionName 
% 
%% Description
% 
% amt typing amt into the command window opens a Help menu with a list of functions available 
% in the Antarctic Mapping Tools package.  
% 
% amt functionName specifies a function name for which you'd like to view documentation.
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at Austin
% Institute for Geophysics (UTIG), February 2017. 
% http://www.chadagreene.com 
% 
% See also: help and showdemo. 

try
   showdemo(strcat(functionName,'_documentation')); 
catch
   if nargin==1
      disp(['Cannot find documentation for an AMT function called ',functionName,'.'])
      disp('Opening the AMT Contents page.')
   end
   showdemo('List_of_Functions') 
end

end
   