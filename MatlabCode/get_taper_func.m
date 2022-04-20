function varargout = get_taper_func(varargin);


%%
%% 'get_taper_func'
%%
%% list available Matlab taper functions
%% that take only taper length as an input 
%% argument;
%% return the function handle for a given 
%% taper function name
%% 
%% syntax:
%% -------
%%   taper_name_list = get_taper_func;
%%   taper_function_handle = get_taper_func(taper_function_name);
%%
%% input:
%% ------
%%   taper_function_name == taper function name string
%%
%% output:
%% -------
%%   taper_function_handle == function handle corresponding to 
%%                            taper name
%%
%% K. A. Cortopassi, April 2004
%%


name_list{1} = 'Bartlett-Hann';
name_list{2} = 'Bartlett';
name_list{3} = 'Blackman' ;
name_list{4} = 'Blackman-Harris';
name_list{5} = 'Bohman';
name_list{6} = 'Flat Top';
name_list{7} = 'Gaussian';
name_list{8} = 'Hamming';
name_list{9} = 'Hann';
name_list{10} = 'Nuttall Blackman-Harris';
name_list{11} = 'Parzen de la Valle-Poussin';
name_list{12} = 'Rectangular';
name_list{13} = 'Triangular';

function_list{1} = @barthannwin;
function_list{2} = @bartlett;
function_list{3} = @blackman; 
function_list{4} = @blackmanharris;
function_list{5} = @bohmanwin;
function_list{6} = @flattopwin;
function_list{7} = @gausswin;
function_list{8} = @hamming;
function_list{9} = @hann;
function_list{10} = @nuttallwin;
function_list{11} = @parzenwin;
function_list{12} = @rectwin;
function_list{13} = @triang;

if ~length(varargin)
  varargout{1} = name_list;
  
else
  window_name = varargin{1};
  varargout{1} = function_list{find(strcmpi(name_list, window_name))};
  
end

return;
