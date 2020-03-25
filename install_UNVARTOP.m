%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% A 200 line topology optimization code using the Unsmooth Variational Topology method %
%   for structural problems                                                            %
%   Paper's Title                                                                      %
%       Educational explanation of Topology Optimization using the UNsmooth VARiational%
%       Topology OPtimization (UnVarTop) method: an efficient implementation in Matlab %
%   Authors                                                                            %
%       Daniel Yago, Juan Carlos Cante, Oriol Lloberas-Valls, Javier Oliver            %
% Versions                                                                             %
%   26/3/2020 - Submission version (install package)                                   %
% GitHub repository																	   %
%	This code and other modifications of it can be downloaded from the Website:		   %
%		https://github.com/DanielYago/UNVARTOP	. Please send your comments to		   %
%		the author: daniel.yago@upc.edu								   				   %
% Copyright																			   %
%	Copyright (C) 2020  Daniel Yago													   %
%																					   %
%	This program is free software: you can redistribute it and/or modify it under the  %
%	terms of the GNU General Public License as published by the Free Software 		   %
%	Foundation, either version 3 of the License, or (at your option) any later version.%
%																					   %
%	This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%	PARTICULAR PURPOSE.  See the GNU General Public License for more details.          %
%																					   %
%	You should have received a copy of the GNU General Public License along with this  %
%	program.  If not, see <https://www.gnu.org/licenses/>.							   %
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%

% Define the path where the codes are stored
dir_path = fullfile(pwd,'Source_codes');

% Add the previous path to Matlab's path
p = genpath(dir_path);
addpath(p);