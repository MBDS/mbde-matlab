classdef mbeModelErrorDef
    % Definition of the kind of errors to introduce in the MB model
    % The effects are implemented in each model's applyErrors() virtual
    % method.
    
	% -----------------------------------------------------------------------------
	% This file is part of MBDE-MATLAB.  See: https://github.com/MBDS/mbde-matlab
	% 
	%     MBDE-MATLAB is free software: you can redistribute it and/or modify
	%     it under the terms of the GNU General Public License as published by
	%     the Free Software Foundation, either version 3 of the License, or
	%     (at your option) any later version.
	% 
	%     MBDE-MATLAB is distributed in the hope that it will be useful,
	%     but WITHOUT ANY WARRANTY; without even the implied warranty of
	%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	%     GNU General Public License for more details.
	% 
	%     You should have received a copy of the GNU General Public License
	%     along with MBDE-MATLAB.  If not, see <http://www.gnu.org/licenses/>.
	% -----------------------------------------------------------------------------
    properties
        % Set error_scale to larger or smaller values to have more or less
        % error. "error_type" selects type of error in the "Wrong model":
        % - 0: No intentionally wrong (for real experiments)
        % - 1: Gravity
        % - 2: Init. pos
        % - 3: Init. vel
        % - 4: Damping "C" param = 0
        % - 5: Damping "C" param += 10
        error_type = 1;
        
        % Magnitude of the error
        error_scale = 1;
        
    end
    
end

