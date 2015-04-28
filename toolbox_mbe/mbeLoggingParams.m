classdef mbeLoggingParams
    % Options to control the level of verbosity of debugging messages,
    % graphics, etc.
    
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

    properties(Access=public)
        % Frequency (FPS) of display graphics of the mechanism animation
        % (0=hide graphs)
        show_mech_anim_fps = 10;
        
        % Verbose level: 0=none; 1:info; 2:verbose
        verbose_level = 1;  
    end

    methods(Access=public)
        % Default constructor:
        function obj = mbeLoggingParams()
            obj.tic_global_id = tic;
        end
        % Display a 'msg' if user's verbose_level>=VERBOSE_LEVEL
        function []=log_level(self, VERBOSE_LEVEL, msg)
            if (VERBOSE_LEVEL<=self.verbose_level)
                fprintf('[%.04f] %s\n',toc(self.tic_global_id), msg);
            end
        end        
    end    
    
    % Private stuff:
    properties(Access=protected)
        tic_global_id; % for tic/toc in log msgs
    end    
    
end

