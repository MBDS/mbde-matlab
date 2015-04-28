classdef mbeLoggingParams
    % Options to control the level of verbosity of debugging messages,
    % graphics, etc.
    
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

