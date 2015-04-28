classdef mbeModelErrorDef
    % Definition of the kind of errors to introduce in the MB model
    % The effects are implemented in each model's applyErrors() virtual
    % method.
    
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

