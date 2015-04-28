classdef mbeMechTypeBase < handle
    % The virtual base class for types of mechanisms: real, simulated.
    
    % --- Public API ---
    methods(Abstract)
        % Generate a Ground Truth solution of the system dynamics, 
        % or loaded the pre-cached one. This method detects if the model 
        % has changed since the last cached file and re-generates only if
        % required.
        % - estim [in]: The estimator object. All required params are
        % gathered from its public properties.
        % - GT [out]: A struct with Nx1 (N=nTimeSteps) q,qp,qpp
        [GT, BadModel] = generate_ground_truth(self,estim );

        % Get the latest observation from the system: either real-time real
        % sensor, or offline simulated. Current time is accessed via
        % properties of mbeEstimatorFilterBase. 
        % - estim [in]: The estimator object
        % - gt_states [in]: (If applicable) GT states from a dyn simulation
        % - obs [out]: The vector of observations (variable length, may be empty)
        [obs] = get_current_observations(self, estim, gt_q,gt_qp,gt_qpp);
        
        % Do whatever is needed to clear the state before starting a fresh
        % simulation from t=0
        [] = resetState(self);
       
        
    end

end

