classdef mbeMechTypeBase < handle
    % The virtual base class for types of mechanisms: real, simulated.
    
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

