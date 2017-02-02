classdef mbeMechTypeSimulated < mbeMechTypeBase
    % A simulated mechanism. See docs for mbeMechTypeBase
    
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
    % --- Public properties  ---
    properties (Access=public)
        % Pick one dyn formulation for the simulations
        dynamic_formulation = mbeDynFormulationMatrixR();
        
        % Enable/disable showing graphs after simulation
        show_graphs_after_GT_sim = 1;
        show_graphs_after_BadModel_sim = 1;
        
        % Enable simulating "multirate sensors": set this to the period
        % between sensor data (e.g. 1.0/f for f in Hz). Leave to 0
        % (default) to simulate sensor data at every time step:
        multirate_sensor_period = 0.0;
        
    end
    
    % --- Public API ---
    methods
        % see docs in mbeMechTypeBase
        function [GT, BadModel] = generate_ground_truth(self,estim)
            self.dynamic_formulation.logging = estim.logging; % Reuse those options
            estim.init_bad_model(); % Make sure the bad model is up-to-date
            
            % 1) Do we need to re-generate datasets?
            simul_params = struct();
            simul_params.mod=estim.mech_phys_model;
            simul_params.mod.installed_sensors = []; % Remove sensors from the list of fields to compare (changes there do not affect to GT)
            simul_params.end_time = estim.end_time;
            simul_params.dt = estim.dt;    
            simul_params.bad_model_errors = estim.bad_model_errors;
            
            do_regenerate=0;
            if (   ~exist('dataset_cache','dir') ...
                || ~exist('dataset_cache/simul_params.mat','file') ...
                || ~exist('dataset_cache/dataset_GT.mat','file') ...
                || ~exist('dataset_cache/dataset_BadModel.mat','file') ...
                )
                do_regenerate=1;
                if (~exist('dataset_cache','dir'))
                    mkdir('dataset_cache');
                end
            else
                tst_simul_params = load('dataset_cache/simul_params.mat');
                if (~ isequaln(simul_params,tst_simul_params.simul_params))
                    estim.log_level(1,'[mbeMechTypeSimulated] ** DETECTED A CHANGE IN THE MODEL/PARAMS: WILL REGENERATE DATASETS**');
                    do_regenerate=1;
                else
                    estim.log_level(1,'[mbeMechTypeSimulated] Cached dataset matches all params. Reusing it.');
                end
            end

            % 2) Re-generate (if needed)
            if (do_regenerate)
                % We are in a simulated model, so: 
                %  - GT: Is the dyn simulation of the actual model.
                %  - BadModel: idem, for a corrupted model.
                
                % 1) Call solver with "perfect" model:
                % -------------------------------------
                estim.log_level(1,'[mbeMechTypeSimulated]  ---- Simulating GT model ---');
                [GT,GT_perf_stats]=self.dynamic_formulation.batch_dyn_solve(...
                    estim.mech_phys_model, ...  % GT model
                    0.0,estim.end_time, estim.dt);  % old was: MBS_solve()
                
                if (self.show_graphs_after_GT_sim)
                    self.dynamic_formulation.plot_perf_stats(GT_perf_stats); % OLD: MBS_plots(...);
                end
                
                % 2) "GenerateDatasets_Imperfect()":
                % Apply a distortion to the model, relaunch the dynamic
                % simulation
                % -------------------------------------
                estim.log_level(1,'[mbeMechTypeSimulated]  ---- Simulating "bad" model ---');
                [BadModel,BadModel_perf_stats]=self.dynamic_formulation.batch_dyn_solve(...
                    estim.bad_mech_phys_model, ... % "Bad" model
                    0.0,estim.end_time, estim.dt);  % old was: MBS_solve()
                
                if (self.show_graphs_after_BadModel_sim)
                    self.dynamic_formulation.plot_perf_stats(BadModel_perf_stats); % OLD: MBS_plots(...);
                end                

                % 3) and save cache:
                % ---------------------
                save('dataset_cache/simul_params.mat','simul_params');
                save('dataset_cache/dataset_GT.mat','GT');                                     
                save('dataset_cache/dataset_BadModel.mat','BadModel');                
            else
                % 3) or, load from files:
                GT=load('dataset_cache/dataset_GT.mat'); GT=GT.GT;
                BadModel=load('dataset_cache/dataset_BadModel.mat');BadModel=BadModel.BadModel;
                
                % Sanity check:
                assert(size(GT.q,1)==ceil(simul_params.end_time/simul_params.dt),'Inconsistent GT data length! Please, delete directory "dataset_cache" and re-run');
            end
            
        end % of generate_ground_truth()

        % see docs in mbeMechTypeBase
        function [obs] = get_current_observations(self, estim, gt_q,gt_qp,gt_qpp)
            if (self.multirate_sensor_period>0)
                time_epsilon = self.multirate_sensor_period * 1e-6; % To avoid missed sensor readings due to numeric accuracy issues
                do_gen_obs = (estim.current_run_time + time_epsilon >= self.sensor_last_simul_time+self.multirate_sensor_period);
            else
                do_gen_obs = 1;
            end
            
            if (do_gen_obs)
                % save time:
                self.sensor_last_simul_time = estim.current_run_time;
                
                % Use the GT model & the GT dyn state to simulate sensors:
                obs = estim.mech_phys_model.sensors_simulate(gt_q,gt_qp,gt_qpp);

                % Add noise:
                obs_stds = estim.mech_phys_model.sensors_std_noise();
                obs = obs + obs_stds.*randn(length(obs),1);           
            else
                obs=[]; % no observations at this time
                % The random number generation is invoqued anyway to keep the same sequence of pseudorandom numbers
                randn(length(estim.mech_phys_model.installed_sensors),1);
            end
        end % of get_current_observations()
        
        % Do whatever is needed to clear the state before starting a fresh
        % simulation from t=0
        function [] = resetState(me)
           me.sensor_last_simul_time = 0.0; 
        end
        
    end
    
    properties (Access=private)
        sensor_last_simul_time = 0.0; % Used for multirate in get_current_observations()
    end
    

end

