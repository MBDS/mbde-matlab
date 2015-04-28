classdef mbeMechTypeReal < mbeMechTypeBase
    % A real mechanism. See docs for mbeMechTypeBase
    
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
        
        OBS; % Observations
        
        dataset; % File to load
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
                % We are in a real model, so:
                %  - GT: Is the kinematic solution for the encoder 1
                %  - BadModel: Is the dyn simulation of the actual model.
                
                % 1) Calculate q,qp and qpp from the encoder in the crank:
                % -------------------------------------
                estim.log_level(1,'[mbeMechTypeSimulated]  ---- Simulating GT model ---');
                addpath('toolbox_mbe','datasets_originales');
                matriz = cargar_outside(self.dataset,estim.end_time);
                q_aprox = estim.mech_phys_model.q_init_aprox;
                % S e n s o r s
                % Just for checking if the IMU axis reference are correct
                migiroAcoplador = mbeSensorGyroscope(...
                    estim.mech_phys_model.installed_sensors{1}.pt_is_fixed,...
                    estim.mech_phys_model.installed_sensors{1}.pt1_idxs,...
                    estim.mech_phys_model.installed_sensors{1}.pt2_idxs, 0);
                migiroBalancin = mbeSensorGyroscope(...
                    estim.mech_phys_model.installed_sensors{2}.pt_is_fixed,...
                    estim.mech_phys_model.installed_sensors{2}.pt1_idxs,...
                    estim.mech_phys_model.installed_sensors{2}.pt2_idxs, 0);
                % Encoder constants
                kBALAN = -(1/(10000*4))*2*pi; % pulsos_contados/(360ppr*x4)*rad/rev ---> rad
                kMANIVELA = -(1/(128*4))*2*pi; % pulsos_contados/(128ppr*x4)*rad/rev ---> rad
                kMANIVELAw = kMANIVELA/8;
                t = matriz.t_enc;
                % se define la variable largo para unificar long. de sensores
                gyro_vel_largoACOPL = -matriz.gyro_velACOPL;
                gyro_velACOPL = gyro_vel_largoACOPL(1:length(t));
                gyro_vel_largoBALAN = -matriz.gyro_velBALAN;
                gyro_velBALAN = gyro_vel_largoBALAN(1:length(t));
                enc0 = matriz.encMANIVELA*kMANIVELA; enc0 = deg2rad(6*(enc0-enc0(1))); % numero magico Maxon: 6?¿?¿
                enc1 = matriz.encBALAN*kBALAN;                
                enc0vel = matriz.encMANIVELAvel*kMANIVELAw;
                enc0acc = zeros(length(t),1); % TODO: CALCULATE REAL ACCELERATIONS FROM ENCODER (TOO NOISY)
                for i=2:length(t)
                    enc0acc(i) = (enc0vel(i)-enc0vel(i-1))/0.005;
                end
                pasos = length(enc0);
                frame_rate = 1/10; % Every frame_rate s of simulation time a frame is plotted
                frame = frame_rate;
                animacion = 1;
                if animacion == 1
                    figure('Name','Simulacion cinematica')
                end
                GT.t = zeros(pasos,1);
                GT.q = zeros(pasos,5);
                GT.qp = zeros(pasos,5);
                GT.qpp = zeros(pasos,5);
                OBS = zeros(pasos,2);
                OBS(:,1) = gyro_velACOPL(1:pasos);
                OBS(:,2) = gyro_velBALAN(1:pasos);
                % auxiliar temporal variables
                qt = mbeKinematicsSolver.pos_problem(estim.mech_phys_model,q_aprox);
                qtp = mbeKinematicsSolver.vel_problem(estim.mech_phys_model,qt,enc0vel(1));
                qtpp = mbeKinematicsSolver.accel_problem(estim.mech_phys_model,qt,qtp,enc0acc(1));
                wAcoplador = zeros(pasos,1);
                wBalancin  = zeros(pasos,1);
                
                for i=1:pasos
                    qt=mbeKinematicsSolver.pos_problem(estim.mech_phys_model,qt);
                    qtp = mbeKinematicsSolver.vel_problem(estim.mech_phys_model,qt,qtp(5));
                    qt(5) = enc0(i);
                    qtp(5) = enc0vel(i);
                    qtpp(5) = enc0acc(i);
                    qtpp = mbeKinematicsSolver.accel_problem(estim.mech_phys_model,qt,qtp,qtpp(5));
                    if t(i) > frame
                        if animacion == 1
                            lineas_outside(qt)
                        end
                        frame = frame + frame_rate;
                    end
                    wAcoplador(i) = migiroAcoplador.sensor_simulate(estim.mech_phys_model,qt,qtp);
                    wBalancin(i) = migiroBalancin.sensor_simulate(estim.mech_phys_model,qt,qtp);
                    GT.t(i) = t(i);
                    GT.q(i,:) = qt;
                    GT.qp(i,:) = qtp;
                    GT.qpp(i,:) = qtpp;
                end
                
                figure('name','comp')
                subplot(2,1,1)
                plot(GT.t,OBS(:,1),GT.t,wAcoplador(1:4000))
                legend('imu','model')
                title('w acoplador')
                subplot(2,1,2)
                plot(GT.t,OBS(:,2),GT.t,wBalancin(1:4000))
                legend('imu','model')
                title('w blancin')
                
                % 2) "GenerateDatasets_Imperfect()":
                % Apply a distortion to the model, relaunch the dynamic
                % simulation
                % -------------------------------------
                estim.log_level(1,'[mbeMechTypeSimulated]  ---- Simulating "bad" model ---');
                [BadModel,BadModel_perf_stats]=self.dynamic_formulation.batch_dyn_solve(...
                    estim.mech_phys_model, ... % "Bad" model
                    0.0,estim.end_time, estim.dt);  % old was: MBS_solve()
                
                if (self.show_graphs_after_BadModel_sim)
                    self.dynamic_formulation.plot_perf_stats(BadModel_perf_stats); % OLD: MBS_plots(...);
                end
                
                % 3) and save cache:
                % ---------------------
                save('dataset_cache/simul_params.mat','simul_params');
                save('dataset_cache/dataset_GT.mat','GT');
                save('dataset_cache/dataset_OBS.mat','OBS');
                self.OBS = OBS;
                save('dataset_cache/dataset_BadModel.mat','BadModel');
            else
                % 3) or, load from files:
                GT=load('dataset_cache/dataset_GT.mat'); GT=GT.GT;
                BadModel=load('dataset_cache/dataset_BadModel.mat');BadModel=BadModel.BadModel;
                load('dataset_cache/dataset_OBS.mat');
                self.OBS = OBS;                
                % Sanity check:
                assert(size(GT.q,1)==ceil(simul_params.end_time/simul_params.dt),'Inconsistent GT data length! Please, delete directory "dataset_cache" and re-run');
            end
            
        end % of generate_ground_truth()
        
        % see docs in mbeMechTypeBase
        function [obs] = get_current_observations(self, estim, gt_q,gt_qp,gt_qpp)
            %load('./dataset_cache/dataset_OBS.mat')
            obs = self.OBS(estim.current_timestep,:)';
        end % of get_current_observations()
        
        % Do whatever is needed to clear the state before starting a fresh
        % simulation from t=0
        function [] = resetState(me)

        end
        
    end
    
    properties (Access=private)

    end
    
end

