classdef mbeEstimatorIncrManifold_errorEKF_shaping < mbeEstimatorFilterBase
    % This function runs a filter which tries to estimate the errors between
    % the real system and the model, so the model is a conventional MB
    % simulation, wich can have its own integrator, formulation, etc, while
    % the model for the errors remains as simple as possible.
    %
    % This filter estimates the states and forces and it is combined with a
    % shaping filter for estimating the plant noise. This is useful for 
    % situations where the plant noise is not white gaussian noise and it
    % is colored.
    %
    % Set the "dynamics_formulation_and_integrator" property to the desired
    % dynamic formulation & integrator. 
    % This property takes preference over the base class property "post_iter_accel_solver", 
    % which is NOT used in this class.
    %
    
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
        % Method to simulate the MBS: Set to any of the dynamics
        % formulation, and configure to the desired integrator.
        %dynamics_formulation_and_integrator = mbeDynFormulationMatrixR();
        dynamics_formulation_and_integrator = mbeDynFormulationI3AL();
        %dynamics_formulation_and_integrator = mbeDynFormulationLagrangeBaumgarte();        
        
        w_noise;
        W1;
        W2;
        w_noise_hist;
        inn_window;
        autocorr_hist;
        W1_hist;
        xplus;
    end
    
    % Private vars
    properties(Access=private)
        CovPlantNoise; 
        CovMeasurementNoise;
        F0; % Transition matrix
        X0; % State from previous step (also used for initialization)
        deltaqpp;
        Iz;
        Oz;
        W; 
        
        %% Variables for the augmented matrix for the noises
        P_aug;
        F_aug;
        F0_aug;
        H_aug;
        w_sigma;
        CovPlant_aug;
        CovMeas_aug;
    end
    
    methods(Access=protected)
        
        % Returns the location (indices) in P for the independent
        % coordinates (z,zp,zpp). Zeros means not estimated by this filter.
        function [P_idx_z, P_idx_zp, P_idx_zpp] = get_indices_in_P(me)
            % Structure state vector (and its cov P) in this filter:
            P_idx_z   =  1:me.lenZ;
            P_idx_zp  = (1:me.lenZ) + me.lenZ;
            P_idx_zpp = (1:me.lenZ) + 2*me.lenZ; 
        end
                
        
        % Init the filter (see docs in mbeEstimatorFilterBase)
        function [] = init_filter(me)
            % 1) q,qp,qpp: already set in base class.
            % 2) Initial covariance: 
            me.P = diag([...
                me.initVar_Z*ones(1,me.lenZ), ...
                me.initVar_Zp*ones(1,me.lenZ),...
                me.initVar_Zpp*ones(1,me.lenZ)]);
            
            me.Iz = eye(me.lenZ);
            me.Oz = zeros(me.lenZ);
            
            me.w_noise = zeros(3*me.lenZ,1);
                                  
            me.P_aug = [ me.P, 0*eye(3*me.lenZ);
                        zeros(3*me.lenZ), 0.001*eye(3*me.lenZ)];
            
            % Terms of the covariance plant noise matrix related to the
            % multibody model
            me.CovPlantNoise = diag([...
                zeros(1,me.lenZ), ...
                zeros(1,me.lenZ),... 
                ones(1,me.lenZ)*me.transitionNoise_Zpp]);
            
            % Terms of the covariance plant noise matrix related to the
            % noise estimation
            me.w_sigma = diag([...
                zeros(1,me.lenZ), ...
                zeros(1,me.lenZ),... 
                ones(1,me.lenZ)]);
            
            % Noise for the 4bar linkage
            if size(me.mech_phys_model.indep_idxs, 2) == 1
                me.w_sigma(3,3) = 2.3428e-06;
            elseif size(me.mech_phys_model.indep_idxs, 2) == 2
                % Noise for the 5bar linkage
                me.w_sigma(5,5) = 7.1105e-04;
                me.w_sigma(6,6) = 2.1196e-04;
            end
            
            % Assembly of the augmented plant covariance matrix
            me.CovPlant_aug = 1 *[ 0*me.CovPlantNoise, zeros(3*me.lenZ);
                               zeros(3*me.lenZ), me.w_sigma];
           
            me.F0 = [ me.Iz, me.Iz*me.dt, 0.5*me.Iz*me.dt^2;
                     me.Oz, me.Iz, me.Iz*me.dt
                     me.Oz, me.Oz, me.Iz];
                 
            % W: parameter which can take any value from 0 to 1 to modulate the
            % influence of the previous value of the noise. It is initially
            % set to 0 and estimated through the algorithm
            me.W1(1:me.lenZ) = 0*0.997; % 0.997 is the expected value
            me.W2 = 0;
            
            weight_matrix = eye(3*me.lenZ);
            for i = 1:me.lenZ
                weight_matrix(3*i-2:3*i,3*i-2:3*i) = me.W1(i)*weight_matrix(3*i-2:3*i,3*i-2:3*i);
            end
            
            me.F0_aug = [ me.F0, eye(3*me.lenZ);
                         zeros(3*me.lenZ), weight_matrix];     
                 
            me.W = 1; 
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter .* me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2); % The noise is already modified with the shaping
            me.X0 = zeros(3*me.lenZ, 1);
            me.deltaqpp = 0;
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            
            % 1) Run dynamics: 
            [qMB,qpMB,qppMB] = me.dynamics_formulation_and_integrator.integrate_one_timestep(...
                me.bad_mech_phys_model,...
                me.current_run_time, ...
                me.dt, ...
                me.q, me.qp, me.qpp);

            % 2) Kalman state propagation
            % the state is DOF position and velocity error            
            [dzpp_dz,dzpp_dzp] = me.bad_mech_phys_model.eval_dzpp_dposvel(me.q,me.qp);

            F = me.F0 + ...
                [0.5*dzpp_dz*me.dt^2,   0.5*dzpp_dzp*me.dt^2,   me.Oz;
                dzpp_dz*me.dt,          dzpp_dzp*me.dt,         me.Oz;
                me.Oz,                       me.Oz,             me.Oz];
            
            weight_matrix = eye(3*me.lenZ);
            for i = 1:me.lenZ
                weight_matrix(3*i-2:3*i,3*i-2:3*i) = me.W1(i)*weight_matrix(3*i-2:3*i,3*i-2:3*i);
            end
            
            me.F_aug = [ F, eye(3*me.lenZ);
                        zeros(3*me.lenZ), weight_matrix];
            
            X_minus = me.F_aug*[me.X0; me.w_noise];
            % time update
            P_minus = me.F_aug*me.P_aug*me.F_aug' + me.CovPlant_aug; %*me.W;
                        
            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp, dh_dzpp] = me.bad_mech_phys_model.sensors_jacob_indep(qMB,qpMB,qppMB); % me.q,me.qp,me.qpp);
                H=[dh_dz, dh_dzp, dh_dzpp]; 
                
                me.H_aug = [H, zeros(size(me.CovMeasurementNoise,1),3*me.lenZ)];
                
                % Kalman gain:
                K = P_minus*me.H_aug'/(me.H_aug*P_minus*me.H_aug'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(qMB,qpMB,qppMB);
                
                me.Innovation = obs-obs_predict;
                X_plus = X_minus + K*me.Innovation;
                me.P_aug = (eye(length(X_minus))-K*me.H_aug)*P_minus;
            else
                % No sensor:
                X_plus = X_minus;
                me.P_aug = P_minus;
                me.Innovation= [];
            end
            
            me.P = me.P_aug(1:3*me.lenZ,1:3*me.lenZ);

            % Recover KF -> Plant and measurement noise
            me.w_noise = X_plus(3*me.lenZ+1:2*3*me.lenZ);
            
            % Data logging
            me.w_noise_hist(:,me.current_timestep) = me.w_noise;
            
            % For studying sequence
            me.xplus = X_plus(1:3*me.lenZ);
            
            %% Estimation of W1
            win_size = 250;
            % Update sliding window
            if me.current_timestep <= 2
                me.inn_window(:,me.current_timestep) = me.xplus(2*me.lenZ+1:end);
                me.autocorr_hist(:,me.current_timestep) = zeros(me.lenZ,1);
            elseif me.current_timestep <= win_size && win_size > 0
                me.inn_window(:,me.current_timestep) = me.xplus(2*me.lenZ+1:end);
                % Estimate W1 as the autocorrelation of the innovation set
                lag = 1:win_size; 
                x = me.inn_window';
                lcorr = 2; % numero de desfases para el que se calcula la autocorrelación (el que interesa es el 2, que es un desfase de una posición)
                autocorrelation = zeros(lcorr,me.lenZ); 
                for j = 1:me.lenZ
                    for i = 1:lcorr
                        x_lag = x(i:end,j);
                        x_ini = x(1:(end+1-i),j);
                        autocorrelation(i,j) = (x_lag-mean(x(:,j)))'*(x_ini-mean(x(:,j)))/sum(((x(:,j)-mean(x(:,j))).^2));
                    end
                    me.W1(j) = me.W1(j) + autocorrelation(2,j)/me.current_timestep;
                    if me.W1(j) < 0; me.W1(j) = 0;
                    elseif me.W1(j) > 1; me.W1(j) = 1;
                    end
                end
                me.autocorr_hist(:,me.current_timestep) = autocorrelation(2,:);
            elseif win_size > 0
                me.inn_window(:,1:win_size-1) = me.inn_window(:,2:win_size);
                me.inn_window(:,win_size) = me.xplus(2*me.lenZ+1:end);
                % Estimate W1 as the autocorrelation of the innovation set
                lag = 1:win_size; 
                x = me.inn_window';
                lcorr = 2; % numero de desfases para el que se calcula la autocorrelación (el que interesa es el 2, que es un desfase de una posición)
                autocorrelation = zeros(lcorr,me.lenZ); 
                for j = 1:me.lenZ
                    for i = 1:lcorr
                        x_lag = x(i:end,j);
                        x_ini = x(1:(end+1-i),j);
                        autocorrelation(i,j) = (x_lag-mean(x(:,j)))'*(x_ini-mean(x(:,j)))/sum(((x(:,j)-mean(x(:,j))).^2));
                    end
                    me.W1(j) = me.W1(j) + autocorrelation(2,j)/win_size;
                    if me.W1(j) < 0; me.W1(j) = 0;
                    elseif me.W1(j) > 1; me.W1(j) = 1;
                    end
                end
                me.autocorr_hist(:,me.current_timestep) = autocorrelation(2,:);
            end
            
            
            % Store the W1 for postprocessing
            me.W1_hist(:,me.current_timestep) = me.W1;
            

            %% Recover KF -> MBS coordinates
            % ------------------------------
            % Position increments fulfill velocity constraints
            deltaq = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, qMB, X_plus(1:me.lenZ) ); 
            me.q = qMB + deltaq;
            
            % Corrected velocities are calculated by applying the velocity
            % problem from the corrected velocities of the degrees of freedom 
            me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, qpMB(me.bad_mech_phys_model.indep_idxs)+X_plus( me.lenZ+(1:me.lenZ) ));

            zpp = qppMB(me.bad_mech_phys_model.indep_idxs)+X_plus(2*me.lenZ+(1:me.lenZ));
            
            R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,me.q);
            M_bar = R'*me.bad_mech_phys_model.M*R;

            deltaQ = M_bar* X_plus(2*me.lenZ+(1:me.lenZ));
            me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) = me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) + deltaQ;
            
            % Acceleration correction
            me.qpp = mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model, me.q, me.qp, zpp);
            me.X0 = zeros(3*me.lenZ,1);
       
        end % run_filter_iter
        
    end
end

