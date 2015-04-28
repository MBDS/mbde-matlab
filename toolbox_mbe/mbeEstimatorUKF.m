classdef mbeEstimatorUKF < mbeEstimatorFilterBase
    % Unscented Kalman filter (UKF) 
    % Configurable numeric integrator & dynamics formulation (see property
    % "dynamics_formulation_and_integrator"). 
    % This property takes preference over the base class property "post_iter_accel_solver", 
    % which is NOT used in this class.
    %
    % As of Roland's paper (+ the optional "plan noise" for a fair
    % comparison with the rest of methods).
    % (old code => impl_UKF.m)
    
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
        % Method to propagate each sample: Set to any of the dynamics
        % formulation, and configure to the desired integrator.
        dynamics_formulation_and_integrator = mbeDynFormulationLagrangeBaumgarte(); %mbeDynFormulationMatrixR();
        
    end
    
    % Private vars
    properties(Access=private)
        CovPlantNoise; 
        CovMeasurementNoise;
    end
    
    methods(Access=protected)
        
        % Returns the location (indices) in P for the independent
        % coordinates (z,zp,zpp). Zeros means not estimated by this filter.
        function [P_idx_z, P_idx_zp, P_idx_zpp] = get_indices_in_P(me)
            % Structure state vector (and its cov P) in this filter:
            P_idx_z   =  1:me.lenZ;
            P_idx_zp  = (1:me.lenZ) + me.lenZ;
            P_idx_zpp = zeros(1,me.lenZ); % Not estimated
        end
                
        
        % Init the filter (see docs in mbeEstimatorFilterBase)
        function [] = init_filter(me)
            % 1) q,qp,qpp: already set in base class.
            % 2) Initial covariance: 
            me.P = diag([...
                me.initVar_Z*ones(1,me.lenZ), ...
                me.initVar_Zp*ones(1,me.lenZ)]);
            
            me.CovPlantNoise = diag([...
                ones(1,me.lenZ)*me.transitionNoise_Z*(me.dt), ...
                ones(1,me.lenZ)*me.transitionNoise_Zp*(me.dt)]);
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            % Kalman state:     
            X = [me.q(me.iidxs) ; me.qp(me.iidxs)];

            L = length(X);  % length of state vector:
            Nsp = 2*L +1;  % Number of samples

            % Generate deterministic samples (aka "sigma points"):
            X_k = mbeUtils.sigma_sampling_generate(X, me.P);

            % Transition model: propagate the samples: X^{-}_{k+1}
            X_minus_kp1   = zeros(Nsp,L);
            q_minus_kp1   = zeros(Nsp,me.lenQ);
            qp_minus_kp1  = zeros(Nsp,me.lenQ);
            qpp_minus_kp1 = zeros(Nsp,me.lenQ);
            for i=1:Nsp,
                % Create the dynamic states for each sample:
                X_ki = X_k(i,:);
                
                q_i  = me.q; q_i(me.iidxs) = X_ki(1:me.lenZ);
                q_i  = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, q_i);
                qp_i = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, q_i, X_ki( me.lenZ+(1:me.lenZ) )');
                [qpp_i,~] = me.dynamics_formulation_and_integrator.solve_for_accelerations(...
                    me.bad_mech_phys_model,...
                    q_i, qp_i, ...
                    struct() ); % simul context
                
                
                %prototype: function [q_next,qp_next,qpp_next] = integrate_one_timestep(me,model,current_time, dt, q,qp,qpp )
                [q_n,qp_n,qpp_n] = ... %me.mbs_propagate_state(me.q, X_k(i,:));
                    me.dynamics_formulation_and_integrator.integrate_one_timestep(...
                        me.bad_mech_phys_model, ...  % model
                        me.current_run_time, me.dt,...
                        q_i,qp_i,qpp_i );
                    
                % Save as rows:
                q_minus_kp1(i,:)=q_n';
                qp_minus_kp1(i,:)=qp_n';
                qpp_minus_kp1(i,:)=qpp_n';
                
                % Extract X (state vector):
                X_minus_kp1(i,:) = [q_minus_kp1(i,me.iidxs)  qp_minus_kp1(i,me.iidxs)];
            end

            % Recover mean & cov before update (NOTE: only for returning P_less, 
            % it could be dropped...)
            [X_minus,P_minus] = mbeUtils.sigma_sampling_mean_and_cov(X_minus_kp1);
            P_minus = P_minus + me.CovPlantNoise; % This wasn't in Roland's paper but we must consider it!!!

            % If there are is no sensor data, we must return this estimate: 
            if (isempty(obs))
                % No sensor:
                % ---------------------------------------------------
                X_plus = X_minus;
                P_plus = P_minus;
            else
                % Yes, there's sensor info: update estimate:
                % ---------------------------------------------------

                % Propagate samples thru observation function (i.e. the sensor)
                nSensors = length(obs);
                Ys = zeros(Nsp, nSensors ); % Each row is a prediction
                for i=1:Nsp,
                    Ys(i,:) = me.bad_mech_phys_model.sensors_simulate(q_minus_kp1(i,:)',qp_minus_kp1(i,:)',qpp_minus_kp1(i,:)');
                end

                % Compute Cov(y) and Cov(X,Y)        
                [y_mean,Pyy, Pxy] = mbeUtils.sigma_sampling_mean_and_cov(Ys, X_minus_kp1,X_minus);

                Pyy = Pyy + me.CovMeasurementNoise; % This wasn't either in Roland's paper!!!

                K = Pxy/Pyy;
                Innovation = obs-y_mean';

                X_plus = X_minus + (K * Innovation)';
                P_plus = P_minus - K * Pyy*K';
            end

            % MBS
            % Convert state vector "X" back to MBS coordinates "q":
            me.P = P_plus;
            
            me.q(me.iidxs) = X_plus(1:me.lenZ);
            me.q = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, me.q);
            me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_plus( me.lenZ+(1:me.lenZ) )');
            % Solve accels:
            [me.qpp,~] = me.dynamics_formulation_and_integrator.solve_for_accelerations(...
                me.bad_mech_phys_model,...
                me.q, me.qp, ...
                struct() ); % simul context
            
        end % run_filter_iter
        
        
    end
end

