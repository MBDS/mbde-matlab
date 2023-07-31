classdef mbeEstimatorIncrManifold_AerrorEKF < mbeEstimatorFilterBase
    % This function runs a filter which tries to estimate the errors between
    % the real system and the model, so the model is a conventional MB
    % simulation, wich can have its own integrator, formulation, etc, while
    % the model for the errors remains as simple as possible. 
    %
    % This filter estimates the states and forces, together with the noise 
    % covariance matrices, in what is known as adaptive Kalman filtering.
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
        
    end
    
    % Private vars
    properties(Access=public)
        CovPlantNoise; 
        CovMeasurementNoise;
        CovPlantNoise_hist;
        F0; % Transition matrix
        X0; % State from previous step (also used for initialization)
        deltaqpp;
        Iz;
        Oz;
        W; 
        % For the adaptive filter
        deltaX;
        deltaInn;
        stack_CovPlantNoise;
        stack_CovMeasurementNoise;
        P_ant;
        cont = 0;
        steps = 0;
        max_steps = 50; % Length of the innovation window
        F_ant;
        H_ant;
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

            % The initial covariance plant noise is set to the default
            % value, which offers a good performance of the errorEKF. It
            % can be modified to see the behaviour of the adaptive filter
            % with different initial covariance plant noises.
            me.CovPlantNoise = diag([...
              zeros(1,me.lenZ), ...
              zeros(1,me.lenZ),... 
              ones(1,me.lenZ)*me.transitionNoise_Zpp]);
            
            me.F0 = [ me.Iz, me.Iz*me.dt, 0.5*me.Iz*me.dt^2;
                     me.Oz, me.Iz, me.Iz*me.dt
                     me.Oz, me.Oz, me.Iz];
            me.W = 1; 
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter .* me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
            me.X0 = zeros(3*me.lenZ, 1);
            me.deltaqpp = 0;

            me.stack_CovPlantNoise = zeros(numel(me.CovPlantNoise), me.max_steps);
            me.stack_CovMeasurementNoise = zeros(numel(me.CovMeasurementNoise), me.max_steps);
            
            me.steps = 0;
            me.CovPlantNoise_hist = [];
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
            X_minus=F*me.X0;  
            % time update
            P_minus = F*me.P*F' + me.CovPlantNoise; 
                       
            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp, dh_dzpp] = me.bad_mech_phys_model.sensors_jacob_indep(qMB,qpMB,qppMB);
                H=[dh_dz, dh_dzp, dh_dzpp]; 
                
                % Kalman gain:
                K = P_minus*H'/(H*P_minus*H'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(qMB,qpMB,qppMB);
                
                me.Innovation = obs-obs_predict;
                                
                X_plus = X_minus + K*me.Innovation;
                me.P = (eye(length(X_minus))-K*H)*P_minus;
                
                % Adaptive equations
                % NOTE: Uncomment the lines referring CovMeasurementNoise 
                % for estimating also the measurement noise covariance 
                % matrix
                if me.steps > 0
                    accum_CovPlantNoise = zeros(size(me.CovPlantNoise));
%                     accum_CovMeasurementNoise = zeros(size(me.CovMeasurementNoise));
                    
                    actual_CovPlantNoise = X_plus*X_plus' + me.P - me.F_ant*me.P_ant*me.F_ant';                    
%                     actual_CovMeasurementNoise = me.Innovation*me.Innovation' - H*P_minus*H';
                    
                    if me.steps < me.max_steps
                        me.stack_CovPlantNoise(:,me.steps) = actual_CovPlantNoise(:);
%                         me.stack_CovMeasurementNoise(:,me.steps) = actual_CovMeasurementNoise(:);
                    else
                        me.stack_CovPlantNoise(:,1:me.max_steps-1)=me.stack_CovPlantNoise(:,2:me.max_steps);
                        me.stack_CovPlantNoise(:,me.max_steps) = actual_CovPlantNoise(:);
%                         me.stack_CovMeasurementNoise(:,1:me.max_steps-1)=me.stack_CovMeasurementNoise(:,2:me.max_steps);
%                         me.stack_CovMeasurementNoise(:,me.max_steps) = actual_CovMeasurementNoise(:);
                    end
                    
                    if me.steps < me.max_steps
                        for i = 1:me.steps
                            accum_CovPlantNoise = accum_CovPlantNoise + reshape(me.stack_CovPlantNoise(:,i),size(me.CovPlantNoise,1),size(me.CovPlantNoise,2));
%                             accum_CovMeasurementNoise = accum_CovMeasurementNoise + reshape(me.stack_CovMeasurementNoise(:,i),size(me.CovMeasurementNoise,1),size(me.CovMeasurementNoise,2));
                        end
                        me.CovPlantNoise = accum_CovPlantNoise/me.steps;
%                         me.CovMeasurementNoise = accum_CovMeasurementNoise/me.steps;
%                         me.CovMeasurementNoise = abs(diag(diag(me.CovMeasurementNoise))); 
                    else
                        for i = 1:me.max_steps
                            accum_CovPlantNoise = accum_CovPlantNoise + reshape(me.stack_CovPlantNoise(:,i),size(me.CovPlantNoise,1),size(me.CovPlantNoise,2));
%                             accum_CovMeasurementNoise = accum_CovMeasurementNoise + reshape(me.stack_CovMeasurementNoise(:,i),size(me.CovMeasurementNoise,1),size(me.CovMeasurementNoise,2));
                        end
                        me.CovPlantNoise = accum_CovPlantNoise/me.max_steps;
%                         me.CovMeasurementNoise = accum_CovMeasurementNoise/me.max_steps;
%                         me.CovMeasurementNoise = abs(diag(diag(me.CovMeasurementNoise)));
                    end
                end
                
                aux = zeros(3*me.lenZ,3*me.lenZ);
                aux(2*me.lenZ+1:end,2*me.lenZ+1:end) = eye(me.lenZ);
                me.CovPlantNoise = abs(me.CovPlantNoise);
                me.CovPlantNoise = aux*me.CovPlantNoise*aux';
                me.CovPlantNoise = diag(diag(me.CovPlantNoise)); 
                
                me.steps = me.steps + 1;
                
                me.P_ant = me.P;
                me.F_ant = F;
                                                                   
            else
                % No sensor:
                X_plus = X_minus;
                me.P = P_minus;
                me.Innovation= [];
            end
            
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

