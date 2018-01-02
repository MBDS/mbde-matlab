classdef mbeEstimatorIncrManifold_full_jac < mbeEstimatorFilterBase
    % This function runs a filter which tries to estimate the errors between
    % the real system and the model, so the model is a conventional MB
    % simulation, wich can have its own integrator, formulation, etc, while
    % the model for the errors remains as simple as possible. Unlike the
    % conventional errorEKF (mbeEstimatorIncrManifold), this version 
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
    properties(Access=private)
        CovPlantNoise; 
        CovMeasurementNoise;
        F0; % Constant part of transition matrix
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
            % These should be constant for all iterations, so eval them once:
            
            Iz = eye(me.lenZ);
            Oz = zeros(me.lenZ);
            lenX = 2*me.lenZ;
            O2z = zeros(lenX);
            
            % 2) Initial covariance: 
            me.P = diag([...
                me.initVar_Z*ones(1,me.lenZ), ...
                me.initVar_Zp*ones(1,me.lenZ)]);
            
%             me.CovPlantNoise = diag([...
%                 ones(1,me.lenZ)*me.transitionNoise_Z*me.dt, ...
%                 ones(1,me.lenZ)*me.transitionNoise_Zp*me.dt]);
%             Discrete covariance from continouos (Van Loan's method)
            ContinuousCovPlantNoise = diag([...
                ones(1,me.lenZ)*me.transitionNoise_Zp, ...
                ones(1,me.lenZ)*me.transitionNoise_Zpp]);
            ContinouosF = [Oz, Iz; 
                            Oz, Oz];
            M = me.dt*[-ContinouosF, ContinuousCovPlantNoise;
                       O2z, ContinouosF' ];
            N = expm(M);
            discreteF = N(lenX+1:2*lenX, lenX+1:2*lenX)';
            me.CovPlantNoise = discreteF*N(1:lenX, lenX+1:2*lenX);
            % Transition matrix
            me.F0 = [ Iz, Iz*me.dt;
                     Oz, Iz ];
                 
            
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter .* me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            
            % 1) Run dynamics: 
            [qMB,qpMB,qppMB] = me.dynamics_formulation_and_integrator.integrate_one_timestep(...
                me.bad_mech_phys_model,...
                me.current_run_time, ...
                me.dt, ...
                me.q, me.qp, me.qpp );

            % 2) Kalman state propagation
            % the state is DOF position and velocity error
            X_less=zeros(2*me.lenZ, 1);
            
            % time update
            [dzpp_dz,dzpp_dzp] = me.bad_mech_phys_model.eval_dzpp_dposvel(me.q,me.qp);
            F = me.F0+[0.5*dzpp_dz*me.dt^2, 0.5*dzpp_dzp*me.dt^2; 
                        dzpp_dz*me.dt, dzpp_dzp*me.dt];
            
            P_less = F*me.P*F' + me.CovPlantNoise;
            
            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp] = me.bad_mech_phys_model.sensors_jacob_indep(qMB,qpMB,qppMB); % me.q,me.qp,me.qpp);
                H=[dh_dz, dh_dzp];
                
                % Kalman gain:
                K = P_less*H'/(H*P_less*H'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(qMB,qpMB,qppMB);
                
                me.Innovation = obs-obs_predict;
                X_plus = X_less + K*me.Innovation;
                me.P = (eye(length(X_less))-K*H)*P_less;
            else
                % No sensor:
                X_plus = X_less;
                me.P = P_less;
                me.Innovation = [];
            end
% % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %             %%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%
% % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             global data_saved
%             data_saved.F = [data_saved.F,me.F];
%             data_saved.H = [data_saved.H,H]; % it works only without multirate
% 
% data_saved.inn{end+1} = me.Innovation;
%                       data_saved.Hcell{end+1} = H;
%                       data_saved.Pcell{end+1} = P_less;
%                       data_saved.Kcell{end+1} = K;
% % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %             %%%%%%%%%%%%%%%% end test %%%%%%%%%%%%%%
% % % % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Recover KF -> MBS coordinates
            % ------------------------------
            %Position increments fulfill velocity constraints
            deltaq = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, qMB, X_plus(1:me.lenZ) ); 
            me.q = qMB + deltaq;
% % % % %             nphi = norm(me.bad_mech_phys_model.phi(me.q));
% % % % %             data_saved.nphi = [data_saved.nphi, nphi];
% % % % %             if nphi > 2e-3
% % % % %                 warning('position_tol')
% % % % % %                 pause
% % % % %             end
% % % % % %             
            
            %Velocity increments also fulfill velocity constraints
            me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, qpMB(me.bad_mech_phys_model.indep_idxs)+X_plus( me.lenZ+(1:me.lenZ) ));
            
            % Accel is already computed:
            
%             [qppMB,zpp,q,qp] =me.dynamics_formulation_and_integrator.solve_for_accelerations( me.bad_mech_phys_model,me.q,me.qp);
            
            me.qpp = qppMB;
       
        end % run_filter_iter
        
    end
end

