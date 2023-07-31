classdef mbeEstimatorDirect_eEKF_FE < mbeEstimatorFilterBase
    % Direct version of the errorEKF_EJ
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
        F0; % Transition matrix
        X0; % State from previous step (also used for initialization)
        deltaqpp;
        Iz;
        Oz;
        W; 
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
%             lenX = 3*me.lenZ;
%             O3z = zeros(lenX);

%             ContinuousCovPlantNoise = diag([...
%                 ones(1,me.lenZ)*me.transitionNoise_Zp, ...
%                 ones(1,me.lenZ)*me.transitionNoise_Zpp*0,... 
%                 ones(1,me.lenZ)*me.transitionNoise_Zpp/me.dt ... 
%                 ]);
%             ContinouosF = [me.Oz, me.Iz, me.Oz; 
%                            me.Oz, me.Oz, me.Iz;
%                            me.Oz, me.Oz, me.Oz];
%             M = me.dt*[-ContinouosF, ContinuousCovPlantNoise;
%                        O3z, ContinouosF' ];
%             N = expm(M);
%             discreteF = N(lenX+1:2*lenX, lenX+1:2*lenX)';
%             me.CovPlantNoise = discreteF*N(1:lenX, lenX+1:2*lenX);
%             me.CovPlantNoise(2*me.lenZ+(1:me.lenZ),2*me.lenZ+(1:me.lenZ)) =  diag(ones(1,me.lenZ)*me.transitionNoise_Zpp);
            me.CovPlantNoise = diag([...
                zeros(1,me.lenZ), ...
                zeros(1,me.lenZ),... 
                ones(1,me.lenZ)*me.transitionNoise_Zpp]);
            
            % These should be constant for all iterations, so eval them once:
%             Iz = eye(me.lenZ);
%             Oz = zeros(me.lenZ);
            
            me.F0 = [ me.Iz, me.Iz*me.dt, 0.5*me.Iz*me.dt^2;
                     me.Oz, me.Iz, me.Iz*me.dt
                     me.Oz, me.Oz, me.Iz];
            me.W = 1; 
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter .* me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
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
%             qppMB = qppMB + mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model, qMB, qpMB, me.X0( 2*me.lenZ+(1:me.lenZ) ));
%             qppMB = qppMB+me.deltaqpp;
            % 2) Kalman state propagation
            % the state is DOF position and velocity error
            
%              me.X0 = zeros(3*me.lenZ, 1);
            X_minus=[...
                qMB(me.iidxs); ...
                qpMB(me.iidxs);...
                qppMB(me.iidxs)];
            
            [dzpp_dz,dzpp_dzp] = me.bad_mech_phys_model.eval_dzpp_dposvel(me.q,me.qp);
%             me.W = 1; 
%             me.W
%             dzpp_dzp = dzpp_dzp *me.W; 
%             dzpp_dz  = dzpp_dz*me.W; 
            F = me.F0 + ...
                [0.5*dzpp_dz*me.dt^2,   0.5*dzpp_dzp*me.dt^2,   me.Oz;
                dzpp_dz*me.dt,          dzpp_dzp*me.dt,         me.Oz;
                me.Oz,                       me.Oz,             me.Oz];
%             X_minus=F*me.X0;  
            % time update
            P_minus = F*me.P*F' + me.CovPlantNoise; %*me.W;
            
            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp, dh_dzpp] = me.bad_mech_phys_model.sensors_jacob_indep(qMB,qpMB,qppMB); % me.q,me.qp,me.qpp);
                H=[dh_dz, dh_dzp, dh_dzpp]; 
                
                % Kalman gain:
                K = P_minus*H'/(H*P_minus*H'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(qMB,qpMB,qppMB);
                
                me.Innovation = obs-obs_predict;
                X_plus = X_minus + K*me.Innovation;
                me.P = (eye(length(X_minus))-K*H)*P_minus;
            else
                % No sensor:
                X_plus = X_minus;
                me.P = P_minus;
                me.Innovation= [];
            end

            
            % Recover KF -> MBS coordinates
            % ------------------------------
            %Position increments fulfill velocity constraints
            deltaz = X_plus(1:me.lenZ)-qMB(me.bad_mech_phys_model.indep_idxs);
            deltaq = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, qMB, deltaz ); 
            me.q = qMB + deltaq;
            
            %Corrected velocities are calculated by applying the velocity
            %problem from the corrected velocities of the degrees of freedom 
            me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_plus( me.lenZ+(1:me.lenZ) ));
            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                         
%             Inverse dynamics
%             R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,me.q); %qp = R*zp_next;
            deltazpp = X_plus(2*me.lenZ+(1:me.lenZ))-qppMB(me.bad_mech_phys_model.indep_idxs);
%             M_bar = R'*me.bad_mech_phys_model.M*R;
%             Rpzp=mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model,me.q,me.qp, zeros(me.lenZ,1) ); 
%             Q=me.bad_mech_phys_model.eval_forces(me.q,me.qp);
%             Q_bar = R'*(Q-me.bad_mech_phys_model.M*Rpzp); 
%             
%             me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) = me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) + M_bar*zpp-Q_bar; % Motor forces at the DOFs
%             Alternative
            
            R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,me.q);
%             R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,qMB);
            M_bar = R'*me.bad_mech_phys_model.M*R;

            %             me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) = me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) + delta_Q_bar; %M_bar* X_plus(2*me.lenZ+(1:me.lenZ));
            deltaQ = M_bar* deltazpp;
            me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) = me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) + deltaQ;
%             norm_deltaQ = norm(deltaQ);
%             norm_Qm = norm(me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs));
%             me.W = 1;%(me.W + norm_deltaQ/(norm_Qm+eps))/2;
      
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             


%                     R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,me.q); %R matrix after position correction
%                     zpp = qppMB(me.bad_mech_phys_model.indep_idxs)+X_plus(2*me.lenZ+(1:me.lenZ)); % zpp corrected
%                     M_bar = R'*me.bad_mech_phys_model.M*R; % Mbar after position correction
%                     Sc=mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model,me.q,me.qp, zeros(me.lenZ,1) ); % Velocity-dependent accelerations after corrections
% %                     Q=me.bad_mech_phys_model.eval_forces(qMB,qpMB); % Forces BEFORE the corrections WATCH OUT!!
%                     Q=me.bad_mech_phys_model.eval_forces(me.q,me.qp); % Forces AFTER the corrections WATCH OUT!!
%                     Q_bar = R'*(Q-me.bad_mech_phys_model.M*Sc);
%                     deltaQ = M_bar*zpp-Q_bar;
%                     norm_deltaQ = norm(deltaQ);
%                     norm_Qm = norm(me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs));
%                     me.W = norm_deltaQ/(norm_deltaQ+norm_Qm);
%                     me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) = me.bad_mech_phys_model.Qm(me.bad_mech_phys_model.indep_idxs) + deltaQ; % Motor forces at the DOFs

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
            

            
%             Acceleration correction
            
            me.qpp = mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model, me.q, me.qp, X_plus(2*me.lenZ+(1:me.lenZ)));
%             me.qpp = qppMB; 
%             me.X0 = zeros(3*me.lenZ,1);
%             me.X0(1:2*me.lenZ) = 0; % = zeros(3*me.lenZ,1);
%             me.X0(2*me.lenZ+(1:me.lenZ)) = X_plus(2*me.lenZ+(1:me.lenZ));
       
        end % run_filter_iter
        
    end
end

