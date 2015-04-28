classdef mbeEstimatorDEKF < mbeEstimatorFilterBase
    % Extended Kalman filter (EKF), forward Euler Integrator
    
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
        % (Is this worth the computational cost? Probably not, but this is how EKF works in theory) 
        % Reflect the new X_less into q,qp,qpp for accurately simulate the
        % sensor & its linearization point:
        eval_sensors_at_X_less = 1;
        dataset;
    end
    
    % Private vars
    properties(Access=private)
        CovPlantNoise; CovMeasurementNoise;
        F; %G; % Transition matrices
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
                ones(1,me.lenZ)*me.transitionNoise_Z*me.dt, ...
                ones(1,me.lenZ)*me.transitionNoise_Zp*me.dt]);
            
            % These should be constant for all iterations, so eval them once:
            me.F = [eye(me.lenZ) eye(me.lenZ)*me.dt;zeros(me.lenZ) eye(me.lenZ)];
            %me.G = [zeros(me.lenZ) ; eye(me.lenZ)*me.dt];
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)

            % time update
            P_less = me.F*(me.P)*me.F' + me.CovPlantNoise;
            
            % Simplified from:  
            %    X_less = me.F*X + me.G*U;
            %    % Kalman state:     
            %    X = [me.q(me.iidxs) ; me.qp(me.iidxs)];
            %    U = me.qpp(me.iidxs); % system input 
            % Transition model (Euler integration)
            X_less=[...
                me.q(me.iidxs) + me.dt*me.qp(me.iidxs); ...
                me.qp(me.iidxs) + me.qpp(me.iidxs)*me.dt ];

            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                if (me.eval_sensors_at_X_less)
                    me.q(me.iidxs) = X_less(1:me.lenZ);
                    me.q = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, me.q);
                    me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_less( me.lenZ+(1:me.lenZ) ));                
                    [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                        me.bad_mech_phys_model,...
                        me.q, me.qp, ...
                        struct() ); % simul context
                end
                
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp] = me.bad_mech_phys_model.sensors_jacob_indep(me.q,me.qp,me.qpp);
                H=[dh_dz, dh_dzp];   % Was: H = [dh_dq(:,me.iidxs) , dh_dqp(:,me.iidxs) ];
                
                % Kalman gain:
                K = P_less*H'/(H*P_less*H'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(me.q,me.qp,me.qpp);
                
                Innovation = obs-obs_predict;
                X_plus = X_less + K*Innovation;
                %     P = (eye(2)-K*H)*P_less*(eye(2)-K*H)'+K*CovMeasurementNoise*K';
                me.P = (eye(length(X_less))-K*H)*P_less;
            else
                % No sensor:
                X_plus = X_less;
                me.P = P_less;
            end
            
            % Recover KF -> MBS coordinates
            % ------------------------------
            me.q(me.iidxs) = X_plus(1:me.lenZ);
            me.q = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, me.q);
            me.qp = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_plus( me.lenZ+(1:me.lenZ) ));
            % Solve accels:
            [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                me.bad_mech_phys_model,...
                me.q, me.qp, ...
                struct() ); % simul context
       
        end % run_filter_iter
        
    end
end

