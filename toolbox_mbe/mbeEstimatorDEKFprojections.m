classdef mbeEstimatorDEKFprojections < mbeEstimatorFilterBase
    % Discrete Iterated Extended Kalman filter (EKF) with "perfect measurements", forward Euler Integrator
    
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
        MAX_INCR_X = 1e-5; % Should it be lower? 1e-10?
%         INCR_X_PLUS = 10*MAX_INCR_X;
        MAX_IEKF_ITERS = 20;
        lenq;
        FULL_COV_SENSORS;
        PM_VAR  = 1e-8;   % Perfect measurement precision
        
    end
    
    methods(Access=protected)
        
        % Returns the location (indices) in P for the independent
        % coordinates (z,zp,zpp). Zeros means not estimated by this filter.
        function [P_idx_z, P_idx_zp, P_idx_zpp] = get_indices_in_P(me)
            % Structure state vector (and its cov P) in this filter:
            P_idx_z   =  me.iidxs;
            P_idx_zp  = me.iidxs + me.lenq;
            P_idx_zpp = zeros(1,me.lenq); % Not estimated
        end
                
        
        % Init the filter (see docs in mbeEstimatorFilterBase)
        function [] = init_filter(me)
            % 1) q,qp,qpp: already set in base class.
            % 2) Initial covariance: the covariances are propagated from
            % the independent coordinates to the dependent coordinates
            % using the R matrix %FIXME: this method is not suitable for
            % the plant noise, because provides different values depending
            % on the initial position
            me.lenq = length(me.q);
%             PZ = diag(me.initVar_Z*ones(1,me.lenZ));
%             PZp = diag(me.initVar_Zp*ones(1,me.lenZ));
%             R = mbeKinematicsSolver.calc_R_matrix(me.mech_phys_model,me.q);
%             Pq = R*PZ*R'; % me.PM_VAR*ones(size(me.q)); 
%             Pqp = R*PZp*R';
%             me.P = zeros(2*me.lenq);
%             % Initial assembly of P. A small number is added to ensure the matrix is positive definite
%             me.P(1:me.lenq,1:me.lenq) = Pq+1e-5; 
%             me.P((me.lenq+1):2*me.lenq,(me.lenq+1):2*me.lenq) = Pqp+1e-5;
            
            me.P = diag(...
                [me.initVar_Z*ones(1,me.lenq),...
                me.initVar_Zp*ones(1,me.lenq)]);
            % Covariance of plant noise in independent coordinates (the
            % covariance of the plant with dependent coordinates will be
            % calculated every time step in the main loop of the filter)
            me.CovPlantNoise = diag([...
                ones(1,me.lenq)*me.transitionNoise_Z*me.dt, ...
                ones(1,me.lenq)*me.transitionNoise_Zp*me.dt]);
            
            
            % These should be constant for all iterations, so eval them once:
            me.F = [eye(me.lenq) eye(me.lenq)*me.dt;zeros(me.lenq) eye(me.lenq)];
            %me.G = [zeros(me.lenq) ; eye(me.lenq)*me.dt];
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)

            % time update. The covariance of the noise of the plant is
            % updated every time step
%             R = mbeKinematicsSolver.calc_R_matrix(me.mech_phys_model,me.q);
%             diagR = [R,zeros(size(R)); 
%                 zeros(size(R)), R];
%             CovPlantNoiseq = diagR*me.CovPlantNoise*diagR'+1e-4;
%             P_less = me.F*(me.P)*me.F' + CovPlantNoiseq;
            P_minus = me.F*(me.P)*me.F' + me.CovPlantNoise;
            
            % Simplified from:  
            %    X_less = me.F*X + me.G*U;
            %    % Kalman state:     
            %    X = [me.q(me.iidxs) ; me.qp(me.iidxs)];
            %    U = me.qpp(me.iidxs); % system input 
            % Transition model (Euler integration)
            X_minus=[...
                me.q + me.dt*me.qp; ...
                me.qp + me.qpp*me.dt ];

           
            me.q = X_minus(1:me.lenq);
            me.qp = X_minus(me.lenq+(1:me.lenq));
            [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                    me.bad_mech_phys_model,...
                    me.q, me.qp, ...
                    struct());
            if (~isempty(obs))
                [dh_dq, dh_dqp]=me.bad_mech_phys_model.sensors_jacob(me.q,me.qp,me.qpp);
                H = [dh_dq, dh_dqp];
                    

                obs_predict = me.bad_mech_phys_model.sensors_simulate(me.q,me.qp,me.qpp);
                Innovation = obs-obs_predict;
                K = P_minus*H'/(H*P_minus*H'+me.CovMeasurementNoise);
                X_plus = X_minus+K*Innovation;
                me.P = (eye(length(X_plus))-K*H)*P_minus;

            else
                X_plus = X_minus;
                me.P = P_minus;
            end
            
            % Projections
            me.q = X_plus(1:me.lenq);
            me.qp = X_plus(me.lenq+(1:me.lenq));
            phi = me.bad_mech_phys_model.phi(me.q);
            error = 1;
%             W = me.P(1:me.lenq,1:me.lenq);
            W_factor = 1e3;
            eps_sum = 1e-4;
            dh_dq = dh_dq+eps_sum;
            dh_dqp = dh_dqp+eps_sum;
            W = inv(me.P(1:me.lenq,1:me.lenq))+dh_dq'*W_factor*dh_dq;
            Sigma = inv(W);
%             W(me.iidx,me.iidx) = 1e3*W(me.iidx,me.iidx);  % Los pesos se
%             calcularan como K*innovation*numero_gordo, de esta manera las
%             variables que son afectadas por los sensores tienen más peso.
%             El inconveniente es que si una variable es afectada por un
%             sensor, pero coincide que está "en el sitio", su innovación
%             será cero, y es posible que se vea afectada por el proceso de
%             proyección, cuando no debería.  (mejor tratar de emplear la
%             H?)
%             
            while error > 1e-5
                D = me.bad_mech_phys_model.jacob_phi_q(me.q);
                d = -phi+D*me.q;
                me.q = me.q-Sigma*D'/(D*Sigma*D')*(D*me.q-d);
                phi = me.bad_mech_phys_model.phi(me.q);
                error = norm(phi);
            end
%             me.P(1:me.lenq,1:me.lenq) = me.P(1:me.lenq,1:me.lenq)-W*D'/(D*W*D')*D*W; %% Where this come from?? Is it correct as is?
            J = Sigma*D'/(D*Sigma*D')*D;
%             I = eye(me.lenq);
%             me.P(1:me.lenq,1:me.lenq) = (I-J)*me.P(1:me.lenq,1:me.lenq)*(I-J)';
            me.P(1:me.lenq,1:me.lenq) = Sigma-J*Sigma;
            
            % Velocities projections
%             jac = me.bad_mech_phys_model.jacob_phi_q(me.q);
%             A = me.bad_mech_phys_model.M + me.dt^2/4*jac'*me.bad_mech_phys_model.alpha*jac;
%             b = me.bad_mech_phys_model.M*me.qp;
%             me.qp = A\b;
            % Covariance for velocities????
            W = inv(me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq))+dh_dqp'*W_factor*dh_dqp;
            Sigma = inv(W);
            D = me.bad_mech_phys_model.jacob_phi_q(me.q);
            me.qp = me.qp-Sigma*D'/(D*Sigma*D')*(D*me.qp);
            J = Sigma*D'/(D*Sigma*D')*D;
%             me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq) = (I-J)*me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq)*(I-J)';
%             me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq) = me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq)-W*D'/(D*W*D')*D*W; %% Where this come from?? Is it correct as is?
            me.P((1:me.lenq)+me.lenq,(1:me.lenq)+me.lenq) = Sigma-J*Sigma;
            % Solve accels:
            [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                me.bad_mech_phys_model,...
                me.q, me.qp, ...
                struct() ); % simul context
       
        end % run_filter_iter
        
    end
end

