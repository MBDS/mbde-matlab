classdef mbeEstimatorSCKF < mbeEstimatorFilterBase
    % Smoothly constrained Kalman filter (EKF) with "perfect measurements", forward Euler Integrator
    
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
        MAX_ERR = 1e-3; % Error in constraints
%         INCR_X_PLUS = 10*MAX_INCR_X;
        MAX_ITERS = 200;
        lenq;
        Iq; % Identity matrix (size of q)
        Ix; % Identity matrix (size of X)
        Oq; % Null matrix

        alpha = 0.1; % this parameter have influence over the fictious noise of constraints sensor.
        beta = 1; %5; % as this parameter is bigger, the fictious noise of the sensors decreases faster.
%         PM_VAR  = 1e-6;   % Perfect measurement precision
        
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
            me.Iq = eye(me.lenq);
            me.Oq = zeros(me.lenq);
            lenX = 2*me.lenq;
            Ox = zeros(lenX);
            me.Ix = eye(lenX);
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
%             me.CovPlantNoise = diag([...
%                 ones(1,me.lenq)*me.transitionNoise_Z*me.dt, ...
%                 ones(1,me.lenq)*me.transitionNoise_Zp*me.dt]);

            ContinuousCovPlantNoise = diag([...
                ones(1,me.lenq)*me.transitionNoise_Zp, ...
                ones(1,me.lenq)*me.transitionNoise_Zpp]);
            factor = 1;
             ContinuousCovPlantNoise([me.iidxs, me.iidxs+me.lenq],[me.iidxs, me.iidxs+me.lenq])= factor*ContinuousCovPlantNoise([me.iidxs, me.iidxs+me.lenq],[me.iidxs, me.iidxs+me.lenq]);
             ContinuousCovPlantNoise= ContinuousCovPlantNoise/factor;
            ContinouosF = [me.Oq, me.Iq; 
                           me.Oq, me.Oq];
            M = me.dt*[-ContinouosF, ContinuousCovPlantNoise;
                       Ox, ContinouosF' ];
            N = expm(M);
            discreteF = N(lenX+1:2*lenX, lenX+1:2*lenX)';
            me.CovPlantNoise = discreteF*N(1:lenX, lenX+1:2*lenX);
            
            
            % These should be constant for all iterations, so eval them once:
            me.F = [eye(me.lenq) eye(me.lenq)*me.dt;zeros(me.lenq) eye(me.lenq)];
            %me.G = [zeros(me.lenq) ; eye(me.lenq)*me.dt];
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            % Transition model (Euler integration)
            X_minus=[...
                me.q + me.dt*me.qp; ...
                me.qp + me.qpp*me.dt ];
            P_minus = me.F*(me.P)*me.F' + me.CovPlantNoise;
            % Correction from sensors measurements
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dq , dh_dqp] = me.bad_mech_phys_model.sensors_jacob(X_minus(1:me.lenq),X_minus(me.lenq+(1:me.lenq)),me.qpp); % me.q,me.qp,me.qpp);
                
                H=[dh_dq, dh_dqp];
                
                % Kalman gain:
                K = P_minus*H'/(H*P_minus*H'+me.CovMeasurementNoise);

                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(X_minus(1:me.lenq),X_minus(me.lenq+(1:me.lenq)),me.qpp);
                
                me.Innovation = obs-obs_predict;
                X_plus = X_minus + K*me.Innovation;
                P_plus = (me.Ix-K*H)*P_minus;
            else
                % No sensor:
                X_plus = X_minus;
                P_plus = P_minus;
                me.Innovation = [];
            end
            
            % Iterative application of the constraints: perfect
            % measurements with weakening matrix (i.e. added virtual noise)
            % to ease the convergence of the nonlinear cocnstraints.
            me.q = X_plus(1:me.lenq);
            me.qp = X_plus(me.lenq+(1:me.lenq));
            [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                        me.bad_mech_phys_model,...
                        me.q, me.qp, ...
                        struct());
            phiq = me.bad_mech_phys_model.jacob_phi_q(me.q);  
            phix11 = phiq;
            phix12 = 0*phix11; % this should be faster than creating a square matrix of 0s.
            phix22 = phix11;
            phix21 = me.bad_mech_phys_model.Phiq_times_qp_q(me.q,me.qp); % Partial derivative of velocity constraints wrt q
            Phi_X = [phix11,phix12; phix21,phix22];
            xi = me.alpha*Phi_X*P_plus*Phi_X';
            phi = me.bad_mech_phys_model.phi(me.q);
            PHI = [phi;
                   phiq*me.qp];
            err= norm(PHI);
            ITERS = 0;
            while (err>me.MAX_ERR && ITERS < me.MAX_ITERS)
                ITERS = ITERS + 1;
                if (ITERS == me.MAX_ITERS)
                    warning('SCKF: Max iters reached')
                end
                K = P_plus*Phi_X'/(Phi_X*P_plus*Phi_X'+xi);
                X_plus = X_plus-K*PHI;
%                 P_plus = (me.Ix-K*Phi_X)*P_plus*(me.Ix-K*Phi_X)'+K*xi*K';
                P_plus = (me.Ix-K*Phi_X)*P_plus;
                xi = xi*exp(-me.beta);
                
                me.q = X_plus(1:me.lenq);
                me.qp = X_plus(me.lenq+(1:me.lenq));
                [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                        me.bad_mech_phys_model,...
                        me.q, me.qp, ...
                        struct());
                    
                phiq = me.bad_mech_phys_model.jacob_phi_q(me.q);  
                phix11 = phiq;
                phix12 = 0*phix11;
                phix22 = phix11;
                phix21 = me.bad_mech_phys_model.Phiq_times_qp_q(me.q,me.qp); % Partial derivative of velocity constraints wrt q
                Phi_X = [phix11,phix12; phix21,phix22];
                phi = me.bad_mech_phys_model.phi(me.q);
                PHI = [phi;
                       phiq*me.qp];
                err= norm(PHI);
            end
            me.P = P_plus;
                
        end % run_filter_iter
        
    end
end

