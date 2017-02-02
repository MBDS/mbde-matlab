classdef mbeEstimatorDIEKF_pm < mbeEstimatorFilterBase
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
        MAX_INCR_X = 1e-8; % Should it be lower? 1e-10?
        MAX_ERR = 1e-3; % Not always achievable with this method
%         INCR_X_PLUS = 10*MAX_INCR_X;
        MAX_IEKF_ITERS = 2000;
        lenq;
        FULL_COV_SENSORS;
        PM_VAR  = 1e-3;   % Perfect measurement variance
        
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
            Iz = eye(me.lenq);
            Oz = zeros(me.lenq);
            lenX = 2*me.lenq;
            O2z = zeros(lenX);
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
            ContinouosF = [Oz, Iz; 
                            Oz, Oz];
            M = me.dt*[-ContinouosF, ContinuousCovPlantNoise;
                       O2z, ContinouosF' ];
            N = expm(M);
            discreteF = N(lenX+1:2*lenX, lenX+1:2*lenX)';
            me.CovPlantNoise = discreteF*N(1:lenX, lenX+1:2*lenX);
            
            
            % These should be constant for all iterations, so eval them once:
            me.F = [eye(me.lenq) eye(me.lenq)*me.dt;zeros(me.lenq) eye(me.lenq)];
            %me.G = [zeros(me.lenq) ; eye(me.lenq)*me.dt];
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
%             me.CovMeasurementNoise = diag(sensors_stds.^2);
            me.FULL_COV_SENSORS = diag([sensors_stds'.^2, me.PM_VAR*ones(1,2*(me.lenq-me.lenZ)) ]);
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
            P_less = me.F*(me.P)*me.F' + me.CovPlantNoise;
            
            % Simplified from:  
            %    X_less = me.F*X + me.G*U;
            %    % Kalman state:     
            %    X = [me.q(me.iidxs) ; me.qp(me.iidxs)];
            %    U = me.qpp(me.iidxs); % system input 
            % Transition model (Euler integration)
            X_less=[...
                me.q + me.dt*me.qp; ...
                me.qp + me.qpp*me.dt ];

           
            incr_X_plus = 1;
            err = 1; 
            IEKF_ITERS = 0;
            X_plus = X_less;
%             while (incr_X_plus>me.MAX_INCR_X && IEKF_ITERS < me.MAX_IEKF_ITERS)
            while (err>me.MAX_ERR && IEKF_ITERS < me.MAX_IEKF_ITERS)
                IEKF_ITERS = IEKF_ITERS + 1;
                me.q = X_plus(1:me.lenq);
                me.qp = X_plus(me.lenq+(1:me.lenq));
                [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                        me.bad_mech_phys_model,...
                        me.q, me.qp, ...
                        struct());
                phiq = me.bad_mech_phys_model.jacob_phi_q(me.q);  
                H21 = phiq;
                H32 = H21;
                H31 = me.bad_mech_phys_model.Phiq_times_qp_q(me.q,me.qp); % Partial derivative of velocity constraints wrt q
                phi = me.bad_mech_phys_model.phi(me.q);
                if (~isempty(obs))
                    [dh_dq, dh_dqp]=me.bad_mech_phys_model.sensors_jacob(me.q,me.qp,me.qpp);
                    H = [dh_dq, dh_dqp;
                        H21, zeros(size(H21));
                        H31, H32];
                    sensors_predict = me.bad_mech_phys_model.sensors_simulate(me.q,me.qp,me.qpp);
                    obs_predict = [sensors_predict;
                                   phi;
                                   phiq*me.qp];
                    full_obs = [obs; zeros(2*length(phi),1)];
                    COV_SENSORS = me.FULL_COV_SENSORS;
                    COV_SENSORS((length(obs)+1):(length(obs)+2*length(phi)),(length(obs)+1):(length(obs)+2*length(phi)))=COV_SENSORS((length(obs)+1):(length(obs)+2*length(phi)),(length(obs)+1):(length(obs)+2*length(phi)))*exp(-IEKF_ITERS);
                    me.Innovation = obs-sensors_predict;
                else
                    H = [H21, zeros(size(H21));
                        H31, H32];
                    
                    obs_predict = [phi;
                                   phiq*me.qp];
                    full_obs = zeros(2*length(phi),1);
                    
                    COV_SENSORS = me.FULL_COV_SENSORS((length(obs)+1):(length(obs)+2*length(phi)),(length(obs)+1):(length(obs)+2*length(phi)))*exp(-IEKF_ITERS);
                    me.Innovation = [];
                end
                err = norm([phi;phiq*me.qp]);
                locInnovation = full_obs-obs_predict;
                K_i = P_less*H'/(H*P_less*H'+COV_SENSORS);
                X_plus_new = X_less+K_i*(locInnovation-H*(X_less-X_plus));
                incr_X_plus = norm(X_plus-X_plus_new);
                X_plus = X_plus_new;
            end
            me.P = (eye(length(X_plus))-K_i*H)*P_less;
%             I_KH = (eye(length(X_plus))-K_i*H);
%             me.P = I_KH*P_less*I_KH'+K_i*COV_SENSORS*K_i'; % Joseph Form
            IEKF_ITERS
            if(norm(phi)>1e-3)
                warning('convergence not achieved')
                phi
                pause();
            end
            
            % Recover KF -> MBS coordinates
            % ------------------------------
            me.q = X_plus(1:me.lenq);
            me.qp = X_plus(me.lenq+(1:me.lenq));
            % Solve accels:
            [me.qpp,~] = me.post_iter_accel_solver.solve_for_accelerations(...
                me.bad_mech_phys_model,...
                me.q, me.qp, ...
                struct() ); % simul context
       
        end % run_filter_iter
        
    end
end

