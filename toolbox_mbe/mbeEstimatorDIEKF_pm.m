classdef mbeEstimatorDIEKF_pm < mbeEstimatorFilterBase
    % Discrete Iterated Extended Kalman filter (EKF) with "perfect measurements", forward Euler Integrator
    
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
            IEKF_ITERS = 0;
            X_plus = X_less;
            while (incr_X_plus>me.MAX_INCR_X && IEKF_ITERS < me.MAX_IEKF_ITERS)
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
                    
                    obs_predict = [me.bad_mech_phys_model.sensors_simulate(me.q,me.qp,me.qpp);
                                   phi;
                                   phiq*me.qp];
                    full_obs = [obs; zeros(2*length(phi),1)];
                    COV_SENSORS = me.FULL_COV_SENSORS;
                else
                    H = [H21, zeros(size(H21));
                        H31, H32];
                    
                    obs_predict = [phi;
                                   phiq*me.qp];
                    full_obs = zeros(2*length(phi),1);
                    
                    COV_SENSORS = me.FULL_COV_SENSORS((length(obs)+1):(length(obs)+2*length(phi)),(length(obs)+1):(length(obs)+2*length(phi)));
                end
                Innovation = full_obs-obs_predict;
                K_i = P_less*H'/(H*P_less*H'+COV_SENSORS);
                X_plus_new = X_less+K_i*(Innovation-H*(X_less-X_plus));
                incr_X_plus = norm(X_plus-X_plus_new);
                X_plus = X_plus_new;
            end
            me.P = (eye(length(X_plus))-K_i*H)*P_less;

            
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

