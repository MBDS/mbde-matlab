classdef mbeEstimatorDIEKFacc < mbeEstimatorFilterBase
    % Iterated Extended Kalman filter (EKF) with virtual acceleration 
    % sensor and forward Euler Integrator
    % Accelerations are predicted using the "post_iter_accel_solver" property
    % of the base class.
    %
    
    properties
        max_xplus_incr_norm = 1e-3;
        max_iekf_iters      = 100;
        
        % Std. deviation of the "virtual observation" that tries to match
        % predicted & real acelerations. Can be a vector to assign
        % different values to each dof.
        assumed_acc_std_noise    = deg2rad(1);        
    end
    
    % Private vars
    properties(Access=private)
        CovPlantNoise; 
        CovMeasurementNoise;
        F; % Transition matrix
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
                me.initVar_Zpp*ones(1,me.lenZ) ]);
            
            me.CovPlantNoise = diag([...
                ones(1,me.lenZ)*me.transitionNoise_Z*me.dt, ...
                ones(1,me.lenZ)*me.transitionNoise_Zp*me.dt, ...
                ones(1,me.lenZ)*me.transitionNoise_Zpp*me.dt ...
                ]);
            
            % These should be constant for all iterations, so eval them once:
            B1 = eye(me.lenZ); B0 = zeros(me.lenZ);
            me.F = [...
                B1, B1*me.dt,       B0; ...
                B0,       B1, B1*me.dt; ...
                B0,       B0,       B1 ];
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            
            % Append uncertainty of "virtual accelerations sensor":
            if (length(me.assumed_acc_std_noise)==1)
                % The same value to all dofs:
                sensors_stds=[sensors_stds; ones(me.lenZ,1)*me.assumed_acc_std_noise ];
            else
                % May have different values for each dof
                sensors_stds=[sensors_stds; me.assumed_acc_std_noise(:)];
            end
            
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            % Kalman state:     
            X = [me.q(me.iidxs) ; me.qp(me.iidxs); me.qpp(me.iidxs)];
            
            % 1) Prediction
            P_minus = me.F * me.P * me.F' + me.CovPlantNoise;
            X_minus = me.F * X;


            % 2) Iterative update
            X_plus = X_minus;
            P_plus = P_minus;
           
            incr_X_plus_norm =10*me.max_xplus_incr_norm;
            IEKF_ITERS = 0;

            % If no sensor data available, P_less should be the output at P
            % variable)
            if (~isempty(obs))
                while(incr_X_plus_norm>me.max_xplus_incr_norm && IEKF_ITERS<me.max_iekf_iters)
                    IEKF_ITERS=IEKF_ITERS+1;

                    % X -> q
                    me.q(me.iidxs) = X_plus(1:me.lenZ);
                    me.q   = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, me.q);
                    me.qp  = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_plus( me.lenZ+(1:me.lenZ) ));
                    me.qpp = mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model, me.q, me.qp, X_plus( 2*me.lenZ+(1:me.lenZ) ));

                    % Eval observation Jacobian wrt the indep. coordinates:
                    [dh_dz , dh_dzp, dh_dzpp] = me.bad_mech_phys_model.sensors_jacob_indep(me.q,me.qp,me.qpp);
                    
                    % "physical" Sensors:
                    sensors_predict = me.bad_mech_phys_model.sensors_simulate(me.q,me.qp,me.qpp);
                    
                    % "virtual sensor" of accs
                    qpp_predict = me.post_iter_accel_solver.solve_for_accelerations(...
                        me.bad_mech_phys_model,...
                        me.q, me.qp, ...
                        struct() ); % simul context
                    
                    
                    % Stack composite observation:
                    obs_predict=[sensors_predict;X_plus( 2*me.lenZ+(1:me.lenZ) )];
                    
                    % H = complete observation jacobian: sensor part + "virtual
                    % observation" of the acceleration.
                    H=[         dh_dz,          dh_dzp,      dh_dzpp; 
                       zeros(me.lenZ),  zeros(me.lenZ),  eye(me.lenZ) ];

                    % Kalman gain:
                    K_i = P_minus*H'/(H*P_minus*H'+me.CovMeasurementNoise);

                    Innovation = [obs;qpp_predict(me.iidxs)] - obs_predict;
                    X_plus_new = X_minus + K_i*( Innovation - H*(X_minus-X_plus) );

                    incr_X_plus_norm = norm(X_plus-X_plus_new);
                    X_plus = X_plus_new;       
                end % while
                
                P_plus = (eye(length(X_plus))-K_i*H)*P_minus;

            else
                % No sensor. Nothing else to do, X_plus & P_plus
                % unmodified.
            end
            
            % Recover KF -> MBS coordinates
            % ------------------------------
            me.P = P_plus;
            
            me.q(me.iidxs) = X_plus(1:me.lenZ);
            me.q   = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, me.q);
            me.qp  = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, me.q, X_plus( me.lenZ+(1:me.lenZ) ));
            me.qpp = mbeKinematicsSolver.accel_problem(me.bad_mech_phys_model, me.q, me.qp, X_plus( 2*me.lenZ+(1:me.lenZ) ));
       
        end % run_filter_iter
        
    end
end

