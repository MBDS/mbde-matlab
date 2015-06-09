classdef mbeEstimatorCEKF < mbeEstimatorFilterBase
    % Continuous-time Extended Kalman filter (EKF), trapezoidal integrator
    
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
        % (Is this worth the computational cost? Probably not) 
        % Reflect the new X_less into q,qp,qpp for accurately simulate the
        % sensor & its linearization point:
        eval_sensors_at_X_less = 0; % 1
        
        tol_dyn = 1e-10;
        iter_max = 100;
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
            
            % Was: CovPlantNoise = diag([0.001 0.01].^2);
            me.CovPlantNoise = diag([...
                ones(1,me.lenZ)*me.transitionNoise_Zp, ...   % cont vs discr: "*(me.dt)^2"
                ones(1,me.lenZ)*me.transitionNoise_Zpp]);
            
            % Sensors noise model:
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();
            me.CovMeasurementNoise = diag(sensors_stds.^2);
        end
        
        % Run one timestep of the estimator (see docs in mbeEstimatorFilterBase)
        function [] = run_filter_iter(me, obs)
            % Aux indices for indep coords in X & X_next
            ind_idxs1 = (1:length(me.iidxs));
            ind_idxs2 = length(me.iidxs) + ind_idxs1;
            
            zpp = me.qpp(me.iidxs);

            % Kalman state:     
            X  = [me.q(me.iidxs) ; me.qp(me.iidxs)];
            Xp = [me.qp(me.iidxs); me.qpp(me.iidxs)];
            
            R = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model,me.q);
            R_q = me.bad_mech_phys_model.jacob_Rq(me.q,R);
            Rp = R; % Allocation (this should be faster than Rp = zeros(size(R))
            for i=1:length(me.iidxs)
                Rp(:,i) = R_q(:,:,i)*me.qp;  % \dot{R}
            end
            
                        
            M = me.bad_mech_phys_model.M;
            
            M_bar = R'*M*R;
            M_Rq_R_zpp = R; % Allocation (this should be faster than Rp = zeros(size(R))
            for i = 1:length(me.iidxs)
                M_Rq_R_zpp(:,i) = M*R_q(:,:,i)*R(:,i)*zpp(i);
            end
            F_21 = -(M_bar)\(R'*(2*M_Rq_R_zpp)); % was F_21 = -(M_bar)\(R'*(2*M*R_q*R*zpp));
            [~, C] = me.bad_mech_phys_model.eval_KC(me.q,me.qp);
            F_22 = -(M_bar)\(R'*(M*Rp + C*R));
            F=[ zeros(me.lenZ)  ,   eye(me.lenZ);
               F_21, F_22 ];
            FP = F*me.P;            
            
            Pp=FP+FP'+me.CovPlantNoise;
            if (~isempty(obs))  % Only reduce uncertainty if we actually have an observation
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp] = me.bad_mech_phys_model.sensors_jacob_indep(me.q,me.qp,me.qpp);
                H=[dh_dz, dh_dzp];   % Was: H = [dh_dq(:,me.iidxs) , dh_dqp(:,me.iidxs) ];
                Pp = Pp - ((me.P*H')/(me.CovMeasurementNoise))*H*me.P;
            end            
            P_next = me.P+Pp*me.dt;
            % Force P_next to be sym (fix tiny roundoff errors)
            P_next = 0.5*(P_next+P_next');
                        
            % Prediction 
            X_next = X + Xp*me.dt;
            Xp_next = Xp;
            q_next = me.q;
            q_next(me.iidxs) = X_next(ind_idxs1);
            q_next  = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, q_next);
            
            
            % R matrix calculation
            R_next = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model, q_next);
            % Velocity problem using R_next
            qp_next = R_next *  X_next(ind_idxs2);
            % Acceleration
            Rpzp = mbeKinematicsSolver.accel_problem( ...  %ace(q_next,qp_next, params,0);
                me.bad_mech_phys_model, ...
                q_next, qp_next, zeros(length(me.iidxs),1) );
            
            M_bar = R_next'*M*R_next;
            
            Q = me.bad_mech_phys_model.eval_forces(q_next,qp_next);
            Q_bar = R_next'*(Q-M*Rpzp);
            
            Xp_next(ind_idxs2) = M_bar\Q_bar;
            
            
            
            % Virtual sensors
            if (~isempty(obs))
                % Eval observation Jacobian wrt the indep. coordinates:
                [dh_dz , dh_dzp] = me.bad_mech_phys_model.sensors_jacob_indep(q_next,qp_next,me.qpp);
                H=[dh_dz, dh_dzp];   % Was: H = [dh_dq(:,me.iidxs) , dh_dqp(:,me.iidxs) ];
				% Kalman gain
				K = (P_next*H')/(me.CovMeasurementNoise);
				
                % measurement update:
                obs_predict = me.bad_mech_phys_model.sensors_simulate(q_next,qp_next,me.qpp); % TODO: acceleration problem should be solved to get qpp_next instaead of me.qpp
                
                Innovation = (obs_predict-obs);
                % Residual
                g1 = Xp_next(ind_idxs1)-X_next(ind_idxs2)+K(ind_idxs1,:)*Innovation;
                g2 = M_bar*Xp_next(ind_idxs2)-Q_bar+M_bar*K(ind_idxs2,:)*Innovation;
                g=[g1;g2];
                error = norm(g); 
                iter = 0; 
            
                while error > me.tol_dyn && iter < me.iter_max
                    iter = iter+1;
                    %Rp = ace(q_next, R_next, l, x, 0); % Derivative of R_next with time 
                    Rp = R;
                    for i = 1:length(me.iidxs)
                        Rp(:,i) = R_q(:,:,i)*qp_next;
                    end
                    % Eval observation Jacobian wrt the indep. coordinates:
                    [dh_dz , dh_dzp] = me.bad_mech_phys_model.sensors_jacob_indep(q_next,qp_next,me.qpp); % TODO: acceleration problem should be solved to get qpp_next instaead of me.qpp
                    H=[dh_dz, dh_dzp];   % Was: H = [dh_dq(:,me.iidxs) , dh_dqp(:,me.iidxs) ];
%                     Rp = R_q*me.qp;
                    eye_matrix = eye(length(me.iidxs));
                    g_q_1=[2/me.dt*eye_matrix, -1*eye_matrix;
                            0*eye_matrix, R_next'*(M*Rp)+2/me.dt*M_bar];
                    g_q_2=[K(ind_idxs1,:)*H(:,ind_idxs1) K(ind_idxs1,:)*H(:,ind_idxs2);
                            M_bar*K(ind_idxs2,:)*H(:,ind_idxs1) M_bar*K(ind_idxs2,:)*H(:,ind_idxs2)];
                        
                    g_q = g_q_1+ g_q_2;  % Dynamical system jacobian        
                    deltaX = -g_q\g;     
                    X_next = X_next + deltaX; 
                    Xp_next = 2/me.dt*(X_next-X)-Xp; % Integrator equation (Trap. rule)
            
                    q_next(me.iidxs) = X_next(ind_idxs1);
                    q_next  = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, q_next);

                    % R matrix calculation
                    R_next = mbeKinematicsSolver.calc_R_matrix(me.bad_mech_phys_model, q_next);
                                        
                    M_bar = R_next'*M*R_next;
                    %Q = ...  (doesn't change)
                    % Acceleration
                    Rpzp = mbeKinematicsSolver.accel_problem( ...  %ace(q_next,qp_next, params,0);
                        me.bad_mech_phys_model, ...
                        q_next, qp_next, zeros(length(me.iidxs),1) );
                    Q_bar = R_next'*(Q-M*Rpzp);
                    qp_next = R_next*X_next(ind_idxs2);
                    %             zpp = M_bar\Q_bar;
                    %             qpp_next = ace(q_next,qp_next, l,x,zpp);
                    %qpp_next = ace(q_next,qp_next, params, Xp_next(2)); % unused??
                    
                    %y = virtual_sensor_eval (q_next,qp_next,params);
                    %Innovation = (y-y_sensor);
                    obs_predict = me.bad_mech_phys_model.sensors_simulate(q_next,qp_next,me.qpp); % TODO: acceleration problem should be solved to get qpp_next instaead of me.qpp
                    Innovation = (obs_predict-obs);                    
                    
                    g1 = Xp_next(ind_idxs1)-X_next(ind_idxs2)+K(ind_idxs1,:)*Innovation;
                    g2 = M_bar*Xp_next(ind_idxs2)-Q_bar+M_bar*K(ind_idxs2,:)*Innovation;
                    g=[g1;g2];
                    error = norm(g);
                end
            end   % End if we have an observation             
			
% 			zpp = Xp_next(ind_idxs2);
            zpp = M_bar\Q_bar;
            
            % Recover KF -> MBS coordinates
            % ------------------------------
            me.P = P_next;
            
            me.q   = mbeKinematicsSolver.pos_problem(me.bad_mech_phys_model, q_next);
            me.qp  = mbeKinematicsSolver.vel_problem(me.bad_mech_phys_model, q_next, X_next(ind_idxs2) );
            me.qpp = mbeKinematicsSolver.accel_problem( ...  %qpp_next = ace(q_next,qp_next, params, zpp);
                me.bad_mech_phys_model, ...
                q_next, qp_next, zpp ...
                );
            
        end % run_filter_iter
        
    end
end

