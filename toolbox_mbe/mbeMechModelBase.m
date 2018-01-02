classdef (Abstract) mbeMechModelBase
    % Virtual base for generic, multibody models. 
    % A derived class represents a physical mechanism and a set of sensors.
        
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

    % Read-only, constant properties of the model
    properties(Abstract,Constant,GetAccess=public)
        % Dependent coordinates (q) count
        dep_coords_count;
        
        % Vector with the indices of the independent coordinates $z \in q$
        indep_idxs;
    
    end
    properties(Abstract,Access=public)
        % List of installed sensors (cell of objects derived from mbeSensorBase)
        installed_sensors;
        
        % Coordinates of fixed points (must be made available for the usage
        % of sensors, etc.)
        fixed_points;
        
        % Initial velocity for independent coords
        zp_init; 
        
        % Global mass matrix
        M;
    end    
    % Read-only properties of the model
    properties(Abstract,GetAccess=public,SetAccess=protected)
        % Initial, approximate position (dep coords) vector
        q_init_approx;
        
        
    end
    
    % Virtual methods to be implemented in all derived classes:
    methods(Abstract,Access=public)
        % Computes the vector of constraints $\Phi(q)$
        val = phi(me,q);
        % Computes the Jacobian $\Phi_q$
        phiq = jacob_phi_q(me,q);
        % Computes the Jacobian $\Phi_{qq}$ (it is a hypermatrix)
        phiqq = eval_phi_q_q(me,q);
        % Computes the product $\Phi_{qq} \dot{q}$
%         phiqq_qp = eval_phiqq_qp(me,q, qp)
        % Computes the time derivative of the Jacobian $\dot{Phi_q}$
        phiqp = eval_phiqp(me,q,qp);
        % Computes the Jacobian $\dot{\Phi_q} \dot{q}$
        phiqpqp = jacob_phiqp_times_qp(me,q,qp);
        % Computes hypermatrix  $\dot{\Phi_q}_q$
        phiqp_q = eval_phiqp_q(me,q,qp);
        % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
        Rq = jacob_Rq(me,q,R);
        % Evaluates the instantaneous forces
        Q = eval_forces(me,q,qp);
        % Evaluates potential energy of the whole system:
        V = eval_potential_energy(me,q);
        % Evaluates the stiffness & damping matrices of the system:
        [K, C] = eval_KC(me,q,dq);
        
        
        
        % Returns a *copy* of "me" after applying the given model errors (of class mbeModelErrorDef)
        [bad_model] = applyErrors(me, error_def); 
        % Plot a simplified version of the mechanism in the current axis.
        % - q: State 
        % - color_code: lines color ('r','b',...)
        % - do_fit: 1/0=whether to call {x,y}lim() 
        [] = plot_model_skeleton(me, q, color_code, do_fit);
        
    end  % end virtual methods
    
    methods(Access=public)
        % Nominal number of degrees of freedom (used in estimators, may be also in
        % dynamic formulations). Returns |z| <= |q| 
        function [n] = dof_count(me)
            n=length(me.indep_idxs);
        end
        
        % Get indices of $z's \in q$
        function [idxs] = get_indep_indxs(me)
            idxs = me.indep_idxs;
        end
        % Get indices of all z not in q
        function [idxs] = get_dep_indxs(me)
            idxs = 1:me.dep_coords_count; % [1,...,n]
            idxs(me.indep_idxs)=[]; % Remove those ones
        end
        % Partial derivative of dependent velocity qp and acceleration qpp
        % with respect to independent coordinates z and velocities zpp
        function [dqp_dz, dqpp_dz, dqpp_dzp] = eval_dDepVelAcel_dzzp(me,q,qp,qpp)
            idxs = me.get_indep_indxs();
            ndof = length(idxs);
            phiq = jacob_phi_q(me,q);
            phiqq = eval_phi_q_q(me,q);
            phiqq_qp = mbeUtils.mat3vec(phiqq,qp);
            phiqp = eval_phiqp(me,q,qp);
            [nrow,ncol] = size(phiq);
            B = zeros(ncol-nrow,ncol);
            for i = 1:ndof
                B(i,idxs(i))=1;
            end
            S = [phiq;B]\[eye(nrow); zeros(ncol-nrow,nrow)];
            R = mbeKinematicsSolver.calc_R_matrix(me,q);
            c_q = - mbeUtils.mat3vec(me.eval_phiqp_q(q,qp),qp); % NOTE: if there are time-dependent constraints, the term -\dot{\Phi_t}_q should be added
            c_qp = -phiqq_qp-phiqp; % NOTE: if there are time-dependent constraints, the term -\Phi_{tq} should be added
            Sc = mbeKinematicsSolver.accel_problem(me,q,qp,zeros(ndof,1));
            Sc_q = S*(-mbeUtils.mat3vec(phiqq, Sc)+c_q);
            Rz = zeros(ncol,ndof, ndof);
            for i = 1:ndof
                for j = 1:ncol
                    Rz(:,:,i) = Rz(:,:,i) - S*phiqq(:,:,j)*R*R(j,i); 
                end
            end
            zpp = qpp(idxs);
            
            dqp_dz   = -S*phiqp*R;
            dqpp_dz  = mbeUtils.mat3vec(Rz,zpp) + Sc_q*R + S*c_qp*dqp_dz;
            dqpp_dzp = S*c_qp*R;
            
        end
        
        % Partial derivative of independent accelerations wrt independent
        % velocities
        function [zpp_z,zpp_zp] = eval_dzpp_dposvel(me,q,qp)
            % Provides the partial derivative of the independent
            % acceleration with respect to the independent positions
            % (zpp_z) and the partial derivative of accelerations with
            % respect to independent velocities (zpp_zp)
            % For further details, see the paper
            % D. Dopico et al., Direct and Adjoint 
            % Sensitivity Analysis of Ordinary Differential Equation 
            % Multibody Formulations, Journal of Computational and 
            % Nonlinear Dynamics, vol 10, 2014 
            R=mbeKinematicsSolver.calc_R_matrix(me,q);
            [K,C] = eval_KC(me,q,qp);
            phiq = jacob_phi_q(me,q);
            phiqq = eval_phi_q_q(me,q);
            phiqq_qp = mbeUtils.mat3vec(phiqq,qp);
            phiqp = eval_phiqp(me,q,qp);
            [nrow,ncol] = size(phiq);
            B = zeros(ncol-nrow,ncol);
            idxs = me.get_indep_indxs();
            ndof = length(idxs);
            for i = 1:ndof
                B(i,idxs(i))=1;
            end
            S = [phiq;B]\[eye(nrow); zeros(ncol-nrow,nrow)];
            c_qp = -phiqq_qp-phiqp; % NOTE: if there are time-dependent constraints, the term -\Phi_{tq} should be added
            Qbar_qp = -R'*(C + me.M *S*c_qp);
            C_bar = -Qbar_qp*R;
            M_bar = R'*me.M*R;
            Sc = mbeKinematicsSolver.accel_problem(me,q,qp,zeros(ndof,1));
            Q = me.eval_forces(q,qp);
            Q_bar = R'*(Q-me.M*Sc);
            invMbar = inv(M_bar);
            zpp_zp = -invMbar*C_bar;     
            Rz = zeros(ncol,ndof, ndof);
            for i = 1:ndof
                for j = 1:ncol
                    Rz(:,:,i) = Rz(:,:,i) - S*phiqq(:,:,j)*R*R(j,i); 
                end
            end
            c_q = - mbeUtils.mat3vec(me.eval_phiqp_q(q,qp),qp); % NOTE: if there are time-dependent constraints, the term -\dot{\Phi_t}_q should be added
            Sc_q = S*(-mbeUtils.mat3vec(phiqq, Sc)+c_q);
            Rqt = -mbeUtils.mat2mat3(R',mbeUtils.mat3mat2(permute(phiqq,[2,1,3]),S'));
            Qbar_q = mbeUtils.mat3vec(Rqt,Q-me.M*Sc)-R'*(K+me.M*Sc_q);
            K_bar = -(Qbar_q-Qbar_qp*S*phiqp)*R; 
            zpp_z = -invMbar* mbeUtils.mat3vec((mbeUtils.mat3mat2(permute(Rz,[2,1,3]),me.M*R)+mbeUtils.mat2mat3(R'*me.M, Rz)),M_bar\Q_bar)...
                -M_bar\K_bar;
        end
        
%         % Partial derivative of independent accelerations wrt independent
%         % coordinates
%         function zpp_z = eval_dzpp_dz(me,q,qp)
%            zpp_z = 0* eval_dzpp_dzp(me,q,qp);
%         end

        % --- Start of sensors API ----

        function [obs] = sensors_simulate(me,q,qp,qpp)
            % Returns an instantaneous snapshop of what would be the sensors
            % readings for the given dynamic state (without noise). 
            % obs: A column with one reading per row. Row count is the number
            % of sensors available.
            
            nSensors = length(me.installed_sensors);
            obs=zeros(nSensors,1);
            for i=1:nSensors
                obs(i) = me.installed_sensors{i}.sensor_simulate(me,q,qp,qpp);
            end             
        end
        
        function [sensors_stds] = sensors_std_noise(me)
            % Returns the standard deviation of each sensor noise model.
            
            nSensors = length(me.installed_sensors);
            sensors_stds = zeros(nSensors,1);
            for i=1:nSensors
                sensors_stds(i) = me.installed_sensors{i}.sensor_std_noise();
            end             
        end
        
        function [dh_dq , dh_dqp, dh_dqpp] = sensors_jacob(me,q,qp,qpp)
            % Evaluates the sensor model jacobian for a given dyn state:
            %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
            % H= [  ...   |   ...   |    ...   ] = [ dh_dq , dh_dqp, dh_dqpp ]
            %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
            
            n = length(q); nSensors = length(me.installed_sensors);
            dh_dq  = zeros(nSensors,n);
            dh_dqp = zeros(nSensors,n);
            dh_dqpp= zeros(nSensors,n);
            for i=1:nSensors
                [dh_dq(i,:), dh_dqp(i,:), dh_dqpp(i,:)] = me.installed_sensors{i}.sensor_jacob(me,q,qp,qpp);
            end             
        end
        
        function [dh_dz , dh_dzp, dh_dzpp] = sensors_jacob_indep(me,q,qp,qpp)
            % Computes the Jacobian of sensors wrt independent coordinates,
            % by evaluating sensors_jacob() and applying the derivative
            % chain rule (below, all are partial derivatives):
            %
            %  Note: h = h(q,qp,qpp)  [See papers, or internal technical report "documentation_ekf" for details]
            %
            %  dh     dh    dq      dh    dqp     dh          dh                   dh
            % ---- = ----  ----  + ----- ----- = ---- R(q) + ----- R_q * R * zp + ---- (...)
            %  dz     dq    dz      dqp    dz     dq          dqp                 dqpp   
            %
            %  dh      dh      d qp      dh          dh         
            % ----- = -----  -------- = ----- * R + ----- * ....
            %  dzp     dqp     dzp       dqp         dq
            %
            %  dh      dh     d qpp      dh                 
            % ----- = -----  -------- = ----- * R
            %  dzpp   dqpp     dzpp      dqpp    
            % 
            
            % dh_dz
            R=mbeKinematicsSolver.calc_R_matrix(me,q);
%             Rq = me.jacob_Rq(q,R);
            [dh_dq , dh_dqp, dh_dqpp] = me.sensors_jacob(q,qp,qpp);
            [dqp_dz, dqpp_dz, dqpp_dzp] = me.eval_dDepVelAcel_dzzp(q,qp,qpp);
%             iindex = me.get_indep_indxs;
%             dh_dz = zeros(length(me.installed_sensors), length(iindex));
%             for i = 1:length(iindex)
%                 dh_dz(:,i) = dh_dq*R(:,i) + dh_dqp*Rq(:,:,i)*R(:,i)*qp(iindex(i)) +dh_dqpp * Rq(:,:,i)*R(:,i)*qpp(iindex(i));
%             end
            dh_dz = dh_dq*R+dh_dqp*dqp_dz+dh_dqpp*dqpp_dz;
			% TODO: Complete partial derivatives!
            
            % dh_dzp: TODO: Complete partial derivatives!
            dh_dzp = dh_dqp *  R + dh_dqpp*dqpp_dzp;
            
            % dh_dzpp:
            dh_dzpp = dh_dqpp *  R;       
        end
        % --- End of sensors API ----
    end
    
    % Static, public, auxiliary methods:
    methods(Static,Access=public)
                
        function [M] = massMatrixBar(m,L)
            % Aux method: builds the 4x4 mass matrix of a planar rod, with
            % mass "m", length between the 2 points "L" and center of
            % gravity cog=[cx cy]
            cog=[0.5*L, 0];
            I=(1/3)*m*(L^2); 
            a=I/(L^2); 
            bx=m*cog(1)/L;
            by=m*cog(2)/L;
            
            M=[m+a-2*bx,        0,  bx-a,   -by; ...
                      0, m+a-2*bx,    by, bx-a; ...
                   bx-a,       by,     a,    0;...
                    -by,     bx-a,     0,    a ];
        end 
        
        function [M] = massMatrixDisc(m,L,radius,cog)
            % Aux method: builds the 4x4 mass matrix of a planar disc, with
            % mass "m", radius between the 2 points "radius" and center of
            % gravity cog=[cx cy] in local coordinates
            I=(1/2)*m*(radius^2); 
            a=I/(L^2); 
            bx=m*cog(1)/L;
            by=m*cog(2)/L;
            
            M=[m+a-2*bx,        0,  bx-a,   -by; ...
                      0, m+a-2*bx,    by, bx-a; ...
                   bx-a,       by,     a,    0;...
                    -by,     bx-a,     0,    a ];
        end
        
        
        

    end % end static methods
    
    
    
end

