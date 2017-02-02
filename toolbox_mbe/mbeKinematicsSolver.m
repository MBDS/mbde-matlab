classdef mbeKinematicsSolver
    % Solves standard position, velocity and acceleration kinematic
    % problems, given the values for the independent coordinates of the
    % multibody model. Generic code, does not assume any particular number 
    % or order of independent variables. 
    % Works with objects of classes derived from mbeMechModelBase.
    
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

    properties(Constant)
        pos_problem_tol = 1e-13;      % Position problem: max |\Phi| value
        pos_problem_max_iters = 1000; % Position problem: max iterations
    end

    methods(Static,Access=public)
        function qpp = accel_problem(model, q, qp, zpp)
            % Solve for dependent accelerations in any mechanism
            % - model: an object of class derived from mbeMechModelBase
            % - q,dq: known pos & vel
            % - zpp: known indep accels            
            assert(size(q,2)==1,'q must be a column vector');
            assert(size(qp,2)==1,'qp must be a column vector');
            assert(size(zpp,2)==1,'zpp must be a column vector');
            
            phiq = model.jacob_phi_q(q); %was: jacob(q,MBK);
            A = phiq(:,model.get_dep_indxs);
            b = -phiq(:,model.get_indep_indxs)*zpp- model.jacob_phiqp_times_qp(q,qp); % was: phiqpqp(q,qp,MBK);
            qpp_dep = A\b;
            % Build complete accel. vector:
            qpp = [];
            qpp(model.get_dep_indxs,1) = qpp_dep;
            qpp(model.get_indep_indxs,1) = zpp;
        end % end accel_problem 

        function qp = vel_problem(model, q, zp)
            % Solve for dependent velocities in any mechanism
            % - model: an object of class derived from mbeMechModelBase
            % - q: known pos
            % - zp: known indep vels
            assert(size(q,2)==1,'q must be a column vector');
            assert(size(zp,2)==1,'zp must be a column vector');
            
            phiq = model.jacob_phi_q(q);
            A = phiq(:,model.get_dep_indxs);
            b = -phiq(:,model.get_indep_indxs)*zp;
            qp_dep = A\b;
            % Build complete vel vector:
            qp = [];
            qp(model.get_dep_indxs,1) = qp_dep;
            qp(model.get_indep_indxs,1) = zp;
        end % end vel_problem
        
        function [q_out, final_phi_norm,iters] = pos_problem(model, q)
            % Solve the iterative position problem for any mechanism
            % - model: an object of class derived from mbeMechModelBase
            % - q: known pos. Only those elements in the indices of independent 
            %      coordinates will remain unchanged.
            % - q_out: output, refined position
            assert(size(q,2)==1,'q must be a column vector');
            assert(size(q,1)==model.dep_coords_count(),'q length mismatch (|q|=%u, %u expected)',length(q),model.dep_coords_count());
            
            phi_0 = model.phi(q); % Was: phi(q,params);
            final_phi_norm = norm(phi_0);
            iters=0;
            didxs = model.get_dep_indxs();
            while final_phi_norm > mbeKinematicsSolver.pos_problem_tol && iters < mbeKinematicsSolver.pos_problem_max_iters
               phiq = model.jacob_phi_q(q); % Was: jacob (q,params);
               deltaq = phiq(:,didxs)\(-phi_0);
               q(didxs) = q(didxs) + deltaq;
               phi_0 = model.phi(q);
               final_phi_norm = norm(phi_0);
               iters = iters+1;
               if iters >= mbeKinematicsSolver.pos_problem_max_iters
                   warning('[mbeKinematicsSolver::pos_problem] Max number of iterations (%u) reached!',mbeKinematicsSolver.pos_problem_max_iters);
               end       
            end
            q_out = q;
        end  % end pos_problem()

        function [R] = calc_R_matrix(model,q)
            % Builds the R matrix of the given mechanism, for the given
            % dynamical state.
            nIndep=length(model.get_indep_indxs());
            R=zeros(model.dep_coords_count,nIndep);
            for k=1:nIndep,
                z_one=zeros(nIndep,1);
                z_one(k)=1.0;
                R(:,k) = mbeKinematicsSolver.vel_problem(model,q, z_one);
            end
        end

        function [dotR] = calc_dotR_matrix(model,q,qp,R)
            % Builds the $\dot{R}$ matrix, as the product of the hypermatrix dR_dq times qp.
            % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
            Rq = model.jacob_Rq(q,R);
            dotR= zeros(size(R,1),size(R,2));
            for k=1:model.dep_coords_count,
                dotR=dotR + Rq(:,:,k)*qp(k);
            end
        end
        
        
    end
end

