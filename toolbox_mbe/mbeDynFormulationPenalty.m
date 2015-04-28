classdef mbeDynFormulationPenalty < mbeDynFormulationBase
    % Dynamic simulator of multibody systems based on the penalty method
    % (based on old code 'MBS_solve')
   
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

    % Penalty parameters:
    properties (Access=public)
        alpha = 1e9; % penalty (it takes influence in conservation of energy)
        omega = 10;
        psi   = 1;
    end
    
    % --- Impl of virtual methods ---
    methods(Access=public)
        % See base class docs
        function [qpp,zpp,q,qp] = solve_for_accelerations(self,model,q,qp, context)
            jac = model.jacob_phi_q(q); % Was: jacob(q,params);
            A = model.M + self.alpha*((jac)')*jac;
            phiqpqp_0 = model.jacob_phiqp_times_qp(q,qp); % Was:  phiqpqp(q,qp,params);
            phi_0 = model.phi(q);
            phip = jac*qp; 
            Q = model.eval_forces(q,qp);
            b = Q-self.alpha*(jac)'*(phiqpqp_0 + 2*self.psi*self.omega*phip + self.omega^2*phi_0);
            qpp = A\b;
            zpp = qpp(model.get_indep_indxs());
        end
    end
    
end

