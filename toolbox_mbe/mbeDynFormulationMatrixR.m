classdef mbeDynFormulationMatrixR < mbeDynFormulationBase
    % Dynamic simulator of multibody systems based on the "Matrix R" method
   
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

    % (No parameters)
    %properties (Access=public)
    %end
    
    % --- Impl of virtual methods ---
    methods(Access=public)
        % See base class docs
        function [qpp,zpp,q,qp] = solve_for_accelerations(self,model,q,qp, context)
            nIndep = length(model.get_indep_indxs());

            R = mbeKinematicsSolver.calc_R_matrix(model,q); %qp = R*zp_next;
            M_bar = R'*model.M*R;
            Rpzp=mbeKinematicsSolver.accel_problem(model,q,qp, zeros(nIndep,1) );  %ace(q_next,qp_next,PARAMS,0);
            Q=model.eval_forces(q,qp);
            Q_bar = R'*(Q-model.M*Rpzp);
            zpp = M_bar\Q_bar;
            qpp = mbeKinematicsSolver.accel_problem(model,q,qp,zpp);
        end
    end
    
end

