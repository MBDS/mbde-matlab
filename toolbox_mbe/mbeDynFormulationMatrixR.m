classdef mbeDynFormulationMatrixR < mbeDynFormulationBase
    % Dynamic simulator of multibody systems based on the "Matrix R" method
   
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

