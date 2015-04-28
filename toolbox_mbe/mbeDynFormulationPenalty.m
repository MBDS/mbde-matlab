classdef mbeDynFormulationPenalty < mbeDynFormulationBase
    % Dynamic simulator of multibody systems based on the penalty method
    % (based on old code 'MBS_solve')
   
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

