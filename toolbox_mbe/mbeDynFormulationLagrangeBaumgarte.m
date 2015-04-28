classdef mbeDynFormulationLagrangeBaumgarte < mbeDynFormulationBase
    % Dynamic simulator of multibody systems based on Lagrange multiplers, 
    % with Baumgarte stabilization.
   
    % Penalty parameters:
    properties (Access=public)
        baum_epsilon = 1;
        baum_omega   = 10;        
    end
    
    % --- Impl of virtual methods ---
    methods(Access=public)
        % See base class docs
        function [qpp,zpp,q,qp] = solve_for_accelerations(self,model,q,qp,context)
            Phi  = model.phi(q);
            Phiq = model.jacob_phi_q(q);

            % right hand side
            % Use Baumgarten Stabilization
            %  c= -\dot{Phi_q} * \dot{q}  - 2*eps*omega*dotPhi - omega^2 * Phi
            %c = dotPhi_q * dq;
            c = model.jacob_phiqp_times_qp(q,qp);
            dotPhi = Phiq*qp;
            cstab = -c - 2*self.baum_epsilon*self.baum_omega*dotPhi - (self.baum_omega.^2) * Phi;
            Q=model.eval_forces(q,qp);
            rhs = [Q; cstab];

            A = [model.M Phiq'; Phiq zeros(size(Phiq,1),size(Phiq,1))];
            sol = A\rhs;
            qpp = sol(1:length(q));
            % Extract indep coords:
            zpp = qpp(model.get_indep_indxs());
        end
    end
    
end

