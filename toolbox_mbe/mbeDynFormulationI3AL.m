classdef mbeDynFormulationI3AL < mbeDynFormulationBase
    % Index 3 augmented lagrangian with mass orthogonal projections
    % TODO: Think about overriding integrator methods in this class, 
    % if it makes sense because better projections can be obtained from
    % I3AL.
   
    % I3AL parameters:
    properties (Access=public)
        alpha = 1e9; 
        tol_dyn = 1e-10; % dynamic tolerance (TODO: Explain!)
        iter_max = 100; % max number of iterations
    end
    
    properties (Access=private)
        % Lagrangian multipliers for I3AL: they are not parameters, but variables, but they
        % are here for convenience.
        lambda = []; %was: LagMul, zeros(4,1);
        % Cached qpp from last step
        prev_qpp = []; 
    end
    
    % --- Impl of virtual methods ---
    methods(Access=public)
        
        % Ctor:
        function [me] = mbeDynFormulationI3AL()
             % This method works integrating the entire "q", not "z".
            me.integrate_indep_coords = 0;
            me.integrator_type = mbeIntegratorTypes.intI3AL;
        end
        
        % See base class docs
        function [qpp,zpp,q_out,qp_out] = solve_for_accelerations(self,model,q,qp,context)
            dt = context.dt;
            % first step?
            if (length(self.prev_qpp)~=length(q))
                self.prev_qpp=zeros(length(q),1);
                if self.integrator_type ~=mbeIntegratorTypes.intI3AL
                    error('I3AL method error: this method cannot be run with any integrator different from mbeIntegratorTypes.intI3AL')
                end
            end
            M=model.M;
            
            %
            qp_g = -(2/dt*q + qp);
            qpp_g = -(4/dt^2*q + 4/dt*qp + self.prev_qpp);

            q = q + dt*qp+0.5*dt^2*self.prev_qpp;
            qp = (2/dt)*q+qp_g;
            qpp = (4/dt^2)*q+qpp_g;
            err = 1;
            iter=0;
            while err > self.tol_dyn && iter < self.iter_max
                iter = iter+1;

                phiq = model.jacob_phi_q(q);

                % first step? resize lambda
                nConstr=size(phiq,1);
                if (length(self.lambda)~=nConstr)
                    self.lambda=zeros(nConstr,1);
                end
                
                
                phi_0 = model.phi(q);
                Q = model.eval_forces(q,qp);
                f = dt^2/4*(M*qpp+phiq'*self.alpha*phi_0+phiq'*self.lambda-Q);
                [K,C]=model.eval_KC();
                f_q = M + 0.5*dt*C+0.25*dt^2*(phiq'*self.alpha*phiq+K);

                deltaq = f_q\-f;

                q = q+deltaq;
                qp = (2/dt)*q+qp_g;
                qpp = (4/dt^2)*q+qpp_g;
                phi_0 = model.phi(q);
                self.lambda = self.lambda + self.alpha*phi_0;
                err = norm(deltaq);
                if iter == self.iter_max
                    warning ('[mbeDynFormulationI3AL] Reached maximum number of iterations!')
                end
            end

            q_out = q;
            %Velocity and acceleration projections (with no time-dependent terms)
            qp_out = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qp);
            phiqpqp_0 = model.jacob_phiqp_times_qp(q_out, qp_out);  %Was: phiqpqp(q_out, qp_out, PARAMS);
            qpp = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qpp - 0.25*dt^2*phiq'*self.alpha*phiqpqp_0);
            self.prev_qpp = qpp;

            % Extract indep coords:
            zpp = qpp(model.get_indep_indxs());
        end
    end
    
end

