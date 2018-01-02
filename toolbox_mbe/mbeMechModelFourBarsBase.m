classdef mbeMechModelFourBarsBase < mbeMechModelBase
    % Mechanism physical model: Generic 4 bars linkage (Abstract base
    % class). Derived classes define particular examples of mechanisms with
    % especific masses, lengths, etc.
    % Modeled in Natural coordinates plus one relative angle coordinate 
    % at the left-hand side fixed end.
    %
    %             2:(q3,q4)
    %               +---------o
    %               |       (xb,yb)
    %               | 
    %               | 
    %     o---------+ 1:(q1,q2)
    %    (xa,ya)       
    %
    %  - q5: Angle (xa,ya)-(q1,q2)
    %
    %
    
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

    % (Abstract) Read-only, constant properties of the model
    properties(Constant,GetAccess=public)
        % Dependent coordinates count
        dep_coords_count = 5;
        
        % A vector with the indices of the independent coordinates in "q":
        indep_idxs = 5;
    end
    % (Abstract) Read-only properties of the model
    properties(GetAccess=public,SetAccess=protected)
        % Initial, approximate position (dep coords) vector
        q_init_approx=zeros(mbeMechModelFourBarsBase.dep_coords_count,1);
        % Generalized vector of constant forces.
        Qconst=zeros(mbeMechModelFourBarsBase.dep_coords_count,1);
        
    end
    
    % Model-specific properties:
    properties(Access=public)
        % Initial velocity for independent coords
        zp_init=[0];
        % Global mass matrix
        M = [0];
        % Fixed point coords:
        xA,yA,xB,yB, fixed_points;
        
        % Gravity:
        g = -9.81;
        
        % lengths:
        bar_lengths; 
        
        % Masses:
        mA1,m12,m2B;
        
        % Force vector (gravity forces only):
        Qg;
        % Force vector: generalized forces vector calculated from estimated
        % input forces
        Qm = zeros(mbeMechModelFourBarsBase.dep_coords_count,1); 
        
        % damping coefficient (TODO: At which joint??)
        C = 0;
    end
    
    methods 
        % Constructor: must be implemented in derived classes to fill-in
        % all mechanical parameters.
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Computes the vector of constraints $\Phi(q)$
        function val = phi(me,q)
            x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4); theta = q(5);
            LA1 = me.bar_lengths(1); L12 = me.bar_lengths(2); L2B = me.bar_lengths(3);
            val = [(me.xA-x1)^2 + (me.yA-y1)^2 - LA1^2;
                        (x1-x2)^2 + (y1-y2)^2 - L12^2;
                        (me.xB-x2)^2 + (me.yB-y2)^2 - L2B^2;
                    mbe_iff(abs(sin(theta)) < 0.7,... 
                        y1-me.yA-LA1*sin(theta), ...
                        x1-me.xA-LA1*cos(theta) ...
                        ) ...
                    ];
        end % of phi()
        
        % Computes the Jacobian $\Phi_q$
        function phiq = jacob_phi_q(me,q)
            % (From old code in jacob.m)
            % q: coordinates
            % l: bar length vector
            % x: fixed points positions
            x1 = q(1); y1 = q(2); x2 = q(3); y2 = q(4); theta = q(5);
            LA1 = me.bar_lengths(1); 
            phiq = [...
                -2*(me.xA-x1), -2*(me.yA-y1),             0,             0,              0;
                    2*(x1-x2),     2*(y1-y2),    -2*(x1-x2),    -2*(y1-y2),              0;
                            0,             0, -2*(me.xB-x2), -2*(me.yB-y2),              0;
                    mbe_iff(abs(sin(theta)) < 0.7,... 
                            [0,             1,             0,             0, -LA1*cos(theta)], ...
                            [1,             0,             0,             0,  LA1*sin(theta)] ...
                        ) ...                            
                    ];
        end % jacob_phi_q()
        
        % Computes the Jacobian $\Phi_{qq} (it is a hypermatrix)$
        function phiqq = eval_phi_q_q(me,q)
            %             x1 = q(1); y1 = q(2); x2 = q(3); y2 = q(4);
            theta = q(5);
            LA1 = me.bar_lengths(1);
            phiqq = zeros(4,5,5);
            phiqq(:,:,1) = [...
                2, 0,  0,  0,  0;
                2,  0,  -2, 0,  0;
                0,  0,  0,  0,  0;
                0,  0,  0,  0,  0;...
                ];
            phiqq(:,:,2) =  [...
                0,  2, 0,  0,  0;
                0,  2,  0,  -2, 0;
                0,  0,  0,  0,  0;
                0,  0,  0,  0,  0;...
                ];
            phiqq(:,:,3) =  [...
                0,  0,  0,  0,  0;
                -2,  0,  2,  0,  0;
                0,  0,  2,  0,  0;
                0,  0,  0,  0,  0;...
                ];
            phiqq(:,:,4) =  [...
                0,  0,  0,  0,  0;
                0,  -2,  0,  2,  0;
                0,  0,  0,  2,  0;
                0,  0,  0,  0,  0;...
                ];
            phiqq(:,:,5) =  [...
                0,  0,  0,  0,  0;
                0,  0,  0,  0,  0;
                0,  0,  0,  0,  0;
                mbe_iff(abs(sin(theta)) < 0.7,...
                [0,             0,             0,             0, LA1*sin(theta)], ...
                [0,             0,             0,             0,  LA1*cos(theta)] ...
                ) ...
                ];
            
        end % eval_phi_q_q

        % Computes the time derivative of the Jacobian $\dot{\Phi_q}$
        function phiqp = eval_phiqp(me,q,qp)
            %x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4);
            theta = q(5);
            x1p = qp(1) ;y1p = qp(2);x2p = qp(3);y2p = qp(4);thetap = qp(5);
            LA1 = me.bar_lengths(1);
            
            phiqp = [...
                2*x1p,        2*y1p,            0,             0,              0;
                2*(x1p-x2p),  2*(y1p-y2p), -2*(x1p-x2p),  -2*(y1p-y2p),              0;
                0,             0,       2*x2p,         2*y2p,              0;
                mbe_iff(abs(sin(theta)) < 0.7,...
                [0,             0,           0,             0, LA1*sin(theta)*thetap ], ...
                [0,             0,           0,             0, LA1*cos(theta)*thetap ]  ...
                ) ...
                ];
        end % eval_phiqp
        
        % Computes hypermatrix  $\dot{\Phi_q}_q$
        function phiqp_q = eval_phiqp_q(me,q,qp)
            %x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4);
            theta = q(5);
            %x1p = qp(1) ;y1p = qp(2);x2p = qp(3);y2p = qp(4);
            thetap = qp(5);
            LA1 = me.bar_lengths(1);
            phiqp_q = zeros(4,5,5);
            phiqp_q(:,:,5) = [...
                0, 0, 0, 0, 0;...
                0, 0, 0, 0, 0;...
                0, 0, 0, 0, 0;...
                mbe_iff(abs(sin(theta)) < 0.7,...
                [0,             0,           0,             0, LA1*cos(theta)*thetap ], ...
                [0,             0,           0,             0, -LA1*sin(theta)*thetap ]  ...
                ) ...
                ];
        end % eval_phiqp_q
        
        % Computes the Jacobian $\dot{\Phi_q} \dot{q}$
        function phiqpqp = jacob_phiqp_times_qp(me,q,qp)
            dotphiq = eval_phiqp(me,q,qp);
            phiqpqp = dotphiq * qp;
        end % jacob_phiqp_times_qp
        
        % Computes the partial derivative of velocity constraints wrt q
        function Phiqqpq = Phiq_times_qp_q(me,q,qp)
%             x1 = q(1); y1 = q(2); x2 = q(3); y2 = q(4); 
            theta = q(5);
            LA1 = me.bar_lengths(1);
            xp1 = qp(1); yp1 = qp(2); xp2 = qp(3); yp2 = qp(4); thetap = qp(5);
            Phiqqpq = [2*xp1, 2*yp1, 0, 0, 0;
                       2*(xp1-xp2), 2*(yp1-yp2), 2*(xp2-xp1), 2*(yp2-yp1),0;
                       0, 0, 2*xp2, 2*yp2, 0;
                       mbe_iff(abs(sin(theta)) < 0.7,... 
                       [0,0,0,0, LA1*sin(theta)*thetap],...
                       [0,0,0,0, LA1*cos(theta)*thetap])];
        end % Phiq_times_qp_q
        
        % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
        function Rq = jacob_Rq(me,q,R)
            theta = q(5); 
            LA1 = me.bar_lengths(1);
            phiq = me.jacob_phi_q(q);
            dep_coord_ind = 1:5;
            dep_coord_ind(me.indep_idxs) = [];
            phiqd = phiq(:,dep_coord_ind); % Jacobian dependent part

            Phi_qqR = [2*R(1), 2*R(2), 0, 0, 0;
                      2*R(1)-2*R(3), 2*(R(2)-R(4)), 2*(-R(1)+R(3)), 2*(-R(2)+R(4)), 0;
                      0, 0, 2*R(3), 2*R(4), 0;
                    mbe_iff(abs(sin(theta)) < 0.7,... 
                      [0, 0, 0, 0, LA1*sin(theta)], ...
                      [0, 0, 0, 0, LA1*cos(theta)] ...
                      ) ...
                      ];
            Rd_q = phiqd\(-Phi_qqR);
            Rq = zeros(5,5);
            Rq(dep_coord_ind,:) = Rd_q;            
        end
        

        % Evaluates the instantaneous forces
        function Q = eval_forces(me,q,qp)
            Q_var = zeros(me.dep_coords_count,1);
            Q_var(5) = -me.C*qp(5);
            Q = me.Qg+Q_var+me.Qconst+me.Qm;
%             Q = zeros(me.dep_coords_count,1)+me.Qconst+me.Qm;
        end % eval_forces
        % sets the value of the constant forces to be applied during all the
        % simulation
        function [estim_wQcnst] = set_Qconst(me, Q_constant)
            estim_wQcnst = me; 
            estim_wQcnst.Qconst = Q_constant;
        end % end of set_Qconst

        % Evaluates the stiffness & damping matrices of the system:
        function [K, C] = eval_KC(me, q,dq)
            K = zeros(me.dep_coords_count,me.dep_coords_count);
            C = zeros(me.dep_coords_count,me.dep_coords_count);
            C(5,5) = me.C;
        end

        % Returns a copy of "me" after applying the given model errors (of
        % class mbeModelErrorDef)
        function [bad_model] = applyErrors(me, error_def)
            bad_model = me; 
             
            % Init with no error:
            ini_vel_error = 0;
            ini_pos_error = 0;
            grav_error = 0;
            damping_coef_error = 0;
            for i = 1:length(error_def.error_type)
                switch error_def.error_type(i)
                    case 0
    %                     ini_vel_error = 0;
    %                     ini_pos_error = 0;
    %                     grav_error = 0;
    %                     damping_coef_error = 0;
                    % 1: Gravity
                    case 1
                        grav_error = 1*error_def.error_scale;
                    % 2: Initial pos error
                    case 2
                        ini_pos_error = error_def.error_scale * pi/16;
                    % 3: Initial vel error
                    case 3
                        ini_vel_error = 10 * error_def.error_scale;
                    % 4: damping (C) param (=0)
                    case 4 
                        damping_coef_error = -1*me.C * error_def.error_scale;
                    % 5: damping (C) param (=10)
                    case 5
%                         ini_vel_error = 0;
%                         ini_pos_error = 0;
%                         grav_error = 0;
                        damping_coef_error = 10 * error_def.error_scale;
                    case 6 
                        bad_model.bar_lengths(1) = bad_model.bar_lengths(1)*(100+error_def.error_scale)/100; % Lenght error in bar 1 (error_scale are 1 % of lenght error)
                    case 7 
                        bad_model.mA1 = bad_model.mA1*(100+10*error_def.error_scale)/100;% Mass error in bar 1(error_scale usits are 10% of mass error)
                    case 8
                        bad_model = bad_model.set_Qconst(bad_model.Qconst*0);
                    otherwise
                        error('Unhandled value!');
                end
            end
            bad_model.g = bad_model.g+grav_error; % gravity error
            bad_model.zp_init = bad_model.zp_init+ini_vel_error; % initial velocity error
            bad_model.q_init_approx(5)=bad_model.q_init_approx(5)+ini_pos_error; %initial position error
            bad_model.C=bad_model.C+damping_coef_error; %damping coefficient error

            % Weight vector 
            % WARNING: This vector MUST be updated here, after modifying the mass and the "g"
            % vector !
            bad_model=bad_model.update_Qg();
            % Update mass matrix (It only need to be updated if the mass error error has
            % been applied)
            bad_model = bad_model.update_M();            
        end % applyErrors
        
        % See docs in base class
        function [] = plot_model_skeleton(me, q, color_code, do_fit)
            plot([me.fixed_points(1),q(1),q(3),me.fixed_points(3)], ...
                 [me.fixed_points(2),q(2),q(4),me.fixed_points(4)] ,color_code,...
                 'LineWidth',5, 'MarkerSize',35);
            if (do_fit)
                axis equal;
                xlim ([me.fixed_points(1)-1.2*me.bar_lengths(1),1.1*me.fixed_points(3)]);
                ylim ([me.fixed_points(2)-1.2*me.bar_lengths(1),me.fixed_points(2)+1.2*me.bar_lengths(3)]);
            end
        end
        
    end % methods
end % class
