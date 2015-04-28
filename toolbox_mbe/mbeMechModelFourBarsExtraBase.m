classdef mbeMechModelFourBarsExtraBase < mbeMechModelBase
    % TODO: Extra points, extra coordinates ....
    % Mechanism physical model: Generic 4 bars linkage (Abstract base
    % class). Derived classes define particular examples of mechanisms with
    % especific masses, lengths, etc.
    % Modeled in Natural coordinates plus one relative angle coordinate
    % at the left-hand side fixed end.
    %
    %             2:(q3,q4)
    %               +---------o    - q6: Angle (xb,yb)-(q3,q4)
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
        dep_coords_count = 8;
        
        % A vector with the indices of the independent coordinates in "q":
        indep_idxs = 5; % dof->theta
        %indep_idxs = 8; % dof->beta
    end
    % (Abstract) Read-only properties of the model
    properties(GetAccess=public,SetAccess=protected)
        % Initial, approximate position (dep coords) vector
        q_init_aprox=zeros(mbeMechModelFourBarsBase.dep_coords_count,1);
        
        % Initial velocity for independent coords
        zp_init=[0];
    end
    
    % Model-specific properties:
    properties(Access=public)
        % Global mass matrix
        M;
        
        % Fixed point coords:
        xA,yA,xB,yB, fixed_points;
        
        % Gravity:
        g = -10;
        
        % lengths:
        bar_lengths;
        
        % Masses:
        mA1,m12,m2B;
        
        % Force vector (gravity forces only):
        Qg;
        
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
            x3 = q(6); y3 = q(7); beta = q(8);
            LA1 = me.bar_lengths(1); L12 = me.bar_lengths(2); L2B = me.bar_lengths(3);
            L13 = me.bar_lengths(4); L32 = me.bar_lengths(5);
            val=zeros(length(q)-1,1);
            val(1) = (x1-me.xA)^2 + (y1-me.yA)^2 - LA1^2;
            val(2) = (x2-x1)^2 + (y2-y1)^2 - L12^2;
            val(3) = (x2-me.xB)^2 + (y2-me.yB)^2 - L2B^2;
            val(4) = mbe_iff(abs(sin(theta)) < 0.7,...
                y1-me.yA-LA1*sin(theta), ...
                x1-me.xA-LA1*cos(theta));
            val(5) = (x3-x1)^2+(y3-y1)^2-L13^2;
            val(6) = (x3-x2)^2+(y3-y2)^2-L32^2;
            val(7) = (x2-me.xB)-L2B*cos(beta);
        end % of phi()
        
        % Computes the Jacobian $\Phi_q$
        function phiq = jacob_phi_q(me,q)
            % (From old code in jacob.m)
            % q: coordinates
            % l: bar length vector
            % x: fixed points positions
            x1 = q(1); y1 = q(2); x2 = q(3); y2 = q(4); theta = q(5);
            x3 = q(6); y3 = q(7); beta = q(8);
            LA1 = me.bar_lengths(1); L2B = me.bar_lengths(3);
            phiq = zeros(7,8);
            phiq(1,1) = 2*(x1-me.xA); phiq(1,2) = 2*(y1-me.yA);
            phiq(2,1) =  -2*(x2-x1); phiq(2,2) =  -2*(y2-y1); phiq(2,3) =  2*(x2-x1); phiq(2,4) = 2*(y2-y1);
            phiq(3,3) = 2*(x2-me.xB); phiq(3,4) = 2*(y2-me.yB);
            phiq(4,:) = mbe_iff(abs(sin(theta)) < 0.7,...
                [0,             1,             0,             0, -LA1*cos(theta),           0,           0,     0], ...
                [1,             0,             0,             0,  LA1*sin(theta),           0,           0,     0]);
            phiq(5,1) = -2*(x3-x1); phiq(5,2) = -2*(y3-y1); phiq(5,6) = 2*(x3-x1); phiq(5,7) = 2*(y3-y1);
            phiq(6,3) = -2*(x3-x2); phiq(6,4) = -2*(y3-y2); phiq(6,6) = 2*(x3-x2); phiq(6,7) = 2*(y3-y2);
            phiq(7,3) = 1; phiq(7,8) = L2B*sin(beta);
            
        end % jacob_phi_q()
        
        % Computes the Jacobian $\dot{\Phi_q} \dot{q}$
        function phiqpqp = jacob_phiqp_times_qp(me,q,qp)
            %x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4);
            theta = q(5);  thetap = qp(5);
            beta  = q(8);  betap = qp(8);
            x1p   = qp(1); y1p = qp(2);
            x2p   = qp(3); y2p = qp(4);
            x3p   = qp(6); y3p = qp(7);
            LA1   = me.bar_lengths(1); L2B = me.bar_lengths(3);
            dotphiq = zeros(7,8);
            dotphiq(1,1) = 2*(x1p); dotphiq(1,2) = 2*(y1p);
            dotphiq(2,1) =  -2*(x2p-x1p); dotphiq(2,2) =  -2*(y2p-y1p); dotphiq(2,3) =  2*(x2p-x1p); dotphiq(2,4) = 2*(y2p-y1p);
            dotphiq(3,3) = 2*(x2p); dotphiq(3,4) = 2*(y2p);
            dotphiq(4,:) = mbe_iff(abs(sin(theta)) < 0.7,...
                [0,             0,             0,             0,  LA1*sin(theta)*thetap,           0,           0,     0], ...
                [0,             0,             0,             0,  LA1*cos(theta)*thetap,           0,           0,     0]);
            dotphiq(5,1) = -2*(x3p-x1p); dotphiq(5,2) = -2*(y3p-y1p); dotphiq(5,6) = 2*(x3p-x1p); dotphiq(5,7) = 2*(y3p-y1p);
            dotphiq(6,3) = -2*(x3p-x2p); dotphiq(6,4) = -2*(y3p-y2p); dotphiq(6,6) = 2*(x3p-x2p); dotphiq(6,7) = 2*(y3p-y2p);
            dotphiq(7,3) = 0; dotphiq(7,8) = L2B*cos(beta)*betap;
            %--------------
%             dotphiq = [...
%                 2*x1p,        2*y1p,            0,             0,              0,              0,  0, 0;
%                 2*(x1p-x2p),  2*(y1p-y2p), -2*(x1p-x2p),  -2*(y1p-y2p),              0,              0,  0, 0;
%                 0,             0,       2*x2p,         2*y2p,              0,              0,  0, 0;
%                 mbe_iff(abs(sin(theta)) < 0.7,...
%                 [0,             0,           0,             0, LA1*sin(theta)*thetap,              0,  0, 0], ...
%                 [0,             0,           0,             0, LA1*cos(theta)*thetap,              0,  0, 0]  ...
%                 ); ...
%                 %0,             0,           0,             0,              0, LA1*sin(theta)*thetap, 0;
%                 -2*(x3p-x1p),             -2*(y3p-y1p),             0,              0,              0,             2*(x3p-x1p),            2*(y3p-y1p),   0;...
%                 0,             0,             -2*(x3p-x2p),              -2*(y3p-y2p),              0,             2*(x3p-x2p),            2*(y3p-y2p),    0;...
%                 0,   0,  0,  0,  0,  0,  0,  L2B*cos(beta)*betap
%                 ];
            
            phiqpqp = dotphiq * qp;
        end % jacob_phiqp_times_qp
        
        % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
        function Rq = jacob_Rq(me,q,R)
            theta = q(5);
            LA1 = me.bar_lengths(1);
            phiq = me.jacob_phi_q(q);
            phiqd = phiq(:,1:4); % Jacobian dependent part
            
            Phi_qqR = [2*R(1), 2*R(2), 0, 0, 0;
                2*R(1)-2*R(3), 2*(R(2)-R(4)), 2*(-R(1)+R(3)), 2*(-R(2)+R(4)), 0;
                0, 0, 2*R(3), 2*R(4), 0;
                mbe_iff(abs(sin(theta)) < 0.7,...
                [0, 0, 0, 0, LA1*sin(theta)], ...
                [0, 0, 0, 0, LA1*cos(theta)] ...
                ) ...
                ];
            Rd_q = phiqd\(-Phi_qqR);
            Rq = [Rd_q;[0,0,0,0,0]];
        end
        
        
        % Evaluates the instantaneous forces
        function Q = eval_forces(me,q,qp)
            Q_var = zeros(me.dep_coords_count,1);
            Q_var(5) = -me.C*qp(5);
            Q = me.Qg+Q_var;
        end % eval_forces
        
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
            
            switch error_def.error_type
                case 0
                    ini_vel_error = 0;
                    ini_pos_error = 0;
                    grav_error = 0;
                    damping_coef_error = 0;
                    % 1: Gravity:
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
                    ini_vel_error = 0;
                    ini_pos_error = 0;
                    grav_error = 0;
                    damping_coef_error = 10 * error_def.error_scale;
                otherwise
                    error('Unhandled value!');
            end
            bad_model.g = bad_model.g+grav_error; % gravity error
            bad_model.zp_init = bad_model.zp_init+ini_vel_error; % initial velocity error
            bad_model.q_init_aprox(5)=bad_model.q_init_aprox(5)+ini_pos_error; %initial position error
            bad_model.C=bad_model.C+damping_coef_error; %initial position error
            
            % Weight vector
            % WARNING: This vector MUST be updated here, after modifying the "g"
            % vector!
            bad_model=bad_model.update_Qg();
        end % applyErrors
        
        % See docs in base class
        function [] = plot_model_skeleton(me, q, color_code, do_fit)
            plot([me.fixed_points(1),q(1),q(3),me.fixed_points(3)], ...
                [me.fixed_points(2),q(2),q(4),me.fixed_points(4)] ,color_code);
            if (do_fit)
                axis equal;
                xlim ([me.fixed_points(1)-1.2*me.bar_lengths(1),1.1*me.fixed_points(3)]);
                ylim ([me.fixed_points(2)-1.2*me.bar_lengths(1),me.fixed_points(2)+1.2*me.bar_lengths(3)]);
            end
        end
        
    end % methods
end % class
