classdef mbeMechModelRodArray1dof < mbeMechModelBase
    % Academic example of a mechanism: 3 identical & parallel rods
    % connected with one solid body.
    %
    %
    %         1:(q1,q2)  2:(q3,q4)  3:(q5,q6)
    %         +----------+----------+ <--- 2*m_bar
    %        /          /          /
    %   L-> /          /          / <-- m_bar
    %      /          /          /  
    %     o          o          o
    %    (0,0)      (D,0)     (2*D,0)
    %
    %  %%%- q7: Angle (0,0)-(q1,q2)
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
    
    % (Abstract) Read-Write properties
    properties(Access=public)
        % List of installed sensors (cell of objects derived from mbeSensorBase)
        installed_sensors = { ...
            	mbeSensorPosIndex(1, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
            };
    end
    
    % (Abstract) Read-only, constant properties of the model
    properties(Constant,GetAccess=public)
        % Dependent coordinates count
        dep_coords_count = 6;
        
        % A vector with the indices of the independent coordinates in "q":
        indep_idxs = 1;
        
    end
    % (Abstract) Read-only properties of the model
    properties(GetAccess=public,SetAccess=public)
        % Initial, approximate position (dep coords) vector
        q_init_approx=zeros(mbeMechModelRodArray1dof.dep_coords_count,1);
        
        
        % Initial velocity for independent coords
        zp_init=[0];
    end
    
    % Model-specific properties:
    properties(Access=public)
        % Global mass matrix
        M;
        
        L; % Length of bars
        D; % Distance between pinned-joints
        
        % Gravity:
        g = -10;
        
        fixed_points;
        
        % Mass of each bar:
        m_bar;
        
        % Force vector (gravity forces only):
        Qg;

    end
    
   
    methods 
        % Constructor
        function [me] = mbeMechModelRodArray1dof()
            me.L = 1;
            me.D = 2;
            me.m_bar = 1.0;
            
            me.fixed_points=[0,0, me.D,0, 2*me.D,0]';
            
            theta = pi/4;
            ct = cos(theta); st=sin(theta);
            
            % Initial approximated position
            me.q_init_approx = [...
                me.L*ct        , me.L*st, ...
                me.D+me.L*ct   , me.L*st, ...
                2*me.D+me.L*ct , me.L*st ...                
                ]'; 
            
            me.zp_init = 0; % Initial DOF velocity
            
            % gravity
            me.g = -10;
            
            % Global Mass Matrix
            me.M = zeros(6,6);
            end1=1:2; end2=3:4;
            p1=1:2; p2=3:4; p3=5:6;
            
            if (0)
                % Way 1:
                MA1  = me.massMatrixBar(me.m_bar,me.L);
                M13  = me.massMatrixBar(2*me.m_bar,2*me.D);

                me.M(p1,p1) = MA1(end2,end2);  % pt: 1
                me.M(p2,p2) = MA1(end2,end2);  % pt: 2
                me.M(p3,p3) = MA1(end2,end2);  % pt: 3
                me.M([p1 p3],[p1 p3]) = me.M([p1 p3],[p1 p3]) + M13;                
            else
                % Way 2:
                MA1  = me.massMatrixBar(me.m_bar,me.L);
                M12  = me.massMatrixBar(me.m_bar,me.D);
                M23  = M12;

                me.M(p1,p1) = MA1(end2,end2);  % pt: 1
                me.M(p2,p2) = MA1(end2,end2);  % pt: 2
                me.M(p3,p3) = MA1(end2,end2);  % pt: 3
                me.M([p1 p2],[p1 p2]) = me.M([p1 p2],[p1 p2]) + M12;
                me.M([p2 p3],[p2 p3]) = me.M([p2 p3],[p2 p3]) + M23;
            end
          
            % Vector of generalized forces
            me=me.update_Qg();
            
        end % end ctor()
        
        function [me] = update_Qg(me)
            % Re-calc Qg from gravity vector & masses.
            %me.Qg = me.g*[0; 0.5*me.m_bar; 0; 0.5*me.m_bar + 2*me.m_bar; 0 ; 0.5*me.m_bar];
            me.Qg = me.g*[0; 0.5*me.m_bar+(1/3)*2*me.m_bar; 0; 0.5*me.m_bar + (1/3)*2*me.m_bar; 0 ; 0.5*me.m_bar+(1/3)*2*me.m_bar];
        end            
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Evaluates potential energy of the whole system:
        function V = eval_potential_energy(me,q)
            y1 = q(2);
            V = -me.g*( ...
                3 * me.m_bar*0.5*y1 + ...
                (2*me.m_bar)*y1 );
        end
        
    end % methods
        
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Computes the vector of constraints $\Phi(q)$
        function val = phi(me,q)
            x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4); x3 = q(5);y3 = q(6); %theta = q(7);
            val = [...
                (x1-0*me.D)^2 + y1^2 - me.L^2;
                (x2-1*me.D)^2 + y2^2 - me.L^2;
                (x3-2*me.D)^2 + y3^2 - me.L^2;
                (x2-x1)^2 + (y2-y1)^2 - me.D^2;
                (x3-x1)-2*(x2-x1);
                (y3-y1)-2*(y2-y1)
                    ];
        end % of phi()
        
        % Computes the Jacobian $\Phi_q$
        function phiq = jacob_phi_q(me,q)
            x1 = q(1) ;y1 = q(2);x2 = q(3);y2 = q(4); x3 = q(5);y3 = q(6);

            phiq = [...
                 2*(x1-0*me.D),      2*y1,             0,           0,             0,             0;
                             0,         0, 2*(x2-1*me.D),        2*y2,             0,             0;
                             0,         0,             0,           0, 2*(x3-2*me.D),          2*y3;
                    -2*(x2-x1),-2*(y2-y1),     2*(x2-x1),   2*(y2-y1),             0,             0;
                             1,         0,            -2,           0,             1,             0;
                             0,         1,            0,           -2,             0,             1;
                    ];
        end % jacob_phi_q()

        % Computes the Jacobian $\dot{\Phi_q} \dot{q}$
        function phiqpqp = jacob_phiqp_times_qp(me,q,qp)
            x1p = qp(1) ;y1p = qp(2);x2p = qp(3);y2p = qp(4); x3p = q(5);y3p = qp(6);

            dotphiq = [...
                 2*(x1p),      2*y1p,             0,           0,             0,             0;
                             0,         0, 2*(x2p),        2*y2p,             0,             0;
                             0,         0,             0,           0, 2*(x3p),          2*y3p;
                    -2*(x2p-x1p),-2*(y2p-y1p),     2*(x2p-x1p),   2*(y2p-y1p),             0,             0;
                             0,         0,            0,           0,             0,             0;
                             0,         0,            0,           0,             0,             0;
                    ];

            phiqpqp = dotphiq * qp;
        end % jacob_phiqp_times_qp
        
        % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
        function Rq = jacob_Rq(me,q,R)
            error('todo');
        end
        

        % Evaluates the instantaneous forces
        function Q = eval_forces(me,q,qp)
            Q = me.Qg;
        end % eval_forces

        % Evaluates the stiffness & damping matrices of the system:
        function [K, C] = eval_KC(me, q,dq)
            K = zeros(me.dep_coords_count,me.dep_coords_count);
            C = zeros(me.dep_coords_count,me.dep_coords_count);
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
            bad_model.q_init_approx(5)=bad_model.q_init_approx(5)+ini_pos_error; %initial position error
            bad_model.C=bad_model.C+damping_coef_error; %initial position error

            % Weight vector 
            % WARNING: This vector MUST be updated here, after modifying the "g"
            % vector!
            bad_model=bad_model.update_Qg();
        end % applyErrors
        
        % See docs in base class
        function [] = plot_model_skeleton(me, q, color_code, do_fit)
            plot([0,q(1)], [0,q(2)], color_code);
            plot([me.D,q(3)], [0,q(4)], color_code);
            plot([2*me.D,q(5)], [0,q(6)], color_code);
            plot([q(1),q(5)], [q(2),q(6)], color_code);
            if (do_fit)
                axis equal;
            end
        end
        
    end % methods
    
end % class
