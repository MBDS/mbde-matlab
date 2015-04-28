classdef mbeMechModelFiveBarsBase < mbeMechModelBase
    % Mechanism physical model: Generic 5 bars linkage (Abstract base
    % class). Derived classes define particular examples of mechanisms with
    % especific masses, lengths, etc.
    % Modeled in Natural coordinates plus two relative angle coordinates 
    % at each fixed end.
    %
    %            2:(q3,q4)
    %               +------------+ 3:(q5,q6)
    %               |            |
    %               |            |
    %               |            |
    %     o---------+            o
    %    (xa,ya) 1:(q1,q2)      (xb,yb)
    %
    % - q7: Angle (xa,ya)-(q1,q2)
    % - q8: Angle (xb,yb)-(q5,q6)
    %
    
    % (Abstract) Read-only, constant properties of the model
    properties(Constant,GetAccess=public)
        % Dependent coordinates count
        dep_coords_count = 8;
        
        % A vector with the indices of the independent coordinates in "q":
        indep_idxs = [7, 8];
    end
    
    % (Abstract) Read-only properties of the model
    properties(GetAccess=public,SetAccess=protected)
        % Initial, approximate position (dep coords) vector
        q_init_aprox=zeros(mbeMechModelFiveBarsBase.dep_coords_count,1);
        
        % Initial velocity for independent coords
        zp_init=[0, 0]';
    end
    
    % Model-specific properties:
    properties(Access=public)
        % Global mass matrix
        M;
        
        % Fixed point coords:
        xA,yA,xB,yB, fixed_points;
        
        % Gravity:
        g = -10;
        
        % lengths: [LA1, L12, L23, L3B]
        bar_lengths; 
        
        % Masses:
        mA1,m12,m23,m3B;
        
        % Force vector (gravity forces only):
        Qg;
    end
    
    methods 
        % Constructor: must be implemented in derived classes to fill-in
        % all mechanical parameters.
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Computes the vector of constraints $\Phi(q)$
        function val = phi(me,q)
            x1 = q(1) ;y1 = q(2); x2 = q(3); y2 = q(4); x3 = q(5); y3 = q(6); ang1 = q(7); ang2 = q(8);
            LA1 = me.bar_lengths(1); L12 = me.bar_lengths(2); L23 = me.bar_lengths(3); L3B = me.bar_lengths(4);
            val = [ (me.xA-x1)^2 + (me.yA-y1)^2 - LA1^2;...
                    (x1-x2)^2 + (y1-y2)^2 - L12^2; ...
                    (x2-x3)^2 + (y2-y3)^2 - L23^2; ...
                    (me.xB-x3)^2 + (me.yB-y3)^2 - L3B^2;...
                  mbe_iff(abs(sin(ang1)) > 0.7,... 
                    x1-me.xA-LA1*cos(ang1), ...
                    y1-me.yA-LA1*sin(ang1)  ...
                    );
                  mbe_iff(abs(sin(ang2)) > 0.7,... 
                    x3-me.xB-L3B*cos(ang2), ...
                    y3-me.yB-L3B*sin(ang2) ...
                    ) ...
                ];
        end % of phi()
        
        % Computes the Jacobian $\Phi_q$
        function phiq = jacob_phi_q(me,q)
            x1 = q(1) ;y1 = q(2); x2 = q(3); y2 = q(4); x3 = q(5); y3 = q(6); ang1 = q(7); ang2 = q(8);
            LA1 = me.bar_lengths(1); L3B = me.bar_lengths(4);
            phiq = [...
                -2*(me.xA-x1), -2*(me.yA-y1),          0,          0,            0,            0,             0,             0; 
                    2*(x1-x2),     2*(y1-y2), -2*(x1-x2), -2*(y1-y2),            0,            0,             0,             0;
                            0,             0,  2*(x2-x3),  2*(y2-y3),   -2*(x2-x3),   -2*(y2-y3),             0,             0;
                            0,             0,          0,          0,-2*(me.xB-x3),-2*(me.yB-y3),             0,             0;
                  mbe_iff(abs(sin(ang1)) > 0.7,... 
                            [1,             0,          0,          0,            0,            0, LA1*sin(ang1),             0],...
                            [0,             1,          0,          0,            0,            0,-LA1*cos(ang1),             0] );
                  mbe_iff(abs(sin(ang2)) > 0.7,... 
                            [0,             0,          0,          0,            1,            0,             0, L3B*sin(ang2)],...
                            [0,             0,          0,          0,            0,            1,             0,-L3B*cos(ang2)] )...
                    ];
                       
        end % jacob_phi_q()

        % Computes the Jacobian $\dot{\Phi_q} \dot{q}$
        function phiqpqp = jacob_phiqp_times_qp(me,q,qp)
            ang1 = q(7); ang2 = q(8);
            x1p = qp(1) ;y1p = qp(2); x2p = qp(3); y2p = qp(4); x3p = qp(5); y3p = qp(6); ang1p = qp(7); ang2p = qp(8);
            LA1 = me.bar_lengths(1); L3B = me.bar_lengths(4);
            dotphiq = [...
                     2*x1p,       2*y1p,            0,            0,            0,            0,             0,             0; 
               2*(x1p-x2p), 2*(y1p-y2p), -2*(x1p-x2p), -2*(y1p-y2p),            0,            0,             0,             0;
                         0,           0,  2*(x2p-x3p),  2*(y2p-y3p), -2*(x2p-x3p), -2*(y2p-y3p),             0,             0;
                         0,           0,            0,            0,        2*x3p,        2*y3p,             0,             0;
                  mbe_iff(abs(sin(ang1)) > 0.7,... 
                         [0,           0,            0,            0,            0,            0, LA1*cos(ang1)*ang1p,             0],...
                         [0,           0,            0,            0,            0,            0, LA1*sin(ang1)*ang1p,             0] ); 
                  mbe_iff(abs(sin(ang2)) > 0.7,... 
                         [0,           0,            0,            0,            0,            0,             0, L3B*cos(ang2)*ang2p],...
                         [0,           0,            0,            0,            0,            0,             0, L3B*sin(ang2)*ang2p])... 
                     ];
                     
            phiqpqp = dotphiq * qp;
        end % jacob_phiqp_times_qp
        
        % Computes the partial derivative of velocity constraints wrt q
        function Phiqqpq = Phiq_times_qp_q(me,q,qp)
%             x1 = q(1); y1 = q(2); x2 = q(3); y2 = q(4); 
            ang1 = q(7); ang2 = q(8);
            LA1 = me.bar_lengths(1); L3B = me.bar_lengths(4);
            xp1 = qp(1); yp1 = qp(2); xp2 = qp(3); yp2 = qp(4); xp3 = qp(5); yp3 = qp(6); ang1p = qp(7); ang2p = qp(8);
            Phiqqpq = [2*xp1, 2*yp1, 0, 0, 0, 0, 0, 0;
                       2*(xp1-xp2), 2*(yp1-yp2), 2*(xp2-xp1), 2*(yp2-yp1),0, 0, 0, 0;
                       0, 0, 2*(xp2-xp3), 2*(yp2-yp3),-2*(xp2-xp3), -2*(yp2-yp3), 0, 0;
                       0, 0, 0, 0, 2*xp3, 2*yp3, 0, 0;
                       mbe_iff(abs(sin(ang1)) < 0.7,... 
                       [0,0,0,0,0,0, LA1*cos(ang1)*ang1p, 0],...
                       [0,0,0,0,0,0, LA1*sin(ang1)*ang1p, 0]);
                       mbe_iff(abs(sin(ang2)) < 0.7,... 
                       [0,0,0,0,0,0,0, L3B*cos(ang2)*ang2p],...
                       [0,0,0,0,0,0,0, L3B*sin(ang2)*ang2p]);];
        end % Phiq_times_qp_q
        
        % Computes the hypermatrix $\frac{\partial R}{\partial q}$, with Rq(:,:,k) the partial derivative of R(q) wrt q(k)
        function Rq = jacob_Rq(me,q,R)
%             error('TO DO!');
            x1 = q(1) ;y1 = q(2); x2 = q(3); y2 = q(4); x3 = q(5); y3 = q(6); ang1 = q(7); ang2 = q(8);
            LA1 = me.bar_lengths(1); L3B = me.bar_lengths(4);
            
            phiq = me.jacob_phi_q(q);
            phiqd = phiq(:,1:6); % Jacobian dependent part
            zeros21 = zeros(1,2);
            Phi_qqR = zeros(6,8,2);
            Phi_qx1R= [2*R(1,1), 2*R(1,2); 
                       2*(R(1,1)-R(3,1)), 2*(R(1,2)-R(3,2));  
                       zeros21;
                       zeros21;
                       zeros21;
                       zeros21];
            
            Phi_qy1R= [2*R(2,1), 2*R(2,2); 
                       2*(R(2,1)-R(4,1)), 2*(R(2,2)-R(4,2));  
                       zeros21;
                       zeros21;
                       zeros21;
                       zeros21];
                   
            Phi_qx2R= [zeros21;
                       2*(R(3,1)-R(1,1)), 2*(R(3,2)-R(1,2)); 
                       2*(R(3,1)-R(5,1)), 2*(R(3,2)-R(5,2));  
                       zeros21;
                       zeros21;
                       zeros21];
            
            Phi_qy2R= [zeros21;
                       2*(R(4,1)-R(2,1)), 2*(R(4,2)-R(2,2)); 
                       2*(R(4,1)-R(6,1)), 2*(R(4,2)-R(6,2));  
                       zeros21;
                       zeros21;
                       zeros21];
            
            Phi_qx3R= [zeros21;
                       zeros21; 
                       2*(R(5,1)-R(3,1)), 2*(R(5,2)-R(3,2)); 
                       2*(R(5,1)), 2*(R(5,2));  
                       zeros21;
                       zeros21];
            
            Phi_qy3R= [zeros21;
                       zeros21;
                       2*(R(6,1)-R(4,1)), 2*(R(6,2)-R(4,2)); 
                       2*(R(6,1)), 2*(R(6,2));  
                       zeros21;
                       zeros21];
                   
                   
            Phi_qang1R= [zeros21;
                       zeros21; 
                       zeros21;
                       zeros21;
                       mbe_iff(abs(sin(ang1)) > 0.7,... 
                         [LA1*cos(ang1), 0],...
                         [LA1*sin(ang1), 0] ); 
                       zeros21 ];
            
            Phi_qang2R= [zeros21;
                       zeros21;
                       zeros21;
                       zeros21;
                       zeros21;
                       mbe_iff(abs(sin(ang2)) > 0.7,... 
                         [0, LA1*cos(ang2)],...
                         [0, LA1*sin(ang2)] ) ];
             
                   
            Phi_qqR(:,1,:) = Phi_qx1R;
            Phi_qqR(:,2,:) = Phi_qy1R;
            Phi_qqR(:,3,:) = Phi_qx2R;
            Phi_qqR(:,4,:) = Phi_qy2R;
            Phi_qqR(:,5,:) = Phi_qx3R;
            Phi_qqR(:,6,:) = Phi_qy3R;
            Phi_qqR(:,7,:) = Phi_qang1R;
            Phi_qqR(:,8,:) = Phi_qang2R;
            
           
            Rd_q1 = phiqd\(-Phi_qqR(:,:,1));
            Rd_q2 = phiqd\(-Phi_qqR(:,:,2));
            Rq = zeros(8,8,2);
            Rq(1:6,:,1) = Rd_q1;
            Rq(1:6,:,2) = Rd_q2;
%             Rd_q = phiqd\(-Phi_qqR);
%             R_q = [Rd_q;[0,0,0,0,0]];            
        end
        

        % Evaluates the instantaneous forces
        function Q = eval_forces(me,q,qp)
            %Q_var = zeros(me.dep_coords_count,1);
            Q = me.Qg; %+Q_var;
        end % eval_forces

        % Evaluates the stiffness & damping matrices of the system:
        function [K, C] = eval_KC(me, q,dq)
            K = zeros(me.dep_coords_count);
            C = zeros(me.dep_coords_count);
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
            if (isprop(bad_model,'C'))
                bad_model.C=bad_model.C+damping_coef_error; %initial position error
            end

            % Weight vector 
            % WARNING: This vector MUST be updated here, after modifying the "g"
            % vector!
            bad_model=bad_model.update_Qg();
        end % applyErrors
        
        % See docs in base class
        function [] = plot_model_skeleton(me, q, color_code, do_fit)
            plot([me.fixed_points(1),q(1),q(3),q(5),me.fixed_points(3)],...
                 [me.fixed_points(2),q(2),q(4),q(6),me.fixed_points(4)],color_code);
            if (do_fit)
                axis equal;
                xlim ([me.fixed_points(1)-1.2*me.bar_lengths(1),1.5*me.fixed_points(3)]);
                ylim ([me.fixed_points(2)-1.8*me.bar_lengths(1),me.fixed_points(2)+1.2*me.bar_lengths(3)]);
            end
        end
        
    end % methods
end % class
