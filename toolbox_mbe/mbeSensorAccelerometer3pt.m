classdef mbeSensorAccelerometer3pt < mbeSensorBase
    % Sensor: acceleration of a 2D body from a mechanism, computed from 
    % two selected points in the solid + a third point where the sensor is placed at. 
    % The constructor must specify those three points, and whether the first two are "fixed" or part of "q".
    % TODO: For now, the acceleration is the direction of +X local.
    
    properties
        pt_is_fixed = [0 0];  % Is point {1,2} fixed? (0:part of q, 1:is fixed)
        pt1_idxs    = [1 3];  % Indices in "q"/"fixed_points" of {x,y} for point #1
        pt2_idxs    = [2 4];  % Indices in "q"/"fixed_points" of {x,y} for point #2
        pt3_idxs    = [2 4];  % Indices in "q"/"fixed_points" of {x,y} for point #2
                
        std_dev;  % Noise stddev (m/s^2)
    end
    
    methods
        function [me]=mbeSensorAccelerometer3pt(pt_is_fixed,pt1_idxs,pt2_idxs,pt3_idxs, noise_std)
            % Constructor. See properties docs of mbeSensorAccelerometer3pt for the
            % meaning of each argument.
            me.pt_is_fixed = pt_is_fixed;
            me.pt1_idxs=pt1_idxs;
            me.pt2_idxs=pt2_idxs;
            me.pt3_idxs=pt3_idxs;
            
            me.std_dev = noise_std;
        end

        function [p1,p2,p3,v1,v2,v3,a1,a2,a3,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp)
            % Build the 2D velocity & acc vector & position of points 1 & 2:
            if (me.pt_is_fixed(1))
                a1=[0, 0];
                v1=[0, 0];
                p1=[model.fixed_points(me.pt1_idxs(1)), model.fixed_points(me.pt1_idxs(2))];
            else
                a1=[qpp(me.pt1_idxs(1)), qpp(me.pt1_idxs(2))];
                v1=[qp(me.pt1_idxs(1)), qp(me.pt1_idxs(2))];
                p1=[ q(me.pt1_idxs(1)),  q(me.pt1_idxs(2))];
            end
            if (me.pt_is_fixed(2))
                a2=[0, 0];
                v2=[0, 0];
                p2=[model.fixed_points(me.pt2_idxs(1)), model.fixed_points(me.pt2_idxs(2))];
            else
                a2=[qpp(me.pt2_idxs(1)), qpp(me.pt2_idxs(2))];
                v2=[qp(me.pt2_idxs(1)), qp(me.pt2_idxs(2))];
                p2=[ q(me.pt2_idxs(1)),  q(me.pt2_idxs(2))];
            end
            % Pt3:
            a3=[qpp(me.pt3_idxs(1)), qpp(me.pt3_idxs(2))];
            v3=[qp(me.pt3_idxs(1)), qp(me.pt3_idxs(2))];
            p3=[ q(me.pt3_idxs(1)),  q(me.pt3_idxs(2))];
            
            L2=sum( (p2-p1).^2);
            assert(L2>0,'Gyro: the two points must be separated!');
        end
        
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        function [obs] = sensor_simulate(me,model,q,qp,qpp)
            % Build the 2D velocity vector & position of points 1 & 2:
            [p1,p2,p3,v1,v2,v3,a1,a2,a3,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);

            % unit vectors (u,v) defining the local coords system: 
            L=sqrt(L2);
            u_vect = [ p2(1)-p1(1); p2(2)-p1(2) ]./L;   % ( Ax, Ay)
            %v_vect = [-(p2(2)-p1(2)) ; p2(1)-p1(1)]./L; % (-Ay, Ax)
            
            % 
            acc_lx = dot(a3',u_vect);
            
            obs = acc_lx;
        end
        
        % The standard deviation of this sensor
        function [sensor_std] = sensor_std_noise(me)
            sensor_std = me.std_dev;
        end
        
        % Evaluates the sensor model jacobian for a given dyn state:
        %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
        % H= [  ...   |   ...   |    ...   ] = [ dh_dq , dh_dqp, dh_dqpp ]
        %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
        function [dh_dq , dh_dqp, dh_dqpp] = sensor_jacob(me,model,q,qp,qpp)
            n=length(q);
            dh_dq=zeros(1,n);
            dh_dqp=zeros(1,n);
            dh_dqpp=zeros(1,n);

            [p1,p2,p3,v1,v2,v3,a1,a2,a3,L2]  = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);
            
            % unit vectors (u,v) defining the local coords system: 
            L=sqrt(L2);
            u_vect = [ p2(1)-p1(1); p2(2)-p1(2) ]./L;   % ( Ax, Ay)
            %v_vect = [-(p2(2)-p1(2)) ; p2(1)-p1(1)]./L; % (-Ay, Ax)
            
            a3x = a3(1);
            a3y = a3(2);
            
            % pt1!=fixed
            if (me.pt_is_fixed(1)==0)
                dh_dq (me.pt1_idxs(1)) = -1*a3x/L2; % + 0*a3y
                dh_dq (me.pt1_idxs(2)) = -1*a3y/L2;
            end
            % pt2!=fixed
            if (me.pt_is_fixed(2)==0)
                dh_dq (me.pt2_idxs(1)) = 1*a3x/L2;
                dh_dq (me.pt2_idxs(2)) = 1*a3y/L2;
            end
            
            % dh_dqpp part:
            dh_dqpp (me.pt3_idxs(1)) = u_vect(1);
            dh_dqpp (me.pt3_idxs(2)) =  u_vect(2);
            
            
        end % sensor_jacob
        
    end
    
end

