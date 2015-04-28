classdef mbeSensorAccelerometer < mbeSensorBase
    % Sensor: acceleration of a 2D body from a mechanism, computed from 
    % two selected points in the solid. The constructor must specify those
    % two points, and whether each of them is "fixed" or part of "q".
    
    properties
        pt_is_fixed = [0 0];  % Is point {1,2} fixed? (0:part of q, 1:is fixed)
        pt1_idxs    = [1 3];  % Indices in "q"/"fixed_points" of {x,y} for point #1
        pt2_idxs    = [2 4];  % Indices in "q"/"fixed_points" of {x,y} for point #2
                
        std_dev;  % Noise stddev (rad/s)
    end
    
    methods
        function [me]=mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std)
            % Constructor. See properties docs of mbeSensorGyroscope for the
            % meaning of each argument.
            me.pt_is_fixed = pt_is_fixed;
            me.pt1_idxs=pt1_idxs;
            me.pt2_idxs=pt2_idxs;
            
            me.std_dev = noise_std;
        end

        function [p1,p2,v1,v2,a1,a2,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp)
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
            L2=sum( (p2-p1).^2);
            assert(L2>0,'Gyro: the two points must be separated!');
        end
        
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        function [obs] = sensor_simulate(me,model,q,qp,qpp)
            % Build the 2D velocity vector & position of points 1 & 2:
            [p1,p2,v1,v2,a1,a2,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);
            
            % Calc angular speed:
            % x1 y1 x2 y2
            %  1  2  3  4
            %obs = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
            w = ((v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)))/L2;
            
            % Calculate the tangential aceleration of pt2 wrt pt1 (in the
            % local coords system defined as +x pt1->pt2, +y orthogonal to +x):
            % acc_tang = a_b_wrt_a \dot v (v=unit vector in local +y)
            %
            % acc_pt2 = acc_pt1 + \alpha \times 12 + w \times w \times 12
            %  =>
            % acc_tang = \alpha \times 12 = (acc_pt2 - acc_pt1 ) - w \times w \times 12
            %
            v_vect = [-(p2(2)-p1(2)) ; p2(1)-p1(1); 0]./norm(p1-p2); % (-Ay, Ax)
            acc_tang = [(a2-a1)';0]; % - times([0;0;w],[(p2-p1)';0]);
            
            % Solve for angular acceleration:
            acc_tang_y_local = dot(v_vect,acc_tang);
            L=sqrt(L2);
            alpha = acc_tang_y_local / L;
            
            
            obs = w;
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

            [p1,p2,v1,v2,a1,a2,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);
            
            % pt1!=fixed
            if (me.pt_is_fixed(1)==0)
                %obs = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
                % dz_dx1 = (v1y - v2y)/L^2
                % dz_dy1 = -(v1x - v2x)/L^2
                % dz_dv1x = -(y1 - y2)/L^2
                % dz_dv1y = (x1 - x2)/L^2
                dh_dq (me.pt1_idxs(1)) = (v1(2) - v2(2))/L2;
                dh_dq (me.pt1_idxs(2)) = -(v1(1) - v2(1))/L2;
                dh_dqp(me.pt1_idxs(1)) = -(p1(2) - p2(2))/L2;
                dh_dqp(me.pt1_idxs(2)) = (p1(1) - p2(1))/L2;
            end
            % pt2!=fixed
            if (me.pt_is_fixed(2)==0)
                % dz_dx2 =  -(v1y - v2y)/L^2
                % dz_dy2 =   (v1x - v2x)/L^2
                % dz_dv2x =     (y1 - y2)/L^2
                % dz_dv2y =    -(x1 - x2)/L^2
                dh_dq (me.pt2_idxs(1)) = -(v1(2) - v2(2))/L2;
                dh_dq (me.pt2_idxs(2)) =  (v1(1) - v2(1))/L2;
                dh_dqp(me.pt2_idxs(1)) =  (p1(2) - p2(2))/L2;
                dh_dqp(me.pt2_idxs(2)) = -(p1(1) - p2(1))/L2;
            end
            
        end % sensor_jacob
        
    end
    
end

