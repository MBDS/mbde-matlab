classdef mbeSensorAngle < mbeSensorBase
    % Sensor: absolute orientation of one 2D body from a mechanism, computed from 
    % two selected points in the solid. The constructor must specify those
    % two points, and whether each of them is "fixed" or part of "q".
    
    properties
        pt_is_fixed = [0 0];  % Is point {1,2} fixed? (0:part of q, 1:is fixed)
        pt1_idxs    = [1 3];  % Indices in "q"/"fixed_points" of {x,y} for point #1
        pt2_idxs    = [2 4];  % Indices in "q"/"fixed_points" of {x,y} for point #2
                
        std_dev;  % Noise stddev (rad)
    end
    
    methods
        function [me]=mbeSensorAngle(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std)
            % Constructor. See properties docs of mbeSensorGyroscope for the
            % meaning of each argument.
            me.pt_is_fixed = pt_is_fixed;
            me.pt1_idxs=pt1_idxs;
            me.pt2_idxs=pt2_idxs;
            
            me.std_dev = noise_std;
        end

        function [p1,p2,L2] = extract_pos_vectors(me,model,q)
            % Build the 2D position of points 1 & 2:
            if (me.pt_is_fixed(1))
                p1=[model.fixed_points(me.pt1_idxs(1)), model.fixed_points(me.pt1_idxs(2))];
            else
                p1=[ q(me.pt1_idxs(1)),  q(me.pt1_idxs(2))];
            end
            if (me.pt_is_fixed(2))
                p2=[model.fixed_points(me.pt2_idxs(1)), model.fixed_points(me.pt2_idxs(2))];
            else
                p2=[ q(me.pt2_idxs(1)),  q(me.pt2_idxs(2))];
            end
            L2=sum( (p2-p1).^2);
            assert(L2>0,'mbeSensorAngle: the two points must be separated!');
        end
        
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        function [obs] = sensor_simulate(me,model,q,qp,qpp)
            % Build the 2D position of points 1 & 2:
            [p1,p2,~] = extract_pos_vectors(me,model,q);
            
            % Calc angle:
            % x1 y1 x2 y2
            %  1  2  3  4
            obs = atan2(p2(2)-p1(2), p2(1)-p1(1));
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

            [p1,p2,L2] = extract_pos_vectors(me,model,q);
            x1 = p1(1); y1 = p1(2);
            x2 = p2(1); y2 = p2(2);
            
            % pt1!=fixed
            if (me.pt_is_fixed(1)==0)
                %obs = atan2(p2(2)-p1(2), p2(1)-p1(1))
                % dz_dx1 = -(y1 - y2)/((x1 - x2)^2 + (y1 - y2)^2)
                % dz_dy1 = (x1 - x2)/((x1 - x2)^2 + (y1 - y2)^2)
                dh_dq (me.pt1_idxs(1)) = -(y1 - y2)/L2;
                dh_dq (me.pt1_idxs(2)) = (x1 - x2)/L2;
            end
            % pt2!=fixed
            if (me.pt_is_fixed(2)==0)
                % dz_dx2 = (y1 - y2)/((x1 - x2)^2 + (y1 - y2)^2)
                % dz_dy2 = -(x1 - x2)/((x1 - x2)^2 + (y1 - y2)^2)
                dh_dq (me.pt2_idxs(1)) = (y1 - y2)/L2;
                dh_dq (me.pt2_idxs(2)) = -(x1 - x2)/L2;
            end
            
        end % sensor_jacob
        
    end
    
end

