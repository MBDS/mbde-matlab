classdef mbeSensorGyroscope < mbeSensorBase
    % Sensor: angular velocity of a 2D body from a mechanism, computed from 
    % two selected points in the solid. The constructor must specify those
    % two points, and whether each of them is "fixed" or part of "q".
    
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
    properties
        pt_is_fixed = [0 0];  % Is point {1,2} fixed? (0:part of q, 1:is fixed)
        pt1_idxs    = [1 3];  % Indices in "q"/"fixed_points" of {x,y} for point #1
        pt2_idxs    = [2 4];  % Indices in "q"/"fixed_points" of {x,y} for point #2
                
        std_dev;  % Noise stddev (rad/s)
    end
    
    methods
        function [me]=mbeSensorGyroscope(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std)
            % Constructor. See properties docs of mbeSensorGyroscope for the
            % meaning of each argument.
            me.pt_is_fixed = pt_is_fixed;
            me.pt1_idxs=pt1_idxs;
            me.pt2_idxs=pt2_idxs;
            
            me.std_dev = noise_std;
        end

        function [p1,p2,v1,v2,L2] = extract_pos_vel_vectors(me,model,q,qp)
            % Build the 2D velocity vector & position of points 1 & 2:
            if (me.pt_is_fixed(1))
                v1=[0, 0];
                p1=[model.fixed_points(me.pt1_idxs(1)), model.fixed_points(me.pt1_idxs(2))];
            else
                v1=[qp(me.pt1_idxs(1)), qp(me.pt1_idxs(2))];
                p1=[ q(me.pt1_idxs(1)),  q(me.pt1_idxs(2))];
            end
            if (me.pt_is_fixed(2))
                v2=[0, 0];
                p2=[model.fixed_points(me.pt2_idxs(1)), model.fixed_points(me.pt2_idxs(2))];
            else
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
            [p1,p2,v1,v2,L2] = extract_pos_vel_vectors(me,model,q,qp);
            
            % Calc angular speed:
            % x1 y1 x2 y2
            %  1  2  3  4
            %obs = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
            obs = ((v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)))/L2;
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

            [p1,p2,v1,v2,L2] = extract_pos_vel_vectors(me,model,q,qp);
            numw = (v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)); %numerator of the angular velocity
            % pt1!=fixed
            if (me.pt_is_fixed(1)==0)
                %obs = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
                % dz_dx1 = (v1y - v2y)/L^2
                % dz_dy1 = -(v1x - v2x)/L^2
                % dz_dv1x = -(y1 - y2)/L^2
                % dz_dv1y = (x1 - x2)/L^2
%                 dh_dq (me.pt1_idxs(1)) = (v1(2) - v2(2))/L2;
%                 dh_dq (me.pt1_idxs(2)) = -(v1(1) - v2(1))/L2;
                dL2_dP1x = -2*(p2(1)-p1(1));
                dL2_dP1y = -2*(p2(2)-p1(2));
                
                dnumw_dP1x = -(v2(2)-v1(2));
                dnumw_dP1y = (v2(1)-v1(1));
                
                dh_dq (me.pt1_idxs(1)) =        (dnumw_dP1x*L2-numw*dL2_dP1x)/L2^2;
                dh_dq (me.pt1_idxs(2)) =        (dnumw_dP1y*L2-numw*dL2_dP1y)/L2^2;
                dh_dqp(me.pt1_idxs(1)) = -(p1(2) - p2(2))/L2;
                dh_dqp(me.pt1_idxs(2)) = (p1(1) - p2(1))/L2;
            end
            % pt2!=fixed
            if (me.pt_is_fixed(2)==0)
                % dz_dx2 =  -(v1y - v2y)/L^2
                % dz_dy2 =   (v1x - v2x)/L^2
                % dz_dv2x =     (y1 - y2)/L^2
                % dz_dv2y =    -(x1 - x2)/L^2
%                 dh_dq (me.pt2_idxs(1)) = -(v1(2) - v2(2))/L2;
%                 dh_dq (me.pt2_idxs(2)) =  (v1(1) - v2(1))/L2;
                dL2_dP2x = 2*(p2(1)-p1(1));
                dL2_dP2y = 2*(p2(2)-p1(2));
                
                dnumw_dP2x = (v2(2)-v1(2));
                dnumw_dP2y = -(v2(1)-v1(1));
                
                dh_dq (me.pt2_idxs(1)) = (dnumw_dP2x*L2-numw*dL2_dP2x)/L2^2;
                dh_dq (me.pt2_idxs(2)) = (dnumw_dP2y*L2-numw*dL2_dP2y)/L2^2;
                dh_dqp(me.pt2_idxs(1)) =  (p1(2) - p2(2))/L2;
                dh_dqp(me.pt2_idxs(2)) = -(p1(1) - p2(1))/L2;
            end
            
        end % sensor_jacob
        
    end
    
end

