classdef mbeSensorAccelerometer < mbeSensorBase
    % Sensor: output from a 1 axis accelerometer. The acceleration is
    % calculated from kinematic magnitudes of two points of a solid.
    % Gravity is 9.81 m/s^2 pointing downwards by default (global y axis). 
    
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
        axisangle   = 0;      % Angle in radians of the sensitive axis of the accelerometer wrt the line defined by pt2-pt1  
        LCC         = [0;0];  % Linear combination coefficients to define the position of the accelerometer.
        dist_acc    = 0; 
        ang_r_acc   = 0;      % dist_acc and ang_r_acc provide the same information than LCC but in polar coordinates
        % The first unit vector points from pt1 to pt2. The second forms a
        % right-handed system if we consider the third vector pointing
        % towards the reader.
        gravity     = [0;-9.81]; % Y axis points upwards by default.
        std_dev;  % Noise stddev (rad/s)
    end
    
    methods
        function [me]=mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
            % Constructor. See properties docs of mbeSensorGyroscope for the
            % meaning of each argument. Particular arguments of this
            % function are axisangle, which sets the angle of the
            % sensitive axis of the accelerometer, and gravity (optional) ,
            % wich is set to [0,-8.81, 0] by defaulty, but can be changed.
            me.pt_is_fixed = pt_is_fixed;
            me.pt1_idxs=pt1_idxs;
            me.pt2_idxs=pt2_idxs;
            if ~isempty(LCC)
                me.LCC = LCC;
            end
            me.dist_acc = norm(me.LCC);
            me.ang_r_acc = atan2(me.LCC(2),me.LCC(1));
            
            if ~isempty(axisangle)
                me.axisangle = axisangle;
            end
            if ~isempty(gravity)
                me.gravity = gravity;
            end
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
            assert(L2>0,'Accelerometer: the two points must be separated!');
        end
        
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        function [obs] = sensor_simulate(me,model,q,qp,qpp)
            % Build the 2D velocity vector & position of points 1 & 2:
            [p1,p2,v1,v2,a1,a2,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);
            
            % Calc angular speed:
            % x1 y1 x2 y2
            %  1  2  3  4
            %w = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
            w = ((v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)))/L2;
            % Angular acceleration
            alpha = ((a2(1)-a1(1))*(p1(2)-p2(2))+(a2(2)-a1(2))*(p2(1)-p1(1)))/L2;
            % Calculate the axes of the local base. UVect1 is a unit vector
            % that points from p2 to p1. UVect2 is a unit vector normal to
            % UVect1 and it forms a right-handed basis if the third vector
            % is pointing towards the reader.
            UVect1 = (p2'-p1')/sqrt(L2);
            UVect2 = [-UVect1(2); UVect1(1)];
            % Calculate the coordinates of the accelerometer in the local
            % base
            unit_r_acc = UVect1*cos(me.ang_r_acc)+UVect2*sin(me.ang_r_acc);
            norm_unit_r_acc = [-unit_r_acc(2);unit_r_acc(1)];
            
            % Acceleration sensed by the accelerometer (acceleration minus
            % acceleration of gravity)
            accelerometer = a1' + alpha*me.dist_acc*norm_unit_r_acc - w^2*me.dist_acc*unit_r_acc - me.gravity;
            
            % Projection over accelerometes axis
            acc_axis = UVect1*cos(me.axisangle)+UVect2*sin(me.axisangle);
            obs = accelerometer'*acc_axis;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obs = w;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            dalpha_dq = zeros(1,n);
            dalpha_dqpp = zeros(1,n);
            dw_dq = zeros(1,n);
            dw_dqp = zeros(1,n);
            dUVect1_dq = zeros(2,n);
            dUVect2_dq = zeros(2,n);
            da1_dqpp = zeros(2,n);
         
            [p1,p2,v1,v2,a1,a2,L2] = extract_pos_vel_acc_vectors(me,model,q,qp,qpp);
            L = sqrt(L2);
            L3 = (sqrt(L2))^3;
            
            %w = ((v2x-v1x)*(y1-y2)+(v2y-v1y)*(x2-x1))/L^2;
            w = ((v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)))/L2;
            % Angular acceleration
            alpha = ((a2(1)-a1(1))*(p1(2)-p2(2))+(a2(2)-a1(2))*(p2(1)-p1(1)))/L2;
            % Calculate the axes of the local base. UVect1 is a unit vector
            % that points from p2 to p1. UVect2 is a unit vector normal to
            % UVect1 and it forms a right-handed basis if the third vector
            % is pointing towards the reader.
            UVect1 = (p2'-p1')/sqrt(L2);
            UVect2 = [-UVect1(2); UVect1(1)];
            % Calculate the coordinates of the accelerometer in the local
            % base
            unit_r_acc = UVect1*cos(me.ang_r_acc)+UVect2*sin(me.ang_r_acc);
            norm_unit_r_acc = [-unit_r_acc(2);unit_r_acc(1)];
            
            % Acceleration sensed by the accelerometer (acceleration minus
            % acceleration of gravity)
            accelerometer = a1' + alpha*me.dist_acc*norm_unit_r_acc - w^2*me.dist_acc*unit_r_acc - me.gravity;
            
            % Projection over accelerometes axis
            acc_axis = UVect1*cos(me.axisangle)+UVect2*sin(me.axisangle);
            
            numw = (v2(1)-v1(1))*(p1(2)-p2(2))+(v2(2)-v1(2))*(p2(1)-p1(1)); %numerator of the angular velocity
            numalpha = (a2(1)-a1(1))*(p1(2)-p2(2))+(a2(2)-a1(2))*(p2(1)-p1(1)); %numerator of the angular acceleration
            % pt1!=fixed
            if (me.pt_is_fixed(1)==0)
                
                dL2_dP1x = -2*(p2(1)-p1(1));
                dL2_dP1y = -2*(p2(2)-p1(2));
                
                dnumw_dP1x = -(v2(2)-v1(2));
                dnumw_dP1y = (v2(1)-v1(1));
                
                dnumalpha_dP1x = -(a2(2)-a1(2));
                dnumalpha_dP1y = (a2(1)-a1(1));
                                
%                 dalpha_dq(me.pt1_idxs(1)) =     (a1(2) - a2(2))/L2;
%                 dalpha_dq(me.pt1_idxs(2)) =     -(a1(1) - a2(1))/L2;
                dalpha_dq(me.pt1_idxs(1)) =     (dnumalpha_dP1x*L2-numalpha*dL2_dP1x)/L2^2;
                dalpha_dq(me.pt1_idxs(2)) =     (dnumalpha_dP1y*L2-numalpha*dL2_dP1y)/L2^2;
                dalpha_dqpp(me.pt1_idxs(1)) =   -(p1(2) - p2(2))/L2;
                dalpha_dqpp(me.pt1_idxs(2)) =   (p1(1) - p2(1))/L2;
                
%                 dw_dq (me.pt1_idxs(1)) =        (v1(2) - v2(2))/L2;
%                 dw_dq (me.pt1_idxs(2)) =        -(v1(1) - v2(1))/L2;
                dw_dq (me.pt1_idxs(1)) =        (dnumw_dP1x*L2-numw*dL2_dP1x)/L2^2;
                dw_dq (me.pt1_idxs(2)) =        (dnumw_dP1y*L2-numw*dL2_dP1y)/L2^2;
                dw_dqp(me.pt1_idxs(1)) =        -(p1(2) - p2(2))/L2;
                dw_dqp(me.pt1_idxs(2)) =        (p1(1) - p2(1))/L2;
                
%                 dUVect1_dq(1,me.pt1_idxs(1)) = -(4*p1(1)*L2+1)/(4*L3);
%                 dUVect1_dq(1,me.pt1_idxs(2)) = -(p2(1)-p1(1))/(4*L3*(p2(2)-p1(2)));
%                 dUVect1_dq(2,me.pt1_idxs(1)) = -(p2(2)-p1(2))/(4*L3*(p2(1)-p1(1)));
%                 dUVect1_dq(2,me.pt1_idxs(2)) = -(4*p1(2)*L2+1)/4*L3;
                dUVect1_dq(1,me.pt1_idxs(1)) = (-L+(p2(1)-p1(1))^2/L)/L2;
                dUVect1_dq(1,me.pt1_idxs(2)) = ((p2(1)-p1(1))*(p2(2)-p1(2)))/L3;
                dUVect1_dq(2,me.pt1_idxs(1)) = dUVect1_dq(1,me.pt1_idxs(2));
                dUVect1_dq(2,me.pt1_idxs(2)) = (-L+(p2(2)-p1(2))^2/L)/L2;
                
                dUVect2_dq(1,me.pt1_idxs(1)) = -dUVect1_dq(2,me.pt1_idxs(1));
                dUVect2_dq(1,me.pt1_idxs(2)) = -dUVect1_dq(2,me.pt1_idxs(2));
                dUVect2_dq(2,me.pt1_idxs(1)) =  dUVect1_dq(1,me.pt1_idxs(1));
                dUVect2_dq(2,me.pt1_idxs(2)) =  dUVect1_dq(1,me.pt1_idxs(2));
                
                
                da1_dqpp(1,me.pt1_idxs(1))  = 1; 
                da1_dqpp(2,me.pt1_idxs(2))  = 1; 
                
            end
            % pt2!=fixed
            if (me.pt_is_fixed(2)==0)
                dL2_dP2x = 2*(p2(1)-p1(1));
                dL2_dP2y = 2*(p2(2)-p1(2));
                
                dnumw_dP2x = (v2(2)-v1(2));
                dnumw_dP2y = -(v2(1)-v1(1));
                
                dnumalpha_dP2x = (a2(2)-a1(2));
                dnumalpha_dP2y = -(a2(1)-a1(1));
                
%                 dalpha_dq (me.pt2_idxs(1)) = -(a1(2) - a2(2))/L2;
%                 dalpha_dq (me.pt2_idxs(2)) =  (a1(1) - a2(1))/L2;
                dalpha_dq (me.pt2_idxs(1)) = (dnumalpha_dP2x*L2-numalpha*dL2_dP2x)/L2^2;
                dalpha_dq (me.pt2_idxs(2)) =  (dnumalpha_dP2y*L2-numalpha*dL2_dP2y)/L2^2;
                dalpha_dqpp(me.pt2_idxs(1)) =  (p1(2) - p2(2))/L2;
                dalpha_dqpp(me.pt2_idxs(2)) = -(p1(1) - p2(1))/L2;
                
%                 dw_dq (me.pt2_idxs(1)) = -(v1(2) - v2(2))/L2;
%                 dw_dq (me.pt2_idxs(2)) =  (v1(1) - v2(1))/L2;
                dw_dq (me.pt2_idxs(1)) = (dnumw_dP2x*L2-numw*dL2_dP2x)/L2^2;
                dw_dq (me.pt2_idxs(2)) = (dnumw_dP2y*L2-numw*dL2_dP2y)/L2^2;
                dw_dqp(me.pt2_idxs(1)) =  (p1(2) - p2(2))/L2;
                dw_dqp(me.pt2_idxs(2)) = -(p1(1) - p2(1))/L2;
                
%                 dUVect1_dq(1,me.pt2_idxs(1)) = (4*p2(2)*L2+1)/(4*L3);
%                 dUVect1_dq(1,me.pt2_idxs(2)) = (p2(2)-p1(2))/(4*L3*(p2(2)-p1(2)));
%                 dUVect1_dq(2,me.pt2_idxs(1)) = (p2(2)-p1(2))/(4*L3*(p2(1)-p1(1)));
%                 dUVect1_dq(2,me.pt2_idxs(2)) = (4*p2(2)*L2+1)/(4*L3);
                
                dUVect1_dq(1,me.pt2_idxs(1)) = (L-(p2(1)-p1(1))^2/L)/L2;
                dUVect1_dq(1,me.pt2_idxs(2)) = -((p2(1)-p1(1))*(p2(2)-p1(2)))/L3;
                dUVect1_dq(2,me.pt2_idxs(1)) = dUVect1_dq(1,me.pt2_idxs(2));
                dUVect1_dq(2,me.pt2_idxs(2)) = (L-(p2(2)-p1(2))^2/L)/L2;
                
                dUVect2_dq(1,me.pt2_idxs(1)) = -dUVect1_dq(2,me.pt2_idxs(1));
                dUVect2_dq(1,me.pt2_idxs(2)) = -dUVect1_dq(2,me.pt2_idxs(2));
                dUVect2_dq(2,me.pt2_idxs(1)) =  dUVect1_dq(1,me.pt2_idxs(1));
                dUVect2_dq(2,me.pt2_idxs(2)) =  dUVect1_dq(1,me.pt2_idxs(2));
            end
            % Partial derivative of h wrt q
            dacc_axis_dq =      cos(me.axisangle)*dUVect1_dq+sin(me.axisangle)*dUVect2_dq;
            dr_acc_dq =         me.LCC(1)*dUVect1_dq+me.LCC(2)*dUVect2_dq;
            dnorm_r_acc_dq =   -me.LCC(2)*dUVect1_dq+me.LCC(1)*dUVect2_dq;
            daccelerometer_dq = me.dist_acc*norm_unit_r_acc*dalpha_dq + alpha*dnorm_r_acc_dq-(2*w*me.dist_acc*unit_r_acc*dw_dq+w^2*dr_acc_dq);
            dh_dq = acc_axis'*daccelerometer_dq+accelerometer'*dacc_axis_dq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             dh_dq = dw_dq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Partial derivative of h wrt qp
            dh_dqp = -2*w*me.dist_acc*unit_r_acc'*acc_axis*dw_dqp;       
            
            % Partial derivative of h wrt qpp
            dh_dqpp = acc_axis'*(da1_dqpp + norm_unit_r_acc*dalpha_dqpp*me.dist_acc);
            
        end % sensor_jacob
        
    end
    
end

