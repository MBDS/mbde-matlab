classdef mbeMechModelFiveBars1 < mbeMechModelFiveBarsBase
    % Mechanism physical model: 5 bars linkage ("Example model #1")
    %  See mbeMechModelFiveBarsBase for a sketch of the mechanism.
    
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
        % List of installed sensors (cell of objects derived from mbeSensorBase)
        installed_sensors = { ...
%             	mbeSensorPosIndex(7, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%             	mbeSensorPosIndex(8, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%                mbeSensorGyroscope([1 0],[1 2],[1 2], deg2rad(1))... % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%                mbeSensorGyroscope([0 1],[5 6],[3 4], deg2rad(1))... % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
               mbeSensorGyroscope([0 0],[1 2],[3 4], deg2rad(1))... % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
               mbeSensorGyroscope([0 0],[3 4],[5 6], deg2rad(1))... % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%             	mbeSensorVelIndex(7, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%             	mbeSensorVelIndex(8, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
               };
    end
   
    methods 
        % Constructor
        function [me] = mbeMechModelFiveBars1()
            % Initial conditions
            x1 = 0.5; y1= 0; 
            x2 =   0; y2= 2;
            x3 = 2.5; y3= 0;
            theta1 = 0; % Initial DOF value
            theta2 = pi; % Initial DOF value
            me.q_init_approx = [x1, y1, x2, y2, x3, y3, theta1, theta2]'; % Initial approximated position
            me.zp_init = [0, 0]'; % Initial DOF velocity
            % gravity
            me.g = -10;

            % Multibody model parameters
            % fixed points:
            me.xA = 0; me.yA = 0; 
            me.xB = 3; me.yB = 0;
            me.fixed_points = [me.xA, me.yA, me.xB, me.yB];
            % geometry
            LA1=norm([me.xA,me.yA]-[x1,y1]); 
            L12=norm([x1,y1]-[x2,y2]); 
            L23=norm([x2,y2]-[x3,y3]); 
            L3B=norm([me.xB,me.yB]-[x3,y3]); 
            me.bar_lengths = [LA1, L12, L23, L3B];
            % mass
            me.mA1 = 3; 
            me.m12 = 1; 
            me.m23 = 2; 
            me.m3B = 3;
            
            % Local Mass matrices. 
            MA1 = me.massMatrixBar(me.mA1,LA1);
            M12 = me.massMatrixBar(me.m12,L12);
            M23 = me.massMatrixBar(me.m23,L23);
            M3B = me.massMatrixBar(me.m3B,L3B);
            
            % Global Mass Matrix
            me.M = zeros(me.dep_coords_count);
            end1=1:2; end2=3:4;
            p1=1:2; p2=3:4; p3=5:6;
            
            me.M(p1,p1) = MA1(end2,end2);  % pt: 1
            me.M(p3,p3) = M3B(end1,end1);  % pt: 3
            me.M([p1 p2],[p1 p2]) = me.M([p1 p2],[p1 p2]) + M12;
            me.M([p2 p3],[p2 p3]) = me.M([p2 p3],[p2 p3]) + M23;
                      
            % Vector of generalized forces
            me=me.update_Qg();
            
        end % end ctor()
        
        function [me] = update_Qg(me)
        % Re-calc Qg from gravity vector & masses.
            me.Qg = 0.5*me.g*[...
                0; me.mA1+me.m12; ... % pt 1
                0; me.m12+me.m23; ... % pt 2
                0; me.m23+me.m3B; ... % pt 3
                0; 0 ...              % rel coords (torques)
                ];
        end
        
        function [me] = update_M(me)
            % Recalculate mass matrix if masses and/or lengths are changed
            LA1 = me.bar_lengths(1);
            L12 = me.bar_lengths(2);
            L23 = me.bar_lengths(3);
            L3B = me.bar_lengths(4);
            MA1 = me.massMatrixBar(me.mA1,LA1);
            M12 = me.massMatrixBar(me.m12,L12);
            M23 = me.massMatrixBar(me.m23,L23);
            M3B = me.massMatrixBar(me.m3B,L3B);

            % Global Mass Matrix
            me.M = zeros(me.dep_coords_count);
            end1=1:2; end2=3:4;
            p1=1:2; p2=3:4; p3=5:6;

            me.M(p1,p1) = MA1(end2,end2);  % pt: 1
            me.M(p3,p3) = M3B(end1,end1);  % pt: 3
            me.M([p1 p2],[p1 p2]) = me.M([p1 p2],[p1 p2]) + M12;
            me.M([p2 p3],[p2 p3]) = me.M([p2 p3],[p2 p3]) + M23;
        end
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Evaluates potential energy of the whole system:
        function V = eval_potential_energy(me,q)
            y1 = q(2); y2 = q(4); y3 = q(6);
            V = -0.5*me.g*( ...
               me.mA1*(y1+me.yA)+ ...  % bar A1
               me.m12*(y1+y2) + ...  % bar 12
               me.m23*(y2+y3) + ...  % bar 23
               me.m3B*(y3+me.yB) ...   % bar 3B
               ); 
        end
        
    end % methods
end % class
