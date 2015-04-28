classdef mbeMechModelPendulum1 < mbeMechModelPendulumBase
    % Mechanism physical model: one-link planar pendulum
    %  See mbeMechModelPendulumBase for a sketch of the mechanism.
    
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
                mbeSensorGyroscope([1 0],[1 2],[1 2], deg2rad(1))... % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
            };
    end
   
    methods 
        % Constructor
        function [me] = mbeMechModelPendulum1()
            % Initial conditions
            x1 = 2; y1= 0; 
            theta = 0; % Initial DOF value
            me.q_init_aprox = [x1, y1, theta]'; % Initial approximated position
            me.zp_init = 0; % Initial DOF velocity
            % gravity
            me.g = -10;

            % Multibody model parameters
            % fixed points:
            me.xA = 0; me.yA = 0; 
            me.fixed_points = [me.xA, me.yA];
            % geometry
            LA1=norm([me.xA,me.yA]-[x1,y1]); 
            me.bar_lengths = [LA1];
            
            % mass
            me.mA1 = 5; 

            % Local Mass matrices. 
            MA1 = me.massMatrixBar(me.mA1,LA1);

            % Global Mass Matrix
            me.M = zeros(3,3);
            me.M(1:2,1:2) = MA1 (3:4,3:4);
          
            % Vector of generalized forces
            me=me.update_Qg();
            
            % damping coefficient
            me.C = 0;             
        end % end ctor()
        
        function [me] = update_Qg(me)
        % Re-calc Qg from gravity vector & masses.
            me.Qg = 0.5*me.g*[0; me.mA1; 0];
        end            
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Evaluates potential energy of the whole system:
        function V = eval_potential_energy(me,q)
           V = -0.5*me.g*(me.mA1*(q(2)+me.yA));
        end
        
    end % methods
end % class
