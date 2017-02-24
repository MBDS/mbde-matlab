classdef mbeMechModelMaqueta < mbeMechModelFourBarsExtraBase
    % Mechanism physical model: 4 bars linkage ("Almeria model, #3")
    %  See mbeMechModelFourBarsBase for a sketch of the mechanism.
    
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
%            	mbeSensorPosIndex(5, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
                mbeSensorGyroscope([1 0],[1 2],[1 2], 0)... % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorGyroscope([0 0],[1 2],[3 4], 0)... % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorGyroscope([0 1],[3 4],[3 4], 0)... % Gyro in third link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorAngle([1 0],[1 2],[1 2], 0)... % Gyro (Angular position) in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorAngle([0 0],[1 2],[3 4], 0)... % Gyro (Angular position) in coupler link: See mbeSensorAngle(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorAngle([0 1],[3 4],[3 4], 0)... % Gyro (Angular position) in third link: See mbeSensorAngle(is_fixed,idxs1,idxs2,noise_std)
            };
    end
   
    methods 
        % Constructor
        function [me] = mbeMechModelMaqueta(LA1,L12,L2B)
            % Initial conditions
            % M A Q U E T A
%             x1 = .120; y1= 0; x2 = .540; y2 = .338; x3 = .330; y3 = .323; theta = 0; beta = deg2rad(127.5); % dof->theta                            
            x1 = .016267; y1= 0.118892; x2 = .489447; y2 = .347018; x3 = .296403; y3 = .381854; theta = deg2rad(82.209); beta = deg2rad(130.2993); % dof->beta                            
            
            
            me.q_init_approx = [x1, y1, x2, y2, theta, x3, y3, beta]'; % Initial approximated position
            wA1k = 2700*(1/150)*2*pi/60; % 2700rpmMot * 1rpmA1/150rpmMOT * 2pirad/rev * 1m/60s ---> rad/s            
            %me.zp_init = -0.4309; % Initial DOF velocity
            me.zp_init = wA1k; % 
            % gravity
            me.g = -9.80665;

            % Multibody model parameters
            % fixed points:
            me.xA = 0; me.yA = 0; me.xB = 0.8; me.yB = 0;
            me.fixed_points = [me.xA, me.yA, me.xB, me.yB];
            % geometry
            if nargin < 3
                disp('default length values')
                LA1=0.126;
                L12 = 0.54;%L12=norm([x2,y2]-[x1,y1]); 
                L2B = 0.455; %L2B=norm([me.xB,me.yB]-[x2,y2]); 
                %L2B = 0.426; %L2B=norm([me.xB,me.yB]-[x2,y2]);
                %L13 = 0.38422;
                %L32 = 0.212191;
            end
            L13 = sqrt(2*(L12/2)^2);
            L32 = L13;
            me.bar_lengths = [LA1, L12, L2B, L13, L32];
            % mass
            me.mA1=1; me.m12 = 2; me.m2B = 4;
            %me.radio = .14; me.cdg = [0 0]';
            
            % Local Mass matrices. 
            %MA1 = me.massMatrixDisc(me.mA1,me.radio,me.cdg);
            MA1 = me.massMatrixDisc(me.mA1,.12,0.14,[0 0]);
            M12 = me.massMatrixBar(me.m12,L12);
            M2B = me.massMatrixBar(me.m2B,L2B);
            
            % Global Mass Matrix
            me.M = zeros(5,5);
            me.M(1:2,1:2) = MA1 (3:4,3:4);
            me.M(3:4,3:4) = M2B (1:2,1:2);
            me.M(1:4,1:4) = me.M(1:4,1:4) + M12;
          
            % Vector of generalized forces
            me=me.update_Qg();
            
            % damping coefficient
            me.C = 0;             
        end % end ctor()
        
        function [me] = update_Qg(me)
        % Re-calc Qg from gravity vector & masses.
            me.Qg = 0.5*me.g*[0; me.mA1+me.m12; 0; me.m12+me.m2B; 0];
        end            
        
    end
    
    % Implementation of Virtual methods from mbeMechModelBase
    methods(Access=public)
        % Evaluates potential energy of the whole system:
        function V = eval_potential_energy(me,q)
           V = -0.5*me.g*(me.mA1*(q(2)+me.yA)+me.m12*(q(2)+q(4))+me.m2B*(q(4)+me.yB)); 
        end
        
    end % methods
end % class
