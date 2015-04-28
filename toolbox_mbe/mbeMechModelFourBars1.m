classdef mbeMechModelFourBars1 < mbeMechModelFourBarsBase
    % Mechanism physical model: 4 bars linkage ("Ferrol model, #1")
    %  See mbeMechModelFourBarsBase for a sketch of the mechanism.
    
    % (Abstract) Read-Write properties
    properties(Access=public)
        % List of installed sensors (cell of objects derived from mbeSensorBase)
        % List of installed sensors (cell of objects derived from mbeSensorBase)
        installed_sensors = { ...
%            	mbeSensorPosIndex(5, deg2rad(1)) ...  % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%                mbeSensorAccelerometer3pt([1 0],[1 2],[1 2], [1 2], 0.1 )... % Accelerometer with 3 points
%            	mbeSensorAngle([1 0],[3 4],[3 4], deg2rad(1))...  % Encoder in an arbitrary bar: See mbeSensorAngle(is_fixed,idxs1,idxs2,noise_std)
%            	mbeSensorAngle([0 0],[1 2],[3 4], deg2rad(1))...  % Encoder in an arbitrary bar: See mbeSensorAngle(is_fixed,idxs1,idxs2,noise_std)
%                mbeSensorGyroscope([1 0],[1 2],[1 2], deg2rad(1))... % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                mbeSensorGyroscope([0 0],[1 2],[3 4], deg2rad(1))... % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
                %mbeSensorGyroscope([0 1],[3 4],[3 4], deg2rad(1)) % Gyro in third link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
            };
    end
   
    methods 
        % Constructor
        function [me] = mbeMechModelFourBars1()
            % Initial conditions
            % F E R R O L
            x1 = 0; y1= 2; x2 = 10; y2 = 5; %
            theta = pi/3; % Initial DOF value
            me.q_init_aprox = [x1, y1, x2, y2, theta]'; % Initial approximated position
            me.zp_init = 0; % Initial DOF velocity
            % gravity
            me.g = -10;

            % Multibody model parameters
            % fixed points:
            me.xA = 0; me.yA = 0; me.xB = 10; me.yB = 0;
            me.fixed_points = [me.xA, me.yA, me.xB, me.yB];
            % geometry
            LA1=2; %norm([me.xA,me.yA]-[x1,y1]); 
            L12=8; %norm([x2,y2]-[x1,y1]); 
            L2B=5; %norm([me.xB,me.yB]-[x2,y2]); 
            me.bar_lengths = [LA1, L12, L2B];
            
            % mass
            me.mA1 = 2; 
            me.m12 = 8; 
            me.m2B = 5;

            % Local Mass matrices. 
            MA1 = me.massMatrixBar(me.mA1,LA1);
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
