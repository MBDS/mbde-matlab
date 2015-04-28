classdef mbeSensorPosIndex < mbeSensorBase
    % Sensor: direct sensor for one of the values in "q" (e.g. encoder
    % angle, linear position, etc.)
    
    properties
        q_idx;    % Index in q
        std_dev;  % Noise stddev
    end
    
    methods
        % Ctor
        function [me]=mbeSensorPosIndex(q_idx,noise_std)
            me.q_idx = q_idx;
            me.std_dev = noise_std;
        end

        
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        function [obs] = sensor_simulate(me,model,q,qp,qpp)
            obs = q(me.q_idx);
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
            dh_dq=zeros(1,n);  dh_dq(me.q_idx)=1;
            dh_dqp=zeros(1,n);
            dh_dqpp=zeros(1,n);
        end
        
        
        
    end
    
end

