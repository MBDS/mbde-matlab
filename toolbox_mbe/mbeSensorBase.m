classdef mbeSensorBase
    % Virtual sensor virtual base class. See docs of derived classes for
    % examples.
    % TODO: Make sensors capable of outputing more than one scalar value.
    
    methods(Abstract)
        % --- virtual methods implementation  ----
        % The (scalar) value of this sensor output for the given dynamical state.
        [obs] = sensor_simulate(me,model,q,qp,qpp);

        % The standard deviation of this sensor
        [sensor_std] = sensor_std_noise(me);
        
        % Evaluates the sensor model jacobian for a given dyn state:
        %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
        % H= [  ...   |   ...   |    ...   ] = [ dh_dq , dh_dqp, dh_dqpp ]
        %    [ dh1_dq | dh1_dqp | dh1_dqpp ]
        [dh_dq , dh_dqp, dh_dqpp] = sensor_jacob(me,model,q,qp,qpp);
    end
end

