classdef mbeSensorBase
    % Virtual sensor virtual base class. See docs of derived classes for
    % examples.
    % TODO: Make sensors capable of outputing more than one scalar value.
    
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

