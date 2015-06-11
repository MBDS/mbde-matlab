% Example of off-line estimation

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

clear;
close all; 

addpath('toolbox_mbe');

% 1) Select the estimator: 
% ---------------------------------
estim = mbeEstimatorDEKF();
% estim = mbeEstimatorIncrManifold(); 
% estim = mbeEstimatorDIEKFacc();
% estim = mbeEstimatorUKF();
% estim = mbeEstimatorCEKF();
% estim = mbeEstimatorDIEKF_pm(); 
% estim = mbeEstimatorBatchMRF();

% Formulation for evaluating accelerations after each time step:
%estim.post_iter_accel_solver = mbeDynFormulationMatrixR(); 
%estim.post_iter_accel_solver = mbeDynFormulationLagrangeBaumgarte(); % Faster

% 2) Pick multibody model:
% ---------------------------------
if (1) 
    % Simulate mechanism:
    estim.mech_phys_model = mbeMechModelFourBars1();
%     estim.mech_phys_model = mbeMechModelFiveBars1();
%     estim.mech_phys_model = mbeMechModelPendulum1();

estim.bad_model_errors.error_type = 2;

estim.mech_phys_model.installed_sensors = {mbeSensorGyroscope([1 0],[1 2],[1 2], deg2rad(1))}; % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%  estim.mech_phys_model.installed_sensors = {mbeSensorGyroscope([0 0],[1 2],[3 4], deg2rad(1))}; % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)    
% estim.mech_phys_model.installed_sensors = {mbeSensorGyroscope([0 1],[3 4],[3 4], deg2rad(1))}; % Gyro in third link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)    

    % Estimator further params:
    estim.dt       = 5e-3;
    estim.mechanism_type.multirate_sensor_period = 10*estim.dt*0; % if == 0 => measurement every time step
    estim.end_time = 10.0; 
   
else
    % Real mechanism. Load dataset as ground-truth:
    % estim.dt       = 20e-3;
    estim.end_time = 20;
    
    estim.mech_phys_model = mbeMechModelFourBars3(); % es la maqueta pero sin el punto 3 ni el segundo angulo
    estim.mechanism_type = mbeMechTypeReal();
    estim.mechanism_type.dataset = 'LOG_2015_02_18_17_06_14_';

    %estim.transitionNoise_Zp = (deg2rad( 500 ))^2;
    %estim.transitionNoise_Zpp = (deg2rad( 10000 ))^2;
    estim.SIGMAp = deg2rad(1);   % Transition noise in "z" (assumed by the estimator)
    estim.SIGMAr = deg2rad(.1);  % Transition noise in "zp" (assumed by the estimator)
    estim.SIGMAa = deg2rad(100);  % Transition noise in "zpp" (assumed by the estimator)
end


% Select dynamics formulation for simulating GT/Bad model:
%estim.mechanism_type.dynamic_formulation = mbeDynFormulationMatrixR();

% Enable simulating "multirate sensors": set this to the period
% between sensor data (e.g. 1.0/f for f in Hz). Leave to 0
% (default) to simulate sensor data at every time step:
%estim.mechanism_type.multirate_sensor_period = 1.0/ 5;


% Change default installed sensors:
% List of installed sensors (cell of objects derived from mbeSensorBase)
%estim.mech_phys_model.installed_sensors = { ...
%        mbeSensorPosIndex(5, deg2rad(0.5)), ...  % mbeSensorPosIndex(qp_index, std_dev_noise)
%    };

% 3) Select kind of error in modeling: 
% -------------------------------------
%estim.bad_model_errors.error_type = 2; 
%estim.bad_model_errors.error_scale = 1; 

% 4) Launch filter offline estimation:
% -------------------------------------

estim.run_offline();

