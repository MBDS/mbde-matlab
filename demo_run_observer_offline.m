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

clear; close all; clear all;

addpath('toolbox_mbe');
% 1) Select the estimator: 
% ---------------------------------
estim = mbeEstimatorDEKF(); % DEKF
% estim = mbeEstimatorIncrManifold(); % AKA errorEKF: by default, I3AL MB formulation with trapezoidal rule integration
% estim = mbeEstimatorIncrManifold(); estim.dynamics_formulation_and_integrator = mbeDynFormulationMatrixR();estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intTrapezoidal; % errorEKF with matrix R formulation and trapezoidal rule intergation
% estim = mbeEstimatorDIEKFacc();
% estim = mbeEstimatorUKF(); % UKF with trapezoidal rule integration
% estim = mbeEstimatorUKF(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intEuler; % UKF with forward Euler integrator
% estim = mbeEstimatorCEKF(); % CEKF with forward Euler integration for covariance matrix of the estimation error
% estim = mbeEstimatorCEKF_TR(); % CEKF with trapezoidal rule integration for covariance matrix of the estimation error
% estim = mbeEstimatorDIEKF_pm(); % Method with "perfect measurements"
% estim = mbeEstimatorSCKF(); % Smoothly constrained Kalman filter

% estim = mbeEstimatorBatchMRF(); % Not an observer, but a batch smoother

% Formulation for evaluating accelerations after each time step:
%estim.post_iter_accel_solver = mbeDynFormulationMatrixR(); 
%estim.post_iter_accel_solver = mbeDynFormulationLagrangeBaumgarte(); % Faster

% 2) Pick multibody model and set of sensors:
% ---------------------------------
if (1) % Normal behaviour is achieved with this condition set to 1. If this condition is chaged to 0, a file with experimental data must be provided
    sen_noise = deg2rad(1);
% Four-bar linkage:
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], sen_noise)}; % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], sen_noise)}; % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 1],[3 4],[3 4], sen_noise)}; % Gyro in rocker link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
    estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(5, sen_noise)}; % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
    
% Five-bar linkage    
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], sen_noise),mbeSensorPosIndex(5, sen_noise)}; % Gyro in coupler and encoder in crank.
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(7, sen_noise) mbeSensorPosIndex(8, sen_noise)}; %Encoders in cranks
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], sen_noise) mbeSensorGyroscope([0 1],[5 6],[3 4], sen_noise)}; %Gyros in cranks
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], sen_noise) mbeSensorGyroscope([0 0],[3 4],[5 6], sen_noise)}; %Gyros in couplers
    
% 3) Select kind of error in modeling: 
% -------------------------------------
    estim.bad_model_errors.error_type = [1,2]; % 1=>Gravity error; 2=>Initial position error;  
    estim.bad_model_errors.error_scale = 1; % Scale of the modeling error
    
        
    % Estimator further params:
    multiple_multirate = 1; % if multiple_multirate equals 0 or 1 => measurement every time step, if multiple multirate ==2 => measurements every 2 time steps, etc.
    estim.dt       = 5e-3; % Time step of the simulations
    estim.mechanism_type.multirate_sensor_period = estim.dt*multiple_multirate; 
    estim.end_time = 10; % Simulation time
    % Set initial covariance matrix and plant covariance matrix (Some of
    % these values might have to be changed depending on the set of
    % sensors, the mechanism, the modeling error, the observer formulation, etc).
    estim.transitionNoise_Z = 0; 
    estim.transitionNoise_Zp = 0;
    if isa(estim.mech_phys_model,'mbeMechModelFourBars1')
        if isa(estim,'mbeEstimatorCEKF')||isa(estim,'mbeEstimatorCEKF_TR')
            estim.transitionNoise_Zpp = deg2rad(3.5);
            estim.initVar_Z = 1.3e-4;
            estim.initVar_Zp = 1.3e-2;
        elseif isa(estim,'mbeEstimatorDIEKF_pm')
            estim.transitionNoise_Zpp = deg2rad(350);
        elseif isa(estim,'mbeEstimatorSCKF')
            estim.transitionNoise_Zpp = deg2rad(525);
        elseif ~isa(estim, 'mbeEstimatorBatchMRF')
            estim.transitionNoise_Zpp = deg2rad(5.25);
        end
    end
    
    if isa(estim.mech_phys_model,'mbeMechModelFiveBars1')
        if estim.bad_model_errors.error_scale==1
            if isa(estim,'mbeEstimatorCEKF')||isa(estim,'mbeEstimatorCEKF_TR')
                % Select one of the following if using CEKF or CEKF_TR
                estim.transitionNoise_Zpp = deg2rad(87.5); % Noise with encoders
%                 estim.transitionNoise_Zpp = deg2rad(5.25); % Noise with gyros
                estim.initVar_Z = 1.3e-5;
                estim.initVar_Zp = 1.3e-3;
            elseif isa(estim,'mbeEstimatorSCKF')
                estim.transitionNoise_Zpp = deg2rad(875);
            else
                estim.transitionNoise_Zpp = deg2rad(31.5);
            end
        elseif estim.bad_model_errors.error_scale==0.5
            if isa(estim,'mbeEstimatorCEKF')||isa(estim,'mbeEstimatorCEKF_TR')
                % Select one of the following if using CEKF or CEKF_TR
                estim.transitionNoise_Zpp = deg2rad(35); % Noise wiht encoders
%                 estim.transitionNoise_Zpp = deg2rad(3.5);% Noise with gyros
                estim.initVar_Z = 1.3e-5;
                estim.initVar_Zp = 1.3e-3;
            elseif isa(estim,'mbeEstimatorSCKF')
                estim.transitionNoise_Zpp = deg2rad(35);
            else
                estim.transitionNoise_Zpp = deg2rad(17.5);
            end
        else
            warning('Noises have not been adjusted for estim.bad_model_errors.error_scale values othr than 0.5 or 1')
            pause; 
        end
        
        
    end
   
else
    % Real mechanism. Load dataset as ground-truth:
    % estim.dt       = 20e-3;
    estim.end_time = 20;
    
    estim.mech_phys_model = mbeMechModelFourBars3(); % es la maqueta pero sin el punto 3 ni el segundo angulo
    estim.mechanism_type = mbeMechTypeReal();
    estim.mechanism_type.dataset = 'LOG_2015_02_18_17_06_14_';

    %estim.transitionNoise_Zp = (deg2rad( 500 ))^2;
    %estim.transitionNoise_Zpp = (deg2rad( 10000 ))^2;
    estim.SIGMAp = deg2rad(1)*0;   % Transition noise in "z" (assumed by the estimator)
    estim.SIGMAr = deg2rad(.1)*0;  % Transition noise in "zp" (assumed by the estimator)
    estim.SIGMAa = deg2rad(100);  % Transition noise in "zpp" (assumed by the estimator)
end


% Select dynamics formulation for simulating GT/Bad model:
estim.mechanism_type.dynamic_formulation = mbeDynFormulationMatrixR();

% 4) Launch filter offline estimation:
% ------------------------------------
rng(10);
estim.run_offline(); 


