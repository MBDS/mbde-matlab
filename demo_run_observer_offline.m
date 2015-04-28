% Example of off-line estimation

clear;  clc, close all 
%clear classes; %clear class;
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
    %estim.mech_phys_model = mbeMechModelFourBars1();
    estim.mech_phys_model = mbeMechModelFiveBars1();
    %estim.mech_phys_model = mbeMechModelPendulum1();

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

