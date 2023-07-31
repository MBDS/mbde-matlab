%% Example of off-line estimation

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
close all;% clear all;

addpath('toolbox_mbe');

%% 1) Select the estimator: 
% ---------------------------------
% estim = mbeEstimatorDEKF();
% estim = mbeEstimatorIncrManifold_full_jac(); % AKA errorEKF_EJ
% estim = mbeEstimatorDirect_eEKF_FJ(); % This is the direct version of the errorEKF_EJ
% estim = mbeEstimatorIncrManifold_full_jac();estim.dynamics_formulation_and_integrator = mbeDynFormulationMatrixR(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intTrapezoidal;
% estim = mbeEstimatorDirect_eEKF_FJ();estim.dynamics_formulation_and_integrator = mbeDynFormulationMatrixR(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intTrapezoidal;
% estim = mbeEstimatorIncrManifold(); % AKA errorEKF
% estim = mbeEstimatorIncrManifold(); estim.dynamics_formulation_and_integrator = mbeDynFormulationMatrixR(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intTrapezoidal;
% estim = mbeEstimatorIncrManifold_ace_full_jac(); % AKA errorEKF_FE
% estim = mbeEstimatorDirect_eEKF_FE(); % Direct version of errorEKF_FE
% estim = mbeEstimatorIncrManifold_AerrorEKF(); % AKA adaptive errorEKF_FE
% estim = mbeEstimatorIncrManifold_errorEKF_shaping(); % AKA errorEKF_FE combined with a shaping filter
% estim = mbeEstimatorIncrManifold_ace_full_jac(); estim.dynamics_formulation_and_integrator = mbeDynFormulationMatrixR(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intTrapezoidal;
% estim = mbeEstimatorDIEKFacc();
% estim = mbeEstimatorUKF();
% estim = mbeEstimatorUKF(); estim.dynamics_formulation_and_integrator.integrator_type = mbeIntegratorTypes.intEuler;
estim = mbeEstimatorCEKF();
% estim = mbeEstimatorDIEKF_pm(); 
% estim = mbeEstimatorBatchMRF();   
% estim = mbeEstimatorDEKFprojections();
%  estim = mbeEstimatorSCKF();
 


% Formulation for evaluating accelerations after each time step:
%estim.post_iter_accel_solver = mbeDynFormulationMatrixR(); 
%estim.post_iter_accel_solver = mbeDynFormulationLagrangeBaumgarte(); % Faster

%% 2) Pick multibody model:
% ---------------------------------
if (1) 
    enc_std_noise = deg2rad(1); 
    gyro_std_noise = 9.839439e-04; 
    acc_std_noise = 5.637583e-02;
%% mbeMechModelPendulum1    
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, [])}; % Accelerometer perpendicular. See mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, [])}; % Accelerometer parallel. See mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise*100,[2;0], pi/2, [])}; 
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[4;0], pi/2, [])}; 
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise)};
%     estim.mech_phys_model = mbeMechModelPendulum1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise)}; 

%% mbeMechModelFourBars1
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise)}; % Gyro in first link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
    estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise)}; % Gyro in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 1],[3 4],[3 4], gyro_std_noise)}; % Gyro in rocker link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(5, enc_std_noise)}; % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, [])};
% %%%% Other possible sensor configurations for the four-bar linkage %%%
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(5, enc_std_noise),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, [])}; % Encoder pos sensor: See mbeSensorPosIndex(q_index, std_dev_noise)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, [])}; % Accelerometer in first link, axis parallel to the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, [])}; estim.sensors_std_magnification4filter = 2; % Accelerometer in first link, axis perpendicular to the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[4;0], pi/2, [])}; % Accelerometer in first link, axis perpendicular to the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], 0, [])}; estim.sensors_std_magnification4filter = [20;4]; % Accelerometer in first link, two axes. see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([0 1],[1 2],[1 2], acc_std_noise,[0;0], pi/2, []),mbeSensorAccelerometer([0 1],[1 2],[1 2], acc_std_noise,[0;0], 0, [])}; % Accelerometer in first link, two axes. see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], -pi/4, [])}; % Accelerometer in first link, axis at 45 deg wrt the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)    
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2+pi/10, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], -pi/2-pi/10, [])}; % Accelerometer in first link, axis at 45 deg wrt the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)    
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[2;0], pi/2, []),mbeSensorAccelerometer([1 0],[1 2],[1 2], acc_std_noise,[3;0], pi/2, [])}; % Accelerometer in first link, axis perpendicular to the bar: see mbeSensorAccelerometer(pt_is_fixed,pt1_idxs,pt2_idxs,noise_std, LCC, axisangle, gravity)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([0 0],[1 2],[3 4], acc_std_noise,[0;0],0, []),mbeSensorAccelerometer([0 0],[1 2],[3 4], acc_std_noise,[0;0],pi/2, [])}; % Accelerometer in coupler link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise),mbeSensorAccelerometer([0 1],[3 4],[3 4], sen_std_noise,[0;0], pi/2, [])}; % Accelerometer in rocker link: See mbeSensorGyroscope(is_fixed,idxs1,idxs2,noise_std)
%     estim.mech_phys_model = mbeMechModelFourBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise),mbeSensorPosIndex(5, enc_std_noise)}; % Gyro in coupler and encoder in crank.

%% mbeMechModelFiveBars1
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(7, enc_std_noise) mbeSensorPosIndex(8, enc_std_noise)}; %Encoders in cranks
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise) mbeSensorGyroscope([0 1],[5 6],[3 4], gyro_std_noise)}; %Gyros in cranks
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise) mbeSensorGyroscope([0 0],[3 4],[5 6], gyro_std_noise)}; %Gyros in couplers
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([0 1],[1 2],[1 2], acc_std_noise,[0;0], pi/2, []), mbeSensorAccelerometer([0 1],[1 2],[1 2], acc_std_noise,[0;0], 0, []), mbeSensorAccelerometer([0 1],[5 6],[3 4], acc_std_noise,[0;0], pi/2, []),mbeSensorAccelerometer([0 1],[5 6],[3 4], acc_std_noise,[0;0], 0, [])}; estim.sensors_std_magnification4filter = [1;1;1;1]; %Accelerometers in cranks
    % %%%% Other possible sensor configurations for the five-bar linkage %%%
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise)}; %Gyros in first crank
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise) }; %Gyros in first coupler
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([0 0],[3 4],[5 6], gyro_std_noise)}; %Gyros in second coupler
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(7, enc_std_noise) mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise) }; %Encoder first crank, Gyros in first coupler
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorGyroscope([1 0],[1 2],[1 2], gyro_std_noise) mbeSensorGyroscope([0 0],[1 2],[3 4], gyro_std_noise)}; %Gyros in first crank, gyro in first coupler
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorPosIndex(7, enc_std_noise)}; %Enoder in first link
%     estim.mech_phys_model = mbeMechModelFiveBars1(); estim.mech_phys_model.installed_sensors =  {mbeSensorAccelerometer([0 1],[1 2],[1 2], acc_std_noise,[0;0], 0, []) mbeSensorAccelerometer([0 1],[5 6],[3 4], acc_std_noise,[0;0], 0, [])}; %Accelerometers in cranks


%% Bad model errors
    estim.bad_model_errors.error_type = [1,2]; % 1=>Gravity error; 2=>Initial pos; 6=>Bar length error; 7=> Mass error, 8=> constant torque error
    estim.bad_model_errors.error_scale = 1; %0.5*2*2;
         
%%  Estimator further params:
    multiple_multirate = 1; % if multiple_multirate == 0 => measurement every time step
    estim.dt       = 5e-3;
    estim.mechanism_type.multirate_sensor_period = estim.dt*multiple_multirate; 
    estim.end_time = 10; 

%% Four-bar noise 
    if isa(estim.mech_phys_model,'mbeMechModelFourBars1')
        if isa(estim,'mbeEstimatorCEKF')||isa(estim,'mbeEstimatorCEKF_TR')
             estim.transitionNoise_Zp = 0;
%              estim.transitionNoise_Zpp = (deg2rad(100))*0.035*1;%*0.4;
%                estim.transitionNoise_Zpp = (deg2rad(3.5))*200;
               estim.transitionNoise_Zpp = (deg2rad(100))*0.0035*15*8;
%                estim.transitionNoise_Zpp = (deg2rad(100))*0.0035*15/estim.dt;
%              estim.initVar_Z = 1.3e-4;
%              estim.initVar_Zp = 1.3e-2;
%              estim.initVar_Z = 0.00018225;
%              estim.initVar_Zp = 0.00297;
%         elseif isa(estim,'mbeEstimatorDIEKF_pm') || isa(estim,'mbeEstimatorDIEKF_pm_test')
%             estim.transitionNoise_Z = (deg2rad(1*100))*0;
%             estim.transitionNoise_Zp = (deg2rad(10*100)*0);
% %             estim.transitionNoise_Zpp = deg2rad(3.5*1.5)/40;
% %             estim.transitionNoise_Zpp = (deg2rad(3.5))*150*0.005;%Errorscale 1, encoders
%             estim.transitionNoise_Zpp = (deg2rad(100))*0.0035*15;
% %             estim.transitionNoise_Zpp = (deg2rad(100))*0.0035*15;%*0.4;
% %             estim.transitionNoise_Zpp = (deg2rad(3.5))*100;%Errorscale 1, encoders
%         elseif isa(estim,'mbeEstimatorSCKF') || isa(estim,'mbeEstimatorSCKF_test')
%             estim.transitionNoise_Z = (deg2rad(1))*0;
%             estim.transitionNoise_Zp = (deg2rad(10))*0;
%             estim.transitionNoise_Zpp = (deg2rad(3.5))*150;%*0.4;
        elseif isa(estim, 'mbeEstimatorIncrManifold_ace_full_jac')||isa(estim, 'mbeEstimatorDirect_eEKF_FE')
            estim.transitionNoise_Z = 0; 
            estim.transitionNoise_Zp = 0; 
            estim.transitionNoise_Zpp = deg2rad(3.5*1.5)/40;
        elseif isa(estim, 'mbeEstimatorIncrManifold_errorEKF_shaping')
            estim.transitionNoise_Z = 0; 
            estim.transitionNoise_Zp = 0; 
            estim.transitionNoise_Zpp = deg2rad(3.5*1.5)/40;
        elseif ~isa(estim, 'mbeEstimatorBatchMRF')
            estim.transitionNoise_Z = (deg2rad(1))*0.1*0; %*0.1*0.2;
    %         estim.transitionNoise_Zp = (deg2rad(100));
            estim.transitionNoise_Zp = (deg2rad(10))*0.000005*0.9*0; %*0.2; %*0.1*0.5;
            estim.transitionNoise_Zpp = (deg2rad(100))*0.0035*15;%*0.4;
       

        end
    end
%% Five-bar noise
    if isa(estim.mech_phys_model,'mbeMechModelFiveBars1')
        estim.transitionNoise_Z = 0; %*0.1*0.2;
        estim.transitionNoise_Zp = 0; %*0.2; %*0.1*0.5;
            if isa(estim,'mbeEstimatorCEKF')
                estim.transitionNoise_Zpp = (deg2rad(100))*0.035*65;%*0.4;
                estim.transitionNoise_Zpp = (deg2rad(100))*0.035*4*0.01;%*0.4;
                estim.transitionNoise_Zpp = (deg2rad(100))*0.035*25;%*0.4; % Noise with encoders
                estim.transitionNoise_Zpp = (deg2rad(100))*0.035*1.5;%*0.4; % Noise with gyros
                estim.initVar_Z = 1.3e-5;
                estim.initVar_Zp = 1.3e-3;
%                 estim.initVar_Z = 1.3e-4;
%                 estim.initVar_Zp = 1.3e-3;
%              estim.initVar_Z = 1.3e-4;
%              estim.initVar_Zp = 1.3e-2;

            elseif isa(estim,'mbeEstimatorSCKF')
                estim.transitionNoise_Zpp = (deg2rad(3.5))*250;%*0.4;
            elseif isa(estim, 'mbeEstimatorIncrManifold_ace_full_jac')
                estim.transitionNoise_Zpp = (deg2rad(3.5))*9/5; %/5e-3;%*0.4;
            else
%                 estim.transitionNoise_Zpp = (deg2rad(100))*0.035*15;%*0.4;
                estim.transitionNoise_Zpp = (deg2rad(3.5))*9;%*0.4;
            end      
        
    end
%% Pendulum noise
if isa(estim.mech_phys_model,'mbeMechModelPendulum1')
    if isa(estim,'mbeEstimatorCEKF')||isa(estim,'mbeEstimatorCEKF_TR')||isa(estim,'mbeEstimatorCEKF_SR')
% % %         estim.transitionNoise_Zp = 0;
% % %         %             estim.transitionNoise_Zpp = (deg2rad(10000*estim.dt))*0;
% % %         %             estim.transitionNoise_Zpp = (deg2rad(100))*sqrt(multiple_multirate)*0.035*1000;%*0.4;
% % %         estim.transitionNoise_Zpp = (deg2rad(100))*0.035*1;%*0.4;
% % %         estim.initVar_Z = 1.3e-4;
% % %         estim.initVar_Zp = 1.3e-2;
% % %         %              estim.initVar_Z = 0.00018225;
% % %         %              estim.initVar_Zp = 0.00297;
    elseif isa(estim,'mbeEstimatorDIEKF_pm')
% %         estim.transitionNoise_Z = (deg2rad(1))*0;
% %         estim.transitionNoise_Zp = (deg2rad(10))*0;
% %         estim.transitionNoise_Zpp = (deg2rad(3.5))*100;%Errorscale 1, encoders
% %         %             estim.transitionNoise_Zpp = (deg2rad(3.5))*100;%Errorscale 1, encoders
    elseif isa(estim,'mbeEstimatorSCKF')
%         estim.transitionNoise_Z = (deg2rad(1))*0;
%         estim.transitionNoise_Zp = (deg2rad(10))*0;
%         estim.transitionNoise_Zpp = (deg2rad(3.5))*150;%*0.4;
    elseif isa(estim, 'mbeEstimatorIncrManifold_ace')
%         estim.transitionNoise_Z = 0; %*0.1*0.2;
%         estim.transitionNoise_Zp = 0; %*0.2; %*0.1*0.5;
%         estim.transitionNoise_Zpp = (deg2rad(100))*0.035;%*0.4;
%         %             estim.initVar_Zpp = (deg2rad(30))^2;
%         %             estim.initVar_Z = (deg2rad(1))^2;
    elseif isa(estim, 'mbeEstimatorIncrManifold_ace_full_jac') || isa(estim, 'mbeEstimatorerrorUKFace') || isa(estim, 'mbeEstimatorUKFace')
%         estim.transitionNoise_Z = 0; %*0.1*0.2;
%         estim.transitionNoise_Zp = 0; %*0.2; %*0.1*0.5;
%         estim.transitionNoise_Zpp = (deg2rad(100))*0.035*0.01; %*1000;%*0.1*0.01;
%         estim.initVar_Zpp = (deg2rad(30))^2;
%         estim.initVar_Z = (deg2rad(10))^2;
    end
end
else
    %% Real mechanism. Load dataset as ground-truth:
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
%% Some other unused options

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

%% 4) Launch filter offline estimation:
% ------------------------------------
rng(10); % Fixes the seed for the pseudorandom noise (all the tests use the same pseudorandom sequence).
estim.mech_phys_model.zp_init = 0*ones(size(estim.mech_phys_model.zp_init)); % Set initial velocities
Qconst_val = 0; 
Qconst = zeros(estim.mech_phys_model.dep_coords_count,1);
Qconst(estim.mech_phys_model.get_indep_indxs) = Qconst_val; 
estim.mech_phys_model = estim.mech_phys_model.set_Qconst(Qconst); % Set constant forces (or torques) applied to the DOFs
estim.debug_show_live_filter_anim_fps = 10; % Frames per second for the animation
estim.run_offline(); 