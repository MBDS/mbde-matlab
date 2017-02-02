classdef (Abstract) mbeEstimatorFilterBase < mbeEstimatorBase
    %mbeEstimatorFilterBase The virtual base class for all FILTER estimators 
    % The entry point for estimation experiments, offline & online.
    
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

    % --- Params shared to ALL filter estimators ---
    % Public read-write params:
    properties (Access=public)
        % =========== Unified parameters for ALL state observers ===========
        % "Initial covariance"
        initVar_Z  = (deg2rad(5))^2;  % Initial variance for "z"  (deg)
        initVar_Zp = (deg2rad(5))^2;  % Initial variance for "z'" (deg/s)
        initVar_Zpp= (deg2rad(5))^2;  % Initial variance for "z''" (deg/s^2)

        % "Plant noise": Variance of the transition model for "z"
        % (In this form, it fits directly continuous-time filters. Must be
        % multiplied by "dt" for discrete-time filters):
        transitionNoise_Z   = (deg2rad( 1 ))^2; % 10  % (variance in deg, per second^(1/2))
        transitionNoise_Zp  = (deg2rad( 50 ))^2; % 100 %  (variance in deg/s, per second^(1/2))
        transitionNoise_Zpp = (deg2rad( 100 ))^2; %  (variance in deg/s^2, per second^(1/2))
        
        % Slower for debug: check for positive definit. of P at each iter.
        debug_always_check_def_pos = 1; 
        
        % (Slow) Show animations during filtering? 
        % Frequency in FPS (of filter time) (0=disable)
        debug_show_live_filter_anim_fps = 10;
        
        % Video
        video_filename = '';
        
        % Formulation for evaluating accelerations after each time step:
        post_iter_accel_solver = mbeDynFormulationMatrixR();
        
    end
    % Public read-only params:
    properties (GetAccess=public, SetAccess=protected)  
        % Online/offline filters may need to know the current "time",
        % in seconds since the beginning of the experiment:
        current_run_time = 0.0;
        % Like current_run_time but counts discrete timesteps.
        current_timestep = 0;
        
        % State vector (position part)
        q = [];
        % State vector (vel part)
        qp = [];
        % State vector (acc part)
        qpp = [];
        
        % Covariance matrix (we assume all filters have one, so it's here)
        P = [];                
        Innovation = [];
    end
    
    % --- Handy vars set-up by methods in this base class for the convenience of derived classes ---
    properties (Access=protected)
        lenQ;  % Number of coords: |q|
        lenZ;  % Number of indep coords: |z|
        iidxs; % Indices of independent coords in q()
        didxs; % Indices of dependent coords in q()
    end
    
    % --- Public API ---
    methods
        function [perfStats] = run_offline(me)
            % Execute one complete offline experiment
            % Optional live animations show 3 mechanism states: GT, Wrong-Model
            % and Estimator.
            
            % Init common vars to all filters:
            me.current_run_time = 0.0;
            me.current_timestep = 0;
            
            me.lenQ = me.mech_phys_model.dep_coords_count();
            me.lenZ = length( me.mech_phys_model.indep_idxs );
            me.iidxs = me.mech_phys_model.get_indep_indxs();
            me.didxs = me.mech_phys_model.get_dep_indxs();
            
            me.init_bad_model(); % Make sure "bad model" is up-to-date
            
            me.mechanism_type.resetState(); % Prepare for simulating from t=0
            
            % Initial state vector mean (from the "bad model") 
            % -------------------------------------------------
            % Make sure it is self-consistent for...
            % ... pos:
            me.q = mbeKinematicsSolver.pos_problem(...
                me.bad_mech_phys_model, ...
                me.bad_mech_phys_model.q_init_approx);
            
            % ... vel:
            me.qp = mbeKinematicsSolver.vel_problem(...
                me.bad_mech_phys_model,...
                me.q,...
                me.bad_mech_phys_model.zp_init);

            % ... initial acceleration (dynamics):
            me.qpp =  me.mechanism_type.dynamic_formulation.solve_for_accelerations(...
                me.bad_mech_phys_model,...
                me.q,...
                me.qp,0);
            
            % Initial covariance: done in each filter's init_filter()

            % Load data: 
            [states_GT, states_BadModel] = me.mechanism_type.generate_ground_truth(me);

            % Experiment duration:
            nTimeSteps = ceil(me.end_time/me.dt);
            assert(size(states_GT.q,1)==nTimeSteps,'Mismatch in GT length!')
            assert(size(states_BadModel.q,1)==nTimeSteps,'Mismatch in BadModel length!')

            % For logging progress:
            logstep_last = 0;
            logstep_nSteps = 20;
            logstep_incr = nTimeSteps/logstep_nSteps;
            
            % For debugging animations:
            liveanim_last = 0;
            if (me.debug_show_live_filter_anim_fps>0)
                liveanim_fig = figure();
%                 set(liveanim_fig,'Position',get( 0, 'ScreenSize' ));
                pause(0.001);
                liveanim_incr = (1.0/me.debug_show_live_filter_anim_fps)/me.dt;
                % Video settings
                if isempty(me.video_filename)
                    temp_name = class(me);
                    me.video_filename = [temp_name(13:end),'_'];
                    temp_name = class(me.mech_phys_model);
                    me.video_filename = [me.video_filename, temp_name(13:end),'_'];
                    temp_name = class(me.bad_mech_phys_model.installed_sensors{1});
                    me.video_filename = [me.video_filename, temp_name(10:end)];
                    clear('temp_name');
                end
%                 videowriter = VideoWriter(me.video_filename,'MPEG-4');
                videowriter = VideoWriter(me.video_filename);
                videowriter.FrameRate = me.debug_show_live_filter_anim_fps;
                open(videowriter);
            else
                liveanim_incr = Inf; % Disabled
            end
            
            % Init filter:
            me.init_filter();
            
            % Reserve stats for benchmark:
            hist=struct();
            if (me.do_benchmark)
                hist.t  = zeros(nTimeSteps,1);  % time of each timestep
                hist.estim_q   = zeros(nTimeSteps,length(me.q));
                hist.estim_qp  = zeros(nTimeSteps,length(me.q));
                hist.estim_qpp = zeros(nTimeSteps,length(me.q));
                hist.P         = cell(nTimeSteps,1); % Cov history
                [hist.Pidxs_z, hist.Pidxs_zp, hist.Pidxs_zpp ] = me.get_indices_in_P();
                hist.sensor_data = cell(nTimeSteps,1); % cell: Variable number of readings in each timestep!
                hist.estim_inn = cell(1,nTimeSteps);
            end


            % ===== Main loop ===== 
            for i = 1:nTimeSteps
                % Set iteration:
                me.current_timestep = i;
                me.current_run_time = me.current_timestep*me.dt;
                
                % Logging:
                if (i>=logstep_last+logstep_incr || i==nTimeSteps)
                    pcdone = 100.0 * i / nTimeSteps;
                    me.log_level(1, sprintf('run_offline: %.02f%% done. iter %u/%u, simul_time=%.03fs/%.03fs',...
                        pcdone,...
                        me.current_timestep,nTimeSteps,...
                        me.current_run_time,me.end_time));
                    logstep_last = i;
                end
                
                % ----------------
                % (1/2) Get observation(s), if ready:
                % ----------------
                obs = me.mechanism_type.get_current_observations( ...
                    me, ...
                    states_GT.q(me.current_timestep,:),...
                    states_GT.qp(me.current_timestep,:),...
                    states_GT.qpp(me.current_timestep,:) );
                
                % ----------------
                % (2/2) Run filter:
                % ----------------
                me.run_filter_iter(obs);
                
                % Sanity check:
                if (me.debug_always_check_def_pos)
                    [~, def_pos]= chol(me.P); % Checking if P is possitive definite
                    if def_pos > 0
                        warning ('P is not positive definite!')
                        pause
                    end
                end
                
                % Debug animations: GT, BadModel, Estim
                if (i>=liveanim_last+liveanim_incr)
                    clf; hold on;
                    % Plot GT, BadModel, Estim
                    me.bad_mech_phys_model.plot_model_skeleton(states_GT.q(i,:), 'r', 0);
                    me.bad_mech_phys_model.plot_model_skeleton(states_BadModel.q(i,:), 'g:', 0);
                    me.bad_mech_phys_model.plot_model_skeleton(me.q, 'b', 1); % the last 1=fit plot axes
                    % Title:
                    title('r=GT; g:=BadModel; b=Estim');
                    drawnow;
%                     frame = getframe;
                    writeVideo(videowriter,getframe);
                    liveanim_last=i;
                end

                % Benchmark stats:
                if (me.do_benchmark)
                    hist.t(i) = me.current_run_time;
                    hist.estim_q(i,:)   = me.q;
                    hist.estim_qp(i,:)  = me.qp;
                    hist.estim_qpp(i,:) = me.qpp;
                    hist.P{i} = me.P;
                    hist.sensor_data{i} = obs;
                    hist.estim_inn{i} = me.Innovation;
                end
                
            end % =====  end of for each time step ===== 
            % close video logging
            if (me.debug_show_live_filter_anim_fps>0)
                close(videowriter);
            end
            % Plots, stats...
            perfStats = struct();
            if (me.do_benchmark)
                perfStats.hist = hist;                
                perfStats.states_GT = states_GT;
                perfStats.states_BadModel = states_BadModel;                
                perfStats.stats = me.plot_estim_performance(states_GT, states_BadModel,hist, me.show_benchmark_graphs);
            end

        end % end run_offline() 
        
    end

    % Private stuff: Abstract methods:
    methods(Abstract,Access=protected)
        % Run one timestep of the estimator:
        % - o [In]: Vector with "observations" (variable length, may be empty)
        [] = run_filter_iter(me, o);
        
        % Initialize the filter state vector and covariance matrices as
        % required:
        [] = init_filter(me);

        % Returns the location (indices) in P for the independent
        % coordinates (z,zp,zpp). Zeros means not estimated by this filter.
        [P_idx_z, P_idx_zp, P_idx_zpp] = get_indices_in_P(me);
    end
    
end

