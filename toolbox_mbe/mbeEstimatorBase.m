classdef (Abstract) mbeEstimatorBase < handle
    %MBEESTIMATORBASE The virtual base class for all estimators 
    % The entry point for estimation experiments, offline & online.
    % Includes a generic method plot_estim_performance() to display 
    % performance graphs from experiments.
    
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

    % --- Params shared to ALL estimators ---
    properties(GetAccess=public, SetAccess=public)
        % Whether to show messages, animations, etc.
        logging = mbeLoggingParams();
        
        % Simulation end time
        end_time = 5.0; 
        
        % time step
        dt = 5e-3;
        
        % The type of mechanism: class derived from mbeMechTypeBase
        mechanism_type = mbeMechTypeSimulated();  
        
        % Object with the physical properties of the mechanism: lengths, masses, etc.
        % This is the "ground-truth" model. A "corrupted" version is
        % automatically generated and stored in bad_mech_phys_model 
        % with the errors given in bad_model_errors
        mech_phys_model = mbeMechModelFourBars1();
        
        % Definition of the kind of errors to introduce in the MB model.
        bad_model_errors = mbeModelErrorDef();
        
        % Factor to multiply std deviations of sensors in estimators' noise model.
        sensors_std_magnification4filter = 1.0; 
        
        % plot_estim_performance(): Confidence interval for plots (default=0.95)
        perf_eval_CI = 0.95;
        
        % Generate performance stats (they are returned in `perfStats` in `run_offline()`);
        do_benchmark = 1;

        % If 'do_benchmark=1', whether to also show all those results
        % graphically.
        show_benchmark_graphs = 1;        
    end

    properties(GetAccess=public, SetAccess=protected)
        % Filled-in by init_bad_model() from mech_phys_model & bad_model_errors
        bad_mech_phys_model;
    end
    
    % --- Public API ---
    methods(Abstract)
        % Run one entire estimation experiment
        [perfStats] = run_offline(me);
    end
    
    methods
        % Display a 'msg' if user's verbose_level>=VERBOSE_LEVEL
        function []=log_level(me, VERBOSE_LEVEL, msg)
            me.logging.log_level(VERBOSE_LEVEL, msg);
        end
        
        function [] = init_bad_model(me)
            % (re)initialize wrong model from the GT model + the errors in "bad_model_errors"
            me.bad_mech_phys_model = me.mech_phys_model.applyErrors(me.bad_model_errors);
        end
        
        function [perfResults] = plot_estim_performance(me,states_GT, states_BadModel,hist,show_benchmark_graphs)
            % Plot estimation performance stats:
            % - states_GT, states_BadModel: structs with fields t,q,qp,qpp
            % - hist: struct with fields: t,estim_q,estim_qp,estim_qpp, P,...
            % (based on old code Kalman_plots.m)
            iidxs = me.mech_phys_model.get_indep_indxs();
            L = length(hist.estim_q);
            lenZ = length(iidxs); % =|z|=|zp|=|zpp|
            
            % Compute quantile for the desired percentile:
            % "r" is the number of dimensions (degrees of freedom)
            nSigmas_1_dof = sqrt(mbeUtils.qchisq(me.perf_eval_CI, 1 )); 
            nSigmas_lenZ_dofs = sqrt(mbeUtils.qchisq(me.perf_eval_CI, lenZ )); 
            
            aux_lgnds_real = repmat({'Real'},[1,lenZ]);
            aux_lgnds_model = repmat({'Model'},[1,lenZ]);
            aux_lgnds_observer = repmat({'Observer'},[1,lenZ]);
            
            x_axis = hist.t';
            % Mahalanobis distance calculation
            has_z_covs = (~isempty(hist.P) && min(hist.Pidxs_z)~=0); % If there is P, AND it contains z
            has_zp_covs = (~isempty(hist.P) && min(hist.Pidxs_zp)~=0); % If there is P, AND it contains zp
            has_zpp_covs = (~isempty(hist.P) && min(hist.Pidxs_zpp)~=0); % If there is P, AND it contains zpp
            % Eval the complete mahalanobis distance at each time step:
            % Build list of indices in P
            idxs=[hist.Pidxs_z, hist.Pidxs_zp, hist.Pidxs_zpp];
            idxs(idxs==0) = []; % Remove non-estimated vars
            
            % Build equivalent indices in q,qp, qpp:
            idxs_in_q=[]; idxs_in_qp=[]; idxs_in_qpp=[];
            if (has_z_covs)
                assert(length(iidxs)==length(hist.Pidxs_z));
                idxs_in_q=iidxs;
            end
            if (has_zp_covs)
                assert(length(iidxs)==length(hist.Pidxs_zp));
                idxs_in_qp=iidxs;
            end
            if (has_zpp_covs)
                assert(length(iidxs)==length(hist.Pidxs_zpp));
                idxs_in_qpp=iidxs;
            end
            
            maha_dist=zeros(L,1);  % Total maha dist.
            maha_dist_pos=zeros(L,1);  % maha dist. in position
            maha_dist_vel=zeros(L,1);  % maha dist. in velocity
            for i=1:L,
                % Build complete P & error vector for [z, zp, zpp],
                % excepting non-estimated variables:
                
                ith_P = hist.P{i}(idxs,idxs);
                ith_estim_state = [hist.estim_q(i,idxs_in_q), hist.estim_qp(i,idxs_in_qp), hist.estim_qpp(i,idxs_in_qpp)];
                ith_GT_state    = [states_GT.q(i,idxs_in_q), states_GT.qp(i,idxs_in_qp), states_GT.qpp(i,idxs_in_qpp)];
                
                error = ith_GT_state - ith_estim_state;
                maha_dist(i,1) = sqrt(error * (ith_P\error'));
                
                % z-only part:
                if (has_z_covs)
                    err_z  = states_GT.q(i,idxs_in_q) - hist.estim_q(i,idxs_in_q);
                    ith_Pz = hist.P{i}(hist.Pidxs_z,hist.Pidxs_z);
                    maha_dist_pos(i,1) = sqrt(err_z * (ith_Pz\err_z'));
                end
                % zp-only part:
                if (has_zp_covs)
                    err_zp  = states_GT.qp(i,idxs_in_qp)- hist.estim_qp(i,idxs_in_qp);
                    ith_Pzp = hist.P{i}(hist.Pidxs_zp,hist.Pidxs_zp);
                    maha_dist_vel(i,1) = sqrt(err_zp * (ith_Pzp\err_zp'));
                end
            end

        if (show_benchmark_graphs)
            figure();
            clf;
            hold on
                % ----------------------------------------------
                title('Independent coordinate position')
                % ----------------------------------------------
                grid on
                
                lgnds={};
                if (has_z_covs)
                    std_th=zeros(L,lenZ)*1e-3;
                    for i=1:L,
                        P_diag = diag(hist.P{i});
                        std_th(i,:) = sqrt( P_diag(hist.Pidxs_z) );
                    end

                    for i=1:lenZ,
                        fill([x_axis fliplr(x_axis)],...
                            [hist.estim_q(:,iidxs(i))+nSigmas_1_dof*std_th(:,i);...
                            flipud(hist.estim_q(:,iidxs(i))-nSigmas_1_dof*std_th(:,i))]', 'c', 'EdgeColor', 'c');
                        lgnds={lgnds{:}, sprintf('%.03f%%CI z(%u)',me.perf_eval_CI*100,i)};
                    end
                end
                plot(hist.t,states_GT.q(:,iidxs),'r--')
                plot(hist.t,states_BadModel.q(:,iidxs),'g')
                plot(hist.t,hist.estim_q(:,iidxs),'b')
                legend({lgnds{:}, aux_lgnds_real{:}, aux_lgnds_model{:}, aux_lgnds_observer{:} } );
            hold off

            % Separate error plots:
            % ----------------------------------------------
            if (has_z_covs)
                for i=1:lenZ,
                    figure();grid on; hold on;
                    title(sprintf('Error in z_%u with %.03f%%CI ',i,me.perf_eval_CI*100));                    
                    fill([x_axis fliplr(x_axis)], [nSigmas_lenZ_dofs*std_th(1:end,i); -nSigmas_lenZ_dofs*std_th(end:-1:1,i)]', 'c', 'EdgeColor', 'c');
                    plot(hist.t,hist.estim_q(:,iidxs(i))-states_GT.q(:,iidxs(i)),'k')
                end
            end          

            
            figure;
            clf;
            hold on;
                % ----------------------------------------------
                title('Independent coordinate speed')
                % ----------------------------------------------
                grid on

                
                lgnds={};
                if (has_zp_covs)
                    std_thp=ones(L,lenZ)*1e-3;
                    for i=1:L,
                        P_diag = diag(hist.P{i});
                        std_thp(i,:) = sqrt( P_diag(hist.Pidxs_zp) );
                    end

                    for i=1:lenZ,
                        fill([x_axis fliplr(x_axis)],...
                            [hist.estim_qp(:,iidxs(i))+nSigmas_1_dof*std_thp(:,i); ...
                            flipud(hist.estim_qp(:,iidxs(i))-nSigmas_1_dof*std_thp(:,i))]', 'c', 'EdgeColor', 'c');
                        lgnds={lgnds{:}, sprintf('%.03f%%CI zp(%u)',me.perf_eval_CI*100,i)};
                    end
                end
                
                plot(hist.t,states_GT.qp(:,iidxs),'r')
                plot(hist.t,states_BadModel.qp(:,iidxs),'g')
                plot(hist.t,hist.estim_qp(:,iidxs),'b')
                legend({lgnds{:}, aux_lgnds_real{:}, aux_lgnds_model{:}, aux_lgnds_observer{:} } );
            hold off
            % Separate error plots:
            % ----------------------------------------------
            if (has_zp_covs)
                for i=1:lenZ,
                    figure();grid on; hold on;
                    title(sprintf('Error in zp_%u with %.03f%%CI ',i,me.perf_eval_CI*100));                    
                    fill([x_axis fliplr(x_axis)], [nSigmas_lenZ_dofs*std_thp(1:end,i); -nSigmas_lenZ_dofs*std_thp(end:-1:1,i)]', 'c', 'EdgeColor', 'c');
                    plot(hist.t,hist.estim_qp(:,iidxs(i))-states_GT.qp(:,iidxs(i)),'k')
                end
            end          
            

            figure;
            clf;
            hold on;
                % ----------------------------------------------
                title('Independent coordinate accelerations')
                % ----------------------------------------------
                grid on

                
                lgnds={};
                if (has_zpp_covs)
                    std_thpp=ones(L,lenZ)*1e-3;
                    for i=1:L,
                        P_diag = diag(hist.P{i});
                        std_thpp(i,:) = sqrt( P_diag(hist.Pidxs_zpp) );
                    end

                    for i=1:lenZ,
                        fill([x_axis fliplr(x_axis)],...
                            [hist.estim_qpp(:,iidxs(i))+nSigmas_1_dof*std_thpp(:,i); ...
                            flipud(hist.estim_qpp(:,iidxs(i))-nSigmas_1_dof*std_thpp(:,i))]', 'c', 'EdgeColor', 'c');
                        lgnds={lgnds{:}, sprintf('%.03f%%CI zpp(%u)',me.perf_eval_CI*100,i)};
                    end
                end
                
                plot(hist.t,states_GT.qpp(:,iidxs),'r')
                plot(hist.t,states_BadModel.qpp(:,iidxs),'g')
                plot(hist.t,hist.estim_qpp(:,iidxs),'b')
                legend({lgnds{:}, aux_lgnds_real{:}, aux_lgnds_model{:}, aux_lgnds_observer{:} } );
            hold off
            % Separate error plots:
            % ----------------------------------------------
            if (has_zpp_covs)
                for i=1:lenZ,
                    figure();grid on; hold on;
                    title(sprintf('Error in zpp_%u with %.03f%%CI ',i,me.perf_eval_CI*100));                    
                    fill([x_axis fliplr(x_axis)], [nSigmas_lenZ_dofs*std_thpp(1:end,i); -nSigmas_lenZ_dofs*std_thpp(end:-1:1,i)]', 'c', 'EdgeColor', 'c');
                    plot(hist.t,hist.estim_qpp(:,iidxs(i))-states_GT.qpp(:,iidxs(i)),'k')
                end
            end          
            
            
            figure()
            clf
            hold on
                % ----------------------------------------------
                title('Ground truth-predicted sensor(s) readings')
                % ----------------------------------------------
                % Keep only the timesteps with observations: (multirate, etc.)
                time_obs = []; % Yes, these vectores grow with the loop, but it's easier like this and we don't worry too much on performance in this method.
                values_obs = [];
                for i=1:L,
                    if (~isempty(hist.sensor_data{i}))
                        time_obs=[time_obs, hist.t(i)];
                        values_obs=[values_obs;hist.sensor_data{i}(:)'];
                    end
                end

                plot(time_obs,values_obs,'.',time_obs,values_obs,'-');
                grid on
            hold off

            figure()
            clf
            hold on
                % ----------------------------------------------
                title('Mahalanobis distance')
                % ----------------------------------------------
                grid on
                if (has_z_covs || has_zp_covs || has_zpp_covs)
                                  
                    nSigmas_lenTotal_dofs = sqrt(mbeUtils.qchisq(me.perf_eval_CI, length(idxs) ));

                    subplot(2,2,[1 2]); hold on;
                    plot(hist.t,maha_dist,'b');
                    plot([hist.t(1),hist.t(end)],[1,1]*nSigmas_lenTotal_dofs^2,'k'); % 3 sigmas
                    legend('MahaDist',sprintf('%.03f%% CI level',me.perf_eval_CI*100));
                    axis([hist.t(1) hist.t(end) 0 10]);  hold off;

                    subplot(2,2,3); hold on;
                    plot(hist.t,maha_dist_pos,'r');
                    plot([hist.t(1) hist.t(end)],[1,1]*nSigmas_lenZ_dofs^2,'k'); % 3 sigmas
                    legend('MahaDist(z only)',sprintf('%.03f%% CI level',me.perf_eval_CI*100));
                    axis([hist.t(1) hist.t(end) 0 10]);  hold off;
                    
                    subplot(2,2,4); hold on;
                    plot(hist.t,maha_dist_vel,'r');
                    plot([hist.t(1),hist.t(end)],[1,1]*nSigmas_lenZ_dofs^2,'k'); % 3 sigmas
                    legend('MahaDist(zp only)',sprintf('%.03f%% CI level',me.perf_eval_CI*100));
                    axis([hist.t(1) hist.t(end) 0 10]);  hold off;
                else
                    title('Mahalanobis distance: Cannot eval, no covariance P available!')
                end
            hold off

        end % if show_benchmark_graphs

            % ----------------------------------------------
            % GT, BadModel, Kalman: RMS errors
            % ----------------------------------------------
            position_error = sqrt( mean( (states_GT.q(:,iidxs)-hist.estim_q(:,iidxs)).^2 ));     
            max_pos=max(abs(states_GT.q(:,iidxs)));
            velocity_error = sqrt( mean( (states_GT.qp(:,iidxs)-hist.estim_qp(:,iidxs)).^2 ));   
            max_vel=max(abs(states_GT.qp(:,iidxs)));
            acceleration_error = sqrt( mean((states_GT.qpp(:,iidxs)-hist.estim_qpp(:,iidxs)).^2 ));
            max_acc=max(abs(states_GT.qpp(:,iidxs)));

            fprintf('RMSE: Pos=%f Vel=%f Acc=%f\n',[position_error' ,velocity_error', acceleration_error']');
            fprintf('RMSE: Pos=%.02f%% Vel=%.02f%% Acc=%.02f%% (%% wrt |max| GT value)\n',position_error/max_pos*100,velocity_error/max_vel*100,acceleration_error/max_acc*100);            
            
            % Save output stats:
            perfResults.pos_rme = position_error;
            perfResults.vel_rme = velocity_error;
            perfResults.acc_rme = acceleration_error;
            perfResults.max_maha_dist = max(maha_dist);
            perfResults.mean_maha_dist = mean(maha_dist);
            perfResults.max_maha_dist_pos = max(maha_dist_pos);
            perfResults.mean_maha_dist_pos = mean(maha_dist_pos);
            perfResults.max_maha_dist_vel = max(maha_dist_vel);
            perfResults.mean_maha_dist_vel = mean(maha_dist_vel);
            
        end % plot_estim_performance
        
        
        % To update stat files with results from run_offline()
        function [] = updateStatsResults(me,perfResults, td)
            % (From old code updateStatsResults.m)
            % INPUTS:
            label = td.label;
            sExp = class(me); %  sExp: observer under study,
            sErrorScale = sprintf('Errorscale%d',td.error_scale);
            sMultirate  = sprintf('Multirate%d',td.multirate_decimation);
%             elapsed_time = perfResults.elapsed_t;%  elapsed_time: (t/real-time)
%             maha_mean = perfResults.stats.mean_maha_dist;
%             maha_mean = 
            %  Estim data: errors.
           % warning('to do');
            %return;
            %sExp = params.observer;
%             if (isfield(params,'experiment'))
%                 sError = params.experiment;
%             else
%                 sError = 'ManuallyRun';
%             end

%             fprintf('Updating stats: sExp="%s" sError="%s"\n',sExp,sError);
            fprintf('Updating stats: "%s" \n',label);

%             % Average Mahalanobis distance
%             [ ~,~,MahaTotal ] = mahalanobis_mean_errors( GT, EstimData, hist_P );
%             maha_mean = mean(MahaTotal);

            % outputs:

            % One separate file for each parallel thread:
            %sStatFile = sprintf('results_stats_%d.mat',labindex());
            sStatFile = 'results_stats.mat';

            % Load:
            if (exist(sStatFile,'file')~=0)
                load(sStatFile);
            else
                results = struct();
            end

%             if (~isfield(results,sExp) || ~isfield(eval(sprintf('results.%s',sExp)),sError))
            if (~isfield(results,sExp) || ~isfield(eval(sprintf('results.%s',sExp)),sErrorScale) ||...
                    ~isfield(eval(sprintf('results.%s.%s',sExp,sErrorScale)),sMultirate))
                % Init all to empty vectors:
                eval(sprintf('results.%s.%s.%s.CPU=[];',sExp,sErrorScale,sMultirate));
                eval(sprintf('results.%s.%s.%s.RMSE_pos=[];',sExp,sErrorScale,sMultirate));
                eval(sprintf('results.%s.%s.%s.RMSE_vel=[];',sExp,sErrorScale,sMultirate));
                eval(sprintf('results.%s.%s.%s.RMSE_acc=[];',sExp,sErrorScale,sMultirate));
                eval(sprintf('results.%s.%s.%s.maha_mean=[];',sExp,sErrorScale,sMultirate));
            end

            % Append acumulative values:
            eval(sprintf('results.%s.%s.%s.CPU(end+1)=perfResults.elapsed_t;',sExp,sErrorScale,sMultirate));
            eval(sprintf('results.%s.%s.%s.RMSE_pos(end+1)=mean(perfResults.stats.pos_rme(:));',sExp,sErrorScale,sMultirate));
            eval(sprintf('results.%s.%s.%s.RMSE_vel(end+1)=mean(perfResults.stats.vel_rme(:));',sExp,sErrorScale,sMultirate));
            eval(sprintf('results.%s.%s.%s.RMSE_acc(end+1)=mean(perfResults.stats.acc_rme(:));',sExp,sErrorScale,sMultirate));
            eval(sprintf('results.%s.%s.%s.maha_mean(end+1)=perfResults.stats.mean_maha_dist;',sExp,sErrorScale,sMultirate));

            % Save back:
            save(sStatFile, 'results');
            
        end % end of updateStatsResults()
        
    end % end public methods
    
end

