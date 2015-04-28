classdef mbeEstimatorBatchMRF < mbeEstimatorBase
    % Batch estimator based on a Markov Random Field (MRF) over the whole history.
    % Trapezoidal integrator.
    %
    % TODO: Allow sensors to be available only at some timesteps!
    
    % NOTE ON THE VARIABLES ORDER IN THE STATE VECTOR (affecting
    % information matrix, Jacobians, etc.): 
    %
    %  x = [ z^1_1...z^g_1 | dz^1_1...dz^g_1 | ddz^1_1...ddz^g_1 | ... | z^1_N...z^g_N | dz^1_N...dz^g_N | ddz^1_N...ddz^g_N ];
    % with: g=number of degrees of freedom; N=number of timesteps
    %

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

    
    properties(Access=public)
        SIGMA_init_pos = deg2rad(20);   % Transition noise in "z" (assumed by the estimator)
        
        SIGMAp = deg2rad(0.1);  % Transition noise in "z" (assumed by the estimator)
        SIGMAr = deg2rad(0.5);  % Transition noise in "zp" (assumed by the estimator)
        SIGMAa = deg2rad(200);  % Transition noise in "zpp" (assumed by the estimator)
        
        
        dynamics_formulation = mbeDynFormulationLagrangeBaumgarte(); %mbeDynFormulationMatrixR();
        
        % Show plots with state vector during optimization (!=0 means true)
        debug_show_live_x_evolution = 1;

        % Minimum |Ax|/len(x) to continue iterating
        epsilon_x = 1e-5;  

        % Maximum optimization iters:
        max_iters = 50;        

        % Set to <1 to enable adaptive re-weighting of constraints
        ADAPTIVE_WEIGHT_INIT_MINIMUM = 1; %0.01;
    end
    
    % Public read-only:
    properties(GetAccess=public,SetAccess=private)
        % Mechanism models  filters may need to know the current "time",
        % in seconds since the beginning of the experiment:
        current_run_time = 0.0;
        % Like current_run_time but counts discrete timesteps.
        current_timestep = 0;
 
    end
    
    % Private vars
    properties(Access=private)
        lenQ;  % Number of coords: |q|
        lenZ;  % Number of indep coords: |z|
        iidxs; % Indices of independent coords in q()
        didxs; % Indices of dependent coords in q()

        % State vector (the reconstructed trajectory): indep coords ONLY
        x_z;    % nTimeSteps x lenZ
        x_zp;   % nTimeSteps x lenZ
        x_zpp;  % nTimeSteps x lenZ
        
        J; % The sparse Jacobian
        
        debug_show_live_x_fig_handle;
        numeric_jacob_params = struct();
    end

    % --- Public API ---
    methods
        function [perfStats] = run_offline(me)
            % Execute one complete offline experiment
            
            % Init vars:
            me.lenQ = me.mech_phys_model.dep_coords_count();
            me.lenZ = length( me.mech_phys_model.indep_idxs );
            me.iidxs = me.mech_phys_model.get_indep_indxs();
            me.didxs = me.mech_phys_model.get_dep_indxs();
            
            me.init_bad_model(); % Make sure "bad model" is up-to-date
            
            me.mechanism_type.resetState(); % Prepare for simulating from t=0
            
            % Load data: 
            [states_GT, states_BadModel] = me.mechanism_type.generate_ground_truth(me);

            % Experiment duration:
            nTimeSteps = ceil(me.end_time/me.dt);
            assert(size(states_GT.q,1)==nTimeSteps,'Mismatch in GT length!')
            assert(size(states_BadModel.q,1)==nTimeSteps,'Mismatch in BadModel length!')
            
            % get all the sensor observations at once:
            % -------------------------------------------
            % nTimeSteps x nObs matrix
            % obs2timestep: observation index -> timestep (because not all
            % timesteps may have an observation)
            [gt_observations, obs2timestep, timestep2obsidx] = me.batch_get_all_real_observations(nTimeSteps, states_GT);
            
            % Temporary buffer with the complete dynamic reconstruction: 
            % dep and indep coords:
            dyn_states = struct();
            dyn_states.q   = zeros(nTimeSteps, me.lenQ );
            dyn_states.qp  = zeros(nTimeSteps, me.lenQ );
            dyn_states.qpp = zeros(nTimeSteps, me.lenQ );
            
            % Initial state vector mean (from the "bad model") 
            % -------------------------------------------------
            % Make sure it is self-consistent for...
            % ... pos:
            dyn_states.q(1,:) = mbeKinematicsSolver.pos_problem(...
                me.bad_mech_phys_model, ...
                me.bad_mech_phys_model.q_init_aprox);
            
            % ... vel:
            dyn_states.qp(1,:) = mbeKinematicsSolver.vel_problem(...
                me.bad_mech_phys_model,...
                dyn_states.q(1,:)',...
                me.bad_mech_phys_model.zp_init);

            % ... accel to zero by default:

            % Init state vector:
            % [ z1,zp1,zp1,  ... , zL,zpL,zpL]
            %       t=1              t=L
%             me.x_z   = zeros( nTimeSteps, me.lenZ); % Use these if we don't have a GT...
%             me.x_zp  = zeros( nTimeSteps,me.lenZ);
%             me.x_zpp = zeros( nTimeSteps,me.lenZ);

            % Initial value for the optimization: 
            %  populate with the "bad model" dynamic simulation
            me.x_z   = states_BadModel.q  (:,me.iidxs);
            me.x_zp  = states_BadModel.qp (:,me.iidxs);
            me.x_zpp = states_BadModel.qpp(:,me.iidxs);

            % Reserve stats for benchmark:
            hist=struct();
            if (me.do_benchmark)
                hist.t  = ((0:(nTimeSteps-1)) * me.dt)';  % time of each timestep
                hist.estim_q   = zeros(nTimeSteps,me.lenQ);
                hist.estim_qp  = zeros(nTimeSteps,me.lenQ);
                hist.estim_qpp = zeros(nTimeSteps,me.lenQ);
                hist.P         = cell(nTimeSteps,1); % Cov history
                hist.Pidxs_z = 0;
                hist.Pidxs_zp= 0;
                hist.Pidxs_zpp = 0;
                hist.sensor_data = cell(nTimeSteps,1); % cell: Variable number of readings in each timestep!
            end

            % Information ("weight") matrix:
            % ------------------------------------
            sensors_stds = me.sensors_std_magnification4filter * me.bad_mech_phys_model.sensors_std_noise();

            % Non-adaptive weights: 
            if (me.ADAPTIVE_WEIGHT_INIT_MINIMUM==1)
                weights= ones(nTimeSteps,1);
                VLim = 1;
                L = me.buildL(nTimeSteps,obs2timestep,sensors_stds, weights );
            end
            
            % Start of iterative estimation
            %----------------------------------------
            iter=0;
            norm_Ax_avr=1;
            lambda = 1e4;  % mix 2nd and 1st order optimization
            
            %tic            
            
            % Generate  first residuals:
            dyn_states= me.batchUpdateDependentPosVelAcc(nTimeSteps,dyn_states);
            r= me.simulAndGenerateResiduals(obs2timestep,timestep2obsidx,dyn_states,gt_observations);
            r_norm = norm(r);
            me.log_level(1, sprintf('Iter: %i, |r|:%.5f',iter,r_norm));

            % ========= GAUSS-NEWTON =============
            while (iter<me.max_iters && norm_Ax_avr>me.epsilon_x)
                
                % Adaptative Information matrix:
                if (me.ADAPTIVE_WEIGHT_INIT_MINIMUM<1)
                    % Constant to make the end-time constraint to be VLim
                    VLim = 1 - (1-me.ADAPTIVE_WEIGHT_INIT_MINIMUM)*exp( -iter/15 ); %1e-2;
                    kweights = -log(VLim)/nTimeSteps;
                    weights = exp(-kweights*(1:nTimeSteps)');
                    L = me.buildL(nTimeSteps,obs2timestep,sensors_stds, weights );
                end

                
                % Debug: live plots
                % ------------------
                if (me.debug_show_live_x_evolution)
                    if (isempty(me.debug_show_live_x_fig_handle))
                        me.debug_show_live_x_fig_handle = figure();
                    else
                        figure(me.debug_show_live_x_fig_handle);
                    end
                    
                    % Pos: 
                    subplot(4,1,1); cla; hold on;
                    plot(dyn_states.q(:,me.iidxs),'r');
                    plot(states_GT.q(:,me.iidxs),'k');
                    title('z'); legend('Estim','GT');
                    
                    % Vel: 
                    subplot(4,1,2); cla; hold on;
                    plot(dyn_states.qp(:,me.iidxs),'r');
                    plot(states_GT.qp(:,me.iidxs),'k');
                    title('zp'); legend('Estim','GT');

                    % Acc: 
                    subplot(4,1,3); cla; hold on;
                    plot(dyn_states.qpp(:,me.iidxs),'r');
                    plot(states_GT.qpp(:,me.iidxs),'k');
                    title('zpp'); legend('Estim','GT');
                    
                    
                    % weights: 
                    subplot(4,1,4); cla; hold on;
                    plot(weights);axis([1 nTimeSteps 0 1]);
                    title('Relative weights');
                    
                    drawnow;                   
                    %pause;
                end % end debug plots
                

                % Solve for increments Ax:
                % Gauss-Newton:
                me.updateJacobian(obs2timestep,timestep2obsidx,dyn_states);
                %me.updateJacobianNumeric(obs2timestep,timestep2obsidx,dyn_states);
                
                
                H = me.J' * L * me.J;
                JtL = me.J'*L;

                g = JtL * r;
                Ax =  - (H+speye(size(H,1))*lambda ) \ g;
                lambda = lambda*0.3; % LevMarq like

                norm_Ax_avr=norm(Ax)/(3*length(me.x_z));

                %x_new = x + Ax;
                % -------------------
                % Reshape as: Ax_mat = [ Az | Azp | Azpp ]
                Ax_mat = reshape(Ax, 3*me.lenZ,nTimeSteps)';  %Ax_mat = reshape(Ax_mat(:, (1:me.lenZ)), nTimeSteps,me.lenZ );
                
                me.x_z   = me.x_z   + Ax_mat(:,(1:me.lenZ));
                me.x_zp  = me.x_zp  + Ax_mat(:,(1:me.lenZ)+me.lenZ);
                me.x_zpp = me.x_zpp + Ax_mat(:,(1:me.lenZ)+2*me.lenZ);
                
                % Re-evaluate errors:
                % Note: batchUpdateDependentPosVelAcc() takes the values from x_z, x_zp, x_zpp
                dyn_states= me.batchUpdateDependentPosVelAcc(nTimeSteps,dyn_states); 
                r_new= me.simulAndGenerateResiduals(obs2timestep,timestep2obsidx,dyn_states,gt_observations);
                r_norm_new = norm(r_new);
                
                % LM: check new error: 
                r_norm_prev = r_norm; 
                r=r_new;
                r_norm=r_norm_new;

                err_improv = r_norm_new/r_norm_prev;

                iter=iter+1;

                me.log_level(1, sprintf('Iter: %i, |Ax|avr=%.03e |r|:%.5f->%5f (improv=%.03f%%); MinWght=%.04f',iter,norm_Ax_avr,r_norm_prev,r_norm_new,100*(1-err_improv),VLim));

            end % End of iterative optimizer

            perfStats = struct();
            % Plots, stats...
            if (me.do_benchmark)
                hist.estim_q   = dyn_states.q;
                hist.estim_qp  = dyn_states.qp;
                hist.estim_qpp = dyn_states.qpp;
                
                % This is *really* costly...
                me.log_level(1, 'Retrieving covariance...');
                t_P= tic;
                global_P = inv(H);
                me.log_level(1, sprintf('Covariance done: %.03f secs',toc(t_P)) );
                
                for i=1:nTimeSteps,
                    idxs = (3*(i-1)*me.lenZ) + (1:(3*me.lenZ));
                    hist.P{i} = full(global_P(idxs,idxs));
                end
                
                % Indices of each coord in each P{i}
                hist.Pidxs_z   = 1:me.lenZ;
                hist.Pidxs_zp  = (1:me.lenZ) + me.lenZ;
                hist.Pidxs_zpp = (1:me.lenZ) + 2*me.lenZ;

                % Save output performance report:
                perfStats.hist = hist;                
                perfStats.states_GT = states_GT;
                perfStats.states_BadModel = states_BadModel;
                
                % Eval and/or plot stats:
                perfStats.stats = me.plot_estim_performance(states_GT, states_BadModel,hist,me.show_benchmark_graphs);
            end

        end % end run_offline() 
       
    end % end public methods
        
    methods(Access=private)
        
        function [r]= simulAndGenerateResiduals(me,obs2timestep,timestep2obsidx,dyn_states,obs)
            % Compute all residuals:
            % obs: t x nObs matrix. Each row=sensors for one timestep
            %
            nTimeSteps = size(dyn_states.q,1);
            nObs = length(obs2timestep);
            lenObs = size(obs,2);
            assert(size(obs,1)==nObs);
            
            % Generate sensor prediction from x,  h(x)
            r_s_mat = zeros(nObs,lenObs);  % Was: %r_s = zeros(nTimeSteps*nObs,1);
            
            % Generate sensor prediction from eom
            r_m_mat = zeros(2*(nTimeSteps-1),me.lenZ); % Was: r_m = zeros(2*(nTimeSteps-1),1);
            r_eom_mat = zeros(nTimeSteps,me.lenZ); % Was: r_eom = zeros(nTimeSteps,1);

            % Update residuals:
            for i=1:nTimeSteps,
                obs_idx = timestep2obsidx(i);
                if (obs_idx~=0)
                    % Equations of errors in sensors:
                    %r_s(i) = virtual_sensor_eval (EstimData(i,1:5)',EstimData(i,6:10)', MBK ) - observations(i);   
                    obs_predict = me.bad_mech_phys_model.sensors_simulate(...
                        dyn_states.q(i,:)', dyn_states.qp(i,:)', dyn_states.qpp(i,:)' );


                    r_s_mat(obs_idx,:) = obs_predict' - obs(obs_idx,:);
                end
                
                % eom:
                [~, zpp_t_pred]= me.dynamics_formulation.solve_for_accelerations(...
                    me.bad_mech_phys_model, ...
                    dyn_states.q(i,:)', dyn_states.qp(i,:)', ...
                    struct() ... % context
                    );
                
                % Errors in predicted accelerations:
                %r_eom(i)= x((i-1)*3+3)-zpp_t_pred;
                r_eom_mat(i,:) = me.x_zpp(i,:) - zpp_t_pred';

                % Integrator:
                if(i<nTimeSteps)
                    z_t   = me.x_z(i,:);
                    zp_t  = me.x_zp(i,:);
                    zpp_t = me.x_zpp(i,:);

                    z_tp1   = me.x_z(i+1,:);
                    zp_tp1  = me.x_zp(i+1,:);
                    zpp_tp1 = me.x_zpp(i+1,:);

                    r_m_mat((i-1)*2+1,:) =  z_t +   me.dt*0.5*(zp_t+zp_tp1)       - z_tp1;
                    r_m_mat((i-1)*2+2,:) = zp_t +   me.dt*0.5*(zpp_t+zpp_tp1)     - zp_tp1;
                end
            end
            
            % r_inipos
            initguess_z0 = me.bad_mech_phys_model.q_init_aprox(me.iidxs);
            r_inipos = me.x_z(1,:)' - initguess_z0;

            % Residuals: r
            r_s = reshape(r_s_mat', nObs*lenObs,1);
            r_m = reshape(r_m_mat', 2*(nTimeSteps-1)*me.lenZ,1);
            r_eom = reshape(r_eom_mat', nTimeSteps*me.lenZ,1);
            r = [ r_s; r_eom; r_m; r_inipos ];
        end % bui

        % Only for debugging the accuracy of the closed-form Jacobian
        function [] = updateJacobianNumeric(me,obs2timestep,timestep2obsidx,dyn_states)
            updateJacobian(me,obs2timestep,timestep2obsidx,dyn_states);
            
            me.numeric_jacob_params.dyn_states = dyn_states;
            %me.numeric_jacob_params.
            Jnum = jacobianest(numeric_error_func,x0,me);
        end
        % Only used in updateJacobianNumeric()
        function [residuals] = numeric_error_func(me,x0)
            me.numeric_jacob_params.dyn_states
            
        end
        
        
        function [] = updateJacobian(me,obs2timestep,timestep2obsidx,dyn_states)
            % Compute J, the complete sparse Jacobian
            
            lZ=me.lenZ;
            dt=me.dt;
            nTimeSteps = size(dyn_states.q,1);
            nObs = length(obs2timestep);
            lenObs = length(me.bad_mech_phys_model.installed_sensors);
            
            % Init empty Jacobian:
            nRows= nObs*lenObs + ...    % Sensors
                   nTimeSteps*lZ + ...   % eom
                   2*(nTimeSteps-1)*lZ + ...  % integrator
                   1*lZ;      % Initial pos guess
               
            nCol= 3*nTimeSteps*lZ;  % (z,zp,zpp) for each timestep
            
            if (isempty(me.J))
                estimNZs = nTimeSteps*lZ^2 + 3*nTimeSteps + 2*lZ*nTimeSteps + 1*nTimeSteps;
                me.J= spalloc(nRows,nCol, estimNZs); %sparse(nRows,nCol);
            end
            
            % Sensors:
            % --------------------
            for state_idx=1:nTimeSteps,
                obs_idx = timestep2obsidx(state_idx);
                if (obs_idx~=0)
                    col=(state_idx-1)*3*lZ;
                    row=(obs_idx-1)*lenObs;
                    % Eval observation Jacobian wrt the indep. coordinates:
                    [dh_dz , dh_dzp, dh_dzpp] = me.bad_mech_phys_model.sensors_jacob_indep(...
                        dyn_states.q(state_idx,:)', dyn_states.qp(state_idx,:)', dyn_states.qpp(state_idx,:)' );

                    me.J( row+(1:lenObs), col+(1:(3*lZ)) ) = [dh_dz , dh_dzp, dh_dzpp];
                end
            end
            last_block_rows = nObs*lenObs;
            
            % equations of motion:
            % --------------------
            %J_eom=[0 0 1];
            for i=(1:nTimeSteps),
                col=(i-1)*3*lZ;
                row= last_block_rows + (i-1)*lZ;

                % Numeric Jacobian estimations:
                % dzpp/dz
            %     Delta_z = 1e-3;
            %     val1p = solve_eom(i,x,EstimData,MBK, +Delta_z, 0);
            %     val1m = solve_eom(i,x,EstimData,MBK, -Delta_z, 0);
            %     dzpp_dz = (val1p-val1m)/2*Delta_z;

                % dzpp/dzp
            %     Delta_zp = 1e-3;
            %     val1p = solve_eom(i,x,EstimData,MBK, 0, +Delta_zp);
            %     val1m = solve_eom(i,x,EstimData,MBK, 0, -Delta_zp);
            %     dzpp_dzp = (val1p-val1m)/2*Delta_zp;
            %     
            % ==> They are almost exactly zero!!
%                 dzpp_dz = 0;
%                 dzpp_dzp = 0;

                % Fill in:
%                 me.J(i+nTimeSteps,col  )= -dzpp_dz;     % Don't fill-in zeros!!
%                 me.J(i+nTimeSteps,col+1)= -dzpp_dzp; 

                me.J(row+ (1:lZ),col+2*lZ+ (1:lZ) )= eye(lZ);  % dzpp/dzpp
            end
            last_block_rows=last_block_rows+  nTimeSteps*lZ;

            % Integrator:
            % --------------------
            dt05 = 0.5*dt;
            for i=(1:nTimeSteps-1),
                row=last_block_rows + 2*lZ*(i-1);
                col=(i-1)*lZ*3;
                
                me.J(row+(1:(2*lZ)),col+(1:(2*3*lZ))) = ...
                    [   eye(lZ), dt05*eye(lZ),    zeros(lZ),  -eye(lZ), dt05*eye(lZ),    zeros(lZ) ;...
                       zeros(lZ),     eye(lZ), dt05*eye(lZ), zeros(lZ),     -eye(lZ), dt05*eye(lZ) ];
%                    [ 1 dt05  0   -1  dt05  0 ;...
%                      0  1    dt05   0   -1    dt05 ];
            end
            last_block_rows=last_block_rows+  2*(nTimeSteps-1)*lZ;
            
            % Initial pos guess
            row=last_block_rows ;
            col=0;
            me.J(row+(1:lZ),col+(1:lZ)) = eye(lZ);
                            
            assert(size(me.J,1)==nRows && size(me.J,2)==nCol); % Detect possible out-of-bounds indices!
            
        end % updateJacobian
        
 
        
        function L = buildL(me,nTimeSteps,obs2timestep,SIGMAsensors, time_reweight )
            % Build the sparse diagonal infomation matrix:
            % - SIGMAsensors: a vector with the std dev of sensor obs.
            % - SIGMAp: sigma of integrator (transition in pos)
            % - SIGMAr: sigma of integrator (transition in velocity)
            % - SIGMAa: sigma of acceleration (eom model vs state)
            % - time_reweight: empty=all timesteps have equal weight.
            % Optionally, a vector of length nTimeSteps with relative
            % weights for the costs of each timestep constraints. Larger
            % values means that the constraints must be fulfilled more
            % closely, values close to zero imply making the constraint more loosy.
            % 
            
            nObservations = length(obs2timestep);
            
            if (isempty(time_reweight))
                time_reweight=ones(nTimeSteps,1);
            end            

            lZ=me.lenZ;
            lenObs = length(me.bad_mech_phys_model.installed_sensors);
            assert(lenObs==length(SIGMAsensors));

            % Convert to "information":
            inf_s = (1./SIGMAsensors).^2;
            inf_p = (1/me.SIGMAp)^2;
            inf_r = (1/me.SIGMAr)^2;
            inf_a = (1/me.SIGMAa)^2;
            inf_init_pos = (1/me.SIGMA_init_pos)^2;

            % Build the diagonal:
            D=[ kron(ones(nObservations,1),inf_s);...   %ones(t,1)*inf_s; ...
                kron(time_reweight, ones(lZ,1))*inf_a; ...
                kron(time_reweight(1:(end-1)),[ones(lZ,1)*inf_p;ones(lZ,1)*inf_r]);...
                ones(1*lZ,1)*inf_init_pos
                ];

            % Convert to sparse:
            %N=nTimeSteps+nTimeSteps+2*(nTimeSteps-1)+1;
            
            % |L| = rows in Jacobian:
            nRows= nObservations*lenObs + ...    % Sensors
                   nTimeSteps*lZ + ...   % eom
                   2*(nTimeSteps-1)*lZ + ...  % integrator
                   1*lZ;      % Initial pos guess            
            
            L = spdiags(D,0,nRows,nRows);
        end % buildL()
        
        function [dynstates_out]= batchUpdateDependentPosVelAcc(me, nTimeSteps, dynstates)
            % Update the list of all q,qp in EstimData from the state vector x
            dynstates_out=dynstates;
            
            % Update q, qp,qpp
            for i=1:nTimeSteps
                % Pos:
                dynstates_out.q(i,me.iidxs) = me.x_z(i,:)';
                dynstates_out.q(i,:)  = mbeKinematicsSolver.pos_problem( ...
                    me.bad_mech_phys_model,dynstates_out.q(i,:)' )';
                % Vel:
                dynstates_out.qp(i,:) = mbeKinematicsSolver.vel_problem( ...
                    me.bad_mech_phys_model, dynstates_out.q(i,:)', me.x_zp(i,:)')';
                % Accel:
                dynstates_out.qpp(i,:) = mbeKinematicsSolver.accel_problem(...
                    me.bad_mech_phys_model, ...
                    dynstates_out.q(i,:)', dynstates_out.qp(i,:)', me.x_zpp(i,:)')';

                % Help the initial position problem in (i+1), which for default 
                % has a q=[0...0]:
                if (i<nTimeSteps),
                    dynstates_out.q(i+1,:) = dynstates_out.q(i,:) + me.dt * dynstates_out.qp(i,:);
                end
            end
        end
        
        function [obs_out,obs2timestep,timestep2obsidx] = batch_get_all_real_observations(me, nTimeSteps, states)
            obs_out=[];
            nActualObservations = 0;
            timestep2obsidx=zeros(nTimeSteps,1);
            for i = 1:nTimeSteps
                % Set iteration:
                me.current_timestep = i;
                me.current_run_time = i*me.dt;

                obs = me.mechanism_type.get_current_observations( ...
                    me, states.q(i,:), states.qp(i,:), states.qpp(i,:) );

                if (~isempty(obs))
                    nSensors = length(obs);
                    if (isempty(obs_out))
                        obs_out = zeros(nTimeSteps, nSensors);
                        obs2timestep=zeros(nTimeSteps,1);
                    end

                    nActualObservations=nActualObservations+1;
                    obs_out(nActualObservations,:) = obs';                 %#ok<AGROW>  % Suppress MATLAB warning here!
                    obs2timestep(nActualObservations)=i;

                    timestep2obsidx(i)=nActualObservations;
                end
            end
            
            assert(nActualObservations>0);
            obs_out = obs_out(1:nActualObservations,:);
            obs2timestep = obs2timestep(1:nActualObservations);
            
        end % batch_get_all_real_observations()
        

    end % end private methods
   
        
end

