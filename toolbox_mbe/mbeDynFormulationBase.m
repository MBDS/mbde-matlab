classdef (Abstract) mbeDynFormulationBase < handle
    % Virtual base for dynamic solvers. Offerts the service of dyn
    % simulations, displaying plots of their results, etc.
    % (portions based on old code 'MBS_solve')
    
    methods(Abstract,Access=public)
        % Abstract method for all dynamic solvers to calculate
        % instantaneous accelerations (both dep & indep accelerations may
        % be returned, depending on the method)
        %
        % [qpp,zpp] = solve_for_accelerations(me,model,q,qp,context)
        %
        % - context: Simulation context struct. Contents are: Current
        % time,dt,...
        % - Returned (qp, q) will normally be the same than passed as
        % input, except in "projection" methods.
        [qpp,zpp, q, qp] = solve_for_accelerations(me,model,q,qp,context);
    end
    
    properties(Access=public)
        % Whether to show messages, animations, etc.
        logging = mbeLoggingParams();
        
        % The integrator scheme for batch_dyn_solve()
        % From enumeration above
        integrator_type = mbeIntegratorTypes.intTrapezoidal;
        
        % Integration policy: 
        %  - "1" (default) -> numerical integration of **independent coords**, then solve
        %  pos, vel problems at each timestep (even if the formulation
        %  produces dependent accelerations, they are ignored).
        %  - "0" -> numerical integration of **dependent coords**, even if
        %  the formulation produces independent accelerations, the
        %  dependent accelerations are computed via a linear "acceleration
        %  problem".
        %
        % **Warning:** Disable indep coords integration only in dynamical methods that
        % have stabilization (e.g. penalizers, baumbarte)
        integrate_indep_coords = 1; 
        
        tol_trapezoidal = 1e-10; % dynamic tolerance (trapezoidal integrator)
        iter_max_trapezoidal = 100; % max iters for integration (trapezoidal integrator)
    end
    
    % --- public API ----
    methods(Access=public)
        function [states,perf_stats] = batch_dyn_solve(me,model,t_ini,t_end,dt)
            % Executes a batch simulation of a given multibody system with the
            % integrator defined in the property "integrator_type". 
            % Accelerations are obtained from the
            % specific method implemented in the derived class.
            % - model: The MB system parameters
            % - t_ini, t_end, dt: Initial, final and step time.
            % - states: A struct with Nx1 (N=nTimeSteps) vectors : q,qp,qpp, t
            % - perf_stats: A struct with energy & residues info
            %
            % See also: plot_perf_stats()
            
            % Initial conditions
            q_approx = model.q_init_aprox;
            assert(length(q_approx)==model.dep_coords_count);
            zp_init = model.zp_init; % Was: thetap = MBS.thetap; % Initial DOF velocity
            
            nTimeSteps = ceil((t_end-t_ini)/dt);
            
            % get indep/dep coord indices:
            iidxs = model.get_indep_indxs();
            
            % Integration context (aux vars,...)
            context = struct();
            context.model = model;
            context.iidxs = iidxs;
            context.dt = dt;
            context.current_run_time = t_ini;            

            % initial position, velocity and acceleration problems
            q = mbeKinematicsSolver.pos_problem(model,q_approx);
            qp = mbeKinematicsSolver.vel_problem(model,q, zp_init);
            [qpp,~,q,qp] = me.solve_for_accelerations(model,q,qp, context);
            

            % TODO: Detect the optional residual output vars, don't waste
            % time computing them if not present.
            do_calc_perf_indicators = 1;

            % Output matrix:
            states = struct();
            states.t   = zeros(nTimeSteps,1);
            states.q   = zeros(nTimeSteps,model.dep_coords_count);
            states.qp  = zeros(nTimeSteps,model.dep_coords_count);
            states.qpp = zeros(nTimeSteps,model.dep_coords_count);
            
            % Verification variables
            perf_stats = struct();
            if (do_calc_perf_indicators)
                perf_stats.T = zeros (1, nTimeSteps);
                perf_stats.Vg = zeros (1, nTimeSteps);
                perf_stats.E_mec = zeros (1, nTimeSteps);
                % Constraints fulfilment
                perf_stats.resid_pos = zeros(1, nTimeSteps);
                perf_stats.resid_vel = zeros(1, nTimeSteps);
                perf_stats.resid_acc = zeros(1, nTimeSteps);
            end
            
            % For logging progress:
            logstep_last = 0;
            logstep_nSteps = 10;
            logstep_incr = nTimeSteps/logstep_nSteps;
            
            % For debugging animations:
            liveanim_last = 0;
            if (me.logging.show_mech_anim_fps>0)
                liveanim_fig = figure();
                liveanim_incr = (1.0/me.logging.show_mech_anim_fps)/dt;
            else
                liveanim_incr = Inf; % Disabled
            end
            
            % MBS loop
            for i = 1:nTimeSteps,
                % Set iteration:
                context.current_run_time = t_ini + i*dt;
                if (i>=logstep_last+logstep_incr || i==nTimeSteps)
                    pcdone = 100.0 * i / nTimeSteps;
                    me.logging.log_level(1, sprintf('batch_dyn_solve: %.02f%% done. iter %u/%u, simul_time=%.03fs/%.03fs',...
                        pcdone, ...
                        i,nTimeSteps, ...
                        context.current_run_time,t_end ));
                    logstep_last = i;
                end

                % Perform numerical integration:
                % ----------------------------------
                switch me.integrator_type
                    case mbeIntegratorTypes.intEuler
                        [q,qp,qpp] = me.ni_step_euler(q,qp,qpp, context);
                    case mbeIntegratorTypes.intTrapezoidal
                        [q,qp,qpp] = me.ni_step_trapezoidal(q,qp,qpp, context);
                    case mbeIntegratorTypes.intI3AL
                        [q,qp,qpp] = me.ni_step_I3AL(q,qp,qpp, context);
                    otherwise
                        error('Unknown integrator scheme!');
                end

                % Store dynamical results:
                states.t(i) = context.current_run_time;
                states.q(i,:)   = q';
                states.qp(i,:)  = qp';
                states.qpp(i,:) = qpp';

                % Store performance data:
                if (do_calc_perf_indicators)
                    perf_stats.T(i) = 0.5*qp'*model.M*qp;
                    perf_stats.Vg(i) = model.eval_potential_energy(q);
                    perf_stats.E_mec(i) = perf_stats.T(i)+perf_stats.Vg(i);
                    phi_q = model.jacob_phi_q(q);
                    perf_stats.resid_pos(i) = norm(model.phi(q));
                    perf_stats.resid_vel(i) = norm(phi_q*qp);
                    perf_stats.resid_acc(i) = norm(model.jacob_phiqp_times_qp(q,qp)+phi_q*qpp);
                end
                
                % Debug animations: GT, BadModel, Estim
                if (i>=liveanim_last+liveanim_incr)
                    clf; hold on;
                    % Plot mechanism state:
                    model.plot_model_skeleton(q, 'b', 1); % the last 1=fit plot axes
                    % Title:
                    title('Mechanism skeleton visualization');
                    drawnow;
                    liveanim_last=i;
                end
                
                
            end % end for each timestep
        end % batch_dyn_solve
        
        function [q_next,qp_next,qpp_next] = integrate_one_timestep(me,model,current_time, dt, q,qp,qpp )
            % perform one single integration timestep.
            
            % Integration context (aux vars,...)
            context = struct();
            context.model = model;
            context.iidxs = model.get_indep_indxs();
            context.dt = dt;
            context.current_run_time = current_time;            
            
            % Perform numerical integration:
            % ----------------------------------
            switch me.integrator_type
                case mbeIntegratorTypes.intEuler
                    [q_next,qp_next,qpp_next] = me.ni_step_euler(q,qp,qpp, context);
                case mbeIntegratorTypes.intTrapezoidal
                    [q_next,qp_next,qpp_next] = me.ni_step_trapezoidal(q,qp,qpp, context);
                case mbeIntegratorTypes.intI3AL
                    [q_next,qp_next,qpp_next] = me.ni_step_I3AL(q,qp,qpp, context);
                otherwise
                    error('Unknown integrator scheme!');
            end            
        end % integrate_one_timestep()
        
        
    end % end of public methods
        
    methods(Static,Access=public)
        function [] = plot_perf_stats(s)
            % Shows the results from batch_dyn_solve: Conservation of
            % energy and constraints fulfilment
            %  [] = plot_perf_stats(s)
            %  - s: perf_stats from batch_dyn_solve()
            % (Old code: MBS_plots)
            long = length(s.E_mec);
            figure
            clf
            hold on
                grid on
                plot (1:long,s.E_mec, 'r')
                plot (1:long,s.T, 'g')
                plot (1:long,s.Vg, 'b')
                EnCon = abs(100 - abs(max(s.E_mec) - min(s.E_mec))*100/s.E_mec(1));
                str = sprintf('Conservation of energy: %.3f %%\n',EnCon);
                Box = uicontrol('style','text');
                set(Box,'String',str)
                set(Box,'Position',[80,-20,200,35])
                legend ('Mechanical', 'Kinetic ', 'Potential')
                title('Conservation of energy')
            hold off

            figure
            clf
            hold on
                grid on
                plot (1:long, s.resid_pos, 'r')
                plot (1:long, s.resid_vel, 'b')
                plot (1:long, s.resid_acc, 'g')
                legend ('Position', 'Velocity', 'Acceleration')
                title('Constraints fulfilment')
            hold off            
        end % plot_perf_stats        
        
    end % end of public (static) methods    
    
    
    methods(Access=protected)
        function [q_next,qp_next,qpp_next] = ni_step_I3AL(me, q,qp,qpp, context)
            % formulation-specific formula:
            [qpp_next,zpp_next,q_next, qp_next] = solve_for_accelerations(me, context.model,q,qp,  context);
        end
        function [q_next,qp_next,qpp_next] = ni_step_trapezoidal(me, q,qp,qpp, context)
        % NUM. INTEGRATION IMPLEMENTATION: Trapezoidal rule
            dt = context.dt;
            if (me.integrate_indep_coords)
                iidxs=context.iidxs;
                z = q(iidxs); 
                zp = qp(iidxs);
                zpp = qpp(iidxs);
                zpp_next=zpp; 
            end
            
            err = 1;  ite = 0;
            
            % Initial values (will be refined in the first loop below):
            q_next = q; qpp_next = qpp; 

            while err>me.tol_trapezoidal && ite < me.iter_max_trapezoidal
                q_old = q_next;

                if (me.integrate_indep_coords)
                    z_next = z + zp*dt + dt^2/4*(zpp+zpp_next);
                    zp_next = zp + dt/2*(zpp+zpp_next);
                    q_next(iidxs) = z_next;
                    q_next = mbeKinematicsSolver.pos_problem( context.model,q_next); %ini_pos(q_next,PARAMS);
                    qp_next = mbeKinematicsSolver.vel_problem( context.model,q_next,zp_next); % Solve for dep. vels
                else
                    q_next = q + qp*dt + dt^2/4*(qpp+qpp_next);
                    qp_next = qp + dt/2*(qpp+qpp_next);
                end
                % New acceleration vector: invoke the
                % formulation-specific formula:
                [qpp_next,zpp_next,q_next, qp_next] = solve_for_accelerations(me, context.model,q_next,qp_next,  context);

                % next iter:
                ite=ite+1;
                err = norm (q_old-q_next);
            end
        end % ni_step_trapezoidal
        
        function [q_next,qp_next,qpp_next] = ni_step_euler(me, q,qp,qpp, context)
        % NUM. INTEGRATION IMPLEMENTATION: Explicit Euler rule
            dt = context.dt;
            
            if (me.integrate_indep_coords)
                iidxs=context.iidxs;
                z = q(iidxs); 
                zp = qp(iidxs);
                zpp = qpp(iidxs);

                z_next  = z  + zp  * dt;
                zp_next = zp + zpp * dt;
                
                q_next = q; % Init value
                q_next(iidxs) = z_next;
                q_next = mbeKinematicsSolver.pos_problem( context.model,q_next); %ini_pos(q_next,PARAMS);
                qp_next = mbeKinematicsSolver.vel_problem( context.model,q_next,zp_next); % Solve for dep. vels
            else
                q_next  = q  + qp  * dt;
                qp_next = qp + qpp * dt;
            end
            
            % New acceleration vector: invoke the formulation-specific formula:
            [qpp_next,~,q_next, qp_next] = solve_for_accelerations(me, context.model,q_next,qp_next,context);
            
        end % ni_step_euler
        
    end
    
end

