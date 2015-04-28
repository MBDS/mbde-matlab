% Example of dynamic simulation

%%
addpath('toolbox_mbe');

%% Instance of multibody model:
mbs = mbeMechModelFourBars1();
%mbs = mbeMechModelFourBars2();
%mbs = mbeMechModelFiveBars1();
%mbs = mbeMechModelRodArray1dof();
%mbs = mbeMechModelPendulum1();

%% Instance of dynamical simulator: 
% sim = mbeDynFormulationMatrixR();
sim = mbeDynFormulationI3AL();
% sim = mbeDynFormulationPenalty();
% sim = mbeDynFormulationLagrangeBaumgarte();

%% Select integrator: 
% See options in mbeIntegratorTypes: intEuler, intTrapezoidal,...
% sim.integrator_type = mbeIntegratorTypes.intTrapezoidal;
sim.integrator_type = mbeIntegratorTypes.intI3AL;

% (Read the docs of mbeDynFormulationBase on this variable!)
%sim.integrate_indep_coords = 1;  % Default
%sim.integrate_indep_coords = 0; 

%% Run batch dyn simulation:
t_ini = 0.0;
t_end = 5;
dt    = 5e-3;
[states,perf_stats]=sim.batch_dyn_solve(mbs,t_ini,t_end,dt);

%% Plot results:
figure, 
ts=states.t(:); % time vector
plot(ts,states.q(:, mbs.indep_idxs(1) )),title('z(t) (1st DOF)');

% Show performance stats:
sim.plot_perf_stats(perf_stats);

