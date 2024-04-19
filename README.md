### MATLAB Toolbox: State Observers for Multibody Dynamic Systems ###

Toolbox features
------------------
  * Entirely based on MATLAB classes.
  * A generic framework for kinematics and dynamics simulation of rigid multibody models of mechanisms.
    * [`mbeKinematicsSolver`](toolbox_mbe/mbeKinematicsSolver.m): Solve kinematic problems (position, velocity, acceleration).
	* [`mbeDynFormulationBase`](toolbox_mbe/mbeDynFormulationBase.m): Abstract base of all dynamics formulations.
  * Mechanical systems can be easily defined by users by inheriting from the abstract base `mbeMechModelBase`. Predefined mechanisms: 
    * `mbeMechModelFourBars1`: 4 bars linkage
    * `mbeMechModelFiveBars1`: 5 bars linkage
  * Any number of sensors can be dynamically attached to any mechanism. See base abstract class `mbeSensorBase`
  * Several implemented state observers (base class `mbeEstimatorBase`).

Licensing
----------
MBDE-MATLAB is licensed by University of A Coruña and University of Almería under GNU GPL v3.

Citation
--------------
This toolbox was introduced in: 
  * E. Sanjurjo, J.L. Blanco, J.L. Torres, M.A. Naya, "Testing the efficiency and accuracy of multibody-based state observers", ECCOMAS Thematic Conference on Multibody Dynamics, 2015. [[PDF](http://ingmec.ual.es/~jlblanco/papers/sanjurjo2015eccomas_mbs_observers.pdf)]
  * E. Sanjurjo, M.Á. Naya, J.L. Blanco-Claraco, J.L. Torres-Moreno, A. Giménez-Fernández. Accuracy and efficiency comparison of various nonlinear Kalman filters applied to multibody models. Nonlinear Dyn (2017). doi:10.1007/s11071-017-3354-z [[LINK](http://rdcu.be/oY1f)]
  * E. Sanjurjo, D. Dopico, A. Luaces, M.Á. Naya. State and force observers based on multibody models and the indirect Kalman filter. Mechanical Systems and Signal Processing. doi:10.1016/j.ymssp.2017.12.041 [[LINK](https://www.sciencedirect.com/science/article/pii/S0888327017306787)]

Installation
--------------
Just add the directory `toolbox_mbe` to your MATLAB `PATH`, e.g. from the menu "Set path" in the "Home" ribbon.

Demos
----------
  * [`demo_dynamic_simulation.m`](demo_dynamic_simulation.m): Runs a dynamic simulation for a given mechanism.
  * [`demo_run_observer_offline.m`](demo_run_observer_offline.m): Runs an offline estimation of a given mechanism, with a given set of sensors and using the especific estimation method. 

Classes diagram
----------
Classes diagram in PNG format. Download the PDF version [`here`](https://github.com/MBDS/mbde-matlab/blob/e254f6c6eb4eb22bc8dea3a2d486e4c08c6dafa9/toolbox_mbe/classes_diagram.pdf) 
![classes_diagram](https://github.com/MBDS/mbde-matlab/assets/3897942/841e5709-f4e1-437b-a7f3-f57d34c9f9d0)
