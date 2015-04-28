### MATLAB Toolbox: State Observers for Multibody Dynamic Systems ###

Toolbox features
------------------
  * Entirely based on MATLAB classes.
  * A generic framework for kinematics and dynamics simulation of rigid multibody models of mechanisms.
    * `mbeKinematicsSolver`: Solve kinematic problems (position, velocity, acceleration).
	* `mbeDynFormulationBase`: Abstract base of all dynamics formulations.
  * Mechanical systems can be easily defined by users by inheriting from the abstract base `mbeMechModelBase`. Predefined mechanisms: 
    * `mbeMechModelFourBars1`: 4 bars linkage
    * `mbeMechModelFiveBars1`: 4 bars linkage
  * Any number of sensors can be dynamically attached to any mechanism. See base abstract class `mbeSensorBase`
  * Several implemented state observers (base class `mbeEstimatorBase`).

Licensing
----------
MBDE-MATLAB is licensed by University of Coruña and University of Almería under GNU GPL v3.

Citation
--------------
This toolbox was introduced in: 
  * J.L. Torres, J.L. Blanco, E. Sanjurjo, A. Gimenez-Fernandez, M.A. Naya, "A testbed for benchmarking state observers in Multibody Dynamics", ECCOMAS Thematic Conference on Multibody Dynamics, 2015.

Installation
--------------
Just add the directory `toolbox_mbe` to your MATLAB `PATH`, e.g. from the menu "Set path" in the "Home" ribbon.

Demos
----------
  * `demo_dynamic_simulation.m`: Runs a dynamic simulation for a given mechanism.
  * `demo_run_observer_offline.m`: Run an offline estimation of a given mechanism, with a given set of sensors and using the especific estimation method. 

This file is part of MBDE-MATLAB.

    MBDE-MATLAB is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MBDE-MATLAB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MBDE-MATLAB.  If not, see <http://www.gnu.org/licenses/>.
