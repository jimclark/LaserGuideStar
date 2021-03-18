# LaserGuideStar
Code for design and simulation of satellite laser guide stars.

The majority of this code was written by Dr. Jim Clark of the Massachusetts Institute of Technology during the course of his doctoral thesis research ("[Space-based laser guide stars for astronomical observatories](http://dspace.mit.edu/handle/1721.1/129146)").  Its copyright belongs to MIT, but MIT has authorized the release of this code under the BSD 2-clause license (case no. 23070).

(The exception is ```hamiltonian.m```, which is copyright 2015 Pramit Biswas and released [via MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/51610-hamiltonian-graph-source-destination) under the BSD 2-clause license, included here as ```license.txt``` -- not to be confused with ```LICENSE.md```, which is the license governing the rest of the code.)

## Usage

To use this software package, download all code and data files into a single folder.  MATLAB is required, as well as the Mapping Toolbox, Optimization Toolbox, and Phased Array System Toolbox.  Run ```LGSmain.m``` to initialize all variables (e.g. propulsion system options, host telescope parameters, observation mission parameters) and run all other code modules.  Several figures will be produced in the course of operation, most of which are (perhaps with some aesthetic edits) included in Clark 2020.

## Function reference

### ```LGSmain```

The main program for setting up the LGS parameters (mass, LUVOIR size, propulsion system options, *etc.*) and invoking the other functions.

### ```linkbudgetG```

A helper function for calculating Gaussian beam link budgets (power received, photons received per second, apparent magnitude, and beamwidth).

### L2 orbit calculations

#### ```OrbitCalcs2dome3```

Calculates the thrust required to hold the Telescope-LGS formation to observe any direction over the course of a full six-month halo orbit period.

#### ```OrbitCalcs2dome2```

Calculates a slice of observations and produces a contour map of thrust vs. time and azimuth (constant elevation).

#### ```OrbitCalcs2dome```

Calculates a slice of observations and produces a contour map of thrust vs. elevation and azimuth (at one time of observation).

#### ```OrbitCalcs3```

Calculates the time and delta-V cost of deploying a LGS from LUVOIR.

#### ```OrbitCalcs3CL3```

Simulates the formation flight of a laser guide star with a telescope, including closed-loop control, for evaluating drift and noise sensitivity.

#### ```cr3bpse```

A differential-equation representation of the circular restricted 3-body problem, for use with ```ode45```.  Used by the various OrbitCalcs scripts.

#### ```cr3bpsepropCLazel3```

An advanced version of ```cr3bpse``` which includes a rudimentary controller, for investigating formation flight.

#### ```cr3bpsepropCLazel3rpt```

A variant of ```cr3bpsepropCLazel3``` which includes additional return variables; it is not to be used with ```ode45```, but can be used on the trajectory coming out of ```ode45``` to see what the controller in ```cr3bpsepropCLazel3``` was doing.

### LGS performance

#### ```NoiseCalcsPropSens```

Calculates the sensitivity of the formation to thruster noise and calculates the required check-in frequency from the telescope during an observation.

#### ```NoiseCalcs```

Calculates the possible negative effects of the LGS spacecraft on the telescope's observation (thermal emission, sun glinting, *etc.*).

#### ```PlancksLaw```

An implementation of Planck's Law.

#### ```PowerCalcs```

Calculations of power and ADCS subsystem performance and requirements during formation flight.

#### ```LGSretro```

Calculations of the brightness of a passive "laser guide star" spacecraft that reflects a laser generated externally (LUVOIR itself, or a ground-based facility).

### Design Reference Mission (DRM)

#### ```StarkSkymap```

Ingests lists of stars and observations for later use, and plots them on a map.

#### ```DRM_prop_options```

Calculates the number of LGS spacecraft required to support a mission, evaluating the different propulsion system options considered in ```LGSmain```.

#### ```DRM_sensitivity```

Conducts sensitivity analyses of the design reference mission with respect to various parameters (range to the telescope, total spacecraft mass, *etc.*).

#### ```DRMfunc```

The main function for evaluating a design reference mission, invoked by ```DRM_prop_options``` and ```DRM_sensitivity```.

#### ```StarkSkymap_TSP```

Uses MATLAB's ```intlinprog``` solver to solve the Traveling Salesman Problem for the list of star targets.  This program takes longer to run than all the rest put together, so it saves its result for later reuse.

#### ```ham_StarkSkymap```

Uses Pramit Biswas's Hamiltonian code to process the graph produced by ```StarkSkymap_TSP``` into a Hamiltonian route, *i.e.* a list of stars visited in order.

#### ```StarkSchedule```

A simple greedy scheduler, for assigning observations to LGS spacecraft.  This code rigidly adheres to the schedule of observations in order of scientific priority (*i.e.* expected exo-Earth yield).

#### ```StarkScheduleAltB```

An attempted evolution of ```StarkSchedule``` that only assigns `nearby' observations to LGS spacecraft, to keep their movements segmented.  Unfortunately it doesn't yet ensure that the segments are butted against each other, so it actually performs worse (*i.e.* assigns more LGS spacecraft) than ```StarkSchedule```.

#### ```StarkScheduleAltD```

A scheduler that ingests the TSP solution and segments it among LGS spacecraft.

#### ```seed_tsp_stars```

Uses Joseph Kirk's [TSP Genetic Algorithm code](https://github.com/rubikscubeguy/matlab-tsp-ga) to solve for an optimal path among the nodes in a Fibonacci spiral, which MATLAB's ```intlinprog``` did not handle very well (due to the narrow "dynamic range" between the optimal solution and similar "good enough" solutions).  Kirk's code also includes routines for solving the multi-salesman problem, which will be investigated in the future for further optimization of the LGS scheduler.  Not currently invoked by the main LGS code.

### Pathfinder

#### ```SkyCalcs```

Calculates the line-of-sight and targets accessible from various ground telescopes through satellites in geostationary orbit.

#### ```SkyCalcsOffGEO```

Calculates the delta-V cost to incline an LGS satellite from GEO to access non-equatorial stars.

#### ```SkyCalcsOffGEO2```

Calculates the line-of-sight and targets accessible from various ground telescopes through a satellite in an inclined geosynchronous orbit.

#### ```SkyCalcsHEO```

Calculates the line-of-sight from a ground telescope through a satellite in a highly elliptical orbit.

#### ```SkyCalcsHEO2```

Calculates the line-of-sight from a telescope in geostationary orbit through a satellite in a super-geostationary orbit.

#### ```ecc_from_mean```

A helper function for calculating the eccentric anomaly from the mean anomaly.

#### ```HEO_LGS```

Calculations of highly-elliptical "sidereal" orbits for Earth-orbiting laser guide stars.

### ```hamiltonian```

Code for computing a Hamiltonian path from a graph, by Pramit Biswas.  Used according to the terms of the 2-clause BSD license.

## Data files 

### ```simbad-trim.csv```

The list of stars which are baselined to be observed by LUVOIR, per [Stark *et al.* 2015](https://iopscience.iop.org/article/10.1088/0004-637X/808/2/149), received from Chris Stark via personal communication and released with his permission.  Each row consists of the star ID (Hipparcos Catalogue), followed by the hour angle and declination, as sourced from [SIMBAD](http://simbad.u-strasbg.fr/simbad/).

### ```bright_stars_simbad_trim.csv```

The celestial coordinates (hour angle and declination) of the stars of apparent magnitude 2 and brighter, as sourced from SIMBAD.

### ```LUVOIR-Architecture_A-NOMINAL_OCCRATES-observations-trim.csv```

The list of observations baselined for LUVOIR, per Stark *et al.* 2015, received from Chris Stark via personal communication and released with his permission.  Each row consists of the star ID (Hipparcos Catalogue), the count of how many times the star has been observed (including the current one, *i.e.* starting at 1), the desired number of years since the first observation of that star, and the duration of the observation in days.
