# MachUp
Machup is a Python 3 library for the design and analysis of 
fixed-wing aircraft. This includes things like calculating lift, 
drag, pitching moments, and stability derivatives. A web-based
version of MachUp can be found at aero.go.usu.edu.

At the heart of machup is a modern numerical lifting-line algorithm that
rapidly predicts flow over multiple lifting surfaces and can 
incorperate viscous effects. For a detailed explanation of the theory 
refer to the following sources.

W. F. Phillips and D. O. Snyder. "Modern Adaptation of Prandtl's
Classic Lifting-Line Theory", Journal of Aircraft, Vol. 37, No. 4
(2000), pp. 662-670.

W. F. Phillips, "Flow over Multiple Lifting Surfaces," Mechanics of
Flight, 2nd ed., Wiley, New Jersey, 2010, pp. 94 -107.


The following code demonstrates how machup might be used in a 
Python script:

```python
import machup.MU as MU
import numpy as np
import matplotlib.pyplot as plt


muAirplane = MU.MachUp("Example_FlyingWing.json")

aero_state = {"V_mag": 10,
              "rho": 0.0023769,
              "alpha": 0,
              "beta": 0}
prop_state = {"electric_throttle": 0.2}
control_state = {"elevator": 0}

muAirplane.create_stl(filename="FlyingWing.stl")

trim_info = muAirplane.pitch_trim(L_target=0.2,
                                  aero_state=aero_state,
                                  control_state=control_state,
                                  prop_state=prop_state)
print(trim_info)

aero_state["alpha"] = trim_info["alpha"]
control_state["elevator"] = trim_info["elevator"]

forcesandmoments = muAirplane.solve(aero_state=aero_state,
                                    control_state=control_state,
                                    prop_state=prop_state)
print(forcesandmoments)

wing_dis, prop_dis = muAirplane.distributions(aero_state=aero_state,
                                              control_state=control_state,
                                              prop_state=prop_state)

dynamic_pressure = 0.5*aero_state["rho"]*aero_state["V_mag"]**2
spanwise_station = wing_dis["control_points"][:, 1][0:160]
CZ = wing_dis["forces_xyz"][0:160, 2]/(dynamic_pressure*wing_dis["area"][0:160])
plt.plot(spanwise_station, CZ, '^')
plt.xlabel("Span Position Y (ft)")
plt.ylabel("Section CZ")
plt.show()

plt.plot(prop_dis["radial_position"], prop_dis["induced_angle"], 'r--')
plt.xlabel("Radial Position (ft)")
plt.ylabel("Induced Angle, ei")
plt.show()
```

## Features

* Scriptable
* Fast
* Incorperates viscous effects
* Handles multiple lifting surfaces that have sweep, dihedral, and twist
* Incorporates effects of prop wash on wings

## Documentation

For documentation, please refer to the docstrings. This can be done either by
using the built in python help command in the interpreter or by refering to the source.
For example, to learn more about a given function or class, type

```python
help(nameoffunctionorclass)
```
.

Note that the modules containing the class or function need to be imported first.

## Installation

The MachUp package can be installed by navigating to the root directory of the project
and using the following command. 

'pip install .'

### Prerequisites

* Python3
* Scipy
* Numpy
* numpy-stl
* Matplotlib (needed for plots used in the examples)

### Getting the Source Code

The source code can be found at [https://github.com/usuaero/MachUp](https://github.com/usuaero/MachUp)

You can either download the source as a ZIP file and extract the contents, or 
clone the MachUp repository using Git. If your system does not already have a 
version of Git installed, you will not be able to use this second option unless 
you first download and install Git. If you are unsure, you can check by typing 
`git --version` into a command prompt.

#### Downloading source as a ZIP file

1. Open a web browser and navigate to [https://github.com/usuaero/MachUp](https://github.com/usuaero/MachUp)
2. Make sure the Branch is set to `Master`
3. Click the `Clone or download` button
4. Select `Download ZIP`
5. Extract the downloaded ZIP file to a local directory on your machine

#### Cloning the Github repository

1. From the command prompt navigate to the directory where MachUp will be installed
2. `git clone https://github.com/usuaero/MachUp`

## Testing
Unit tests are implemented using the pytest module and are run using the following command.

'python3 -m pytest test/'

##Support
Contact doug.hunsaker@usu.edu with any questions.

##License
This project is licensed under the MIT license. See LICENSE file for more information. 
