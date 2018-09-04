# MachUp
Machup is a Python 3 library for the design and analysis of 
fixed-wing aircraft. This includes things like calculating lift, 
drag, pitching moments, and stability derivatives. 

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
import machup.geometry as geom
from machup import LLModel


myplane = geom.Airplane()
myplane.cg_location(-0.29, 0., 0.25)

# add main wing and aileron control surface
mainwing = myplane.add_wing("main_wing",
                            position=[0., 0., 0.],
                            semispan=4.,
                            root_chord=1.,
                            tip_chord=0.5,
                            dihedral=3.)

mainwing.airfoil("NACA2412",
                 alpha_L0=-0.036899751,
                 CL_alpha=6.283185307,
                 Cm_L0=-0.0527,
                 Cm_alpha=-0.08,
                 CD0=0.0055,
                 CD0_L=-0.0045,
                 CD0_L2=0.01,
                 CL_max=1.4)

mainwing.control_surface(percent_span=(0.4, 0.9),
                         percent_chord=0.25,
                         mix={"aileron": 1.})

# generate lifting-line model for airplane
myllmodel = LLModel(myplane)

# solve the lifting-line model for the given condition
controls = {
    "aileron": 5.,
}
aero_state = {
    "V_mag": 65.,
    "alpha": 5.,
    "rho": 1.225
}
results = myllmodel.solve(stype="linear",
                          control_state=controls,
                          aero_state=aero_state)

print("Lift", results["FL"])
print("Drag", results["FD"])
print("Rolling Moment", results["l"])
print("Pitch Moment", results["m"])
print("Yaw Moment", results["n"])
```

## Features

* Easy to couple with other python codes (flight simulators, optimization tools, etc...)
* Fast
* Incorperates viscous effects
* Handles multiple lifting surfaces that have sweep, dihedral, and twist

## Documentation

Refer to docstrings for documentation of public methods.

## Installation

Once the package is cloned or downloaded, it can be installed by running one of the
following commands in the base directory of the package.

`pip install .'

or

`python setup.py install`

### Prerequisites

* Numpy 1.15.1 (at least, that is what it has been tested with)

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

## Support
Contact doug.hunsaker@usu.edu with any questions.

## License
This project is licensed under the MIT license. See LICENSE file for more information.
