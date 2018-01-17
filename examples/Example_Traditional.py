import machup.MU as MU
import numpy as np
import matplotlib.pyplot as plt


muAirplane = MU.MachUp("Example_Traditional.json")

aero_state = {"V_mag":10,
              "rho":0.0023769,
              "alpha": 0,
              "beta":0}
prop_state = {"J": 0.5}
control_state = {"elevator":0}

muAirplane.create_stl()

trim_info = muAirplane.pitch_trim(L_target = 0.5,aero_state = aero_state, control_state = control_state, prop_state = prop_state)
print(trim_info)

aero_state["alpha"] = trim_info["alpha"]
control_state["elevator"] = trim_info["elevator"]

FM = muAirplane.solve(aero_state = aero_state, control_state = control_state, prop_state = prop_state)
print(FM)

derivs = muAirplane.derivatives(aero_state = aero_state, control_state = control_state, prop_state = prop_state)
for key in derivs:
    print(key, derivs[key])
