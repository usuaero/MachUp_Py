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
