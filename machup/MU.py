from machup import Airplane
from machup import LLModel
from machup import PropsSolver
import machup.outputstl as stl

import numpy as np
import matplotlib.pyplot as plt

class MachUp:
    """The MachUp model.

    The MachUp model combines a model for calculating the aerodynamics
    on propellers with a separate model for calculating the aerodynamics
    on lifting surfaces such as wings. The integration of these two 
    models allows a user to construct an aircraft with multiple propellers
    and wings and then analyze its aerodynamics as a whole.

    Parameters
    ----------
    inputfile : string
        The filename for a .json file input that provides all of the
        necessary information about the airplane.

    inputAirplane : machup.plane
        Plane object which provides necessary geometry information to
        lifting line algorithm.

    Returns
    -------
    MachUp
        Returns the newly created MachUp object.

    Examples
    --------
    A simple use case is shown below

    import machup.MU as MU

    filename = "SimplePlane.json"
    muAirplane = MU.MachUp(filename)
    

    control_state = {
        "aileron": 10.,
        "elevator": 5.,
        "rudder": 0.,
        "flap": 0.,
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 0.,
        "rho": 0.0023769
    }
    prop_state = {"J": 0.25}

    results = muAirplane.solve(aero_state = aero_state,
                               control_state = control_state,
                               prop_state = prop_state)

    """
    def __init__(self,inputfile=None, inputAirplane = None):
        self.myairplane = None
        if inputfile:
            self.myairplane = Airplane(inputfile=inputfile)
        elif inputAirplane:
            self.myairplane = inputAirplane

        if self.myairplane:
            if self.myairplane._wings:
                self.myllmodel = LLModel(self.myairplane)
                self.Sw = self.myllmodel._grid.get_reference_area()
                self.lat_ref = self.myairplane.get_lat_ref()
                self.long_ref = self.myairplane.get_long_ref()

                if self.lat_ref is None or self.long_ref is None:
                    print("Reference lengths not given. Assuming following reference lengths:")

                    if self.lat_ref == None and self.long_ref == None:
                        self.lat_ref = 0
                        segs = self.myairplane.get_wingsegments()
                        for seg in segs:
                            span = 2*seg.get_span()
                            if span>self.lat_ref:
                                self.lat_ref = span
                        self.long_ref = self.Sw/self.lat_ref

                    elif self.long_ref == None:
                        self.long_ref = self.Sw/self.lat_ref

                    elif self.lat_ref == None:
                        self.lat_ref = self.Sw/self.long_ref
                    print("Longitudinal Reference: ", self.long_ref)
                    print("Lateral Reference: ", self.lat_ref)
                

            if self.myairplane._props:
                self.myprops = PropsSolver(self.myairplane)

    def solve(self, aero_state = None, control_state = None, prop_state = None):
        """Solve the combined propeller and wing models to determine the
        aerodynamic forces for the provided Plane.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        results : dict
            Python dictionary containing the resulting forces and
            moments about the X, Y, and Z axis in the standard
            body-fixed coordinate system as well as the lift and 
            drag forces. Dictionary keys are "FL", "FD", "FX",
            "FY", "FZ", "MX", "MY", and "MZ".

        """
        if aero_state:
            prop_results = {}
            wing_results = {}
            total_results = {}
            if self.myairplane._props:
                if self.myairplane._wings:
                    control_points = self.myllmodel._grid.get_control_point_pos()
                else: 
                    control_points = None
                prop_results,velocity_array = self.myprops.solve(aero_state,prop_state,control_points)

                if velocity_array is not None:
                    aero_state["local_state"] = velocity_array

            if self.myairplane._wings:
                wing_results = self.myllmodel.solve("linear",aero_state = aero_state,control_state = control_state)

            total_results['FL'] = prop_results.get('FL',0)+wing_results.get('FL',0)
            total_results['FD'] = prop_results.get('FD',0)+wing_results.get('FD',0)
            total_results['FS'] = prop_results.get('FS',0)+wing_results.get('FS',0)

            total_results['FX'] = prop_results.get('FX',0)+wing_results.get('FX',0)
            total_results['FY'] = prop_results.get('FY',0)+wing_results.get('FY',0)
            total_results['FZ'] = prop_results.get('FZ',0)+wing_results.get('FZ',0)

            total_results['MX'] = prop_results.get('MX',0)+wing_results.get('MX',0)
            total_results['MY'] = prop_results.get('MY',0)+wing_results.get('MY',0)
            total_results['MZ'] = prop_results.get('MZ',0)+wing_results.get('MZ',0)

        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769,
                          "alpha": 0}
            total_results = self.solve(aero_state = aero_state,control_state = control_state, prop_state = prop_state)

        return total_results

    def create_stl(self,filename = "plane.stl"):
        """Create a .stl file showing the geometry of the MachUp object

        Parameters
        ----------
        filename : string
            Filename for the .stl file that will be created. Should end in '.stl'
            If no filename is specified, .stl file is saved under default name

        Returns
        -------
        None
        """
        self.solve()

        if self.myairplane._wings:
            if self.myairplane._props:
                stl.create_from_grid(llgrid = self.myllmodel._grid, prop_models = self.myprops.get_prop_models(),filename = filename)
            else:
                stl.create_from_grid(llgrid = self.myllmodel._grid,filename = filename)
        elif self.myairplane._props:
            stl.create_from_grid(prop_models = self.myprops.get_prop_models(),filename = filename)
        

    def distributions(self,filename = None, aero_state = None, control_state = None, prop_state = None):
        """Obtain aerodynamic information about model. Values are returned for each
           section along the wings and propellers.

        Parameters
        ----------
        filename : string
            Filename for an output .txt file containing all distributions info.
            If no filename is specified, no output file is created
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        wing_distributions : dict
            Python dictionary containing information about the analytical 
            sections of the wings. Dictionary keys are "name", "control_points", "chord",
            "twist", "dihedral", "sweep", "area", "alpha", "section_CL", 
            "section_Cm", "section_CD_parasitic", and "section_alpha_L0".
        prop_distributions : dict
            Python dictionary containing information about the analytical 
            sections of the props. Dictionary keys are "name", "radial_position", "chord",
            "pitch", "advance_angle", "induced_angle", "alpha", "section_CL", 
            "section_CD", "Vi", "Vi_x",  and "Vi_t".

        """
        self.solve(aero_state, control_state, prop_state)
        if self.myairplane._wings:
            wing_distributions = self.myllmodel._distributions()
            wing_distributions["section_CL_ref"] = wing_distributions["section_CL"]*wing_distributions["chord"]/self.long_ref
        else: 
            wing_distributions = {}
            
        if self.myairplane._props:
            prop_distributions = self.myprops._distributions()
        else:
            prop_distributions = {}

        if filename:
            print("Distributions file output not currently supported.")
            #do stuff to print out distributions info to .txt file

        return wing_distributions, prop_distributions

    def find_aero_center(self,aero_state = None, control_state = None, prop_state = None):
        """Determine the location of the aerodynamic center for the given state.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        ac_loc : numpy array
            Array containing the location of the aerodynamic center in the 
            MachUp coordinate frame

        """

        if aero_state:
            ao = aero_state.get("alpha",0.)
            DP = 0.5*aero_state.get("rho",0.0023769)*(aero_state.get("V_mag",10.)**2)*self.Sw
            Normal = []
            Axial = []
            Moment = []
            a = np.arange(-10,10,0.05)
            for alpha in a:
                aero_state["alpha"] = alpha
                FM = self.solve(aero_state = aero_state, control_state = control_state, prop_state = prop_state)
                Normal.append(-FM["FZ"])
                Axial.append(-FM["FX"])
                Moment.append(FM["MY"])

            alpha = np.radians(a)
            
            CN = np.array(Normal)/DP
            CA = np.array(Axial)/DP
            Cm = np.array(Moment)/DP/self.long_ref

            ca = np.cos(alpha)
            sa = np.sin(alpha)

            CL = CN*ca-CA*sa
            CD = CA*ca+CN*sa


            sCLca = np.sum(CL*ca)
            ss2a = np.sum(sa*sa)
            sCLsa = np.sum(CL*sa)
            ssaca = np.sum(sa*ca)
            sc2a = np.sum(ca*ca)

            aL0 = np.arctan((sCLca*ss2a-sCLsa*ssaca)/(sCLca*ssaca-sCLsa*sc2a))

            ta0 = np.tan(aL0)
            CL0_a = np.sum(CL*(sa-ta0*ca))/np.sum((sa-ta0*ca)**2)


            sCL = np.sum(CL)
            sCL2 = np.sum(CL*CL)
            sCL3 = np.sum(CL*CL*CL)
            sCL4 = np.sum(CL*CL*CL*CL)
    
            sCD = np.sum(CD)
            sCDCL = np.sum(CD*CL)
            sCDCL2 = np.sum(CD*CL*CL)

            LHS = np.array([[CL.size,sCL,sCL2],
                            [  sCL ,sCL2,sCL3],
                            [  sCL2,sCL3,sCL4]])
            LHSinv = np.linalg.inv(LHS)
            RHS = np.array([sCD,sCDCL,sCDCL2])
            CD0,CD0_L,CD0_L2 = np.matmul(LHSinv,RHS)

            s2a = np.sin(2*alpha)

            ss22a = np.sum(s2a*s2a)
            ss2aCN = np.sum(s2a*CN)
            ss2aCA = np.sum(s2a*CA)
            sCN2 = np.sum(CN*CN)
            sCNCA = np.sum(CN*CA)
            sCA2 = np.sum(CA*CA)

            sCms2a = np.sum(Cm*s2a)
            sCmCN = np.sum(Cm*CN)
            sCmCA = np.sum(Cm*CA)

            LHS = np.array([[ss22a,ss2aCN,ss2aCA],
                            [ss2aCN,sCN2,sCNCA],
                            [ss2aCA,sCNCA,sCA2]])
            LHSinv = np.linalg.inv(LHS)
            RHS = np.array([sCms2a,sCmCN,sCmCA])

            Cm0_a,Cm_N,Cm_A = np.matmul(LHSinv,RHS)

            alphao = np.radians(ao)

            Xac = -self.long_ref*(-((4+(2*CD0/CL0_a)+(12*alphao*alphao-4)*CD0_L2*CL0_a)/(2*CL0_a+3*CD0-2*CD0_L2*CL0_a*CL0_a))*Cm0_a-Cm_N)
            Zac = -self.long_ref*(((4*aL0+4*CD0_L+12*alphao*CD0_L2*CL0_a)/(2*CL0_a+3*CD0-2*CD0_L2*CL0_a*CL0_a))*Cm0_a+Cm_A)

            cg_loc = self.myllmodel._grid.get_cg_location()
            ac_loc = np.array([Xac+cg_loc[0],0.,Zac+cg_loc[2]]) 
            
            
        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769,
                          "alpha": 0}
            ac_loc = self.find_aero_center(aero_state,control_state,prop_state)
        return ac_loc


    def stall_onset(self,aero_state = None, control_state = None, prop_state = None):
        """Determine the approximate angle of attack at which the aircraft will stall
           for the given operating conditions.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        stall_info : dict
            Python dictionary containing the angle of attack at which the aircraft
            stalled and the amount of lift generated at that angle of attack. 
            Dictionary keys are "alpha" and "lift".

        """
        if self.myairplane._wings:
            if aero_state:
                a0 = aero_state.get("alpha",0.) 
                a = 0.
                while True:
                    aero_state["alpha"] = a
                    results = self.solve(aero_state,control_state,prop_state)
                    stall_info = self.myllmodel._find_stall_location()
                    if stall_info:
                        stall_info["alpha"] = round(a,2)
                        stall_info["lift"] = results["FL"]
                        return stall_info
                    a+=0.1
            else:
                aero_state = {"V_mag":10,
                              "rho":0.0023769,
                              "alpha": 0}
                stall_info = self.stall_onset(aero_state,control_state,prop_state)
                return stall_info
        else:
            raise RuntimeError("Current model does not have wings to stall")

    def stall_airspeed(self, weight, aero_state = None, control_state = None, prop_state = None, error = 10E-10):
        """Determine the minimum airspeed at which the aircraft can fly without 
           stalling. Requires that the weight of the aircraft be specified.

        Parameters
        ----------
        weight : float
            The weight of the aircraft in the correct units. (lbs in english 
            units or Newtons in metric units)
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.
        error : float
            max acceptable error for the secant method convergence. Default value
            is 10E-10.

        Returns
        -------
        min_airspeed : float
            Minimum airspeed at which the aircraft can fly without stalling.

        """
        if aero_state:
            aero0 = aero_state.copy()
            aero1 = aero_state.copy()
            x0 = aero0["V_mag"] #airspeed guess
            x1 = x0*1.5
            aero0["V_mag"] = x0
            aero1["V_mag"] = x1
            y0 = self.stall_onset(aero0,control_state,prop_state)["lift"]-weight #difference between max lift and weight
            y1 = self.stall_onset(aero1,control_state,prop_state)["lift"]-weight

            while True:
                #new guess based on secant method
                x2 = x1-((y1*(x1-x0))/(y1-y0))
                #if solution has converged, return results
                if abs(x2-x1)<error:
                    min_airspeed = x2
                    return min_airspeed
                #set variables for next iteration
                x0=x1
                y0=y1
                x1=x2
                aero1["V_mag"] = x1
                y1=self.stall_onset(aero1,control_state,prop_state)["lift"]-weight
        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769,
                          "alpha": 0}
            return self.stall_airspeed(self, weight, aero_state, control_state, prop_state)

    def derivatives(self,aero_state = None, control_state = None, prop_state = None):
        """Compute the stability, damping, and control derivatives at given state.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        results : dict
            Python dictionary containing the stability, damping, and control 
            derivatives. 

        """

        stab_deriv = self.stability_derivatives(aero_state,control_state,prop_state)
        cont_deriv = self.control_derivatives(aero_state,control_state,prop_state)
        damp_deriv = self.damping_derivatives(aero_state,control_state,prop_state)
        all_derivatives = {**stab_deriv,**cont_deriv,**damp_deriv}
        all_derivatives["static_margin"] = -all_derivatives["Cm_a"]/all_derivatives["CL_a"]
        return all_derivatives


    def stability_derivatives(self, aero_state = None, control_state = None, prop_state = None):
        """Compute the stability derivatives at given state.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        stab_deriv : dict
            Python dictionary containing the stability derivatives. 

        """
        stab_deriv = {}

        if aero_state:
            DP = 0.5*aero_state.get("rho",0.0023769)*(aero_state.get("V_mag",10.)**2)*self.Sw

            delta = 0.5
            delta_rad = np.radians(delta)
            #alpha
            a0 = aero_state.get("alpha",0.)
            aero_state["alpha"] = a0-delta
            back = self.solve(aero_state,control_state, prop_state)
            aero_state["alpha"] = a0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            stab_deriv["CL_a"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            stab_deriv["CD_a"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            stab_deriv["CS_a"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            stab_deriv["CX_a"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            stab_deriv["CY_a"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            stab_deriv["CZ_a"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            stab_deriv["Cl_a"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            stab_deriv["Cm_a"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            stab_deriv["Cn_a"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            aero_state["alpha"] = a0

            #beta
            b0 = aero_state.get("beta",0.)
            aero_state["beta"] = b0-delta
            back = self.solve(aero_state,control_state, prop_state)
            aero_state["beta"] = b0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            stab_deriv["CL_b"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            stab_deriv["CD_b"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            stab_deriv["CS_b"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            stab_deriv["CX_b"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            stab_deriv["CY_b"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            stab_deriv["CZ_b"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            stab_deriv["Cl_b"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            stab_deriv["Cm_b"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            stab_deriv["Cn_b"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            aero_state["beta"] = b0

        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769,
                          "alpha": 0}
            stab_deriv = self.stability_derivatives(aero_state,control_state,prop_state)

        return stab_deriv

    def damping_derivatives(self, aero_state = None, control_state = None, prop_state = None):
        """Compute the damping derivatives at given state.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        damp_deriv : dict
            Python dictionary containing the damping derivatives. 

        """
        damp_deriv = {}
        
        if aero_state:
            V = aero_state.get("V_mag",10.)
            DP = 0.5*aero_state.get("rho",0.0023769)*V*V*self.Sw

            delta = 0.005
            dpbar = delta
            dp = dpbar*2*V/self.lat_ref
            non_dim = 2*dpbar*DP

            p0 = aero_state.get("roll_rate",0.)
            aero_state["roll_rate"] = p0-dp
            back = self.solve(aero_state,control_state, prop_state)
            aero_state["roll_rate"] = p0+dp
            forw = self.solve(aero_state,control_state, prop_state)
            damp_deriv["CL_pbar"] = (forw["FL"]-back["FL"])/non_dim
            damp_deriv["CD_pbar"] = (forw["FD"]-back["FD"])/non_dim
            damp_deriv["CS_pbar"] = (forw["FS"]-back["FS"])/non_dim
            damp_deriv["CX_pbar"] = (forw["FX"]-back["FX"])/non_dim
            damp_deriv["CY_pbar"] = (forw["FY"]-back["FY"])/non_dim
            damp_deriv["CZ_pbar"] = (forw["FZ"]-back["FZ"])/non_dim
            damp_deriv["Cl_pbar"] = (forw["MX"]-back["MX"])/non_dim/self.lat_ref
            damp_deriv["Cm_pbar"] = (forw["MY"]-back["MY"])/non_dim/self.long_ref
            damp_deriv["Cn_pbar"] = (forw["MZ"]-back["MZ"])/non_dim/self.lat_ref
            aero_state["roll_rate"] = p0

            dqbar = delta
            dq = dqbar*2*V/self.long_ref

            q0 = aero_state.get("pitch_rate",0.)
            aero_state["pitch_rate"] = q0-dq
            back = self.solve(aero_state,control_state, prop_state)
            aero_state["pitch_rate"] = q0+dq
            forw = self.solve(aero_state,control_state, prop_state)
            damp_deriv["CL_qbar"] = (forw["FL"]-back["FL"])/non_dim
            damp_deriv["CD_qbar"] = (forw["FD"]-back["FD"])/non_dim
            damp_deriv["CS_qbar"] = (forw["FS"]-back["FS"])/non_dim
            damp_deriv["CX_qbar"] = (forw["FX"]-back["FX"])/non_dim
            damp_deriv["CY_qbar"] = (forw["FY"]-back["FY"])/non_dim
            damp_deriv["CZ_qbar"] = (forw["FZ"]-back["FZ"])/non_dim
            damp_deriv["Cl_qbar"] = (forw["MX"]-back["MX"])/non_dim/self.lat_ref
            damp_deriv["Cm_qbar"] = (forw["MY"]-back["MY"])/non_dim/self.long_ref
            damp_deriv["Cn_qbar"] = (forw["MZ"]-back["MZ"])/non_dim/self.lat_ref
            aero_state["pitch_rate"] = q0

            drbar = delta
            dr = drbar*2*V/self.lat_ref
            r0 = aero_state.get("yaw_rate",0.)
            aero_state["yaw_rate"] = r0-dr
            back = self.solve(aero_state,control_state, prop_state)
            aero_state["yaw_rate"] = r0+dr
            forw = self.solve(aero_state,control_state, prop_state)
            damp_deriv["CL_rbar"] = (forw["FL"]-back["FL"])/non_dim
            damp_deriv["CD_rbar"] = (forw["FD"]-back["FD"])/non_dim
            damp_deriv["CS_rbar"] = (forw["FS"]-back["FS"])/non_dim
            damp_deriv["CX_rbar"] = (forw["FX"]-back["FX"])/non_dim
            damp_deriv["CY_rbar"] = (forw["FY"]-back["FY"])/non_dim
            damp_deriv["CZ_rbar"] = (forw["FZ"]-back["FZ"])/non_dim
            damp_deriv["Cl_rbar"] = (forw["MX"]-back["MX"])/non_dim/self.lat_ref
            damp_deriv["Cm_rbar"] = (forw["MY"]-back["MY"])/non_dim/self.long_ref
            damp_deriv["Cn_rbar"] = (forw["MZ"]-back["MZ"])/non_dim/self.lat_ref
            aero_state["yaw_rate"] = q0

        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769}
            damp_deriv = self.damping_derivatives(aero_state,control_state,prop_state)

        return damp_deriv


    
    def control_derivatives(self, aero_state = None, control_state = None, prop_state = None):
        """Compute the control derivatives at given state.

        Parameters
        ----------
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.

        Returns
        -------
        cont_deriv : dict
            Python dictionary containing the control derivatives. 

        """
        cont_deriv = {}

        if control_state:
            DP = 0.5*aero_state.get("rho",0.0023769)*(aero_state.get("V_mag",10.)**2)*self.Sw

            delta = 0.25
            delta_rad = np.radians(delta)

            #aileron
            a0 = control_state.get("aileron",0.)
            control_state["aileron"] = a0-delta
            back = self.solve(aero_state,control_state, prop_state)
            control_state["aileron"] = a0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            cont_deriv["CL_da"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            cont_deriv["CD_da"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            cont_deriv["CS_da"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            cont_deriv["CX_da"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            cont_deriv["CY_da"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            cont_deriv["CZ_da"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            cont_deriv["Cl_da"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            cont_deriv["Cm_da"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            cont_deriv["Cn_da"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            control_state["aileron"] = a0

            #elevator
            e0 = control_state.get("elevator",0.)
            control_state["elevator"] = e0-delta
            back = self.solve(aero_state,control_state, prop_state)
            control_state["elevator"] = e0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            cont_deriv["CL_de"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            cont_deriv["CD_de"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            cont_deriv["CS_de"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            cont_deriv["CX_de"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            cont_deriv["CY_de"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            cont_deriv["CZ_de"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            cont_deriv["Cl_de"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            cont_deriv["Cm_de"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            cont_deriv["Cn_de"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            control_state["elevator"] = e0

            #rudder
            r0 = control_state.get("rudder",0.)
            control_state["rudder"] = r0-delta
            back = self.solve(aero_state,control_state, prop_state)
            control_state["rudder"] = r0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            cont_deriv["CL_dr"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            cont_deriv["CD_dr"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            cont_deriv["CS_dr"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            cont_deriv["CX_dr"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            cont_deriv["CY_dr"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            cont_deriv["CZ_dr"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            cont_deriv["Cl_dr"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            cont_deriv["Cm_dr"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            cont_deriv["Cn_dr"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            control_state["rudder"] = r0

            #flap
            f0 = control_state.get("flap",0.)
            control_state["flap"] = f0-delta
            back = self.solve(aero_state,control_state, prop_state)
            control_state["flap"] = f0+delta
            forw = self.solve(aero_state,control_state, prop_state)
            cont_deriv["CL_df"] = (forw["FL"]-back["FL"])/2/delta_rad/DP
            cont_deriv["CD_df"] = (forw["FD"]-back["FD"])/2/delta_rad/DP
            cont_deriv["CS_df"] = (forw["FS"]-back["FS"])/2/delta_rad/DP
            cont_deriv["CX_df"] = (forw["FX"]-back["FX"])/2/delta_rad/DP
            cont_deriv["CY_df"] = (forw["FY"]-back["FY"])/2/delta_rad/DP
            cont_deriv["CZ_df"] = (forw["FZ"]-back["FZ"])/2/delta_rad/DP
            cont_deriv["Cl_df"] = (forw["MX"]-back["MX"])/2/delta_rad/DP/self.lat_ref
            cont_deriv["Cm_df"] = (forw["MY"]-back["MY"])/2/delta_rad/DP/self.long_ref
            cont_deriv["Cn_df"] = (forw["MZ"]-back["MZ"])/2/delta_rad/DP/self.lat_ref
            control_state["flap"] = f0

        else: 
            control_state = {"aileron":0.,
                             "elevator":0.,
                             "rudder": 0.,
                             "flap": 0.}
            cont_deriv = self.control_derivatives(aero_state,control_state,prop_state)

        return cont_deriv

    def pitch_trim(self,L_target,m_target=0.,aero_state = None, control_state = None, prop_state = None, error = 10E-10):
        """Determine the angle of attack and elevator deflection to 
           trim the aircraft in pitch at given airspeed.

        Parameters
        ----------
        L_target : float
            Desired lift force to be generated by trimmed aircraft
        m_target : float
            Desired pitching moment to be generated by trimmed aircraft.
            If no m_target is specified, m_target is set to 0.
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            aero state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, rudder, and flap deflections in
            degrees. Dictionary keys are "aileron", "elevator",
            "rudder", and "flap". If no control_state is specified, all control
            surfaces are assumed to be at zero deflection.
        prop_state : dict
            Only necessary if MachUp model has at least one propeller.
            Contains information for setting state of propeller. This information 
            can include a rotation rate (in rotations per second), advance ratio, 
            or braking power. If one is given, the others should not be included 
            in the dictionary. 
            Dictionary keys are "J", "rot_per_sec", and
            "brake_power". If no prop_state is specified, an advance ratio of 
            J = 0.5 is assumed.
        error : float
            max acceptable error for the secant method convergence. Default value
            is 10E-10.

        Returns
        -------
        results : dict
            Python dictionary containing the resulting forces and
            moments about the X, Y, and Z axis in the standard
            body-fixed coordinate system as well as the lift and 
            drag forces. Dictionary keys are "FL", "FD", "FX",
            "FY", "FZ", "MX", "MY", and "MZ".

        """
        if aero_state and control_state:

            a, de = 0.,0.
            iters = 0.
            while True:
                aero_state["alpha"] = a
                control_state["elevator"] = de
                current = self.solve(aero_state, control_state, prop_state)
                F = np.array([current["FL"]-L_target,current["MY"]-m_target])
                if abs(F[0])<error and abs(F[1])<error:
                    trim_state = {"alpha": a,
                                  "elevator": de}
                    return trim_state
                L_a,L_de = self._lift_slope(aero_state, control_state, prop_state)
                m_a,m_de = self._moment_slope(aero_state, control_state, prop_state)
                J = np.array([[L_a, L_de],
                              [m_a, m_de]])
                Ji = np.linalg.inv(J)
                dP = np.degrees(np.matmul(Ji,-F))
                a += dP[0]
                de += dP[1]
                iters+=1
                if iters>200:
                    raise RuntimeError("Unable to trim aircraft. Consider increasing elevator size or moving center of gravity")

        elif aero_state:
            control_state = {"aileron":0.,
                             "elevator":0.,
                             "rudder": 0.,
                             "flap": 0.}
            return self.pitch_trim(L_target,m_target,aero_state,control_state,prop_state)

        elif control_state:
            aero_state = {"V_mag":10,
                          "rho":0.0023769}
            return self.pitch_trim(L_target,m_target,aero_state,control_state,prop_state)

        else:
            aero_state = {"V_mag":10,
                          "rho":0.0023769}
            control_state = {"aileron":0.,
                             "elevator":0.,
                             "rudder": 0.,
                             "flap": 0.}
            return self.pitch_trim(L_target,m_target,aero_state,control_state,prop_state)

    def _lift_slope(self,aero_state = None, control_state = None, prop_state = None):
        delta = 0.25

        #alpha
        a0 = aero_state.get("alpha",0.)
        aero_state["alpha"] = a0-delta
        back = self.solve(aero_state,control_state, prop_state)
        aero_state["alpha"] = a0+delta
        forw = self.solve(aero_state,control_state, prop_state)
        L_a = (forw["FL"]-back["FL"])/np.radians(2*delta)
        aero_state["alpha"] = a0

        #elevator
        e0 = control_state.get("elevator",0.)
        control_state["elevator"] = e0-delta
        back = self.solve(aero_state,control_state, prop_state)
        control_state["elevator"] = e0+delta
        forw = self.solve(aero_state,control_state, prop_state)
        L_de = (forw["FL"]-back["FL"])/np.radians(2*delta)
        control_state["elevator"] = e0

        return L_a,L_de

    def _moment_slope(self,aero_state = None, control_state = None, prop_state = None):
        delta = 0.25

        #alpha
        a0 = aero_state.get("alpha",0.)
        aero_state["alpha"] = a0-delta
        back = self.solve(aero_state,control_state, prop_state)
        aero_state["alpha"] = a0+delta
        forw = self.solve(aero_state,control_state, prop_state)
        m_a = (forw["MY"]-back["MY"])/np.radians(2*delta)
        aero_state["alpha"] = a0

        #elevator
        e0 = control_state.get("elevator",0.)
        control_state["elevator"] = e0-delta
        back = self.solve(aero_state,control_state, prop_state)
        control_state["elevator"] = e0+delta
        forw = self.solve(aero_state,control_state, prop_state)
        m_de = (forw["MY"]-back["MY"])/np.radians(2*delta)
        control_state["elevator"] = e0

        return m_a,m_de

