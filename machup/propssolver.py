import numpy as np
import math as m
from machup import PropModel

class PropsSolver:

    def __init__(self, plane,wing_points):
        self._prop_models = {}
        self._total_nodes = 0
        props = plane.get_props()
        for prop in props:
            self._prop_models[prop.name] = PropModel(prop,plane.get_cg_location(),plane.get_units(),wing_points)
            self._total_nodes+= prop._nodes
        self._wing_points = wing_points
        

    def solve(self,aero_state = None,prop_state = None):
        """Solve the blade element algorithm for the provided Plane.

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
        prop_state : dict
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
            body-fixed coordinate system. Dictionary keys are "FL",
            "FD", "FS", "FX", "FY", "FZ", "MX", "MY", and "MZ".
        velocity_array : array
            If prop solver has wing_points array, returns an array containing
            containing the density and induced velocity at each point in 
            the order: [rho, Vx, Vy, Vz]. If no wing_points are specified, 
            returns None.


        """
        results = self._calc_forces(aero_state, prop_state)
        if self._wing_points is not None:
            velocity_array = self._initialize_velocity_field(aero_state,self._wing_points)
            for prop in self._prop_models:
                pm = self._prop_models[prop]
                velocities = pm._calc_prop_wing_vel()
                velocity_array[:,1:]+=velocities
            return results,velocity_array
        else:
            return results,None


    def _calc_forces(self,aero_state = None,prop_state = None):
        self._F = np.zeros(3)
        self._M = np.zeros(3)
        if aero_state:
            a = m.radians(aero_state.get("alpha",0.))
        else:
            a = 0.
        self._results = {
            "total":{}
            }
        for prop in self._prop_models:
            #check for individual propeller controls
            if prop_state:
                if prop in prop_state:
                    local_prop_state = prop_state[prop]
                else:
                    local_prop_state = prop_state
            else:
                local_prop_state = None

            pm = self._prop_models[prop]
            pm._set_aero_state(aero_state)
            pm._set_prop_state(local_prop_state)
            F,M = pm._calc_forces()
            self._F += F
            self._M += M

        lift,drag,side = self._rotate_aero_forces(self._F,aero_state)

        self._results["total"]['FL'] = lift
        self._results["total"]['FD'] = drag
        self._results["total"]['FS'] = side
        self._results["total"]['FX'] = self._F[0]
        self._results["total"]['FY'] = self._F[1]
        self._results["total"]['FZ'] = self._F[2]
        self._results["total"]['MX'] = self._M[0]
        self._results["total"]['MY'] = self._M[1]
        self._results["total"]['MZ'] = self._M[2]

        return self._results

    def _rotate_aero_forces(self,F,aero_state):
        a = m.radians(aero_state.get("alpha",0.))
        b = m.radians(aero_state.get("beta",0.))
        ca = m.cos(a)
        sa = m.sin(a)
        cb = m.cos(b)
        sb = m.sin(b)
        v_xyz = np.zeros(3)
        v_xyz[:] = 1/m.sqrt(1.-sa*sa*sb*sb)
        v_xyz[0] *= -ca*cb
        v_xyz[1] *= -ca*sb
        v_xyz[2] *= -sa*cb
        u_inf = v_xyz

        L_xyz = np.cross(u_inf,[0.,1.,0.])
        L_xyz /= np.linalg.norm(L_xyz)
        S_xyz = np.cross(L_xyz,u_inf)
        S_xyz /= np.linalg.norm(S_xyz)

        D = np.dot(F,u_inf)
        D_vec = D*u_inf
        
        L = np.dot(F,L_xyz)
        L_vec = L*L_xyz
        
        S_vec = F-L_vec-D_vec

        S = np.dot(S_vec,S_xyz)

        return L,D,S


    def _initialize_velocity_field(self,aero_state,wing_points):
        if aero_state:
            a = m.radians(aero_state.get("alpha",0.))
            b = m.radians(aero_state.get("beta",0.))
            ca = m.cos(a)
            sa = m.sin(a)
            cb = m.cos(b)
            sb = m.sin(b)
            v_xyz = np.zeros(3)
            v_xyz[:] = 1/m.sqrt(1.-sa*sa*sb*sb)
            v_xyz[0] *= -ca*cb
            v_xyz[1] *= -ca*sb
            v_xyz[2] *= -sa*cb
            u_inf = v_xyz

            if "V_mag" not in aero_state:
                raise RuntimeError("Must supply 'V_mag' key and value")
            if "rho" not in aero_state:
                raise RuntimeError("Must supply 'rho' key and value")

            rho = aero_state["rho"]
            Vinf = aero_state["V_mag"]
        else:
            print("'aero_state' not provided. Using default values: V=10 ft/s, rho = 0.0023769 slugs/ft^3")
            rho = 0.0023769
            Vinf = 10
            u_inf = np.array([-1,0,0])

        velocity_array = np.empty([wing_points.shape[0],4])
        velocity_array[:,0] = rho
        velocity_array[:,1:] = Vinf*u_inf
        return velocity_array

    def _distributions(self):
        
        self._initialize_arrays()
        index = 0
        for prop in self._prop_models:
            pm = self._prop_models[prop]
            num_nodes = pm._nodes
            cur_slice = slice(index,index+num_nodes)
            index+=num_nodes
            pm._find_Vi()

            self._data['name'][cur_slice] = pm.get_name()
            self._data['radial_position'][cur_slice] = pm.get_radial_position()
            self._data['chord'][cur_slice] = pm.get_chord()
            self._data['pitch'][cur_slice] = np.degrees(pm.get_aero_pitch_angle())
            self._data['advance_angle'][cur_slice] = np.degrees(pm.get_advance_angle())
            self._data['induced_angle'][cur_slice] = np.degrees(pm.get_induced_angle())
            self._data['alpha'][cur_slice] = np.degrees(pm.get_alpha())
            self._data['section_CL'][cur_slice] = pm.get_section_CL()
            self._data['section_CD'][cur_slice] = pm.get_section_CD()
            self._data['Vi'][cur_slice] = pm.get_Vi()
            self._data['Vi_x'][cur_slice] = pm.get_axial_Vi()
            self._data['Vi_t'][cur_slice] = pm.get_tangential_Vi()

        return self._data

    def _initialize_arrays(self):
        self._data = {
            'name': np.empty(self._total_nodes,dtype = object),
            'radial_position': np.zeros(self._total_nodes),
            'chord': np.zeros(self._total_nodes),
            'pitch': np.zeros(self._total_nodes),
            'advance_angle': np.zeros(self._total_nodes),
            'induced_angle': np.zeros(self._total_nodes),
            'alpha': np.zeros(self._total_nodes),
            'section_CL': np.zeros(self._total_nodes),
            'section_CD': np.zeros(self._total_nodes),
            'Vi': np.zeros(self._total_nodes),
            'Vi_x': np.zeros(self._total_nodes),
            'Vi_t': np.zeros(self._total_nodes)}


    def get_prop_models(self):
        """Get a list of all of the prop models.

        Returns
        -------
        list
            List of all PropModel objects in airplane.

        """
        props = []
        for prop in self._prop_models.values():
            props.append(prop)
        return props

