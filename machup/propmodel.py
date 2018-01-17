import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate
import math as m


class PropModel:
    """Descritizes prop geometry information and performs blade element theory analysis.
       
    Parameters
    ----------
    machup.Plane.Prop
        Prop object that contains all of the necessary information
        about propeller geometry.

    Returns
    -------
    PropGrid
        The newly created PropModel object.

    """
    def __init__(self, prop, cg_location,units):
        self._prop = prop
        self._name = prop.name
        self._position = prop.get_position()
        self._calc_spacing()
        self._calc_prop_forward()
        self._calc_coefficients()
        self._calc_pitch()
        self._calc_chord()
        self._calc_motor()
        self._init_speed_vars()
        self._cg_location = cg_location
        self._units = units

    #---------------------------------------Functions for setting up variables based on geometry-----------------------------------------

    def _calc_spacing(self):
        self._d = self._prop.get_diameter()
        hd = self._prop.get_hub_diameter()
        self._nodes = self._prop.get_number_of_nodes()
        self._zeta = np.linspace(hd/self._d,0.999,self._nodes)
        self._r = self._zeta*self._d*0.5

    def _calc_prop_forward(self):
        self._theta,self._gamma = np.radians(self._prop.get_orientation())
        st = m.sin(self._theta)
        ct = m.cos(self._theta)
        sg = m.sin(self._gamma)
        cg = m.cos(self._gamma)
        self._prop_forward = np.array([ct*cg,ct*sg,-st])
        self._rot_dir = self._prop.get_rot_dir()

    def _calc_coefficients(self):
        r_airfoil, t_airfoil = self._prop.get_airfoils()
        self._alpha_L0 = self._linearInterp(r_airfoil.get_zero_lift_alpha(),t_airfoil.get_zero_lift_alpha(),self._zeta)
        self._CL_alpha = self._linearInterp(r_airfoil.get_lift_slope(),t_airfoil.get_lift_slope(),self._zeta)
        self._CL_max = self._linearInterp(r_airfoil.get_max_lift(),t_airfoil.get_max_lift(),self._zeta)
        self._stall_deg = np.degrees(self._CL_max/self._CL_alpha)

        rCD0,rCDL,rCDL2 = r_airfoil.get_drag_coefficients()
        tCD0,tCDL,tCDL2 = t_airfoil.get_drag_coefficients()
        self._CD_0 = self._linearInterp(rCD0,tCD0,self._zeta)
        self._CD_L = self._linearInterp(rCDL,tCDL,self._zeta)
        self._CD_L2 = self._linearInterp(rCDL2,tCDL2,self._zeta)

    def _calc_pitch(self):
        pitch_info = self._prop.get_pitch_info()
        
        if pitch_info["type"] == "default":
            Kc = 0.5
        elif pitch_info["type"] == "from_file":
            filename = pitch_info["filename"]
            r,p = self._import_data(filename)
            p = np.radians(p)
            pitch_func = interp1d(r,p,kind='cubic')
            self._beta = pitch_func(self._zeta)
            return
        elif pitch_info["type"] == "ratio":
            Kc = pitch_info["value"]
        elif pitch_info["type"] == "value":
            Kc = pitch_info["value"]/self._d
        else:
            raise RuntimeError("Invalid pitch type")

        K = np.pi*self._zeta*((Kc-np.pi*self._zeta*np.tan(self._alpha_L0))/(np.pi*self._zeta+Kc*np.tan(self._alpha_L0)))
        self._beta = np.arctan2(K,np.pi*self._zeta)

    def _calc_chord(self):
        chord_info = self._prop.get_chord_info()
        
        if chord_info["type"] == "default":
            self._cb = 0.075*self._d*np.sqrt(1-self._zeta**2)

        elif chord_info["type"] == "from_file":
            filename = chord_info["filename"]
            r,c = self._import_data(filename)
            chord_func = interp1d(r,c,kind='cubic')
            self._cb = chord_func(self._zeta)

        elif chord_info["type"] == "linear":
            root = chord_info["root"]
            tip = chord_info.get("tip",root)
            self._cb = self._linearInterp(root,tip,self._zeta)

        elif chord_info["type"] == "elliptical":
            self._cb = chord_info["root"]*np.sqrt(1-self._zeta**2)
        else:
            raise RuntimeError("Invalid chord type")

        self._k = self._prop.get_num_of_blades()
        self._chatb = (self._k*self._cb)/self._d

    def _calc_motor(self):
        motor_info = self._prop.get_motor_info()
        self._has_motor = motor_info.get("has_motor",False)
        if self._has_motor:
            self._Kv = motor_info["motor_kv"]
            self._Rm = motor_info["motor_resistance"]
            self._Io = motor_info["motor_no_load_current"]
            self._Gr = motor_info["gear_reduction"]
            self._Rb = motor_info["battery_resistance"]
            self._Eo = motor_info["battery_voltage"]
            self._Rc = motor_info["speed_control_resistance"]

    def _init_speed_vars(self):
        self._ei = np.zeros(self._nodes)
        self._rot_per_sec = 15


    #---------------------------------------Functions for setting prop operating state-----------------------------------------
    def _set_aero_state(self, aero_state):
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
            self._v_xyz = v_xyz


            if "V_mag" not in aero_state:
                raise RuntimeError("Must supply 'V_mag' key and value")
            if "rho" not in aero_state:
                raise RuntimeError("Must supply 'rho' key and value")
            self._rho = aero_state["rho"]
            self._Vinf = aero_state["V_mag"]
        else:
            print("'aero_state' not provided. Using default values: V=10 ft/s, rho = 0.0023769 slugs/ft^3")
            self._rho = 0.0023769
            self._Vinf = 10
            self._v_xyz = np.array([-1,0,0])

    def _set_prop_state(self, prop_state):
        if prop_state:
            self._solve_from_power = False
            self._solve_from_motor = False
            if "J" in prop_state:
                self._J = prop_state["J"]
                if self._J == 0. or self._Vinf == 0.:
                    if "rot_per_sec" in prop_state:
                        self._rot_per_sec = prop_state["rot_per_sec"]
                    else:
                        raise RuntimeError("For advance ratio of zero, must also supply 'rot_per_sec' key and value")
                else:
                    self._rot_per_sec = self._Vinf/(self._d*self._J)
            elif "rot_per_sec" in prop_state:
                self._rot_per_sec = prop_state["rot_per_sec"]
                self._J = self._Vinf/(self._rot_per_sec*self._d)
            elif "brake_power" in prop_state:
                self._brake_power = prop_state["brake_power"]
                self._solve_from_power = True
            elif "electric_throttle" in prop_state:
                self._throttle = prop_state["electric_throttle"]
                self._solve_from_power = True
                self._solve_from_motor = True

                
            else:
                raise RuntimeError("Must supply 'J', 'rot_per_sec', or 'brake_power' key and value")
        else:
            #print("'prop_state' not provided. Using default values: J=0.5")
            self._J = 0.5
            self._rot_per_sec = self._Vinf/(self._d*self._J)
            self._solve_from_power = False
            self._solve_from_motor = False

    #---------------------------------------Functions for calculating prop performance-----------------------------------------

    def _find_coefficients(self):
        self._einf = np.arctan2(self._J, (m.pi*self._zeta))

        
        #initial guess for iterative solution of ei if not previously solved
        if not np.any(self._ei):
            self._ei = self._beta-self._einf

        #set constant values for ei solver
        self._ei_func_a = (self._chatb/(8*self._zeta))
        self._ei_func_b = np.arccos(np.exp(-(self._k*(1-self._zeta))/(2*np.sin(self._beta[-1]))))

        #solve for ei
        self._ei = self._find_ei()
        self._alpha = self._beta-self._einf-self._ei


        CLift,CLift_a,CDrag = self.get_coefficients_alpha(self._alpha)	
        self._CL = CLift
        self._CD = CDrag

        cei = np.cos(self._ei)
        ceinf = np.cos(self._einf)
        cee = np.cos(self._einf+self._ei)
        see = np.sin(self._einf+self._ei)
        teinf = np.tan(self._einf)

        z2 = self._zeta*self._zeta
        z3 = z2*self._zeta
        cei2 = cei*cei
        ceinf2 = ceinf*ceinf
        teinf2 = teinf*teinf
        pi2 = m.pi*m.pi

        CT_z = z2*self._chatb*((cei2)/(ceinf2))*(CLift*cee-CDrag*see)
        Cl_z = z3*self._chatb*((cei2)/(ceinf2))*(CDrag*cee+CLift*see)
        CN_z = (z2)*self._chatb*(cei2)*(2*teinf*(CDrag*cee+CLift*see)-(teinf2)*(CLift*cee-(CLift_a+CDrag)*see))
        Cn_z = (z3)*self._chatb*(cei2)*(2*teinf*(CLift*cee-CDrag*see)+(teinf2)*((CLift_a+CDrag)*cee+CLift*see))

        self._CT = ((pi2)/4)*integrate.simps(CT_z,self._zeta)
        self._Cl = ((pi2)/8)*integrate.simps(Cl_z,self._zeta)
        self._CN_a = ((pi2)/8)*integrate.simps(CN_z,self._zeta)
        self._Cn_a = ((pi2)/16)*integrate.simps(Cn_z,self._zeta)

        self._CP = self._Cl*2*m.pi

        self._efficiency = (self._CT*self._J)/self._CP
        if self._CT<=0.:
            self._efficiency = 0.

    def _find_ei(self,error = 1E-11):
        window = int(self._zeta.size/7)
        x0 = np.copy(self._ei)
        x1 = np.copy(x0)*(1+1e-4)-1e-4
        y0 = self._ei_function(x0)
        y1 = self._ei_function(x1)

        x2 = np.copy(x0)
        i = 0
        while True:
            loc = np.where(abs(y1-y0)>error)
            x2[loc] = x1[loc]-((y1[loc]*(x1[loc]-x0[loc]))/(y1[loc]-y0[loc]))
            i+=1

            x2 = np.where(x2<np.radians(60),x2,0.1)
            x2 = np.where(x2>np.radians(-30),x2,0.1)
            ma = self._moving_average(x2,window)
            ma[:int(window/2)+1] = x2[:int(window/2)+1]
            ma[-int(window/2):] = x2[-int(window/2):]
            x2 = np.where(abs(x2-ma)>0.05,ma,x2)
            if loc[0].size == 0:
                return x2

            x0 = np.copy(x1)
            x1 = np.copy(x2)
            y0 = np.copy(y1)
            y1 = self._ei_function(x1)

    def _ei_function(self,ei):
        alpha = self._beta-self._einf-ei
        CLift = self._get_CL(alpha)
        zero = (self._ei_func_a*CLift)-self._ei_func_b*np.tan(ei)*np.sin(self._einf+ei)
        return zero

    def _moving_average(self,values,window):
        weights = np.repeat(1.0,window)/window
        sma = np.convolve(values,weights,'same')
        return sma
    #--------------------Functions for finding prop_state based on input power------------------------
    def _prop_state_from_power(self):
        if self._solve_from_motor:
            if self._has_motor:
                self._rot_per_sec = self._find_rot_rate_motor()
            else:
                raise RuntimeError("Propeller does not have electric motor information")
        else:
            self._rot_per_sec = self._find_rot_rate_power()
        self._J = self._Vinf/(self._rot_per_sec*self._d)

    def _find_rot_rate_power(self, error = 1E-14):
        x0 = self._rot_per_sec #get a better guess
        x1 = x0*1.5
        y0 = self._zero_power_difference(x0)
        if abs(y0)<error:
            return x0
        y1 = self._zero_power_difference(x1)
        while True:
            #new guess based on secant method
            x2 = x1-((y1*(x1-x0))/(y1-y0))
            #if solution has converged, return results
            if abs(x2-x1)<error:
                return x2
            #set variables for next iteration
            x0=x1
            y0=y1
            x1=x2
            y1=self._zero_power_difference(x1)

    def _zero_power_difference(self,rot_per_sec_est):
        self._rot_per_sec = rot_per_sec_est
        self._J = self._Vinf/(self._rot_per_sec*self._d)
        self._find_coefficients()
        prop_power = (self._rho*(self._rot_per_sec**3)*(self._d**5)*self._CP)
        return prop_power-self._brake_power

    def _find_rot_rate_motor(self,error = 1E-14):
        x0 = self._rot_per_sec #get a better guess
        x1 = x0*1.5
        y0 = self._zero_torque_difference(x0)
        if abs(y0)<error:
            return x0
        y1 = self._zero_torque_difference(x1)
        while True:
            #new guess based on secant method
            x2 = x1-((y1*(x1-x0))/(y1-y0))
            #if solution has converged, return results
            if abs(x2-x1)<error:
                return x2
            #set variables for next iteration
            x0=x1
            y0=y1
            x1=x2
            y1=self._zero_torque_difference(x1)

    def _zero_torque_difference(self,rot_per_sec_est):
        self._rot_per_sec = rot_per_sec_est
        self._J = self._Vinf/(self._rot_per_sec*self._d)
        self._find_coefficients()
        prop_torque = (self._rho*(self._rot_per_sec**2)*(self._d**5)*self._Cl)
        motor_torque = self._electric_motor_torque()
        return prop_torque-motor_torque

    def _electric_motor_torque(self):
        N = self._rot_per_sec*60
        ns = 1-0.079*(1-self._throttle)
        Im = (ns*self._throttle*self._Eo-(self._Gr/self._Kv)*N)/(ns*self._throttle*self._Rb+self._Rc+self._Rm)
        Tq = (7.0432*self._Gr/self._Kv)*(Im-self._Io)
        if self._units == None:
            raise RuntimeError("Units must be provided for electric motor analysis. Please specify units as either 'English' or 'SI'")
        if self._units == "English":
            return Tq
        else: 
            return Tq*1.3558179483314

    #---------------------------------------Functions for calculating prop forces-----------------------------------------

    def _calc_forces(self):
        if self._solve_from_power == True:
            self._prop_state_from_power()
        else:
            self._find_coefficients()

        angle = self._vec_angle(self._prop_forward,-self._v_xyz)
        if angle == 0.:
            N,n=0.,0.
            N_vec = np.zeros(3)
        else:
            N = self._rho*(self._rot_per_sec*self._rot_per_sec)*(self._d**4)*self._CN_a*angle
            n = self._rho*(self._rot_per_sec*self._rot_per_sec)*(self._d**5)*self._Cn_a*angle
            N_dir = np.cross(np.cross(-self._v_xyz, self._prop_forward),self._prop_forward)
            N_vec = N_dir/m.sqrt(N_dir[0]**2+N_dir[1]**2+N_dir[2]**2)

        self._T = self._rho*(self._rot_per_sec*self._rot_per_sec)*(self._d**4)*self._CT
        self._l = self._rho*(self._rot_per_sec*self._rot_per_sec)*(self._d**5)*self._Cl*(-self._rot_dir)
        self._F = (self._T*self._prop_forward)+(N*N_vec)
        x,y,z = self._position-self._cg_location
        thrust_offset_M = np.array([self._F[2]*y-self._F[1]*z,
                                    self._F[0]*z-self._F[2]*x,
                                    self._F[1]*x-self._F[0]*y])
        self._M = thrust_offset_M+(self._rot_dir*n*N_vec)+self._l*self._prop_forward
        return self._F, self._M


    #---------------------------------------Functions for calculating prop slipstream effects on control points-----------------------------------------

    def _find_points_in_stream(self,control_points):
        #determine which control points fall in cylinder behind prop and their respective distances along the streamtube and from the center of the streamtube
        prop_pos = self._position
        self._t = np.dot((control_points-prop_pos),self._prop_forward)#distance down centerline to position of min distance to point. negative is behind prop, positive is in front of prop
        stream_vec = self._prop_forward*self._t[:, np.newaxis] #vector down centerline to position of min distance to point
        P = prop_pos+stream_vec #position on center line of min distance to point
        center2point_vec = control_points-P #vector from point to centerline
        self._center2point_dist = self._vec_mag_mult(center2point_vec) #distance from point to center line
        tangent = np.cross(self._prop_forward,center2point_vec)
        self._tangent_norm = tangent/self._vec_mag_mult(tangent)[:,None]

        self._points_in_stream = np.where((abs(self._center2point_dist)<(self._d/2))&(self._t<0.))

    def _calc_velocities(self, control_points):
        Vx_mag, Vt_mag = self._find_slipstream_mult(self._points_in_stream)

        tangential_vel = np.zeros_like(control_points)
        axial_vel = np.zeros_like(control_points)
        tangential_vel[self._points_in_stream] = self._tangent_norm[self._points_in_stream]*Vt_mag[:,None]
        axial_vel[self._points_in_stream] = -self._prop_forward*Vx_mag[:,None]
        total_vel = tangential_vel+axial_vel
        return total_vel

    def _find_slipstream_mult(self,loc):
        s = -self._t[loc]
        Vxi,Vti = self._find_Vi()
        if s.size == 0:
            return s,s

        kd = self._find_ss_dev_mult(s)
        r_prime = np.empty([self._r.size,s.size])
        r_prime[0] = 0.00 #nacelle radius

        Vx2 = Vxi[1:,None]+Vxi[0:-1,None]
        Kv = (2*self._Vinf+Vx2)/(2*self._Vinf+kd[None,:]*Vx2)

        for i in range(1,self._r.size):
            r_prime[i] = np.sqrt((r_prime[i-1]*r_prime[i-1])+((self._r[i]*self._r[i])-(self._r[i-1]*self._r[i-1]))*Kv[i-1])


        Vxim = Vxi[:,None]*kd
        Vtim = np.empty_like(Vxim)
        Vtim[int(self._r.size/10):][:] = 2*Vti[int(self._r.size/10):,None]*self._r[int(self._r.size/10):,None]/r_prime[int(self._r.size/10):][:]
        Vtim[0:int(self._r.size/10)][:] = (Vtim[int(self._r.size/10)][:]/r_prime[int(self._r.size/10)][:])*r_prime[0:int(self._r.size/10)][:]

        Vx_point = np.empty_like(s)
        Vt_point = np.empty_like(s)

        in_straight_tube = self._center2point_dist[loc]

        for i in range(0,s.size):
            Vx_point[i] = np.interp(in_straight_tube[i],r_prime[:,i],Vxim[:,i], left=0,right=0)
            Vt_point[i] = np.interp(in_straight_tube[i],r_prime[:,i],Vtim[:,i], left=0,right=0)*self._rot_dir

        return Vx_point, Vt_point

    def _find_Vi(self):
        self._Vi = self._rot_per_sec*2*m.pi*self._r*np.sin(self._ei)/np.cos(self._einf)
        eb = self._ei+self._einf
        self._Vxi = self._Vi*np.cos(eb)
        self._Vti = self._Vi*np.sin(eb)
        return self._Vxi,self._Vti

    def _find_ss_dev_mult(self,s):
        R = self._d/2
        return 1+(s/np.sqrt((s*s)+(R*R)))

    #---------------------------------------Functions for calculating lift and drag coefficients-----------------------------------------

    def get_coefficients_alpha(self, alpha):
        CL = self._get_CL(alpha)
        CL_a = self._get_CL_a(alpha)
        CD = self._get_CD(CL)
        return CL,CL_a,CD

    def _get_CL(self,alpha):
        neg = np.where(alpha<0.)
        pos = np.where(alpha>0.)
        CL = np.empty_like(alpha)
        CL[neg] = -self._lift_curve(-alpha,neg)
        CL[pos] = self._lift_curve(alpha,pos)
        return CL

    def _lift_curve(self,alpha,loc):
        step = (np.arctan(np.degrees(alpha[loc])-self._stall_deg[loc]-3.0)/np.pi+0.5)
        linear = self._CL_alpha[loc]*np.sin(alpha[loc])
        pressure = 2.2*np.sin(alpha[loc])*np.cos(alpha[loc])
        cL = linear*(1-step)+pressure*step
        return cL

    def _get_CL_a(self,alpha):
        neg = np.where(alpha<0.)
        pos = np.where(alpha>0.)
        CL_a = np.empty_like(alpha)
        CL_a[neg] = self._lift_slope_curve(-alpha,neg)
        CL_a[pos] = self._lift_slope_curve(alpha,pos)
        return CL_a

    def _lift_slope_curve(self,alpha,loc):
        step = (np.arctan(np.degrees(alpha[loc])-self._stall_deg[loc]-3.0)/np.pi+0.5)
        linear = self._CL_alpha[loc]*np.sin(alpha[loc])
        pressure = 2.2*np.sin(alpha[loc])*np.cos(alpha[loc])
        s_a = (180/(np.pi*np.pi))*(1/(1+(np.degrees(alpha[loc])-self._stall_deg[loc]-3.0)**2))
        l_a = self._CL_alpha[loc]
        p_a = 2.2*np.cos(2*alpha[loc])

        cL_a = -linear*s_a+(1-step)*l_a+pressure*s_a+step*p_a
        return cL_a

    def _get_CD(self,CL):
        return self._CD_0+self._CD_L*CL+self._CD_L2*CL*CL

    #-------------------------------get characteristics functions---------------------------------

    def get_name(self):
        """Get the name of the propeller. 

        Returns
        -------
        string
            name of the propeller

        """
        return self._name

    def get_num_nodes(self):
        """Get the number of nodes in the propeller. 

        Returns
        -------
        int
            number of nodes used in blade element analysis

        """
        return self._nodes

    def get_num_blades(self):
        """Get the number of blades in the propeller. 

        Returns
        -------
        int
            number of blades in the propeller

        """
        return self._k
    
    def get_position(self):
        """Get the position of the propeller. 

        Returns
        -------
        Array
            Array containing the position [x,y,z] of the propeller

        """
        return self._position

    def get_orientation(self):
        """Get the orientation of the propeller. 

        Returns
        -------
        float
            elevation angle of propeller
        float
            heading angle of propeller

        """
        return self._theta, self._gamma

    def get_rotation_direction(self):
        """Get the rotation direction of the propeller. 

        Returns
        -------
        int
            Rotation direction of the propeller. 1 indicates CCW rotation, 
            0 indicates CW rotation.

        """
        return self._rot_dir

    def get_radial_position(self):
        """Get the radial distance to the control points of the propeller. 

        Returns
        -------
        array
            Array containing the radial distance to each of the control points

        """
        return self._r
    
    def get_chord(self):
        """Get chord of the propeller. 

        Returns
        -------
        array
            Array containing the chord of the propeller at each of the control points

        """
        return self._cb

    def get_aero_pitch_angle(self):
        """Get pitch angle of the propeller. 

        Returns
        -------
        array
            Array containing the pitch angle of the propeller at each of the control points

        """
        return self._beta

    def get_advance_angle(self):
        """Get downwash angle due to forward movement of the propeller. 

        Returns
        -------
        array
            Array containing the advance angle of the propeller at each of the control points

        """
        return self._einf

    def get_induced_angle(self):
        """Get induced downwash angle due to the other blades of the propeller. 

        Returns
        -------
        array
            Array containing the induced angle of the propeller at each of the control points

        """
        return self._ei

    def get_alpha(self):
        """Get aerodynamic angle of attack of the propeller. 

        Returns
        -------
        array
            Array containing the aerodynamic angle of attack of the propeller at each of the control points

        """
        return self._alpha

    def get_section_CL(self):
        """Get CL for each section of the propeller. 

        Returns
        -------
        array
            Array containing the section lift coefficient at each of the control points

        """
        return self._CL

    def get_section_CD(self):
        """Get CD for each section of the propeller. 

        Returns
        -------
        array
            Array containing the section drag coefficient at each of the control points

        """
        return self._CD

    def get_zero_lift_alpha(self):
        """Get zero lift angle of attack of propeller. 

        Returns
        -------
        array
            Array containing the zero lift angle of attack at each of the control points

        """
        return self._alpha_L0

    def get_Vi(self):
        """Get total induced velocity of the propeller. 

        Returns
        -------
        array
            Array containing the total induced velocity at each of the control points

        """
        return self._Vi

    def get_axial_Vi(self):
        """Get axial induced velocity of the propeller. 

        Returns
        -------
        array
            Array containing the axial induced velocity at each of the control points

        """
        return self._Vxi

    def get_tangential_Vi(self):
        """Get tangential induced velocity of the propeller. 

        Returns
        -------
        array
            Array containing the tangential induced velocity at each of the control points

        """
        return self._Vti


    #---------------------------------------helpful assist functions-----------------------------------------

    def _import_data(self, filename):
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]
        return x,y

    def _linearInterp(self, left, right, array):
        slope = (right-left)/(array[-1]-array[0])
        inter = left+((array-array[0])*slope)
        return inter

    def _vec_mag_mult(self,vec):
        return np.sqrt(vec[:,0]**2+vec[:,1]**2+vec[:,2]**2)

    def _vec_angle(self, v1,v2):
        a,b,c = v1
        d,e,f = v2
        dot = a*d+b*e+c*f
        mag = m.sqrt(a*a+b*b+c*c)*m.sqrt(d*d+e*e+f*f)
        return m.acos(dot/mag)
