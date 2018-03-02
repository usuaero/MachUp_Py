"""Collection of classes that describe aircraft geometry.

Each geometry class groups together all of the physical information
needed to fully define a given type of geometry. This information
is used in constructing the aerodynamic models.

Geometry list
-------------
Airplane
    Defines the geometry for a single airplane.

Wing
    Defines the geometry for a wing.

WingSegment
    Defines the geometry for a straight section of wing.

Airfoil
    Defines the properties of a 2D aifoil section.

ControlSurface
    Defines the properties of a control surface.

"""

from collections import OrderedDict
import json
import numpy as np
import machup.helpers



class Airplane:
    """Defines the geometry for an airplane.

    A plane object fully describes a single airplane and basically
    consists of a set of wings that describe any wing and tail
    surfaces and the location of the center of gravity.

    Parameters
    ----------
    name : str
        The name of the airplane object.

    inputfile : str
        The filename for a .json file input that provides all of the
        necessary information about the airplane.

    Returns
    -------
    Airplane
        Returns the newly created Airplane object.

    Raises
    ------
    IOError
        If input filepath or filename is invalid.

    Examples
    --------

    """

    def __init__(self, name="", position=None, inputfile=None):
        self.name = name
        self._wings = {}
        self._props = {}
        self._cg_loc = np.zeros(3)
        self._long_ref = None
        self._lat_ref = None
        self._position = position
        self._connections = {}

        if inputfile:
            self._buildfrominputfile(inputfile)

    def _buildfrominputfile(self, filepath):
        # Construct airplane from information in .json file
        machup.helpers.check_valid_filepath(filepath)
        with open(filepath) as file:
            data = json.load(file, object_pairs_hook=OrderedDict)
            self.name = data["plane"]["name"]
            self.cg_location(data["plane"]["CGx"],
                             data["plane"]["CGy"],
                             data["plane"]["CGz"])
            self._long_ref = data["reference"]["longitudinal_length"]
            self._lat_ref = data["reference"]["lateral_length"]
            self._units = self._check_units(data["condition"].get("units",None))

            if "wings" in data:
                self._map_wing_connections(data["wings"])
                for wing_key, wing_dict in data["wings"].items():
                    self._add_wing_from_inputfile(wing_key, wing_dict)
            if "propellers" in data:
                for prop_key, prop_dict in data["propellers"].items():
                    self._add_prop_from_inputfile(prop_key, prop_dict)

    def _check_units(self, units):
        if units == "English" or units == "english": 
            return "English"
        elif units == "SI" or units == "si" or units == "Si":
            return "SI"
        elif units == None:
            return None
        else:
            raise RuntimeError(units+"is not a valid units descriptor")

    def _map_wing_connections(self, wings_dict):
        # Create a map of wing_name to parent name and connect location.
        id_map = {}
        connections = self._connections
        for wing_name, wing_dict in wings_dict.items():
            id_map[wing_dict["ID"]] = wing_name

            parent_id = wing_dict["connect"]["ID"]
            if parent_id == 0:
                connections[wing_name] = (None, 'tip')
            else:
                connections[wing_name] = (id_map[parent_id],
                                          wing_dict["connect"]["location"])

    def _add_wing_from_inputfile(self, wing_name, wing_dict):
        # parses inputfile and adds corresponding wing
        wing = self.add_wing(wing_name,
                             connect_to=self._connections[wing_name],
                             side=wing_dict["side"],
                             position=[wing_dict["connect"]["dx"],
                                       wing_dict["connect"]["dy"],
                                       wing_dict["connect"]["dz"]],
                             is_main = wing_dict["is_main"],
                             yoffset=wing_dict["connect"]["yoffset"],
                             semispan=wing_dict["span"],
                             sweep=wing_dict["sweep"],
                             dihedral=wing_dict["dihedral"],
                             mount_angle=wing_dict["mounting_angle"],
                             washout=wing_dict["washout"],
                             root_chord=wing_dict["root_chord"],
                             tip_chord=wing_dict["tip_chord"],
                             grid=wing_dict["grid"])

        # add airfoils to wing
        airfoils = list(wing_dict["airfoils"])
        if len(airfoils) > 1:
            root = wing_dict["airfoils"][airfoils[0]]["properties"]
            tip = wing_dict["airfoils"][airfoils[1]]["properties"]
            wing.airfoil(airfoils[0],
                         "root",
                         alpha_L0=root["alpha_L0"],
                         CL_alpha=root["CL_alpha"],
                         Cm_L0=root["Cm_L0"],
                         Cm_alpha=root["Cm_alpha"],
                         CD0=root["CD0"],
                         CD0_L=root["CD0_L"],
                         CD0_L2=root["CD0_L2"],
                         CL_max=root["CL_max"])
            wing.airfoil(airfoils[1],
                         "tip",
                         alpha_L0=tip["alpha_L0"],
                         CL_alpha=tip["CL_alpha"],
                         Cm_L0=tip["Cm_L0"],
                         Cm_alpha=tip["Cm_alpha"],
                         CD0=tip["CD0"],
                         CD0_L=tip["CD0_L"],
                         CD0_L2=tip["CD0_L2"],
                         CL_max=tip["CL_max"])
        else:
            root = wing_dict["airfoils"][airfoils[0]]["properties"]
            wing.airfoil(airfoils[0],
                         "both",
                         alpha_L0=root["alpha_L0"],
                         CL_alpha=root["CL_alpha"],
                         Cm_L0=root["Cm_L0"],
                         Cm_alpha=root["Cm_alpha"],
                         CD0=root["CD0"],
                         CD0_L=root["CD0_L"],
                         CD0_L2=root["CD0_L2"],
                         CL_max=root["CL_max"])

        control_dict = wing_dict["control"]

        # add control surface to wing
        if control_dict:
            wing.control_surface(percent_span=(control_dict["span_root"],
                                               control_dict["span_tip"]),
                                 percent_chord=(control_dict["chord_root"],
                                                control_dict["chord_tip"]),
                                 mix=control_dict["mix"],
                                 sealed=control_dict["is_sealed"])

    def _add_prop_from_inputfile(self, prop_name, prop_dict):
        # parses inputfile and adds corresponding wing
        prop = self.add_prop(prop_name,
                             position=[prop_dict["position"]["dx"],
                                       prop_dict["position"]["dy"],
                                       prop_dict["position"]["dz"]],
                             orientation = [prop_dict["orientation"]["elevation_angle"],
                                            prop_dict["orientation"]["heading_angle"]],
                             nodes=prop_dict["radial_nodes"],
                             diameter=prop_dict["diameter"],
                             hub_diameter=prop_dict["hub_diameter"],
                             num_of_blades=prop_dict["num_of_blades"],
                             rot_dir=prop_dict["rotation_direction"],
                             pitch_info=prop_dict["pitch"],
                             chord_info=prop_dict["chord"],
                             electric_motor = prop_dict.get("electric_motor",{}),
                             airfoils=prop_dict["airfoils"])

        
        # add airfoils to prop
        airfoils = list(prop_dict["airfoils"])
        if len(airfoils) > 1:
            root = prop_dict["airfoils"][airfoils[0]]["properties"]
            tip = prop_dict["airfoils"][airfoils[1]]["properties"]
            prop.airfoil(airfoils[0],
                         span_position = "root",
                         alpha_L0=root["alpha_L0"],
                         CL_alpha=root["CL_alpha"],
                         Cm_L0=root["Cm_L0"],
                         Cm_alpha=root["Cm_alpha"],
                         CD0=root["CD0"],
                         CD0_L=root["CD0_L"],
                         CD0_L2=root["CD0_L2"],
                         CL_max=root["CL_max"])
            prop.airfoil(airfoils[1],
                         span_position = "tip",
                         alpha_L0=tip["alpha_L0"],
                         CL_alpha=tip["CL_alpha"],
                         Cm_L0=tip["Cm_L0"],
                         Cm_alpha=tip["Cm_alpha"],
                         CD0=tip["CD0"],
                         CD0_L=tip["CD0_L"],
                         CD0_L2=tip["CD0_L2"],
                         CL_max=tip["CL_max"])
        else:
            root = prop_dict["airfoils"][airfoils[0]]["properties"]
            prop.airfoil(airfoils[0],
                         span_position = None,
                         alpha_L0=root["alpha_L0"],
                         CL_alpha=root["CL_alpha"],
                         Cm_L0=root["Cm_L0"],
                         Cm_alpha=root["Cm_alpha"],
                         CD0=root["CD0"],
                         CD0_L=root["CD0_L"],
                         CD0_L2=root["CD0_L2"],
                         CL_max=root["CL_max"])

    def get_num_sections(self):
        """Get the total number of sections of all of the wings.

        Returns
        -------
        int
            Total number of wing sections in airplane.

        """
        total_sections = 0
        for wing in self._wings.values():
            total_sections += wing.get_num_sections()
        return total_sections

    def get_wingsegments(self):
        """Get a list of all of the WingSegments in the airplane.

        Returns
        -------
        list
            List of all WingSegments in airplane.

        """
        segments = []
        for wing in self._wings.values():
            segments.extend(wing.get_wingsegments())
        return segments

    def get_props(self):
        """Get a list of all of the props in the airplane.

        Returns
        -------
        list
            List of all Props in airplane.

        """
        props = []
        for prop in self._props.values():
            props.append(prop)
        return props

    def add_wing(self, name, connect_to=(None, 'tip'), side='both', **dims):
        """Add wing to airplane.

        Parameters
        ----------
        wing : numpy.Wing
            Wing object that describes the wing to be added.

        Returns
        -------
        None

        """
        self._wings[name] = Wing(name, side, dims)

        parent_name = connect_to[0]
        connect_at = connect_to[1]

        if parent_name:
            parent = self._wings[parent_name]
        else:
            parent = self

        self._wings[name].connect_to(parent, connect_at)

        return self._wings[name]
    
    def add_prop(self, name, **dims):
        """Add prop to airplane.

        Parameters
        ----------
        name: name of the prop to be added
        **dims:
            All of properties need to specify a Propeller can be passed
            in as key word arguments. These include the following:
            position - list containing x,y,z position of prop,
            orientation - list containing elevation angle and heading angle of prop in degrees,
            nodes - number of points along radius used for calculating prop performace,
            diameter - diameter of prop,
            hub_diameter - diameter of prop hub,
            num_of_blades - number of propeller blades,
            rot_dir - direction of rotation of prop. Can be given as string ("CCW" or "CW") or int ("1" or "0"),
            pitch_type - method for defining pitch of prop. Values can be "ratio" for pitch to diameter ratio, "unit" for a given pitch value with unit of measurement, or "from_file" for .txt file containing blade angle distribution in degrees.
            pitch_info - dictionary containing necessary info for pitch_type selected
            chord_type - method for defining chord distribution of prop. Values can be "linear" for a linear interpolation between given root and tip chord values, "Elliptical" for a simple elliptical chord distribution, or "from_file" for .txt file containing chord distribution.
            chord_info - dictionary containing necessary info for chord_type selected
            electric_motor - dictionary containing necessary info for electric motor analysis
            airfoils - dicitonary containing airfoil data

        Returns
        -------
        None

        """
        self._props[name] = Prop(name, dims)

        return self._props[name]

    def get_cg_location(self):
        """Get the location of the center of gravity.

        Returns
        -------
        Numpy Array
            Cartesian coordinates of center of gravity.

        """
        return self._cg_loc

    def get_units(self):
        """Get the type of units used in calculation.
           This information is necessary for determining
           the rotation speed of an electric motor based 
           on throttle input. 

        Returns
        -------
        string
            string describing units used in calculations. Can be either "English" or "SI"

        """
        return self._units

    def cg_location(self, x_coord=None, y_coord=None, z_coord=None):
        """Set the location of the center of gravity.

        Parameters
        -------
        x_coord
            The updated x-coordinate of the center of gravity.
        y_coord
            The updated y-coordinate of the center of gravity.
        z_coord
            The updated z-coordinate of the center of gravity.

        """
        if x_coord:
            self._cg_loc[0] = x_coord
        if y_coord:
            self._cg_loc[1] = y_coord
        if z_coord:
            self._cg_loc[2] = z_coord

    def set_reference_lengths(self, lateral = None, longitudinal = None):
        """Set the reference lengths used to nondimensionalize the coefficients.

        Parameters
        -------
        lateral
            lateral reference length (typically the wingspan of the main wing)
        longitudinal
            longitudinal reference length (typically the average chord of the main wing)
            

        """
        self._lat_ref = lateral
        self._long_ref = longitudinal

    def get_position(self, at):
        """Get the position of the Airplane in the simulation.

        This allows the user to specify a location for the plane in
        the case that it is being used in some larger simulation.

        Parameters
        ----------
        at
            Required because get_position is a common method between
            all geometry objects but it is currently ignored in the
            case of an Airplane object.

        Returns
        -------
        Position
            Position of Airplane as specified by the user.
        """
        # pylint: disable=unused-argument
        if self._position:
            return self._position

        return np.array([0., 0., 0.])

    def get_long_ref(self):
        """Get the longitudinal reference length.

        Returns
        -------
        float
            longitudinal reference length.

        """
        return self._long_ref

    def get_lat_ref(self):
        """Get the lateral reference length.

        Returns
        -------
        float
            lateral reference length.

        """
        return self._lat_ref

class Wing:
    """Defines the geometry for a wing.

    Note that the main wing and horizontal and vertical stabilizers can
    all be described by the same set of information and thus do not have
    seperate classes.

    Parameters
    ----------
    wing_name : str
        Name of wing to be added. This will be used if the user wants to
        access the wing object at a later time.
    wing_dict: dict
        Python dictionary that contains all of the necessary information
        to build the wing.

    Returns
    -------
    Wing
        Returns the newly created Wing object.

    Examples:
    --------

    """

    def __init__(self, wing_name, side, dims):
        self.name = wing_name
        self._side = side
        self._left_segment = None
        self._right_segment = None

        if self._side == "both":
            self._left_segment = WingSegment("left_"+self.name,
                                             "left",
                                             dims)
            self._right_segment = WingSegment("right_"+self.name,
                                              "right",
                                              dims)
        elif self._side == "right":
            self._right_segment = WingSegment("right_"+self.name,
                                              "right",
                                              dims)
        elif self._side == "left":
            self._left_segment = WingSegment("left_"+self.name,
                                             "left",
                                             dims)
        else:
            raise RuntimeError(self._side+"side specification not recognized")

    def get_num_sections(self):
        """Get the number of sections of the wings.

        Returns
        -------
        int
            Total number of wing sections in wing.

        """
        total_sections = 0
        if self._side == "both":
            total_sections += self._left_segment.get_num_sections()
            total_sections += self._right_segment.get_num_sections()
        elif self._side == "right":
            total_sections += self._right_segment.get_num_sections()
        elif self._side == "left":
            total_sections += self._left_segment.get_num_sections()

        return total_sections

    def get_wingsegments(self):
        """Get a list of all of the WingSegments in the wing.

        Returns
        -------
        list
            List of all WingSegments in wing.

        """
        if self._side == "both":
            return [self._left_segment, self._right_segment]
        elif self._side == "right":
            return [self._right_segment]
        elif self._side == "left":
            return [self._left_segment]

    def connect_to(self, parent, at):
        """Define a connection between wing and a "parent" geometry.

        This allows for a wing to be position relative to another geometry
        to facilitate grouping of geometry together. This is done for
        convenience purposes. For example, to group tail surfaces together
        so that they can all be moved by changing the position of just one
        object instead of each individual surface.

        Parameters
        ----------
        Parent
            This can be another wing or an airplane object to connect to.

        at
            This is the location to attach to on the parent object. The
            only options currently available are "tip" which signifies a
            connection to wing tips of a parent wing.
        """
        if self._side == "both":
            self._left_segment.connect_to(parent, at)
            self._right_segment.connect_to(parent, at)
        elif self._side == "left":
            self._left_segment.connect_to(parent, at)
        elif self._side == "right":
            self._right_segment.connect_to(parent, at)

    def get_position(self, at):
        """Get the position connection location (at).

        Parameters
        ----------
        at
            The point where position is being queried. The only points
            currently available are the "left_tip" and "right_tip"
            connection points.

        Returns
        -------
        Position
            Position of connection point specified by 'at'.
        """
        if at == "left_tip":
            return self._left_segment.get_position(at)
        elif at == "right_tip":
            return self._right_segment.get_position(at)
        else:
            raise RuntimeError(at+" is not a valid connection point")

    def airfoil(self, name, end="both", **properties):
        """Update the properties of the root and/or tip airfoils.

        Parameters
        ----------
        end
            Can be set as 'both', 'root', or 'tip' to set both to be
            the same or modify each individually.
        name
            The name of the airfoil.

        **properties
            All of properties need to specify an Airfoil can be passed
            in as key word arguments. These include the following:
            alpha_L0 - zero-lift angle of attack,
            CL_alpha - The lift coefficient slope,
            Cm_alpha - The moment coefficient slope,
            Cm_L0 - The zero-lift coefficient of moment,
            CD0 - The zero-lift drag,
            CD0_L - The drag coefficient proportional to lift,
            CD0_L2 - The drag coefficient proportional to lift squared.
            CL_max - The max lift coefficient.

        Returns
        -------
        Airfoil
            The newly created Airfoil object.
        """
        airfoil = Airfoil(name, properties)

        if self._side == "both":
            self._left_segment.airfoil(end, airfoil)
            self._right_segment.airfoil(end, airfoil)
        elif self._side == "left":
            self._left_segment.airfoil(end, airfoil)
        elif self._side == "right":
            self._right_segment.airfoil(end, airfoil)

        return airfoil

    def control_surface(self, **properties):
        """Set the properties of the control surface.

        Parameters
        ----------
        mix : Dict
            Contains the ratio of control surface deflection to control input.
            For example, inputing {"elevator": 1.} would mix the control
            surface 100 percent with elevator control input.

        percent_span : tuple
            The spanwise location of the control surface root and tip
            respectively, normalized by semispan.

        percent_chord : tuple or float
            The control surface width at the control surface root and tip
            respectively, normalized by chord length. If a single value is
            passed in then width is assumed to be the same at both the root
            and tip.

        is_sealed : bool
            Boolean flag for whether the control surface hinge is sealed.

        Returns
        -------
        ControlSurface
            The newly created ControlSurface object.
        """
        control_surface = ControlSurface(properties)

        if self._side == "both":
            self._left_segment.control_surface(control_surface)
            self._right_segment.control_surface(control_surface)
        elif self._side == "left":
            self._left_segment.control_surface(control_surface)
        elif self._side == "right":
            self._right_segment.control_surface(control_surface)

        return control_surface


class WingSegment:
    """Defines the geometry for a wing segment.

    A wing segment is a subdivision of a wing and describes a straight
    portion of the wing. For example, a standard main wing might have two
    segments, one for the left hand portion of the wing and another for
    right hand portion. Another example might be a wing that has an
    inboard and an outboard segment each with a different sweep and
    dihedral.

    Parameters
    ----------
    name : str
        Name of WingSegment to be added. This will be used if the user wants to
        access the WingSegment object at a later time.
    wing_dict: dict
        Python dictionary that contains all of the necessary information
        to build the WingSegment.

    Returns
    -------
    WingSegment
        Returns the newly created WingSegment object.

    """

    def __init__(self, name, side, dims):
        self.name = name
        self._parent = None
        self._is_main = 1
        self._connect_at = "tip"
        self._side = side
        self._delta_pos = np.array([0., 0., 0.])
        self._dimensions = {
            "yoffset": 0.,
            "span": 4.,
            "root_chord": 1.,
            "tip_chord": 1.,
            "sweep": 0.,
            "dihedral": 0.,
            "mounting_angle": 0.,
            "washout": 0
        }
        self._root_airfoil = None
        self._tip_airfoil = None
        self._control_surface = ControlSurface()
        self._num_sections = 40
        self._unpack(dims)

    def _unpack(self, dims):
        delta_pos = dims.get("position", [0., 0., 0.])
        self._is_main = dims.get("is_main",1)
        self._delta_pos[0] = delta_pos[0]
        self._delta_pos[1] = delta_pos[1]
        self._delta_pos[2] = delta_pos[2]
        self._dimensions["yoffset"] = dims.get("yoffset", 0.)
        self._dimensions["span"] = dims.get("semispan",4)
        self._dimensions["root_chord"] = dims.get("root_chord",1)
        self._dimensions["tip_chord"] = dims.get("tip_chord",
                                                 self._dimensions["root_chord"])
        self._dimensions["sweep"] = dims.get("sweep", 0.)
        self._dimensions["dihedral"] = dims.get("dihedral", 0.)
        self._dimensions["mounting_angle"] = dims.get("mount_angle", 0.)
        self._dimensions["washout"] = dims.get("washout", 0.)
        self._num_sections = dims.get("grid", 40)
        '''
        if "airfoils" in dims:
            airfoils = list(dims["airfoils"].keys())
            if len(airfoils) > 1:
                self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
                self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[1]])
            else:
                self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
                self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
        else:
        '''
        self._root_airfoil = Airfoil()
        self._tip_airfoil = Airfoil()


    def connect_to(self, parent, at):
        """Define a connection between WingSegment and a "parent" geometry.

        This allows for a WingSegment to be position relative to another
        geometry to facilitate grouping of geometry together.

        Parameters
        ----------
        Parent
            This can be another wing or an airplane object to connect to.

        at
            This is the location to attach to on the parent object. The
            only options currently available are "tip" which signifies a
            connection to wing tips of a parent wing.
        """
        self._parent = parent
        self._connect_at = self._side+"_"+at

    def airfoil(self, end, airfoil):
        """Set the airfoil(s) of the wing segment.

        Parameters
        ----------
        end : str
            Whether the airfoil is being set as the root or tip airfoil.
        airfoil : Airfoil
            The airfoil to be used.
        """
        if end == "both":
            self._root_airfoil = airfoil
            self._tip_airfoil = airfoil
        elif end == "root":
            self._root_airfoil = airfoil
        elif end == "tip":
            self._tip_airfoil = airfoil

    def control_surface(self, control_surface):
        """Set the control surface for the wing segment.

        Parameters
        ----------
        ControlSurface
            The control surface to be used.
        """
        self._control_surface = control_surface

    def get_airfoils(self):
        """Get the root and tip Airfoil of the WingSegment.

        Returns
        -------
        tuple
            The root and tip Airfoils of the segment respectively.

        """
        return self._root_airfoil, self._tip_airfoil

    def get_num_sections(self):
        """Get the number of sections of the WingSegment.

        Returns
        -------
        int
            Total number of wing sections in segment.

        """
        return self._num_sections

    def get_span(self):
        """Get the span of the wing segment.

        Returns
        -------
        float
            The span of the wing segment.

        """
        return self._dimensions["span"]

    def get_position(self, side):
        """Get the position vector of a side of the WingSegment.

        Specifically, this returns the position at right or left
        wing tip at the quarter chord.

        Parameters
        ----------
        side : str
            The side of interest. ("left" or "right")

        Returns
        -------
        np.array
            The position vector of the left or right tip at the quarter chord.

        """
        if side == "left_tip":
            left_pos = np.copy(self._delta_pos)

            if self._side == "left":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                left_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                left_pos[1] -= span*np.cos(dihedral)
                left_pos[1] -= self._dimensions["yoffset"]
                left_pos[2] -= span*np.sin(dihedral)
            else:
                left_pos[1] += self._dimensions["yoffset"]

            return left_pos + self._get_parent_position()
        elif side == "right_tip":
            right_pos = np.copy(self._delta_pos)

            if self._side == "right":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                right_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                right_pos[1] += span*np.cos(dihedral)
                right_pos[1] += self._dimensions["yoffset"]
                right_pos[2] -= span*np.sin(dihedral)
            else:
                right_pos[1] -= self._dimensions["yoffset"]

            return right_pos + self._get_parent_position()
        else:
            raise RuntimeError("side specified incorrectly")

    def _get_parent_position(self):
        if self._parent:
            return self._parent.get_position(self._connect_at)
        return np.array([0., 0., 0.])

    def get_chord(self):
        """Get the root and tip chords of the wing segment.

        Returns
        -------
        tuple
            The root and tip chords of the segment respectively.

        """
        root_chord = self._dimensions["root_chord"]
        tip_chord = self._dimensions["tip_chord"]

        return root_chord, tip_chord

    def get_mounting_angle(self):
        """Get the mounting angle of the wing segment.

        Returns
        -------
        float
            The mounting angle of the segment (degrees).

        """
        return self._dimensions["mounting_angle"]

    def get_washout(self):
        """Get the washout of the wing segment.

        Returns
        -------
        float
            The washout of the segment (degrees).

        """
        return self._dimensions["washout"]

    def get_dihedral(self):
        """Get the dihedral of the wing segment.

        Returns
        -------
        float
            The dihedral of the segment (degrees).

        """
        return self._dimensions["dihedral"]

    def get_sweep(self):
        """Get the dihedral of the wing segment.

        Returns
        -------
        float
            The dihedral of the segment (degrees).

        """
        return self._dimensions["sweep"]

    def get_side(self):
        """Get the side of the airplane that the wing segment is on.

        Returns
        -------
        float
            The side of the segment.

        """
        return self._side

    def get_control_surface_span(self):
        """Get location of control surface along WingSegment.

        Position of start and end of control surface are given as a
        percentage of the WingSegment.

        Returns
        -------
        tuple
            The percent span that the control surface starts at and
            ends at starting from the left side of the WingSegment.

        """
        if self._side == "left":
            end, start = self._control_surface.get_control_span()
            end = 1. - end
            start = 1. - start
        elif self._side == "right":
            start, end = self._control_surface.get_control_span()

        return start, end

    def get_control_surface_chord(self):
        """Get the chord of the control surface.

        Chord given as a percent of the WingSegment chord.

        Returns
        -------
        tuple
            The chord on the left side of the control surface and the
            chord on the right side of the control surface.

        """
        if self._side == "left":
            right_chord, left_chord = self._control_surface.get_control_chord()
        elif self._side == "right":
            left_chord, right_chord = self._control_surface.get_control_chord()

        return left_chord, right_chord

    def get_control_mix(self):
        """Get the mixing parameters of the WingSegment.

        Returns
        -------
        tuple
            The aileron, elevator, rudder, and flap mixing respectively.

        """
        control_mix = self._control_surface.get_control_mix(self._side)
        mix_aileron = control_mix.get("aileron", 0.)
        mix_elevator = control_mix.get("elevator", 0.)
        mix_rudder = control_mix.get("rudder", 0.)
        mix_flap = control_mix.get("flap", 0.)

        return mix_aileron, mix_elevator, mix_rudder, mix_flap

    def is_control_surface_sealed(self):
        """Check if control surface is sealed.

        Returns
        -------
        boolean
            True if control surface has been specified as sealed.

        """
        return self._control_surface.get_is_control_sealed()

class Prop:
    """Defines the geometry for a propeller.

    Parameters
    ----------
    prop_name : str
        Name of wing to be added. This will be used if the user wants to
        access the wing object at a later time.
    prop_dict: dict
        Python dictionary that contains all of the necessary information
        to build the wing.

    Returns
    -------
    Prop
        Returns the newly created prop object.

    Examples:
    --------

    """
    def __init__(self, prop_name, prop_dict):
        self.name = prop_name
        self._unpack(prop_dict)



    def _unpack(self, prop_dict):
        self._position = prop_dict.get("position", [0.,0.,0.])
        self._orientation = prop_dict.get("orientation", [0.,0.])
        self._nodes = prop_dict.get("nodes", 100)
        self._diameter = prop_dict.get("diameter", 1.0)
        self._hub_diameter = prop_dict.get("hub_diameter", 0.1*self._diameter)
        self._blades = prop_dict.get("num_of_blades", 2)
        self._rot_dir = prop_dict.get("rot_dir", 1.)
        self._pitch_info = prop_dict.get("pitch_info",{"type":"default"})
        self._chord_info = prop_dict.get("chord_info",{"type":"default"})
        self._motor_info = prop_dict.get("electric_motor",{})
        '''
        self._root_airfoil = Airfoil()
        self._tip_airfoil = Airfoil()
        '''
        self._airfoils = {}

    def airfoil(self, name = 'NACA2412', **properties):
        """Update the properties of the root and/or tip airfoils.

        Parameters
        ----------
        end
            Can be set as 'both', 'root', or 'tip' to set both to be
            the same or modify each individually.
        name
            The name of the airfoil.

        **properties
            All of properties need to specify an Airfoil can be passed
            in as key word arguments. These include the following:
            alpha_L0 - zero-lift angle of attack,
            CL_alpha - The lift coefficient slope,
            Cm_alpha - The moment coefficient slope,
            Cm_L0 - The zero-lift coefficient of moment,
            CD0 - The zero-lift drag,
            CD0_L - The drag coefficient proportional to lift,
            CD0_L2 - The drag coefficient proportional to lift squared.
            CL_max - The max lift coefficient.

        Returns
        -------
        Airfoil
            The newly created Airfoil object.
        """
        airfoil = Airfoil(name, properties)
        '''
        if span_position == "both":
            self._root_airfoil = airfoil
            self._tip_airfoil = airfoil
        elif span_position == "root":
            self._root_airfoil = airfoil
        elif span_position == "tip":
            self._tip_airfoil = airfoil

        '''
        self._airfoils[name] = airfoil

        return airfoil

    def get_position(self):
        """Get the position of the prop.

        Returns
        -------
        numpy array
            array containing the x,y,z position of the propeller.

        """
        return self._position

    def get_orientation(self):
        """Get the orientation of the prop.

        Returns
        -------
        list
            list containing elevation angle and heading angle of prop (degrees).

        """
        return self._orientation
        
    def get_number_of_nodes(self):
        """Get the number of radial analysis nodes along the prop.

        Returns
        -------
        int
            number of radial nodes used for analysis.

        """
        return self._nodes

    def get_diameter(self):
        """Get the diameter of the prop.

        Returns
        -------
        float
            diameter of the propeller.

        """
        return self._diameter

    def get_hub_diameter(self):
        """Get the diameter of the prop hub.

        Returns
        -------
        float
            diameter of the propeller hub.

        """
        return self._hub_diameter

    def get_num_of_blades(self):
        """Get the number of blades on the propeller.

        Returns
        -------
        float
            number of blades on the propeller.

        """
        return self._blades

    def get_rot_dir(self):
        """Get the direction of rotation of the prop.

        Returns
        -------
        int
            1 denotes a ccw rotation, -1 denotes a cw rotation.

        """
        if self._rot_dir == "CCW" or self._rot_dir =="ccw" or self._rot_dir == 1.:
            return 1.
        elif self._rot_dir == "CW" or self._rot_dir =="cw" or self._rot_dir == -1.:
            return -1
        else: 
            return 1.

    def get_airfoils(self):
        """Get the root and tip Airfoil of the propeller.

        Returns
        -------
        tuple
            The root and tip Airfoils of the propeller respectively.

        """
        #return self._root_airfoil, self._tip_airfoil
        return self._airfoils

    def get_pitch_info(self):
        """Get the pitch info for the prop.

        Returns
        -------
        dict
            The pitch info dictionary.

        """
        return self._pitch_info

    def get_chord_info(self):
        """Get the chord info for the prop.

        Returns
        -------
        dict
            The chord info dictionary.

        """
        return self._chord_info

    def get_motor_info(self):
        """Get the electric motor info for the prop.

        Returns
        -------
        dict
            Dictionary containing necessary info for electric motor analysis.

        """
        return self._motor_info

class Airfoil:
    """Defines the aerodynamic properties of an airfoil.

    Parameters
    ----------
    airfoil_data : dict
        A python dictionary that contains all fo the properties of the
        Airfoil.

    Returns
    -------
    Airfoil
        Returns the newly created Airfoil object.

    Examples:
    --------

    """

    def __init__(self, airfoil_name='NACA2412', airfoil_data=None):
        self.name = airfoil_name
        if airfoil_data:
            self._properties = {
                "span_position": airfoil_data.get("span_position",None),
                "alpha_L0": airfoil_data["alpha_L0"],
                "CL_alpha": airfoil_data["CL_alpha"],
                "CL_max": airfoil_data["CL_max"],
                "Cm_L0": airfoil_data["Cm_L0"],
                "Cm_alpha": airfoil_data["Cm_alpha"],
                "CD_0": airfoil_data["CD0"],
                "CD_L": airfoil_data["CD0_L"],
                "CD_L2": airfoil_data["CD0_L2"]
            }
        else:
            self._properties = {
                "span_position": None,
                "alpha_L0": -0.0369,
                "CL_alpha": 6.2832,
                "CL_max": 1.4,
                "Cm_L0": -0.0527,
                "Cm_alpha": -0.08,
                "CD_0": 0.0055,
                "CD_L": -0.0045,
                "CD_L2": 0.01
            }
    def get_name(self):
        """Get the name of the airfoil.

        Returns
        -------
        string
            The name of the Airfoil.

        """
        return self.name

    def get_span_position(self):
        """Get the spanwise position of the airfoil.

        Returns
        -------
        float
            The spanwise position of the airfoil. Should be a value from 0 to 1

        """
        return self._properties["span_position"]

    def get_lift_slope(self):
        """Get the lift coefficient slope of the Airfoil.

        Returns
        -------
        float
            The lift slope of the Airfoil (per radian).

        """
        return self._properties["CL_alpha"]

    def get_zero_lift_alpha(self):
        """Get the zero-lift angle of attack of the Airfoil.

        Returns
        -------
        float
            The zero-lift angle of attack of the Airfoil (radians).

        """
        return self._properties["alpha_L0"]

    def get_max_lift(self):
        """Get the max lift coefficient of the Airfoil.

        Returns
        -------
        float
            The max lift coefficient of the Airfoil.

        """
        return self._properties["CL_max"]

    def get_moment_slope(self):
        """Get the moment coefficient slope of the Airfoil.

        Returns
        -------
        float
            The moment slope of the Airfoil (per radians).

        """
        return self._properties["Cm_alpha"]

    def get_zero_lift_moment(self):
        """Get the zero-lift moment coefficient of the Airfoil.

        Returns
        -------
        float
            The zero-lift moment coefficient of the Airfoil.

        """
        return self._properties["Cm_L0"]

    def get_drag_coefficients(self):
        """Get the drag coefficients of the Airfoil.

        Returns
        -------
        tuple
            The drag coefficients of the Airfoil.

        """
        return (self._properties["CD_0"],
                self._properties["CD_L"],
                self._properties["CD_L2"])


class ControlSurface:
    """Defines the dimensions and properties of wing control surface.

    Parameters
    ----------

    Returns
    -------
    ControlSurface
        Returns the newly created ControlSurface object.

    Examples
    --------
    """
    def __init__(self, dimensions=None):
        self._span_root = 0.
        self._span_tip = 1.
        self._chord_root = 0.
        self._chord_tip = 0.
        self._is_sealed = True
        self._mix = {}

        if dimensions:
            self._span_root = dimensions["percent_span"][0]
            self._span_tip = dimensions["percent_span"][1]
            try:
                self._chord_root = dimensions["percent_chord"][0]
                self._chord_tip = dimensions["percent_chord"][1]
            except TypeError:
                self._chord_root = dimensions["percent_chord"]
                self._chord_tip = dimensions["percent_chord"]
            self._is_sealed = dimensions.get("sealed", 1)

            mix_dict = dimensions["mix"]
            if "aileron" in mix_dict:
                self._mix["aileron"] = mix_dict["aileron"]
            if "elevator" in mix_dict:
                self._mix["elevator"] = mix_dict["elevator"]
            if "rudder" in mix_dict:
                self._mix["rudder"] = mix_dict["rudder"]
            if "flap" in mix_dict:
                self._mix["flap"] = mix_dict["flap"]

    def get_control_span(self):
        """Get the spanwise location of control surface,

        Returns
        -------
        tuple
            The spanwise location of the root and tip of the control
            surface, normalized by semispan.
        """
        return self._span_root, self._span_tip

    def get_control_chord(self):
        """Get the control surface width.

        Returns
        -------
        tuple
            The control surface width at the root and tip of the control
            surface, normalized by chord.
        """
        return self._chord_root, self._chord_tip

    def get_control_mix(self, side):
        """Get the control surface mixing parameters.

        Parameters
        ----------
        side
            The side of the airplane that the WingSegment is on. This
            allows the method to return the proper mixing for control
            surfaces that are deflected in an asymmetric manner.

        Returns
        -------
        control_mix : Dict
        """
        segment_mix = {}
        if "aileron" in self._mix:
            if side == "right":
                segment_mix["aileron"] = self._mix["aileron"]
            else:
                segment_mix["aileron"] = -self._mix["aileron"]
        if "elevator" in self._mix:
            segment_mix["elevator"] = self._mix["elevator"]
        if "rudder" in self._mix:
            if side == "right":
                segment_mix["rudder"] = self._mix["rudder"]
            else:
                segment_mix["rudder"] = -self._mix["rudder"]
        if "flap" in self._mix:
            segment_mix["flap"] = self._mix["flap"]

        return segment_mix

    def get_is_control_sealed(self):
        """Check if control surface is sealed.

        Returns
        -------
        boolean
            True if control surface has been specified as sealed.

        """
        return self._is_sealed
