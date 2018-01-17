import numpy as np
import matplotlib.pyplot as plt
from numpy import cos as c
from numpy import sin as s
from stl import mesh


def create_from_grid(llgrid = None, prop_models = None,filename = None):

    num_faces = 0
    n = 50

    if llgrid:
        wing_segments = llgrid._plane.get_wingsegments()
        total_sections = llgrid._num_sections
        num_faces += (total_sections-len(wing_segments))*(n-1)*2

    if prop_models:
        for prop in prop_models:
            num_faces += prop._k*(prop._nodes-1)*(n-1)*2

    plane_stl = mesh.Mesh(np.zeros(num_faces, dtype = mesh.Mesh.dtype))

    index=0

    if llgrid:
        control_points = llgrid.get_control_point_pos()
        chord = llgrid.get_chord_length()
        twist = llgrid.get_twist()
        dihedral = llgrid.get_dihedral()
        percent = llgrid.get_cp_spacing()
        slices = llgrid._segment_slices

        segment_index = 0
        for counter, seg in enumerate(wing_segments):
            num_sections = seg.get_num_sections()
            if seg.get_side() == "left":
                right_airfoil, left_airfoil = seg.get_airfoils()
            else:
                left_airfoil, right_airfoil = seg.get_airfoils()

            left_name = left_airfoil.get_name()
            right_name = right_airfoil.get_name()
        
            if left_name[:4] == 'NACA' or left_name[:4] == 'naca':
                if len(left_name) == 8:
                    left_NACA = left_name[4:]
            else: left_NACA = '2412'
                    
            if right_name[:4] == 'NACA' or right_name[:4] == 'naca':
                if len(right_name) == 8:
                    right_NACA = right_name[4:]
            else: right_NACA = '2412'

            left_X,left_Z = plotNACA(left_NACA,n)
            right_X,right_Z = plotNACA(right_NACA,n)
            X_coords = linear_inter(left_X, right_X, percent[slices[counter]])
            Z_coords = linear_inter(left_Z, right_Z, percent[slices[counter]])

            X_coords *=chord[slices[counter],None]
            Z_coords *=chord[slices[counter],None]

            Coords = np.zeros((num_sections,n,3))
            Coords[:,:,0] = X_coords
            Coords[:,:,2] = Z_coords

            if seg.get_side() == "left":
                e1 = Euler2Quat((np.radians(dihedral[slices[counter]]),0., 0.))
            else:
                e1 = Euler2Quat((np.radians(-dihedral[slices[counter]]),0., 0.))
            e2 = Euler2Quat((0.,np.radians(twist[slices[counter]]),0.))
            e = np.transpose(QuatMult(e1,e2))
            
            for i in range(num_sections):
                for j in range(n):
                    Coords[i][j] = Body2Fixed(Coords[i][j], e[i])

            Coords = Coords+control_points[slices[counter],None,:]

            for i in range(num_sections-1):
                section_index = 2*(n-1)*i
                for j in range(n-1):
                    slice_index = 2*j
                    plane_stl.vectors[segment_index+section_index+slice_index][0] = Coords[i][j]
                    plane_stl.vectors[segment_index+section_index+slice_index][1] = Coords[i][j+1]
                    plane_stl.vectors[segment_index+section_index+slice_index][2] = Coords[i+1][j]
                    plane_stl.vectors[segment_index+section_index+slice_index+1][0] = Coords[i][j+1]
                    plane_stl.vectors[segment_index+section_index+slice_index+1][1] = Coords[i+1][j+1]
                    plane_stl.vectors[segment_index+section_index+slice_index+1][2] = Coords[i+1][j]

            segment_index+=(num_sections-1)*(n-1)*2
        index+=segment_index

    if prop_models:
        for prop in prop_models:
            prop_geom = prop._prop
            num_nodes = prop.get_num_nodes()
            num_blades = prop.get_num_blades()
            position = prop.get_position()
            rot_dir = prop.get_rotation_direction()

            root_airfoil,tip_airfoil = prop_geom.get_airfoils()

            root_name = root_airfoil.get_name()
            tip_name = tip_airfoil.get_name()
        
            if root_name[:4] == 'NACA' or root_name[:4] == 'naca':
                if len(root_name) == 8:
                    root_NACA = root_name[4:]
            else: root_NACA = '2412'
                    
            if tip_name[:4] == 'NACA' or tip_name[:4] == 'naca':
                if len(tip_name) == 8:
                    tip_NACA = tip_name[4:]
            else: tip_NACA = '2412'


            root_X,root_Z = plotNACA(root_NACA,n)
            tip_X,tip_Z = plotNACA(tip_NACA,n)
            X_coords = linear_inter(root_X,tip_X,prop._zeta)
            Z_coords = linear_inter(root_Z,tip_Z,prop._zeta)
            
            chord = prop.get_chord()
            X_coords *= chord[:,None]
            Z_coords *= chord[:,None]

            Coords = np.zeros((num_nodes,n,3))
            Coords[:,:,0] = -Z_coords
            Coords[:,:,1] =  X_coords*rot_dir
            twist = prop.get_aero_pitch_angle()-prop.get_zero_lift_alpha()
            e = Euler2Quat((0.,0.,-rot_dir*twist))
            e = np.transpose(e)


            for i in range(num_nodes):
                for j in range(n):
                    Coords[i][j] = Body2Fixed(Coords[i][j], e[i])
            Coords[:,:,2] = -prop.get_radial_position()[:,None]

            blade_rotation = 2*np.pi/num_blades
            e_blade = Euler2Quat((blade_rotation,0.,0.))

            theta,gamma = prop.get_orientation()
            e_full_prop = Euler2Quat((0.,theta,gamma))
            rotated_coords = np.zeros_like(Coords)

            for k in range(prop._k):
                for i in range(num_nodes):
                    for j in range(n):
                        rotated_coords[i][j] = Body2Fixed(Coords[i][j], e_full_prop)

                rotated_coords+=position
                blade_index = (num_nodes-1)*(n-1)*2*k
                for i in range(num_nodes-1):
                    node_index = 2*(n-1)*i
                    for j in range(n-1):
                        slice_index = 2*j
                        if rot_dir == 1:
                            plane_stl.vectors[index+blade_index+node_index+slice_index][0] = rotated_coords[i][j]
                            plane_stl.vectors[index+blade_index+node_index+slice_index][1] = rotated_coords[i][j+1]
                            plane_stl.vectors[index+blade_index+node_index+slice_index][2] = rotated_coords[i+1][j]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][0] = rotated_coords[i][j+1]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][1] = rotated_coords[i+1][j+1]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][2] = rotated_coords[i+1][j]
                        else:
                            plane_stl.vectors[index+blade_index+node_index+slice_index][0] = rotated_coords[i][j]
                            plane_stl.vectors[index+blade_index+node_index+slice_index][1] = rotated_coords[i+1][j]
                            plane_stl.vectors[index+blade_index+node_index+slice_index][2] = rotated_coords[i][j+1]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][0] = rotated_coords[i][j+1]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][1] = rotated_coords[i+1][j]
                            plane_stl.vectors[index+blade_index+node_index+slice_index+1][2] = rotated_coords[i+1][j+1]

                for i in range(num_nodes):
                    for j in range(n):
                        Coords[i][j] = Body2Fixed(Coords[i][j], e_blade)

            index+=num_blades*(num_nodes-1)*(n-1)*2

    plane_stl.save(filename)


def plotNACA(airfoil,n):
    max_camb = float(airfoil[0])/100
    max_camb_loc = float(airfoil[1])/10
    max_thickness = float(airfoil[2:4])/100
    m = max_camb
    p = max_camb_loc
    t = max_thickness

    #cosine spacing
    theta = np.linspace(-np.pi,np.pi,n)
    Xc = 0.5*(1-np.cos(theta))
    Yt = 5*t*(0.2969*np.sqrt(Xc)-0.1260*Xc-0.3516*Xc**2+0.2843*Xc**3-0.1015*Xc**4)
    if m == 0 or p == 0:
        Yc = np.zeros_like(Yt)
        Yc_x = np.copy(Yc)
    else:
        Yc = np.where(Xc<p,(m/(p*p))*(2*p*Xc-Xc*Xc),(m/((1-p)**2))*((1-2*p)+2*p*Xc-Xc*Xc))
        Yc_x = np.where(Xc<p,(2*m/(p*p))*(p-Xc),(m/(1-p)**2)*(p-Xc))

    Yul = Yc+(Yt*np.cos(np.arctan(Yc_x))*np.sign(theta))
    Xul = Xc-(Yt*np.sin(np.arctan(Yc_x))*np.sign(theta))

    X = -Xul+0.25
    Z = -Yul
    return X,Z

def linear_inter(left, right, spacing):
    dif = right-left
    interpolated = left[None,:]+dif[None,:]*spacing[:,None]
    return interpolated

def Euler2Quat(euler):
    p = euler[0]/2.
    t = euler[1]/2.
    g = euler[2]/2.

    cp = c(p)
    ct = c(t)
    cg = c(g)
    sp = s(p)
    st = s(t)
    sg = s(g)
    return [cp*ct*cg+sp*st*sg,
            sp*ct*cg-cp*st*sg,
            cp*st*cg+sp*ct*sg,
            cp*ct*sg-sp*st*cg]

def Body2Fixed(vec,e):
	x,y,z = vec
	e0,ex,ey,ez = e

	To =  x*ex+y*ey+z*ez
	Tx =  x*e0-y*ez+z*ey
	Ty =  x*ez+y*e0-z*ex
	Tz = -x*ey+y*ex+z*e0

	a = e0*Tx+ex*To+ey*Tz-ez*Ty
	b = e0*Ty-ex*Tz+ey*To+ez*Tx
	c = e0*Tz+ex*Ty-ey*Tx+ez*To


	return [a,b,c]

def QuatMult(A, B):
	w1,x1,y1,z1 = A
	w2,x2,y2,z2 = B
	return np.array([w1*w2-x1*x2-y1*y2-z1*z2,
			w1*x2+x1*w2+y1*z2-z1*y2,
			w1*y2-x1*z2+y1*w2+z1*x2,
			w1*z2+x1*y2-y1*x2+z1*w2])

