#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Generates a 2D finite element domain to be meshed with GMSH (a geo file).
    This script generates a cross section of a spiral coil.
"""

import pygmsh
import geoutils as guti
import numpy as np
from itertools import chain

# Iterate an iterable pairwise so that adjacent are paired up, and
# last and first element are paired up.
def iter_loop(it):
    return chain(zip(it[:-1], it[1:]), [(it[-1], it[0])])    


def show_in_gmsh(file):
    try:
        Popen(["gmsh", file])
    except Exception as e:
        print("Gmsh failed to open:")
        print(e)


#%% Geometry of a coil
c_geom = pygmsh.built_in.Geometry()

a_coar = 0.05
c_coar = 0.003

cdomain_height = 0.4
cdomain_width = 1.0

coil_radius = 0.75
conductor_radius = 0.01

corners = np.array([(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)])

# boundary corners
bcn = 0.5*np.array([cdomain_width, cdomain_height, 0])*corners

# points for boundary corners
bps = [c_geom.add_point(p, a_coar) for p in bcn]

# Adding lines for boundaries
oblines = [c_geom.add_line(p0, p1) for p0, p1 in iter_loop(bps)]

# Adding physical line bounding the domain
obline = c_geom.add_physical_line(oblines, label=1)

obll = c_geom.add_line_loop(oblines)

create_coil = guti.create_spiral_winding_2d

wind_lineloops = create_coil(c_geom, (-0.4, 0, 0), (0.4, 0, 0), 18, 
                             conductor_radius, c_coar, 3, 4)

# Adding physical domain with holes
dompls = c_geom.add_plane_surface(obll, holes=wind_lineloops)
c_geom.add_physical_surface(dompls, label=2)

with open('spiral_coil.geo', 'w') as file:
    file.write(c_geom.get_code())

show_in_gmsh("spiral_coil.geo")

#%% Geometry of the medium domain
m_geom = pygmsh.built_in.Geometry()

far_coar = 0.2
near_coar = 0.06

coil_dist = 0.6

# Nodes in external boundary and coupling boundary corners
bns = 0.5*np.array([2, 2, 0])*corners
ucdomns = 0.5*np.array([cdomain_width, cdomain_height, 0])*corners+np.array([0, coil_dist*0.5,0])
lcdomns = 0.5*np.array([cdomain_width, cdomain_height, 0])*corners-np.array([0, coil_dist*0.5,0])

# Adding points for boundary corners
bps = [m_geom.add_point(p, far_coar) for p in bns]
uc1ibps = [m_geom.add_point(p, near_coar) for p in ucdomns]
lc1ibps = [m_geom.add_point(p, near_coar) for p in lcdomns]

# Adding lines for boundaries
oblines = [m_geom.add_line(p0, p1) for p0, p1 in iter_loop(bps)]
uiblines = [m_geom.add_line(p0, p1) for p0, p1 in iter_loop(uc1ibps)]
liblines = [m_geom.add_line(p0, p1) for p0, p1 in iter_loop(lc1ibps)]

# Adding physical lines bounding the domain and coupling domains
obline = m_geom.add_physical_line(oblines, label=1)
uibline = m_geom.add_physical_line(uiblines, label=2)
libline = m_geom.add_physical_line(liblines, label=3)

obll = m_geom.add_line_loop(oblines)
uibll = m_geom.add_line_loop(uiblines)
libll = m_geom.add_line_loop(liblines)

# Adding physical domain with holes
dompls = m_geom.add_plane_surface(obll, holes=[uibll,libll])
m_geom.add_physical_surface(dompls, label=4)

#print(geom.get_code())

with open('spiral_medium.geo', 'w') as file:
    file.write(m_geom.get_code())
    
from subprocess import Popen


show_in_gmsh("spiral_medium.geo")


