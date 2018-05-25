#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pygmsh
import numpy as np
from itertools import chain

from subprocess import Popen

# Iterate an iterable pairwise so that adjacent are paired up, and
# last and first element are paired up.
def iter_loop(it):
    return chain(zip(it[:-1], it[1:]), [(it[-1], it[0])])    

def iterp(it):
    return zip(it[:-1], it[1:])

def show_in_gmsh(file):
    try:
        Popen(["gmsh", file])
    except Exception as e:
        print("Gmsh failed to open:")
        print(e)


#%% Geometry of a coil
c_geom = m_geom = pygmsh.built_in.Geometry()

a_coar = 0.08
c_coar = 0.06

hcdh = 0.6 * 0.5
cdw = 1.0

coil_radius = 0.4
wind_r = 0.1
#wind_w = 0.2

corners = np.array([(0,-1,0),(1,-1,0),(1,1,0),(0,1,0)])

# boundary corners
bcn = np.array([cdw, hcdh, 0])*corners
# winding center and nodes for the arcs
wns = np.array([coil_radius, 0, 0])
wns = wns + np.array([
        [0, 0, 0],
        [wind_r, 0, 0],
        [0, wind_r, 0],
        [-wind_r, 0, 0],
        [0, -wind_r, 0]])
    
    
# points for boundary corners
bps = [c_geom.add_point(p, a_coar) for p in bcn]
coillps = [c_geom.add_point(p, c_coar) for p in wns]

# Adding lines for boundaries
oblines = [m_geom.add_line(p0, p1) for p0, p1 in iter_loop(bps)]
# Symmetry axel
#symmline = m_geom.add_line(bps[-1], bps[0])

cp = coillps[0]
windarcs = [m_geom.add_circle_arc(p0, cp, p1) 
            for p0, p1 in iter_loop(coillps[1:])]

# Adding physical line bounding the domain
obline = m_geom.add_physical_line(oblines[:-1], label=1)
symmline = m_geom.add_physical_line(oblines[-1], label=2)


obll = m_geom.add_line_loop(oblines)
windll = m_geom.add_line_loop(windarcs)


# Adding physical domain with holes
dompls = m_geom.add_plane_surface(obll, holes=[windll])
m_geom.add_physical_surface(dompls, label=3)

windpl = m_geom.add_plane_surface(windll)
m_geom.add_physical_surface(windpl, label=4)

with open('coil.geo', 'w') as file:
    file.write(c_geom.get_code())
    
show_in_gmsh('coil.geo')
#%% Geometry of the medium domain
m_geom = pygmsh.built_in.Geometry()

far_coar = 0.2
near_coar = 0.06

hcoildist = 1.0 * 0.5

medw = 2
hmedh = 1.5

# Simple and straightforward list of node coordinates and
# the coarsity which will be associated
ncords = [
        ([0,    -hmedh,            0], far_coar),  #corners
        ([medw, -hmedh,            0], far_coar),  
        ([medw,  hmedh,            0], far_coar),  
        ([0,     hmedh,            0], far_coar),  
        ([0,     hcoildist + hcdh, 0], near_coar), #upper coil
        ([cdw,   hcoildist + hcdh, 0], near_coar),
        ([cdw,   hcoildist - hcdh, 0], near_coar),
        ([0,     hcoildist - hcdh, 0], near_coar),
        ([0,    -hcoildist + hcdh, 0], near_coar), # lower coil
        ([cdw,  -hcoildist + hcdh, 0], near_coar),
        ([cdw,  -hcoildist - hcdh, 0], near_coar),
        ([0,    -hcoildist - hcdh, 0], near_coar)
        ]


# Nodes in external boundary and coupling boundary corners
points = [m_geom.add_point(p, c) for p, c in ncords]
    
#bns = np.array([2, 2, 0])*corners
#ucdomns = np.array([cdw, hcdh, 0])*corners+np.array([0, hcoildist,0])
#lcdomns = np.array([cdw, hcdh, 0])*corners-np.array([0, hcoildist,0])

# Adding points for boundary corners
bps = points[0:4]
uc1ibps = points[4:8]
lc1ibps = points[8:12]

# Adding lines for infty boundary
inflines = [m_geom.add_line(p0, p1) 
            for p0, p1 in iterp(bps)]

# upper slot boundary
uiblines = [m_geom.add_line(p0, p1) 
            for p0, p1 in iterp(uc1ibps)] 

# lower slot boundary
liblines = [m_geom.add_line(p0, p1) 
            for p0, p1 in iterp(lc1ibps)]

symmlines = [m_geom.add_line(points[3], points[4]),
             m_geom.add_line(points[7], points[8]),
             m_geom.add_line(points[11], points[0])]

# Adding physical lines bounding the domain and coupling domains
infline = m_geom.add_physical_line(inflines, label=1)
symmline = m_geom.add_physical_line(symmlines, label=2)
uibline = m_geom.add_physical_line(uiblines, label=3)
libline = m_geom.add_physical_line(liblines, label=4)

# Line loop for the medium air mesh
obll = m_geom.add_line_loop(inflines + 
                            symmlines[0:1] + 
                            uiblines + 
                            symmlines[1:2] +
                            liblines +
                            symmlines[2:3])

# Adding physical domain with holes
dompls = m_geom.add_plane_surface(obll)
m_geom.add_physical_surface(dompls, label=5)

#print(geom.get_code())

with open('medium.geo', 'w') as file:
    file.write(m_geom.get_code())

show_in_gmsh('medium.geo')

#%% Do the meshing

handles = [Popen(["gmsh", "-2", *mshf]) 
           for mshf in [("coil.geo", "-o" , "coil.msh"), 
                        ("medium.geo", "-o", "medium.msh")]]

for h in handles:
    print(h.wait())

