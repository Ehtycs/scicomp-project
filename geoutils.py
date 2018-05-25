
#import pygmsh
import numpy as np
from itertools import islice


def create_spiral_winding_2d(geom, fst_point, lst_point, nturns, 
                          cond_radius, ccoar, label_one , label_two):
    """ Create a spiral coil 2d slice.
    Inputs: 
        geom - pygmsh geometry where entities will get added. Will be mutated!
        fst_point - first coil conductor center point
        lst_point - last coil conductor center point
        nturns - number of turns
        cond_radius - conductor radius
        ccoar - mesh coarsity near generated points
        kwargs are forwarded to pygmsh add_physical_surfaces (e.g label)
    
    Outputs: 
        list of generated line loops for physical surface generation (holes)
        
    """
        
    xs, ys, zs = [np.linspace(fc, lc, nturns) 
                  for fc,lc in zip(fst_point,lst_point)]
    
    offsets = [(0,0), (1,0), (0,1), (-1,0), (0,-1)]
    allps = []
    lines = []
    surfs = []
    line_loops = []
    
    # Create points for arcs which form the conductor circles. Each circle
    # consists of four arcs.
    allps = [geom.add_point((x+ox*cond_radius, y+oy*cond_radius, z), ccoar)
             for x,y,z in zip(xs, ys, zs)
             for ox,oy in offsets]
    
    for i, _ in islice(enumerate(allps), 0, None, 5):
        # Iterate every fifth, add four arcs to create a circle
        cp = allps[i]
        lines.append(geom.add_circle_arc(allps[i+1], cp, allps[i+2]))
        lines.append(geom.add_circle_arc(allps[i+2], cp, allps[i+3]))
        lines.append(geom.add_circle_arc(allps[i+3], cp, allps[i+4]))
        lines.append(geom.add_circle_arc(allps[i+4], cp, allps[i+1]))
    
    # Create line loops and surfaces from circle arcs.
    for i, _ in islice(enumerate(lines), 0, None, 4):
        ll = geom.add_line_loop(lines[i:i+4])
        line_loops.append(ll)
        surfs.append(geom.add_plane_surface(ll))
    
    halfway = int(len(surfs)/2)
    geom.add_physical_surface(surfs[:halfway], label = label_one)
    geom.add_physical_surface(surfs[halfway:], label = label_two)
    
    return line_loops