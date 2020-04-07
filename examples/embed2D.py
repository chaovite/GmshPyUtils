from GmshPyUtils import gmsh_utils
import numpy as np

geofile='EdgePtEmbed2D.geo'

L = 10.0          
H = 10.0       
size_b = 1.0   
size_l = 0.1   
size_p = 0.05

#define the x,y coordinate the line end points
p1xy = [6, 4]
p2xy = [3, 7]

#define the x,y coordinate the embedded point (injector)
p3xy = [6, 6]


objs = []
rect = gmsh_utils.Rectangle(0, 0, L, H, size_b, z0 = 0)
rect.create_rectangle_geometry(createSurface=True)
objs.append(rect)

surface = rect.Surfaces[0]

p1   = gmsh_utils.Point(p1xy[0], p1xy[1], 0, size_l) 
p2   = gmsh_utils.Point(p2xy[0], p2xy[1], 0, size_l)
p3   = gmsh_utils.Point(p3xy[0], p3xy[1], 0, size_p)

objs +=[p1, p2, p3]

l1   = gmsh_utils.Line(p1,p2)
objs.append(l1)

embed1 = gmsh_utils.Embed([l1], 'Line', surface, 'Surface')
embed2 = gmsh_utils.Embed([p3], 'Point', surface, 'Surface')
objs  +=[embed1, embed2]

with open(geofile,'w') as f:
        for obj in objs:
                    f.write(obj.write_txt())
