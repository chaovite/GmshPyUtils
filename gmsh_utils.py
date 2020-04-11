"""
gmsh_utils.py is python module that contains utility functions 
that helps to generate gmsh geo files.
"""
import numpy as np

class Point:
    # use static variable to track the counts
    PointIndex = 0
    def __init__(self, x, y, z, size):
        Point.PointIndex += 1
        self.index= Point.PointIndex
        self.x = x
        self.y = y
        self.z = z
        self.size = size
        
    def reset_index(index=0):
        Point.PointIndex  = index
        
    def write_txt(self):
        txt = "Point(%d) = {%.5f, %.5f, %.5f, %.5f};\n\n"%(
                self.index, self.x, self.y, self.z, self.size)
        return txt

class Line:
    # use static variable to track the counts
    LineIndex = 0
    def __init__(self, start, end, center=None, major=None):
        """
        index: line index
        start, end: start and end points
        center: None (straight line), point for circular arc
        major: only needed for elliptical arc
        """
        Line.LineIndex  += 1
        self.index  = Line.LineIndex
        self.start  = start
        self.end    = end
        self.center = center
        self.major  = major
    
    #@classmethod
    def reset_index(index=0):
        Line.LineIndex  = index
    
    def invert(self):
        self.index = -self.index
        
    def write_txt(self):
        if not self.center:
            # straight line
            txt = "Line(%d) = {%d, %d};\n\n"%(
                  self.index, self.start.index, 
                  self.end.index)
        else:
            if not self.major:
                # circle arc
                txt = "Circle(%d) = {%d, %d, %d};\n\n"%(
                      self.index, self.start.index, 
                      self.center.index, self.end.index)
            else:
                # ellipse arc
                txt = "Ellipse(%d) = {%d, %d, %d, %d};\n\n"%(
                      self.index, self.start.index, 
                      elf.center.index, self.major.index)
                
        return txt

class LineLoop:
    # use static variable to track the counts
    LineLoopIndex = 0
    def __init__(self, lines, direction=1):
        LineLoop.LineLoopIndex += 1
        self.index = LineLoop.LineLoopIndex
        self.lines = lines
        self.lines_index = []
        self.direction = direction
        self.compute_LineLoop()
        
    def reset_index(index=0):
        LineLoop.LineLoopIndex = index
        
    def invert(self):
        self.index = -self.index
        
    def compute_LineLoop(self):
        """
        reorder the lines so that they for a loop
        
        direction is set by the first line and direction
        """
        
        line_index = []
        prevline   = self.lines[0]
        
        # add the first line into the loop
        line_index.append(prevline.index)
        end_pt_index= prevline.end.index
        
        countsMax = 100
        count  = 0
        while len(line_index)<len(self.lines):
            count += 1
            if count>countsMax:
                print('error in defining line loop!')
                break
            for nextline in self.lines:
                ind = nextline.index
                if (ind in line_index) or (-ind in line_index):
                    # has already been added
                    continue
                
                if nextline.start.index==end_pt_index:
                    line_index.append(ind)
                    end_pt_index = nextline.end.index
                elif nextline.end.index==end_pt_index:
                    line_index.append(-ind)
                    end_pt_index = nextline.start.index
                else:
                    pass
        self.lines_index = [i*self.direction for i in line_index]
    def write_txt(self):
        line_index = self.lines_index
        txt = "Line Loop(%d) = {%s};\n\n"%(
              self.index, str(line_index).strip('[]'))
        return txt

class Surface:
    # use static variable to track the counts
    SurfaceIndex = 0
    def __init__(self, lineLoops, isPlaneSurface=False):
        """
        lineLoops[0] defines the boundary
        lineLoops[1:] defines the holes
        """
        Surface.SurfaceIndex  += 1
        self.index     = Surface.SurfaceIndex
        self.lineLoops = lineLoops
        self.isPlaneSurface = isPlaneSurface

    def invert(self):
        self.index = -self.index
    def reset_index(index=0):
        Surface.SurfaceIndex  = index
    def write_txt(self):
        index_lineloops = [lineloop.index 
                           for lineloop in self.lineLoops]
        if not self.isPlaneSurface:
            txt = "Surface(%d)={%s};\n\n"%(
               self.index, str(index_lineloops).strip('[]'))
        else:
            txt = "Plane Surface(%d)={%s};\n\n"%(
               self.index, str(index_lineloops).strip('[]'))

        return txt

class SurfaceLoop:
    # use static variable to track the counts
    SurfaceLoopIndex = 0
    def __init__(self, Surfaces, normals=None):
        SurfaceLoop.SurfaceLoopIndex += 1
        self.index = SurfaceLoop.SurfaceLoopIndex
        self.Surfaces = Surfaces
        if not normals:
            normals = [1]*len(Surfaces)
        self.normals  = normals
        
    def reset_index(index=0):
        SurfaceLoop.SurfaceLoopIndex = index
    
    def invert():
        self.index = -self.index
        
    def write_txt(self):
        Surface_index = [i.index*j for (i,j) in zip(self.Surfaces, self.normals)]
        txt = "Surface Loop(%d) = {%s};\n\n"%(
              self.index, str(Surface_index).strip('[]'))
        return txt

class Volume:
    # use static variable to track the counts
    VolumeIndex = 0
    def __init__(self, SurfaceLoops):
        """
        SurfaceLoops[0] defines the boundary
        SurfaceLoops[1:] defines the holes
        """
        Volume.VolumeIndex  += 1
        self.index     = Volume.VolumeIndex
        self.SurfaceLoops = SurfaceLoops
    
    def reset_index(index=0):
        Volume.VolumeIndex = index
        
    def write_txt(self):
        index_surfaceloops = [surfaceloop.index 
                              for surfaceloop in self.SurfaceLoops]
        txt = "Volume(%d)={%s};\n\n"%(
               self.index, str(index_surfaceloops).strip('[]'))
        return txt

class PhysicalGroup:
    """
    create a physical group:
    Physical line, surface, or volume
    """
    def __init__(self, objType, objList, tag):
        assert objType in ['Line','Surface','Volume']
        self.tag     = tag        
        self.objType = objType
        self.objList = objList
    
    def write_txt(self):
        inds = [o.index for o in self.objList]
        if type(self.tag) == str:
            txt  = "Physical %s(%s) = {%s}; \n\n"%(self.objType, self.tag, 
                                                   str(inds).strip('[]'))
        else:
            txt  = "Physical %s(%d) = {%s}; \n\n"%(self.objType, self.tag, 
                                                   str(inds).strip('[]'))
        return txt

class Field:
    """
    create a size field
    
    A list of types supported:
    Threshold, Min, MathEval
    
    Min: Take the minimum value of a list of fields.
        options:
            -FieldsList
    
    MathEval: Evaluate a mathematical expression. 
              The expression can contain x, y, z for spatial coordinates, 
              F0, F1, ... for field values, and and mathematical functions.
        options:
            -F: a math expression string.
    
    Box: The value of this field is VIn inside the box, VOut outside the box. The box is given by
            Xmin <= x <= XMax &&
            YMin <= y <= YMax &&
            ZMin <= z <= ZMax
         options:
            -VIn (float): Value inside the box
            -VOut (float): Value outside the box
            -XMax (float): Maximum X coordinate of the box
            -XMin (float): Minimum X coordinate of the box
            -YMax (float): Maximum Y coordinate of the box
            -YMin (float): Minimum Y coordinate of the box
            -ZMax (float): Maximum Z coordinate of the box
            -ZMin (float): Minimum Z coordinate of the box
    
    Threshold: F = LCMin if Field[IField] <= DistMin,
               F = LCMax if Field[IField] >= DistMax,
               F = interpolation between LcMin and LcMax if DistMin < Field[IField] < DistMax
        options:
            -DistMax (float):  Distance from entity after which element size will be LcMax
            -DistMin (float):  Distance from entity up to which element size will be LcMin
            -IField (integer): Index of the field to evaluate
            -LcMax  (float):    Element size outside DistMax
            -LcMin  (float):    Element size inside DistMin
            -Sigmoid (boolean): True to interpolate between LcMin and LcMax using a sigmoid, 
                                false to interpolate linearly
            -StopAtDistMax (boolean): True to not impose element size outside DistMax 
                                    (i.e., F = a very big value if Field[IField] > DistMax)
    """
    
    FieldIndex = 0
    def __init__(self, ftype, options, setbackground=False):
        Field.FieldIndex  += 1
        self.index         = Field.FieldIndex
        self.ftype         = ftype
        self.options       = options
        self.setbackground = setbackground
        
    def reset_index(index=0):
        Field.FieldIndex = index
    
    def write_txt(self):
        txt = "Field[%d] = %s;\n\n"%(self.index, self.ftype)
        
        txt_options = ""
        
        ListFields  = ["FieldsList"]
        
        IntFields   = ["IField","Sigmoid","StopAtDistMax"]
        
        StrFields   = ['F']
        
        FloatFields = ['DistMax','DistMin','LcMax',"LcMin",
                      'VIn','VOut','XMax','YMax','ZMax',
                      'XMin','YMin','ZMin']
        
        txt_option  = ""
        for k in self.options:
            if k in ListFields:
                str_list = str(self.options[k]).strip('[]')
                txt_option    += "Field[%d].%s = {%s};\n\n"%(self.index, k, str_list)
            if k in IntFields:
                txt_option    += "Field[%d].%s = %d;\n\n"%(self.index, k, self.options[k])
            if k in StrFields:
                txt_option    += """Field[%d].%s = "%s";\n\n"""%(self.index, k, self.options[k])
            if k in FloatFields:
                txt_option    += "Field[%d].%s = %.5f;\n\n"%(self.index, k, self.options[k])
        
        txt += txt_option
        if self.setbackground:
            txt += "Background Field = %d;\n\n"%(self.index)
        
        return txt

class Embed:
    """
    A class that embed points, lines, surfaces into higher order objects
    
    Support types: 
        points/lines/surfaces in volume
        points/lines in in surface
    """
    def __init__(self, objsLowDim, typeLowDim, objHighDim, typeHighDim):
        assert typeLowDim in ['Point','Line','Surface']
        assert typeHighDim in ['Surface','Volume']
        if typeLowDim=='Surface':
            assert typeHighDim=='Volume'

        self.objsLowDim = objsLowDim
        self.typeLowDim = typeLowDim
        self.objHighDim = objHighDim
        self.typeHighDim= typeHighDim
    
    def write_txt(self):
        objsLowIndexList = [str(i.index) for i in self.objsLowDim]
        objsLowIndexListStr = ','.join(objsLowIndexList)
        txt = "%s{%s} In %s{%d};\n\n"%(self.typeLowDim, objsLowIndexListStr,
                                      self.typeHighDim, self.objHighDim.index)
        return txt
    
class Sphere:
    def __init__(self, xc, yc, zc, r, size):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.r  = r
        self.size = size
        self.Points        = []
        self.Lines         = []
        self.LineLoops     = []
        self.Surfaces      = []
        self.SurfaceLoops  = []
        self.Volumes       = []
        
    def create_sphere_geometry(self):
        xc = self.xc
        yc = self.yc
        zc = self.zc
        r  = self.r
        size = self.size
        
        # create all the points
        p0 = Point(xc, yc, zc, size)
        p1 = Point(xc+r, yc, zc, size)
        p2 = Point(xc, yc+r, zc, size)
        p3 = Point(xc-r, yc, zc, size)
        p4 = Point(xc, yc-r, zc, size)
        p5 = Point(xc, yc, zc+r, size)
        p6 = Point(xc, yc, zc-r, size)
        self.Points=[p0,p1,p2,p3,p4,p5,p6]
        
        # create lines: circular arcs
        arc1 = Line(p1, p2, p0)
        arc2 = Line(p2, p3, p0)
        arc3 = Line(p3, p4, p0)
        arc4 = Line(p4, p1, p0)
        
        arct1 = Line(p5, p1, p0)
        arct2 = Line(p5, p2, p0)
        arct3 = Line(p5, p3, p0)
        arct4 = Line(p5, p4, p0)
        
        arcb1 = Line(p6, p1, p0)
        arcb2 = Line(p6, p2, p0)
        arcb3 = Line(p6, p3, p0)
        arcb4 = Line(p6, p4, p0)
        
        self.Lines =[arc1,arc2,arc3,arc4,
                     arct1,arct2,arct3,arct4,
                     arcb1,arcb2,arcb3,arcb4]
        
        # create LineLoops top hemisphere
        
        # create the lineloop such that 
        # the normal (right hand rule) points outside of the sphere
        #
        
        LP1 = LineLoop([arct1, arc1, arct2]) 
        LP2 = LineLoop([arct2, arc2, arct3]) 
        LP3 = LineLoop([arct3, arc3, arct4]) 
        LP4 = LineLoop([arct4, arc4, arct1]) 
        
        # create LineLoops bottom hemisphere
        LP5 = LineLoop([arcb2, arc1, arcb1]) 
        LP6 = LineLoop([arcb3, arc2, arcb2]) 
        LP7 = LineLoop([arcb4, arc3, arcb3]) 
        LP8 = LineLoop([arcb1, arc4, arcb4])
        
        self.LineLoops = [LP1, LP2, LP3, LP4,
                          LP5, LP6, LP7, LP8]
        
        # create 8 Surfaces
        S1 = Surface([LP1])
        S2 = Surface([LP2])
        S3 = Surface([LP3])
        S4 = Surface([LP4])
        S5 = Surface([LP5])
        S6 = Surface([LP6])
        S7 = Surface([LP7])
        S8 = Surface([LP8])
        
        self.Surfaces = [S1, S2, S3, S4, 
                         S5, S6, S7, S8]
        # create SurfaceLoops
        SL = SurfaceLoop(self.Surfaces)
        self.SurfaceLoops = [SL]
        
        # create volumes
        V  = Volume(self.SurfaceLoops)
        self.Volumes = [V]
        
    def write_txt(self):
        
        txt = "// Creating a sphere\n"
        
        # concatenate all the objects
        # Points, Lines, LineLoops, Surfaces, SurfaceLoops, Volumes
        #
        
        objs = self.Points + self.Lines + self.LineLoops + \
                  self.Surfaces + self.SurfaceLoops + self.Volumes
        
        # write all the objects
        for obj in objs:
            txt+=obj.write_txt()
        
        txt += "// Done creating a sphere\n"
            
        return txt

class EllipseCylinderVertical:
    """
    build a cylinder with elliptical cross-section, the top and base layers can 
    have slight offsets in center and orientation
    """
    def __init__(self, xyzBase, rxBase, ryBase, angleBase, sizeBase, 
                      xyzTop, rxTop, ryTop, angleTop, sizeTop, 
                      createSurfaceBase=False, createSurfaceTop=False,
                      createVolume=False):
        self.xyzBase=xyzBase
        self.rxBase = rxBase
        self.ryBase = ryBase
        self.angleBase = angleBase
        self.sizeBase  = sizeBase
        self.createSurfaceBase=createSurfaceBase

        self.xyzTop=xyzTop
        self.rxTop = rxTop
        self.ryTop = ryTop
        self.angleTop = angleTop
        self.sizeTop  = sizeTop
        self.createSurfaceTop=createSurfaceTop
        
        self.Points         = []
        self.Lines          = []
        self.LineLoops      = []
        self.Surfaces       = []
        self.LineLoopBase   = None
        self.LineLoopTop    = None

    def createGeometry(self):

        # create base geometry
        xc,yc,zc = self.xyzBase
        rx = self.rxBase
        ry = self.ryBase
        
        size  = self.sizeBase
        theta = self.angleBase*np.pi/180.
        
        # create all the points
        p0Base   = Point(xc,   yc, zc, size)
        p1x, p1y  =  rx*np.cos(theta) + xc,  rx*np.sin(theta) + yc
        p2x, p2y  = -ry*np.sin(theta) + xc,  ry*np.cos(theta) + yc
        p3x, p3y  = -rx*np.cos(theta) + xc, -rx*np.sin(theta) + yc
        p4x, p4y  = ry*np.sin(theta)  + xc, -ry*np.cos(theta) + yc
        
        p1Base = Point(p1x, p1y, zc, size)
        p2Base = Point(p2x, p2y, zc, size)
        p3Base = Point(p3x, p3y, zc, size)
        p4Base = Point(p4x, p4y, zc, size)
        
        self.Points = [p0Base, p1Base, p2Base, p3Base, p4Base]
        
        # create arc
        line1Base = Line(p1Base, p2Base, center=p0Base, major=p1Base)
        line2Base = Line(p2Base, p3Base, center=p0Base, major=p1Base)
        line3Base = Line(p3Base, p4Base, center=p0Base, major=p1Base)
        line4Base = Line(p4Base, p1Base, center=p0Base, major=p1Base)
        
        self.Lines = [line1Base, 
                      line2Base, 
                      line3Base, 
                      line4Base]
        
        # create line loop
        LP1   =  LineLoop([line1Base, line2Base, line3Base, line4Base], direction = -1)
        
        self.LineLoops = [LP1]
        
        # create Top geometry
        xc,yc,zc = self.xyzTop
        rx = self.rxTop
        ry = self.ryTop
        
        size  = self.sizeTop
        theta = self.angleTop*np.pi/180.
        
        # create all the points
        p0Top   = Point(xc,   yc, zc, size)
        p1x, p1y  =  rx*np.cos(theta) + xc,  rx*np.sin(theta) + yc
        p2x, p2y  = -ry*np.sin(theta) + xc,  ry*np.cos(theta) + yc
        p3x, p3y  = -rx*np.cos(theta) + xc, -rx*np.sin(theta) + yc
        p4x, p4y  = ry*np.sin(theta)  + xc, -ry*np.cos(theta) + yc
        
        p1Top = Point(p1x, p1y, zc, size)
        p2Top = Point(p2x, p2y, zc, size)
        p3Top = Point(p3x, p3y, zc, size)
        p4Top = Point(p4x, p4y, zc, size)
        
        self.Points += [p0Top, p1Top, p2Top, p3Top, p4Top]
        
        # create arc
        line1Top = Line(p1Top, p2Top, center=p0Top, major=p1Top)
        line2Top = Line(p2Top, p3Top, center=p0Top, major=p1Top)
        line3Top = Line(p3Top, p4Top, center=p0Top, major=p1Top)
        line4Top = Line(p4Top, p1Top, center=p0Top, major=p1Top)
        
        self.Lines += [line1Top, 
                      line2Top, 
                      line3Top, 
                      line4Top]
        
        # create line loop
        LP2   =  LineLoop([line1Top, line2Top, line3Top, line4Top], direction = -1)
        
        self.LineLoops += [LP2]

        # create lines connecting Base and Top


class Ellipse:
    """
    A class to build a Ellipse
    """
    def __int__(self, xc, yc, rx, ry, size, theta=0, zc=0):
        self.xc    = xc
        self.yc    = yc
        self.zc    = zc
        self.rx    = rx
        self.ry    = ry 
        self.size  = size
        self.theta = theta # angle from x direction, counter-clock wise
        
        self.Points        = []
        self.Lines         = []
        self.LineLoops     = []
        self.Surfaces      = []
        
    def create_ellipse_geometry(self):
        xc = self.xc
        yc = self.yc
        zc = self.zc
        rx = self.rx
        ry = self.ry
        
        size  = self.size
        theta = self.theta*np.pi/180.
        
        # create all the points
        p0   = Point(xc,   yc, zc, size)
        p1x, p1y  =  rx*np.cos(theta) + xc,  rx*np.sin(theta) + yc
        p2x, p2y  = -ry*np.sin(theta) + xc,  ry*np.cos(theta) + yc
        p3x, p3y  = -rx*np.cos(theta) + xc, -rx*np.sin(theta) + yc
        p4x, p4y  = ry*np.sin(theta) + xc,  -ry*np.cos(theta) + yc
        
        p1 = Point(p1x, p1y, zc, size)
        p2 = Point(p2x, p2y, zc, size)
        p3 = Point(p3x, p3y, zc, size)
        p4 = Point(p4x, p4y, zc, size)
        
        self.Points = [p0, p1, p2, p3, p4]
        
        # create arc
        line1 = Line(p1, p2, center=p0, major=p1)
        line2 = Line(p2, p3, center=p0, major=p1)
        line3 = Line(p3, p4, center=p0, major=p1)
        line4 = Line(p4, p1, center=p0, major=p1)
        
        self.Lines = [line1, line2, line3, line4]
        
        # create line loop
        LP1   =  LineLoop([line1, line2, line3, line4], direction = 1)
        
        self.LineLoops = [LP1]
        
        # create surface
        SF1   = Surface([LP1])
        self.Surfaces = [SF1]
        
    def write_txt(self):
        
        txt = "// Creating an ellipse \n"
        
        # concatenate all the objects
        # Points, Lines, LineLoops, Surfaces
        #
        
        objs = self.Points + self.Lines + \
               self.LineLoops + self.Surfaces
               
        
        # write all the objects
        for obj in objs:
            txt+=obj.write_txt()
        
        txt += "// Done creating an ellipse \n"
            
        return txt
    
class Circle:
    """
    A class to build a circle
    """
    def __init__(self, xc, yc, r, size, zc=0):
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.r  = r
        self.size = size
        self.Points        = []
        self.Lines         = []
        self.LineLoops     = []
        self.Surfaces      = []
        
    def create_circle_geometry(self):
        xc = self.xc
        yc = self.yc
        zc = self.zc
        r  = self.r
        size = self.size
        
        # create all the points
        p0 = Point(xc,   yc, zc, size)
        p1 = Point(xc+r, yc, zc, size)
        p2 = Point(xc, yc+r, zc, size)
        p3 = Point(xc-r, yc, zc, size)
        p4 = Point(xc, yc-r, zc, size)
        
        self.Points = [p0, p1, p2, p3, p4]
        
        # create arc
        line1 = Line(p1, p2, center=p0, major=None)
        line2 = Line(p2, p3, center=p0, major=None)
        line3 = Line(p3, p4, center=p0, major=None)
        line4 = Line(p4, p1, center=p0, major=None)
        
        self.Lines = [line1, line2, line3, line4]
        
        # create line loop
        LP1   =  LineLoop([line1, line2, line3, line4], direction = 1)
        
        self.LineLoops = [LP1]
        
        # create surface
        SF1   = Surface([LP1])
        self.Surfaces = [SF1]
        
    def write_txt(self):
        
        txt = "// Creating a circle \n"
        
        # concatenate all the objects
        # Points, Lines, LineLoops, Surfaces, SurfaceLoops, Volumes
        #
        
        objs = self.Points + self.Lines + \
               self.LineLoops + self.Surfaces
               
        
        # write all the objects
        for obj in objs:
            txt+=obj.write_txt()
        
        txt += "// Done creating a circle \n"
            
        return txt

class Rectangle:
    """
    A class to build a rectangular surface
    """
    def __init__(self, x0, y0, dx, dy, size, z0 = 0):
        self.x0   = x0
        self.y0   = y0
        self.z0   = z0
        self.dx   = dx
        self.dy   = dy
        self.size = size
        self.Points     =[]
        self.Lines      = []
        self.LineLoops  = []
        self.Surfaces   = []
        
    def create_rectangle_geometry(self, createSurface=True):
        x    = self.x0
        y    = self.y0
        z    = self.z0
        dx   = self.dx
        dy   = self.dy
        size = self.size
        
        # creating points
        p1   = Point(x, y, z, size)
        p2   = Point(x+dx, y, z, size)
        p3   = Point(x+dx, y+dx, z, size)
        p4   = Point(x, y+dx, z, size)
        
        self.Points = [p1, p2, p3, p4]
        
        # creating lines
        line1 = Line(p1, p2)
        line2 = Line(p2, p3)
        line3 = Line(p3, p4)
        line4 = Line(p4, p1)
        
        self.Lines = [line1, line2, line3, line4]
        
        # creating line loops
        LP1   = LineLoop([line1, line2, line3, line4], direction = 1)
        self.LineLoops = [LP1]
        
        if createSurface:
            SF1  = Surface([LP1])
            self.Surfaces = [SF1]
        
    def write_txt(self):
        
        txt = "// Creating a Rectangle\n"
        
        # concatenate all the objects
        # Points, Lines, LineLoops, Surfaces, SurfaceLoops, Volumes
        #
        
        objs = self.Points + self.Lines + \
               self.LineLoops + self.Surfaces
               
        
        # write all the objects
        for obj in objs:
            txt+=obj.write_txt()
        
        txt += "// Done creating a Rectangle\n"
            
        return txt
        
    
class Box:
    """
    A class to build a box from points, lines, lineloops, surfaces, 
    """
    def __init__(self, x0, y0, z0, dx, dy, dz, size_tp, size_bt):
        """
        offer two meshes, one at bottom face, one at top face.
        """
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dx = dx
        self.dy = dy
        self.dz = dz
        
        self.size_tp = size_tp
        self.size_bt = size_bt
        
        self.Points        = []
        self.Lines         = []
        self.LineLoops     = []
        self.Surfaces      = []
        self.SurfaceLoops  = []
        self.Volumes       = []
        
    def create_box_geometry(self, createVolume=True):
        x  = self.x0
        y  = self.y0
        z  = self.z0
        dx  = self.dx
        dy  = self.dy
        dz  = self.dz
        size_tp = self.size_tp
        size_bt = self.size_bt
        
        # create points in tbe bottom face(bt)
        p1 = Point(x,    y,    z, size_bt)
        p2 = Point(x+dx, y,    z, size_bt)
        p3 = Point(x+dx, y+dy, z, size_bt)
        p4 = Point(x,    y+dy, z, size_bt)
        
        p5 = Point(x,    y,    z+dz, size_tp)
        p6 = Point(x+dx, y,    z+dz, size_tp)
        p7 = Point(x+dx, y+dy, z+dz, size_tp)
        p8 = Point(x,    y+dy, z+dz, size_tp)
        
        self.Points=[p1, p2, p3, p4, p5, p6, p7, p8]
        
        # create lines:
        
        # bottom face
        lineb1 = Line(p1, p2)
        lineb2 = Line(p2, p3)
        lineb3 = Line(p3, p4)
        lineb4 = Line(p4, p1)

        # top face
        linet1 = Line(p5, p6)
        linet2 = Line(p6, p7)
        linet3 = Line(p7, p8)
        linet4 = Line(p8, p5)

        # 4 vertical lines
        linev1 = Line(p1, p5)
        linev2 = Line(p2, p6)
        linev3 = Line(p3, p7)
        linev4 = Line(p4, p8)
        
        self.Lines =[lineb1,lineb2,lineb3,lineb4,
                     linet1,linet2,linet3,linet4,
                     linev1,linev2,linev3,linev4]
        
        # create LineLoops
        
        # create the lineloop such that 
        # the normal (right hand rule) points outside of the box
        
        # botoom
        LP_bt = LineLoop([lineb1, lineb2, lineb3, lineb4], direction= -1)
        
        # top
        LP_tp = LineLoop([linet1, linet2, linet3, linet4], direction=  1)
        
        # left
        LP_lt = LineLoop([lineb1, linev2, linet1, linev1], direction=  1)
        
        # front
        LP_ft = LineLoop([lineb2, linev3, linet2, linev2], direction=  1)
        
        # right
        LP_rt = LineLoop([lineb3, linev4, linet3, linev3], direction=  1)
        
        # back
        LP_bk = LineLoop([lineb4, linev1, linet4, linev4], direction=  1)
        
        self.LineLoops = [LP_bt, LP_tp, LP_lt, LP_ft,LP_rt,LP_bk]
        
        # create 6 Surfaces
        S_bt = Surface([LP_bt])
        S_tp = Surface([LP_tp])
        S_lt = Surface([LP_lt])
        S_ft = Surface([LP_ft])
        S_rt = Surface([LP_rt])
        S_bk = Surface([LP_bk])
        
        self.Surfaces = [S_bt, S_tp, S_lt, S_ft, S_rt, S_bk]
        normals = [1]*6
        
        # create SurfaceLoops
        SL = SurfaceLoop(self.Surfaces, normals = normals)
        self.SurfaceLoops = [SL]
        
        # create volumes
        if createVolume:
            V  = Volume(self.SurfaceLoops)
            self.Volumes = [V]
    
    def write_txt(self):
        
        txt = "// Creating a box\n"
        
        # concatenate all the objects
        # Points, Lines, LineLoops, Surfaces, SurfaceLoops, Volumes
        #
        
        objs = self.Points + self.Lines + self.LineLoops + \
                  self.Surfaces + self.SurfaceLoops + self.Volumes
        
        # write all the objects
        for obj in objs:
            txt+=obj.write_txt()
        
        txt += "// Done creating a box\n"
            
        return txt

def createCircularSizeFields(cs, rs, sizes, DistMax, DistMin, LcMax, LcMin):
    """
    create a size field near each circle
    so that mesh is refined in the region near circular sizes.
    """
    objs = []
    sizeFieldsEvals = []
    sizeFieldsThres = []
    for c, r, size in zip(cs, rs, sizes):
        # create a MathEvalField
        option_matheval = {'F':"( (x-(%.5f))^2 + (y-(%.5f))^2 )^0.5/%.5f"%(c[0],c[1], r)}
        field_matheval = Field('MathEval', option_matheval)
        sizeFieldsEvals.append(field_matheval)

        option_threshold = {'DistMax': DistMax, 
                            'DistMin': DistMin, 
                            'IField': field_matheval.index, 
                            'LcMax': size*LcMax,
                            'LcMin': size*LcMin,
                            'StopAtDistMax': 1}

        field_threhold = Field('Threshold', option_threshold)
        sizeFieldsThres.append(field_threhold)

    objs +=  sizeFieldsEvals +  sizeFieldsThres
    return objs, sizeFieldsThres
    
def createSphereShellSizeFields(cs, rs, sizes, DistMax, DistMin, LcMax, LcMin):
    """
    create a size field near each sphere 
    so that mesh is refined in the region near sphere sizes.
    """
    objs = []
    sizeFieldsEvals = []
    sizeFieldsThres = []
    for c, r, size in zip(cs, rs, sizes):
        # create a MathEvalField
        option_matheval = {'F':"((x-(%.5f))^2 + (y-(%.5f))^2 + (z-(%.5f))^2)^0.5/%.5f"%(c[0],c[1],c[2], r)}
        field_matheval = Field('MathEval', option_matheval)
        sizeFieldsEvals.append(field_matheval)

        option_threshold = {'DistMax': DistMax, 
                            'DistMin': DistMin, 
                            'IField': field_matheval.index, 
                            'LcMax': size*LcMax,
                            'LcMin': size*LcMin,
                            'StopAtDistMax': 1}

        field_threhold = Field('Threshold', option_threshold)
        sizeFieldsThres.append(field_threhold)

    # create a Min field as set it as background field
    # option_min = {'FieldsList': [i.index for i in sizeFieldsThres]}
    # field_min  = Field('Min', option_min, setbackground=True)

    objs +=  sizeFieldsEvals +  sizeFieldsThres
    return objs, sizeFieldsThres

def createEllipseShellSizeFields2D(cs, rs, axes, sizes, DistMax, DistMin, LcMax, LcMin, StopAtDistMax=1):
    """
    create a size field near ellipsoids 
    so that mesh is refined in the region near ellipsoids.
    """
    objs = []
    sizeFieldsEvals = []
    sizeFieldsThres = []
    
    for c, r, axe, size in zip(cs, rs, axes, sizes):
        # create a MathEvalField
        ellipseExpr = writeMathEvalEllipse2D(c, r, axe)
        
        option_matheval = {'F': ellipseExpr}
        field_matheval = Field('MathEval', option_matheval)
        sizeFieldsEvals.append(field_matheval)

        option_threshold = {'DistMax': DistMax, 
                            'DistMin': DistMin, 
                            'IField' : field_matheval.index, 
                            'LcMax'  : size*LcMax,
                            'LcMin'  : size*LcMin,
                            'StopAtDistMax': StopAtDistMax}

        field_threhold = Field('Threshold', option_threshold)
        sizeFieldsThres.append(field_threhold)

    objs +=  sizeFieldsEvals +  sizeFieldsThres
    return objs, sizeFieldsThres

def createEllipseShellSizeFields(cs, rs, axes, sizes, DistMax, DistMin, LcMax, LcMin, StopAtDistMax=1):
    """
    create a size field near ellipsoids 
    so that mesh is refined in the region near ellipsoids.
    """
    objs = []
    sizeFieldsEvals = []
    sizeFieldsThres = []
    
    for c, r, axe, size in zip(cs, rs, axes, sizes):
        # create a MathEvalField
        ellipseExpr = writeMathEvalEllipse(c, r, axe)
        
        option_matheval = {'F': ellipseExpr}
        field_matheval = Field('MathEval', option_matheval)
        sizeFieldsEvals.append(field_matheval)

        option_threshold = {'DistMax': DistMax, 
                            'DistMin': DistMin, 
                            'IField' : field_matheval.index, 
                            'LcMax'  : size*LcMax,
                            'LcMin'  : size*LcMin,
                            'StopAtDistMax': StopAtDistMax}

        field_threhold = Field('Threshold', option_threshold)
        sizeFieldsThres.append(field_threhold)

    objs +=  sizeFieldsEvals +  sizeFieldsThres
    return objs, sizeFieldsThres

def writeMathEvalEllipse(cs, rs, axes):
    """
    write a MathEval string expression 
    for a transformed distance to ellipsoid centroid
    normalized to a unit sphere
    cs: a list of 3 values, xc, yc, zc
    rs: a list of 3 values, rx, ry, rz
    axes: a list of 3 vectors, ax1, ax2, ax3, normal vectors
    """
    xc, yc, zc   = cs
    rx, ry, rz   = rs
    
    # normalize the axe vectors
    for ax in axes:
        norm = sum([v**2 for v in ax])
        norm = np.sqrt(norm)
        ax   = [v/norm for v in ax]
    
    ax1,ax2, ax3 = axes
    
    txt_x = "((x-(%.6f))*(%.6f) + (y-(%.6f))*(%.6f) + (z-(%.6f))*(%.6f))^2/%.6f^2"%(
                  xc,    ax1[0],       yc,   ax1[1],       zc,   ax1[2],    rx)
    txt_y = "((x-(%.6f))*(%.6f) + (y-(%.6f))*(%.6f) + (z-(%.6f))*(%.6f))^2/%.6f^2"%(
                  xc,    ax2[0],       yc,   ax2[1],       zc,   ax2[2],    ry)
    txt_z = "((x-(%.6f))*(%.6f) + (y-(%.6f))*(%.6f) + (z-(%.6f))*(%.6f))^2/%.6f^2"%(
                  xc,    ax3[0],       yc,   ax3[1],       zc,   ax3[2],    rz)
    expr   = "(" + txt_x + "+" + txt_y + "+" + txt_z + ')^0.5'
    
    return expr


def writeMathEvalEllipse2D(cs, rs, axes):
    """
    write a MathEval string expression 
    for a transformed distance to 2D ellipse centroid
    normalized to a unit sphere
    cs: a list of 2 values, xc, yc
    rs: a list of 2 values, rx, ry
    axes: a list of 2 vectors, ax1, ax2, normal vectors
    """
    xc, yc   = cs
    rx, ry   = rs
    
    # normalize the axe vectors
    for ax in axes:
        norm = sum([v**2 for v in ax])
        norm = np.sqrt(norm)
        ax   = [v/norm for v in ax]
    
    ax1,ax2 = axes
    
    txt_x = "((x-(%.6f))*(%.6f) + (y-(%.6f))*(%.6f))^2/%.6f^2"%(
                  xc,    ax1[0],       yc,   ax1[1],    rx)
    txt_y = "((x-(%.6f))*(%.6f) + (y-(%.6f))*(%.6f) )^2/%.6f^2"%(
                  xc,    ax2[0],       yc,   ax2[1],    ry)
    expr   = "(" + txt_x + "+" + txt_y +")^0.5"
    
    return expr

def createRectangleCircles(L, H, size_bc, rs, cs, sizes, 
                           DistMax=3, DistMin=1, LcMax=3, 
                           LcMin=1, geofile=None):
    """
    
    create a mesh with rectangle boundary and circular inclusion
    
    """
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    objs = []
    
    # create a rectangle without creating the surface
    x0   = -L/2.
    y0   = -H/2.
    dx   = L
    dy   = H
    
    Rect = Rectangle(x0, y0, dx, dy, size_bc, z0 = 0)
    Rect.create_rectangle_geometry(createSurface=False)
    
    objs += [Rect]
    
    # create multiple circular inclusions
    LP   = [Rect.LineLoops[0]]
    
    for c, r, size in zip(cs, rs, sizes):
        CRC_i = Circle(c[0], c[1], r, size, zc=0)
        CRC_i.create_circle_geometry()
        objs.append(CRC_i)
        LP.append(CRC_i.LineLoops[0])
    
    # create surface with rectangular outboundary and multile inclusions
    SF = Surface(LP)
    objs += [SF]
    
    # create size fields
    fields, fieldsThres = createCircularSizeFields(cs, rs, sizes, DistMax, DistMin, LcMax, LcMin)
    objs += fields
        
    # create a min field
    option_min = {'FieldsList': [i.index for i in fieldsThres]}
    field_min  = Field('Min', option_min, setbackground=True)
    
    objs += [field_min]
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs
    

def createTwoLayers_Cylinder(L, H, Z0, size_bt, size_md, size_tp, 
                    DistMax, DistMin, LcMax, LcMin, StopAtDistMax=1, 
                    geofile=None):
    """
    Create a mesh with two layers conforming at the interface Z0
    """

    # reset all the indexes
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    dZ1 = H + Z0
    dZ2 = -Z0

    X,  Y, Z = -L/2.0, -L/2.0, -H
    dX,  dY  = L, L

    objs = []

    # create bottom layer without creating volume
    BX1    = Box(X, Y, Z, dX, dY, dZ1, size_md, size_bt)
    BX1.create_box_geometry(createVolume = True)

    objs.append(BX1)
    
    # create a Box from a bottom box.
    objs += createBoxFromBtBox(BX1, dZ2, size_tp)
    
    # 
    meval_opt = {'F': '(x^2+y^2)^0.5'}
    sf_meval = Field('MathEval', meval_opt)
    
    option_threshold = {'DistMax': DistMax, 
                        'DistMin': DistMin, 
                        'IField' : sf_meval.index, 
                        'LcMax'  : LcMax,
                        'LcMin'  : LcMin,
                        'StopAtDistMax': StopAtDistMax}
    
    sf_thre = Field('Threshold', option_threshold, setbackground=True)
    
    objs += [sf_meval, sf_thre]
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs

def createMultiLayersBands(L, H, Z0, dZ1, dZ2, size_bt, size_md, size_tp, thetas,
                           DistMax_up, DistMin_up, LcMax_up, LcMin_up,
                           DistMax_lo, DistMin_lo, LcMax_lo, LcMin_lo,
                           DistMax, DistMin, LcMax, LcMin,
                           geofile = None):
    """
    Create a mesh with two layers conforming at the interface Z0
    """

    # reset all the indexes
    objs = []

    # create multiple layer geometry
    zs   = list(np.arange(0, Z0, -dZ1)) + list(np.arange(Z0, -H-dZ2, -dZ2))
    
    # create all the points
    pts  = []

    for z in zs:
        if z>Z0:
            pt_size = size_tp
        elif z<Z0:
            pt_size = size_bt
        else:
            pt_size = size_md

        p1 = Point(-L/2., -L/2., z, pt_size) 
        p2 = Point(+L/2., -L/2., z, pt_size) 
        p3 = Point(+L/2., +L/2., z, pt_size) 
        p4 = Point(-L/2., +L/2., z, pt_size) 
        pts.append([p1, p2, p3, p4])
        objs += [p1, p2, p3, p4]
    
    lines_h = []
    lines_v = []
    lineloops_h = []
    surfaces_h  = []
    N = len(pts)

    for i in range(N):
        # create horizontal lines
        lines_h_i = []
        pts_i = pts[i]
        npt_i  = len(pts_i)
        for j in range(npt_i):
            if j<npt_i-1:
                lj = Line(pts_i[j], pts_i[j+1])
            else:
                lj = Line(pts_i[j], pts_i[0])
            lines_h_i.append(lj)
        lines_h.append(lines_h_i)
        lloop = LineLoop(lines_h_i)
        lineloops_h.append(lloop)
        s_h   = Surface([lloop], isPlaneSurface=True)
        surfaces_h.append(s_h)
        objs += lines_h_i
        objs.append(lloop)
        objs.append(s_h)
        
        if i==N-1:
            continue
        lines_v_i = []
        pts_i_p   = pts[i+1]
        for j in range(npt_i):
            lines_v_i.append(Line(pts_i[j], pts_i_p[j]))
        lines_v.append(lines_v_i)
        objs += lines_v_i

    # create vertical lineloops, surfaces and volumes
    for i in range(N-1):
        npt = len(pts[i])
        surface_list = [surfaces_h[i], surfaces_h[i+1]]
        for j in range(npt):
            if j<npt-1:
                llooplist= [lines_h[i][j], lines_v[i][j+1], lines_h[i+1][j], lines_v[i][j]]
            else:
                llooplist= [lines_h[i][j], lines_v[i][0], lines_h[i+1][j], lines_v[i][j]]
            lloop = LineLoop(llooplist)
            objs.append(lloop)
            s_v   = Surface([lloop], isPlaneSurface=True)
            objs.append(s_v)
            surface_list.append(s_v)
        # create surface loop
        sloop = SurfaceLoop(surface_list)
        objs.append(sloop)
        # create a volume
        vol = Volume([sloop])
        objs.append(vol)

    # create shear band refinement.
    sfs_thres = []
    for theta in thetas:
        
        F_up, F_lo = writeDistToPlaneTwoLayers(theta, Z0)
        
        # upper layer
        meval_opt_up = {'F':F_up}
        sf_meval_up = Field('MathEval', meval_opt_up)
        
        objs.append(sf_meval_up)
        option_threshold_up = {'DistMax': DistMax_up, 
                            'DistMin': DistMin_up, 
                            'IField' : sf_meval_up.index, 
                            'LcMax'  : LcMax_up,
                            'LcMin'  : LcMin_up,
                            'StopAtDistMax': 1}
        sf_thre_up = Field('Threshold', option_threshold_up)
        sfs_thres.append(sf_thre_up)
        
        # lower layer
        meval_opt_lo = {'F':F_lo}
        sf_meval_lo = Field('MathEval', meval_opt_lo)
        
        objs.append(sf_meval_lo)
        
        option_threshold_lo = {'DistMax': DistMax_lo, 
                            'DistMin': DistMin_lo, 
                            'IField' : sf_meval_lo.index, 
                            'LcMax'  : LcMax_lo,
                            'LcMin'  : LcMin_lo,
                            'StopAtDistMax': 1}
        sf_thre_lo = Field('Threshold', option_threshold_lo)
        sfs_thres.append(sf_thre_lo)
    
    # add a cylinder
    meval_opt = {'F': '(x^2+y^2)^0.5'}
    sf_meval = Field('MathEval', meval_opt)
    
    option_threshold = {'DistMax': DistMax, 
                        'DistMin': DistMin, 
                        'IField' : sf_meval.index, 
                        'LcMax'  : LcMax,
                        'LcMin'  : LcMin,
                        'StopAtDistMax': 1}
    
    sf_thre = Field('Threshold', option_threshold)
    objs.append(sf_meval)
    
    sfs_thres.append(sf_thre)
    
    objs += sfs_thres
    # compute a min of all these
    option_min = {'FieldsList': [i.index for i in sfs_thres]}
    field_min  = Field('Min', option_min, setbackground=True)
    
    objs += [field_min]
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs

def createTwoLayersBands(L, H, Z0, size_bt, size_md, size_tp, thetas,
                         DistMax_up, DistMin_up, LcMax_up, LcMin_up,
                         DistMax_lo, DistMin_lo, LcMax_lo, LcMin_lo,
                         DistMax, DistMin, LcMax, LcMin,
                         geofile = None):
    """
    Create a mesh with two layers conforming at the interface Z0
    """

    # reset all the indexes
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    dZ1 = H + Z0
    dZ2 = -Z0

    X,  Y, Z = -L/2.0, -L/2.0, -H
    dX,  dY  = L, L

    objs = []

    # create bottom layer without creating volume
    BX1    = Box(X, Y, Z, dX, dY, dZ1, size_md, size_bt)
    BX1.create_box_geometry(createVolume = True)

    objs.append(BX1)
    
    # create a Box from a bottom box.
    objs += createBoxFromBtBox(BX1, dZ2, size_tp)
    
    # create shear band refinement.
    sfs_thres = []
    for theta in thetas:
        
        F_up, F_lo = writeDistToPlaneTwoLayers(theta, Z0)
        
        # upper layer
        meval_opt_up = {'F':F_up}
        sf_meval_up = Field('MathEval', meval_opt_up)
        
        objs.append(sf_meval_up)
        option_threshold_up = {'DistMax': DistMax_up, 
                            'DistMin': DistMin_up, 
                            'IField' : sf_meval_up.index, 
                            'LcMax'  : LcMax_up,
                            'LcMin'  : LcMin_up,
                            'StopAtDistMax': 1}
        sf_thre_up = Field('Threshold', option_threshold_up)
        sfs_thres.append(sf_thre_up)
        
        # lower layer
        meval_opt_lo = {'F':F_lo}
        sf_meval_lo = Field('MathEval', meval_opt_lo)
        
        objs.append(sf_meval_lo)
        
        option_threshold_lo = {'DistMax': DistMax_lo, 
                            'DistMin': DistMin_lo, 
                            'IField' : sf_meval_lo.index, 
                            'LcMax'  : LcMax_lo,
                            'LcMin'  : LcMin_lo,
                            'StopAtDistMax': 1}
        sf_thre_lo = Field('Threshold', option_threshold_lo)
        sfs_thres.append(sf_thre_lo)
    
    # add a cylinder
    meval_opt = {'F': '(x^2+y^2)^0.5'}
    sf_meval = Field('MathEval', meval_opt)
    
    option_threshold = {'DistMax': DistMax, 
                        'DistMin': DistMin, 
                        'IField' : sf_meval.index, 
                        'LcMax'  : LcMax,
                        'LcMin'  : LcMin,
                        'StopAtDistMax': 1}
    
    sf_thre = Field('Threshold', option_threshold)
    objs.append(sf_meval)
    
    sfs_thres.append(sf_thre)
    
    objs += sfs_thres
    # compute a min of all these
    option_min = {'FieldsList': [i.index for i in sfs_thres]}
    field_min  = Field('Min', option_min, setbackground=True)
    
    objs += [field_min]
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs


def writeDistToPlaneTwoLayers(theta, z0, c=[0,0]):
    """
    write distance to a plane but distinguish the upper/lower layer 
    
    for meshing purposes.
    distance is seperated defined for upper and lower layer
    
    """
    theta  = np.pi*theta/180.
    nx     = -np.sin(theta)
    ny     =  np.cos(theta)
    x0, y0 = c[0], c[1]
    
    dist       = "Abs((x-%.6f)*(%.6f) + (y-%.6f)*(%.6f))"%(x0, nx, y0, ny)
    
    # select the region 
    # if z>z0 return 99999999
    # if z<0 return dist
    
    ind_upper   = "(((z  - (%.6f) + 0.000001)/Abs(z  - (%.6f) + 0.000001) + 1)/2)"%(z0, z0)
    ind_lower   = "((((%.6f)  - z + 0.000001)/Abs((%.6f) - z  + 0.000001) + 1)/2)"%(z0, z0)
    
    dist_z     = "Abs(z-(%.6f))"%(z0)
    
    F_up  = "%s + %s * %s"%(dist, dist_z, ind_lower)
    F_lo  = "%s + %s * %s"%(dist, dist_z, ind_upper)
    
    return F_up, F_lo

def writeDistToPlane(theta, c=[0,0]):
    
    """
    write a distance to a plane with azimuth theta in degree
    counter clockwise from positive x
    """
    
    theta  = np.pi*theta/180.
    nx     = -np.sin(theta)
    ny     =  np.cos(theta)
    x0, y0 = c[0], c[1]
    F      = "Abs((x-%.6f)*(%.6f) + (y-%.6f)*(%.6f))"%(x0, nx, y0, ny)
    
    return F

def AxesVerticalPenny(theta):
    """
    return axes for vertical pennyshape crack.
    
    theta: angle from x in degrees
    """
    pi = 3.1415926535897932384626433
    theta = theta*pi/180.
    ax1   = [ np.cos(theta), np.sin(theta), 0]
    ax2   = [-np.sin(theta), np.cos(theta), 0]
    ax3   = [0, 0, 1]
    return [ax1, ax2, ax3]

def AxesEllipse2D(theta):
    """
    return axes for 2D ellipse.
    
    theta: angle from x in degrees of the major axis
    """
    theta = theta*np.pi/180.
    ax1   = [ np.cos(theta), np.sin(theta)]
    ax2   = [-np.sin(theta), np.cos(theta)]
    return [ax1, ax2]

def createBoxSpheres(L, H, size_bt, size_tp, cs, rs, sizes, 
                     DistMax=3, DistMin=1, LcMax=3, LcMin=1, geofile=None):
    """
    Create geometry: One box, multiple spheres
    Example:
    # two layer geometry
    L   =  100.
    H   =  40.

    # size of bottom, middle and top layer
    size_bt  = 2
    size_md  = 2
    size_tp  = 2

    # spherical inclusions
    cs     = [[0.5*L, 0.5*L, Z0-7],[0.75*L, 0.75*L, Z0-10]]
    rs     = [5, 5]
    sizes  = [0.5, 0.5]
    """
    
    # reset all the indexes
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    dZ1 = H

    X,  Y, Z = -L/2.0, -L/2.0, -H
    dX,  dY  = L, L

    objs = []

    # create bottom layer without creating volume
    BX1    = Box(X, Y, Z, dX, dY, dZ1, size_tp, size_bt)
    BX1.create_box_geometry(createVolume = False)

    objs.append(BX1)

    SurfaceLoops1 = [BX1.SurfaceLoops[0]]

    SPHs   = []
    for (c, r, size) in zip(cs, rs, sizes):
        SPHs.append(Sphere(c[0], c[1], c[2], r, size))
        SPHs[-1].create_sphere_geometry()
        SurfaceLoops1.append(SPHs[-1].SurfaceLoops[0])

    # create a volume with a spherical hole
    V  = Volume(SurfaceLoops1)

    objs += SPHs
    objs += [V]
    
    # write size fields.
    SFs, SFThres= createSphereShellSizeFields(cs, rs, sizes, DistMax, DistMin, LcMax, LcMin)
    
    objs += SFs
    
    # create a min field
    option_min = {'FieldsList': [i.index for i in SFThres]}
    field_min  = Field('Min', option_min, setbackground=True)
    
    objs += [field_min]

    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs


def createTwoLayersSpheres(L, H, Z0, size_bt, size_md, size_tp, cs, rs, sizes, 
                           DistMax=3, DistMin=1, LcMax=3, LcMin=1, geofile=None):
    """
    Create geometry: two layers, multiple spheres in the lower layers
    Example:
    # two layer geometry
    L   =  100.
    H   =  40.
    Z0  = -15.

    # size of bottom, middle and top layer
    size_bt  = 2
    size_md  = 2
    size_tp  = 2

    # spherical inclusions
    cs     = [[0.5*L, 0.5*L, Z0-7],[0.75*L, 0.75*L, Z0-10]]
    rs     = [5, 5]
    sizes  = [0.5, 0.5]
    """
    
    # reset all the indexes
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    dZ1 = H + Z0
    dZ2 = -Z0

    X,  Y, Z = -L/2.0, -L/2.0, -H
    dX,  dY  = L, L

    objs = []

    # create bottom layer without creating volume
    BX1    = Box(X, Y, Z, dX, dY, dZ1, size_md, size_bt)
    BX1.create_box_geometry(createVolume = False)

    objs.append(BX1)

    SurfaceLoops1 = [BX1.SurfaceLoops[0]]

    SPHs   = []
    for (c, r, size) in zip(cs, rs, sizes):
        SPHs.append(Sphere(c[0], c[1], c[2], r, size))
        SPHs[-1].create_sphere_geometry()
        SurfaceLoops1.append(SPHs[-1].SurfaceLoops[0])

    # create a volume with a spherical hole
    V  = Volume(SurfaceLoops1)

    objs += SPHs
    objs += [V]

    # create a Box from a bottom box.
    objs_tp = createBoxFromBtBox(BX1, dZ2, size_tp)
    objs += objs_tp
    
    # write size fields.
    fields= createSphereShellSizeFields(cs, rs, sizes, DistMax, DistMin, LcMax, LcMin)
    
    objs += fields
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs, fields

def createTwoLayersSpheresNoFields(L, H, Z0, size_bt, size_md, size_tp, cs, rs, sizes, geofile=None):
    """
    Create geometry: two layers, multiple spheres in the lower layers
    Example:
    # two layer geometry
    L   =  100.
    H   =  40.
    Z0  = -15.

    # size of bottom, middle and top layer
    size_bt  = 2
    size_md  = 2
    size_tp  = 2

    # spherical inclusions
    cs     = [[0.5*L, 0.5*L, Z0-7],[0.75*L, 0.75*L, Z0-10]]
    rs     = [5, 5]
    sizes  = [0.5, 0.5]
    """
    
    # reset all the indexes
    Point.reset_index()
    Line.reset_index()
    LineLoop.reset_index()
    Surface.reset_index()
    SurfaceLoop.reset_index()
    Volume.reset_index()
    Field.reset_index()
    
    dZ1 = H + Z0
    dZ2 = -Z0

    X,  Y, Z = -L/2.0, -L/2.0, -H
    dX,  dY  = L, L

    objs = []

    # create bottom layer without creating volume
    BX1    = Box(X, Y, Z, dX, dY, dZ1, size_md, size_bt)
    BX1.create_box_geometry(createVolume = False)

    objs.append(BX1)

    SurfaceLoops1 = [BX1.SurfaceLoops[0]]

    SPHs   = []
    for (c, r, size) in zip(cs, rs, sizes):
        SPHs.append(Sphere(c[0], c[1], c[2], r, size))
        SPHs[-1].create_sphere_geometry()
        SurfaceLoops1.append(SPHs[-1].SurfaceLoops[0])

    # create a volume with a spherical hole
    V  = Volume(SurfaceLoops1)

    objs += SPHs
    objs += [V]

    # create a Box from a bottom box.
    objs_tp = createBoxFromBtBox(BX1, dZ2, size_tp)
    objs += objs_tp
    
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
    
    return objs

def createTwoLayerEllipses(L, H, Z0, size_bt, size_md, size_tp, cs, rs, sizes, angles, 
                           DistMax, DistMin, LcMax, LcMin, 
                           R_in, DistMax_LD, DistMin_LD, LcMax_LD, LcMin_LD, sizes_LD, 
                           geofile=None):
    axes   = [AxesVerticalPenny(angle) for angle in angles]
    
    # one sphere is added to enable the meshes at ellipsoidal inclusion.
    cs_sphere    = [[0, 0, cs[0][2]]]
    rs_sphere    = [rs[0][0]]
    sizes_sphere = [size_bt]
    
    objs = createTwoLayersSpheresNoFields(L, H, Z0, size_bt, 
                                          size_md, size_tp, cs_sphere, 
                                          rs_sphere, sizes_sphere, geofile=None)
    SFs, SFThres = createEllipseShellSizeFields(cs, rs, axes, sizes, 
                                                DistMax, DistMin, LcMax, LcMin)
    objs += SFs
    
    rs_LD      = [[R_in]*3]
    cs_LD      = [[0, 0, -H/2.0]]
    axe3       = AxesVerticalPenny(0)
    axes_LD    = [axe3]
    
    SFs_LD, SFThres_LD = createEllipseShellSizeFields(cs_LD, rs_LD, 
                                                      axes_LD, sizes_LD, 
                                                      DistMax_LD, DistMin_LD, 
                                                      LcMax_LD, LcMin_LD)
    objs += SFs_LD
    # create a min field
    option_min = {'FieldsList': [i.index for i in SFThres+SFThres_LD]}
    field_min  = Field('Min', option_min, setbackground=True)
    objs += [field_min]
    if geofile:
        f = open(geofile,'w')
        for obj in objs:
            f.write(obj.write_txt())
        f.close()
        
    return objs
    
    
def createBoxFromBtBox(BX1, H, size_tp):
    """
    create new a box from a bottom box given the fact that the 
    bottom face of the new box has already been created by the lower box.
    
    return a list of new objects
    """
    
    objs = []
    
    # create another box on top of the first block
    p1, p2, p3, p4 = BX1.Points[4:]

    # create another 4 points 
    p5 = Point(p1.x, p1.y, p1.z+H, size_tp)
    p6 = Point(p2.x, p2.y, p2.z+H, size_tp)
    p7 = Point(p3.x, p3.y, p3.z+H, size_tp)
    p8 = Point(p4.x, p4.y, p4.z+H, size_tp)

    newpts = [p5, p6, p7, p8]
    objs  += newpts
    
    # read lines for the bottom face
    lineb1, lineb2, lineb3, lineb4 = BX1.Lines[4:8]

    # top face
    linet1 = Line(p5, p6)
    linet2 = Line(p6, p7)
    linet3 = Line(p7, p8)
    linet4 = Line(p8, p5)

    # 4 vertical lines
    linev1 = Line(p1, p5)
    linev2 = Line(p2, p6)
    linev3 = Line(p3, p7)
    linev4 = Line(p4, p8)

    newlines = [linet1, linet2, linet3, linet4, 
                linev1, linev2, linev3, linev4] 

    objs  += newlines

    # create line loops
    # create the lineloop such that 
    # the normal (right hand rule) points outside of the box

    # botoom was created, the top surface of BX1, the second Surface
    LP_bt = BX1.LineLoops[1]

    # top
    LP_tp = LineLoop([linet1, linet2, linet3, linet4], direction=  1)

    # left
    LP_lt = LineLoop([lineb1, linev2, linet1, linev1], direction=  1)

    # front
    LP_ft = LineLoop([lineb2, linev3, linet2, linev2], direction=  1)

    # right
    LP_rt = LineLoop([lineb3, linev4, linet3, linev3], direction=  1)

    # back
    LP_bk = LineLoop([lineb4, linev1, linet4, linev4], direction=  1)

    newLineLoops = [LP_tp, LP_lt, LP_ft,LP_rt,LP_bk]
    
    objs  += newLineLoops

    # create new surfaces

    # bottom surface was created by BX1, the second Surface
    S_bt = BX1.Surfaces[1]
    
    S_tp = Surface([LP_tp])
    S_lt = Surface([LP_lt])
    S_ft = Surface([LP_ft])
    S_rt = Surface([LP_rt])
    S_bk = Surface([LP_bk])

    Surfaces    = [S_bt, S_tp, S_lt, S_ft, S_rt, S_bk]
    newSurfaces = [S_tp, S_lt, S_ft, S_rt, S_bk]
    objs  += newSurfaces

    # define a surface loop
    normals    = [1]*6
    normals[0] = -1    # revert the normal of the bottom face
    SL_BOX = SurfaceLoop(Surfaces, normals = normals)
    
    newSurfaceLoops = [SL_BOX]
    objs += newSurfaceLoops

    # Create a new volume for the top layer
    V  = Volume([SL_BOX])
    objs += [V]
    
    return objs
