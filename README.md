# GmshPyUtils
A python module that contains helper functions and classes to write Gmsh code programmatically. Although pygmsh does a better job, this module is a tool I developed before I knew pygmsh existed.

This module automatically track the numbering of basic gmsh objects, points, lines, lineloops, surfaces, volumes, fields. Each object contains the same method "write_txt" to generate relevant gmsh code piece.

It also contains simple geometric objects such as rectangle, box, ellipse, sphere, ellipsoids, built using a bottom-up approach from the default kernel.
