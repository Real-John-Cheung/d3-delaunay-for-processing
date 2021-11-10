# d3-delaunay-for-processing
To port the [d3-delaunay library](https://github.com/d3/d3-delaunay) to Java
using the [delaunator-java library by waveware4ai](https://github.com/waveware4ai/delaunator-java)

<hr>

## Usage
You can download the packaged <code>.jar</code> file from the release, then import them into processing as a library. (you might need the newest Processing 4 to use it)
Or you can also build it from source. As you can see, I am using gradle to package it. To do so:<br>
clone the repo to your machine:<br><code> git clone https://github.com/Real-John-Cheung/d3-delaunay-for-processing.git </code><br>
after that, you can make some change according to your needs and run
<code> gradle clean shadowJar</code> 
to build the library.
It will be in <code>./lib/build/libs</code>

<hr>

## APIs
This library, just like the JS one, provides two new class.

## Delaunay

### <i>static</i> Delaunay.from(double[][] points)

This method return a new Delaunay instant from a double[][] array like this [ [x0, y0] , [ x1, y1] , ...]

<b>Return type: Delaunay</b>

### <i>constructor</i> Delaunay(double[] points) 

This is the recommanded way to create a new Delaunay instant. It accepts a 2D double array like this[x0, y0, x1, y1 , ...]

<b>Return type: Delaunay</b>

### <i>attribute</i> delaunay.points

An 1D double array of all the points input. Like this [ x0, y0, x1, y1, ...]

<b>Type: double[]</b>

### <i>attribute</i> delaunay.halfedges

The halfedge indexes as an int array: [ j0, j1, ...]

<b>Type: int[]</b>

### <i>attribute</i> delaunay.hull

An int array of point indexes that form the convex hull in counterclockwise order. If the points are collinear, returns them ordered.

<b>Type: int[]</b>

### <i>attribute</i> delaunay.triangles

The triangle vertex indexes as an int array [ i0, j0, k0, i1, j1, k1, ...]

<b>Type: int[]</b>

### <i>attribute</i> delaunay.inedges

The incoming halfedge indexes as a int array [ e0, e1, e2, ...]. For each point i, inedges[i] is the halfedge index e of an incoming halfedge. For coincident points, the halfedge index is -1; for points on the convex hull, the incoming halfedge is on the convex hull; for other points, the choice of incoming halfedge is arbitrary.

<b>Type: int[]</b>

### <i>method</i> delaunay.find(double x, double y, int i)

Returns the index of the input point that is closest to the specified point (x, y). The search is started at the specified point i. If i is not specified, it defaults to zero.

<b>Return type: int</b>

### <i>method</i> delaunay.neighbors(int i)

Returns an int array of the indexes of the neighboring points to the specified point i. 

<b>Return type: int[]</b>

### <i>method</i> delaunay.hullPolygon()

Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], ... , [ x0, y0 ]] representing the convex hull.

<b>Return type: double[][]</b>

### <i>method</i> delaunay.trianglePolygon(int i)

Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ], [ x0, y0 ] ] representing the triangle i.

<b>Return type: double[][]</b>

### <i>method</i> delaunay.trianglePolygons()

Returns all triangles in an array.

<b>Return type: double[][][]</b>

### <i>method</i> delaunay.voronoi(double[] bounds)

Returns the corresponding Voronoi instant.The diagram will be clipped to the specified bounds = [ xmin, ymin, xmax, ymax ]

<b>Return type: Voronoi</b>

<hr>

## Voronoi

### <i>attribute</i> voronoi.delaunay

Reference back to the corresponding delaunay instant.

<b>Type: Delaunay</b>

### <i>attribute</i> voronoi.circumcenters

The circumcenters of the Delaunay triangles as a double array [ cx0, cy0, cx1, cy1, ... ]

<b>Type: double[]</b>

### <i>attribute</i> voronoi.vectors

A double array [ vx0, vy0, wx0, wy0, ... ] where each non-zero quadruple describes an open (infinite) cell on the outer hull, giving the directions of two open half-lines.

<b>Type: double[]</b>

### <i>method</i> voronoi.contains(int i, double x, double y)

Returns true if the cell with the specified index i contains the specified point (x, y)

<b>Return type: boolean</b>

### <i>method</i> voronoi.neighbors(int i)

Returns an int array of the indexes of the cells that share a common edge with the specified cell i. 

<b>Return type: int[]</b>

### <i>method</i> voronoi.cellPolygon(int i)

Returns the convex, closed polygon [ [ x0, y0 ], [ x1, y1 ], ..., [ x0, y0 ] ] representing the cell for the specified point i.

<b>Return type: double[][]</b>





















