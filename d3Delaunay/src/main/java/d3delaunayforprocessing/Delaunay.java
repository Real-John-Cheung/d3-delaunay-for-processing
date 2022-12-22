
package d3delaunayforprocessing;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Comparator;
import org.waveware.delaunator.Delaunator;
import org.waveware.delaunator.DTriangle;
import org.waveware.delaunator.DPoint;
import org.checkerframework.checker.signedness.qual.Unsigned;
import org.waveware.delaunator.DEdge;
import d3delaunayforprocessing.Polygon;
import d3delaunayforprocessing.Voronoi;

/**
 * d3-delaunay-for-processing
 * To port the d3-delaunay library to Java using the delaunator-java library by
 * waveware4ai
 * This is the Delaunay class
 * 
 * @author John C
 * @version 1.1
 * @since 2022-12-22
 */

public class Delaunay {
    public static double tau = (double) (2 * Math.PI);
    public Delaunator _delaunator;
    /**
     * The incoming halfedge indexes as a int array [ e0, e1, e2, ...]. For each
     * point i, inedges[i] is the halfedge index e of an incoming halfedge. For
     * coincident points, the halfedge index is -1; for points on the convex hull,
     * the incoming halfedge is on the convex hull; for other points, the choice of
     * incoming halfedge is arbitrary.
     */
    public int[] indeges;
    public int[] _hullIndex;
    /**
     * An 1D double array of all the points input. Like this [ x0, y0, x1, y1, ...]
     */
    public double[] points;
    public int[] collinear;
    /**
     * The halfedge indexes as an int array: [ j0, j1, ...]
     */
    public int[] halfedges;
    /**
     * An int array of point indexes that form the convex hull in counterclockwise
     * order. If the points are collinear, returns them ordered.
     */
    public int[] hull;
    /**
     * The triangle vertex indexes as an int array [ i0, j0, k0, i1, j1, k1, ...]
     */
    public int[] triangles;

    private double pointX(double[] p) {
        return p[0];
    }

    private double pointY(double[] p) {
        return p[1];
    }

    private boolean collinearF(Delaunator d) {
        int[] triangles = d.triangles;
        DPoint[] points = d.points;
        double[] coords = new double[2 * points.length];
        for (int i = 0; i < points.length; i++) {
            DPoint p = points[i];
            coords[2 * i] = p.x;
            coords[2 * i + 1] = p.y;
        }
        for (int i = 0; i < triangles.length; i += 3) {
            int a = 2 * triangles[i];
            int b = 2 * triangles[i + 1];
            int c = 2 * triangles[i + 2];
            double cross = (coords[c] - coords[a]) * (coords[b + 1] - coords[a + 1])
                    - (coords[b] - coords[a]) * (coords[c + 1] - coords[a + 1]);
            if (cross > 1e-10)
                return false;
        }
        return true;
    }

    private double[] jitter(double x, double y, double r) {
        return new double[] { x + Math.sin(x + y), y + Math.cos(x - y) * r };
    }

    /**
     * This method return a new Delaunay instant from a double[][] array like this [
     * [x0, y0] , [ x1, y1] , ...].
     * 
     * @param points the input points
     * @return a Delaunay Object
     */

    static public Delaunay from(double[][] points) {
        double[] flatArr = new double[points.length * 2];
        for (int i = 0; i < points.length; i++) {
            flatArr[i * 2] = points[i][0];
            flatArr[i * 2 + 1] = points[i][1];
        }
        return new Delaunay(flatArr);
    }

    /**
     * This method return a new Delaunay instant from a double[][] array like this [
     * [x0, y0] , [ x1, y1] , ...]. Input with float.
     * 
     * @param points the input points
     * @return a Delaunay Object
     */

    static public Delaunay fromFloat(float[][] points) {
        double[] flatArr = new double[points.length * 2];
        for (int i = 0; i < points.length; i++) {
            flatArr[i * 2] = (double) points[i][0];
            flatArr[i * 2 + 1] = (double) points[i][1];
        }
        return new Delaunay(flatArr);
    }

    /**
     * This is the recommanded way to create a new Delaunay instant. It accepts a 2D
     * double array like this[x0, y0, x1, y1 , ...]
     * 
     * @param points the input points
     */

    public Delaunay(double[] points) {
        DPoint[] dps = new DPoint[points.length / 2];
        for (int i = 0; i < dps.length; i++) {
            dps[i] = new DPoint(points[i * 2], points[i * 2 + 1]);
        }
        this._delaunator = new Delaunator(dps);
        this.indeges = new int[points.length / 2];
        this._hullIndex = new int[points.length / 2];
        DPoint[] dPoints = this._delaunator.points;
        this.points = new double[dPoints.length * 2];
        for (int i = 0; i < dPoints.length; i++) {
            this.points[2 * i] = dPoints[i].x;
            this.points[2 * i + 1] = dPoints[i].y;
        }
        this._init();
    }

    /**
     * This is the recommanded way to create a new Delaunay instant. It accepts a 2D
     * double array like this[x0, y0, x1, y1 , ...]. Float version
     * 
     * @param points the input points
     */

    public Delaunay(float[] points) {
        double[] da = new double[points.length];
        for (int i = 0; i < da.length; i++) {
            da[i] = (double) points[i];
        }
        DPoint[] dps = new DPoint[da.length / 2];
        for (int i = 0; i < dps.length; i++) {
            dps[i] = new DPoint(da[i * 2], da[i * 2 + 1]);
        }
        this._delaunator = new Delaunator(dps);
        this.indeges = new int[da.length / 2];
        this._hullIndex = new int[da.length / 2];
        DPoint[] dPoints = this._delaunator.points;
        this.points = new double[dPoints.length * 2];
        for (int i = 0; i < dPoints.length; i++) {
            this.points[2 * i] = dPoints[i].x;
            this.points[2 * i + 1] = dPoints[i].y;
        }
        this._init();
    }

    /*
     * public Delaunay update() {
     * //not avaliable now
     * return this;
     * }
     */

    public void _init() {
        Delaunator d = this._delaunator;
        double[] points = this.points;
        // check for collinear
        if (d.hull != null && d.hull.length > 2 && this.collinearF(d)) {
            Integer[] collinearWIP = new Integer[points.length / 2];
            for (int i = 0; i < collinearWIP.length; i++) {
                collinearWIP[i] = i;
            }
            Arrays.sort(collinearWIP, new Comparator<Integer>() {
                @Override
                public int compare(Integer i, Integer j) {
                    if (points[2 * i] - points[2 * j] != 0)
                        return points[2 * i] > points[2 * j] ? 1 : 0;
                    return points[2 * i + 1] > points[2 * j + 1] ? 1 : 0;
                }
            });
            this.collinear = Arrays.asList(collinearWIP).stream().mapToInt(i -> i).toArray();
            int e = this.collinear[0];
            int f = this.collinear[this.collinear.length - 1];
            double[] bounds = new double[] { points[2 * e], points[2 * e + 1], points[2 * f], points[2 * f + 1] };
            double r = 1e-8 * Math.hypot(bounds[3] - bounds[1], bounds[2] - bounds[0]);
            for (int i = 0; i < points.length / 2; ++i) {
                double[] p = this.jitter(points[2 * i], points[2 * i + 1], r);
                points[2 * i] = p[0];
                points[2 * i + 1] = p[1];
            }
            DPoint[] newDPoints = new DPoint[points.length / 2];
            for (int i = 0; i < newDPoints.length; i++) {
                newDPoints[i] = new DPoint(points[2 * i], points[2 * i + 1]);
            }
            this._delaunator = new Delaunator(newDPoints);
        } else {
            this.collinear = null;
        }
        this.halfedges = this._delaunator.halfedges;
        int[] halfedges = this.halfedges;
        this.hull = this._delaunator.hull;
        int[] hull = this.hull;
        this.triangles = this._delaunator.triangles;
        int[] triangles = this.triangles;
        for (int i = 0; i < this.indeges.length; i++) {
            this.indeges[i] = -1;
        }
        int[] indeges = this.indeges;
        for (int i = 0; i < this._hullIndex.length; i++) {
            this._hullIndex[i] = -1;
        }
        int[] hullIndex = this._hullIndex;
        int n = halfedges.length;
        for (int e = 0; e < n; ++e) {
            int p = triangles[e % 3 == 2 ? e - 2 : e + 1];
            if (halfedges[e] == -1 || indeges[p] == -1)
                indeges[p] = e;
        }
        for (int i = 0; i < hull.length; i++) {
            hullIndex[hull[i]] = i;
        }
        if (hull.length <= 2 && hull.length > 0) {
            this.triangles = new int[] { -1, -1, -1 };
            this.halfedges = new int[] { -1, -1, -1 };
            this.triangles[0] = hull[0];
            indeges[hull[0]] = 1;
            if (hull.length == 2) {
                indeges[hull[1]] = 0;
                this.triangles[1] = hull[1];
                this.triangles[2] = hull[2];
            }
        }
    }

    /**
     * Returns the corresponding Voronoi instant.The diagram will be clipped to the
     * specified bounds = [ 0, 0, 960, 500 ]
     * 
     * @return a Voronoi object
     */

    public Voronoi voronoi() {
        return new Voronoi(this, new double[] { 0, 0, 960, 500 });
    }

    /**
     * Returns the corresponding Voronoi instant.The diagram will be clipped to the
     * specified bounds = [ xmin, ymin, xmax, ymax ].
     * 
     * @param bounds the bounds of the Voronoi graph
     * @return a Voronoi object
     */

    public Voronoi voronoi(double[] bounds) {
        return new Voronoi(this, bounds);
    }

    /**
     * Returns the corresponding Voronoi instant.The diagram will be clipped to the
     * specified bounds = [ xmin, ymin, xmax, ymax ]. Float version
     * 
     * @param bounds the bounds of the Voronoi graph
     * @return a Voronoi object
     */

    public Voronoi voronoiFloat(float[] bounds) {
        double[] da = new double[bounds.length];
        for (int i = 0; i < da.length; i++) {
            da[i] = (double) bounds[i];
        }
        return this.voronoi(da);
    }

    /**
     * Returns an int array of the indexes of the neighboring points to the
     * specified point i.
     * 
     * @param i the index of the point
     * @return an int array of the indexes of the neighboring points
     */

    public int[] neighbors(int i) {
        int[] indeges = this.indeges;
        int[] hull = this.hull;
        int[] hullIndex = this._hullIndex;
        int[] halfedges = this.halfedges;
        int[] triangles = this.triangles;
        int[] collinear = this.collinear;
        List<Integer> result = new ArrayList<Integer>();

        if (collinear != null) {
            int l = indexOfIntArray(collinear, i);
            if (l > 0) {
                result.add(collinear[l - 1]);
                System.out.println("Delaunay.neighbors() BP1");
            } else if (l < collinear.length - 1) {
                result.add(collinear[l + 1]);
                System.out.println("Delaunay.neighbors() BP2");
            } else {
                return result.stream().mapToInt(ii -> ii).toArray();
            }
        }

        int e0 = indeges[i];
        if (e0 == -1)
            return result.stream().mapToInt(ii -> ii).toArray();
        int e = e0;
        int p0 = -1;
        do {
            p0 = triangles[e];
            result.add(p0);
            e = e % 3 == 2 ? e - 2 : e + 1;
            if (triangles[e] != i)
                return result.stream().mapToInt(ii -> ii).toArray();
            e = halfedges[e];
            if (e == -1) {
                int p = hull[(_hullIndex[i] + 1 + hull.length) % hull.length];
                if (p != p0) {
                    result.add(p);
                }
                return result.stream().mapToInt(ii -> ii).toArray();
            }
        } while (e != e0);
        return result.stream().mapToInt(ii -> ii).toArray();
    }

    /**
     * Returns the index of the input point that is closest to the specified point
     * (x, y). The search is started at the point 0.
     * 
     * @param x point.x
     * @param y point.y
     * @return the found point index
     */

    public int find(double x, double y) {
        return this.find(x, y, 0);
    }

    /**
     * Returns the index of the input point that is closest to the specified point
     * (x, y). The search is started at the specified point i.
     * 
     * @param x point.x
     * @param y point.y
     * @param i the index to start search
     * @return the found point index
     */

    public int find(double x, double y, int i) {
        int i0 = i;
        int c;
        while ((c = this._step(i, x, y)) >= 0 && c != i && c != i0)
            i = c;
        return c;
    }

    /**
     * Returns the index of the input point that is closest to the specified point
     * (x, y). The search is started at the point 0. Float version
     * 
     * @param x point.x
     * @param y point.y
     * @return the found point index
     */

    public int findFloat(float x, float y) {
        return this.findFloat(x, y, 0);
    }

    /**
     * Returns the index of the input point that is closest to the specified point
     * (x, y). The search is started at the specified point i. Float version
     * 
     * @param x point.x
     * @param y point.y
     * @param i the index to start search
     * @return the found point index
     */

    public int findFloat(float x, float y, int i) {
        return this.find((double) x, (double) y, i);
    }

    /**
     * @param i
     * @param x
     * @param y
     * @return int
     */
    public int _step(int i, double x, double y) {
        int[] indeges = this.indeges;
        int[] hull = this.hull;
        int[] _hullIndex = this._hullIndex;
        int[] halfedges = this.halfedges;
        int[] triangles = this.triangles;
        double points[] = this.points;
        if (indeges[i] == -1 || points.length == 0)
            return (i + 1) % (points.length >> 1);
        int c = i;
        double dc = Math.pow(x - points[i * 2], 2) + Math.pow(y - points[i * 2 + 1], 2);
        int e0 = indeges[i];
        int e = e0;
        do {
            int t = triangles[e];
            double dt = Math.pow(x - points[t * 2], 2) + Math.pow(y - points[t * 2 + 1], 2);
            if (dt < dc) {
                dc = dt;
                c = t;
            }
            e = e % 3 == 2 ? e - 2 : e + 1;
            if (triangles[e] != i)
                break;
            e = halfedges[e];
            if (e == -1) {
                e = hull[(_hullIndex[i] + 1) % hull.length];
                if (e != t) {
                    if (Math.pow(x - points[e * 2], 2) + Math.pow(y - points[e * 2 + 1], 2) < dc)
                        return e;
                }
                break;
            }
        } while (e != e0);
        return c;
    }

    /*
     * public void render(){
     * 
     * }
     */

    /*
     * public void renderPoint(){
     * 
     * }
     */

    public void renderHull(Polygon context) {
        int[] hull = this.hull;
        double[] points = this.points;
        int h = hull[0] * 2;
        int n = hull.length;
        context.moveTo(points[h], points[h + 1]);
        for (int i = 1; i < n; i++) {
            int hh = 2 * hull[i];
            context.moveTo(points[hh], points[hh + 1]);
        }
        context.closePath();
    }

    /**
     * Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], ... , [ x0, y0 ]]
     * representing the convex hull.
     * 
     * @return an array representing the convex hull
     */

    public double[][] hullPolygon() {
        Polygon polygon = new Polygon();
        this.renderHull(polygon);
        return polygon.value();
    }

    /**
     * Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], ... , [ x0, y0 ]]
     * representing the convex hull. Float version
     * 
     * @return an array representing the convex hull
     */

    public float[][] hullPolygonFloat() {
        Polygon polygon = new Polygon();
        this.renderHull(polygon);
        double[][] da = polygon.value();
        float[][] fa = new float[da.length][2];
        for (int i = 0; i < da.length; i++) {
            fa[i] = new float[] { (float) da[i][0], (float) da[i][1] };
        }
        return fa;
    }

    /**
     * @param i
     * @param context
     */
    public void renderTriangle(int i, Polygon context) {
        int[] triangles = this.triangles;
        double[] points = this.points;
        int t0 = triangles[i *= 3] * 2;
        int t1 = triangles[i + 1] * 2;
        int t2 = triangles[i + 2] * 2;
        context.moveTo(points[t0], points[t0 + 1]);
        context.lineTo(points[t1], points[t1 + 1]);
        context.lineTo(points[t2], points[t2 + 1]);
        context.closePath();
    }

    /**
     * Returns an array of the closed polygon [ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ],
     * [ x0, y0 ] ]
     * representing the each triangle.
     * 
     * @return an array containing all triangles
     */

    public double[][][] trianglePolygons() {
        List<double[][]> res = new ArrayList<double[][]>();
        int[] triangles = this.triangles;
        for (int i = 0; i < triangles.length / 3; i++) {
            res.add(this.trianglePolygon(i));
        }
        return res.toArray(new double[0][0][2]);
    }

    /**
     * Returns an array of the closed polygon [ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ],
     * [ x0, y0 ] ]
     * representing the each triangle. Float version
     * 
     * @return an array containing all triangles
     */

    public float[][][] trianglePolygonsFloat() {
        double[][][] da = this.trianglePolygons();
        float[][][] fa = new float[da.length][4][2];
        for (int i = 0; i < da.length; i++) {
            double[][] dtri = da[i];
            for (int j = 0; j < dtri.length; j++) {
                double[] dp = dtri[j];
                fa[i][j][0] = (float) dp[0];
                fa[i][j][1] = (float) dp[1];
            }
        }
        return fa;
    }

    /**
     * Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ], [ x0, y0 ] ]
     * representing the triangle i.
     * 
     * @param i the index of the triangle
     * @return an array representing the triangle
     */

    public double[][] trianglePolygon(int i) {
        Polygon polygon = new Polygon();
        this.renderTriangle(i, polygon);
        return polygon.value();
    }

    /**
     * Returns the closed polygon [ [ x0, y0 ], [ x1, y1 ], [ x2, y2 ], [ x0, y0 ] ]
     * representing the triangle i. Float version
     * 
     * @param i the index of the triangle
     * @return an array representing the triangle
     */

    public float[][] trianglePolygonFloat(int i) {
        Polygon polygon = new Polygon();
        this.renderTriangle(i, polygon);
        double[][] da = polygon.value();
        float[][] fa = new float[da.length][2];
        for (int j = 0; j < fa.length; j++) {
            fa[i] = new float[] { (float) da[i][0], (float) da[i][1] };
        }
        return fa;
    }

    public static int indexOfIntArray(int[] array, int key) {
        int returnvalue = -1;
        for (int i = 0; i < array.length; ++i) {
            if (key == array[i]) {
                returnvalue = i;
                break;
            }
        }
        return returnvalue;
    }
}
