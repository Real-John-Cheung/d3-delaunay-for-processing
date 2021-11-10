
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


public class Delaunay {
    public static double tau = (double) (2 * Math.PI);
    public Delaunator _delaunator;
    public int[] indeges;
    public int[] _hullIndex;
    public double[] points;
    public int[] collinear;
    public int[] halfedges;
    public int[] hull;
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

    static public Delaunay from(double[][] points) {
        double[] flatArr = new double[points.length * 2];
        for (int i = 0; i < points.length; i ++) {
            flatArr[i * 2] = points[i][0];
            flatArr[i * 2 + 1] = points[i][1];
        }
        return new Delaunay(flatArr);
    }

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

    /*
    public Delaunay update() {
        //not avaliable now
        return this;
    }
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

    public Voronoi voronoi() {
        return new Voronoi(this, new double[] { 0, 0, 960, 500 });
    }

    public Voronoi voronoi(double[] bounds) {
        return new Voronoi(this, bounds);
    }

    public int[] neighbors(int i) {
        int[] indeges = this.indeges;
        int[] hull = this.hull;
        int[] hullIndex = this._hullIndex;
        int[] halfedges = this.halfedges;
        int[] triangles = this.triangles;
        int[] collinear = this.collinear;
        List<Integer> result = new ArrayList<Integer>();

        if (collinear != null) {
            int l = Arrays.asList(collinear).indexOf(i);
            if (l > 0){
                result.add(collinear[l - 1]);
            } else if (l < collinear.length - 1) {
                result.add(collinear[l + 1]);
            } else {
                return result.stream().mapToInt(ii -> ii).toArray();
            }
        }
        
        int e0 = indeges[i];
        if (e0 == -1) return result.stream().mapToInt(ii -> ii).toArray();
        int e = e0;
        int p0 = -1;
        do {
            p0 = triangles[e];
            result.add(p0);
            e = e % 3 == 2 ? e - 2 : e + 1;
            if (triangles[e] != i) return result.stream().mapToInt(ii -> ii).toArray();
            e = halfedges[e];
            if (e == -1) {
                int p = hull[(_hullIndex[i] + 1) % hull.length];
                if (p != p0) {
                    result.add(p);
                } else {
                    return result.stream().mapToInt(ii -> ii).toArray();
                }
            }
        } while (e != e0);
        return result.stream().mapToInt(ii -> ii).toArray();
    }
    
    public int find(int x, int y) {
        return this.find(x, y, 0);
    }

    public int find(int x, int y, int i) {
        int i0 = i;
        int c;
        while ((c = this._step(i, x, y)) >= 0 && c != i && c != i0)
            i = c;
        return c;
    }

    public int _step(int i, double x, double y) {
        int[] indeges = this.indeges;
        int[] hull = this.hull;
        int[] _hullIndex = this._hullIndex;
        int[] halfedges = this.halfedges;
        int[] triangles = this.triangles;
        double points[] = this.points;
        if (indeges[i] == -1 || points.length == 0) return (i + 1) % (points.length >> 1);
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
            if (triangles[e] != i)  break;
            e = halfedges[e];
            if (e == -1) {
                e = hull[(_hullIndex[i] + 1) % hull.length];
                if (e != t) {
                    if (Math.pow(x - points[e * 2], 2) + Math.pow(y - points[e * 2 + 1], 2) < dc) return e;
                }
                break;
            }
        } while (e != e0);
        return c;
    }

    /*
    public void render(){
    
    }
    */

    /*
    public void renderPoint(){
    
    }
    */

    
    public void renderHull(Polygon context){
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

    public double[][] hullPolygon() {
        Polygon polygon = new Polygon();
        this.renderHull(polygon);
        return polygon.value();
    }

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
    
    public double[][][] trianglePolygons() {
        List<double[][]> res = new ArrayList<double[][]>();
        int[] triangles = this.triangles;
        for (int i = 0; i < triangles.length / 3; i++) {
            res.add(this.trianglePolygon(i));
        }
        return res.toArray(new double[0][0][2]);
    }
    
    public double[][] trianglePolygon(int i) {
        Polygon polygon = new Polygon();
        this.renderTriangle(i, polygon);
        return polygon.value();
    }
}

