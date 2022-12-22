
package d3delaunayforprocessing;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import javax.swing.text.html.StyleSheet;

import java.util.Comparator;
import org.waveware.delaunator.Delaunator;
import org.waveware.delaunator.DTriangle;
import org.waveware.delaunator.DPoint;
import org.checkerframework.checker.signedness.qual.Unsigned;
import org.waveware.delaunator.DEdge;
import d3delaunayforprocessing.Polygon;

/**
 * d3-delaunay-for-processing
 * To port the d3-delaunay library to Java using the delaunator-java library by
 * waveware4ai
 * This is the Voronoi class
 * 
 * @author John C
 * @version 1.1
 * @since 2022-12-22
 */

public class Voronoi {
    /**
     * Reference Reference back to the corresponding delaunay instant.
     */
    public Delaunay delaunay;
    private double[] _circumcenters;
    /**
     * The circumcenters of the Delaunay triangles as a double array [ cx0, cy0,
     * cx1, cy1, ... ]
     */
    public double[] circumcenters;
    /**
     * A double array [ vx0, vy0, wx0, wy0, ... ] where each non-zero quadruple
     * describes an open (infinite) cell on the outer hull, giving the directions of
     * two open half-lines.
     */
    public double[] vectors;
    public double xmin;
    public double xmax;
    public double ymin;
    public double ymax;

    /**
     * Voronoi class constructor.
     * 
     * @param delaunay A Delaunay object
     * @param bounds   the boundary for the Voronoi graph, [xmin, ymin, xmax, ymax]
     * @see Delaunay
     */
    public Voronoi(Delaunay delaunay, double[] bounds) {
        double xmin = bounds[0];
        double ymin = bounds[1];
        double xmax = bounds[2];
        double ymax = bounds[3];
        if (xmin > xmax || ymin > ymax) {
            throw new Error("invaild bounds");
        }
        this.delaunay = delaunay;
        this._circumcenters = new double[delaunay.points.length * 2];
        this.vectors = new double[delaunay.points.length * 2];
        this.xmin = xmin;
        this.xmax = xmax;
        this.ymin = ymin;
        this.ymax = ymax;
        this._init();
    }

    /**
     * Voronoi class constructor, but with float.
     * 
     * @param delaunay A Delaunay object
     * @param bounds   the boundary for the Voronoi graph, [xmin, ymin, xmax, ymax]
     * @see Delaunay
     */

    public Voronoi(Delaunay delaunay, float[] bounds) {
        double[] da = new double[bounds.length];
        for (int i = 0; i < da.length; i++) {
            da[i] = (double) bounds[i];
        }
        double xmin = da[0];
        double ymin = da[1];
        double xmax = da[2];
        double ymax = da[3];
        if (xmin > xmax || ymin > ymax) {
            throw new Error("invaild bounds");
        }
        this.delaunay = delaunay;
        this._circumcenters = new double[delaunay.points.length * 2];
        this.vectors = new double[delaunay.points.length * 2];
        this.xmin = xmin;
        this.xmax = xmax;
        this.ymin = ymin;
        this.ymax = ymax;
        this._init();
    }

    /*
     * public Voronoi update(){
     * 
     * }
     */

    public void _init() {
        double[] points = this.delaunay.points;
        int[] hull = this.delaunay.hull;
        int[] triangles = this.delaunay.triangles;
        double[] vectors = this.vectors;

        // draw circles
        this.circumcenters = Arrays.copyOfRange(this._circumcenters, 0, triangles.length / 3 * 2);
        double[] circumcenters = this.circumcenters;
        int j = 0;
        double x;
        double y;
        for (int i = 0; i < triangles.length; i += 3) {
            int t1 = triangles[i] * 2;
            int t2 = triangles[i + 1] * 2;
            int t3 = triangles[i + 2] * 2;
            double x1 = points[t1];
            double y1 = points[t1 + 1];
            double x2 = points[t2];
            double y2 = points[t2 + 1];
            double x3 = points[t3];
            double y3 = points[t3 + 1];

            double dx = x2 - x1;
            double dy = y2 - y1;
            double ex = x3 - x1;
            double ey = y3 - y1;
            double ab = (dx * ey - dy * ex) * 2;
            if (Math.abs(ab) < 1e-9) {
                double a = 1e9;
                int r = triangles[0] * 2;
                a *= ((points[r] - x1) * ey - (points[r + 1] - y1) * ex) > 0 ? 1
                        : (((points[r] - x1) * ey - (points[r + 1] - y1) * ex) == 0 ? 0 : -1);
                x = (x1 + x3) / 2 - a * ey;
                y = (y1 + y3) / 2 + a * ex;
            } else {
                double d = 1 / ab;
                double bl = dx * dx + dy * dy;
                double cl = ex * ex + ey * ey;
                x = x1 + (ey * bl - dy * cl) * d;
                y = y1 + (dx * cl - ex * bl) * d;
            }
            circumcenters[j] = x;
            circumcenters[j + 1] = y;
            // at the end change j
            j += 2;
        }

        // rays
        int h = hull[hull.length - 1];
        int p0;
        int p1 = h * 4;
        double x0;
        double x1 = points[2 * h];
        double y0;
        double y1 = points[2 * h + 1];
        for (int i = 0; i < vectors.length; i++) {
            vectors[i] = 0;
        }
        for (int i = 0; i < hull.length; i++) {
            h = hull[i];
            p0 = p1;
            x0 = x1;
            y0 = y1;
            p1 = h * 4;
            x1 = points[2 * h];
            y1 = points[2 * h + 1];
            vectors[p0 + 2] = vectors[p1] = y0 - y1;
            vectors[p0 + 3] = vectors[p1 + 1] = x1 - x0;
        }
    }

    /*
     * public void render(Polygon context){
     * 
     * }
     */

    /*
     * public void renderBounds(Polygon context){
     * 
     * }
     */
    public void renderCell(int i, Polygon context) {
        double[] points = this._clip(i);
        if (points == null || points.length < 2)
            return;
        context.moveTo(points[0], points[1]);
        int n = points.length;
        while (points[0] == points[n - 2] && points[1] == points[n - 1] && n > 1)
            n -= 2;
        for (int j = 2; j < n; j += 2) {
            if (points[j] != points[j - 2] || points[j + 1] != points[j - 1])
                context.lineTo(points[j], points[j + 1]);
        }
        context.closePath();
    }

    /*
     * public int[] cellPolygons(){
     * 
     * }
     */

    /**
     * Returns the convex, closed polygon [ [ x0, y0 ], [ x1, y1 ], ..., [ x0, y0 ]
     * ] representing the cell for the specified point i.
     * 
     * @param i the index of the cell
     * @return the points in an array [[x0, y0],[x1, y1], ...]
     */

    public double[][] cellPolygon(int i) {
        Polygon polygon = new Polygon();
        this.renderCell(i, polygon);
        return polygon.value();
    }

    /**
     * Returns the convex, closed polygon [ [ x0, y0 ], [ x1, y1 ], ..., [ x0, y0 ]
     * ] representing the cell for the specified point i. Return in float.
     * 
     * @param i the index of the cell
     * @return the points in an array [[x0, y0],[x1, y1], ...]
     */

    public float[][] cellPolygonFloat(int i) {
        Polygon polygon = new Polygon();
        this.renderCell(i, polygon);
        double[][] da = polygon.value();
        float[][] re = new float[da.length][2];
        for (int j = 0; j < re.length; j++) {
            re[i] = new float[] { (float) da[i][0], (float) da[i][1] };
        }
        return re;
    }

    /*
     * public void _renderSegment(){
     * 
     * }
     */

    /**
     * Returns true if the cell with the specified index i contains the specified
     * point (x, y).
     * 
     * @param i the index of the cell
     * @param x x coordinate value of the point
     * @param y y coordinate value of the point
     * @return true or false
     */

    public boolean contains(int i, double x, double y) {
        return this.delaunay._step(i, x, y) == i;
    }

    /**
     * Returns true if the cell with the specified index i contains the specified
     * point (x, y).
     * 
     * @param i the index of the cell
     * @param x x coordinate value of the point
     * @param y y coordinate value of the point
     * @return true or false
     */

    public boolean contains(int i, float x, float y) {
        return this.contains(i, (double) x, (double) y);
    }

    /**
     * Returns an int array of the indexes of the cells that share a common edge
     * with the specified cell i.
     * 
     * @param i the index of the cell
     * @return an int array of the indexes of the cell's neighbors
     */

    public int[] neighbors(int i) {
        List<Integer> result = new ArrayList<Integer>();
        double[] ci = this._clip(i);
        if (ci != null) {
            for (int ind = 0; ind < this.delaunay.neighbors(i).length; ind++) {
                int j = this.delaunay.neighbors(i)[ind];
                double[] cj = this._clip(j);
                if (cj != null) {
                    int li = ci.length;
                    int lj = cj.length;
                    for (int ai = 0; ai < ci.length; ai += 2) {
                        for (int aj = 0; aj < cj.length; aj += 2) {
                            if (ci[ai] == cj[aj]
                                    && ci[ai + 1] == cj[aj + 1]
                                    && ci[(ai + 2) % li] == cj[(aj + lj - 2) % lj]
                                    && ci[(ai + 3) % li] == cj[(aj + lj - 1) % lj]) {
                                result.add(j);
                            }
                        }
                    }
                }
            }
        }
        return result.stream().mapToInt(ii -> ii).toArray();
    }

    public double[] _cell(int i) {
        double[] circumcenters = this.circumcenters;
        int[] indeges = this.delaunay.indeges;
        int[] halfedges = this.delaunay.halfedges;
        int[] triangles = this.delaunay.triangles;
        int e0 = indeges[i];
        if (e0 == -1)
            return null;
        List<Double> points = new ArrayList<Double>();
        int e = e0;
        do {
            int t = (int) Math.floor(e / 3);
            points.add(circumcenters[t * 2]);
            points.add(circumcenters[t * 2 + 1]);
            e = e % 3 == 2 ? e - 2 : e + 1;
            if (triangles[e] != i)
                break;
            e = halfedges[e];
        } while (e != e0 && e != -1);
        return points.stream().mapToDouble(d -> d).toArray();
    }

    public double[] _clip(int i) {
        if (i == 0 && this.delaunay.hull.length == 1) {
            return new double[] { this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin,
                    this.ymin };
        }
        double[] points = this._cell(i);
        if (points == null || points.length < 1)
            return null;
        double[] V = this.vectors;
        int v = i * 4;
        if (((v > 0 && v < V.length) && V[v] != 0) || ((v + 1 > 0 && v + 1 < V.length) && V[v + 1] != 0)) {
            return this._clipInfinite(i, points, V[v], V[v + 1], V[v + 2], V[v + 3]);
        }
        return this._clipFinite(i, points);
    }

    public double[] _clipFinite(int i, double[] points) {
        int n = points.length;
        List<Double> P = null;
        double x0;
        double y0;
        double x1 = points[n - 2];
        double y1 = points[n - 1];
        int c0;
        int c1 = this._regioncode(x1, y1);
        int e0;
        int e1 = 0;
        for (int j = 0; j < n; j += 2) {
            x0 = x1;
            y0 = y1;
            x1 = points[j];
            y1 = points[j + 1];
            c0 = c1;
            c1 = this._regioncode(x1, y1);
            if (c0 == 0 && c1 == 0) {
                e0 = e1;
                e1 = 0;
                if (P != null && P.size() > 0) {
                    P.add(x1);
                    P.add(y1);
                } else {
                    P = new ArrayList<Double>();
                    P.add(x1);
                    P.add(y1);
                }
            } else {
                double[] S;
                double sx0;
                double sy0;
                double sx1;
                double sy1;
                if (c0 == 0) {
                    S = this._clipSegment(x0, y0, x1, y1, c0, c1);
                    if (S == null)
                        continue;
                    sx0 = S[0];
                    sy0 = S[1];
                    sx1 = S[2];
                    sy1 = S[3];
                } else {
                    S = this._clipSegment(x1, y1, x0, y0, c1, c0);
                    if (S == null)
                        continue;
                    sx1 = S[0];
                    sy1 = S[1];
                    sx0 = S[2];
                    sy0 = S[3];
                    e0 = e1;
                    e1 = this._edgecode(sx0, sy0);
                    if (e0 != 0 && e1 != 0)
                        this._edge(i, e0, e1, P, P.size());
                    if (P != null && P.size() > 0) {
                        P.add(sx0);
                        P.add(sy0);
                    } else {
                        P = new ArrayList<Double>();
                        P.add(sx0);
                        P.add(sy0);
                    }
                }
                e0 = e1;
                e1 = this._edgecode(sx1, sy1);
                if (e0 != 0 && e1 != 0)
                    this._edge(i, e0, e1, P, P.size());
                if (P != null && P.size() > 0) {
                    P.add(sx1);
                    P.add(sy1);
                } else {
                    P = new ArrayList<Double>();
                    P.add(sx1);
                    P.add(sy1);
                }
            }
        }
        if (P != null && P.size() > 1) {
            e0 = e1;
            e1 = this._edgecode(P.get(0), P.get(1));
            if (e0 != 0 && e1 != 0)
                this._edge(i, e0, e1, P, P.size());
        } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
            return new double[] { this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin,
                    this.ymin };
        }
        return P.stream().mapToDouble(d -> d).toArray();
    }

    public double[] _clipSegment(double x0, double y0, double x1, double y1, int c0, int c1) {
        while (true) {
            if (c0 == 0 && c1 == 0)
                return new double[] { x0, y0, x1, y1 };
            if ((c0 & c1) != 0)
                return null;
            double x;
            double y;
            int c = c0 == 0 ? c1 : c0;
            if ((c & 0b1000) != 0) {
                x = x0 + (x1 - x0) * (this.ymax - y0) / (y1 - y0);
                y = this.ymax;
            } else if ((c & 0b0100) != 0) {
                x = x0 + (x1 - x0) * (this.ymin - y0) / (y1 - y0);
                y = this.ymin;
            } else if ((c & 0b0010) != 0) {
                y = y0 + (y1 - y0) * (this.xmax - x0) / (x1 - x0);
                x = this.xmax;
            } else {
                y = y0 + (y1 - y0) * (this.xmin - x0) / (x1 - x0);
                x = this.xmin;
            }
            if (c0 != 0) {
                x0 = x;
                y0 = y;
                c0 = this._regioncode(x0, y0);
            } else {
                x1 = x;
                y1 = y;
                c1 = this._regioncode(x1, y1);
            }
        }
    }

    public double[] _clipInfinite(int i, double[] points, double vx0, double vy0, double vxn, double vyn) {
        List<Double> P = DoubleStream.of(points.clone()).boxed().collect(Collectors.toList());
        double[] p = null;
        p = this._project(P.get(0), P.get(1), vx0, vy0);
        if (p != null && p.length > 1) {
            P.add(0, p[1]);
            P.add(0, p[0]);
        }
        p = this._project(P.get(P.size() - 2), P.get(P.size() - 1), vxn, vyn);
        if (p != null && p.length > 1) {
            P.add(p[0]);
            P.add(p[1]);
        }
        P = DoubleStream.of(this._clipFinite(i, P.stream().mapToDouble(d -> d).toArray())).boxed()
                .collect(Collectors.toList());
        if (P != null && P.size() > 1) {
            int n = P.size();
            int c0;
            int c1 = this._edgecode(P.get(n - 2), P.get(n - 1));
            for (int j = 0; j < n; j += 2) {
                c0 = c1;
                c1 = this._edgecode(P.get(j), P.get(j + 1));
                if (c0 != 0 && c1 != 0) {
                    j = this._edge(i, c0, c1, P, j);
                    n = P.size();
                }
            }
        } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
            return new double[] { this.xmin, this.ymin, this.xmax, this.ymin, this.xmax, this.ymax, this.xmin,
                    this.ymax };
        }
        return P.stream().mapToDouble(d -> d).toArray();
    }

    public int _edge(int i, int e0, int e1, List<Double> P, int j) {
        while (e0 != e1) {
            double x;
            double y;
            switch (e0) {
                case 0b0101:
                    e0 = 0b0100;
                    continue;
                case 0b0100:
                    e0 = 0b0110;
                    x = this.xmax;
                    y = this.ymin;
                    break;
                case 0b0110:
                    e0 = 0b0010;
                    continue;
                case 0b0010:
                    e0 = 0b1010;
                    x = this.xmax;
                    y = this.ymax;
                    break;
                case 0b1010:
                    e0 = 0b1000;
                    continue;
                case 0b1000:
                    e0 = 0b1001;
                    x = this.xmin;
                    y = this.ymax;
                    break;
                case 0b1001:
                    e0 = 0b0001;
                    continue;
                case 0b0001:
                    e0 = 0b0101;
                    x = this.xmin;
                    y = this.ymin;
                    break;
                default:
                    x = 0;
                    y = 0;
                    break;
            }
            // if P[j] or P[j+1] are undefined, the conditional statement should be
            // executed.
            if ((j < 0 || j >= P.size()) || (j + 1 < 0 || j + 1 >= P.size())) {
                if (this.contains(i, x, y)) {
                    P.add(j, y);
                    P.add(j, x);
                    j += 2;
                }
            } else if ((P.get(j) != x || P.get(j + 1) != y) && this.contains(i, x, y)) {
                P.add(j, y);
                P.add(j, x);
                j += 2;
            }
        }
        if (P.size() > 4) {
            for (int j2 = 0; j2 < P.size(); j2 += 2) {
                int j3 = (i + 2) % P.size();
                int k = (i + 4) % P.size();
                if ((j2 < 0 || j2 > P.size() - 1) || (j3 < 0 || j3 > P.size() - 1) || (k < 0 || k > P.size() - 1))
                    continue;
                if ((j2 + 1 < 0 || j2 + 1 > P.size() - 1) || (j3 + 1 < 0 || j3 + 1 > P.size() - 1)
                        || (k + 1 < 0 || k + 1 > P.size() - 1))
                    continue;
                if ((P.get(j2) == P.get(j3) && P.get(j3) == P.get(k))
                        || (P.get(j2 + 1) == P.get(j3 + 1) && P.get(j3 + 1) == P.get(k + 1))) {
                    P.remove(j3);
                    P.remove(j3);
                    j2 -= 2;
                }
            }
        }
        return j;
    }

    public double[] _project(double x0, double y0, double vx, double vy) {
        double t = Double.POSITIVE_INFINITY;
        double c;
        double x = 0;
        double y = 0;
        if (vy < 0) {
            if (y0 <= this.ymin)
                return null;
            c = (this.ymin - y0) / vy;
            if (c < t) {
                y = this.ymin;
                x = x0 + (t = c) * vx;
            }
        } else if (vy > 0) {
            if (y0 >= this.ymax)
                return null;
            c = (this.ymax - y0) / vy;
            if (c < t) {
                y = this.ymax;
                x = x0 + (t = c) * vx;
            }
        }
        if (vx > 0) {
            if (x0 >= this.xmax)
                return null;
            c = (this.xmax - x0) / vx;
            if (c < t) {
                x = this.xmax;
                y = y0 + (t = c) * vy;
            }
        } else if (vx < 0) {
            if (x0 <= this.xmin)
                return null;
            c = (this.xmin - x0) / vx;
            if (c < t) {
                x = this.xmin;
                y = y0 + (t = c) * vy;
            }
        }
        return new double[] { x, y };
    }

    public int _edgecode(double x, double y) {
        return (x == this.xmin ? 0b0001 : x == this.xmax ? 0b0010 : 0b0000)
                | (y == this.ymin ? 0b0100 : y == this.ymax ? 0b1000 : 0b0000);
    }

    public int _regioncode(double x, double y) {
        return (x < this.xmin ? 0b0001 : x > this.xmax ? 0b0010 : 0b0000)
                | (y < this.ymin ? 0b0100 : y > this.ymax ? 0b1000 : 0b0000);
    }
}
