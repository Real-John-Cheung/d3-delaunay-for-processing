package d3delaunayforprocessing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.crypto.interfaces.PBEKey;

/**
 * d3-delaunay-for-processing
 * To port the d3-delaunay library to Java using the delaunator-java library by
 * waveware4ai
 * This is the Polygon class
 * 
 * @author John C
 * @version 1.1
 * @since 2022-12-22
 */

public class Polygon {
    public List<double[]> dash;

    public Polygon() {
        this.dash = new ArrayList<double[]>();
    }

    public void moveTo(double x, double y) {
        double[] point = new double[] { x, y };
        this.dash.add(point);
    }
    
    public void closePath() {
        double[] firstPoint = this.dash.get(0).clone();
        this.dash.add(firstPoint);
    }

    public void lineTo(double x, double y) {
        double[] point = new double[] { x, y };
        this.dash.add(point);
    }

    public double[][] value() {
        if (this.dash.size() < 1) return null;
        return this.dash.toArray(new double[0][2]);
    }
}