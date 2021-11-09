package d3delaunayforprocessing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.crypto.interfaces.PBEKey;

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