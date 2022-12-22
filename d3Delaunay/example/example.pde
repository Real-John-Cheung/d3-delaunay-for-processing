import d3delaunayforprocessing.*;

double[] points = new double[40];


void setup() {
  size(600,600);
  for (int i = 0; i < points.length; i ++){
    points[i] = random(10,590);
  }
  Delaunay d = new Delaunay(points);
  background(255);
  noFill();
  stroke(0);
  strokeWeight(3);
  for (int i = 0; i < points.length; i +=2){
    point((float)points[i], (float)points[i+1]);
  }
  strokeWeight(0.5);
  Voronoi v = d.voronoi(new double[] {0,0,width,height});
  for (int i = 0; i < points.length/2; i ++){
    double[][] cell = v.cellPolygon(i);
    if (cell == null) continue;
    beginShape();
    for (int j = 0; j < cell.length; j ++){
      vertex((float)cell[j][0],(float)cell[j][1]);
    }
    endShape();
  }
  stroke(255,0,0);
  double[][][] triangles = d.trianglePolygons();
  for(int i = 0; i < triangles.length; i ++) {
    double[][] tri = triangles[i];
    beginShape();
    for (int j = 0; j < tri.length; j ++) {
      double[] p = tri[j];
      vertex((float)p[0], (float)p[1]);
    }
    endShape();
  }
}
