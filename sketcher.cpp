/*
 * (c) 2017 the mullican group
 * kevin mullican
 *
 * sketcher.cpp
 *
 * tool to extract 'nodes' from a collada file exported from sketchup, then
 * write the nodes and beams to a beamng.drive jbeam file
 */

#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

struct vect {
  double x;
  double y;
  double z;
  vect() { x = 0; y = 0; z = 0; }
  vect(double _x, double _y, double _z) {
    x = _x; y = _y; z = _z;
  }
  vect (const vect &v) {
    x = v.x; y = v.y; z = v.z;
  }
  vect &set (double _x, double _y, double _z) {
    x = _x; y = _y; z = _z;
    return *this;
  }
  vect cross(const vect &v) {
    vect c;
    c.x = (v.y * z) - (v.z * y);
    c.y = (v.z * x) - (v.x * z);
    c.z = (v.x * y) - (v.y * x);
    return c;
  }
  double mag_sq() {
    return x*x + y*y + z*z;
  }
  double mag() {
    return sqrt(mag_sq());
  }
  vect norm() {
    double l = mag();
    return vect (x/l, y/l, z/l);
  }
  vect neg() {
    return vect(-x, -y, -z);
  }
  vect& operator+ (const vect &v) {
    x += v.x; y += v.y; z += v.z;
    return *this;
  }
  vect& operator- (const vect &v) {
    x -= v.x; y -= v.y; z -= v.z;
    return *this;
  }
  vect& operator= (const vect &v) {
    x = v.x; y = v.y; z = v.z;
    return *this;
  }
  bool operator== (const vect &v) const {
    if (x == v.x && y == v.y && z == v.z)
      return true;
    return false;
  }
};

struct beam {
  vect p1;
  vect p2;
  beam() {}
  beam(vect _p1, vect _p2) {
    p1 = _p1; p2 = _p2;
  }
  beam(const beam &b) {
    p1 = b.p1; p2 = b.p2;
  }
  beam &set(const vect &_p1, const vect &_p2) {
    p1 = _p1;
    p2 = _p2;
    return *this;
  }
  beam &operator= (const beam &b) {
    p1 = b.p1;
    p2 = b.p2;
    return *this;
  }
  bool operator== (const beam &b) const {
    if (p1 == b.p1 && p2 == b.p2) return true;
    if (p1 == b.p2 && p2 == b.p1) return true;
    return false;
  }
};

struct triangle {
  vect p1;
  vect p2;
  vect p3;
  triangle() {}
  triangle(vect _p1, vect _p2, vect _p3) {
    p1 = _p1;
    p2 = _p2;
    p3 = _p3;
  }
  triangle(const triangle &t) {
    p1 = t.p1;
    p2 = t.p2;
    p3 = t.p3;
  }
  triangle &set (vect _p1, vect _p2, vect _p3) {
    p1 = _p1;
    p2 = _p2;
    p3 = _p3;
    return *this;
  }
  int sharedPoints(const triangle &t) const {
    int same = 0;
    if (p1 == t.p1 || p1 == t.p2 || p1 == t.p3)
      same++;
    if (p2 == t.p1 || p2 == t.p2 || p2 == t.p3)
      same++;
    if (p3 == t.p1 || p3 == t.p2 || p3 == t.p3)
      same++;
    return same;
  }
  vect normal() {
    vect pp2 = p2 - p1;
    vect pp3 = p3 - p1;
    vect cr = pp2 . cross (pp3);
    return cr.norm();
  }
  bool touching(const triangle &t) {
    if (sharedPoints(t) == 1) return true;
    return false;
  }
  bool adjacent(const triangle &t) {
    if (sharedPoints(t) == 2) return true;
    return false;
  }
  bool same(const triangle &t) {
    if (sharedPoints(t) == 3) return true;
    return false;
  }
  bool coplaner(triangle &t) {
    vect tn = t.normal();
    vect n = normal();
    if (tn == n) return true;
    if (tn == n.neg()) return true;
    return false;
  }
  bool contains (const vect &pt) const {
    if (pt == p1 || pt == p2 || pt == p3)
      return true;
    return false;
  }
  triangle& operator= (const triangle &t) {
    p1 = t.p1; p2 = t.p2; p3 = t.p3;
    return *this;
  }
  bool operator== (const triangle &t) const {
    if (sharedPoints(t) == 3)
      return true;
    return false;
  }
};

vect VectFromNodesAt (vector <double> &nodes, const int at) {
  int idx = at * 3;
  vect v (nodes [idx], nodes [idx+1], nodes [idx+2]);
  return v;
}

bool containsBeam(const vector <beam> &beams, const beam &theBeam) {
  for (int i = 0; i < beams . size (); ++i) {
    if (beams[i] == theBeam)
      return true;
  }
  return false;
}

void addUniqueBeam(vector <beam> &beams, const beam &theBeam) {
  if (! containsBeam(beams, theBeam))
    beams . push_back (theBeam);
}

bool oppositePoints(beam &theBeam, const triangle &t1, const triangle &t2) {
  int shared = t1 . sharedPoints (t2);
  if (shared != 2) {
    printf ("shared: %d\n", shared);
    return false;
  }

  vect p1;
  if (! t2.contains(t1.p1))
    p1 = t1.p1;
  else if (! t2.contains(t1.p2))
    p1 = t1.p2;
  else if (! t2.contains(t1.p3))
    p1 = t1.p3;

  vect p2;
  if (! t1.contains(t2.p1))
    p2 = t2.p1;
  else if (! t1.contains(t2.p1))
    p2 = t2.p1;
  else if (! t1.contains(t2.p1))
    p2 = t2.p1;

  theBeam.set(p1, p2);
  return true;
}

vector <unsigned int> UintSplit (string str, const char *delim) {
  vector <unsigned int> tokens;
  int len = str.length ();
  char *cstr = new char[len + 1];
  if (! cstr) return tokens;
  strncpy (cstr, str.c_str(), len);
  cstr[len] = 0;
  char *buf = cstr;
  while (char *got = strtok (buf, delim))
   { tokens . push_back ((unsigned int)strtoul (got, NULL, 0));
     buf = NULL;
   }
  delete [] (cstr);
  return tokens;
}

vector <double> DoubleSplit (string str, const char *delim) {
  vector <double> tokens;
  int len = str.length ();
  char *cstr = new char[len + 1];
  if (! cstr) return tokens;
  strncpy (cstr, str.c_str(), len);
  cstr[len] = 0;
  char *buf = cstr;
  while (char *got = strtok (buf, delim))
   { tokens . push_back (strtod (got, NULL));
     buf = NULL;
   }
  delete [] (cstr);
  return tokens;
}

XMLElement *FindElement (XMLElement *parent,
                         vector <string> hierarchy,
                         const char *attr = NULL,
                         const char *named = NULL) {
  if (! parent) return NULL;
  int elems = hierarchy . size ();
  if (! elems) return NULL;
  XMLElement *it = parent;
  for (int i = 0; i < elems; ++i) {
    XMLElement *buf = it -> FirstChildElement (hierarchy [i] . c_str());
    if (! buf)
      return NULL;
    if (attr && named && strcmp (buf -> Attribute(attr), named))
      return NULL;
    it = buf;
  }
  return it;
}

int main (int argc, char **argv) {
  string fname;
  int acm1 = argc - 1;
  for (int i = 1; i < acm1; ++i) {
    if (! strncmp ("-f", argv[i], 2))
      fname = argv[i+1];
  }

  if (fname.empty()) {
    printf ("please specify a file with '-f <filename>'\n");
    return 1;
  }

  if (access (fname.c_str(), R_OK)) {
    printf ("unable to read: %s\n", fname.c_str());
    return 2;
  }

  XMLDocument doc;
  printf ("loading %s\n", fname.c_str());
  XMLError ok = doc.LoadFile (fname.c_str());
  if (ok != XML_SUCCESS) {
    printf ("unable to parse the xml file: %s\n", fname.c_str());
    return 3;
  }

  XMLElement *collada = doc.FirstChildElement("COLLADA");
  if (! collada) {
    printf ("unable to find the COLLADA XML element\n");
    return 4;
  }

  vector <string> mh;
  mh . push_back ("library_geometries");
  mh . push_back ("geometry");
  mh . push_back ("mesh");
  XMLElement *mesh = FindElement (collada, mh);
  if (! mesh) {
    printf ("unable to find the mesh\n");
    return 5;
  }

  vector <string> fh;
  fh . push_back ("source");
  fh . push_back ("float_array");
  XMLElement *fa = FindElement (mesh, fh);
  if (! fa) {
    printf ("unable to find the float_array XML element\n");
    return 6;
  }

  vector <string> th;
  th . push_back ("triangles");
  XMLElement *tri = FindElement (mesh, th);
  if (! tri) {
    printf ("unable to find the triangles XML element\n");
    return 7;
  }

  vector <string> tvh;
  tvh . push_back ("p");
  XMLElement *tri_vert = FindElement (tri, tvh);
  if (! tri_vert) {
    printf ("unable to find the triangle vertex XML element\n");
    return 8;
  }

  int want = fa -> IntAttribute ("count");
  string node_text = fa -> GetText ();
  vector <double> nodes = DoubleSplit (node_text, " ");
  int node_points = (int)nodes . size ();
  int node_count = (int)node_points / 3;
  if (node_points != want)
    printf ("node want count %d not equal got count %d\n", want, node_points);
  printf ("found %d node points\n", node_points);

  want = tri -> IntAttribute ("count");
  string tri_text = tri_vert -> GetText ();
  vector <unsigned int> tris = UintSplit (tri_text, " ");
  int tri_points = (int)tris . size ();
  int tri_count = (int)tri_points / 3;
  if (tri_count != want)
    printf ("triangle want count %d not equal got count %d\n", want, tri_count);
  printf ("found %d triangles\n", tri_count);

  vector <triangle> triangles;
  for (int i = 0; i < tri_count; ++i) {
    int idx = i * 3;
    unsigned int pidx1 = tris [idx];
    unsigned int pidx2 = tris [idx+1];
    unsigned int pidx3 = tris [idx+2];

    if (pidx1*3 >= node_points ||
        pidx2*3 >= node_points ||
        pidx3*3 >= node_points) {
      printf ("triangle node reference out of range: %u, %u, %u\n",
               pidx1, pidx2, pidx3);
      continue;
    }

    triangle t (VectFromNodesAt (nodes, pidx1),
                VectFromNodesAt (nodes, pidx2),
                VectFromNodesAt (nodes, pidx3));
    triangles . push_back (t);
  }

  printf ("constructed %d triangle objects\n", (int)triangles . size ());
  
  vector <beam> beams;
  for (int i = 0; i < triangles . size (); ++i) {
    triangle t1 = triangles [i];
    // all triangle edges are beams
    addUniqueBeam(beams, beam(t1.p1, t1.p2));
    addUniqueBeam(beams, beam(t1.p2, t1.p3));
    addUniqueBeam(beams, beam(t1.p3, t1.p1));
    // for adjacent, co-planar triangles,
    // add the beam between the opposing points
    for (int j = 0; j < triangles . size (); ++j) {
      // continue if same triangle
      if (i == j) continue;
      triangle t2 = triangles [j];
      // continue if not co-planer
      if (! t2.coplaner(t1)) {
        //printf ("not coplanar\n");
        continue;
      }
      beam theBeam;
      if (oppositePoints(theBeam, t1, t2)) {
          printf ("opposite beam\n");
          addUniqueBeam(beams, theBeam);
      } else {
        printf ("not opposite\n");
      }
    }
  }

  printf ("extracted %d beams\n", (int)beams . size ());

  return 0;
}

// Any two triangles that have the same normal and share a line,
// should have a beam between the vertex opposite the shared line.
