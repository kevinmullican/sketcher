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
#include <algorithm>
#include <sstream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "tinyxml2.h"


unsigned int node_weight = 10;
double coef_friction = 0.7;
unsigned int spring = 2000000;
unsigned int damp = 200;
unsigned int deform = 80000;
unsigned int strength = 800000;
double wheel_factor = 0.05;
unsigned int wheel_lock = 460;
unsigned int wheel_degrees = 25;

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
  // x = vy * z - vz * y
  // y = vz * x - vx * z
  // z = vx * y - vy * x
  vect cross(const vect &v) const {
    vect c;
    c.x = (v.y * z) - (v.z * y);
    c.y = (v.z * x) - (v.x * z);
    c.z = (v.x * y) - (v.y * x);
    return c;
  }
  double mag_sq() const {
    return x*x + y*y + z*z;
  }
  double mag() const {
    return sqrt(mag_sq());
  }
  vect norm() const {
    double l = mag();
    if (l == 0.0) return vect();
    return vect (x/l, y/l, z/l);
  }
  vect neg() const {
    return vect(-x, -y, -z);
  }
  vect &neg() {
    x = -x; y = -y; z = -z;
    return *this;
  }
  vect& operator+ (const vect &v) {
    x += v.x; y += v.y; z += v.z;
    return *this;
  }
  vect& operator- (const vect &v) {
    x -= v.x; y -= v.y; z -= v.z;
    return *this;
  }
  vect operator+ (const vect &v) const {
    return vect (x + v.x, y + v.y, z + v.z);
  }
  vect operator- (const vect &v) const {
    return vect (x - v.x, y - v.y, z - v.z);
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
  void print(bool newline = false) const {
    printf ("%0.2f, %0.2f, %0.2f", x, y, z);
    if (newline) printf("\n");
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
  void print(bool newline = false) const {
    printf ("p1: [");
    p1.print();
    printf("] p2: [");
    p2.print();
    printf ("]");
    if (newline) printf ("\n");
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
  vect normal() const {
    vect pp2 = p2 - p1;
    vect pp3 = p3 - p1;
    vect cr = pp2 . cross (pp3);
    return cr.norm();
  }
  bool touching(const triangle &t) const {
    if (sharedPoints(t) == 1) return true;
    return false;
  }
  bool adjacent(const triangle &t) const {
    if (sharedPoints(t) == 2) return true;
    return false;
  }
  bool same(const triangle &t) const {
    if (sharedPoints(t) == 3) return true;
    return false;
  }
  bool sameOrientation(const triangle &t) const {
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
  bool contains (const beam &b) const {
    beam b12 (p1, p2);
    if (b12 == b) return true;
    beam b13 (p1, p3);
    if (b13 == b) return true;
    beam b23 (p2, p3);
    if (b23 == b) return true;
    return false;
  }
  bool isLongest (const beam &b) const {
    if (! contains (b)) return false;
    vect vb = b.p1 - b.p2;
    double mb = vb . mag_sq ();
    vect v12 = p1 - p2;
    double m12 = v12 . mag_sq ();
    vect v13 = p1 - p3;
    double m13 = v13 . mag_sq ();
    vect v23 = p2 - p3;
    double m23 = v23 . mag_sq ();
    if (mb >= m12 && mb >= m13 && mb >= m23)
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
  void print (bool newline = false) const {
    printf("p1: [");
    p1.print();
    printf("] p2: [");
    p2.print();
    printf("] p3: [");
    p3.print();
    printf("]");
    if (newline) printf ("\n");
  }
};

bool containsBeam(const vector <beam> &beams, const beam &theBeam) {
  for (int i = 0; i < beams . size (); ++i) {
    if (beams[i] == theBeam)
      return true;
  }
  return false;
}

bool addUniqueBeam(vector <beam> &beams, const beam &theBeam) {
  if (containsBeam(beams, theBeam))
    return false;
  beams . push_back (theBeam);
  return true;
}

bool squarePoints(beam &opposite, beam &shared, const triangle &t1, const triangle &t2) {
  if (t1 . sharedPoints (t2) != 2)
    return false;

  vect o1;
  if (! t2.contains(t1.p1))
    o1 = t1.p1;
  else if (! t2.contains(t1.p2))
    o1 = t1.p2;
  else if (! t2.contains(t1.p3))
    o1 = t1.p3;

  vect o2;
  if (! t1.contains(t2.p1))
    o2 = t2.p1;
  else if (! t1.contains(t2.p2))
    o2 = t2.p2;
  else if (! t1.contains(t2.p3))
    o2 = t2.p3;

  vect s1, s2;
  if (t1.p1 == o1) {
      s1 = t1.p2;
      s2 = t1.p3;
  } else if (t1.p2 == o1) {
      s1 = t1.p1;
      s2 = t1.p3;
  } else if (t1.p3 == o1) {
      s1 = t1.p1;
      s2 = t1.p2;
  }

  opposite.set(o1, o2);
  shared.set(s1, s2);

  return true;
}

vector <beam> extractBeams (const vector <triangle> &triangles) {
  int tsize = (int)triangles . size ();

  vector <beam> beams;
  for (int i = 0; i < tsize; ++i) {
    triangle t1 = triangles [i];
    // all triangle edges are beams
    addUniqueBeam(beams, beam(t1.p1, t1.p2));
    addUniqueBeam(beams, beam(t1.p2, t1.p3));
    addUniqueBeam(beams, beam(t1.p3, t1.p1));

    // for adjacent, co-planar triangles,
    // add the beam between the opposing points
    for (int j = 0; j < tsize; ++j) {
      // continue if same triangle
      if (i == j) continue;
      triangle t2 = triangles [j];

      // continue if not co-planer
      if (! t2.sameOrientation(t1)) continue;

      // create a beam, if it is indeed a beam
      beam opposite;
      beam shared;
      if (squarePoints(opposite, shared, t1, t2)) {
          if (t1.isLongest(shared) && t2.isLongest(shared))
            if (addUniqueBeam(beams, opposite)) {
              printf("found cross beam: ");
              opposite.print(true);
            }
      }
    }
  }

  printf ("extracted %d beams\n", (int)beams . size ());
  return beams;
}

vector <triangle> extractTriangles (const vector <unsigned int> &tridx, const vector <vect> &nodes) {
  int node_count = (int)nodes . size ();

  int tri_points = (int)tridx . size ();
  int tri_count = tri_points / 3;
  if (tri_count * 3 != tri_points)
    printf ("incomplete triangle count: %d", tri_count);

  vector <triangle> triangles;
  for (int i = 0; i < tri_count; ++i) {
    int idx = i * 3;
    int vidx1 = tridx [idx];
    int vidx2 = tridx [idx+1];
    int vidx3 = tridx [idx+2];

    if (vidx1 >= node_count ||
        vidx2 >= node_count ||
        vidx3 >= node_count) {
      printf ("triangle vertex index out of node range: %d, %d, %d > %d\n",
               vidx1, vidx2, vidx3, node_count);
      continue;
    }

    triangle t (nodes[vidx1],
                nodes[vidx2],
                nodes[vidx3]);
    triangles . push_back (t);
  }

  printf ("extracted %d triangles\n", (int)triangles . size ());
  return triangles;
}

vector <vect> extractNodes (const vector <double> &node_dims) {
  int node_elems = (int)node_dims . size ();
  int node_count = node_elems / 3;
  if (node_count * 3 != node_elems)
    printf ("incomplete node count: %d", node_count);

  vector <vect> nodes;
  for (int i = 0; i < node_count; ++i) {
    int idx = i * 3;
    vect v (node_dims [idx],
            node_dims [idx+1],
            node_dims [idx+2]);
    nodes . push_back (v);
  }

  printf ("extracted %d nodes\n", (int)nodes . size ());
  return nodes;
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

string lower (string s) {
  transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

void writeNodes (FILE *fp, const vector <vect> &nodes, const string &group, const char pfx) {
  if (! fp) return;
  int ns = (int)nodes . size ();
  if (! ns) return;
  fprintf (fp, "        {\"group\":\"%s\"},\n",
                group.c_str());

  for (int i = 0; i < ns; ++i) {
      vect n = nodes [i];
      fprintf (fp, "        [\"%c%d\",%0.3f,%0.3f,%0.3f],\n", pfx, i, n.x, n.y, n.z);
  }

}

void writeBeams (FILE *fp, const vector <beam> &beams,
                 const vector <vect> &first, const char first_char,
                 const vector <vect> &second, const char second_char,
                 unsigned int spring, unsigned int damp,
                 unsigned int deform = 0, unsigned int strength = 0)
{ if (! fp) return;
  int bs = (int)beams . size ();
  if (! bs) return;

  string def = "FLT_MAX";
  if (deform) {
      ostringstream convert;
      convert << deform;
      def = convert.str();
  }

  string strn = "FLT_MAX";
  if (strength) {
      ostringstream convert;
      convert << strength;
      strn = convert.str();
  }

  fprintf (fp, "        {\"beamSpring\":%u,\"beamDamp\":%u},\n"
               "        {\"beamDeform\":\"%s\",\"beamStrength\":\"%s\"},\n",
               spring,
               damp,
               def.c_str(),
               strn.c_str()
          );

  for (int i = 0; i < bs; ++i) {
      beam b = beams [i];
      vect &p1 = b.p1;
      vect &p2 = b.p2;

      int i1 = 0;
      char i1_set = 0;

      int i2 = 0;
      char i2_set = 0;

      for (int j = 0; j < first . size (); ++j) {
        vect n = first [j];
        if (n == p1) { i1 = j; i1_set = first_char; }
        if (n == p2) { i2 = j; i2_set = first_char; }
      }

      for (int j = 0; j < second . size (); ++j) {
        vect n = second [j];
        if (! i1_set && n == p1) { i1 = j; i1_set = second_char; }
        if (! i2_set && n == p2) { i2 = j; i2_set = second_char; }
      }

      fprintf(fp, "        [\"%c%d\",\"%c%d\"],\n", i1_set, i1, i2_set, i2);
  }
}

void writeMaterial (FILE *mat, const string &body) {
  if (! mat) return;

  string pic = body + ".png";

  fprintf (mat, "singleton Material(%s)\n"
                "{\n"
                "    mapTo = \"%s\";\n"
                "    diffuseMap[0] = \"%s\";\n"
                "    specularPower[0] = \"15\";\n"
                "    useAnisotropic[0] = \"1\";\n"
                "    castShadows = \"1\";\n"
                "    translucent = \"0\";\n"
                "    alphaTest = \"0\";\n"
                "    alphaRef = \"0\";\n"
                "}\n",
                body . c_str(),
                body . c_str(),
                pic . c_str()
          );
}

bool exportJBeam (const string &author, const string &model,
                  const vector <vect> &nodes, const vector <beam> &beams,
                  const vector <vect> &axle_nodes, const vector <beam> &axle_beams,
                  const vector <beam> &steering_beams) {

  string jbeam = model + ".jbeam";

  // jbeam file
  FILE *fp = fopen (jbeam . c_str (), "w");
  if (! fp) return false;

  // body naming
  string body = model + "_body";
  string body_group = lower (body + "_g");
  char body_char = 'b';

  // axle naming
  string axles = model + "_axles";
  string axles_group = lower (axles + "_g");
  char axle_char = 'a';

  // wheel naming
  string wheel = model + "_wheel";

  string wbase = "Wheel";

  string wheel_fl = wbase + "_FL";
  string wheel_fr = wbase + "_FR";
  string wheel_rl = wbase + "_RL";
  string wheel_rr = wbase + "_RR";

  string wheel_fl_group = lower (wheel_fl + "_g");
  string wheel_fr_group = lower (wheel_fr + "_g");
  string wheel_rl_group = lower (wheel_rl + "_g");
  string wheel_rr_group = lower (wheel_rr + "_g");

  // header
  fprintf (fp, "{\"%s\":\n"
               "\n"
               "{\n"
               "    \"information\":{\n"
               "         \"authors\":\"%s\",\n"
               "         \"name\":\"%s\",\n"
               "    }\n"
               "\n"
               "    \"slotType\" : \"main\",\n"
               "\n"
               "    \"flexbodies\": [\n"
               "        [\"mesh\", \"[group]:\", \"nonFlexMaterials\"],\n"
               "        [\"%s\", [\"%s\"]],\n"
               "        [\"%s\", [\"%s\"]],\n"
               "        [\"%s\", [\"%s\"]],\n"
               "        [\"%s\", [\"%s\"]],\n"
               "        [\"%s\", [\"%s\"]],\n"
               "    ],\n"
               "\n",
                model.c_str(),
                author.c_str(),
                model.c_str(),
                body.c_str(),
                body_group.c_str(),
                wheel_fl.c_str(),
                wheel_fl_group.c_str(),
                wheel_fr.c_str(),
                wheel_fr_group.c_str(),
                wheel_rl.c_str(),
                wheel_rl_group.c_str(),
                wheel_rr.c_str(),
                wheel_rr_group.c_str()
          );

  // nodes
  fprintf (fp, "    \"nodes\": [\n"
               "        [\"id\", \"posX\", \"posY\", \"posZ\"],\n"
               "        {\"nodeWeight\":%u},\n"
               "        {\"frictionCoef\":%0.2f},\n"
               "        {\"nodeMaterial\":\"|NM_METAL\"},\n"
               "        {\"collision\":true},\n"
               "        {\"selfCollision\":true},\n",
               node_weight,
               coef_friction);

  writeNodes (fp, nodes, body_group, 'b');
  writeNodes (fp, axle_nodes, axles_group, 'a');

  fprintf (fp, "    ],\n"
               "\n");

  // beams
  fprintf (fp, "    \"beams\": [\n"
               "        [\"id1:\", \"id2:\"],\n");

  vector <vect> mt;
  writeBeams (fp, beams, nodes, body_char, mt, 0, spring, damp, deform, strength);
  writeBeams (fp, axle_beams, nodes, body_char, axle_nodes, axle_char, spring, damp);

  fprintf (fp, "    ],\n"
               "\n");

  // steering hydros
  fprintf (fp, "    \"hydros\": [\n"
               "        [\"id1:\", \"id2:\"],\n");

  for (int i = 0; i < steering_beams . size (); ++i) {
    beam b = steering_beams [i];
    vect &p1 = b.p1;
    vect &p2 = b.p2;
    int i1 = 0;
    int i2 = 0;
    for (int j = 0; j < axle_nodes . size (); ++j) {
      vect n = axle_nodes [j];
      if (n == p1) i1 = j;
      if (n == p2) i2 = j;
    }
    fprintf (fp, "        [\"%c%d\",\"%c%d\",{\"factor\":%0.2f,\"steeringWheelLock\":%u,\"lockDegrees\":%u}],\n",
                 axle_char,
                 i1,
                 axle_char,
                 i2,
                 wheel_factor,
                 wheel_lock,
                 wheel_degrees);
  }

  fprintf (fp, "    ],\n"
               "\n");

  // TODO: write "hubWheels"
  // TODO: write "enginetorque"
  // TODO: write "engine"

  // footer
  fprintf (fp, "}\n");
  fprintf (fp, "}\n");
  fclose (fp);

  // info file
  FILE *info = fopen ("info.json", "w");
  if (! info) return false;

  fprintf (info, "{\n"
                 "    \"Name\":\"%s\",\n"
                 "    \"Author\":\"%s\",\n"
                 "    \"Type\":\"Car\",\n"
                 "    \"default_pc\":\"default\",\n"
                 "    \"colors\":{\n"
                 "        \"Pearl White\": \"1 1 1 1\"\n"
                 "    }\n"
                 "}\n",
                 model . c_str(),
                 author . c_str()
          );

  fclose (info);

  // material file
  FILE *mat = fopen ("material.cs", "w");
  if (! mat) return false;

  writeMaterial (mat, body);
  fprintf (mat, "\n");
  writeMaterial (mat, wheel);

  fclose (mat);

  return true;
}

int main (int argc, char **argv) {
  string fname;
  string model;
  string author;

  int acm1 = argc - 1;
  for (int i = 1; i < acm1; ++i) {
    if (! strncmp ("-f", argv[i], 2))
      fname = argv[i+1];
    else if (! strncmp ("-m", argv[i], 2))
      model = argv[i+1];
    else if (! strncmp ("-n", argv[i], 2))
      author = argv[i+1];
  }

  if (fname . empty () ||
      model . empty () ||
      author . empty ()) {
    printf ("usage: %s -f <input_filename> -m <model_name> -n <author_name>\n",
             argv[0]);
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

  // search for the mesh
  vector <string> mh;
  mh . push_back ("library_geometries");
  mh . push_back ("geometry");
  mh . push_back ("mesh");
  XMLElement *mesh = FindElement (collada, mh);
  if (! mesh) {
    printf ("unable to find the mesh\n");
    return 5;
  }

  // find the first set of sources, which are the node dimensions
  vector <string> fh;
  fh . push_back ("source");
  fh . push_back ("float_array");
  XMLElement *fa = FindElement (mesh, fh);
  if (! fa) {
    printf ("unable to find the float_array XML element\n");
    return 6;
  }

  // find the first set of triangles, which are the triangle
  // vertex indices in the nodes (node dimensions mod 3)
  vector <string> th;
  th . push_back ("triangles");
  XMLElement *tri = FindElement (mesh, th);
  if (! tri) {
    printf ("unable to find the triangles XML element\n");
    return 7;
  }

  // the actual triangle points are here
  vector <string> tvh;
  tvh . push_back ("p");
  XMLElement *tri_vert = FindElement (tri, tvh);
  if (! tri_vert) {
    printf ("unable to find the triangle vertex XML element\n");
    return 8;
  }

  // check the node dim count and parse the node dims
  int want = fa -> IntAttribute ("count");
  string node_text = fa -> GetText ();
  vector <double> node_dims = DoubleSplit (node_text, " ");
  int node_elems = (int)node_dims . size ();
  if (node_elems != want)
    printf ("node element want count %d not equal to got count %d\n", want, node_elems);
  printf ("found %d node elements\n", node_elems);

  // check the triangle count and parse the triangle indices
  want = tri -> IntAttribute ("count");
  string tri_text = tri_vert -> GetText ();
  vector <unsigned int> tridx = UintSplit (tri_text, " ");
  int tri_points = (int)tridx . size ();
  if (tri_points / 3 != want)
    printf ("triangle index want count %d not equal to got count %d\n", want, tri_points / 3);
  printf ("found %d triangle indices\n", tri_points);

  vector <vect> nodes = extractNodes (node_dims);
  vector <triangle> triangles = extractTriangles (tridx, nodes);
  vector <beam> beams = extractBeams (triangles);

  if (! mkdir (model . c_str(), 0755))
      chdir (model . c_str());

  vector <vect> mt_vect;
  vector <beam> mt_beam;

  if (exportJBeam (author, model, nodes, beams, mt_vect, mt_beam, mt_beam))
    printf ("successfully exported model %s\n", model . c_str());
  else
    printf ("error exporting %s\n", model . c_str());

  return 0;
}

// Any two triangles that have the same normal and share a line,
// should have a beam between the vertex opposite the shared line.
