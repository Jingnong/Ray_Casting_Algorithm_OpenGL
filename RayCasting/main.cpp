//
//  main.cpp
//  Assignment3_290_1124
//
//  Created by Jingnong Wang on 11/24/17. Student ID: 1281672
//  Copyright © 2017 Jingnong Wang. All rights reserved.
//

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<GLUT/glut.h>

using namespace std;

#define  PI  2*acos(0)

/* data structures */
typedef struct{
  float x;
  float y;
  float z;
} Point; // 3D point

typedef struct{
  float x;
  float y;
  float z;
} Vector; // 3D vector

typedef struct{
  float r;
  float g;
  float b;
} RGB_float; // color

typedef struct sphere{
  Point center;
  float radius;
  float color[3];
} Sphere; // Sphere

/*  Viewing parameters.  */
Point from, at, up;
float  VXR, VXL, VYB, VYT;
float  ax, ay, az, bx, by, bz, cx, cy, cz;
float  viewangle, angle, tanv2;
float  xinterval, yinterval;

/*  Illumination parameters.  */
Point  light;
RGB_float  il, ia;
RGB_float  ka1, kd1, ks1;
RGB_float  ka2, kd2, ks2;
int  phong1, phong2;

/*  Image parameters.  */
int    win_width, win_height;

float *texture_R;
float *texture_G;
float *texture_B;

/*  Object parameters.  */
Sphere spheres[2];
int numSphere;

/* Utility functions */

/* Normalizes the given vector. */
void Normalize(float *x,float *y,float *z)
{
  float  norm;
  
  norm = sqrt( *x * *x + *y * *y + *z * *z );
  if (norm != 0.0) {
    *x = *x / norm;
    *y = *y / norm;
    *z = *z / norm;
  }
}

/* Computes dot product of two vectors */
float dot(float v0[3], float v1[3])
{
  return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] );
}

/* Computes the power of the given base and exponent. */
float Power(float base,int exp)
{
  int  i;
  float  value;
  
  value = 1.0;
  
  for (i=1; i<=exp; i++) {
    value *= base;
  }
  
  return( value );
}

/* Computes the transformation matrix to be used in the perspective viewing model. */
void Compute_M()
{
  
  /*  Compute the line-of-sight vector, c.  */
  cx = at.x - from.x;
  cy = at.y - from.y;
  cz = at.z - from.z;
  Normalize(&cx, &cy, &cz);
  
  /*  Compute the cross product of vector c and the up vector.  */
  ax = cy*up.z - up.y*cz;
  ay = up.x*cz - cx*up.z;
  az = cx*up.y - up.x*cy;
  Normalize(&ax, &ay, &az);
  
  /*  Compute the cross product of vector a and c.  */
  bx = ay*cz - cy*az;
  by = cx*az - ax*cz;
  bz = ax*cy - cx*ay;
}

/* Computes the transformation matrix for the perspective viewing model, and sets up the default illumination parameters. */
void Setup_Parameters()
{
  /*  Viewing parameters.  */
  from.x = 12.0;
  from.y = 7.0;
  from.z = 25.0;
  
  at.x = 6.0;
  at.y = 2.0;
  at.z = 0.0;
  
  up.x = 0.0;
  up.y = 1.0;
  up.z = 0.0;
  
  /* Veiwing angle */
  viewangle = 100.0;
  angle = viewangle * PI/180.0;
  tanv2 = tan(angle/2.0);
  
  /* Veiwing port */
  VXL = 0.5;
  VXR = -0.5;
  VYB = 0.5;
  VYT = -0.5;
  
  /* Light */
  light.x = 5.0;
  light.y = 2.0;
  light.z = 2.0;
  
  /*  Sphere parameters.  */
  spheres[0].center.x = 7.0;
  spheres[0].center.y = 6.0;
  spheres[0].center.z = 7.0;
  spheres[0].radius = 2.0;
  spheres[1].center.x = 2.0;
  spheres[1].center.y = 2.0;
  spheres[1].center.z = 7.0;
  spheres[1].radius = 3.0;
  numSphere = 2;
  
  /* Window size */
  win_width = 500;
  win_height = 500;
  
  /* Texture */
  texture_R = new float [win_width * win_height];
  texture_G = new float [win_width * win_height];
  texture_B = new float [win_width * win_height];

  /*  Compute the transformation matrix for converting world coordinates to eye coordinates.  */
  Compute_M();

  /*  Normalized the given directional light vector. */
  Normalize(&light.x, &light.y, &light.z);

  /*  Set up the conversion factors for converting from pixel coordinates to view port coordinates.  */
  xinterval = (VXR - VXL) / win_width;
  yinterval = (VYT - VYB) / win_height;

  /*  Set up default illumination (Phong lighting) parameters.  */
  il.r = 1.0;  il.g = 1.0;  il.b = 1.0;
  ia.r = 1.0;  ia.g = 1.0;  ia.b = 1.0;

  /*  Phong lighting parameters for the spheres.  */
  ka1.r = 0.0;  ka1.g = 0.2;  ka1.b = 1.0;
  kd1.r = 0.5;  kd1.g = 0.4;  kd1.b = 1.0;
  ks1.r = 0.5;  ks1.g = 1.0;  ks1.b = 1.0;
  phong1 = 100;
  ka2.r = 0.0;  ka2.g = 0.2;  ka2.b = 1.0;
  kd2.r = 0.5;  kd2.g = 0.4;  kd2.b = 1.0;
  ks2.r = 0.5;  ks2.g = 1.0;  ks2.b = 1.0;
  phong2 = 100;
}

/* Check if the give ray intercepts the given sphere. */
void Check_Sphere(float px,float py,float pz,float dx,float dy,float dz,float xc,float yc,float zc,float r,
                  float *t1, float *t2)
{
  float  a, b, c, xdiff, ydiff, zdiff, discr;
  
  xdiff = px-xc;
  ydiff = py-yc;
  zdiff = pz-zc;
  a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
  b = 2.0*( dx*xdiff + dy*ydiff + dz*zdiff );
  c = dx*dx + dy*dy + dz*dz;
  
  /*  Check if there are any intersections.  */
  
  discr = b*b - 4.0*a*c;
  if (discr < 0.0) {
    *t1 = -1.0;
    *t2 = -1.0;
  }
  else if (discr == 0.0) {
    *t1 = -b / (2.0*c);
    *t2 = -1.0;
  }
  else {
    discr = sqrt(discr);
    *t1 = (-b + discr) / (2.0*c);
    *t2 = (-b - discr) / (2.0*c);
  }
}

/* Checks if there is any object between light and  the point */
int Check_Shadow(float px, float py, float pz){
  
  float t1, t2;
  
  for(int i = 0; i<numSphere; i++){
    Check_Sphere(px, py, pz, light.x, light.y, light.z, spheres[i].center.x, spheres[i].center.y, spheres[i].center.z , spheres[i].radius, &t1, &t2);
    if(t1>= 0.0001 || t2>=0.0001) {
      return 1;
    }
  }
  return 0;
}

/* Computes the intersection of ray with an object. */
void Compute_Intersection(float px,float py,float pz,float dx,float dy, float dz,float t,float *newx,float *newy,float *newz)
{
  *newx = px + t*dx;
  *newy = py + t*dy;
  *newz = pz + t*dz;
}

/* Computes the intensity of the color for the given location based on the Phong lighting model. */
void Compute_Color(int shadow_flag, float ipx,float ipy,float  ipz,float  nx,float  ny,float  nz,
                   RGB_float ia,RGB_float ka,RGB_float kd, RGB_float ks,int n,float *r,float *g, float *b)
{
  float  vx, vy, vz, rx, ry, rz;
  float  ndotl, vdotr, cosalphapower;
  
  /*  Compute the view vector.  */
  
  vx = from.x - ipx;
  vy = from.y - ipy;
  vz = from.z - ipz;
  Normalize(&vx, &vy, &vz);
  
  /*  Compute the R (reflection) vector.  */
  
  ndotl = nx*light.x + ny*light.y + nz*light.z;
  rx = 2.0*ndotl*nx - light.x;
  ry = 2.0*ndotl*ny - light.y;
  rz = 2.0*ndotl*nz - light.z;
  
  /* Compute the V (view) vector. */
  
  vdotr = vx*rx + vy*ry + vz*rz;
  
  /* Compute Ia * Ka.  */
  
  *r = ia.r * ka.r;
  *g = ia.g * ka.g;
  *b = ia.b * ka.b;
  
  /* Compute diffuse reflection. */
  if (ndotl >= 0.0 && shadow_flag==0) {
    
    /*  diffuse reflection = kd * N dot L * Il  */
    *r = *r + kd.r*ndotl*il.r;
    *g = *g + kd.g*ndotl*il.g;
    *b = *b + kd.b*ndotl*il.b;
    
    if (vdotr >= 0.0) {
      /*  specular reflection = ks * cos(alpha)**K^n * Il */
      cosalphapower = Power(vdotr, n);
      *r = *r + ks.r*cosalphapower*il.r;
      *g = *g + ks.g*cosalphapower*il.g;
      *b = *b + ks.b*cosalphapower*il.b;
    }
  }
  
  /*  Make sure that the color is within range.  */
  if (*r > 1.0) *r = 1.0;
  if (*g > 1.0) *g = 1.0;
  if (*b > 1.0) *b = 1.0;
}

/* Compute color for each pixel. */
RGB_float Ray_Tracer_recur(int level, Point frm, float dx, float dy, float dz, float r_ind, int obj_type, int obj_id){
  
  int    obj, obj_num, shadow_flag, texture, i;
  float  nx, ny, nz;
  float  t_min, t1, t2, ipx, ipy, ipz;
  float  r, g, b;
  RGB_float newcolor;
  
  t_min = 999.0;
  obj_num = 0;
  obj = 0;
  texture = 0;
  
  /*  Check if the current ray intercepts spheres  */
  
  for(i = 0; i < numSphere; i++){
    if(obj_type == 1 && obj_id == i)continue;
    Check_Sphere(frm.x, frm.y, frm.z, dx, dy, dz, spheres[i].center.x, spheres[i].center.y, spheres[i].center.z ,spheres[i].radius, &t1, &t2);
    
    if (t1 >= 0.0) {
      t_min = t1;
      obj = 1;
      obj_num = i;
      Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t1, &ipx, &ipy, &ipz);
    }
    
    if (t2 >= 0.0 && t2 < t_min) {
      t_min = t2;
      obj = 1;
      obj_num = i;
      Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t2, &ipx, &ipy, &ipz);
    }
  }
  
  if(fabs(frm.x - ipx) < 0.00001 && fabs(frm.y - ipy) < 0.00001 && fabs(frm.z - ipz) < 0.00001)
  {
    obj = 0;
  }
  
  /*  Compute the intensity to use at the current pixel.  */
  
  /*  The current ray does not intersect any of the objects.  */
  if(obj == 0)
  {
    r = 0.0;
    g = 0.0;
    b = 0.0;
  }
  
  /*  The current ray intercept spheres.  */
  else if(obj == 1)
  {
    nx = ipx - spheres[obj_num].center.x;
    ny = ipy - spheres[obj_num].center.y;
    nz = ipz - spheres[obj_num].center.z;
    Normalize(&nx, &ny, &nz);
    
    shadow_flag = 0;
    shadow_flag = Check_Shadow(ipx, ipy, ipz );
    texture = 0;
    
    if(obj_num == 0){
      Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b);
    }
    if(obj_num == 1){
      Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b);
    }
  }
  
  Point newFrom;
  newFrom.x = ipx;
  newFrom.y = ipy;
  newFrom.z = ipz;
  
  newcolor.r = r ;
  newcolor.g = g ;
  newcolor.b = b ;
  return newcolor;
}

// Generates the primary ray from camera to every pixel of the window and call  ray tracer
void Ray_Generate(){
  
  int    i, j;
  int    buf_ptr;
  float  xv, yv, dx, dy, dz;
  float  u, v;
  RGB_float  newcolor;
  
  
  /*  Generate a ray for each pixel in the desired image.  */
  buf_ptr = 0;
  for (i = 0; i < win_width; i++) {
    u = (float)i/win_width;
    
    for (j = 0; j < win_height; j++) {
      v = (float)j/win_height;
      
      /*  Compute the corresponding view port coordinates.  */
      xv = VXL + i * xinterval;
      yv = VYB + j * yinterval;
      
      /*  Compute the direction of the current ray from the "From" point to the current position on the image.  */
      dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
      dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
      dz = az*xv*tanv2 + bz*yv*tanv2 + cz;
      
      newcolor = Ray_Tracer_recur(1, from, dx, dy, dz, 1, -1, -1);
      
      /* Save the computed color intensity to the image buffer. */
      texture_R[i + win_width * j] = newcolor.r;
      texture_G[i + win_width * j] = newcolor.g;
      texture_B[i + win_width * j] = newcolor.b;
      
    }
  }
}

/* Initialize the projection matrix.  */
void myinit(void)
{
  /* attributes */
  
  glClearColor(0.0, 0.0, 0.0, 1.0);
  
  /* set up viewing */
  /* 512 x 512 window with origin lower left */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 512.0, 0.0, 512.0);
  glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image.   A more efficient method is to use glDrawPixels(). */
void display( void )
{
  int s, t;
  float  r, g, b;
  
  glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */
  
  for(t = 0; t < win_height; t++) {
    for(s = 0; s < win_width; s++) {
      
      r = texture_R[s + win_width * t];
      g = texture_G[s + win_width * t];
      b = texture_B[s + win_width * t];
      
      glColor3f(r, g, b);
      glBegin(GL_POINTS);
      glVertex2f(s,t);
      glEnd();
    }
  }
  glFlush(); /* clear buffers */
}


int main(int argc, char**argv)
{
  /* GLUT initialization */
  glutInit(&argc,argv);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(500,500);
  glutInitWindowPosition(0,0);
  glutCreateWindow("Ray Casting");
  glutDisplayFunc(display);
  
  Setup_Parameters();
  Ray_Generate();
    
  myinit();
  glutMainLoop();
  
  return(0);
}
