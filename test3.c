/* Output from p2c, the Pascal-to-C translator */
/* From input file "test3.pas" */

/****************************************************************************/
/* program : hidden lines                                                   */
/* purpose : This does hidden lines on the three dimentional package.       */
/****************************************************************************/

/*       NOTICE  NOTICE  NOTICE  NOTICE  NOTICE  NOTICE  NOTICE             */

/* This has been a collaborative attempt (of debugging) by  Dan Fields and  */
/* Fred Meczywor.  Fred did the original typing and debugging of the code.  */
/* Subsequently, Dan Fields did the (real difficult debugging) of the math  */
/* stuff. Therefore, this is a dual handing in of the project.              */

#include "p2c.h"

#define MAXSEGS 50
#define MAX_VERTS 250
#define MAX_LINES 350
#define MAX_FACETS 150
#define MAX_SIDES 6

#define KLUDGE 0.00001

typedef enum
{
  X,
  Y,
  Z,
  W
} XYZW;
typedef double POINT[4];
typedef double MATRIX[4][4];
enum
{
  LEFT,
  RIGHT,
  BOTTOM,
  TOP,
  FRONT,
  BACK
};
typedef long DIRECTION;

typedef POINT POINT_TYPE[MAX_VERTS];
typedef double RM_TYPE[MAXSEGS][2];
typedef double COORD_TYPE[MAX_VERTS];
typedef uchar LINV_TYPE[MAX_LINES][2];
typedef uchar LINF_TYPE[MAX_FACETS][MAX_SIDES];
typedef uchar INDEX_TYPE[MAX_FACETS];

typedef struct OBJ_TYPE
{
  uchar NOV;
  POINT_TYPE POINTS;
  COORD_TYPE XP, YP;
  uchar NOL;
  LINV_TYPE LINV;
  uchar NOF;
  INDEX_TYPE INDEXF;
  LINF_TYPE LINF;
} OBJ_TYPE;

Static FILE *NIB;
Static POINT VRP, VPN, VUP, COP, P1, P2, VRPPRIME;
Static double D, F, B, UMIN, VMIN, UMAX, VMAX, MINDIST;
Static MATRIX TRANSMATRIX, THEMATRIX, VIEWMATRIX;
Static long NOWX, NOWY, COLOR;
Static double WX, WY, WL, WH, VH, VL;
Static OBJ_TYPE OBJ;
Static RM_TYPE RM;

/********************************************************************/

static double
stretch(double x1, double x2)
{
  return (x1 + (x2 - x1) * 319 / 199);
}

/********************************************************************/

static long
tox(double x)
{
  return ((long)floor(VL * (x - WX) / WL + 0.5));
}

/********************************************************************/

static long
toy(double y)
{
  return ((long)floor(VH * (WY - y) / WH + VH + 0.5));
}

/********************************************************************/

static void
window(double left, double right, double bottom, double top)
{
  WX = left;
  WL = right - left;
  WY = bottom;
  WH = top - bottom;
}

/********************************************************************/

static void
viewport(double left, double right, double bottom, double top)
{
  long X1, Y1, X2, Y2;

  VL = (right - left) * 319;
  VH = (top - bottom) * 199;
  X1 = (long)floor(left * 319 + 0.5);
  Y1 = (long)floor((1 - top) * 199 + 0.5);
  X2 = (long)floor(right * 319 + 0.5);
  Y2 = (long)floor((1 - bottom) * 199 + 0.5);
  /*  GRAPHWINDOW(X1, Y1, X2, Y2); */
  /* p2c: test3.pas, line 103:
   * Warning: Symbol 'GRAPHWINDOW' is not defined [221] */
}

/********************************************************************/

static void
drawon()
{
  COLOR = 3;
}

/********************************************************************/

static void
drawoff()
{
  COLOR = 0;
}

/********************************************************************/

static void
initialize()
{
  /* PURPOSE : INITIALIZES THE MATRIX AND COLORS */
  /*  GRAPHMODE(); */
  /* p2c: test3.pas, line 125:
   * Warning: Symbol 'GRAPHMODE' is not defined [221] */
  /*  PALETTE(1); */
  /* p2c: test3.pas, line 126:
   * Warning: Symbol 'PALETTE' is not defined [221] */
  drawon();
  window(0.0, 319.0, 0.0, 199.0);
  viewport(0.0, 1.0, 0.0, 1.0);
}

/********************************************************************/

static void
plotpoint(double x, double y)
{
  /* PURPOSE : PLOTS A SINGLE POINT ONTO THE SCREEN */
  NOWX = tox(x);
  NOWY = toy(y);
  /*  PLOT(NOWX, NOWY, COLOR); */
  /* p2c: test3.pas, line 139: Warning: Symbol 'PLOT' is not defined [221] */
}

/********************************************************************/

static void
moves(double x, double y)
{
  /* PURPOSE : MOVES TO <X,Y> WITHOUT DRAWING */
  NOWX = tox(x);
  NOWY = toy(y);
}

/********************************************************************/

static void
draws(double x, double y)
{
  /* PURPOSE : DRAWS TO <X,Y> AND MAKES IT THE CURRENT POINT */
  long temp;

  temp = toy(y);
  /*  DRAW(NOWX, NOWY, tox(X_), TEMP, COLOR); */
  /* p2c: test3.pas, line 160: Warning: Symbol 'DRAW' is not defined [221] */
  NOWX = tox(x);
  NOWY = temp;
}

/********************************************************************/

static double
dot_prod(double *v1, double *v2)
{
  /* PURPOSE : TO PERFORM THE DOT PRODUCT OF V1 AND V2 */
  return (v1[(long)X] * v2[(long)X] + v1[(long)Y] * v2[(long)Y] +
          v1[(long)Z] * v2[(long)Z]);
}

/*******************************************************************/

static void
cross_prod(double *v1, double *v2, double *result)
{
  /* PURPOSE : TO PERFORM THE CROSS PROD WITH V1 AND V2 */
  /*            WITH THE RESULT IN RESULT */
  result[(long)X] = v1[(long)Y] * v2[(long)Z] - v1[(long)Z] * v2[(long)Y];
  result[(long)Y] = v1[(long)Z] * v2[(long)X] - v1[(long)X] * v2[(long)Z];
  result[(long)Z] = v1[(long)X] * v2[(long)Y] - v1[(long)Y] * v2[(long)X];
}

/********************************************************************/

Static Void NORMALIZE(V1)
double *V1;
{
  /* PURPOSE : TO NORMALIZE A VECTOR */
  double LENGTH;
  XYZW INDEX;

  LENGTH = sqrt(dot_prod(V1, V1));
  for (INDEX = X; (long)INDEX <= (long)Z; INDEX = (XYZW)((long)INDEX + 1))
    V1[(long)INDEX] /= LENGTH;
}

/********************************************************************/

Static Void MULTIPLY_MAT(A, B, C) double (*A)[4], (*B)[4], (*C)[4];
{
  /* PURPOSE : TO MULTIPLY TWO 4X4 MATRIX'S WITH RESULT IN C */
  long I, J, K;

  for (I = 0; I <= 3; I++)
  {
    for (J = 0; J <= 3; J++)
    {
      C[I][J] = 0.0;
      for (K = 0; K <= 3; K++)
        C[I][J] += A[I][K] * B[K][J];
    }
  }
}

/********************************************************************/

Static Void MULTIPLY_VEC(A, V1, V2) double (*A)[4];
double *V1, *V2;
{
  /* PURPOSE : TO MULTIPLY TWO VECTORS GIVING A POINT */
  XYZW I;
  long J;

  for (I = X; (long)I <= (long)W; I = (XYZW)((long)I + 1))
  {
    J = (int)I + 1;
    V2[(long)I] = V1[(long)X] * A[0][J - 1] + V1[(long)Y] * A[1][J - 1] +
                  V1[(long)Z] * A[2][J - 1] + V1[(long)W] * A[3][J - 1];
  }
}

/********************************************************************/

Static Void SUBTRACT_VEC(V1, V2, RESULT)
double *V1, *V2, *RESULT;
{
  /* PURPOSE : TO SUBTRACT TWO VECTORS GIVING A VECTOR */
  RESULT[(long)X] = V1[(long)X] - V2[(long)X];
  RESULT[(long)Y] = V1[(long)Y] - V2[(long)Y];
  RESULT[(long)Z] = V1[(long)Z] - V2[(long)Z];
}

/********************************************************************/

Static Void GET_AXIS(V)
double *V;
{
  /* PURPOSE : TO CALCULATE THE AXIS'S FROM THE SETUP PROCEDURE */
  double SCALAR;
  POINT TEMPVECT;
  XYZW I;

  NORMALIZE(VPN);
  SCALAR = dot_prod(VPN, VUP);
  for (I = X; (long)I <= (long)Z; I = (XYZW)((long)I + 1))
    TEMPVECT[(long)I] = SCALAR * VPN[(long)I];
  TEMPVECT[(long)W] = 1.0;
  SUBTRACT_VEC(VUP, TEMPVECT, V);
}

/********************************************************************/

Static Void FILL_MATRIX(A, A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4, D1,
                        D2, D3, D4) double (*A)[4];
double A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4, D1, D2, D3, D4;
{
  /* PURPOSE : TO FILL A MATRIX WITH THE PROPER VALUES */
  A[0][0] = A1;
  A[0][1] = A2;
  A[0][2] = A3;
  A[0][3] = A4;
  A[1][0] = B1;
  A[1][1] = B2;
  A[1][2] = B3;
  A[1][3] = B4;
  A[2][0] = C1;
  A[2][1] = C2;
  A[2][2] = C3;
  A[2][3] = C4;
  A[3][0] = D1;
  A[3][1] = D2;
  A[3][2] = D3;
  A[3][3] = D4;
}

/********************************************************************/

Static Void RESET_TRANS()
{
  /* PURPOSE : TO INITIALIZE THE MATRIX BACK TO AN IDENTITY MATRIX */
  FILL_MATRIX(TRANSMATRIX, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
              1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
}

/********************************************************************/

Static Void SETUP()
{
  /* PURPOSE : SETS UP THE MATRIXES AND INITIALIZES EVERYTHING */
  POINT U, V, CENT;
  MATRIX T, ROT, SH, SC;
  double SX, SY, SZ;

  GET_AXIS(V);
  cross_prod(VPN, V, U);
  FILL_MATRIX(T, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
              -(VRP[(long)X] + COP[(long)X]), -(VRP[(long)Y] + COP[(long)Y]),
              -(VRP[(long)Z] + COP[(long)Z]), 1.0);
  FILL_MATRIX(ROT, U[(long)X], V[(long)X], VPN[(long)X], 0.0, U[(long)Y],
              V[(long)Y], VPN[(long)Y], 0.0, U[(long)Z], V[(long)Z],
              VPN[(long)Z], 0.0, 0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(T, ROT, VIEWMATRIX);
  MULTIPLY_VEC(VIEWMATRIX, VRP, VRPPRIME);
  CENT[(long)Z] = VRPPRIME[(long)Z];
  CENT[(long)X] = VRPPRIME[(long)X] + (UMAX + UMIN) / 2;
  CENT[(long)Y] = VRPPRIME[(long)Y] + (VMAX + VMIN) / 2;
  FILL_MATRIX(SH, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
              -(CENT[(long)X] / CENT[(long)Z]),
              -(CENT[(long)Y] / CENT[(long)Z]), 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(VIEWMATRIX, SH, VIEWMATRIX);
  SZ = 1 / (VRPPRIME[(long)Z] + B);
  SX = 2 * VRPPRIME[(long)Z] * SZ / (UMAX - UMIN);
  SY = 2 * VRPPRIME[(long)Z] * SZ / (VMAX - VMIN);
  FILL_MATRIX(SC, SX, 0.0, 0.0, 0.0, 0.0, SY, 0.0, 0.0, 0.0, 0.0, SZ, 0.0,
              0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(VIEWMATRIX, SC, VIEWMATRIX);
  D = VRPPRIME[(long)Z] * SZ;
  window(-D, D, -D, D);
  MINDIST = (VRPPRIME[(long)Z] + F) * SZ;
}

/********************************************************************/

Static Void REGION(P, ENDPOINT)
double *P;
DIRECTION *ENDPOINT;
{
  /* PURPOSE : TO DETERMINE IF A POINT IS IN THE WINDOW */
  *ENDPOINT = 0;
  if (P[(long)X] < -P[(long)Z] - 0.001)
    *ENDPOINT = 1L << ((long)LEFT);
  else
  {
    if (P[(long)X] > P[(long)Z] + 0.001)
      *ENDPOINT = 1L << ((long)RIGHT);
  }
  if (P[(long)Y] < -P[(long)Z] - 0.001)
    *ENDPOINT |= 1L << ((long)BOTTOM);
  else
  {
    if (P[(long)Y] > P[(long)Z] + 0.001)
      *ENDPOINT |= 1L << ((long)TOP);
  }
  if (P[(long)Z] < MINDIST - 0.001)
    *ENDPOINT |= 1L << ((long)FRONT);
  else
  {
    if (P[(long)Z] > 1.001)
      *ENDPOINT |= 1L << ((long)BACK);
  }
}

/********************************************************************/

Static Void GET_P(P, T)
double *P;
double T;
{
  P[(long)X] = (P2[(long)X] - P1[(long)X]) * T + P1[(long)X];
  P[(long)Y] = (P2[(long)Y] - P1[(long)Y]) * T + P1[(long)Y];
  P[(long)Z] = (P2[(long)Z] - P1[(long)Z]) * T + P1[(long)Z];
}

/********************************************************************/

Static Void CLIP_LEFT(P)
double *P;
{
  double T;

  T = (-P1[(long)Z] - P1[(long)X]) /
      (P2[(long)X] - P1[(long)X] + P2[(long)Z] - P1[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void CLIP_RIGHT(P)
double *P;
{
  double T;

  T = (P1[(long)Z] - P1[(long)X]) /
      (P2[(long)X] - P1[(long)X] + P1[(long)Z] - P2[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void CLIP_TOP(P)
double *P;
{
  double T;

  T = (P1[(long)Z] - P1[(long)Y]) /
      (P2[(long)Y] - P1[(long)Y] - P2[(long)Z] + P1[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void CLIP_BOTTOM(P)
double *P;
{
  double T;

  T = (P1[(long)Y] + P1[(long)Z]) /
      (P1[(long)Y] - P2[(long)Y] - P2[(long)Z] + P1[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void CLIP_BACK(P)
double *P;
{
  double T;

  T = (1 - P1[(long)Z]) / (P2[(long)Z] - P1[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void CLIP_FRONT(P)
double *P;
{
  double T;

  T = (MINDIST - P1[(long)Z]) / (P2[(long)Z] - P1[(long)Z]);
  GET_P(P, T);
}

/********************************************************************/

Static Void MAKE(P, OLDP)
double *P, *OLDP;
{
  P[(long)X] = OLDP[(long)X];
  P[(long)Y] = OLDP[(long)Y];
  P[(long)Z] = OLDP[(long)Z];
  P[(long)W] = 1.0;
}

/********************************************************************/

Static Void CLIPPING(P1_, P2_)
double *P1_, *P2_;
{
  /* PURPOSE : DOES THE ACTUAL CLIPPING OF A LINE IF IT IS OUTSIDE */
  /*           THE WINDOW AND THE DRAWING TO THE SCREEN */
  POINT P1, P2, P;
  DIRECTION DIR, DIR1, DIR2;
  boolean INVISIBLE;
  double XP, YP;

  memcpy(P1, P1_, sizeof(POINT));
  memcpy(P2, P2_, sizeof(POINT));
  INVISIBLE = false;
  REGION(P1, &DIR1);
  REGION(P2, &DIR2);
  while (DIR1 != 0 || DIR2 != 0)
  {
    if ((DIR1 & DIR2) != 0)
    {
      DIR1 = 0;
      DIR2 = 0;
      INVISIBLE = true;
      continue;
    }
    if (DIR1 == 0)
    {
      DIR = DIR2;
      MAKE(P, P2);
    }
    else
    {
      DIR = DIR1;
      MAKE(P, P1);
    }
    if (((1L << ((long)LEFT)) & DIR) != 0)
      CLIP_LEFT(P);
    else
    {
      if (((1L << ((long)RIGHT)) & DIR) != 0)
        CLIP_RIGHT(P);
      else
      {
        if (((1L << ((long)TOP)) & DIR) != 0)
          CLIP_TOP(P);
        else
        {
          if (((1L << ((long)BOTTOM)) & DIR) != 0)
            CLIP_BOTTOM(P);
          else
          {
            if (((1L << ((long)BACK)) & DIR) != 0)
              CLIP_BACK(P);
            else
            {
              if (((1L << ((long)FRONT)) & DIR) != 0)
                CLIP_FRONT(P);
            }
          }
        }
      }
    }
    if (DIR == DIR1)
    {
      MAKE(P1, P);
      REGION(P1, &DIR1);
    }
    else
    {
      MAKE(P2, P);
      REGION(P2, &DIR2);
    }
  }
  if (INVISIBLE)
    return;
  XP = D * P1[(long)X] / P1[(long)Z];
  YP = D * P1[(long)Y] / P1[(long)Z];
  moves(XP, YP);
  XP = D * P2[(long)X] / P2[(long)Z];
  YP = D * P2[(long)Y] / P2[(long)Z];
  draws(XP, YP);
}

/****************************************************************************/
/* procedure : update                                                       */
/* purpose   : This is used to update the visible line segment array with   */
/*             the visible line segments for that particular line under     */
/*             study.                                                       */
/*             It leaves the RM array filled with the Mus of the endpoints  */
/*             of each visible segment.                                     */
/****************************************************************************/

Static Void UPDATE(RMIN, RMAX, NRL, RM)
double RMIN, RMAX;
uchar *NRL;
double (*RM)[2];
{
  uchar MORERL, K;
  double R1, R2;
  uchar FORLIM;

  MORERL = *NRL;
  FORLIM = *NRL;
  for (K = 0; K < FORLIM; K++)
  {
    R1 = RM[K][0];
    R2 = RM[K][1];
    if (R1 <= RMAX && R2 >= RMIN)
    {
      if (R1 >= RMIN && R2 <= RMAX)
        RM[K][0] = -1.0;
      else
      {
        if (R1 < RMIN && R2 > RMAX)
        {
          MORERL++;
          RM[MORERL - 1][0] = RMAX;
          RM[MORERL - 1][1] = R2;
          RM[K][1] = RMIN;
        }
        else
        {
          if (R1 < RMIN)
            RM[K][1] = RMIN;
          else
            RM[K][0] = RMAX;
        }
      }
    }
  } /* FOR */
  *NRL = 0;
  for (K = 0; K < MORERL; K++)
  {
    if (RM[K][0] >= -KLUDGE)
    {
      (*NRL)++;
      RM[*NRL - 1][0] = RM[K][0];
      RM[*NRL - 1][1] = RM[K][1];
    }
  }
}

/****************************************************************************/
/* function : eval                                                          */
/* synopsis : Believe it or not, this will convert a real number to units   */
/*            of integers.                                                  */
/****************************************************************************/

Static long EVAL(VAL)
double VAL;
{
  if (VAL > KLUDGE)
    return 1;
  else
  {
    if (VAL < -KLUDGE)
      return -1;
    else
      return 0;
  }
}

/***************************************************************************/
/* function : G (what a pnemonic name)                                     */
/* purpose  : This guy returns the functional for a line that is perpen-   */
/*            dicular to the line under study.  This is used by the funct- */
/*            ion offend .                                                 */
/***************************************************************************/

Static long G(XX, YY, X4, Y4, XD, YD)
double XX, YY, X4, Y4, XD, YD;
{
  double VAL;

  VAL = (YY - Y4) * YD + (XX - X4) * XD;
  return (EVAL(VAL));
}

/****************************************************************************/
/* function : visible                                                       */
/* purpose  : This one will find out of the reverse projection of the RMID  */
/*            (that is X,Y,Z hat) is on the same or different size of the   */
/*            plane that the line is in intersection with.  Essentially, it */
/*            does this by plugging the (origin) and the X,Y,Z hat point    */
/*            into the functional for the plane.  and testing the result.   */
/*            That is -1, 0, or 1.                                          */
/****************************************************************************/

Static boolean VISIBLE(J, L1, L2, X1, Y1, X2, Y2, RMIN, RMAX)
uchar J, L1, L2;
double X1, Y1, X2, Y2, RMIN, RMAX;
{
  boolean Result;
  double RMID, RXX, XMID, YMID, DENOM, PHI, ZHAT, DDD, XHAT, YHAT, D1, F1, F2;
  uchar I1, I2, V1, V2, V3;
  POINT P, P1, P3;

  Result = true;
  RMID = (RMAX + RMIN) * 0.5;
  RXX = 1.0 - RMID;
  XMID = RXX * X1 + RMID * X2;
  YMID = RXX * Y1 + RMID * Y2;
  DENOM = D * (OBJ.POINTS[L2 - 1][(long)X] - OBJ.POINTS[L1 - 1][(long)X]) -
          XMID * (OBJ.POINTS[L2 - 1][(long)Z] - OBJ.POINTS[L1 - 1][(long)Z]);
  if (fabs(DENOM) > KLUDGE)
    PHI = (XMID * OBJ.POINTS[L1 - 1][(long)Z] - D * OBJ.POINTS[L1 - 1]
                                                              [(long)X]) /
          DENOM;
  else
  {
    DENOM = D * (OBJ.POINTS[L2 - 1][(long)Y] - OBJ.POINTS[L1 - 1][(long)Y]) -
            YMID * (OBJ.POINTS[L2 - 1][(long)Z] - OBJ.POINTS[L1 - 1][(long)Z]);
    PHI = YMID * OBJ.POINTS[L1 - 1][(long)Z] - D * OBJ.POINTS[L1 - 1][(long)Y] / DENOM;
  } /* ELSE */
  ZHAT = (1.0 - PHI) * OBJ.POINTS[L1 - 1][(long)Z] + PHI * OBJ.POINTS[L2 - 1]
                                                                     [(long)Z];

  DDD = ZHAT / D;
  XHAT = XMID * DDD;
  YHAT = YMID * DDD;

  I1 = OBJ.LINF[J - 1][0];
  I2 = OBJ.LINF[J - 1][1];
  V1 = OBJ.LINV[I1 - 1][0];
  V2 = OBJ.LINV[I1 - 1][1];
  V3 = OBJ.LINV[I2 - 1][0];
  if (V1 == V3 || V2 == V3)
    V3 = OBJ.LINV[I2 - 1][1];
  SUBTRACT_VEC(OBJ.POINTS[V2 - 1], OBJ.POINTS[V1 - 1], P1);
  SUBTRACT_VEC(OBJ.POINTS[V2 - 1], OBJ.POINTS[V3 - 1], P3);
  cross_prod(P1, P3, P);
  D1 = P[(long)X] * OBJ.POINTS[V2 - 1][(long)X] + P[(long)Y] * OBJ.POINTS[V2 - 1][(long)Y] + P[(long)Z] * OBJ.POINTS[V2 - 1][(long)Z];
  F1 = P[(long)X] * XHAT + P[(long)Y] * YHAT + P[(long)Z] * ZHAT - D1;
  F2 = -D1;
  if (fabs(F1) < KLUDGE)
    Result = false;
  else
    fprintf(NIB, "eval is %12ld\n", labs(EVAL(F1) - EVAL(F2)));
  if (labs(EVAL(F1) - EVAL(F2)) >= 2) /* WITH */
    return false;

  return Result;
}

/****************************************************************************/
/* function : intersects                                                    */
/* purpose  : This will tell if the particular line under study is in       */
/*            intersection a facet.  It does this by simultaneously solving */
/*            the equasions of the line under study and each side of the    */
/*            facet under study. If there is an intersection, the endpoint  */
/*            RMIN and RMAX are changed to the points on the line of the    */
/*            intersection with the edges.                                  */
/****************************************************************************/

Static boolean INTERSECTS(J, X1, Y1, XD, YD, INS, RMIN, RMAX)
uchar J;
double X1, Y1, XD, YD, INS, *RMIN, *RMAX;
{
  boolean ANS;
  uchar K, NN, V1, V2;
  double XE, YE, XF, YF, DISK, XSI, MU;

  ANS = true;
  *RMAX = 0.0;
  *RMIN = 1.0;
  K = 1;
  while (K <= INS && ANS)
  {
    NN = OBJ.LINF[J - 1][K - 1];
    V1 = OBJ.LINV[NN - 1][0];
    V2 = OBJ.LINV[NN - 1][1];
    XE = OBJ.XP[V1 - 1] - OBJ.XP[V2 - 1];
    YE = OBJ.YP[V1 - 1] - OBJ.YP[V2 - 1];
    XF = OBJ.XP[V1 - 1] - X1;
    YF = OBJ.YP[V1 - 1] - Y1;
    DISK = XD * YE - XE * YD;
    if (fabs(DISK) > KLUDGE)
    {                                   /* THE LINES DO NOT OVERLAP */
      XSI = (XD * YF - YD * XF) / DISK; /* XSI ON FACET EDGE */
      if (XSI > -KLUDGE && XSI < 1 + KLUDGE)
      {
        MU = (YE * XF - XE * YF) / DISK; /* MU OF SEARCH LINE */
        if (*RMAX < MU)
          *RMAX = MU;
        if (*RMIN > MU)
          *RMIN = MU;
      }
    }
    else
    {
      if (fabs(XD) > KLUDGE)
      {
        XSI = XF / XD;
        if (fabs(YF - XSI * YD) < KLUDGE)
          ANS = false;
      }
      else
      {
        if (fabs(XF) < KLUDGE)
          ANS = false;
      }
    } /* ELSE */
    K++;
  } /* WHILE */

  /*zzz*/
  if (*RMIN > 1 + KLUDGE || *RMAX < KLUDGE)
  {
    ANS = false;
    return ANS;
  }
  if (*RMAX > 1.0)
    *RMAX = 1.0;
  if (*RMIN < 0.0)
    *RMIN = 0.0;
  if (*RMAX - *RMIN < KLUDGE)
    ANS = false;
  return ANS; /* WITH */

  /* ELSE */
}

/****************************************************************************/
/* function : off_end                                                       */
/* purpose  : This guy will tell if the points of the facet are off one end */
/*            of the line of study (or not). It does so by finding a line   */
/*            that is perpendicular to the line of study and then seeing if */
/*            all the endpoints of the facet are on one side of this per-   */
/*            pendicular.                                                   */
/****************************************************************************/

Static boolean OFF_END(J, INS, X3, Y3, X4, Y4, XD, YD)
uchar J, INS;
double X3, Y3, X4, Y4, XD, YD;
{
  boolean INSEG;
  uchar K, NN, VERTEX;
  long GVAL;

  INSEG = false;
  GVAL = G(X3, Y3, X4, Y4, XD, YD);
  if (GVAL == 0)
  {
    INSEG = true;
    return (!INSEG);
  } /* IF */
  K = 1;
  while (K <= INS && !INSEG)
  {
    NN = OBJ.LINF[J - 1][K - 1];
    VERTEX = OBJ.LINV[NN - 1][0];
    if (labs(GVAL - G(OBJ.XP[VERTEX - 1], OBJ.YP[VERTEX - 1], X4, Y4, XD, YD)) < 2)
      INSEG = true;
    else
    {
      VERTEX = OBJ.LINV[NN - 1][1];
      if (labs(GVAL - G(OBJ.XP[VERTEX - 1], OBJ.YP[VERTEX - 1], X4, Y4, XD,
                        YD)) < 2)
        INSEG = true;
    }
    K++;
  }                /* WHILE */
  return (!INSEG); /* WITH */
}

/****************************************************************************/
/* function : fit                                                           */
/* purpose  : This returns the functional of the points that are passed in  */
/*            through parameters.  This tells wether (or not) the points are*/
/*            on the line.                                                  */
/****************************************************************************/

Static long FIT(XX, YY, X1, Y1, XD, YD)
double XX, YY, X1, Y1, XD, YD;
{
  double VAL;

  VAL = (YY - Y1) * XD - (XX - X1) * YD;
  return (EVAL(VAL));
}

/****************************************************************************/
/* function : diff sides                                                    */
/* purpose  : This tells wether the points of a facet are on either side of */
/*            a line.  If they turn out to be on both sides (that is the    */
/*            somehow intersects the facet) the situation must be studied   */
/*            further.                                                      */
/*            Essentially used to elliminated facets and line that absol-   */
/*            utly do not intersect.                                        */
/****************************************************************************/

Static boolean DIFF_SIDE(J, INS, X1, Y1, XD, YD)
uchar J, INS;
double X1, Y1, XD, YD;
{
  boolean SAME_SIDE;
  uchar K, NN, VERTEX;
  long FVAL;

  SAME_SIDE = true;
  NN = OBJ.LINF[J - 1][0];
  VERTEX = OBJ.LINV[NN - 1][0];
  FVAL = FIT(OBJ.XP[VERTEX - 1], OBJ.YP[VERTEX - 1], X1, Y1, XD, YD);
  if (FVAL == 0)
  {
    SAME_SIDE = false;
    return (!SAME_SIDE);
  } /* IF */
  K = 2;
  while (K <= INS && SAME_SIDE)
  {
    NN = OBJ.LINF[J - 1][K - 1];
    VERTEX = OBJ.LINV[NN - 1][0];
    if (FVAL != FIT(OBJ.XP[VERTEX - 1], OBJ.YP[VERTEX - 1], X1, Y1, XD, YD))
      SAME_SIDE = false;
    else
    {
      VERTEX = OBJ.LINV[NN - 1][1];
      if (FVAL != FIT(OBJ.XP[VERTEX - 1], OBJ.YP[VERTEX - 1], X1, Y1, XD, YD))
        SAME_SIDE = false;
    }
    K++;
  }                    /* WHILE */
  return (!SAME_SIDE); /* WITH */
}

/****************************************************************************/
/* function : in facet                                                      */
/* purpose  : This tells wether the line in a member of one of the sides of */
/*            facet.  It does this by checking the indexzing numbers of the */
/*            line against the indexzing numbers of the lines that make up  */
/*            the facet.                                                    */
/****************************************************************************/

Static boolean IN_FACET(INS, I, J)
uchar INS, I, J;
{
  uchar K;
  boolean IN_IT;

  K = 1;
  IN_IT = false;
  while (K <= INS && !IN_IT)
  {
    if (OBJ.LINF[J - 1][K - 1] == I)
      IN_IT = true;
    K++;
  }
  return IN_IT; /* WITH */
}

/*****************************************************************************/
/* procedure : findsegs                                                      */
/* purpose   : The main driver for the search algorithim.  This checks the   */
/*             interaction between the I th line and J th facet.             */
/*****************************************************************************/

Static Void FIND_SEGS(I, J, L1, L2, X1, Y1, X2, Y2, XD, YD, NRL, RM)
uchar I, J, L1, L2;
double X1, Y1, X2, Y2, XD, YD;
uchar *NRL;
double (*RM)[2];
{
  uchar INS;
  double RMIN, RMAX;

  INS = OBJ.INDEXF[J - 1];
  if (IN_FACET(INS, I, J)) /* WITH */
    return;
  if (!DIFF_SIDE(J, INS, X1, Y1, XD, YD))
    return;
  if (OFF_END(J, INS, X2, Y2, X1, Y1, XD, YD))
    return;
  if (OFF_END(J, INS, X1, Y1, X2, Y2, XD, YD))
    return;
  if (INTERSECTS(J, X1, Y1, XD, YD, (double)INS, &RMIN, &RMAX))
  {
    if (VISIBLE(J, L1, L2, X1, Y1, X2, Y2, RMIN, RMAX))
      UPDATE(RMIN, RMAX, NRL, RM);
  }
}

/****************************************************************************/
/* procedure : draw segs                                                    */
/* purpose   : This guy will draw the visible segements of a line from the  */
/*             Mus in the RM array.                                         */
/****************************************************************************/

static void
    draw_segs(X1, Y1, X2, Y2, MU1, MU2) double X1,
    Y1, X2, Y2, MU1, MU2;
{
  double XP1, YP1, XP2, YP2, ONE_MU1, ONE_MU2;

  ONE_MU1 = 1 - MU1;
  ONE_MU2 = 1 - MU2;
  XP1 = X1 * ONE_MU1 + X2 * MU1;
  YP1 = Y1 * ONE_MU1 + Y2 * MU1;
  XP2 = X1 * ONE_MU2 + X2 * MU2;
  YP2 = Y1 * ONE_MU2 + Y2 * MU2;
  if (fabs(XP1 - XP2) > KLUDGE || fabs(YP1 - YP2) > KLUDGE)
  {
    moves(XP1, YP1);
    draws(XP2, YP2);
  }
}

/****************************************************************************/
/* procedure : hidden                                                       */
/* purpose   : This is the main routine of the hidden lines thing.          */
/****************************************************************************/

static void
hidden()
{
  uchar I, J, NRL, L1, L2;
  double X1, X2, Y1, Y2, XD, YD;
  uchar FORLIM;

  FORLIM = OBJ.NOL;
  for (I = 1; I <= FORLIM; I++)
  { /* WITH */
    L1 = OBJ.LINV[I - 1][0];
    L2 = OBJ.LINV[I - 1][1];
    NRL = 1;
    RM[0][0] = 0.0;
    RM[0][1] = 1.0;
    X1 = OBJ.XP[L1 - 1];
    Y1 = OBJ.YP[L1 - 1];
    X2 = OBJ.XP[L2 - 1];
    Y2 = OBJ.YP[L2 - 1];
    XD = X2 - X1;
    YD = Y2 - Y1;
    J = 1;
    while (J <= OBJ.NOF && NRL != 0)
    {
      FIND_SEGS(I, J, L1, L2, X1, Y1, X2, Y2, XD, YD, &NRL, RM);
      J++;
    }
    for (J = 0; J < NRL; J++)
      draw_segs(X1, Y1, X2, Y2, RM[J][0], RM[J][1]);
  }
}

/****************************************************************************/

Static Void GET_VIEW()
{
  /* PURPOSE : TO INPUT THE VRP AND PREPARE THE ANGLE OF VIEW */
  printf("ENTER X Y Z OF VIEW REFERENCE POINT\n");
  scanf("%lg%lg%lg%*[^\n]", VRP, &VRP[(long)Y], &VRP[(long)Z]);
  getchar();
  VRP[(long)W] = 1.0;
  VPN[(long)X] = -VRP[(long)X];
  VPN[(long)Y] = -VRP[(long)Y];
  VPN[(long)Z] = -VRP[(long)Z];
  VPN[(long)W] = 1.0;
  VUP[(long)X] = 0.0;
  VUP[(long)Y] = 1.0;
  VUP[(long)Z] = 0.0;
  VUP[(long)W] = 1.0;
  COP[(long)X] = VRP[(long)X] + 10;
  COP[(long)Y] = VRP[(long)Y] + 10;
  COP[(long)Z] = VRP[(long)Z] + 10;
  COP[(long)W] = 1.0;
  UMIN = -50.0;
  VMIN = -50.0;
  UMAX = 50.0;
  VMAX = 50.0;
  F = -1000.0;
  B = 1000.0;
  SETUP();
}

/********************************************************************/

static void
move3(double x1, double y1, double z1)
{
  /* PURPOSE : MOVES TO A THREE DIMENSIONAL COORD WITHOUT DRAWWINF */
  POINT p3;

  p3[(long)X] = x1;
  p3[(long)Y] = y1;
  p3[(long)Z] = z1;
  p3[(long)W] = 1.0;
  MULTIPLY_VEC(THEMATRIX, p3, P1);
}

/********************************************************************/

static void
draw3(double x2, double y2, double z2)
{
  /* PURPOSE : DRAWS FROM THE CURRENT POINT TO X2,Y2,Z2 AND MAKES */
  /*           THE SECOND POINT THE CURRENT POINT */
  POINT p3;

  p3[(long)X] = x2;
  p3[(long)Y] = y2;
  p3[(long)Z] = z2;
  p3[(long)W] = 1.0;
  MULTIPLY_VEC(THEMATRIX, p3, P2);
  CLIPPING(P1, P2);
  MAKE(P1, P2);
}

/********************************************************************/

Static Void READY_DRAW()
{
  /* PURPOSE : TO PREPARE THE MATRIX WITH THE VIEW MATRIX */
  MULTIPLY_MAT(TRANSMATRIX, VIEWMATRIX, THEMATRIX);
}

/********************************************************************/

Static Void TRANSLATE(TX, TY, TZ)
double TX, TY, TZ;
{
  MATRIX T;

  FILL_MATRIX(T, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
              -TX, -TY, -TZ, 1.0);
  MULTIPLY_MAT(TRANSMATRIX, T, TRANSMATRIX);
}

/********************************************************************/

Static Void SCALE(SX, SY, SZ)
double SX, SY, SZ;
{
  MATRIX S;

  FILL_MATRIX(S, SX, 0.0, 0.0, 0.0, 0.0, SY, 0.0, 0.0, 0.0, 0.0, SZ, 0.0, 0.0,
              0.0, 0.0, 1.0);
  MULTIPLY_MAT(TRANSMATRIX, S, TRANSMATRIX);
}

/********************************************************************/

static void
rotate_x(double theta)
{
  MATRIX R;
  double SINT, COST;

  SINT = sin(theta);
  COST = cos(theta);
  FILL_MATRIX(R, 1.0, 0.0, 0.0, 0.0, 0.0, COST, SINT, 0.0, 0.0, -SINT, COST,
              0.0, 0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(TRANSMATRIX, R, TRANSMATRIX);
}

/********************************************************************/

static void
rotate_y(double theta)
{
  MATRIX R;
  double SINT, COST;

  SINT = sin(theta);
  COST = cos(theta);
  FILL_MATRIX(R, COST, 0.0, -SINT, 0.0, 0.0, 1.0, 0.0, 0.0, SINT, 0.0, COST,
              0.0, 0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(TRANSMATRIX, R, TRANSMATRIX);
}

/********************************************************************/

static void
rotate_z(theta)
{
  MATRIX R;
  double SINT, COST;

  SINT = sin(theta);
  COST = cos(theta);
  FILL_MATRIX(R, COST, SINT, 0.0, 0.0, -SINT, COST, 0.0, 0.0, 0.0, 0.0, 1.0,
              0.0, 0.0, 0.0, 0.0, 1.0);
  MULTIPLY_MAT(TRANSMATRIX, R, TRANSMATRIX);
}

/********************************************************************/

Static Void CUBE()
{
  static double XX[8] = {
      4.0, 4.0, -4.0, -4.0, 4.0, 4.0, -4.0, -4.0};

  static double YY[8] = {
      4.0, -4.0, -4.0, 4.0, 4.0, -4.0, -4.0, 4.0};

  static double ZZ[8] = {
      4.0, 4.0, 4.0, 4.0, -4.0, -4.0, -4.0, -4.0};

  static uchar LV[12][2] = {
      {1, 2},
      {2, 3},
      {3, 4},
      {4, 1},
      {5, 6},
      {6, 7},
      {7, 8},
      {8, 5},
      {1, 5},
      {2, 6},
      {4, 8},
      {3, 7}};

  static uchar LF[6][4] = {
      {1, 2, 3, 4},
      {1, 10, 5, 9},
      {5, 6, 7, 8},
      {7, 12, 3, 11},
      {4, 9, 8, 11},
      {2, 10, 6, 12}};

  uchar I, J, INOV, INOL, INOF;
  POINT P;
  double DD;

  for (I = 1; I <= 8; I++)
  {
    INOV = I + OBJ.NOV;
    P[(long)X] = XX[I - 1];
    P[(long)Y] = YY[I - 1];
    P[(long)Z] = ZZ[I - 1];
    P[(long)W] = 1.0;
    MULTIPLY_VEC(THEMATRIX, P, OBJ.POINTS[INOV - 1]);
    DD = OBJ.POINTS[INOV - 1][(long)Z] / D;
    OBJ.XP[INOV - 1] = OBJ.POINTS[INOV - 1][(long)X] * DD;
    OBJ.YP[INOV - 1] = OBJ.POINTS[INOV - 1][(long)Y] * DD;
  }
  for (I = 1; I <= 12; I++)
  {
    INOL = I + OBJ.NOL;
    OBJ.LINV[INOL - 1][0] = LV[I - 1][0] + OBJ.NOV;
    OBJ.LINV[INOL - 1][1] = LV[I - 1][1] + OBJ.NOV;
  }
  for (I = 1; I <= 6; I++)
  {
    INOF = I + OBJ.NOF;
    OBJ.INDEXF[INOF - 1] = 4; /* NUMBER OF SIDES/FACETS */
    for (J = 0; J <= 3; J++)
      OBJ.LINF[INOF - 1][J] = LF[I - 1][J] + OBJ.NOL;
  }
  OBJ.NOV += 8;
  OBJ.NOL += 12;
  OBJ.NOF += 6; /* WITH */
}

/********************************************************************/

static void
pyramid(double x)
{
  /* THE FIXED POINT IS THE LOWER LEFT FRONT OF THE OBJECT */
  double y;

  y = x * 2;
  move3(0.0, 0.0, 0.0);
  draw3(y, 0.0, 0.0);
  draw3(y, 0.0, -y);
  draw3(0.0, 0.0, -y);
  draw3(0.0, 0.0, 0.0);
  draw3(x, y, -x);
  draw3(0.0, 0.0, -y);
  move3(y, 0.0, 0.0);
  draw3(x, y, -x);
  draw3(y, 0.0, -y);
}

/********************************************************************/

Static Void SETDRAW()
{
  OBJ.NOV = 0;
  OBJ.NOL = 0;
  OBJ.NOF = 0;
  RESET_TRANS();
  SCALE(2.0, 2.0, 2.0);
  TRANSLATE(1.0, 1.0, 1.0);
  READY_DRAW();
  CUBE();

  RESET_TRANS();
  TRANSLATE(-20.0, -30.0, -5.0);
  READY_DRAW();
  CUBE();
}

/***************    MAIN  PROGRAM   *********************************/

main(argc, argv) int argc;
Char *argv[];
{
  /*PASCAL_MAIN(argc, argv);*/
  NIB = NULL;
  /*ASSIGN (NIB,'A:FRED.XX');*/
  /*REWRITE (NIB);*/
  GET_VIEW();
  initialize();
  /*  GRAPHCOLORMODE(); */
  /* p2c: test3.pas, line 1112:
   * Warning: Symbol 'GRAPHCOLORMODE' is not defined [221] */
  SETUP();
  SETDRAW();
  hidden();
  /*CLOSE (NIB);*/
  if (NIB != NULL)
    fclose(NIB);
  exit(EXIT_SUCCESS);
}

/* End. */
