
/*Distance Subroutine - Spherical Earth*/
int delaz5_(thei, alei, thsi, alsi, delt, deltdg, deltkm, 
            azes, azesdg, azse, azsedg, i)

float *thei, *alei, *thsi, *alsi, *delt, *deltdg, *deltkm, *azes;
float *azesdg, *azse, *azsedg;
int *i;
{
    double d__1, d__2, d__3;

    double tan(), atan(), sin(), cos(), acos(), atan2(), sqrt(), asin();

    /* Local variables */
    static double a, b, c, d, e, g, h;
    static float c1, c3, c4, c5, c6;
    static double ak, ap, bp, cp, dp, ep, gp, hp;
    static float aaa, ale;
    static double akp;
    static float als, the, ths;

    if (*i <= 0) {
	goto L50;
    } else {
	goto L51;
    }
/*     IF  COORDINATES ARE GEOGRAPH DEG I=0 */
/*     IF COORDINATES ARE GEOCENT RADIAN  I=1 */
L50:
    the = *thei * (float).01745329252;
    ale = *alei * (float).01745329252;
    ths = *thsi * (float).01745329252;
    als = *alsi * (float).01745329252;
    aaa = tan(the) * (float).9931177;
    the = atan(aaa);
    aaa = tan(ths) * (float).9931177;
    ths = atan(aaa);
    goto L32;
L51:
    the = *thei;
    ale = *alei;
    ths = *thsi;
    als = *alsi;
L32:
    c = sin(the);
    ak = -(double)cos(the);
    d = sin(ale);
    e = -(double)cos(ale);
    a = ak * e;
    b = -ak * d;
    g = -c * e;
    h = c * d;
    cp = sin(ths);
    akp = -(double)cos(ths);
    dp = sin(als);
    ep = -(double)cos(als);
    ap = akp * ep;
    bp = -akp * dp;
    gp = -cp * ep;
    hp = cp * dp;
    c1 = a * ap + b * bp + c * cp;
    if (c1 - (float).94 >= (float)0.) {
	goto L31;
    } else {
	goto L30;
    }
L30:
    if (c1 + (float).94 <= (float)0.) {
	goto L28;
    } else {
	goto L29;
    }
L29:
    *delt = acos(c1);
L33:
    *deltkm = *delt * (float)6371.;
/* Computing 2nd power */
    d__1 = ap - d;
/* Computing 2nd power */
    d__2 = bp - e;
/* Computing 2nd power */
    d__3 = cp;
    c3 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - (float)2.;
/* Computing 2nd power */
    d__1 = ap - g;
/* Computing 2nd power */
    d__2 = bp - h;
/* Computing 2nd power */
    d__3 = cp - ak;
    c4 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - (float)2.;
/* Computing 2nd power */
    d__1 = a - dp;
/* Computing 2nd power */
    d__2 = b - ep;
/* Computing 2nd power */
    d__3 = c;
    c5 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - (float)2.;
/* Computing 2nd power */
    d__1 = a - gp;
/* Computing 2nd power */
    d__2 = b - hp;
/* Computing 2nd power */
    d__3 = c - akp;
    c6 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - (float)2.;
    *deltdg = *delt * (float)57.29577951;
    *azes = atan2(c3, c4);
    if (*azes >= (float)0.) {
	goto L81;
    } else {
	goto L80;
    }
L80:
    *azes += (float)6.283185308;
L81:
    *azse = atan2(c5, c6);
    if (*azse >= (float)0.) {
	goto L71;
    } else {
	goto L70;
    }
L70:
    *azse += (float)6.283185308;
L71:
    *azesdg = *azes * (float)57.29577951;
    *azsedg = *azse * (float)57.29577951;
    return 0;
L31:
/* Computing 2nd power */
    d__1 = a - ap;
/* Computing 2nd power */
    d__2 = b - bp;
/* Computing 2nd power */
    d__3 = c - cp;
    c1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    c1 = sqrt(c1);
    c1 /= (float)2.;
    *delt = asin(c1);
    *delt *= (float)2.;
    goto L33;
L28:
/* Computing 2nd power */
    d__1 = a + ap;
/* Computing 2nd power */
    d__2 = b + bp;
/* Computing 2nd power */
    d__3 = c + cp;
    c1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    c1 = sqrt(c1);
    c1 /= (float)2.;
    *delt = acos(c1);
    *delt *= (float)2.;
    goto L33;
} /* delaz5_ */

