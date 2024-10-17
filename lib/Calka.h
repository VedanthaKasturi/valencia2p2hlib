#ifndef Calka_H
#define Calka_H
#include <math.h>

template <class F, class U, class V, class X, class Y, class W, class Q, class P, class E>
class Fix12345678
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;
    W _w;
    Q _q;
    P _p;
    E _e;
    
    
public:
    Fix12345678(F  f, U u,V v,X x,Y y, W w, Q q, P p, E e) :_f(f), _u(u), _v(v), _x(x),_y(y), _w(w), _q(q), _p(p), _e(e){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _u, _v,_x,_y,_w, _q, _p,_e);
    }
};

template <class F, class U, class V, class X, class Y, class W, class Q, class P, class E>
inline Fix12345678<F,U,V,X,Y,W,Q,P,E> fix12345678(F  f, U u, V v, X x, Y y, W w, Q q, P p, E e)
{
    return Fix12345678<F,U,V,X,Y,W,Q,P,E>(f,u,v,x,y,w,q,p,e);
}

template <class F, class U, class V, class X, class Y, class W, class Q, class P>
class Fix1234567
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;
    W _w;
    Q _q;
    P _p;
    
    
public:
    Fix1234567(F  f, U u,V v,X x,Y y, W w, Q q, P p) :_f(f), _u(u), _v(v), _x(x),_y(y), _w(w), _q(q), _p(p){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _u, _v,_x,_y,_w, _q, _p);
    }
};

template <class F, class U, class V, class X, class Y, class W, class Q, class P>
inline Fix1234567<F,U,V,X,Y,W,Q,P> fix1234567(F  f, U u, V v, X x, Y y, W w, Q q, P p)
{
    return Fix1234567<F,U,V,X,Y,W,Q,P>(f,u,v,x,y,w,q,p);
}

template <class F, class U, class V, class X, class Y, class W, class Q>
class Fix123456
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;
    W _w;
    Q _q;


public:
    Fix123456(F  f, U u,V v,X x,Y y, W w, Q q) :_f(f), _u(u), _v(v), _x(x),_y(y), _w(w), _q(q){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _u, _v,_x,_y,_w, _q);
    }
};

template <class F, class U, class V, class X, class Y, class W, class Q>
inline Fix123456<F,U,V,X,Y,W,Q> fix123456(F  f, U u, V v, X x, Y y, W w, Q q)
{
    return Fix123456<F,U,V,X,Y,W,Q>(f,u,v,x,y,w,q);
}


template <class F, class U, class V, class X, class Y, class W>
class Fix12345
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;
    W _w;


public:
    Fix12345(F  f, U u,V v,X x,Y y, W w) :_f(f), _u(u), _v(v), _x(x),_y(y), _w(w){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _u, _v,_x,_y,_w);
    }
};

template <class F, class U, class V, class X, class Y, class W>
inline Fix12345<F,U,V,X,Y,W> fix12345(F  f, U u, V v, X x, Y y, W w)
{
    return Fix12345<F,U,V,X,Y,W>(f,u,v,x,y,w);
}

template <class F, class U, class V, class X, class Y>
class Fix1345
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;
    
    
public:
    Fix1345(F  f, U u,V v,X x,Y y) :_f(f), _u(u), _v(v), _x(x),_y(y){}
    template <class Z>
    double operator()(Z z)
    { return _f(_u, z, _v,_x,_y);
    }
};

template <class F, class U, class V, class X, class Y>
inline Fix1345<F,U,V,X,Y> fix1345(F  f, U u, V v, X x, Y y)
{
    return Fix1345<F,U,V,X,Y>(f,u,v,x,y);
}


template <class F, class U, class V, class X, class Y>
class Fix1234
{ F _f;
    U _u;
    V _v;
    X _x;
    Y _y;


public:
    Fix1234(F  f, U u,V v,X x,Y y) :_f(f), _u(u), _v(v), _x(x),_y(y){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _u, _v,_x,_y);
    }
};

template <class F, class U, class V, class X, class Y>
inline Fix1234<F,U,V,X,Y> fix1234(F  f, U u, V v, X x, Y y)
{
    return Fix1234<F,U,V,X,Y>(f,u,v,x,y);
}





template <class F, class V, class X, class Y>
class Fix123
{ F _f;
    V _v;
    X _x;
    Y _y;


public:
    Fix123(F  f,V v,X x,Y y) :_f(f), _v(v), _x(x),_y(y){}
    template <class Z>
    double operator()(Z z)
    { return _f(z, _v,_x,_y);
    }
};

template <class F, class V, class X, class Y>
inline Fix123<F,V,X,Y> fix123(F  f, V v, X x, Y y)
{
    return Fix123<F,V,X,Y>(f,v,x,y);
}


template <class F, class X, class Y>
class Fix12
{ F _f;
    X _x;
    Y _y;
public:
    Fix12(F  f, X x,Y y) :_f(f),_x(x),_y(y){}
    template <class Z>
    double operator()(Z z)
    { return _f(z,_x,_y);
    }
};

template <class F, class X, class Y>
inline Fix12<F,X,Y> fix12(F  f, X x, Y y)
{
    return Fix12<F,X,Y>(f,x,y);
}


template <class F, class Y>
class Fix2nd
{ F _f;
       Y _y;
       public:
       Fix2nd(F  f, Y y) :_f(f),_y(y){}
     template <class X>
     double operator()(X x)
     { return _f(x,_y);
          }
     } ;
template <class F, class Y>
inline Fix2nd<F,Y> fix2nd(F  f, Y y)
 {
       return Fix2nd<F,Y>(f,y);
    }

template <class F, class X, class Y, class Z, class T, class U, class W, class A, class B, class C>
class Fix9
{ F _f;
    X _x;
    Y _y;
    Z _z;
    T _t;
    U _u;
    W _w;
    A _a;
    B _b;
    C _c;
    public:
    Fix9(F  f, X x, Y y, Z z, T t, U u, W w, A a, B b, C c) :_f(f),_x(x),_y(y),_z(z),_t(t), _u(u), _w(w), _a(a), _b(b), _c(c){}
    double operator()(double xx)
                    { return _f(xx,_x,_y,_z,_t,_u, _w, _a, _b, _c);
                    }
} ;


template <class F, class X, class Y, class Z, class T, class U, class W, class A, class B, class C>
inline Fix9<F,X,Y,Z,T,U,W,A,B,C> fix9(F  f, X x, Y y, Z z, T t, U u, W w, A a, B b, C c)
{
    return Fix9<F,X,Y,Z,T,U,W,A,B,C>(f,x,y,z,t,u,w,a,b,c);
}

template <class F, class X, class Y, class Z, class T, class U, class W, class A, class B>
class Fix8_5_special
{ F _f;
	X _x;
	Y _y;
	Z _z;
	T _t;
	U _u;
	W _w;
	A _a;
	B _b;
public:
	Fix8_5_special(F  f, X x, Y y, Z z, T t, U u, W w, A a, B b) :_f(f),_x(x),_y(y),_z(z),_t(t), _u(u), _w(w), _a(a), _b(b) {}
	void operator()(double xx, double& w1, double& w2, double& w3, double& w4, double& w5)
	{ return _f(xx, w1, w2, w3, w4, w5, _x,_y,_z,_t,_u, _w, _a, _b);
	}
} ;


template <class F, class X, class Y, class Z, class T, class U, class W, class A, class B>
inline Fix8_5_special<F,X,Y,Z,T,U,W,A,B> fix8_5_special(F  f, X x, Y y, Z z, T t, U u, W w, A a, B b)
{
	return Fix8_5_special<F,X,Y,Z,T,U,W,A,B>(f,x,y,z,t,u,w,a,b);
}


template <class _T,class F>
class Makefun
{ _T & _t;
  F _f;
 public:
 Makefun(_T& t, F f):_t(t),_f(f){}
  template <class X>
 double operator()(X x)
 {return (_t.*_f)(x);
 }


template <class X,class Y>
 double operator()(X x, Y y)
 {return (_t.*_f)(x,y);
 }

template <class X,class Y,class Z>
double operator()(X x, Y y,Z z)
 {return (_t.*_f)(x,y,z);
 }


template <class X,class Y,class Z, class T >
double operator()(X x, Y y,Z z,T t)
 {return (_t.*_f)(x,y,z,t);
 }

template <class X,class Y,class Z, class T, class U >
double operator()(X x, Y y,Z z,T t,U u)
 {return (_t.*_f)(x,y,z,t,u);
 }

template <class X,class Y,class Z, class T, class U,class V >
double operator()(X x, Y y,Z z,T t,U u,V v)
 {return (_t.*_f)(x,y,z,t,u,v);
 }

template <class X,class Y,class Z, class T, class U,class V, class Q >
double operator()(X x, Y y,Z z,T t,U u,V v, Q q)
{return (_t.*_f)(x,y,z,t,u,v,q);
    }

template <class X,class Y,class Z, class T, class U,class V, class Q, class P >
double operator()(X x, Y y,Z z,T t,U u,V v, Q q, P p)
{return (_t.*_f)(x,y,z,t,u,v,q,p);
}

};


template <class T,class F>
Makefun<T,F>  makefun(T& t, F f)
 {return Makefun<T,F>(t,f);
 }


template <class F>
void calg20_special(F f,double x1,double x2, double& wyn1, double& wyn2, double& wyn3, double& wyn4, double& wyn5, int n=3)
{
	static const double w20[10]={
		0.017614007139152118311861962351853,
		0.040601429800386941331039952274932,
		0.062672048334109063569506535187042,
		0.083276741576704748724758143222046,
		0.101930119817240435036750135480350,
		0.118194531961518417312377377711382,
		0.131688638449176626898494499748163,
		0.142096109318382051329298325067165,
		0.149172986472603746787828737001969,
		0.152753387130725850698084331955098
	};
	
	static const double g20[10]={
		0.993128599185094924786122388471320,
		0.963971927277913791267666131197277,
		0.912234428251325905867752441203298,
		0.839116971822218823394529061701521,
		0.746331906460150792614305070355642,
		0.636053680726515025452836696226286,
		0.510867001950827098004364050955251,
		0.373706088715419560672548177024927,
		0.227785851141645078080496195368575,
		0.076526521133497333754640409398838
	};
	
	//if (n>1) return calg20(f,x1,x1+(x2-x1)*(n/2)/n,n/2)+calg20(f,x1+(x2-x1)*(n/2)/n,x2,n-n/2);
	
	double wyn11,wyn12,wyn21,wyn22,wyn31,wyn32,wyn41,wyn42,wyn51,wyn52;
	wyn11 = wyn12 = wyn21 = wyn22 = wyn31 = wyn32 = wyn41 = wyn42 = wyn51 = wyn52= 0;
	
	double A=(x2-x1)/2/n;
	for (int i=0;i<n;i++)
	{
		double B=(x1+x2)/2+(2*i+1-n)*A;
		for (int j=0; j<10; j++) {
			f(B+A*g20[j], wyn11, wyn21, wyn31, wyn41, wyn51);
			f(B-A*g20[j], wyn12, wyn22, wyn32, wyn42, wyn52);
			wyn1+=(wyn11+wyn12)*w20[j];
			wyn2+=(wyn21+wyn22)*w20[j];
			wyn3+=(wyn31+wyn32)*w20[j];
			wyn4+=(wyn41+wyn42)*w20[j];
			wyn5+=(wyn51+wyn52)*w20[j];
			
		}
	}
	wyn1*=A;
	wyn2*=A;
	wyn3*=A;
	wyn4*=A;
	wyn5*=A;
	
}


template <class F> double calg20(F f,double x1,double x2, int n=3)
{
	static const double w20[10]={
	  0.017614007139152118311861962351853,
	  0.040601429800386941331039952274932,
	  0.062672048334109063569506535187042,
	  0.083276741576704748724758143222046,
	  0.101930119817240435036750135480350,
	  0.118194531961518417312377377711382,
	  0.131688638449176626898494499748163,
	  0.142096109318382051329298325067165,
	  0.149172986472603746787828737001969,
	  0.152753387130725850698084331955098
	};

	static const double g20[10]={
	  0.993128599185094924786122388471320,
	  0.963971927277913791267666131197277,
	  0.912234428251325905867752441203298,
	  0.839116971822218823394529061701521,
	  0.746331906460150792614305070355642,
	  0.636053680726515025452836696226286,
	  0.510867001950827098004364050955251,
	  0.373706088715419560672548177024927,
	  0.227785851141645078080496195368575,
	  0.076526521133497333754640409398838
	};

	//if (n>1) return calg20(f,x1,x1+(x2-x1)*(n/2)/n,n/2)+calg20(f,x1+(x2-x1)*(n/2)/n,x2,n-n/2);

	double wynik=0;
	double A=(x2-x1)/2/n;
	for (int i=0;i<n;i++)
	{
		double B=(x1+x2)/2+(2*i+1-n)*A;
		for (int j=0; j<10; j++)
			wynik+=(f(B+A*g20[j])+f(B-A*g20[j]))*w20[j];
	}
	return wynik*A;
}

template <class F> double calg8(F f,double x1,double x2)
{
    static const double w8[4]={
        0.1012285362,
        0.2223810344,
        0.3137066458,
        0.3626837833
    };
    
    static const double g8[4]={
        -.9602898564,
        -.7966664774,
        -.5255324099,
        -.1834346424
    };
    
    //if (n>1) return calg20(f,x1,x1+(x2-x1)*(n/2)/n,n/2)+calg20(f,x1+(x2-x1)*(n/2)/n,x2,n-n/2);
    
    double wynik=0;
    double A=(x2-x1)/2;
    double B=(x1+x2)/2.;
    for (int j=0; j<4; j++)
        wynik+=(f(B+A*g8[j])+f(B-A*g8[j]))*w8[j];
    
    return wynik*A;
}

template <class F> double supercalka(F f, double x1, double x2, double precision, double renorm)
 {
    double results[2] = { calg20(f,x1,x2,1), calg20(f,x1,x1+(x2-x1)/2.0,1) + calg20(f,x1+(x2-x1)/2.0,x2,1)};

    if (results[0]==0 && results [1]==0) {
        return 0;}
    else {
        if ( 2.0*fabs(results[0]-results[1])< precision*fabs(results[0]+ results[1]) ) {
            return results[1];
        }
        else {
            return supercalka(f,x1,x1+(x2-x1)/2.0 ,precision*renorm, renorm) + supercalka(f,x1+(x2-x1)/2.0,x2,precision*renorm, renorm);
        }
    }
}





#endif
