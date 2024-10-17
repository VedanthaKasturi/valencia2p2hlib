#ifndef MatrixDir_h_
#define MatrixDir_h_
#include <cmath>
#include <complex>
#include <iostream>

static unsigned short int dim=4;

class Spinor;

//4x4 Complex matrix class
class DM
{
        friend class Spinor;
public:
    DM(unsigned short i=0)
    {
        
        //default 4x4 matrix of zeroes
        for (unsigned short a=0; a<4; a++)
            for (unsigned short b=0; b<4; b++)
            {
                matrix[a][b]= std::complex<double>(0,0);
                matrix[a][b]= std::complex<double>(0,0);
            }
        //gamma0
        if (i==0)
        {
            matrix[0][0]= std::complex<double>(1,0);
            matrix[1][1]= std::complex<double>(1,0);
            matrix[2][2]= std::complex<double>(-1,0);
            matrix[3][3]= std::complex<double>(-1,0);
        }
        //gamma1
        if (i==1)
        {
            matrix[0][3]= std::complex<double>(1,0);
            matrix[1][2]= std::complex<double>(1,0);
            matrix[3][0]= std::complex<double>(-1,0);
            matrix[2][1]= std::complex<double>(-1,0);
        }
        //gamma2
        if (i==2)
        {
            matrix[0][3]= std::complex<double>(0,-1);
            matrix[1][2]= std::complex<double>(0,1);
            matrix[3][0]= std::complex<double>(0,-1);
            matrix[2][1]= std::complex<double>(0, 1);
        }
        //gamma3
        if (i==3)
        {
            matrix[0][2]= std::complex<double>(1,0);
            matrix[1][3]= std::complex<double>(-1,0);
            matrix[3][1]= std::complex<double>(1,0);
            matrix[2][0]= std::complex<double>(-1,0);
        }
        // identity 4x4
        if (i==4)
        {
            matrix[0][0]= std::complex<double>(1,0);
            matrix[1][1]= std::complex<double>(1,0);
            matrix[2][2]= std::complex<double>(1,0);
            matrix[3][3]= std::complex<double>(1,0);
        }
        //gamma5
        if (i==5)
        {
            matrix[0][2]= std::complex<double>(1,0);
            matrix[1][3]= std::complex<double>(1,0);
            matrix[2][0]= std::complex<double>(1,0);
            matrix[3][1]= std::complex<double>(1,0);
        }
        //wrt();
    }
	DM(const double &a, const double &b, const double &c, const double &d);
	DM(const double &a, const double &b, const double &c, const double &d, const double & e);
	~DM(){;};
	
	//bunch of overloaded operators
	inline DM& operator=(const DM& mat)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]=mat.matrix[a][b];
	
		}
		
		return *this;
	}
	
	inline friend DM operator+ (const DM  &m1,const DM &m2)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]+m2.matrix[a][b];
	
		}
		return result;

	}

	//add number to the diagonal elements. Stupid implementation,
	//for some reason this bullshit it is 20%
	//faster, than accessing just the diagonals directly
	//blame the processor architecture. Fixing requires Assembler!
	inline friend DM operator+ (const DM  &m1,const double &d)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]+diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator+ (const double &d, const DM  &m1)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]+diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator+ (const DM  &m1,const std::complex<double> &d)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]+diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator+ (const std::complex<double> &d, const DM  &m1)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]+diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline DM& operator+= (const DM& m)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]+=m.matrix[a][b];
	
		}
		return *this;

	}

	//add number to the diagonal elements. Stupid implementation,
	//for some reason this bullshit it is 20%
	//faster, than accessing just the diagonals directly
	//blame the processor architecture. Fixing requires Assembler!
	inline DM& operator+= (const double &d)
	{
		DM identi(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]+=identi.matrix[a][b]*d;
	
		} //Why The Hell is this the fastest implementation??????????
		//for (unsigned short a=0; a<dim; a++)
		//	matrix[a][a]+=d;//->how it should look!
		return *this;
	}

	inline DM& operator-= (const double &d)
	{
		for (unsigned short a=0; a<dim; a++)
			matrix[a][a]-=d;
		return *this;
	}

	inline DM& operator+= (const std::complex<double> &d)
	{
		for (unsigned short a=0; a<dim; a++)
			matrix[a][a]+=d;
		return *this;
	}

	inline DM& operator-= (const std::complex<double> &d)
	{
		for (unsigned short a=0; a<dim; a++)
			matrix[a][a]-=d;
		return *this;
	}


	inline DM& operator-= (const DM& m)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]-=m.matrix[a][b];
	
		}
		return *this;

	}


	inline friend DM operator- (const DM  &m1,const DM& m2)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]-m2.matrix[a][b];
	
		}
		return result;

	}

	inline friend DM operator- (const DM  &m1,const double &d)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]-diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator- (const double &d, const DM  &m1)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]-diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator- (const DM  &m1,const std::complex<double> &d)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]-diag.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator- (const std::complex<double> &d, const DM  &m1)
	{
		DM result(11);
		DM diag(4);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]-diag.matrix[a][b]*d;
	
		}
		return result;

	}


	inline friend DM operator* (const DM  &m1,const DM &m2)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		for (unsigned short c=0; c<dim; c++)
		{
			result.matrix[a][b]+=m1.matrix[a][c]*m2.matrix[c][b];
	
		}
		return result;

	}
	
	//multiplication A=A*B
	inline DM& operator*= (const DM&  m)
	{
//MOJA ZMIANA dim->4
		std::complex<double> tmp[4][4];
			
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			tmp[a][b]=std::complex<double>(0,0);
			for (unsigned short c=0; c<dim; c++)
			{
				tmp[a][b]+=matrix[a][c]*m.matrix[c][b];
	
			}
			
		}

		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		matrix[a][b]=tmp[a][b];
			
		return *this;

	}

	//multiplication A=B*A
	inline DM& operator&= (const DM&  m)
	{
        //MOJA ZMIANA dim->4
		std::complex<double> tmp[4][4];
			
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			tmp[a][b]=std::complex<double>(0,0);
			for (unsigned short c=0; c<dim; c++)
			{
				tmp[a][b]+=matrix[c][b]*m.matrix[a][c];
	
			}
			
		}

		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		matrix[a][b]=tmp[a][b];
			
		return *this;
	}


	inline friend DM operator* (const DM  &m1,const double &d)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]*d;
	
		}
		return result;

	}
	
	inline DM& operator*= (const double& d)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]*=d;
	
		}
		return *this;

	}

	inline friend DM operator* (const double &d,const DM &m1)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]*d;
	
		}
		return result;

	}
 
	inline friend DM operator* (const DM  &m1,const std::complex<double> &d)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]*d;
	
		}
		return result;

	}

	inline DM& operator*= (const std::complex<double>& d)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]*=d;
	
		}
		return *this;

	}

	inline friend DM operator* (const std::complex<double> &d,const DM &m1)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]*d;
	
		}
		return result;

	}

	inline friend DM operator/ (const DM  &m1,const double &d)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]/d;
	
		}
		return result;

	}

	inline DM& operator/= (const double &d)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]/=d;
	
		}
		return *this;

	}

	
 
	inline friend DM operator/ (const DM  &m1,const std::complex<double> &d)
	{
		DM result(11);
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			result.matrix[a][b]=m1.matrix[a][b]/d;
	
		}
		return result;

	}

	inline DM& operator/= (const std::complex<double> &d)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
			matrix[a][b]/=d;
	
		}
		return *this;

	}

	void hermit();
	void transp();


	inline std::complex<double> operator() (unsigned short row, unsigned short col) const
	{
//		if ((row>dim)||(col>dim))
//			return std::complex<double>(0,0);
		return matrix[row%dim][col%dim];
	}

	inline std::complex<double>& operator() (unsigned short row, unsigned short col)
	{
		
		//if ((row>4)||(col>4))
		//	return (*this).matrix[0][0];
		//else
			return (*this).matrix[row%dim][col%dim];
	}
  
  	//real part of trace
	inline double Trace()
	{
		return real(matrix[0][0]+matrix[1][1]+matrix[2][2]+matrix[3][3]);

	}
	//complex trace
	inline std::complex<double> CTrace()
	{
		return matrix[0][0]+matrix[1][1]+matrix[2][2]+matrix[3][3];

	}
	
	inline std::complex<double> BTrace(const DM & m2)
	{
		return matrix[0][0]*m2.matrix[0][0] + matrix[0][1]*m2.matrix[1][0]
		     //+ matrix[0][2]*m2.matrix[2][0] + matrix[0][3]*m2.matrix[3][0]
		     + matrix[1][0]*m2.matrix[0][1] + matrix[1][1]*m2.matrix[1][1]
		    // + matrix[1][2]*m2.matrix[2][1] + matrix[1][3]*m2.matrix[3][1]
		     + matrix[2][0]*m2.matrix[0][2] + matrix[2][1]*m2.matrix[1][2]
		    // + matrix[2][2]*m2.matrix[2][2] + matrix[2][3]*m2.matrix[3][2]
		     + matrix[3][0]*m2.matrix[0][3] + matrix[3][1]*m2.matrix[1][3];
		    // + matrix[3][2]*m2.matrix[2][3] + matrix[3][3]*m2.matrix[3][3];

	}
//inactive! writing down the matrix contents
    void wrt();
	//friend ostream& operator<<(ostream& os, const DM& mat);
private:

	std::complex<double> matrix[4][4];
};

/********
 
 SPINOR
 
 ********/

class Spinor{

public:
    std::complex<double> matrix[4];
    
    Spinor(const std::complex<double> &a, const std::complex<double> &b, const std::complex<double> &c, const std::complex<double> &d){
        matrix[0] = a;
        matrix[1]=b;
        matrix[2]=c;
        matrix[3]=d;
        
    }
    inline std::complex<double> operator() (unsigned short j) const
    {
            return matrix[j];
    }
    
    inline friend std::complex<double> operator* (const Spinor  &m1,const Spinor &m2)
    {
        std::complex<double> result = m1.matrix[0]*m2.matrix[0]+ m1.matrix[1]*m2.matrix[1]+ m1.matrix[2]*m2.matrix[2]+ m1.matrix[3]*m2.matrix[3];
        return result;
        
    }
    
    inline friend Spinor operator* (const Spinor  &m1,const DM &m2)
    {
        Spinor result(0,0,0,0);
        for (unsigned short a=0; a<dim; a++)
            for (unsigned short b=0; b<dim; b++)
            {
                result.matrix[a]+=m1.matrix[b]*m2.matrix[b][a];
                
            }
        return result;
        
    }
    
    

    
};


//generic template of a Lorentz 4-vector
template< class T>
class D4V
{
	public:
    
    D4V(){}
    
	D4V(T t0,T t1, T t2, T t3): t0_(t0), t1_(t1), t2_(t2), t3_(t3)
	{}
	
	inline D4V& operator=(const D4V& mat)
	{
		t0_=mat.t0_;
		t1_=mat.t1_;
		t2_=mat.t2_;
		t3_=mat.t3_;
		return *this;
	}
	inline T operator() (unsigned short j) const
	{
        
		if (j==0)
			return t0_;
		else if (j==1)
			return t1_;
		else if (j==2)
			return t2_;
		else if (j==3)
			return t3_;
		else
		{
			double failsafe=0;
			return t1_*failsafe;
		}
	}

	inline T& operator() (unsigned short j)
	{
		if (j==0)
			return (*this).t0_;
		else if (j==1)
			return (*this).t1_;
		else if (j==2)
			return (*this).t2_;
		else if (j==3)
			return (*this).t3_;
		else
		{
			
			return (*this).t1_;
		} 
	}


	//again, bunch of overloaded operators
	inline friend D4V operator+ (const D4V<T> &v1,const D4V<T> &v2)
	{
		D4V<T> result(v1.t0_ + v2.t0_, v1.t1_ + v2.t1_, v1.t2_ + v2.t2_, v1.t3_ + v2.t3_);
		return result;

	}

	
	inline D4V& operator+= (const D4V<T> &v1)
	{
		t0_+=v1.t0_;
		t1_+=v1.t1_;
		t2_+=v1.t2_;
		t3_+=v1.t3_;
		return *this;
	}

	inline friend D4V operator- (const D4V<T> &v1,const D4V<T> &v2)
	{
		D4V<T> result(v1.t0_ - v2.t0_, v1.t1_ - v2.t1_, v1.t2_ - v2.t2_, v1.t3_ - v2.t3_);
		return result;

	}

	inline D4V& operator-= (const D4V<T> &v1)
	{
		t0_-=v1.t0_;
		t1_-=v1.t1_;
		t2_-=v1.t2_;
		t3_-=v1.t3_;
		return *this;
	}

	inline friend D4V<DM> operator* (const D4V<T> &v1, const  DM &d)
	{
		D4V<DM> result(v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
		
		return result;

	}
	
	inline friend D4V<DM> operator* (const  DM &d, const D4V<T> &v1)
	{
		D4V<DM> result(d*v1.t0_, d*v1.t1_, d*v1.t2_, d*v1.t3_);
		
		return result;

	}
    
   //DODANA
    inline friend D4V<Spinor> operator* (const  Spinor &d, const D4V<T> &v1)
    {
        D4V<Spinor> result(d*v1.t0_, d*v1.t1_, d*v1.t2_, d*v1.t3_);
        
        return result;
        
    }
	
	inline friend D4V<T> operator* (const D4V<T> &v1, const std::complex<double> &d)
	{
		D4V<T> result(v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
		
		return result;

	}

		
	inline friend D4V<T> operator* (const D4V<T> &v1, const double &d)
	{
		D4V<T> result(v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
		
		return result;

	}
    
    inline friend D4V<T> operator* (const double &d, const D4V<T> &v1)
    {
        D4V<T> result(v1.t0_*d, v1.t1_*d, v1.t2_*d, v1.t3_*d);
        
        return result;
        
    }
	//vector-vector covariant multiplication. notice metric= (+,-,-,-)
	inline friend T operator* (const D4V<T> &v1,const D4V<double> &v2)
	{
		T result=v1.t0_*v2.t0_-v1.t1_*v2.t1_-v1.t2_*v2.t2_-v1.t3_*v2.t3_;
		return result;

	}

	inline friend T operator* (D4V<T> v1,D4V< std::complex<double> > v2)
	{
		T result=v1.t0_*v2.t0_-v1.t1_*v2.t1_-v1.t2_*v2.t2_-v1.t3_*v2.t3_;
		return result;

	}

	inline friend DM operator* (D4V<T> v1,D4V< DM > v2)
	{
		DM result=v1.t0_*v2.t0_-v1.t1_*v2.t1_-v1.t2_*v2.t2_-v1.t3_*v2.t3_;
		return result;

	}

	inline D4V& operator*= (const double &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}
	
	inline D4V& operator*= (const std::complex<double> &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}

	inline D4V& operator*= (const DM &d)
	{
		t0_*=d;
		t1_*=d;
		t2_*=d;
		t3_*=d;
		
		return *this;

	}
//left-handed multiplication by a matrix
	inline D4V& operator&= (const DM &d)
	{
		t0_=d*t0_;
		t1_=d*t1_;
		t2_=d*t2_;
		t3_=d*t3_;
		
		return *this;

	}


	inline void hermit()
	{
		(*this)(0).hermit();
		(*this)(1).hermit();
		(*this)(2).hermit();
		(*this)(3).hermit();
	}
	
	//returns the trace as a vector
	inline D4V< std::complex<double> > CT(const D4V<DM> &D )
	{
		D4V< std::complex<double> > result(std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0));
		
		for (unsigned short a=0; a<dim; a++)
				result(a)=D(a).CTrace();
		return result;
	}
    
	
	//private:
		T t0_;
		T t1_;
		T t2_;
		T t3_;
};

//generic 4x4 Lorentz tensor

template< class T >
class D4T
{
	public:
    
        D4T(){}
		//constructor as a tensor product of two vectors
		D4T(const D4V<T> &t,const D4V<T> &u):   t0_(t(0)*u(0), t(0)*u(1), t(0)*u(2), t(0)*u(3)),
							t1_(t(1)*u(0), t(1)*u(1), t(1)*u(2), t(1)*u(3)),
							t2_(t(2)*u(0), t(2)*u(1), t(2)*u(2), t(2)*u(3)),
							t3_(t(3)*u(0), t(3)*u(1), t(3)*u(2), t(3)*u(3)){}
		//element-by element constructor
		D4T(const T &t00,const T &t01, const T &t02, const T &t03,
		    const T &t10,const T &t11, const T &t12, const T &t13,
		    const T &t20,const T &t21, const T &t22, const T &t23,
		    const T &t30,const T &t31, const T &t32, const T &t33):
				t0_( t00, t01, t02, t03),
				t1_( t10, t11, t12, t13),
				t2_( t20, t21, t22, t23),
				t3_( t30, t31, t32, t33) {}
	

	//overloaded operators	
	inline D4T& operator=(const D4T& t)
	{
		t0_=t.t0_;
		t1_=t.t1_;
		t2_=t.t2_;
		t3_=t.t3_;
		return *this;
	}

	inline T operator() (unsigned short i, unsigned short j) const
	{
		if(i==0)
			return t0_(j);
		else if (i==1)
			return t1_(j);
		else if (i==2)
			return t2_(j);
		else if (i==3)
			return t3_(j);
		else
		{
			double failsafe=0;
			return t1_(0)*failsafe;
		}
	}

	inline T& operator() (unsigned short i, unsigned short j)
	{
		if(i==0)
			return (*this).t0_(j);
		else if (i==1)
			return (*this).t1_(j);
		else if (i==2)
			return (*this).t2_(j);
		else if (i==3)
			return (*this).t3_(j);
		else
		{
			return (*this).t1_(0);
		}
	}


	inline D4V<T> operator() (unsigned short i)
	{
		if(i==0)
			return t0_;
		else if (i==1)
			return t1_;
		else if (i==2)
			return t2_;
		else if (i==3)
			return t3_;
		else
		{
			double failsafe=0;
			return t1_*failsafe;
		}
	}
	
//transposition
	void transp()
	{
		T tmp=(*this)(0,0);
		for (int a=1; a<4; a++)
			for (int b=0; b<a; b++)
			{
				tmp=(*this)(a,b);
				(*this)(a,b)=(*this)(b,a);
				(*this)(b,a)=tmp;
			}
		
	}

inline void hermit()
	{
		T tmp=(*this)(0,0);
		for (int a=1; a<4; a++)
			for (int b=0; b<a; b++)
			{
				tmp=(*this)(a,b);
				tmp.hermit();
				(*this)(a,b)=(*this)(b,a);
				(*this)(a,b).hermit();
				(*this)(b,a)=tmp;
			}
		(*this)(0,0).hermit();
		(*this)(1,1).hermit();
		(*this)(2,2).hermit();
		(*this)(3,3).hermit();
		/*for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b).hermit();
			}*/
	
}


	//funny addition procedures: imagine you have a tensor ~g^{\mu\nu} and want to add a piece ~p^\mu p^\nu without
	//creating a temporary tensor. you .Add(p^\mu*something,p^\nu)	Saves processor time
	inline void Add(const D4V<double> &t,const D4V<double> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	inline void Add(const D4V<std::complex<double> > &t,const D4V<std::complex<double> > &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	inline void Add(const D4V<std::complex<double> > &t,const D4V<double> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	inline void Add(const D4V<double> &t,const D4V<std::complex<double> > &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	//funny addition procedures: imagine you have a tensor ~g^{\mu\nu} and want to add a piece ~\gamma^\mu \gamma^\nu without
	//creating a temporary tensor. you .Add(p^\mu*something,p^\nu)	
	inline void Add(const D4V<DM> &t,const D4V<DM> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}
	
	inline void Replace(const D4V<DM> &t,const D4V<DM> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)=t(a)*u(b);
		}
	}

	//funny addition procedures: imagine you have a tensor ~g^{\mu\nu} and want to add a piece ~\gamma^\mu p^\nu without
	//creating a temporary tensor. you .Add(\gamma^\mu*something,p^\nu)	
	
	inline void Add(const D4V<DM> &t,const D4V<double> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}
	
	inline void Replace(const D4V<DM> &t,const D4V<double> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)=t(a)*u(b);
		}
	}

	inline void Add(const D4V<double> &t,const D4V<DM> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}
	
	inline void Replace(const D4V<double> &t,const D4V<DM> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)=t(a)*u(b);
		}
	}


	inline void Replace(const D4V<std::complex<double>> &t,const D4V<std::complex<double>> &u)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)=t(a)*u(b);
			}
	}
	
	inline void Add(const D4V<DM> &t,const D4V<std::complex<double> > &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}

	inline void Add(const D4V<std::complex<double> > &t,const D4V<DM> &u)
	{
		for (unsigned short a=0; a<dim; a++)
		for (unsigned short b=0; b<dim; b++)
		{
				(*this)(a,b)+=t(a)*u(b);
		}
	}



	inline friend D4T operator+ (const D4T<T> &t1,const D4T<T> &t2)
	{
		D4T<T> result(t1.t0_,t2.t0_);

		result.t0_ = t1.t0_ +t2.t0_;
		result.t1_ = t1.t1_ +t2.t1_;
		result.t2_ = t1.t2_ +t2.t2_;
		result.t3_ = t1.t3_ +t2.t3_;
		return result;

	}


	inline D4T& operator+= (const D4T<T> &t1)
	{
		t0_ += t1.t0_;
		t1_ += t1.t1_;
		t2_ += t1.t2_;
		t3_ += t1.t3_;
		return *this;

	}

	inline friend D4T operator- (const D4T<T> &t1,const D4T<T> &t2)
	{
		D4T<T> result(t1.t0_,t2.t0_);

		result.t0_ = t1.t0_ -t2.t0_;
		result.t1_ = t1.t1_ -t2.t1_;
		result.t2_ = t1.t2_ -t2.t2_;
		result.t3_ = t1.t3_ -t2.t3_;
		return result;

	}

	inline D4T& operator-= (const D4T<T> &t1)
	{
		t0_ -= t1.t0_;
		t1_ -= t1.t1_;
		t2_ -= t1.t2_;
		t3_ -= t1.t3_;
		return *this;

	}

	//multiplication T1^\mu_\alpha T2^{\alpha \nu}
	inline friend D4T operator* (D4T<T> &T1,const D4T<T> &T2)
	{
		D4T<T> result(T1.t0_*0,T1.t0_*0);
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		//T2.transp();
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)=(T1(a)*tmptransp(b));
				
			}


		/*result.t0_.t0_=T1.t0_.t0_*T2.t0_.t0_ - T1.t0_.t1_*T2.t1_.t0_ - T1.t0_.t2_*T2.t2_.t0_ - T1.t0_.t3_*T2.t3_.t0_;
		result.t0_.t1_=T1.t0_.t0_*T2.t0_.t1_ - T1.t0_.t1_*T2.t1_.t1_ - T1.t0_.t2_*T2.t2_.t1_ - T1.t0_.t3_*T2.t3_.t1_;
		result.t0_.t2_=T1.t0_.t0_*T2.t0_.t2_ - T1.t0_.t1_*T2.t1_.t2_ - T1.t0_.t2_*T2.t2_.t2_ - T1.t0_.t3_*T2.t3_.t2_;
		result.t0_.t3_=T1.t0_.t0_*T2.t0_.t3_ - T1.t0_.t1_*T2.t1_.t3_ - T1.t0_.t2_*T2.t2_.t3_ - T1.t0_.t3_*T2.t3_.t3_;

		result.t1_.t0_=T1.t1_.t0_*T2.t0_.t0_ - T1.t1_.t1_*T2.t1_.t0_ - T1.t1_.t2_*T2.t2_.t0_ - T1.t1_.t3_*T2.t3_.t0_;
		result.t1_.t1_=T1.t1_.t0_*T2.t0_.t1_ - T1.t1_.t1_*T2.t1_.t1_ - T1.t1_.t2_*T2.t2_.t1_ - T1.t1_.t3_*T2.t3_.t1_;
		result.t1_.t2_=T1.t1_.t0_*T2.t0_.t2_ - T1.t1_.t1_*T2.t1_.t2_ - T1.t1_.t2_*T2.t2_.t2_ - T1.t1_.t3_*T2.t3_.t2_;
		result.t1_.t3_=T1.t1_.t0_*T2.t0_.t3_ - T1.t1_.t1_*T2.t1_.t3_ - T1.t1_.t2_*T2.t2_.t3_ - T1.t1_.t3_*T2.t3_.t3_;

		result.t2_.t0_=T1.t2_.t0_*T2.t0_.t0_ - T1.t2_.t1_*T2.t1_.t0_ - T1.t2_.t2_*T2.t2_.t0_ - T1.t2_.t3_*T2.t3_.t0_;
		result.t2_.t1_=T1.t2_.t0_*T2.t0_.t1_ - T1.t2_.t1_*T2.t1_.t1_ - T1.t2_.t2_*T2.t2_.t1_ - T1.t2_.t3_*T2.t3_.t1_;
		result.t2_.t2_=T1.t2_.t0_*T2.t0_.t2_ - T1.t2_.t1_*T2.t1_.t2_ - T1.t2_.t2_*T2.t2_.t2_ - T1.t2_.t3_*T2.t3_.t2_;
		result.t2_.t3_=T1.t2_.t0_*T2.t0_.t3_ - T1.t2_.t1_*T2.t1_.t3_ - T1.t2_.t2_*T2.t2_.t3_ - T1.t2_.t3_*T2.t3_.t3_;

		result.t3_.t0_=T1.t3_.t0_*T2.t0_.t0_ - T1.t3_.t1_*T2.t1_.t0_ - T1.t3_.t2_*T2.t2_.t0_ - T1.t3_.t3_*T2.t3_.t0_;
		result.t3_.t1_=T1.t3_.t0_*T2.t0_.t1_ - T1.t3_.t1_*T2.t1_.t1_ - T1.t3_.t2_*T2.t2_.t1_ - T1.t3_.t3_*T2.t3_.t1_;
		result.t3_.t2_=T1.t3_.t0_*T2.t0_.t2_ - T1.t3_.t1_*T2.t1_.t2_ - T1.t3_.t2_*T2.t2_.t2_ - T1.t3_.t3_*T2.t3_.t2_;
		result.t3_.t3_=T1.t3_.t0_*T2.t0_.t3_ - T1.t3_.t1_*T2.t1_.t3_ - T1.t3_.t2_*T2.t2_.t3_ - T1.t3_.t3_*T2.t3_.t3_;*/
		
		return result;
		
	}
	
	inline friend D4T operator* (const D4T<T> &T1,const double &d)
	{
		D4T<T> result=T1;
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)*=d;
			}
		return result;
		
	}

    inline friend D4T operator* (const D4T<T> &T1,const std::complex<double> &d)
    {
        D4T<T> result=T1;
        for (unsigned short a=0; a<dim; a++)
            for (unsigned short b=0; b<dim; b++)
            {
                result(a,b)*=d;
            }
        return result;
    }

	
	inline friend D4T operator* (const std::complex<double> &d, const D4T<T> &T1)
	{
		D4T<T> result=T1;
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)*=d;
			}
		return result;
		
	}
	
	inline friend D4T operator* (const D4T<T> &T1,const DM &d)
	{
		D4T<T> result=T1;
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				result(a,b)*=d;
			}
		return result;
		
	}
	
	//multiplication V_\alpha T^{\alpha \nu}
	
	friend D4V<T> operator* (const D4V<double> &T1,const D4T<T> &T2)
	{
		D4V<T> result(T2(0,0),T2(0,0),T2(0,0),T2(0,0));
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1*tmptransp(a));
		return result;
		
	}
///Warning!!! this does not recognize right types. I do not know, why. Yet
	friend D4V<T> operator* (D4T<T> &T1, const D4V<double> &T2)
	{
		DM zero(222);
		D4V<DM> result(zero,zero,zero,zero);
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1(a)*T2);
		return result;
		
	}

	friend D4V<DM> operator* (const D4V<DM> &T1,const D4T<T> &T2)
	{
		DM zero(222);
		D4V<DM> result(zero,zero,zero,zero);
		D4T<T> tmptransp =T2;
		tmptransp.transp();
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1*tmptransp(a));
		return result;
		
	}
///Warning!!!
	friend D4V<DM> operator* (D4T<T> &T1,const D4V<DM> &T2)
	{
		DM zero(222);
		D4V<DM> result(zero,zero,zero,zero);
		for (unsigned short a=0; a<dim; a++)
			result(a)=(T1(a)*T2);
		return result;
		
	}

	//multiplication T1^\mu_\alpha T2^{\alpha \nu}
	inline D4T& operator*= (const D4T<T> &T2)
	{
		D4T<T> result(t0_.t0_*T2.t0_.t0_ -  t0_.t1_*T2.t1_.t0_ -  t0_.t2_*T2.t2_.t0_ -  t0_.t3_*T2.t3_.t0_,
		t0_.t0_*T2.t0_.t1_ -  t0_.t1_*T2.t1_.t1_ -  t0_.t2_*T2.t2_.t1_ -  t0_.t3_*T2.t3_.t1_,
		t0_.t0_*T2.t0_.t2_ -  t0_.t1_*T2.t1_.t2_ -  t0_.t2_*T2.t2_.t2_ -  t0_.t3_*T2.t3_.t2_,
		t0_.t0_*T2.t0_.t3_ -  t0_.t1_*T2.t1_.t3_ -  t0_.t2_*T2.t2_.t3_ -  t0_.t3_*T2.t3_.t3_,
		t1_.t0_*T2.t0_.t0_ -  t1_.t1_*T2.t1_.t0_ -  t1_.t2_*T2.t2_.t0_ -  t1_.t3_*T2.t3_.t0_,
		t1_.t0_*T2.t0_.t1_ -  t1_.t1_*T2.t1_.t1_ -  t1_.t2_*T2.t2_.t1_ -  t1_.t3_*T2.t3_.t1_,
		t1_.t0_*T2.t0_.t2_ -  t1_.t1_*T2.t1_.t2_ -  t1_.t2_*T2.t2_.t2_ -  t1_.t3_*T2.t3_.t2_,
		t1_.t0_*T2.t0_.t3_ -  t1_.t1_*T2.t1_.t3_ -  t1_.t2_*T2.t2_.t3_ -  t1_.t3_*T2.t3_.t3_,
		t2_.t0_*T2.t0_.t0_ -  t2_.t1_*T2.t1_.t0_ -  t2_.t2_*T2.t2_.t0_ -  t2_.t3_*T2.t3_.t0_,
		t2_.t0_*T2.t0_.t1_ -  t2_.t1_*T2.t1_.t1_ -  t2_.t2_*T2.t2_.t1_ -  t2_.t3_*T2.t3_.t1_,
		t2_.t0_*T2.t0_.t2_ -  t2_.t1_*T2.t1_.t2_ -  t2_.t2_*T2.t2_.t2_ -  t2_.t3_*T2.t3_.t2_,
		t2_.t0_*T2.t0_.t3_ -  t2_.t1_*T2.t1_.t3_ -  t2_.t2_*T2.t2_.t3_ -  t2_.t3_*T2.t3_.t3_,
		t3_.t0_*T2.t0_.t0_ -  t3_.t1_*T2.t1_.t0_ -  t3_.t2_*T2.t2_.t0_ -  t3_.t3_*T2.t3_.t0_,
		t3_.t0_*T2.t0_.t1_ -  t3_.t1_*T2.t1_.t1_ -  t3_.t2_*T2.t2_.t1_ -  t3_.t3_*T2.t3_.t1_,
		t3_.t0_*T2.t0_.t2_ -  t3_.t1_*T2.t1_.t2_ -  t3_.t2_*T2.t2_.t2_ -  t3_.t3_*T2.t3_.t2_,
		t3_.t0_*T2.t0_.t3_ -  t3_.t1_*T2.t1_.t3_ -  t3_.t2_*T2.t2_.t3_ -  t3_.t3_*T2.t3_.t3_);
		*this=result;
		return *this;
	}
	
	inline D4T& operator&= (const D4T<T> &T1)
	{
		D4T<T> resu(T1.t0_.t0_*t0_.t0_ - T1.t0_.t1_*t1_.t0_ - T1.t0_.t2_*t2_.t0_ - T1.t0_.t3_*t3_.t0_,
					T1.t0_.t0_*t0_.t1_ - T1.t0_.t1_*t1_.t1_ - T1.t0_.t2_*t2_.t1_ - T1.t0_.t3_*t3_.t1_,
					T1.t0_.t0_*t0_.t2_ - T1.t0_.t1_*t1_.t2_ - T1.t0_.t2_*t2_.t2_ - T1.t0_.t3_*t3_.t2_,
					T1.t0_.t0_*t0_.t3_ - T1.t0_.t1_*t1_.t3_ - T1.t0_.t2_*t2_.t3_ - T1.t0_.t3_*t3_.t3_,
					T1.t1_.t0_*t0_.t0_ - T1.t1_.t1_*t1_.t0_ - T1.t1_.t2_*t2_.t0_ - T1.t1_.t3_*t3_.t0_,
					T1.t1_.t0_*t0_.t1_ - T1.t1_.t1_*t1_.t1_ - T1.t1_.t2_*t2_.t1_ - T1.t1_.t3_*t3_.t1_,
					T1.t1_.t0_*t0_.t2_ - T1.t1_.t1_*t1_.t2_ - T1.t1_.t2_*t2_.t2_ - T1.t1_.t3_*t3_.t2_,
					T1.t1_.t0_*t0_.t3_ - T1.t1_.t1_*t1_.t3_ - T1.t1_.t2_*t2_.t3_ - T1.t1_.t3_*t3_.t3_,
					T1.t2_.t0_*t0_.t0_ - T1.t2_.t1_*t1_.t0_ - T1.t2_.t2_*t2_.t0_ - T1.t2_.t3_*t3_.t0_,
					T1.t2_.t0_*t0_.t1_ - T1.t2_.t1_*t1_.t1_ - T1.t2_.t2_*t2_.t1_ - T1.t2_.t3_*t3_.t1_,
					T1.t2_.t0_*t0_.t2_ - T1.t2_.t1_*t1_.t2_ - T1.t2_.t2_*t2_.t2_ - T1.t2_.t3_*t3_.t2_,
					T1.t2_.t0_*t0_.t3_ - T1.t2_.t1_*t1_.t3_ - T1.t2_.t2_*t2_.t3_ - T1.t2_.t3_*t3_.t3_,
					T1.t3_.t0_*t0_.t0_ - T1.t3_.t1_*t1_.t0_ - T1.t3_.t2_*t2_.t0_ - T1.t3_.t3_*t3_.t0_,
					T1.t3_.t0_*t0_.t1_ - T1.t3_.t1_*t1_.t1_ - T1.t3_.t2_*t2_.t1_ - T1.t3_.t3_*t3_.t1_,
					T1.t3_.t0_*t0_.t2_ - T1.t3_.t1_*t1_.t2_ - T1.t3_.t2_*t2_.t2_ - T1.t3_.t3_*t3_.t2_,
					T1.t3_.t0_*t0_.t3_ - T1.t3_.t1_*t1_.t3_ - T1.t3_.t2_*t2_.t3_ - T1.t3_.t3_*t3_.t3_);
		*this=resu;
		return *this;
		
	}

	inline D4T& operator*= (const std::complex<double> &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}

	inline D4T& operator*= (const double &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}


	inline D4T& operator*= (const DM &d)
	{
		for (unsigned short a=0; a<dim; a++)
			for (unsigned short b=0; b<dim; b++)
			{
				(*this)(a,b)*=d;
			}
		return *this;
	}

	inline D4T& operator&= (const DM &d)
	{
		for (unsigned short a=0; a<4; a++)
			for (unsigned short b=0; b<4; b++)
			{
				(*this)(a,b)=d*(*this)(a,b);
			}
		return *this;
	}

    //returns the Real part
    inline D4T< double > Re(const D4T<std::complex<double>> &D )
    {
        D4T<double> result;
        
        for (unsigned short a=0; a<dim; a++)
            for(unsigned short b=0; b<dim;b++)
                result(a,b)=D(a,b).real();
        return result;
    }

    
	//returns the trace as a matrix
	inline D4T< std::complex<double> > CT(const D4T<DM> &D )
	{
		D4V< std::complex<double> > dummy(std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(0,0));
		D4T< std::complex<double> > result(dummy,dummy);

		for (unsigned short a=0; a<dim; a++)
			for(unsigned short b=0; b<dim;b++)
			result(a,b)=D(a,b).CTrace();
		return result;
	}

    //returns the Real part trace as a matrix
    inline D4T< double > ReT(const D4T<DM> &D )
    {
         D4T<double> result;
        
        for (unsigned short a=0; a<dim; a++)
            for(unsigned short b=0; b<dim;b++)
                result(a,b)=D(a,b).Trace();
        return result;
    }

    //returns the Imaginary part trace as a matrix
    inline D4T< double > ImT(const D4T<DM> &D )
    {
        D4T<double> result;
        
        for (unsigned short a=0; a<dim; a++)
            for(unsigned short b=0; b<dim;b++)
                result(a,b)=(D(a,b).CTrace()).imag();
        return result;
    }
    
	/*inline D4T<std::complex<double> > CT1(const D4T<DM> &D )
	{
		D4T< std::complex<double> > result((t(0)*u(0)).CTrace(), (t(0)*u(1)).CTrace(), (t(0)*u(2)).CTrace(), t(0)*u(3),
											t(1)*u(0), t(1)*u(1), t(1)*u(2), t(1)*u(3),
											t(2)*u(0), t(2)*u(1), t(2)*u(2), t(2)*u(3),
											t(3)*u(0), t(3)*u(1), t(3)*u(2), t(3)*u(3))
		return result;
	}*/


	//antisymmetrization by -\epsilon^{\mu\nu\alpha\beta}T_{\alpha\beta}
	inline D4T Ants(const D4T<T> &t)
	{
		D4T<T> result=t;
		double zero=0;
		double minus=-1.0;
		result.t0_.t0_*=zero;
		result.t1_.t1_*=zero;
		result.t2_.t2_*=zero;
		result.t3_.t3_*=zero;

		result.t0_.t1_=t.t2_.t3_-t.t3_.t2_;
		result.t1_.t0_=(result.t0_.t1_*minus);

		result.t0_.t2_=t.t3_.t1_-t.t1_.t3_;
		result.t2_.t0_=(result.t0_.t2_*minus);

		result.t0_.t3_=t.t1_.t2_-t.t2_.t1_;
		result.t3_.t0_=(result.t0_.t3_*minus);

		result.t1_.t2_=t.t3_.t0_-t.t0_.t3_;
		result.t2_.t1_=(result.t1_.t2_*minus);

		result.t2_.t3_=t.t1_.t0_-t.t0_.t1_;
		result.t3_.t2_=(result.t2_.t3_*minus);

		result.t3_.t1_=t.t2_.t0_-t.t0_.t2_;
		result.t1_.t3_=(result.t3_.t1_*minus);
		//result*=minus;
		return result;
	}
	//adds antisymmetric part -\epsilon^{\mu\nu\alpha\beta}k_\alpha k'_\beta
	inline void AddAnts(const D4V<double> &u,const D4V<double> &v)
	{
		double aaa=(u(2)*v(3)-u(3)*v(2));
		
		(*this)(0,1)+=aaa;
		(*this)(1,0)-=aaa;
		
		aaa=(u(3)*v(1)-u(1)*v(3));
		(*this)(0,2)+=aaa;
		(*this)(2,0)-=aaa;

		aaa=(u(1)*v(2)-u(2)*v(1));
		(*this)(0,3)+=aaa;
		(*this)(3,0)-=aaa;

		aaa=(u(3)*v(0)-u(0)*v(3));
		(*this)(1,2)+=aaa;
		(*this)(2,1)-=aaa;

		aaa=(u(1)*v(0)-u(0)*v(1));
		(*this)(2,3)+=aaa;
		(*this)(3,2)-=aaa;

		aaa=(u(2)*v(0)-u(0)*v(2));
		(*this)(3,1)+=aaa;
		(*this)(1,3)-=aaa;
	}

    //adds antisymmetric part -i\epsilon^{\mu\nu\alpha\beta}k_\alpha k'_\beta
    inline void AddAntsImag(const D4V<double> &u,const D4V<double> &v)
    {
        std::complex<double> im(0,1);
        std::complex<double> aaa=(u(2)*v(3)-u(3)*v(2))*im;
        
        (*this)(0,1)+=aaa;
        (*this)(1,0)-=aaa;
        
        aaa=(u(3)*v(1)-u(1)*v(3));
        (*this)(0,2)+=aaa;
        (*this)(2,0)-=aaa;
        
        aaa=(u(1)*v(2)-u(2)*v(1));
        (*this)(0,3)+=aaa;
        (*this)(3,0)-=aaa;
        
        aaa=(u(3)*v(0)-u(0)*v(3));
        (*this)(1,2)+=aaa;
        (*this)(2,1)-=aaa;
        
        aaa=(u(1)*v(0)-u(0)*v(1));
        (*this)(2,3)+=aaa;
        (*this)(3,2)-=aaa;
        
        aaa=(u(2)*v(0)-u(0)*v(2));
        (*this)(3,1)+=aaa;
        (*this)(1,3)-=aaa;
    }
    
	//T1^{\mu\nu}T2_{\mu\nu} for the same type
	inline T contraction(const D4T<T> &t1,const D4T<T> &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}	
	
	//type of the second tensor is <double>
	inline T contrd (const D4T<T> &t1,const D4T<double> &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}
	//type of the second tensor is complex <double>
	
	inline T contrc (const D4T<T> &t1,const D4T<std::complex<double> > &t2)
	{
		T result =(t1(0,0)*t2(0,0) + t1(1,1)*t2(1,1) + t1(2,2)*t2(2,2) + t1(3,3)*t2(3,3));
		  result+=(t1(1,2)*t2(1,2) + t1(2,1)*t2(2,1) + t1(2,3)*t2(2,3) + t1(3,2)*t2(3,2) + t1(3,1)*t2(3,1) + t1(1,3)*t2(1,3));
		  result-=(t1(1,0)*t2(1,0) + t1(0,1)*t2(0,1) + t1(2,0)*t2(2,0) + t1(0,2)*t2(0,2) + t1(3,0)*t2(3,0) + t1(0,3)*t2(0,3));
		return result;
	}
		D4V<T> t0_;
		D4V<T> t1_;
		D4V<T> t2_;
		D4V<T> t3_;
};



#endif
