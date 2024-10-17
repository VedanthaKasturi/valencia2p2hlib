#include "MatrixDir.h"
//#include <iostream>
//4x4 complex matrix constructor//Dirac algebra in the Dirac representation!



DM::DM(const double &a, const double &b, const double &c, const double &d)
{
		matrix[0][0]= std::complex<double>(a,0);
		matrix[0][1]= std::complex<double>(0,0);
		matrix[0][2]= std::complex<double>(-d,0);
		matrix[0][3]= std::complex<double>(-b,c);//xx
		matrix[1][0]= std::complex<double>(0,0);
		matrix[1][1]= std::complex<double>(a,0);
		matrix[1][2]= std::complex<double>(-b,-c);//xx
		matrix[1][3]= std::complex<double>(d,0);
		matrix[2][0]= std::complex<double>(d,0);
		matrix[2][1]= std::complex<double>(b, -c);
		matrix[2][2]= std::complex<double>(-a,0);
		matrix[2][3]= std::complex<double>(0,0);
		matrix[3][0]= std::complex<double>(b,c);//xx
		matrix[3][1]= std::complex<double>(-d,0);
		matrix[3][2]= std::complex<double>(0,0);
		matrix[3][3]= std::complex<double>(-a,0);		
}

DM::DM(const double &a, const double &b, const double &c, const double &d, const double &e)
{
		matrix[0][0]= std::complex<double>(a+e,0);
		matrix[0][1]= std::complex<double>(0,0);
		matrix[0][2]= std::complex<double>(-d,0);
		matrix[0][3]= std::complex<double>(-b,c);//xx
		matrix[1][0]= std::complex<double>(0,0);
		matrix[1][1]= std::complex<double>(a+e,0);
		matrix[1][2]= std::complex<double>(-b,-c);//xx
		matrix[1][3]= std::complex<double>(d,0);
		matrix[2][0]= std::complex<double>(d,0);
		matrix[2][1]= std::complex<double>(b, -c);
		matrix[2][2]= std::complex<double>(e-a,0);
		matrix[2][3]= std::complex<double>(0,0);
		matrix[3][0]= std::complex<double>(b,c);//xx
		matrix[3][1]= std::complex<double>(-d,0);
		matrix[3][2]= std::complex<double>(0,0);
		matrix[3][3]= std::complex<double>(e-a,0);		
}


//Hermitian adjounsigned short
void DM::hermit()
{
		std::complex<double> tmp;
		matrix[0][0]=conj(matrix[0][0]);
		matrix[1][1]=conj(matrix[1][1]);
		matrix[2][2]=conj(matrix[2][2]);
		matrix[3][3]=conj(matrix[3][3]);
		for (unsigned short a=1; a<4; a++)
			for (unsigned short b=0; b<a; b++)
			{
				tmp=conj(matrix[a][b]);
				matrix[a][b]=conj(matrix[b][a]);
				matrix[b][a]=tmp;
			}
		/*for (unsigned short a=0; a<4; a++)
		for (unsigned short b=0; b<4; b++)
		{
			tmp[a][b]=conj(matrix[b][a]);
	
		}
		for (unsigned short a=0; a<4; a++)
		for (unsigned short b=0; b<4; b++)
		{
			matrix[a][b]=tmp[a][b];
	
		}*/
}

//transposition
void DM::transp()
{
		std::complex<double> tmp;
		
		/*for (unsigned short a=0; a<4; a++)
		for (unsigned short b=0; b<4; b++)
		{
			tmp[a][b]=matrix[b][a];
	
		}
		for (unsigned short a=0; a<4; a++)
		for (unsigned short b=0; b<4; b++)
		{
			matrix[a][b]=tmp[a][b];
	
		}*/
		
		for (unsigned short a=1; a<4; a++)
			for (unsigned short b=0; b<a; b++)
			{
				tmp=matrix[a][b];
				matrix[a][b]=matrix[b][a];
				matrix[b][a]=tmp;
			}
}

void DM::wrt()
{
	std::cout<<"\n\n";
	for (unsigned short a=0; a<4; a++)
	{
		for (unsigned short b=0; b<4; b++)
		{
			std::cout<<matrix[a][b]<<"  ";
	
		}
		std::cout<<"\n\n";
	}
	std::cout<<"\n\n";
}



