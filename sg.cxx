#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;

//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);
void step(cmplx *u0, cmplx* u1, const double dx, const double dt, const int Nx,
		  const double omega, const double xmin);

void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	
	const int Nx =300 ;
	
	const double xmin = -40;
    const double xmax =40 ;
	
	const double Tend = 10*PI;
	const double dx =(xmax-xmin)/(Nx-1) ;
	const double dt =dx*0.01;
    double t = 0.;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
    const double omega = 0.2;
	const double k=pow(omega,2.);
	const cmplx  ii=cmplx(0.0, 1.0);
	const double alpha=sqrt(omega);

    stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* u1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0.txt", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		 
		 step(psi0, u1, dx, dt, Nx, omega, xmin);
		 h=psi0;
		 psi0=u1;
		 u1=h;

         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i<<".txt";
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
    cout << "t = " << t << endl;
	delete[] psi0;
	delete[] u1;
	return 0;
}
//-----------------------------------
void step(cmplx *u0, cmplx* u1, const double dx, const double dt, const int Nx, 
		  const double omega, double xmin)
{
	double x;
	cmplx* d= new cmplx[Nx];
	cmplx* dcon= new cmplx[Nx];  // konjugierte wert
	//cmplx a=-ii*dt/(4*dx*dx);

	cmplx* u=new cmplx[Nx];
	cmplx* ucon= new cmplx[Nx];

	cmplx* l=new cmplx[Nx]; 
	cmplx* lcon= new cmplx[Nx];

	for(int i=0; i<Nx; i++)
	{
		x=xmin+i*dx;
		d[i]=1.+cmplx(0.0, dt/(2*dx*dx)+dt*0.25*pow(omega,2.)*pow(x,2.));
		dcon[i]=1.-cmplx(0.0, dt/(2*dx*dx)+dt*0.25*pow(omega,2.)*pow(x,2.));
		u[i]=cmplx(0.0, -dt/(4*dx*dx));
		ucon[i]=cmplx(0.0, dt/(4*dx*dx));
		l[i]=cmplx(0.0, -dt/(4*dx*dx));
		lcon[i]=cmplx(0.0, dt/(4*dx*dx));

	}

	for(int i=0; i<Nx-1; i++)
	{
		d[i+1]=d[i+1]-l[i+1]/d[i]*u[i];
		

		lcon[i+1]=lcon[i+1]-l[i+1]/d[i]*dcon[i];
		dcon[i+1]=dcon[i+1]-l[i+1]/d[i]*ucon[i];
	}

	//backward substitution

	u1[Nx-1]=(lcon[Nx-1]*u0[Nx-2]+dcon[Nx-1]*u0[Nx-1])/d[Nx-1]; //letzte zeile extra

	for(int i=Nx-2; i>=1; i--)
		u1[i]=(lcon[i]*u0[i-1]
			+dcon[i]*u0[i]
			+ucon[i]*u0[i+1]
			-u[i]*u1[i+1])/d[i];
	
	u1[0]=(dcon[0]*u0[0]+ucon[0]*u0[1]-u[0]*u1[1])/d[0]; //erste zeile extra

	delete[] d;
	delete[] dcon;
	delete[] l;
	delete[] lcon;
	delete[] u;
	delete[] ucon;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
 
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
    out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
        << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	//out<<x <<"\t "<< norm(v[i]) << "\t "<< norm(ana)<<endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
