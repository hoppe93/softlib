#ifndef _RKDP45_H
#define _RKDP45_H

#include <softlib/RungeKutta.h>
#include <softlib/Vector.h>

template<unsigned int N>
class RKDP45 : public RungeKutta<N> {
	private:
		Vector<N> k2, k3, k4, k5, k6;
		slibreal_t *rcont1=NULL, *rcont2=NULL, *rcont3=NULL, *rcont4=NULL, *rcont5=NULL;
		unsigned int nrcontallocated=0;
	
		/* Dormand-Prince coefficients */
		const slibreal_t
			c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0,
			a21=0.2,
			a31=3.0/40.0, a32=9.0/40.0,
			a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
			a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
			a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
			a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0,
			e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;

		/* Coefficients for generating "dense" output */
		const slibreal_t
			d1=-12715105075.0/11282082432.0,
			d3=87487479700.0/32700410799.0,
			d4=-10690763975.0/1880347072.0,
			d5=701980252875.0/199316789632.0,
			d6=-1453857185.0/822651844.0,
			d7=69997945.0/29380423.0;
	public:
        RKDP45<N>(const slibreal_t tol=RK_DEFAULT_TOLERANCE) : RungeKutta<N>(tol) {}
		Vector<N>& InnerStep(const slibreal_t, const slibreal_t, const Vector<N>&, const Vector<N>&);
		void PrepareDense(const slibreal_t, const Vector<N>&, const Vector<N>&, const Vector<N>&, const Vector<N>&);
		void OutputDense(const unsigned int, slibreal_t*, slibreal_t*);

        void OutputDense(const unsigned int nt, const slibreal_t tmin, const slibreal_t tmax, slibreal_t* out, slibreal_t* outt) {
            Integrator<N>::OutputDense(nt, tmin, tmax, out, outt);
        }
};

// Implementation
#include <softlib/Integrator/RKDP45.tcc>

#endif/*_RKDP45_H*/
