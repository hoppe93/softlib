/**
 * Implementation of the Dormand-Prince
 * fifth order Runge-Kutta method.
 */

#include <softlib/RKDP45.h>
#include <softlib/RungeKutta.h>
#include <softlib/SOFTLibException.h>
#include <softlib/Vector.h>

/**
 * Advance the integrator a step
 * 'timestep' in time.
 */
template<unsigned int N>
Vector<N>& RKDP45<N>::InnerStep(const slibreal_t T, const slibreal_t h, const Vector<N>& y, const Vector<N>& dydx) {
	Vector<N> ytemp;

	/* Get old solution */
	ytemp = y + h*a21*dydx;
	
	this->eqn->Evaluate(T+c2*h, ytemp, k2);
	ytemp = y + (h*a31)*dydx + (h*a32)*k2;

	this->eqn->Evaluate(T+c3*h, ytemp, k3);
	ytemp = y + (h*a41)*dydx + (h*a42)*k2 + (h*a43)*k3;

	this->eqn->Evaluate(T+c4*h, ytemp, k4);
	ytemp = y + (h*a51)*dydx + (h*a52)*k2 + (h*a53)*k3 + (h*a54)*k4;

	this->eqn->Evaluate(T+c5*h, ytemp, k5);
	ytemp = y + (h*a61)*dydx + (h*a62)*k2 + (h*a63)*k3 + (h*a64)*k4 + (h*a65)*k5;

	this->eqn->Evaluate(T+h, ytemp, k6);
	this->yout = y + (h*a71)*dydx + (h*a73)*k3 + (h*a74)*k4 + (h*a75)*k5 + (h*a76)*k6;

	this->eqn->Evaluate(T+h, this->yout, this->dydxout);
	this->yerr = (h*e1)*dydx + (h*e3)*k3 + (h*e4)*k4 + (h*e5)*k5 + (h*e6)*k6 + (h*e7)*(this->dydxout);

	return this->yout;
}

/**
 * Calculate coefficients needed to
 * generate dense output at the end
 * of the run.
 *
 * h:     Time-step used to reach current point.
 * y1:    Previous solution (y1 = y(T-h)).
 * y2:    Current solution (y2 = y(T)).
 * dydx1: Derivative in previous time-step.
 * dydx2: Derivative in current time-step.
 */
template<unsigned int N>
void RKDP45<N>::PrepareDense(
	const slibreal_t h,
	const Vector<N>& y1, const Vector<N>& y2,
	const Vector<N>& dydx1, const Vector<N>& dydx2
) {
	slibreal_t ydiff, bspl;
	unsigned int i, t = (this->ntimesteps-2)*N;

	if (nrcontallocated < this->ntsallocated) {
		rcont1 = (slibreal_t*)realloc(rcont1, sizeof(slibreal_t) * N * this->ntsallocated);
		rcont2 = (slibreal_t*)realloc(rcont2, sizeof(slibreal_t) * N * this->ntsallocated);
		rcont3 = (slibreal_t*)realloc(rcont3, sizeof(slibreal_t) * N * this->ntsallocated);
		rcont4 = (slibreal_t*)realloc(rcont4, sizeof(slibreal_t) * N * this->ntsallocated);
		rcont5 = (slibreal_t*)realloc(rcont5, sizeof(slibreal_t) * N * this->ntsallocated);

		nrcontallocated = this->ntsallocated;
	}

	for (i = 0; i < N; i++) {
		rcont1[t+i] = y1[i];
		ydiff = y2[i]-y1[i];
		rcont2[t+i] = ydiff;
		bspl = h*dydx1[i]-ydiff;
		rcont3[t+i] = bspl;
		rcont4[t+i] = ydiff - h*dydx2[i] - bspl;
		rcont5[t+i] = h*(d1*dydx1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*dydx2[i]);
	}
}

/**
 * Generate dense output.
 *
 * nt:   Number of elements in 'outt' (i.e. number of timesteps).
 * time: Array containing times to return solution for.
 * out:  Output vector. Contains dense output on return.
 *       Must be of size nt * N.
 */
template<unsigned int N>
//void RKDP45<N>::OutputDense(const unsigned int nt, const slibreal_t tmin, const slibreal_t tmax, slibreal_t *out, slibreal_t *outt) {
void RKDP45<N>::OutputDense(const unsigned int nt, slibreal_t *time, slibreal_t *out) {
	unsigned int it = 0, i, j, indx;
	slibreal_t T, s, s1, stepsize;

	for (i = 0; i < nt; i++) {
        T = time[i];

		/* Find next time */
		while (it+2 < this->ntimesteps && this->time[it+1] < T)
			it++;

		if (this->time[it+1] < T)
			throw SOFTLibException("The requested tmax requires extrapolation.");

		stepsize = this->time[it+1]-this->time[it];
		
		s = (T-this->time[it]) / stepsize;
		s1 = 1.0 - s;
		for (j = 0; j < N; j++) {
			indx = it*N + j;
			out[i*N+j] = rcont1[indx] + s*(rcont2[indx] + s1*(rcont3[indx] + s*(rcont4[indx] + s1*rcont5[indx])));
		}
	}
}

