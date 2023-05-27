
#ifndef BEZIER2_H
#define BEZIER2_H

template <typename Real>
void bezier2_point(const Real *vi, const Real *cij, const Real *vj, Real t, Real *p)
{
	p[0] = (1 - t) * (1 - t) * vi[0] +
		2 * t * (1 - t) * cij[0] +
		t * t * vj[0];
	p[1] = (1 - t) * (1 - t) * vi[1] +
		2 * t * (1 - t) * cij[1] +
		t * t * vj[1];
}

#endif
