
#ifndef DOUGLAS_PEUCKER_H
#define DOUGLAS_PEUCKER_H

template <typename Real>
Real perp_distance(const Real *src, const Real *end, const Real *query)
{
	Real dx = end[0] - src[0];
	Real dy = end[1] - src[1];

	//Normalize
	Real mag = sqrt(dx * dx + dy * dy);
	if (mag > 0.0)
	{
		dx /= mag;
		dy /= mag;
	}

	Real pvx = query[0] - src[0];
	Real pvy = query[1] - src[1];

	//Get dot product (project pv onto normalized direction)
	Real pvdot = dx * pvx + dy * pvy;

	//Scale line direction vector
	Real dsx = pvdot * dx;
	Real dsy = pvdot * dy;

	//Subtract this from pv
	Real ax = pvx - dsx;
	Real ay = pvy - dsy;

	return sqrt(ax * ax + ay * ay);
}

template <typename Real>
void Douglas_Peucker(const Real *points, int a, int b, double distT, int nT, bool *picked)
{
	double pa[2] = { double(points[2 * a]), double(points[2 * a + 1]) };
	double pb[2] = { double(points[2 * b]), double(points[2 * b + 1]) };

	double maxdist = 0.0;
	int maxk = -1;

	for (int k = a + 1; k < b; ++k)
	{
		double pc[2] = { double(points[2 * k]), double(points[2 * k + 1]) };

		double dist = perp_distance(pa, pb, pc);
		if (dist > maxdist)
		{
			maxdist = dist;
			maxk = k;
		}
	}

	if (maxdist < distT)
	{
		assert(nT > 1);
		if (b - a + 1 <= nT)
			return;

		maxk = (a + b) / 2;
	}

	picked[a] = true;
	picked[b] = true;
	
	if (maxk > -1)
	{
		Douglas_Peucker(points, a, maxk, distT, nT, picked);
		Douglas_Peucker(points, maxk, b, distT, nT, picked);
	}
}

#endif
