
#ifndef SVG_UTILS_H
#define SVG_UTILS_H

template <typename Real>
void svg_linestops(const Real *pts, int nb, const Real *dir, Real *stopa, Real *stopb)
{
	Real norm = dir[0] * dir[0] + dir[1] * dir[1];
	assert(norm > 0.0);

	norm = std::sqrt(norm);
	Real udir[2] = { dir[0] / norm, dir[1] / norm };

	const Real *p = &pts[0];
	Real tmin = Real(0.0), tmax = Real(0.0);

	for (int k = 1; k < nb; ++k)
	{
		const Real *q = &pts[2 * k];
		Real pq[2] = { q[0] - p[0], q[1] - p[1] };
		
		double dot = pq[0] * udir[0] + pq[1] * udir[1];
		if (dot < tmin)
			tmin = dot;
		if (dot > tmax)
			tmax = dot;
	}

	stopa[0] = p[0] + tmin * udir[0];
	stopa[1] = p[1] + tmin * udir[1];
	stopb[0] = p[0] + tmax * udir[0];
	stopb[1] = p[1] + tmax * udir[1];
}

#endif
