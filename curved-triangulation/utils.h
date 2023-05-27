
#ifndef UTILS_H
#define UTILS_H

template <typename Real>
bool loop_path(const Real *path, int nb, Real t)
{
	const Real *first = path;
	const Real *last = &path[2 * nb - 2];

	Real dx = abs(last[0] - first[0]);
	Real dy = abs(last[1] - first[1]);

	if (dx <= t || dy <= t)
		return true;

	return false;
}

// distance from a to bc
template <typename Real>
Real point_line_distance(const Real *a, const Real *b, const Real *c)
{
	Real dir[2] = { c[0] - b[0], c[1] - b[1] };
	Real norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
	dir[0] /= norm;
	dir[1] /= norm;

	Real dist = (a[0] - b[0]) * dir[0] + (a[1] - b[1]) * dir[1];
	return dist;
}

template <typename Real>
bool inside_triangle(const Real *tri, const Real *query)
{
	Real vec[3][2];
	for (int k = 0; k < 3; ++k)
	{
		vec[k][0] = tri[2 * k] - query[0];
		vec[k][1] = tri[2 * k + 1] - query[1];
	}

	Real cross[3];
	for (int k = 0; k < 3; ++k)
	{
		int next = (k + 1) % 3;
		cross[k] = vec[k][0] * vec[next][1] - vec[k][1] * vec[next][0];
	}

	for (int k = 0; k < 3; ++k)
	{
		int next = (k + 1) % 3;
		if (cross[k] * cross[next] < 0)
			return false;
	}

	return true;
}

template <typename Real>
bool inside_convex(const Real *polygon, int size, const Real *query)
{
	Real a[2], b[2];
	
	a[0] = polygon[2 * (size - 1)] - query[0];
	a[1] = polygon[2 * size - 1] - query[1];

	b[0] = polygon[0] - query[0];
	b[1] = polygon[1] - query[1];

	double base = a[0] * b[1] - a[1] * b[0];

	a[0] = b[0];
	a[1] = b[1];

	for (int k = 1; k < size; ++k)
	{
		b[0] = polygon[2 * k] - query[0];
		b[1] = polygon[2 * k + 1] - query[1];

		double cross = a[0] * b[1] - a[1] * b[0];
		if (base * cross < 0)
			return false;

		a[0] = b[0];
		a[1] = b[1];
	}

	return true;
}

template <typename Real>
Real compute_step(const Real *query, const Real *dir, const Real *start, const Real *stop)
{
	Real line[2] = { stop[0] - start[0], stop[1] - start[1] };
	Real across = dir[0] * line[1] - dir[1] * line[0];
	if (abs(across) < 1e-16) // parallel
		return Real(1e10); // infinite safe length

	Real temp[2] = { start[0] - query[0], start[1] - query[1] };
	Real bcross = temp[0] * line[1] - temp[1] * line[0];
	
	Real u = bcross / across;
	return u;
}

template <typename Real>
Real compute_safestep(const Real *polygon, int nb, const Real *query, const Real *dir)
{
	Real minsafe = Real(1e10);

	for (int k = 0; k < nb; ++k)
	{
		int next = (k + 1) % nb;
		
		const Real *p = &polygon[2 * k];
		const Real *q = &polygon[2 * next];

		Real step = compute_step(query, dir, p, q);
		if (step >= 0 && step < minsafe)
			minsafe = step;
	}

	return minsafe;
}

template <typename Real>
void triangle_barycoord(const Real *triangle, const Real *query, Real *coord)
{
	Real ab[2] = { triangle[2] - triangle[0], triangle[3] - triangle[1] };
	Real ac[2] = { triangle[4] - triangle[0], triangle[5] - triangle[1] };
	Real tarea = ab[0] * ac[1] - ab[1] * ac[0];

	Real pc[2] = { triangle[4] - query[0], triangle[5] - query[1] };
	Real pa[2] = { triangle[0] - query[0], triangle[1] - query[1] };
	Real parea = pc[0] * pa[1] - pc[1] * pa[0];
	coord[1] = parea / tarea;

	Real pb[2] = { triangle[2] - query[0], triangle[3] - query[1] };
	parea = pa[0] * pb[1] - pa[1] * pb[0];
	coord[2] = parea / tarea;

	coord[0] = 1.0 - coord[1] - coord[2];
}

template <typename Real, int DIM>
void barycoord_interpolation(const Real *attribs, const Real *coord, int nb, Real *result)
{
	for (int d = 0; d < DIM; ++d)
		result[d] = 0;

	for (int k = 0; k < nb; ++k)
	{
		for (int d = 0; d < DIM; ++d)
			result[d] += attribs[k * DIM + d] * coord[k];
	}
}


#endif
