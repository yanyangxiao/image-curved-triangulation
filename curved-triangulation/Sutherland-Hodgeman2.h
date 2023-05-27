
#ifndef SUTHERLAND_HODGEMAN2_H
#define SUTHERLAND_HODGEMAN2_H

#include <vector>

template <typename Real>
Real plane_side2(const Real *point, const Real *normal, const Real *query)
{
	Real result = (query[0] - point[0]) * normal[0] + (query[1] - point[1]) * normal[1];
	return result;
}

template <typename Real>
bool Sutherland_Hodgeman2(
	const Real *point, const Real *normal, 
	const std::vector<Real> &polygon, std::vector<Real> &result)
{
	result.clear();	

	int nb = (int)polygon.size() / 2;
	
	const Real *p = &polygon[2 * (nb - 1)];
	Real pside = plane_side2(point, normal, p);

	bool clipped = false;
	for (int i = 0; i < nb; ++i)
	{
		const Real *q = &polygon[2 * i];
		Real qside = plane_side2(point, normal, q);

		Real sign = pside * qside;
		if (sign < Real(0.0)) // intersect!
		{
			Real ratio = abs(pside) / (abs(pside) + abs(qside));

			Real x[2];
			x[0] = p[0] + ratio * (q[0] - p[0]);
			x[1] = p[1] + ratio * (q[1] - p[1]);

			if (pside < Real(0.0))
			{
				result.push_back(p[0]);
				result.push_back(p[1]);

				result.push_back(x[0]);
				result.push_back(x[1]);
			}
			else
			{
				result.push_back(x[0]);
				result.push_back(x[1]);
			}

			clipped = true;
		}
		else if (sign > Real(0.0))
		{
			if (pside < Real(0.0))
			{
				result.push_back(p[0]);
				result.push_back(p[1]);
			}
		}
		else // sign == 0.0
		{
			if (pside <= Real(0.0))
			{
				result.push_back(p[0]);
				result.push_back(p[1]);
			}
		}

		p = q;
		pside = qside;
	}

	return clipped;
}

#endif
