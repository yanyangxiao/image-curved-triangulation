
#ifndef PIXELSET_H
#define PIXELSET_H

#include <vector>

struct PixelSet
{
	int ymin, ymax;
	std::vector<int> left, right;

	int xmin, xmax;
	std::vector<int> bottom, top;

	PixelSet() {}

	PixelSet& operator= (const PixelSet &rhs)
	{
		ymin = rhs.ymin;
		ymax = rhs.ymax;
		left = rhs.left;
		right = rhs.right;

		xmin = rhs.xmin;
		xmax = rhs.xmax;
		bottom = rhs.bottom;
		top = rhs.top;

		return *this;
	}
};

#endif
