
#ifndef RASTERIZER_H
#define RASTERIZER_H

#include "pixelset.h"

template <typename Real>
class Rasterizer
{
protected:
	int      _width;
	int      _height;

public:
	Rasterizer(int w, int h) : _width(w), _height(h) { }
	void compute(const Real *polygon, int nb, PixelSet &pixels);
	static void rasterize(const int *xs, const int *ys, int nb, PixelSet &pixels);
};

template <typename Real>
void Rasterizer<Real>::compute(const Real *polygon, int nb, PixelSet &pixels)
{
	pixels.left.clear();
	pixels.right.clear();
	pixels.bottom.clear();
	pixels.top.clear();

	pixels.left.resize(_height, _width);
	pixels.right.resize(_height, -1);
	pixels.bottom.resize(_width, _height);
	pixels.top.resize(_width, -1);

	pixels.ymin = _height;
	pixels.ymax = -1;
	pixels.xmin = _width;
	pixels.xmax = -1;

	std::vector<int> xs(nb, -1);
	std::vector<int> ys(nb, -1);
	for (int k = 0; k < nb; ++k)
	{
		xs[k] = int(polygon[2 * k] * _width);
		ys[k] = int(polygon[2 * k + 1] * _width);

		if (xs[k] < 0) xs[k] = 0;
		if (xs[k] >= _width) xs[k] = _width - 1;
		if (ys[k] < 0) ys[k] = 0;
		if (ys[k] >= _height) ys[k] = _height - 1;

		pixels.ymin = pixels.ymin < ys[k] ? pixels.ymin : ys[k];
		pixels.ymax = pixels.ymax > ys[k] ? pixels.ymax : ys[k];

		pixels.xmin = pixels.xmin < xs[k] ? pixels.xmin : xs[k];
		pixels.xmax = pixels.xmax > xs[k] ? pixels.xmax : xs[k];
	}

	rasterize(&xs[0], &ys[0], nb, pixels);
}

template <typename Real>
void Rasterizer<Real>::rasterize(const int *xs, const int *ys, int nb, PixelSet &pixels)
{
	for (int k = 0; k < nb; ++k)
	{
		int next = (k + 1) % nb;
		int x1 = xs[k];
		int y1 = ys[k];
		int x2 = xs[next];
		int y2 = ys[next];

		if (y1 == y2)
		{
			int xmin = x1 < x2 ? x1 : x2;
			int xmax = x1 > x2 ? x1 : x2;

			pixels.left[y1] = pixels.left[y1] < xmin ? pixels.left[y1] : xmin;
			pixels.right[y1] = pixels.right[y1] > xmax ? pixels.right[y1] : xmax;

			for (int i = xmin; i <= xmax; ++i)
			{
				pixels.bottom[i] = pixels.bottom[i] < y1 ? pixels.bottom[i] : y1;
				pixels.top[i] = pixels.top[i] > y1 ? pixels.top[i] : y1;
			}

			continue;
		}

		// Bresenham algo.
		int dx = x2 - x1;
		int dy = y2 - y1;
		int sx = dx > 0 ? 1 : -1;
		int sy = dy > 0 ? 1 : -1;
		dx *= sx;
		dy *= sy;
		int x = x1;
		int y = y1;

		pixels.left[y] = pixels.left[y] < x ? pixels.left[y] : x;
		pixels.right[y] = pixels.right[y] > x ? pixels.right[y] : x;
		pixels.bottom[x] = pixels.bottom[x] < y ? pixels.bottom[x] : y;
		pixels.top[x] = pixels.top[x] > y ? pixels.top[x] : y;

		int e = dy - dx;
		while ((sy > 0 && y < y2) || (sy < 0 && y > y2))
		{
			while (e <= 0)
			{
				x += sx;
				e += 2 * dy;

				pixels.bottom[x] = pixels.bottom[x] < y ? pixels.bottom[x] : y;
				pixels.top[x] = pixels.top[x] > y ? pixels.top[x] : y;
			}

			y += sy;
			e -= 2 * dx;

			pixels.left[y] = pixels.left[y] < x ? pixels.left[y] : x;
			pixels.right[y] = pixels.right[y] > x ? pixels.right[y] : x;
		}

		while ((sx > 0 && x < x2) || (sx < 0 && x > x2))
		{
			x += sx;
			pixels.bottom[x] = pixels.bottom[x] < y ? pixels.bottom[x] : y;
			pixels.top[x] = pixels.top[x] > y ? pixels.top[x] : y;
		}

		pixels.left[y2] = pixels.left[y2] < x2 ? pixels.left[y2] : x2;
		pixels.right[y2] = pixels.right[y2] > x2 ? pixels.right[y2] : x2;
		pixels.bottom[x2] = pixels.bottom[x2] < y2 ? pixels.bottom[x2] : y2;
		pixels.top[x2] = pixels.top[x2] > y2 ? pixels.top[x2] : y2;
	}
}

#endif
