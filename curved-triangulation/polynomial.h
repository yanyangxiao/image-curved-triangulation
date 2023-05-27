
/**
* author: Yanyang Xiao
* email : yanyangxiaoxyy@gmail.com
*/

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <math.h>
#include "../extern/Eigen/Eigen"

#include "pixelset.h"

template <typename Real>
class Polynomial
{
private:
	int                 _degree;
	std::vector<Real>   _coeff;

public:
	Polynomial(int d = 1)
	{
		set_degree(d);
	}

	Polynomial(const Polynomial &rhs)
	{
		_degree = rhs._degree;
		_coeff = rhs._coeff;
	}

	Polynomial& operator= (const Polynomial &rhs)
	{
		_degree = rhs._degree;
		_coeff = rhs._coeff;

		return *this;
	}

	int degree() const { return _degree; }
	void set_degree(int d) { _degree = d; }
	const std::vector<Real>* coefficients() const { return &_coeff; }

	template <typename T>
	void compute_factors(
		const T *image,
		int width,
		int height,
		int channel,
		const PixelSet* pixels);

	Real evaluate(int c, Real x, Real y) const;
	Real evaluate(int c, const Real* p) const
	{
		return evaluate(c, p[0], p[1]);
	}

	template <typename T>
	Real compute_energy(
		const T *image,
		int width,
		int height,
		int channel,
		const PixelSet* pixels,
		int Lp = 2) const;

protected:
	template <typename T>
	void compute_constant_factors(
		const T *image,
		int width,
		int height,
		int channel,
		const PixelSet* pixels);

	template <typename T>
	void compute_linear_factors(
		const T *image,
		int width,
		int height,
		int channel,
		const PixelSet* pixels);

	template <typename T>
	void compute_quadratic_factors(
		const T *image,
		int width,
		int height,
		int channel,
		const PixelSet* pixels);
};

template <typename Real>
template <typename T>
void Polynomial<Real>::compute_factors(
	const T *image,
	int width,
	int height,
	int channel,
	const PixelSet* pixels)
{
	if (!image || !pixels)
		return;

	switch (_degree)
	{
	case 1:
		compute_linear_factors(image, width, height, channel, pixels);
		break;
	case 2:
		compute_quadratic_factors(image, width, height, channel, pixels);
		break;
	default:
		compute_constant_factors(image, width, height, channel, pixels);
		break;
	}
}

template <typename Real>
template <typename T>
void Polynomial<Real>::compute_constant_factors(
	const T *image,
	int width,
	int height,
	int channel,
	const PixelSet* pixels)
{
	_coeff.clear();
	_coeff.resize(channel, Real(0.0));

	int count = 0;
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		int lineStart = j * width;

		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			int pixID = lineStart + i;
			const T* pixColor = &image[channel * pixID];

			for (int c = 0; c < channel; ++c)
			{
				_coeff[c] += Real(pixColor[c]);
			}

			++count;
		}
	}

	if (count > 0)
	{
		for (int c = 0; c < channel; ++c)
		{
			_coeff[c] /= count;
		}
	}
}

template <typename Real>
template <typename T>
void Polynomial<Real>::compute_linear_factors(
	const T *image,
	int width,
	int height,
	int channel,
	const PixelSet* pixels)
{
	_coeff.clear();
	_coeff.resize(channel * 3, Real(0.0));

	int pnb = 0;
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			++pnb;
		}
	}
	
	Eigen::MatrixXd matX(pnb, 3);
	std::vector<Eigen::VectorXd> vecX(channel, Eigen::VectorXd(pnb));
	
	pnb = 0;
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		Real y = (Real(j) + 0.5) / width;
		int lineStart = j * width;
		
		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			Real x = (Real(i) + 0.5) / width;

			matX(pnb, 0) = x;
			matX(pnb, 1) = y;
			matX(pnb, 2) = 1;

			int pixID = lineStart + i;
			const T* pixColor = &image[channel * pixID];
			
			for (int c = 0; c < channel; ++c)
			{
				vecX[c](pnb) = Real(pixColor[c]);
			}

			++pnb;
		}
	}

	Eigen::Matrix3d matA = matX.transpose() * matX;
	if (matA.determinant() != Real(0.0))
	{
		for (int c = 0; c < channel; ++c)
		{
			Eigen::Vector3d vecB = matX.transpose() * vecX[c];
			Eigen::Vector3d result = matA.colPivHouseholderQr().solve(vecB);

			_coeff[c * 3] = result(0);
			_coeff[c * 3 + 1] = result(1);
			_coeff[c * 3 + 2] = result(2);
		}
	}
	else
	{
		// back to constant
		for (int c = 0; c < channel; ++c)
		{
			Eigen::Vector3d vecB = matX.transpose() * vecX[c];

			_coeff[c * 3] = Real(0.0);
			_coeff[c * 3 + 1] = Real(0.0);
			_coeff[c * 3 + 2] = Real(0.0);
			if (matA(2, 2) != Real(0.0))
			{
				_coeff[c * 3 + 2] = vecB(2) / matA(2, 2);
			}
		}
	}
}

template <typename Real>
template <typename T>
void Polynomial<Real>::compute_quadratic_factors(
	const T *image,
	int width,
	int height,
	int channel,
	const PixelSet* pixels)
{
	_coeff.clear();
	_coeff.resize(channel * 6, Real(0.0));
	
	int pnb = 0;
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			++pnb;
		}
	}

	Eigen::MatrixXd matX(pnb, 6);
	std::vector<Eigen::VectorXd> vecX(channel, Eigen::VectorXd(pnb));

	pnb = 0;
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		Real y = (Real(j) + 0.5) / width;
		int lineStart = j * width;

		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			Real x = (Real(i) + 0.5) / width;

			matX(pnb, 0) = x * x;
			matX(pnb, 1) = x * y;
			matX(pnb, 2) = y * y;
			matX(pnb, 3) = x;
			matX(pnb, 4) = y;
			matX(pnb, 5) = 1;

			int pixID = lineStart + i;
			const T* pixColor = &image[channel * pixID];
			
			for (int c = 0; c < channel; ++c)
			{
				vecX[c](pnb) = Real(pixColor[c]);
			}

			++pnb;
		}
	}

	Eigen::MatrixXd matA = matX.transpose() * matX;
	if (matA.determinant() != Real(0.0))
	{
		for (int c = 0; c < channel; ++c)
		{
			Eigen::VectorXd vecB = matX.transpose() * vecX[c];
			Eigen::VectorXd result = matA.colPivHouseholderQr().solve(vecB);
			for (int k = 0; k < 6; ++k)
			{
				_coeff[c * 6 + k] = result(k);
			}
		}
	}
	else
	{
		// back to linear
		Eigen::Matrix3d matC;
		matC(0, 0) = matA(0, 5);
		matC(0, 1) = matA(1, 5);
		matC(0, 2) = matA(3, 5);
		matC(1, 1) = matA(2, 5);
		matC(0, 2) = matA(4, 5);
		matC(2, 2) = matA(5, 5);
		matC(1, 0) = matC(0, 1);
		matC(2, 0) = matC(0, 2);
		matC(2, 1) = matC(1, 2);

		if (matC.determinant() != Real(0.0))
		{
			for (int c = 0; c < channel; ++c)
			{
				_coeff[c * 6] = Real(0.0);
				_coeff[c * 6 + 1] = Real(0.0);
				_coeff[c * 6 + 2] = Real(0.0);

				Eigen::VectorXd vecB = matX.transpose() * vecX[c];
				Eigen::Vector3d vecD;
				for (int k = 0; k < 3; ++k)
					vecD(k) = vecB(k + 3);

				Eigen::Vector3d result = matC.colPivHouseholderQr().solve(vecD);

				for (int k = 0; k < 3; ++k)
					_coeff[c * 6 + k + 3] = result(k);
			}
		}
		else
		{
			// back to constant
			for (int c = 0; c < channel; ++c)
			{
				_coeff[c * 6] = Real(0.0);
				_coeff[c * 6 + 1] = Real(0.0);
				_coeff[c * 6 + 2] = Real(0.0);
				_coeff[c * 6 + 3] = Real(0.0);
				_coeff[c * 6 + 4] = Real(0.0);
				_coeff[c * 6 + 5] = Real(0.0);
				if (matA(5, 5) != Real(0.0))
				{
					Eigen::VectorXd vecB = matX.transpose() * vecX[c];
					_coeff[c * 6 + 5] = vecB(5) / matA(5, 5);
				}
			}
		}
	}
}

template <typename Real>
Real Polynomial<Real>::evaluate(int c, Real x, Real y) const
{
	Real result = Real(0.0);

	switch (_degree)
	{
	case 0:
		result = _coeff[c];
		break;
	case 1:
		result = _coeff[c * 3] * x + _coeff[c * 3 + 1] * y + _coeff[c * 3 + 2];
		break;
	case 2:
		result = _coeff[c * 6] * x * x + _coeff[c * 6 + 1] * x * y + _coeff[c * 6 + 2] * y * y
			+ _coeff[c * 6 + 3] * x + _coeff[c * 6 + 4] * y + _coeff[c * 6 + 5];
		break;
	}

	return result;
}

template <typename Real>
template <typename T>
Real Polynomial<Real>::compute_energy(
	const T *image,
	int width,
	int height,
	int channel,
	const PixelSet* pixels,
	int Lp = 2) const
{
	if (!image || !pixels)
		return Real(0.0);

	Real result = Real(0.0);
	for (int j = pixels->ymin; j <= pixels->ymax; ++j)
	{
		Real y = (Real(j) + 0.5) / width;
		int lineStart = j * width;

		for (int i = pixels->left[j]; i <= pixels->right[j]; ++i)
		{
			if (j < pixels->bottom[i] || j > pixels->top[i])
				continue;

			Real x = (Real(i) + 0.5) / width;

			int pixID = lineStart + i;
			const T* pixColor = &image[channel * pixID];

			Real tempEnergy = 0.0;
			for (int c = 0; c < channel; ++c)
			{
				Real approxVal = evaluate(c, x, y);
				Real pixVal = (Real)pixColor[c];
				Real absError = abs(pixVal - approxVal);
				tempEnergy += std::pow(absError, Lp);
			}

			result += tempEnergy;
		}
	}

	return result;
}

#endif