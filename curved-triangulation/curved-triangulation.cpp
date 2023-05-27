
#include "curved-triangulation.h"

CurvedTriangulation::CurvedTriangulation()
{
	_image = NULL;
	_width = 0;
	_height = 0;
	_channel = 0;
	_pixwidth = 0;
	
	_mesh.add_property(FPropData, "fprop_data");
	_mesh.add_property(EPropControl, "eprop_control");
	_mesh.add_property(EPropFeature, "eprop_feature");
	_mesh.add_property(VPropFeature, "vprop_feature");

	_mesh.property(EPropControl).set_persistent(true);
	_mesh.property(EPropFeature).set_persistent(true);
	_mesh.property(VPropFeature).set_persistent(true);

	set_degree(1);
	set_regweight(0.001);
}

CurvedTriangulation::~CurvedTriangulation()
{
	_mesh.remove_property(FPropData);
	_mesh.remove_property(EPropControl);
	_mesh.remove_property(EPropFeature);
	_mesh.remove_property(VPropFeature);
	_mesh.clear();

	if (_image)
	{
		free(_image);
		_image = NULL;
	}
}

void CurvedTriangulation::set_image(const unsigned char *image, int width, int height, int channel)
{
	if (_image)
	{
		free(_image);
		_image = NULL;
	}

	_width = width;
	_height = height;
	_channel = channel;

	_pixwidth = 1.0 / _width;

	_image = (float *)malloc(sizeof(float) * width * height * channel);
	for (int j = 0; j < height; ++j)
	{
		for (int i = 0; i < width; ++i)
		{
			int idx = j * width + i;
			for (int k = 0; k < channel; ++k)
				_image[channel * idx + k] = float(image[channel * idx + k]) / 255.0f;
		}
	}
}

