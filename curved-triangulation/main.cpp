
/**
* Copyright Xiao, Yanyang (yanyangxiaoxyy@gmail.com)
* This file is part of the project "Curved Optimal Triangulation"
*
* author: Xiao, Yanyang
* email : yanyangxiaoxyy@gmail.com
*
* the main program
*/

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "../extern/stb/stb_image.h"
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../extern/stb/stb_image_write.h"
#endif

#include <fstream>
#include "curved-triangulation.h"
#include "../xlog.h"
#include "../timer.h"

bool load_features(const char *filename, std::vector<std::vector<int>> &features)
{
	features.clear();

	xlog("%s", filename);

	std::ifstream file(filename);
	if (!file)
	{
		xlog("failed");
		return false;
	}

	int nb;
	file >> nb;

	for (int k = 0; k < nb; ++k)
	{
		std::vector<int> path;
		int len = 0;
		int p[2];

		file >> len;
		for (int i = 0; i < len; ++i)
		{
			file >> p[0] >> p[1];
			path.push_back(p[0]);
			path.push_back(p[1]);
		}

		if (!path.empty())
			features.push_back(path);
	}

	file.close();
	return true;
}

int main(int argc, char *argv[])
{
	// demo.exe degree regweight input.jpg init.om features.txt vnb iters
	
	// 1 ----- set parameters
	xlog("1 ----- parameters");

	int degree = atoi(argv[1]);
	double regweight = atof(argv[2]);

	CurvedTriangulation cot;
	cot.set_degree(degree);
	cot.set_regweight(regweight);

	xlog("degree = %d, regularization weight = %f", degree, regweight);

	// 2 ----- load image
	xlog("2 ----- load image");

	int width = -1;
	int height = -1;
	int channel = -1;
	
	stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load(argv[3], &width, &height, &channel, 0);
	xlog("input image: width = %d, height = %d, channel = %d", width, height, channel);	
	
	cot.set_image(image, width, height, channel);

	// 3 ----- load mesh
	xlog("3 ----- load mesh");

	if (!cot.load_mesh(argv[4]))
	{
		// load failed, let us generate

		// 3 ----- init triangulation
		xlog("3 ----- load failed, init mesh");

		std::vector<std::vector<int>> features;
		load_features(argv[5], features);
		xlog("input features: %d paths", (int)features.size());

		int vnb = atoi(argv[6]);

		//cot.random_init(vnb);
		//cot.greedy_init(vnb);
		cot.features_init(features, vnb);

		cot.save_mesh("3-init-mesh.om");
	}

	cot.svg_mesh("3-init-mesh.svg");

	xlog("self-intersection edges: %d", (int)cot.check_self_intersection());

	// 4 ----- approx data of init
	xlog("4 ----- approx data of init");

	cot.compute_approx();
	cot.svg_approx("4-init-approx.svg");

	// 5 ----- optimize mesh
	xlog("5 ----- optimize mesh");

	int iter = atoi(argv[7]);
	cot.optimize(iter, true);
	cot.save_mesh("5-result-mesh.om");
	cot.svg_mesh("5-result-mesh.svg");

	xlog("self-intersection edges: %d", (int)cot.check_self_intersection());

	// 6 ----- save results
	xlog("6 ----- approx data of result");
	cot.svg_approx("6-result-approx.svg");
	
	// 7 ----- free resources
	xlog("7 ----- free resources");
	stbi_image_free(image);

	return 0;
}

