#include <CL/cl.hpp>
#include <string>
#include <stdio.h>
#include <iostream>
using namespace std;

extern "C" {
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
}

#include "systems.hpp"

int merge_test() {
	Face f1, f2;
	f1.size = 3;
	f1.subp = new IntVector3[f1.size] { IntVector3(-2, 2, 3), IntVector3(-4, 3, 1), IntVector3(0, 0, 0)};
	f1.point = Vector3(-1, 0.5, 2);
	f1.normal = Normalize(Vector3(-5, 3, -4));
	f1.box = IntAlignedBox3(IntVector3(-5, 0, 0), IntVector3(0, 3, 3));
	f1.certainty = 1.0F;

	//f2.size = 5;


	return 0;
}

int depth_test() {
    int x, y, n;
    unsigned int* image1 = reinterpret_cast<unsigned int*>(stbi_load("MCImg1.png", &x, &y, &n, 4));
    unsigned int* image2 = reinterpret_cast<unsigned int*>(stbi_load("MCImg2.png", &x, &y, &n, 4));

	Gen_State gen = Gen_State();
	Cam_State cam(x, y - 128 * 4, 0.57735F, 1);
	DP_State dps(&gen, &cam, 16, 1, 1);
	IFE_State ife(&gen, &dps);
	PMU_State pmu(&ife);
	std::cout << GetErrors(&dps) << std::endl;
	std::cout << GetErrors(&ife) << std::endl;
	std::cout << Stringify(&gen) << std::endl;
	std::cout << Stringify(&dps) << std::endl;
	std::cout << Stringify(&ife) << std::endl;
	std::cout << "Image halves: " + std::to_string(dps.image_halves) << std::endl;

	dps.SetImages(image1, image2);
	std::cout << GetErrors(&dps) << std::endl;
	
	dps.Execute();
	std::cout << GetErrors(&dps) << std::endl;
	
	ife.Execute();
	std::cout << GetErrors(&ife) << std::endl;

	pmu.Execute();

	//pmu.SetPose(Transform3(MakeRotation(1, -0.75), Vector3(3, -20, -7)));

	/*
	Frustum frust;
	float fhw = 3;
	float fz = 5;
	frust.p[0] = Vector3(fhw, fhw, fz);
	frust.p[1] = Vector3(-fhw, fhw, fz);
	frust.p[2] = Vector3(-fhw, -fhw, fz);
	frust.p[3] = Vector3(fhw, -fhw, fz);
	frust.p[4] = Vector3(0, 0, 0);
	pmu.frustum = frust;

	Vector3 plane_point(0.5, 0.5, 0.5);
	//pmu.plane_map->planes = new PlaneCube;
	//pmu.plane_map->planes->c = 1.0F;
	for (int i = 0; i < 20; i++) Insert(pmu.plane_map, pmu.log2_pmwidth, pmu.pm_disp, Plane3(Normalize(Vector3(-5, 3, -4)), plane_point),
		plane_point, Vector3(1, 1, 1), 0.5);
	*/

	std::cout << Stringify(*(pmu.plane_map), "") << std::endl;
	std::cout << pmu.pose.R << std::endl;
	std::cout << pmu.pose.t << std::endl;
	std::cout << Stringify(pmu.pose * pmu.frustum) << std::endl;
	std::cout << (1 << pmu.log2_pmwidth) << std::endl;
	std::cout << pmu.pm_disp << std::endl;
	std::cout << GetPlaneMapBox(&pmu).minV << std::endl;
	std::cout << GetPlaneMapBox(&pmu).maxV << std::endl;

	/*
	pmu.Shrink();

	std::cout << Stringify(*(pmu.plane_map)) << std::endl;
	std::cout << Stringify(pmu.pose * pmu.frustum) << std::endl;
	std::cout << (1 << pmu.log2_pmwidth) << std::endl;
	std::cout << GetPlaneMapBox(&pmu).minV << std::endl;
	std::cout << GetPlaneMapBox(&pmu).maxV << std::endl;
	*/
	
	stbi_image_free(image1);
	stbi_image_free(image2);
	return 0;
}

int main() {
	return merge_test();
}