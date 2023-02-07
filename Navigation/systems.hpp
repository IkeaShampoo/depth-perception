#pragma once
#include <CL/cl.hpp>
#include <string>
#include <stdio.h>
#include <iostream>
using namespace std;

#include "kernels.hpp"
#include "geometry.hpp"
using namespace gmtry3D;
using namespace gmtry3DInts;
using namespace gmtry2D;

#define SQRT3 (1.7320508756)

#define PASSTHROUGHCOEFFICIENT (0x1.p-6)
#define NEWPLANECERTAINTYSCALE (1) //0.0625
#define COLOROPTIMIZATIONRATE (1)
#define PRUNERATE (0.125)

unsigned int log2(unsigned int n, bool ceil) { //sign of exponent is positive
	unsigned int i = 0; while (1 << i < n) i++; if ((1 << i == n) || ceil) return i; else return i - 1; }

unsigned int sqrt(unsigned int n, bool ceil) {
	unsigned int i = 0; while (i * i < n) i++; if ((i * i == n) || ceil) return i; else return i - 1; }

inline float square(float n) { return n * n; }

int fit_units(int unit, int step, int space) { return 1 + (space - unit) / step; }

int ceil_div(int a, int b) { return a / b + (a % b != 0); }

int make_even_div(int a, int b) { return ceil_div(a, b) * b; }

template<typename T> int sign(T value) { return (value > 0) - (value < 0); }

const char* getErrorString(cl_int error)
{
	switch (error) {
		// run-time and JIT compiler errors
	case 0: return "CL_SUCCESS";
	case -1: return "CL_DEVICE_NOT_FOUND";
	case -2: return "CL_DEVICE_NOT_AVAILABLE";
	case -3: return "CL_COMPILER_NOT_AVAILABLE";
	case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case -5: return "CL_OUT_OF_RESOURCES";
	case -6: return "CL_OUT_OF_HOST_MEMORY";
	case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
	case -8: return "CL_MEM_COPY_OVERLAP";
	case -9: return "CL_IMAGE_FORMAT_MISMATCH";
	case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case -11: return "CL_BUILD_PROGRAM_FAILURE";
	case -12: return "CL_MAP_FAILURE";
	case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
	case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
	case -15: return "CL_COMPILE_PROGRAM_FAILURE";
	case -16: return "CL_LINKER_NOT_AVAILABLE";
	case -17: return "CL_LINK_PROGRAM_FAILURE";
	case -18: return "CL_DEVICE_PARTITION_FAILED";
	case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

		// compile-time errors
	case -30: return "CL_INVALID_VALUE";
	case -31: return "CL_INVALID_DEVICE_TYPE";
	case -32: return "CL_INVALID_PLATFORM";
	case -33: return "CL_INVALID_DEVICE";
	case -34: return "CL_INVALID_CONTEXT";
	case -35: return "CL_INVALID_QUEUE_PROPERTIES";
	case -36: return "CL_INVALID_COMMAND_QUEUE";
	case -37: return "CL_INVALID_HOST_PTR";
	case -38: return "CL_INVALID_MEM_OBJECT";
	case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case -40: return "CL_INVALID_IMAGE_SIZE";
	case -41: return "CL_INVALID_SAMPLER";
	case -42: return "CL_INVALID_BINARY";
	case -43: return "CL_INVALID_BUILD_OPTIONS";
	case -44: return "CL_INVALID_PROGRAM";
	case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
	case -46: return "CL_INVALID_KERNEL_NAME";
	case -47: return "CL_INVALID_KERNEL_DEFINITION";
	case -48: return "CL_INVALID_KERNEL";
	case -49: return "CL_INVALID_ARG_INDEX";
	case -50: return "CL_INVALID_ARG_VALUE";
	case -51: return "CL_INVALID_ARG_SIZE";
	case -52: return "CL_INVALID_KERNEL_ARGS";
	case -53: return "CL_INVALID_WORK_DIMENSION";
	case -54: return "CL_INVALID_WORK_GROUP_SIZE";
	case -55: return "CL_INVALID_WORK_ITEM_SIZE";
	case -56: return "CL_INVALID_GLOBAL_OFFSET";
	case -57: return "CL_INVALID_EVENT_WAIT_LIST";
	case -58: return "CL_INVALID_EVENT";
	case -59: return "CL_INVALID_OPERATION";
	case -60: return "CL_INVALID_GL_OBJECT";
	case -61: return "CL_INVALID_BUFFER_SIZE";
	case -62: return "CL_INVALID_MIP_LEVEL";
	case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
	case -64: return "CL_INVALID_PROPERTY";
	case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
	case -66: return "CL_INVALID_COMPILER_OPTIONS";
	case -67: return "CL_INVALID_LINKER_OPTIONS";
	case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

		// extension errors
	case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
	case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
	case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
	case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
	case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
	case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
	default: return "Unknown OpenCL error";
	}
}

void PrintError(cl_int error) {
	cout << getErrorString(error) << endl;
}

int* Visualize(const int* shifts, const float* certainties, int length, float scalar) {
	float color0[3] = { 1.0F, 0.0F, 0.0F };
	float color1[3] = { 0.0F, 1.0F, 0.0F };
	float coldisp[3] = { 
		scalar * (color0[0] - color1[0]), 
		scalar * (color0[1] - color1[1]), 
		scalar * (color0[2] - color1[2])
	};
	int* image = new int[length];
	int r, g, b;
	for (int i = 0; i < length; i++) {
		r = certainties[i] * (color0[0] + static_cast<float>(shifts[i]) * coldisp[0]);
		r += (r > 255) * (255 - r);
		g = certainties[i] * (color0[1] + static_cast<float>(shifts[i]) * coldisp[1]);
		g += (g > 255) * (255 - g);
		b = certainties[i] * (color0[2] + static_cast<float>(shifts[i]) * coldisp[2]);
		b += (b > 255) * (255 - b);
		image[i] = (b << 16) + (g << 8) + r;
	}
	return image;
}

struct Gen_State {
	cl_int* err;
	int err_len;
	cl_device_id device;
	cl_platform_id platform;
	cl_context context;
	cl_command_queue queue;
	size_t wg_max_dims[3];
	Gen_State() {
		int err_idx = 0;
		err = new cl_int[err_len = 5];
		err[err_idx++] = clGetPlatformIDs(1, &platform, NULL);
		err[err_idx++] = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, NULL);
		context = clCreateContext(NULL, 1, &device, NULL, NULL, &(err[err_idx++]));
		queue = clCreateCommandQueueWithProperties(context, device, NULL, &(err[err_idx++]));
		err[err_idx++] = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, &wg_max_dims, NULL);
	}
	~Gen_State() {
		delete[] err;
		err = new cl_int[1];
		err[0] = clFlush(queue);
		err[0] = clFinish(queue);
		err[0] = clReleaseCommandQueue(queue);
		err[0] = clReleaseContext(context);
		err[0] = clReleaseDevice(device);
		delete[] err;
	}
};

std::string Stringify(Gen_State* gs) {
	char data[1024] = "";
	delete[] gs->err;
	gs->err = new cl_int[gs->err_len = 8];
	int err_idx = 0;

	std::string str = "General State: ";
	gs->err[err_idx++] = clGetDeviceInfo(gs->device, CL_DEVICE_NAME, 1024, data, NULL);
	str += "Device name: " + std::string(data) + '\n';
	str += "Device work group max dims: " + 
		std::to_string(gs->wg_max_dims[0]) + ", " + 
		std::to_string(gs->wg_max_dims[1]) + ", " +
		std::to_string(gs->wg_max_dims[2]) + '\n';
	for (int i = 0; i < 1024; i++) data[i] = 0;
	gs->err[err_idx++] = clGetPlatformInfo(gs->platform, CL_PLATFORM_NAME, 1024, data, NULL);
	str += "Platform name: " + std::string(data) + '\n';
	return str;
}

struct Cam_State {
	int res_dims[2];
	float tangent, spacing;
	Cam_State(int width, int height, float tan_fov, float cam_space) {
		res_dims[0] = width;
		res_dims[1] = height;
		tangent = tan_fov;
		spacing = cam_space;
	}
};

struct DP_State {
	const Gen_State* gen_state;
	const Cam_State* cam_state;

	cl_int* err;
	int err_len;
	cl_int build_err;
	char build_log[1024] = "";
	cl_program program;
	cl_kernel halver_kernel, loss_kernel, shift_kernel, vertex_kernel;
	size_t halver_lsize[2], loss_lsize[2], vertex_lsize[2], vertex_gsize[2];
	size_t halver_wg_max_size, loss_wg_max_size, shift_wg_max_size, vertex_wg_max_size;
	cl_mem image1, image2, vcolors, losses, shifts, certainties, vertices;
	int image_dims[2], losses_dims[2], shifts_dims[2], vcolors_dims[2];
	int image_buf_size, vcolors_buf_size, losses_buf_size, shifts_buf_size, certainties_buf_size, vertices_buf_size;
	int log2_bwidth, step_size, image_halves, cushion;
	float vertex_constants[4];
	DP_State(const Gen_State* general_state, const Cam_State* camera_state, int batch_width, int batch_step_size, int end_cushion) {
		int err_idx = 0;
		err = new cl_int[err_len = 39 - 1];
		gen_state = general_state;
		cam_state = camera_state;

		//PROGRAM & KERNEL
		string programString = 
			img_halver_kernel_string + " " + 
			dp_loss_kernel_string + " " + 
			dp_shift_kernel_string + " " +
			dp_vertex_kernel_string;
		const char* programChars = programString.c_str();
		size_t programStringSize = programString.size();
		program = clCreateProgramWithSource(gen_state->context, 1, &programChars, &programStringSize, &(err[err_idx++]));
		build_err = err[err_idx++] = clBuildProgram(program, 1, &(gen_state->device), NULL, NULL, NULL);
		//err[err_idx++] = clGetProgramBuildInfo(program, gen_state->device, CL_PROGRAM_BUILD_LOG, 1024, build_log, NULL);
		if (build_err) return;
		halver_kernel = clCreateKernel(program, "halve_img", &(err[err_idx++]));
		loss_kernel = clCreateKernel(program, "get_losses", &(err[err_idx++]));
		shift_kernel = clCreateKernel(program, "get_shifts", &(err[err_idx++]));
		vertex_kernel = clCreateKernel(program, "get_vertices", &(err[err_idx++]));
		err[err_idx++] = clGetKernelWorkGroupInfo(halver_kernel, gen_state->device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &halver_wg_max_size, NULL);
		err[err_idx++] = clGetKernelWorkGroupInfo(loss_kernel, gen_state->device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &loss_wg_max_size, NULL);
		err[err_idx++] = clGetKernelWorkGroupInfo(shift_kernel, gen_state->device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &shift_wg_max_size, NULL);
		err[err_idx++] = clGetKernelWorkGroupInfo(vertex_kernel, gen_state->device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &vertex_wg_max_size, NULL);

		//HALVER LOCAL SIZE
		halver_lsize[1] = sqrt(halver_wg_max_size, false);
		halver_lsize[0] = halver_wg_max_size / halver_lsize[1];

		//BATCH OPTIMIZATION
		if (batch_width * batch_width > loss_wg_max_size)
			log2_bwidth = log2(loss_wg_max_size, false) / 2;
		else log2_bwidth = log2(batch_width, false);
		cushion = end_cushion;

		//LOSS LOCAL SIZE
		loss_lsize[1] = sqrt(loss_wg_max_size, false);
		loss_lsize[0] = loss_wg_max_size / loss_lsize[1];

		//IMAGE CONFIGURATION
		step_size = batch_step_size;
		image_halves = 0;
		while (fit_units(1 << log2_bwidth, step_size, cam_state->res_dims[0] >> image_halves) > shift_wg_max_size) image_halves++;
		if (fit_units(1 << log2_bwidth, step_size, cam_state->res_dims[0] >> image_halves) < shift_wg_max_size && image_halves > 0) image_halves--;
		image_dims[0] = cam_state->res_dims[0] >> image_halves; image_dims[1] = cam_state->res_dims[1] >> image_halves;
		image_buf_size = image_dims[0] * image_dims[1] * sizeof(cl_uchar4);
		image1 = clCreateBuffer(gen_state->context, CL_MEM_READ_ONLY, image_buf_size, NULL, &(err[err_idx++]));
		image2 = clCreateBuffer(gen_state->context, CL_MEM_READ_ONLY, image_buf_size, NULL, &(err[err_idx++]));
		
		//BUFFER CONFIGURATION
		losses_dims[0] = fit_units(1 << log2_bwidth, step_size, image_dims[0]); 
		if (losses_dims[0] > shift_wg_max_size * 2) losses_dims[0] = shift_wg_max_size * 2;
		losses_dims[1] = image_dims[1] >> log2_bwidth;
		losses_buf_size = losses_dims[0] * losses_dims[1] * sizeof(cl_float);
		shifts_dims[0] = (image_dims[0] >> log2_bwidth) - cushion; shifts_dims[1] = losses_dims[1];
		vcolors_dims[0] = image_dims[0] >> log2_bwidth; vcolors_dims[1] = shifts_dims[1];
		shifts_buf_size = shifts_dims[0] * shifts_dims[1] * sizeof(cl_int);
		certainties_buf_size = shifts_dims[0] * shifts_dims[1] * sizeof(cl_float);
		vertices_buf_size = shifts_dims[0] * shifts_dims[1] * sizeof(cl_float4);
		vcolors_buf_size = vcolors_dims[0] * vcolors_dims[1] * sizeof(cl_uchar4);
		losses = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, losses_buf_size, NULL, &(err[err_idx++]));
		shifts = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, shifts_buf_size, NULL, &(err[err_idx++]));
		certainties = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, certainties_buf_size, NULL, &(err[err_idx++]));
		vertices = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, vertices_buf_size, NULL, &(err[err_idx++]));
		vcolors = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, vcolors_buf_size, NULL, &(err[err_idx++]));

		//VERTEX LOCAL & GLOBAL SIZE
		vertex_lsize[0] = sqrt(vertex_wg_max_size, false);
		vertex_lsize[1] = vertex_wg_max_size / vertex_lsize[0];
		vertex_gsize[0] = make_even_div(shifts_dims[0], vertex_lsize[0]);
		vertex_gsize[1] = make_even_div(shifts_dims[1], vertex_lsize[1]);

		//VERTEX CONSTANTS
		const float uni_vertex_constant = (2.0F * cam_state->tangent) /
			(image_dims[0] + (image_dims[1] > image_dims[0]) * (image_dims[1] - image_dims[0]));
		const int actual_batch_width = 1 << log2_bwidth;
		vertex_constants[0] = cam_state->spacing / (uni_vertex_constant * step_size);
		vertex_constants[1] = uni_vertex_constant * actual_batch_width;
		vertex_constants[2] = 0.5F - shifts_dims[0] / 2.0F;
		vertex_constants[3] = 0.5F - shifts_dims[1] / 2.0F;

		//LOSS KERNEL ARGUMENT
		err[err_idx++] = clSetKernelArg(loss_kernel, 0, sizeof(&image1), &image1);
		err[err_idx++] = clSetKernelArg(loss_kernel, 1, sizeof(&image2), &image2);
		err[err_idx++] = clSetKernelArg(loss_kernel, 2, sizeof(cl_int), &(image_dims[0]));
			//SKIP 3 ("losses_width")
		err[err_idx++] = clSetKernelArg(loss_kernel, 4, sizeof(cl_int), &(losses_dims[1]));
		err[err_idx++] = clSetKernelArg(loss_kernel, 5, sizeof(cl_int), &log2_bwidth);
			//SKIP 6 ("image2_col")
		err[err_idx++] = clSetKernelArg(loss_kernel, 7, sizeof(cl_int), &step_size);
		err[err_idx++] = clSetKernelArg(loss_kernel, 8, 
			loss_lsize[0] * loss_lsize[1] * sizeof(cl_float), NULL);
		err[err_idx++] = clSetKernelArg(loss_kernel, 9, sizeof(&losses), &losses);

		//SHIFT KERNEL ARGUMENT
		err[err_idx++] = clSetKernelArg(shift_kernel, 0, sizeof(&losses), &losses);
		err[err_idx++] = clSetKernelArg(shift_kernel, 1, sizeof(cl_int), &(shifts_dims[0]));
		err[err_idx++] = clSetKernelArg(shift_kernel, 2, sizeof(cl_int), &(shifts_dims[1]));
			//SKIP 3 ("losses_width")
			//SKIP 4 ("image2_col")
			//SKIP 5 ("loss_min")
			//SKIP 6 ("loss_min_idx")
			//SKIP 7 ("loss_sum")
		err[err_idx++] = clSetKernelArg(shift_kernel, 8, sizeof(&shifts), &shifts);
		err[err_idx++] = clSetKernelArg(shift_kernel, 9, sizeof(&certainties), &certainties);

		//VERTEX KERNEL ARGUMENT
		err[err_idx++] = clSetKernelArg(vertex_kernel, 0, sizeof(&shifts), &shifts);
		err[err_idx++] = clSetKernelArg(vertex_kernel, 1, sizeof(cl_int), &(shifts_dims[0]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 2, sizeof(cl_int), &(shifts_dims[1]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 3, sizeof(cl_float), &(vertex_constants[0]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 4, sizeof(cl_float), &(vertex_constants[1]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 5, sizeof(cl_float), &(vertex_constants[2]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 6, sizeof(cl_float), &(vertex_constants[3]));
		err[err_idx++] = clSetKernelArg(vertex_kernel, 7, sizeof(&vertices), &vertices);
	}
	void SetImages(const unsigned int* img_og1, const unsigned int* img_og2) {
		delete[] err;
		int err_idx = 0;

		//STRAIGHT WRITE SOURCE IMAGE
		if (image_halves == 0) {
			err = new cl_int[err_len = 3];
			cl_event write_events[2];
			err[err_idx++] = clEnqueueWriteBuffer(gen_state->queue, image1, false, NULL, image_buf_size, img_og1, 0, NULL, &(write_events[0]));
			err[err_idx++] = clEnqueueWriteBuffer(gen_state->queue, image2, false, NULL, image_buf_size, img_og2, 0, NULL, &(write_events[1]));
			err[err_idx++] = clWaitForEvents(2, write_events);
			return;
		}
		else err = new cl_int[err_len = 6 + log2_bwidth * 6 + 2 * (4 + image_halves * 6)];

		//TEMP IMAGE BUFFER INITIALIZATION
		int larger_img_dims[2], smaller_img_dims[2];
		const int img_buf_sizes[2] = {
			cam_state->res_dims[0] * cam_state->res_dims[1] * sizeof(cl_uchar4),
			(cam_state->res_dims[0] >> 1) * (cam_state->res_dims[1] >> 1) * sizeof(cl_uchar4)
		};
		cl_mem img_bufs[2] = {
			clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, img_buf_sizes[0], NULL, &(err[err_idx++])),
			clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, img_buf_sizes[1], NULL, &(err[err_idx++]))
		};

		//GROUPING BUFFERS
		cl_mem* dst_img_bufs[2] = { &image1, &image2 };
		const unsigned int* src_img_bufs[2] = { img_og1, img_og2 };

		//MISC INITIALIZATION
		cl_event write_complete, halver_complete, copy_complete;
		size_t halver_gsize[2];
		bool larger_img_idx;

		for (int i = 0; i < 2; i++) {
			larger_img_idx = 0;

			//IMAGE DIMENSION INITIALIZATION
			larger_img_dims[0] = cam_state->res_dims[0]; larger_img_dims[1] = cam_state->res_dims[1];
			smaller_img_dims[0] = larger_img_dims[0] >> 1; smaller_img_dims[1] = larger_img_dims[1] >> 1;

			//WRITE SOURCE IMAGE TO BUFFER
			err[err_idx++] = clEnqueueWriteBuffer(gen_state->queue, img_bufs[larger_img_idx], false, 
				0, img_buf_sizes[larger_img_idx], src_img_bufs[i], 0, NULL, &write_complete);
			err[err_idx++] = clWaitForEvents(1, &write_complete);

			for (int j = 0; j < image_halves; j++) {
				//HALVER GLOBAL SIZE
				halver_gsize[0] = make_even_div(smaller_img_dims[0], halver_lsize[0]);
				halver_gsize[1] = make_even_div(smaller_img_dims[1], halver_lsize[1]);

				//HALVER KERNEL LAUNCH
				err[err_idx++] = clSetKernelArg(halver_kernel, 0, sizeof(&(img_bufs[larger_img_idx])), &(img_bufs[larger_img_idx]));
				err[err_idx++] = clSetKernelArg(halver_kernel, 1, sizeof(&(img_bufs[!larger_img_idx])), &(img_bufs[!larger_img_idx]));
				err[err_idx++] = clSetKernelArg(halver_kernel, 2, sizeof(cl_int), &(larger_img_dims[0]));
				err[err_idx++] = clSetKernelArg(halver_kernel, 3, sizeof(cl_int), &(larger_img_dims[1]));
				err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, halver_kernel, 2, NULL, halver_gsize, halver_lsize, 0, NULL, &halver_complete);
				err[err_idx++] = clWaitForEvents(1, &halver_complete);

				//BUFFER SWAP
				larger_img_idx = !larger_img_idx;

				//IMAGE HALF DIMENSION
				larger_img_dims[0] = smaller_img_dims[0]; larger_img_dims[1] = smaller_img_dims[1];
				smaller_img_dims[0] = larger_img_dims[0] >> 1; smaller_img_dims[1] = larger_img_dims[1] >> 1;
			}
			//COPY REDUCED IMAGE TO DESTINATION
			err[err_idx++] = clEnqueueCopyBuffer(gen_state->queue, img_bufs[larger_img_idx], *(dst_img_bufs[i]), 
				0, 0, image_buf_size, 0, NULL, &copy_complete);
			err[err_idx++] = clWaitForEvents(1, &copy_complete);
		}

		for (int i = 0; i < log2_bwidth; i++) {
			//HALVER GLOBAL SIZE
			halver_gsize[0] = make_even_div(smaller_img_dims[0], halver_lsize[0]);
			halver_gsize[1] = make_even_div(smaller_img_dims[1], halver_lsize[1]);

			//HALVER KERNEL LAUNCH
			err[err_idx++] = clSetKernelArg(halver_kernel, 0, sizeof(&(img_bufs[larger_img_idx])), &(img_bufs[larger_img_idx]));
			err[err_idx++] = clSetKernelArg(halver_kernel, 1, sizeof(&(img_bufs[!larger_img_idx])), &(img_bufs[!larger_img_idx]));
			err[err_idx++] = clSetKernelArg(halver_kernel, 2, sizeof(cl_int), &(larger_img_dims[0]));
			err[err_idx++] = clSetKernelArg(halver_kernel, 3, sizeof(cl_int), &(larger_img_dims[1]));
			err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, halver_kernel, 2, NULL, halver_gsize, halver_lsize, 0, NULL, &halver_complete);
			err[err_idx++] = clWaitForEvents(1, &halver_complete);

			//BUFFER SWAP
			larger_img_idx = !larger_img_idx;

			//IMAGE HALF DIMENSION
			larger_img_dims[0] = smaller_img_dims[0]; larger_img_dims[1] = smaller_img_dims[1];
			smaller_img_dims[0] = larger_img_dims[0] >> 1; smaller_img_dims[1] = larger_img_dims[1] >> 1;
		}
		err[err_idx++] = clEnqueueCopyBuffer(gen_state->queue, img_bufs[larger_img_idx], vcolors, 0, 0, vcolors_buf_size, 0, NULL, &copy_complete);
		err[err_idx++] = clWaitForEvents(1, &copy_complete);

		err[err_idx++] = clReleaseMemObject(img_bufs[0]);
		err[err_idx++] = clReleaseMemObject(img_bufs[1]);
	}
	void Execute() {
		delete[] err;
		err = new cl_int[err_len = 2 + shifts_dims[0] * 11];
		int err_idx = 0;

		cl_event loss_complete, shift_complete, vertex_complete;

		//LOSS GLOBAL HEIGHT & DEPTH
		size_t loss_gsize[2], shift_lsize[2], shift_gsize[2];
		loss_gsize[1] = make_even_div(losses_dims[1] << log2_bwidth, loss_lsize[1]);

		for (int i = 0; i < shifts_dims[0]; i++) {
			//LOSS GLOBAL WIDTH
			int losses_width = fit_units(1 << log2_bwidth, step_size, image_dims[0] - (i << log2_bwidth));
			if (losses_width > shift_wg_max_size * 2) losses_width = 2 * shift_wg_max_size;
			loss_gsize[0] = make_even_div(losses_width << log2_bwidth, loss_lsize[0]);
			std::cout << "Loss kernel global width (" << i << "th): " << loss_gsize[0] << std::endl;

			//SHIFT LOCAL & GLOBAL SIZE
			shift_lsize[0] = shift_gsize[0] = ceil_div(losses_width, 2);
			shift_lsize[1] = shift_wg_max_size / shift_lsize[0];
			shift_gsize[1] = make_even_div(shifts_dims[1], shift_lsize[1]);

			//LOSS KERNEL LAUNCH
			err[err_idx++] = clSetKernelArg(loss_kernel, 3, sizeof(cl_int), &losses_width);
			err[err_idx++] = clSetKernelArg(loss_kernel, 6, sizeof(cl_int), &i);
			err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, loss_kernel, 2, NULL, loss_gsize, loss_lsize, 0, NULL, &loss_complete);
			if (err[err_idx - 1]) std::cout << getErrorString(err[err_idx - 1]) << std::endl;
			err[err_idx++] = clWaitForEvents(1, &loss_complete);

			//SHIFT KERNEL LAUNCH
			err[err_idx++] = clSetKernelArg(shift_kernel, 3, sizeof(cl_int), &losses_width);
			err[err_idx++] = clSetKernelArg(shift_kernel, 4, sizeof(cl_int), &i);
			err[err_idx++] = clSetKernelArg(shift_kernel, 5, losses_width * shift_lsize[1] * sizeof(cl_float), NULL);
			err[err_idx++] = clSetKernelArg(shift_kernel, 6, losses_width * shift_lsize[1] * sizeof(cl_int), NULL);
			err[err_idx++] = clSetKernelArg(shift_kernel, 7, losses_width * shift_lsize[1] * sizeof(cl_float), NULL);
			err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, shift_kernel, 2, NULL, shift_gsize, shift_lsize, 0, NULL, &shift_complete);
			err[err_idx++] = clWaitForEvents(1, &shift_complete);
		}
		err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, vertex_kernel, 2, NULL, vertex_gsize, vertex_lsize, 0, NULL, &vertex_complete);
		err[err_idx++] = clWaitForEvents(1, &vertex_complete);
	}
	~DP_State() {
		delete[] err;
		err = new cl_int[1];
		err[0] = clReleaseKernel(halver_kernel);
		err[0] = clReleaseKernel(loss_kernel);
		err[0] = clReleaseKernel(shift_kernel);
		err[0] = clReleaseKernel(vertex_kernel);
		err[0] = clReleaseProgram(program);
		err[0] = clReleaseMemObject(image1);
		err[0] = clReleaseMemObject(image2);
		err[0] = clReleaseMemObject(losses);
		err[0] = clReleaseMemObject(shifts);
		err[0] = clReleaseMemObject(certainties);
		err[0] = clReleaseMemObject(vertices);
		err[0] = clReleaseMemObject(vcolors);
		delete[] err;
	}
};

int* GetShifts(const DP_State* dps) {
	int* shifts_cpy = new int[dps->shifts_dims[0] * dps->shifts_dims[1]];
	cl_event copy_complete;
	cl_int err = clEnqueueReadBuffer(dps->gen_state->queue, dps->shifts, false, 
		0, dps->shifts_buf_size, shifts_cpy, 0, NULL, &copy_complete);
	err = clWaitForEvents(1, &copy_complete);
	return shifts_cpy;
}

float* GetCertainties(const DP_State* dps) {
	float* certainties_cpy = new float[dps->shifts_dims[0] * dps->shifts_dims[1]];
	cl_event copy_complete;
	cl_int err = clEnqueueReadBuffer(dps->gen_state->queue, dps->certainties, false,
		0, dps->certainties_buf_size, certainties_cpy, 0, NULL, &copy_complete);
	err = clWaitForEvents(1, &copy_complete);
	return certainties_cpy;
}

Vector3* GetVertices(const DP_State* dps) {
	const int num_vertices = dps->shifts_dims[0] * dps->shifts_dims[1];
	float* vertices_cpy = new float[num_vertices * 4];
	cl_event copy_complete;
	cl_int err = clEnqueueReadBuffer(dps->gen_state->queue, dps->vertices, false,
		0, dps->vertices_buf_size, vertices_cpy, 0, NULL, &copy_complete);
	err = clWaitForEvents(1, &copy_complete);
	Vector3* vertices_cpy_cpy = new Vector3[num_vertices];
	for (int i = 0; i < num_vertices; i++) 
		vertices_cpy_cpy[i] = Vector3(vertices_cpy[4 * i], vertices_cpy[4 * i + 1], vertices_cpy[4 * i + 2]);
	delete[] vertices_cpy;
	return vertices_cpy_cpy;
}

unsigned int* GetVColors(const DP_State* dps) {
	unsigned int* vcolors_cpy = new unsigned int[(dps->shifts_dims[0] + dps->cushion) * dps->shifts_dims[1]];
	cl_event copy_complete;
	cl_int err = clEnqueueReadBuffer(dps->gen_state->queue, dps->vcolors, false,
		0, dps->vcolors_buf_size, vcolors_cpy, 0, NULL, &copy_complete);
	err = clWaitForEvents(1, &copy_complete);
	return vcolors_cpy;
}

unsigned int* GetInternalImg(const DP_State* dps, bool idx) {
	unsigned int* img_cpy = new unsigned int[dps->image_buf_size];
	cl_event copy_complete;
	const cl_mem* img_bufs[2] = { &(dps->image1), &(dps->image2) };
	cl_int err = clEnqueueReadBuffer(dps->gen_state->queue, *(img_bufs[idx]), false,
		0, dps->image_buf_size, img_cpy, 0, NULL, &copy_complete);
	err = clWaitForEvents(1, &copy_complete);
	return img_cpy;
}

unsigned int* GetReconstructedImg(const DP_State* dps) {
	unsigned int* image2 = GetInternalImg(dps, 1);
	int* shifts = GetShifts(dps);
	int batch_width = 1 << dps->log2_bwidth;
	int step_size = dps->step_size;

	int image1_size = dps->image_dims[0] * dps->image_dims[1];
	unsigned int* image1 = new unsigned int[image1_size];
	for (int i = 0; i < image1_size; i++) image1[i] = 0;

	int shift, px_idx;
	for (int bx = 0; bx < dps->shifts_dims[0]; bx++) {
		for (int by = 0; by < dps->shifts_dims[1]; by++) {
			shift = shifts[bx + by * dps->shifts_dims[0]] * step_size;
			for (int x = 0; x < batch_width; x++) {
				if (x + shift + bx * batch_width >= dps->image_dims[0]) continue;
				for (int y = 0; y < batch_width; y++) {
					px_idx = x + bx * batch_width + (y + by * batch_width) * dps->image_dims[0];
					image1[px_idx + shift] = image2[px_idx];
				}
			}
		}
	};

	delete[] image2;
	delete[] shifts;
	return image1;
}

unsigned int* GetDepthVis(const DP_State* dps, float shift_scalar, float certainty_scalar) {
	int* shifts = GetShifts(dps);
	float* certainties = GetCertainties(dps);
	int batch_width = 1 << dps->log2_bwidth;
	int step_size = dps->step_size;

	int visualization_size = dps->image_dims[0] * dps->image_dims[1];
	unsigned int* visualization = new unsigned int[visualization_size];
	for (int i = 0; i < visualization_size; i++) visualization[i] = 0;

	float color_disp_scalar, color_scalar;
	float color0[3] = { 1.0F, 0.0F, 0.0F };
	float color1[3] = { 1.0F, 1.0F, 1.0F };
	float coldisp[3] = {
		color0[0] - color1[0],
		color0[1] - color1[1],
		color0[2] - color1[2]
	};
	int r, g, b;
	for (int bx = 0; bx < dps->shifts_dims[0]; bx++) {
		for (int by = 0; by < dps->shifts_dims[1]; by++) {
			color_disp_scalar = shift_scalar * shifts[bx + by * dps->shifts_dims[0]] * step_size;
			if (color_disp_scalar > 1.0F) color_disp_scalar = 1.0F;
			if (color_disp_scalar < 0.0F) color_disp_scalar = 0.0F;
			color_scalar = certainty_scalar * certainties[bx + by * dps->shifts_dims[0]];
			if (color_scalar > 1.0F) color_scalar = 1.0F;
			for (int x = 0; x < batch_width; x++) {
				for (int y = 0; y < batch_width; y++) {
					r = (color1[0] + color_disp_scalar * coldisp[0]) * 255.0F * color_scalar;
					if (r > 255) r = 255;
					g = (color1[1] + color_disp_scalar * coldisp[1]) * 255.0F * color_scalar;
					if (g > 255) g = 255;
					b = (color1[2] + color_disp_scalar * coldisp[2]) * 255.0F * color_scalar;
					if (b > 255) b = 255;
					visualization[x + bx * batch_width + (y + by * batch_width) * dps->image_dims[0]] = (b << 16) + (g << 8) + r;
				}
			}
		}
	};

	delete[] shifts;
	delete[] certainties;
	return visualization;
}

struct Frustum {
	Vector3 p[5];
	Frustum() = default;
	Frustum(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4, const Vector3& p5) {
		p[0] = p1;
		p[1] = p2;
		p[2] = p3;
		p[3] = p4;
		p[4] = p5;
	}
};

Frustum operator *(const Transform3& T, const Frustum& F) {
	return Frustum(T * F.p[0], T * F.p[1], T * F.p[2], T * F.p[3], T * F.p[4]);
}

string Stringify(const Vector3& v) {
	return "{ " + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + " }";
}

string Stringify(const IntVector3& v) {
	return "{ " + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + " }";
}

string Stringify(const Frustum& frustum) {
	return "{" +
		Stringify(frustum.p[0]) + ", \n " +
		Stringify(frustum.p[1]) + ", \n " +
		Stringify(frustum.p[2]) + ", \n " +
		Stringify(frustum.p[3]) + ", \n " +
		Stringify(frustum.p[4]) + "}";
}

Frustum GetFrustum(const DP_State* dps) {
	Frustum frustum;
	float far_plane = 2.0F * dps->vertex_constants[0];
	float xy_dim = dps->cam_state->tangent * far_plane;
	frustum.p[0] = Vector3(xy_dim, xy_dim, far_plane);
	frustum.p[1] = Vector3(-xy_dim, xy_dim, far_plane);
	frustum.p[2] = Vector3(-xy_dim, -xy_dim, far_plane);
	frustum.p[3] = Vector3(xy_dim, -xy_dim, far_plane);
	frustum.p[4] = Vector3(0, 0, 0);
	return frustum;
}

std::string Stringify(const DP_State* dps) {
	std::string str = std::string("Depth Perception State: ") + '\n';
	str += "Halver wg dims: " + std::to_string(dps->halver_lsize[0]) + ", " + std::to_string(dps->halver_lsize[1]) + '\n';
	str += "Loss wg dims: " + std::to_string(dps->loss_lsize[0]) + ", " + std::to_string(dps->loss_lsize[1]) + '\n';
	str += "Halver max wg size: " + std::to_string(dps->halver_wg_max_size) + '\n';
	str += "Loss max wg size: " + std::to_string(dps->loss_wg_max_size) + '\n';
	str += "Shift max wg size: " + std::to_string(dps->shift_wg_max_size) + '\n';
	str += "New image dims: " + std::to_string(dps->image_dims[0]) + ", " + std::to_string(dps->image_dims[1]) + '\n';
	str += "Losses buffer dims: " + std::to_string(dps->losses_dims[0]) + ", " + std::to_string(dps->losses_dims[1]) + '\n';
	str += "Shifts buffer dims: " + std::to_string(dps->shifts_dims[0]) + ", " + std::to_string(dps->shifts_dims[1]) + '\n';
	str += "Image buffer size: " + std::to_string(dps->image_buf_size) + '\n';
	str += "Losses buffer size: " + std::to_string(dps->losses_buf_size) + '\n';
	str += "Shifts buffer size: " + std::to_string(dps->shifts_buf_size) + '\n';
	str += "Certainties buffer size: " + std::to_string(dps->certainties_buf_size) + '\n';
	str += "log_2(width of batch): " + std::to_string(dps->log2_bwidth) + '\n';
	str += "Step size: " + std::to_string(dps->step_size) + '\n';
	str += "Image halves: " + std::to_string(dps->image_halves) + '\n';
	return str;
}

std::string GetErrors(const Gen_State& gen) {
	std::string err_string = "Depth Perception Errors: \n";
	for (int i = 0; i < gen.err_len; i++) if (gen.err[i]) err_string += std::string(getErrorString(gen.err[i])) + "\n";
	return err_string;
}

std::string GetErrors(const DP_State* dps) {
	std::string err_string = "Depth Perception Errors: \n";
	if (dps->build_err) err_string += "Program build log: \n" + std::string(dps->build_log) + "\n";
	for (int i = 0; i < dps->err_len; i++)
		if (dps->err[i]) err_string += std::to_string(i) + ": " + std::string(getErrorString(dps->err[i])) + "\n";
	return err_string;
}

struct Face {
	IntVector3* subp = 0;
	int size = 0;
	IntAlignedBox3 box;
	Vector3 normal;
	Vector3 point;
	float certainty;
	Face() {
		subp = 0;
		size = 0;
		box = IntAlignedBox3();
		normal = Vector3();
		point = Vector3();
		certainty = 0;
	}
	Face(const Face& f) {
		size = f.size;
		subp = new IntVector3[size];
		for (int i = 0; i < size; i++) subp[i] = f.subp[i];
		box = f.box;
		normal = f.normal;
		point = f.point;
		certainty = f.certainty;
	}
	Face(const IntVector3& idx, const Vector3& n, const Vector3& p, float c) {
		subp = new IntVector3[1];
		subp[0] = idx;
		size = 1;
		box = IntAlignedBox3(idx, idx + IntVector3(1, 1, 1));
		normal = n;
		point = p;
		certainty = c;
	}
	Face& operator =(const Face& f) {
		delete[] subp;
		size = f.size;
		subp = new IntVector3[size];
		for (int i = 0; i < size; i++) subp[i] = f.subp[i];
		box = f.box;
		normal = f.normal;
		point = f.point;
		certainty = f.certainty;
		return *this;
	}
	~Face() {
		delete[] subp;
	}
};

std::string Stringify(const Face& f, std::string indent) {
	std::string str = indent + "Face: \n";
	str += indent + "  Subplane indices: { ";
	for (int i = 0; i < f.size; i++) str += Stringify(f.subp[i]) + ", ";
	str += "} (" + std::to_string(f.size) + ")\n";
	str += indent + "  Box: { " + Stringify(f.box.minV) + ", " + Stringify(f.box.maxV) + " }\n";
	str += indent + "  Average normal: { " + 
		std::to_string(f.normal.x) + ", " + std::to_string(f.normal.y) + ", " + std::to_string(f.normal.z) + 
		" }\n";
	str += indent + "  Average point: " + Stringify(f.point) + '\n';
	str += indent + "  Average certainty: " + std::to_string(f.certainty);
	return str;
}

Face Combine(const Face& f1, const Face& f2, const IntAlignedBox3& box, const Vector3& normal, const Vector3& point, float certainty) {
	Face f;
	f.size = f1.size + f2.size;
	f.subp = new IntVector3[f.size];
	int fidx = 0;
	for (int i = 0; i < f1.size; i++) f.subp[fidx++] = f1.subp[i];
	for (int i = 0; i < f2.size; i++) f.subp[fidx++] = f2.subp[i];
	f.box = box;
	f.normal = normal;
	f.point = point;
	f.certainty = certainty;
	return f;
}

struct FaceContainer {
	Face* ifaces = 0;
	Face* ofaces = 0;
	int isize; // can be 0
	int osize; // can be 0
	IntAlignedBox3 box;
	~FaceContainer() {
		delete[] ifaces;
		delete[] ofaces;
	}
};

std::string Stringify(const FaceContainer& fc) {
	std::string str = "Face Container:\n";
	str += "  Inner faces : {\n";
	for (int i = 0; i < fc.isize; i++) str += Stringify(fc.ifaces[i], "  ") + '\n';
	str += " } (" + std::to_string(fc.isize) + ")\n";
	str += "  Outer faces : {\n";
	for (int i = 0; i < fc.osize; i++) str += Stringify(fc.ofaces[i], "  ") + '\n';
	str += " } (" + std::to_string(fc.osize) + ")";
	str += "  Box: {" + Stringify(fc.box.minV) + ", " + Stringify(fc.box.maxV) + " }";
	return str;
}

bool FaceIsOuter(const IntAlignedBox3& facebox, const IntAlignedBox3& fcbox) {
	return facebox.minV[0] == fcbox.minV[0] || facebox.maxV[0] == fcbox.maxV[0] ||
		   facebox.minV[1] == fcbox.minV[1] || facebox.maxV[1] == fcbox.maxV[1] ||
		   facebox.minV[2] == fcbox.minV[2] || facebox.maxV[2] == fcbox.maxV[2];
}

FaceContainer* MergeFaces(const FaceContainer* fc1, const FaceContainer* fc2, 
	int linger_threshold, float* normal_thresholds, float* point_thresholds, int num_slevels) {
	// NEW FACE CONTAINER INITIALIZATION
	int merge_dim = 0 * (fc1->box.maxV.x == fc2->box.minV.x) + 
					1 * (fc1->box.maxV.y == fc2->box.minV.y) + 
					2 * (fc1->box.maxV.z == fc2->box.minV.z);
	FaceContainer* fc = new FaceContainer();
	fc->box = Combine(fc1->box, fc2->box);
	fc->ifaces = new Face[fc1->isize + fc1->osize + fc2->isize + fc2->osize + 1];
	fc->ofaces = new Face[fc1->osize + fc2->osize + 1];
	fc->isize = 0; fc->osize = 0;
	for (int i = 0; i < fc1->isize; i++) fc->ifaces[fc->isize++] = fc1->ifaces[i];
	for (int i = 0; i < fc2->isize; i++) fc->ifaces[fc->isize++] = fc2->ifaces[i];

	// FC2 MANAGEMENT ARRAYS
	int* fc2_outer_contacts = new int[fc2->osize];
	bool* fc2_outer_inclusions = new bool[fc2->osize];
	int fc2_outer_contacts_len = 0;
	for (int i = 0; i < fc2->osize; i++) {
		if (fc2->ofaces[i].box.minV[merge_dim] == fc2->box.minV[merge_dim]) fc2_outer_contacts[fc2_outer_contacts_len++] = i;
		fc2_outer_inclusions[i] = false;
	}

	// INCLUDING FC1-FC2 MERGED FACES AND UNMERGED FC1 FACES
	float nthreshold, pthreshold; int slevel;
	bool fc1_outer_inclusion;
	Vector3 fc1f_norm, fc2f_norm, avgnormal, avgpoint; 
	float avgptw, fc1f_weighed_cert, fc2f_weighed_cert, inv_total_cert; 
	IntAlignedBox3 avgbox;
	Face* fc1f;
	Face* fc2f;
	for (int i = 0; i < fc1->osize; i++) {
		fc1f = &(fc1->ofaces[i]);
		if (fc1f->box.maxV[merge_dim] == fc1->box.maxV[merge_dim]) {
			fc1_outer_inclusion = false;
			for (int j = 0; j < fc2_outer_contacts_len; j++) {
				if (fc2_outer_inclusions[fc2_outer_contacts[j]]) continue;

				fc1f_weighed_cert = fc1f->certainty * fc1f->size;
				fc1f_norm = fc1f->normal;
				fc2f = &(fc2->ofaces[fc2_outer_contacts[j]]);
				fc2f_weighed_cert = fc2f->certainty * fc2f->size;
				fc2f_norm = fc2f->normal;

				slevel = fc1f->size + fc2f->size;
				if (slevel < num_slevels) {
					nthreshold = normal_thresholds[slevel];
					pthreshold = point_thresholds[slevel];
				}
				else {
					nthreshold = normal_thresholds[num_slevels - 1];
					pthreshold = point_thresholds[num_slevels - 1];
				}

				inv_total_cert = 1.0F / (fc1f_weighed_cert + fc2f_weighed_cert);
				avgnormal = Normalize(fc1f_norm * fc1f_weighed_cert + fc2f_norm * fc2f_weighed_cert);
				avgpoint = (fc1f->point * fc1f_weighed_cert + fc2f->point * fc2f_weighed_cert) * inv_total_cert;
				avgptw = -Dot(avgnormal, avgpoint);
				if (square(Dot(avgnormal, fc1f_norm)) > nthreshold &&
					square(Dot(avgnormal, fc2f_norm)) > nthreshold &&
					square(Dot(avgnormal, fc1f->point) + avgptw) < pthreshold &&
					square(Dot(avgnormal, fc2f->point) + avgptw) < pthreshold) {
					avgbox = Combine(fc1f->box, fc2f->box);
					if (FaceIsOuter(avgbox, fc->box))
						fc->ofaces[fc->osize++] = Combine(*fc1f, *fc2f, avgbox, avgnormal, avgpoint,
							(fc1f->certainty * fc1f_weighed_cert + fc2f->certainty * fc2f_weighed_cert) * inv_total_cert);
					else if ((fc1f->size + fc2f->size) >= linger_threshold)
						fc->ifaces[fc->isize++] = Combine(*fc1f, *fc2f, avgbox, avgnormal, avgpoint,
							(fc1f->certainty * fc1f_weighed_cert + fc2f->certainty * fc2f_weighed_cert) * inv_total_cert);
					fc1_outer_inclusion = true;
					fc2_outer_inclusions[fc2_outer_contacts[j]] = true;
				}
			}
			if (!fc1_outer_inclusion) {
				if (FaceIsOuter(fc1f->box, fc->box)) fc->ofaces[fc->osize++] = *fc1f;
				else if (fc1f->size >= linger_threshold) fc->ifaces[fc->isize++] = *fc1f;
			}
		}
		else fc->ofaces[fc->osize++] = *fc1f;
	}

	// INCLUDING UNMATCHED FC2 FACES
	for (int i = 0; i < fc2->osize; i++) if (!fc2_outer_inclusions[i]) {
		if (FaceIsOuter(fc2->ofaces[i].box, fc->box)) fc->ofaces[fc->osize++] = fc2->ofaces[i];
		else if (fc2->ofaces[i].size >= linger_threshold) fc->ifaces[fc->isize++] = fc2->ofaces[i];
	}

	delete[] fc2_outer_contacts;
	delete[] fc2_outer_inclusions;
	return fc;
}

FaceContainer* ExtractImageFaces(const Plane3* planes, const Vector3* points, const float* certainties, 
	int pwidth, int pheight, int xdisp, int ydisp, int log2_width) {
	int halfwidth = 1 << (log2_width - 1);
	if (halfwidth == 0) {
		FaceContainer* single_fc = new FaceContainer;
		single_fc->isize = 0; single_fc->osize = 1;
		single_fc->ifaces = 0; single_fc->ofaces = new Face[1]{ 
			Face(IntVector3(xdisp, ydisp, 0), 
			planes[xdisp + pwidth * ydisp].GetNormal(), 
			points[xdisp + pwidth * ydisp], 
			certainties[xdisp + pwidth * ydisp])
		};
		return single_fc;
	}
	bool right = xdisp + halfwidth < pwidth;
	bool bottom = ydisp + halfwidth < pheight;
	return 0;
}

struct IFE_State {
	const Gen_State* gen_state;
	const DP_State* dp_state;

	cl_int* err;
	int err_len;
	cl_int build_err;
	char build_log[1024] = "";
	cl_program program;
	cl_kernel plane_kernel;
	size_t plane_lsize[2], plane_gsize[2];
	size_t plane_wg_max_size;
	cl_mem planes_buf, points_buf, colors_buf, certainties_buf;
	Plane3* planes; Vector3* points; Vector3* colors; float* certainties;
	int planes_dims[2];
	int planes_buf_size, points_buf_size, colors_buf_size, certainties_buf_size, num_planes;

	IFE_State(const Gen_State* general_state, const DP_State* depth_perception_state) {
		err = new cl_int[err_len = 13];
		int err_idx = 0;
		gen_state = general_state;
		dp_state = depth_perception_state;

		//PROGRAM & KERNEL
		string programString = ife_plane_kernel_string;
		const char* programChars = programString.c_str();
		size_t programStringSize = programString.size();
		program = clCreateProgramWithSource(gen_state->context, 1, &programChars, &programStringSize, &(err[err_idx++]));
		build_err = err[err_idx++] = clBuildProgram(program, 1, &(gen_state->device), NULL, NULL, NULL);
		if (build_err) return;
		plane_kernel = clCreateKernel(program, "get_planes", &(err[err_idx++]));
		err[err_idx++] = clGetKernelWorkGroupInfo(plane_kernel, gen_state->device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &plane_wg_max_size, NULL);

		//BUFFER CONFIGURATION
		planes_dims[0] = dp_state->shifts_dims[0] - 1; planes_dims[1] = dp_state->shifts_dims[1] - 1;
		num_planes = planes_dims[0] * planes_dims[1];
		planes_buf_size = num_planes * sizeof(cl_float4);
		points_buf_size = planes_buf_size;
		colors_buf_size = planes_buf_size;
		certainties_buf_size = num_planes * sizeof(cl_float);
		planes = new Plane3[num_planes];
		points = new Vector3[num_planes];
		colors = new Vector3[num_planes];
		certainties = new float[num_planes];

		//PLANE WORK SIZE
		plane_lsize[1] = 1 << (log2(plane_wg_max_size, false) / 2);
		plane_lsize[0] = plane_wg_max_size / plane_lsize[1];
		plane_gsize[0] = make_even_div(planes_dims[0], plane_lsize[0]);
		plane_gsize[1] = make_even_div(planes_dims[1], plane_lsize[1]);

		//PLANE KERNEL ARGUMENT
		err[err_idx++] = clSetKernelArg(plane_kernel, 0, sizeof(&(dp_state->vertices)), &(dp_state->vertices));
		err[err_idx++] = clSetKernelArg(plane_kernel, 1, sizeof(&(dp_state->vcolors)), &(dp_state->vcolors));
		err[err_idx++] = clSetKernelArg(plane_kernel, 2, sizeof(&(dp_state->certainties)), &(dp_state->certainties));
		err[err_idx++] = clSetKernelArg(plane_kernel, 3, sizeof(cl_int), &(planes_dims[0]));
		err[err_idx++] = clSetKernelArg(plane_kernel, 4, sizeof(cl_int), &(planes_dims[1]));
		err[err_idx++] = clSetKernelArg(plane_kernel, 5, sizeof(cl_int), &(dp_state->vcolors_dims[0]));
		err[err_idx++] = clSetKernelArg(plane_kernel, 6, (plane_lsize[0] + 1) * (plane_lsize[1] + 1) * sizeof(cl_float4), NULL);
		err[err_idx++] = clSetKernelArg(plane_kernel, 7, (plane_lsize[0] + 1) * (plane_lsize[1] + 1) * sizeof(cl_float4), NULL);
		err[err_idx++] = clSetKernelArg(plane_kernel, 8, (plane_lsize[0] + 1) * (plane_lsize[1] + 1) * sizeof(cl_float), NULL);
			//SKIP 9 ("planes")
			//SKIP 10 ("ppoints")
			//SKIP 11 ("pcolors")
			//SKIP 12 ("pcertainties")
	}

	void Execute() {
		delete[] err;
		err = new cl_int[err_len = 22];
		int err_idx = 0;
		cl_event plane_complete, copy_complete;

		// CREATE BUFFERS
		planes_buf = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, planes_buf_size, NULL, &(err[err_idx++]));
		points_buf = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, points_buf_size, NULL, &(err[err_idx++]));
		colors_buf = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, colors_buf_size, NULL, &(err[err_idx++]));
		certainties_buf = clCreateBuffer(gen_state->context, CL_MEM_READ_WRITE, certainties_buf_size, NULL, &(err[err_idx++]));

		// PLANE KERNEL LAUNCH
		err[err_idx++] = clSetKernelArg(plane_kernel, 9, sizeof(&planes_buf), &planes_buf);
		err[err_idx++] = clSetKernelArg(plane_kernel, 10, sizeof(&points_buf), &points_buf);
		err[err_idx++] = clSetKernelArg(plane_kernel, 11, sizeof(&colors_buf), &colors_buf);
		err[err_idx++] = clSetKernelArg(plane_kernel, 12, sizeof(&certainties_buf), &certainties_buf);
		err[err_idx++] = clEnqueueNDRangeKernel(gen_state->queue, plane_kernel, 2, 0, plane_gsize, plane_lsize, 0, NULL, &plane_complete);
		err[err_idx++] = clWaitForEvents(1, &plane_complete);

		// READ PLANES BUFFER
		err[err_idx++] = clEnqueueReadBuffer(gen_state->queue, planes_buf, false, 0, planes_buf_size, planes, 0, NULL, &copy_complete);
		err[err_idx++] = clWaitForEvents(1, &copy_complete);

		// READ POINTS BUFFER
		float* temp = new float[num_planes * 4];
		err[err_idx++] = clEnqueueReadBuffer(gen_state->queue, points_buf, false, 0, points_buf_size, temp, 0, NULL, &copy_complete);
		err[err_idx++] = clWaitForEvents(1, &copy_complete);
		for (int i = 0; i < num_planes; i++) points[i] = Vector3(temp[i * 4], temp[i * 4 + 1], temp[i * 4 + 2]);

		// READ COLORS BUFFER
		err[err_idx++] = clEnqueueReadBuffer(gen_state->queue, colors_buf, false, 0, colors_buf_size, temp, 0, NULL, &copy_complete);
		err[err_idx++] = clWaitForEvents(1, &copy_complete);
		for (int i = 0; i < num_planes; i++) colors[i] = Vector3(temp[i * 4], temp[i * 4 + 1], temp[i * 4 + 2]);
		delete[] temp;

		// READ CERTAINTIES BUFFER
		err[err_idx++] = clEnqueueReadBuffer(gen_state->queue, certainties_buf, false, 0, certainties_buf_size, certainties, 0, NULL, &copy_complete);
		err[err_idx++] = clWaitForEvents(1, &copy_complete);

		// RELEASE BUFFERS
		err[err_idx++] = clReleaseMemObject(planes_buf);
		err[err_idx++] = clReleaseMemObject(points_buf);
		err[err_idx++] = clReleaseMemObject(colors_buf);
		err[err_idx++] = clReleaseMemObject(certainties_buf);
	}

	~IFE_State() {
		delete[] planes;
		delete[] points;
		delete[] colors;
		delete[] certainties;
		delete[] err;
		err = new cl_int[1];
		err[0] = clReleaseKernel(plane_kernel);
		err[0] = clReleaseProgram(program);
		delete[] err;
	}
};

const Plane3* GetPlanes(const IFE_State* ife) {
	return ife->planes; }

const Vector3* GetPoints(const IFE_State* ife) {
	return ife->points; }

const Vector3* GetColors(const IFE_State* ife) {
	return ife->colors; }

const float* GetCertainties(const IFE_State* ife) {
	return ife->certainties; }

std::string Stringify(const IFE_State* ife) {
	std::string str = std::string("Image Face Extraction State: ") + '\n';
	str += "Plane kernel wg dims: " + std::to_string(ife->plane_lsize[0]) + ", " + std::to_string(ife->plane_lsize[1]) + '\n';
	str += "Plane kernel gbl dims: " + std::to_string(ife->plane_gsize[0]) + ", " + std::to_string(ife->plane_gsize[1]) + '\n';
	str += "Plane kernel wg max size: " + std::to_string(ife->plane_wg_max_size) + '\n';
	str += "Number of planes: " + std::to_string(ife->num_planes) + '\n';
	str += "Planes buffer dims: " + std::to_string(ife->planes_dims[0]) + ", " + std::to_string(ife->planes_dims[1]) + '\n';
	str += "Planes buffer size: " + std::to_string(ife->planes_buf_size) + '\n';
	str += "Points buffer size: " + std::to_string(ife->points_buf_size) + '\n';
	str += "Colors buffer size: " + std::to_string(ife->colors_buf_size) + '\n';
	str += "Certainties buffer size: " + std::to_string(ife->certainties_buf_size) + '\n';
	return str;
}

std::string GetErrors(const IFE_State* ife) {
	std::string err_string = "Image Face Extraction Errors: \n";
	if (ife->build_err) err_string += "Program build log: \n" + std::string(ife->build_log) + "\n";
	for (int i = 0; i < ife->err_len; i++)
		if (ife->err[i]) err_string += std::to_string(i) + ": " + std::string(getErrorString(ife->err[i])) + "\n";
	return err_string;
}

struct PlaneCube {
	Vector3 n[2][2][2];
	Vector3 p[2][2][2];
	Vector3 h[3][3][3];
	float c[2][2][2];
	PlaneCube() {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++) {
			n[i][j][k] = Vector3();
			p[i][j][k] = Vector3();
			c[i][j][k] = 0;
		}
		for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
		for (int k = 0; k < 3; k++) 
		h[i][j][k] = Vector3();
	}
	PlaneCube(const PlaneCube& pc) {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++) {
			n[i][j][k] = pc.n[i][j][k];
			p[i][j][k] = pc.p[i][j][k];
			c[i][j][k] = pc.c[i][j][k];
		}
		for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
		for (int k = 0; k < 3; k++) 
		h[i][j][k] = pc.h[i][j][k];
	}
};

string Stringify(const PlaneCube& pc, const std::string& indent) {
	string str = indent + "Planes: { \n";
	for (int i = 0; i < 2; i++) 
	for (int j = 0; j < 2; j++) 
	for (int k = 0; k < 2; k++) 
	if (pc.c[i][j][k])
	str += indent + " [ " + to_string(i) + ", " + to_string(j) + ", " + to_string(k) + " ]: { " + 
		to_string(pc.n[i][j][k].x) + ", " + to_string(pc.n[i][j][k].y) + ", " + to_string(pc.n[i][j][k].z) + 
		" } \n";
	// Add points
	str += indent + "}\n";
	str += indent + "Hues: {\n";
	for (int i = 0; i < 3; i++) 
	for (int j = 0; j < 3; j++) 
	for (int k = 0; k < 3; k++)
	if (Dot(pc.h[i][j][k], pc.h[i][j][k]))
	str += indent + " [ " + to_string(i) + ", " + to_string(j) + ", " + to_string(k) + " ]: { " +
		to_string(pc.h[i][j][k].x) + ", " + 
		to_string(pc.h[i][j][k].y) + ", " + 
		to_string(pc.h[i][j][k].z) + " } \n";
	str += indent + "}\n";
	//str += indent + "Certainty: " + to_string(pc.c);
	return str;
}

bool pcHasContents (const PlaneCube& pc) {
	return pc.c[0][0][0] || pc.c[0][0][1] || pc.c[0][1][0] || pc.c[0][1][1] ||
		   pc.c[1][0][0] || pc.c[1][0][1] || pc.c[1][1][0] || pc.c[1][1][1];
}

struct HueCoefficients {
	float c[2][2][2];
	HueCoefficients(const Vector3& p) {
		Vector3 minusp = Vector3(1, 1, 1) - p;
		c[0][0][0] = minusp.x * minusp.y * minusp.z;
		c[0][0][1] = minusp.x * minusp.y * p.z;
		c[0][1][0] = minusp.x * p.y * minusp.z;
		c[0][1][1] = minusp.x * p.y * p.z;
		c[1][0][0] = p.x * minusp.y * minusp.z;
		c[1][0][1] = p.x * minusp.y * p.z;
		c[1][1][0] = p.x * p.y * minusp.z;
		c[1][1][1] = p.x * p.y * p.z;
	}
};

Vector3 InterpolateHue(const HueCoefficients& hc, const PlaneCube& pc, int x, int y, int z) {
	Vector3 interpolation = Vector3();
	for (int i = 0; i < 2; i++) 
	for (int j = 0; j < 2; j++) 
	for (int k = 0; k < 2; k++) 
	interpolation += pc.h[x + i][y + j][z + k] * hc.c[i][j][k];
	return interpolation;
}

void OptimizeHues(PlaneCube* pc, const HueCoefficients& hc, const Vector3& dL_dh, int x, int y, int z) {
	for (int i = 0; i < 2; i++) 
	for (int j = 0; j < 2; j++) 
	for (int k = 0; k < 2; k++) 
	pc->h[x + i][y + j][z + k] -= dL_dh * hc.c[i][j][k];
}

struct PlaneContainer {
	PlaneContainer* sub[2][2][2];
	PlaneCube* planes;
	PlaneContainer() {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++) 
		sub[i][j][k] = 0;
		planes = 0;
	}
	PlaneContainer(const PlaneContainer& PC) {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++) {
			if (PC.sub[i][j][k]) sub[i][j][k] = new PlaneContainer(*(PC.sub[i][j][k]));
			else sub[i][j][k] = 0;
		}
		if (PC.planes) planes = new PlaneCube(*(PC.planes));
		else planes = 0;
	}
	~PlaneContainer() {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++) 
		delete sub[i][j][k];
		delete planes;
	}
};

bool PCHasContents(const PlaneContainer& PC) {
	return PC.sub[0][0][0] || PC.sub[0][0][1] ||
		   PC.sub[0][1][0] || PC.sub[0][1][1] ||
		   PC.sub[1][0][0] || PC.sub[1][0][1] ||
		   PC.sub[1][1][0] || PC.sub[1][1][1] ||
		   PC.planes;
}

std::string Stringify(const PlaneContainer& PC, const std::string& indent) {
	std::string str = "";
	for (int i = 0; i < 2; i++) 
	for (int j = 0; j < 2; j++) 
	for (int k = 0; k < 2; k++) 
	if (PC.sub[i][j][k])
	if (PCHasContents(*(PC.sub[i][j][k])))
	str += indent + "[ " + to_string(i) + ", " + to_string(j) + ", " + to_string(k) + " ]: {\n" + 
		Stringify(*(PC.sub[i][j][k]), indent + ' ') + "\n" + 
		indent + "}\n";
	if (PC.planes) if (pcHasContents(*(PC.planes))) str += indent + "Plane Cube:\n" + Stringify(*(PC.planes), indent + ' ');
	return str;
}

unsigned int GetSoloSubPC(const PlaneContainer& PC) {
	unsigned int solo_idx = 8;
	int num_sub_pcs = 0;
	for (int i = 0; i < 2; i++)
	for (int j = 0; j < 2; j++)
	for (int k = 0; k < 2; k++)
	if (PC.sub[i][j][k]) { 
		solo_idx = i + j * 2 + k * 4;
		num_sub_pcs++;
	}
	if ((num_sub_pcs > 1) || PC.planes) solo_idx = 8;
	return solo_idx;
}

void Insert(PlaneContainer* PC, int log2_width, Vector3 PC_disp, const Vector3& n, const Vector3& p, const Vector3& hue, float certainty) {
	PlaneContainer* current_PC = PC;
	int i, j, k;
	for (int halfwidth = 1 << (log2_width - 1); halfwidth > 1; halfwidth >>= 1) {
		i = (p.x - PC_disp.x) > halfwidth; j = (p.y - PC_disp.y) > halfwidth; k = (p.z - PC_disp.z) > halfwidth;
		if (current_PC->sub[i][j][k] == 0) current_PC->sub[i][j][k] = new PlaneContainer();
		current_PC = current_PC->sub[i][j][k];
		PC_disp += Vector3(i * halfwidth, j * halfwidth, k * halfwidth);
	}
	i = (p.x - PC_disp.x) > 1; j = (p.y - PC_disp.y) > 1; k = (p.z - PC_disp.z) > 1; PC_disp += Vector3(i, j, k);
	if (current_PC->planes == 0) current_PC->planes = new PlaneCube();
	PlaneCube* pc = current_PC->planes;
	Vector3 avg_normal = Normalize(n * certainty * NEWPLANECERTAINTYSCALE + pc->n[i][j][k] * pc->c[i][j][k]);
	Vector3 avg_point = p * certainty * NEWPLANECERTAINTYSCALE;
	float certainty_sum = certainty * NEWPLANECERTAINTYSCALE;
	if (pc->c[i][j][k]) {
		avg_point += pc->p[i][j][k] * pc->c[i][j][k];
		certainty_sum += pc->c[i][j][k];
	}
	avg_point /= certainty_sum;
	pc->n[i][j][k] = avg_normal;
	pc->p[i][j][k] = avg_point;
	HueCoefficients hc(p - PC_disp);
	OptimizeHues(pc, hc, (InterpolateHue(hc, *pc, i, j, k) - hue) * COLOROPTIMIZATIONRATE, i, j, k);
	pc->c[i][j][k] = (pc->c[i][j][k] * pc->c[i][j][k] + certainty * certainty) / (pc->c[i][j][k] + certainty);
}

void Blast(PlaneContainer* PC, int log2_width, const Vector3& PC_disp, const Vector3& z_dir, float z, float certainty, const Line3& ray) {
	int log2_halfwidth = log2_width - 1;
	int halfwidth = 1 << log2_halfwidth;
	Vector3 PC_center = Vector3(halfwidth, halfwidth, halfwidth) + PC_disp;
	Vector3 nearest_disp = PC_center - Project(PC_center, ray); // or - (ray.t + Project(PC_center - ray.t, ray.d)) for unnormalized ray
	float center_z = Dot(z_dir, PC_center - ray.t);
	if (Dot(nearest_disp, nearest_disp) > (3 * halfwidth * halfwidth)) return;
	if ((center_z + SQRT3 * halfwidth) <= 0) return;
	if ((center_z - SQRT3 * halfwidth) > z) return;
	if (log2_halfwidth) {
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++)
			if (PC->sub[i][j][k]) {
				Blast(PC->sub[i][j][k], log2_halfwidth, PC_disp + Vector3(halfwidth * i, halfwidth * j, halfwidth * k), 
					z_dir, z, certainty, ray);
				if (!PCHasContents(*(PC->sub[i][j][k]))) { delete PC->sub[i][j][k]; PC->sub[i][j][k] = 0; }
			}
	}
	else {
		Vector3 collision;
		float collision_z;
		if (PC->planes == 0) return;
		for (int i = 0; i < 2; i++) 
		for (int j = 0; j < 2; j++) 
		for (int k = 0; k < 2; k++)
		if (PC->planes->c[i][j][k]) {
			collision = Meet(ray, Plane3(PC->planes->n[i][j][k], PC->planes->p[i][j][k]));
			collision_z = Dot(collision - ray.t, z_dir);
			collision -= PC_disp + Vector3(i, j, k);
			if ((collision_z > 0) && (collision_z < z) && ( 
				(collision.x < i) || (collision.x > i + 1) ||
				(collision.y < j) || (collision.y > j + 1) ||
				(collision.z < k) || (collision.z > k + 1)))
				PC->planes->c[i][j][k] = max(PC->planes->c[i][j][k] - PASSTHROUGHCOEFFICIENT * certainty, 0);
		}
		if (!pcHasContents(*(PC->planes))) { delete PC->planes; PC->planes = 0; }
	}
}

void Prune(PlaneContainer* PC, int log2_width, float subtractor) {
    if (log2_width > 1) {
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
        if (PC->sub[i][j][k]) {
			Prune(PC->sub[i][j][k], log2_width - 1, subtractor);
			if (!PCHasContents(*(PC->sub[i][j][k]))) {
				delete PC->sub[i][j][k];
				PC->sub[i][j][k] = 0;
			}
		}
    }
    else if (PC->planes) {
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
		PC->planes->c[i][j][k] = max(PC->planes->c[i][j][k] - subtractor, 0);
		if (!pcHasContents(*(PC->planes))) { delete PC->planes; PC->planes = 0; }
    }
}

struct PMU_State {
	const IFE_State* ife_state;

	PlaneContainer* plane_map;
	int log2_pmwidth;
	Vector3 pm_disp;
	Transform3 pose;
	Frustum frustum;

	PMU_State(const IFE_State* image_face_extraction_state) {
		ife_state = image_face_extraction_state;
		plane_map = new PlaneContainer();
		log2_pmwidth = 1;
		pm_disp = -Vector3(1, 1, 1);
		pose = Transform3();
		frustum = GetFrustum(ife_state->dp_state);
	}

	void SetPose(const Transform3& T) {
		pose = T;
	}

	void Fit() {
		int width;
		Vector3 pm_point;
		for (int f = 0; f < 5; f++) {
			width = 1 << log2_pmwidth;
			pm_point = pose * frustum.p[f] - pm_disp;
			while ((pm_point.x < 0) || (pm_point.x > width) ||
				   (pm_point.y < 0) || (pm_point.y > width) ||
				   (pm_point.z < 0) || (pm_point.z > width)) {
				PlaneContainer* new_plane_map = new PlaneContainer();
				int i = pm_point.x < 0; int j = pm_point.y < 0; int k = pm_point.z < 0;
				new_plane_map->sub[i][j][k] = plane_map;
				plane_map = new_plane_map;
				pm_disp -= Vector3(width * i, width * j, width * k);

				width = 1 << (++log2_pmwidth);
				pm_point = pose * frustum.p[f] - pm_disp;
			}
		}
	}

	void Shrink() {
		unsigned int sub_pc_idx;
		PlaneContainer* sub_pc;
		int i, j, k, halfwidth;
		while ((sub_pc_idx = GetSoloSubPC(*plane_map)) < 8) {
			i = sub_pc_idx & 1; j = (sub_pc_idx & 2) >> 1; k = (sub_pc_idx & 4) >> 2;
			sub_pc = plane_map->sub[i][j][k];
			plane_map->sub[i][j][k] = 0;
			delete plane_map;
			plane_map = sub_pc;

			halfwidth = 1 << (--log2_pmwidth);
			pm_disp += Vector3(halfwidth * i, halfwidth * j, halfwidth * k);
		}
	}

	void Execute() {
		const Plane3* planes = GetPlanes(ife_state);
		const Vector3* points = GetPoints(ife_state);
		const Vector3* colors = GetColors(ife_state);
		const float* certainties = GetCertainties(ife_state);
		Fit();
		int idx;
		for (int y = 0; y < ife_state->planes_dims[1]; y++) {
			for (int x = 0; x < ife_state->planes_dims[0]; x++) {
				idx = x + y * ife_state->planes_dims[0];
				Blast(plane_map, log2_pmwidth, pm_disp, pose.R(2), points[idx].z, certainties[idx],
					Line3(Normalize(pose.R * points[idx]), pose.t));
				Insert(plane_map, log2_pmwidth, pm_disp, pose * planes[idx].GetNormal(), pose * points[idx], colors[idx], certainties[idx]);
			}
		}
		Shrink();
		Prune(plane_map, log2_pmwidth, PRUNERATE);
	}

	~PMU_State() {
		delete plane_map;
	}
};

AlignedBox3 GetPlaneMapBox(const PMU_State* pmu) {
	float width = 1 << pmu->log2_pmwidth;
	return AlignedBox3(pmu->pm_disp, Vector3(width, width, width) + pmu->pm_disp);
}

