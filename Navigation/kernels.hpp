#pragma once

#include <string>

std::string img_halver_kernel_string =
"kernel void halve_img(global const uchar4* img, global uchar4* hlf, const int width, const int height) {" // 0, 1, 2, 3
"	const int x = get_global_id(0);"
"	const int y = get_global_id(1);"
"	if ((x >= (width >> 1)) || (y >= (height >> 1))) return;"
"	const int top_idx = x * 2 + y * 2 * width;"
"	const int bot_idx = top_idx + width;"
"	uchar8 top_sum = (uchar8)(img[top_idx], img[top_idx + 1]);"
"	uchar8 bot_sum = (uchar8)(img[bot_idx], img[bot_idx + 1]);"
"	hlf[x + y * (width >> 1)] = convert_uchar4((convert_float4(top_sum.hi) + convert_float4(top_sum.lo) + "
"								convert_float4(bot_sum.hi) + convert_float4(bot_sum.lo)) / 4.0F);"
"}";

std::string dp_loss_kernel_string =
"kernel void get_losses(global const uchar4* image1, global const uchar4* image2, " // 0, 1
"const int width, const int losses_width, const int losses_height, const int batch_size, " // 2, 3, 4, 5
"const int image2_col, const int step_size, local float* pixel_loss, global float* batch_loss) {" // 6, 7, 8, 9

"	const int gid0 = get_global_id(0);"
"	const int gid1 = get_global_id(1);"
"	const int gid0b = gid0 >> batch_size;"
"	const int gid1b = gid1 >> batch_size;"

"	if (gid0b >= losses_width || gid1b >= losses_height) return;"

"	const int lid0 = get_local_id(0);"
"	const int lid1 = get_local_id(1);"
"	const int lid0b = lid0 >> batch_size;"
"	const int lid1b = lid1 >> batch_size;"

"	const int lsize0 = get_local_size(0);"
"	const int lsize1 = get_local_size(1);"
"	const int lsize0b = lsize0 >> batch_size;"
"	const int lsize1b = lsize1 >> batch_size;"
"	const int barea = 1 << (2 * batch_size);"

"	const int sbid0 = gid0 - (gid0b << batch_size);"
"	const int sbid1 = gid1 - (gid1b << batch_size);"
"	const int sbidx = sbid0 + (sbid1 << batch_size);"

"	const int pixel2_idx = gid1 * width + (image2_col << batch_size) + sbid0;"
"	const int pixel1_idx = pixel2_idx + gid0b * step_size;"

"	const uchar4 pixel1_color = image1[pixel1_idx];"
"	const uchar4 pixel2_color = image2[pixel2_idx];"
"	float4 color_disp = convert_float4(pixel2_color) - convert_float4(pixel1_color);"
"	color_disp *= color_disp;"
//"	const int pixel_disp = pixel2_idx - pixel1_idx;"

//"	const float batch_disp_weight = 0;" //(0x1.p-32); // should be equal to weight / squared width

"	const int lidx = (lid1b * lsize0b + lid0b) * barea + sbidx;"
"	pixel_loss[lidx] = (0x1.p-16) * (color_disp.x + color_disp.y + color_disp.z) + 0x1.p-16;" 
// + batch_disp_weight * pixel_disp * pixel_disp

"	for (int array_width = barea >> 1; sbidx < array_width; array_width >>= 1) {"
"		barrier(CLK_LOCAL_MEM_FENCE);"
"		pixel_loss[lidx] += pixel_loss[lidx + array_width];"
"	}"

"	if (sbidx == 0) batch_loss[gid0b + gid1b * losses_width] = pixel_loss[lidx];"
"}";

// global id 0 = image1 batch column#
// global id 1 = image batch row#
// global id 2 = batch pixel#
// global size 2 = batch area in square pixels
// 
// local id 0 = batch column# within work group
// local size 0 = work group width in batches
// local id 1 = batch row# within work group
// local size 1 = work group height in batches
// local id 2 = batch pixel# (same as global)
// local size 2 = batch area in square pixels (same as global)
// 
// images, pixel_loss, and batch_loss are all stored in row_major order
// 
// pixel_loss = float[local size 0][local size 1][local size 2 = 1<<(2*batch_size)]

std::string dp_shift_kernel_string =
"kernel void get_shifts(global const float* batch_losses, const int widthb, const int heightb, " // 0, 1, 2
"const int losses_width, const int image2_col, local float* loss_min, " // 3, 4, 5
"local int* loss_min_idx, local float* loss_sum, global int* shift, global float* certainty) {" // 6, 7, 8, 9

"	const int batchx = get_local_id(0);"
"	const int batchy = get_global_id(1);"
"	if (!(batchy < heightb)) return;"
"	const int row_idx = batchy * losses_width;"

"	const int wg_row_idx = get_local_id(1) * losses_width;"
"	const int wg_idx = wg_row_idx + batchx;"

"	int len_in = losses_width;"
"	int len_out = get_local_size(0);"
"	bool extra_inclusion = (batchx + len_out) < len_in;"
"	const int outer_batchx1 = batchx + len_out * extra_inclusion;"
"	int outer_wg_idx = outer_batchx1 + wg_row_idx;"
"	bool comparison = false;"

"	loss_min_idx[wg_idx] = batchx;"
"	loss_min_idx[outer_wg_idx] = outer_batchx1;"
"	loss_sum[wg_idx] = batch_losses[row_idx + batchx];"
"	loss_sum[outer_wg_idx] = batch_losses[row_idx + outer_batchx1];"
"	loss_min[wg_idx] = loss_sum[wg_idx];"
"	loss_min[outer_wg_idx] = loss_sum[outer_wg_idx];"

"	while ((batchx < len_out) && (len_in - 1)) {"
"		barrier(CLK_LOCAL_MEM_FENCE);"
"		loss_sum[wg_idx] += extra_inclusion * loss_sum[outer_wg_idx];"
"		comparison = loss_min[wg_idx] > loss_min[outer_wg_idx];"
"		loss_min_idx[wg_idx] += comparison * (loss_min_idx[outer_wg_idx] - loss_min_idx[wg_idx]);"
"		loss_min[wg_idx] += comparison * (loss_min[outer_wg_idx] - loss_min[wg_idx]);"

"		len_in = len_out;"
"		len_out = (len_in >> 1) + (len_in & 1);" // ceilinged division of len_in by 2
"		extra_inclusion = (batchx + len_out) < len_in;"
"		outer_wg_idx = wg_idx + len_out * extra_inclusion;"
"	}"
"	const gb_idx = image2_col + batchy * widthb;"
"	if (batchx == 0) {"
"		shift[gb_idx] = loss_min_idx[wg_row_idx];"
"		certainty[gb_idx] = (loss_sum[wg_row_idx] / (losses_width * loss_min[wg_row_idx])) * "
"			(1.0 - 0.5 * (loss_min_idx[wg_row_idx] == 0));"
"	}"
"}";

// certainty = avg loss : min loss ratio
// all matrices here are stored in row-major order

std::string dp_vertex_kernel_string =
"kernel void get_vertices(global const int* shifts, const int shifts_width, const int shifts_height, " // 0, 1, 2
"const float zconstant, const float xyconstant, " // 3, 4
"const float xoffset, const float yoffset, global float4* vertices) {" // 5, 6, 7
"	const int gid0 = get_global_id(0);"
"	const int gid1 = get_global_id(1);"
"	if (gid0 >= shifts_width || gid1 >= shifts_height) return;"
"	const int shift_idx = gid0 + gid1 * shifts_width;"
"	float4 vertex = (float4)(0, 0, 0, 0);"
"	const int shift = shifts[shift_idx];"
"	const float xyconstantz = xyconstant * (vertex.z = zconstant/(shift + (shift == 0) * 0x1.p-1));"
"	vertex.x = (gid0 + xoffset) * xyconstantz;"
"	vertex.y = -(gid1 + yoffset) * xyconstantz;"
"	vertices[shift_idx] = vertex;"
"}";

std::string ife_plane_kernel_string =
"kernel void get_planes(global const float4* vertices, global const uchar4* vcolors, global const float* certainties, " // 0, 1, 2
"const int planes_width, const int planes_height, const int vcolors_width, " // 3, 4, 5
"local float4* local_vertices, local float4* local_vcolors, local float* local_certainties, " // 6, 7, 8
"global float4* planes, global float4* ppoints, global float4* pcolors, global float* pcertainties) {" // 9, 10, 11, 12
"	const int gid0 = get_global_id(0);"
"	const int gid1 = get_global_id(1);"
"	if (gid0 >= planes_width || gid1 >= planes_height) return;"
"	const int lid0 = get_local_id(0);"
"	const int lid1 = get_local_id(1);"
"	const int lsize0 = get_local_size(0);"
"	const int lsize1 = get_local_size(1);"
"	const int vertices_width = planes_width + 1;"
"	const int lvertices_width = lsize0 + 1;"
"	const int vcolors_idx = gid0 + gid1 * vcolors_width;"
"	const int vertices_idx = gid0 + gid1 * vertices_width;"
"	const int lvertices_idx = lid0 + lid1 * lvertices_width;"

"	local_vertices[lvertices_idx] = vertices[vertices_idx];"
"	local_certainties[lvertices_idx] = certainties[vertices_idx];"
"	local_vcolors[lvertices_idx] = convert_float4(vcolors[vcolors_idx]);"

"	bool right_edge = (lid0 == lsize0 - 1) || (gid0 == planes_width - 1);"
"	bool bottom_edge = (lid1 == lsize1 - 1) || (gid1 == planes_height - 1);"
"	if (right_edge && bottom_edge) {"
"		local_vertices[lvertices_idx + 1 + lvertices_width] = vertices[vertices_idx + 1 + vertices_width];"
"		local_certainties[lvertices_idx + 1 + lvertices_width] = certainties[vertices_idx + 1 + vertices_width];"
"		local_vcolors[lvertices_idx + 1 + lvertices_width] = convert_float4(vcolors[vcolors_idx + 1 + vcolors_width]);"
"	}"
"	if (right_edge) { "
"		local_vertices[lvertices_idx + 1] = vertices[vertices_idx + 1];"
"		local_certainties[lvertices_idx + 1] = certainties[vertices_idx + 1];"
"		local_vcolors[lvertices_idx + 1] = convert_float4(vcolors[vcolors_idx + 1]);"
"	}"
"	if (bottom_edge) { "
"		local_vertices[lvertices_idx + lvertices_width] = vertices[vertices_idx + vertices_width];"
"		local_certainties[lvertices_idx + lvertices_width] = certainties[vertices_idx + vertices_width];"
"		local_vcolors[lvertices_idx + lvertices_width] = convert_float4(vcolors[vcolors_idx + vcolors_width]);"
"	}"
"	const float4 points[4] = { "
"		local_vertices[lvertices_idx], local_vertices[lvertices_idx + 1], "
"		local_vertices[lvertices_idx + lvertices_width], local_vertices[lvertices_idx + 1 + lvertices_width]"
"	};"
"	float4 plane = normalize(cross(points[0] - points[1], points[3] - points[1]) + "
"		cross(points[3] - points[2], points[0] - points[2]));"
"	const float4 ppoint = (points[0] + points[1] + points[2] + points[3]) / 4.0F;"
"	plane.w = -dot(plane, ppoint);"

"	const int plane_idx = gid0 + gid1 * planes_width;"
"	planes[plane_idx] = plane;"
"	ppoints[plane_idx] = ppoint;"
"	pcolors[plane_idx] = (local_vcolors[lvertices_idx] + local_vcolors[lvertices_idx + 1] + "
"		local_vcolors[lvertices_idx + lvertices_width] + local_vcolors[lvertices_idx + 1 + lvertices_width]) / 1020.0F;"
"	pcertainties[plane_idx] = min(0.25F * (local_certainties[lvertices_idx] + local_certainties[lvertices_idx + 1] + "
"		local_certainties[lvertices_idx + lvertices_width] + local_certainties[lvertices_idx + 1 + lvertices_width]) / 4.0F, "
"		1.0F);"
"}";