__kernel void updateEFields(__global float* Hx, __global float* Hy, 
		__global float* Ez, __global float* Epsilon, __global float* Mu, 
		float dt, float dy, float dx, int width, int height) {
	int x = get_global_id(0);
	int y = get_global_id(1);
	int index = y * width + x;

	Ez[index] += (dt / Epsilon[index]) * ((Hy[index] - Hy[index - 1])
			/ dx - (Hx[index] - Hx[(y - 1) * width + x]) / dy);
}

__kernel void updateHFields(__global float* Hx, __global float* Hy, 
		__global float* Ez, __global float* Epsilon, __global float* Mu, 
		float dt, float dy, float dx, int width, int height) {
	int x = get_global_id(0);
	int y = get_global_id(1);
	int index = y * width + x;

	Hx[index] -= dt / (Mu[index] * dy) * (Ez[index + width] - Ez[index]);
	Hy[index] += dt / (Mu[index] * dx) * (Ez[index + 1] - Ez[index]);
}

__kernel void visualizeTE2(__global float* image, __global float* Hx, 
		__global float* Hy, __global float* Ez, float minField, 
		float maxField, int width) {
	int x = get_global_id(0);
	int y = get_global_id(1);
	int index = y * width + x;

	image[3 * index] = (Ez[index]*Ez[index] - minField) 
			/ (maxField - minField);
	image[3 * index + 1] = (Hx[index]*Hx[index] - minField) 
			/ (maxField - minField);
	image[3 * index + 2] = (Hy[index]*Hy[index] - minField)
			/ (maxField - minField);
}

__kernel void drawMaterialBoundaries(__global float* image, 
		__global float* boundMask, int width) {
	int x = get_global_id(0);
	int y = get_global_id(1);
	int index = y * width + x;

	float maskColor = 0;
	float L = (image[3 * index] + image[3 * index + 1] 
			+ image[3 * index + 2]) / 3;
	if (L < 0.5) maskColor = 1;

	if (boundMask[index] == 1) {
		image[3 * index] = maskColor;
		image[3 * index + 1] = maskColor;
		image[3 * index + 2] = maskColor;
	}
}
