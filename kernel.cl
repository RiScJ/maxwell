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
