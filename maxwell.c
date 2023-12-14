#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DT 1e-2
#define MAX_FIELD 1e3
#define MIN_FIELD 0
#define SPEED_OF_LIGHT 299792458.0

void initFields(float* Ez, float* Hx, float* Hy, int width, int height) {
	int index;
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			index = y * width + x;
			Ez[index] = 0;
			Hx[index] = 0;
			Hy[index] = 0;
		}
	}
}

float randNormalFloat(void) {
	return (float)rand() / (float)RAND_MAX;
}

void glfw_error_callback(int error, const char* desc) {
	fprintf(stderr, "GLFW error %d: %s\n", error, desc);
}

float gaussianPulse(float t, float t0, float spread) {
	return exp(-pow(t / (t0 * spread), 2));
}

void updateFields(float* Ez, float* Hx, float* Hy, int width, int height, 
		float* time) {
	const float dx = 1e1;
	const float dy = 1e1;
	const float eps = 8.854e-12;
	const float mu = 1.2566e-6;
	const float dt = 0.9e0 / (SPEED_OF_LIGHT * sqrt(1 / (dx * dx) + 1 / (dy * dy)));

	*time += dt;

	float t0 = 1e-9;
	float spread = 1e-2;
	int sourceX = width / 2;
	int sourceY = height / 2;

	float frequency = 1.5e6;
	float omega = 2 * M_PI * frequency;

//	Ez[sourceY * width + sourceX] += gaussianPulse(time, t0, spread);

	float phi = sin(*time * 1e4);

	Ez[500 * width + 460] += sin(omega * *time + (1 * phi));
	Ez[500 * width + 465] += sin(omega * *time + (2 * phi));
	Ez[500 * width + 470] += sin(omega * *time + (3 * phi));
	Ez[500 * width + 475] += sin(omega * *time + (4 * phi));
	Ez[500 * width + 480] += sin(omega * *time + (5 * phi));
	Ez[500 * width + 485] += sin(omega * *time + (6 * phi));
	Ez[500 * width + 490] += sin(omega * *time + (7 * phi));
	Ez[500 * width + 495] += sin(omega * *time + (8 * phi));
	Ez[500 * width + 500] += sin(omega * *time + (9 * phi));
	Ez[500 * width + 505] += sin(omega * *time + (10 * phi));
	Ez[500 * width + 510] += sin(omega * *time + (11 * phi));
	Ez[500 * width + 515] += sin(omega * *time + (12 * phi));
	Ez[500 * width + 520] += sin(omega * *time + (13 * phi));
	Ez[500 * width + 525] += sin(omega * *time + (14 * phi));
	Ez[500 * width + 530] += sin(omega * *time + (15 * phi));
	Ez[500 * width + 535] += sin(omega * *time + (16 * phi));
	Ez[500 * width + 540] += sin(omega * *time + (17 * phi));

//	Ez[sourceY * width + sourceX + 100] += sin(omega * time);

	for (int j = 0; j < height - 1; j++) {
		for (int i = 0; i < width; i++) {
			int index = j * width + i;
			if (i < width - 1) {
				Hx[index] -= dt / (mu * dy) * (Ez[index + width] - Ez[index]);
			}
			if (j < height - 1) {
				Hy[index] += dt / (mu * dx) * (Ez[index + 1] - Ez[index]);
			}
		}
	}
/*
	for (int j = 1; j < height; j++) {
		for (int i = 1; i < width; i++) {
			int index = j * width + i;
			Ez[index] += dt / eps * ((Hy[index] - Hy[index - 1]) / dx - (Hx[index] - Hx[index - width]) / dy);
		}
	}

	for (int j = 0; j < height - 1; j++) {
        for (int i = 0; i < width - 1; i++) {
            Hx[j * width + i] += (dt / mu) * (Ez[j * width + i + 1] - Ez[j * width + i]) / dy;
            Hy[j * width + i] += (dt / mu) * (Ez[(j + 1) * width + i] - Ez[j * width + i]) / dx;
        }
   	}
*/
   	// Update E field
    for (int j = 1; j < height - 1; j++) {
       	for (int i = 1; i < width - 1; i++) {
       	    Ez[j * width + i] += (dt / eps) * ((Hy[j * width + i] - Hy[j * width + i - 1]) / dx - (Hx[j * width + i] - Hx[(j - 1) * width + i]) / dy);
       	}
   	}
    // Apply PEC boundary conditions: E field is zero at all boundaries
    // Top and bottom boundaries
    for (int i = 0; i < width; i++) {
        Ez[i] = 0; // Top boundary
        Ez[(height - 1) * width + i] = 0; // Bottom boundary
    }

    // Left and right boundaries
    for (int j = 0; j < height; j++) {
        Ez[j * width] = 0; // Left boundary
        Ez[j * width + (width - 1)] = 0; // Right boundary
    }
}

void updateImage(float* Ez, float* Hx, float* Hy, 
		int width, int height, float* time) {
	updateFields(Ez, Hx, Hy, width, height, time);	
	
	printf("%f\n", Ez[512 * width + 512]);
	float* visualData = (float*)malloc(3 * width * height * sizeof(float));
	int index;

	float logMax = log10(MAX_FIELD);
	float logMin = log10(MIN_FIELD);

	float ezMin, hxMin, hyMin;
	float ezMax, hxMax, hyMax;
	ezMin = MAX_FIELD;
	hxMin = MAX_FIELD;
	hyMin = MAX_FIELD;
	ezMax = MIN_FIELD;
	hxMax = MIN_FIELD;
	hyMax = MIN_FIELD;

	for (int i = 0; i < width * height; ++i) {
		if (Ez[i] < ezMin) ezMin = Ez[i];
		if (Ez[i] < ezMin) ezMin = Ez[i];
		if (Hx[i] < hxMin) hxMin = Hx[i];
		if (Hx[i] > hxMax) hxMax = Hx[i];
		if (Hy[i] > hyMax) hyMax = Hy[i];
		if (Hy[i] > hyMax) hyMax = Hy[i];
	}

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			index = y * width + x;
			float ezVal = Ez[index];
			float hxVal = Hx[index];
			float hyVal = Hy[index];
		
			
			visualData[3 * index + 2] = (ezVal - ezMin) / (ezMax - ezMin);
			visualData[3 * index + 1] = (ezVal - ezMin) / (ezMax - ezMin);
			visualData[3 * index] = (hxVal - hxMin) / (hxMax - hxMin);
			
			//visualData[3 * index + 2] = (hyVal*hyVal - hyMin*hyMin) / (hyMax*hyMax - hyMin*hyMin); 
			//visualData[3 * index] = (hxVal*hxVal - hxMin*hxMin) / (hxMax*hxMax - hxMin*hxMin);
			//visualData[3 * index + 2] = (hyVal*hyVal - hyMin*hyMin) / (hyMax*hyMax - hyMin*hyMin); 
			
			//visualData[3 * index] = (ezVal*ezVal - MIN_FIELD) / (MAX_FIELD - MIN_FIELD);
			//visualData[3 * index + 1] = (hxVal*hxVal - MIN_FIELD) / (MAX_FIELD - MIN_FIELD);
			//visualData[3 * index + 2] = (hyVal*hyVal - MIN_FIELD) / (MAX_FIELD - MIN_FIELD); 
			
			//visualData[3 * index] = (log10(ezVal) - logMin) / (logMax - logMin);
			//visualData[3 * index + 1] = (log10(hxVal) - logMin) / (logMax - logMin);
			//visualData[3 * index + 2] = (log10(hyVal) - logMin) / (logMax - logMin);
		}
	}
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_FLOAT, 
			visualData);
	free(visualData);
}

int main(void) {
	glfwSetErrorCallback(glfw_error_callback);

	GLFWwindow* window;
	const int width = 1024, height = 1024;
	
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize glfw\n");
		return -1;
	}

	window = glfwCreateWindow(width, height, "Example", NULL, NULL);
	if (!window) {
		fprintf(stderr, "Failed to create glfw window\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	float* Ez = (float*)malloc(width * height * sizeof(float));
	float* Hx = (float*)malloc(width * height * sizeof(float));
	float* Hy = (float*)malloc(width * height * sizeof(float));
	initFields(Ez, Hx, Hy, width, height);

	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	float* time = (float*)malloc(sizeof(float));
	*time = 0.0f;
	while (!glfwWindowShouldClose(window)) {
		updateImage(Ez, Hx, Hy, width, height, time);
		glClear(GL_COLOR_BUFFER_BIT);
		
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture);
		glBegin(GL_QUADS);
			glTexCoord2f(0, 0); glVertex2f(-1, -1);
			glTexCoord2f(1, 0); glVertex2f(1, -1);
			glTexCoord2f(1, 1); glVertex2f(1, 1);
			glTexCoord2f(0, 1); glVertex2f(-1, 1);
		glEnd();
		glDisable(GL_TEXTURE_2D);

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	free(Ez);
	free(Hx);
	free(Hy);
	glfwTerminate();
	return 0;
}
