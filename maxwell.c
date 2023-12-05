#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>

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

void updateFields(float* Ez, float* Hx, float* Hy, int width, int height) {
	const float dt = 1e-9;
	const float dx = 1e-3;
	const float dy = 1e-3;
	const float eps = 8.854e-12;
	const float mu = 1.2566e-6;

 	for (int j = 0; j < height - 1; j++) {
        for (int i = 0; i < width - 1; i++) {
            Hx[j * width + i] -= (dt / mu) * (Ez[j * width + i + 1] - Ez[j * width + i]) / dy;
            Hy[j * width + i] += (dt / mu) * (Ez[(j + 1) * width + i] - Ez[j * width + i]) / dx;
        }
    }

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
		int width, int height) {
	updateFields(Ez, Hx, Hy, width, height);	
	
	float* visualData = (float*)malloc(width * height * sizeof(float));
	int index;
	float minVisualDatum = 0, maxVisualDatum = 0;
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			index = y * width + x;
			visualData[index] = Ez[index]*Ez[index];
			if (visualData[index] < minVisualDatum) 
					minVisualDatum = visualData[index];
			if (visualData[index] > maxVisualDatum)
					maxVisualDatum = visualData[index];
		}
	}
	for (int i = 0; i < width * height; ++i) {
		visualData[i] -= minVisualDatum;
		visualData[i] /= maxVisualDatum;
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width, height, 0, GL_RED, GL_FLOAT, 
			visualData);
	free(visualData);
}

int main(void) {
	glfwSetErrorCallback(glfw_error_callback);

	GLFWwindow* window;
	const int width = 512, height = 512;
	
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

	while (!glfwWindowShouldClose(window)) {
		updateImage(Ez, Hx, Hy, width, height);
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
