#include "maxwell.h"

bool sim_running = true;
bool reset_sim = false;
bool cycle_vis = false;
bool draw_material_boundaries = true;
bool report_framerate = false;
bool just_resumed = false;
bool gpu_support = true;
bool trying_gpu = true;

int min(int a, int b) {
	return b ^ ((a ^ b) & -(a < b));
}

int max (int a, int b) {
	return a ^ ((a ^ b) & -(a < b));
}

void key_callback(GLFWwindow* window, int key, int __attribute__((unused)) 
		scancode, int action, int mods) {
	// Handle Ctrl+C to exit the program
	if (key == GLFW_KEY_C && mods == GLFW_MOD_CONTROL && (action == GLFW_PRESS 
			|| action == GLFW_REPEAT)) {
		printf("Caught interrupt - exiting...\n");
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}

	// Toggle simulation running state with spacebar
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		if (sim_running) {
			printf("Pausing simulation.\n");
		} else {
			printf("Resuming simulation.\n");
			just_resumed = true;
		}
		sim_running = !sim_running;
	}

	// Toggle material boundary rendering with 'B' key
	if (key == GLFW_KEY_B && action == GLFW_PRESS) {
		if (draw_material_boundaries) {
			printf("Disabling material boundary rendering.\n");
		} else {
			printf("Enabling material boundary rendering.\n");
		}
		draw_material_boundaries = !draw_material_boundaries;
	}

	// Report framerate using 'F' key
	if (key == GLFW_KEY_F && action == GLFW_PRESS) {
		report_framerate = true;
	}

	// Reset simulation with 'R' key
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		printf("Resetting simulation.\n");
		sim_running = false;
		reset_sim = true;
	}

	// Cycle through visualization functions with 'V' key
	if (key == GLFW_KEY_V && action == GLFW_PRESS) {
		printf("Advancing to next visualization function.\n");
		cycle_vis = true;
	}
}

int PMLayer(int x, int y, int layers, int width, int height) {
	int min_x = x < (width - x) ? x : (width - x - 1);
	int min_y = y < (height - y) ? y : (height - y - 1);
	int dist_to_edge = min_x < min_y ? min_x : min_y;

	if (dist_to_edge < layers) {
		return layers - 1 - dist_to_edge;
	} else {
		return -1;
	}
}

float conductivityPML(Simulation* simulation, int layer) {
	float norm_layer = (float)layer / (simulation->pml_layers - 1);
	return simulation->pml_conductivity * pow(norm_layer, 
			simulation->pml_sigma_polyorder);
}

void initFields(Field* field, Simulation* simulation) {
	int index, layer;
	for (int y = 0; y < simulation->height; ++y) {
		for (int x = 0; x < simulation->width; ++x) {
			index = y * simulation->width + x;
			field->Epsilon[index] = VACUUM_PERMITTIVITY;
			field->Mu[index] = VACUUM_PERMEABILITY;
			field->Ex[index] = 0;
			field->Ey[index] = 0;
			field->Ez[index] = 0;
			field->Hx[index] = 0;
			field->Hy[index] = 0;
			field->Hz[index] = 0;

			if (simulation->boundary_condition == BC_PML
					&& (layer = PMLayer(x, y, simulation->pml_layers, 
					simulation->width, simulation->height)) >= 0) {
				field->Sigma[index] = conductivityPML(simulation, layer);
			} else {
				field->Sigma[index] = 0;
			}
		}
	}
}

void addMaterials(Field* field, Simulation* simulation, Material* materials) {
	// For each user-specified material
	for (int m = 0; m < simulation->materialc; m++) {
		printf("\rApplying material characteristics... (%d/%d)", m, 
			simulation->materialc);
	
	int index;
	float rel_eps, rel_mu;
	
		// We will implement the properties in a region determined by the 
		// geometry
		switch (materials[m].geom) {
			case MG_UNKNOWN:
				break;
			case MG_TRIANGLE:
				// Extract the relative permittivity and permeability for the
				// triangular region
				rel_eps = materials[m].argv[0].value.floatVal;
				rel_mu = materials[m].argv[1].value.floatVal;

				// Extract the bounding vertices of the triangle
				int x1, y1, x2, y2, x3, y3;
				x1 = materials[m].argv[2].value.intVal;	
				y1 = materials[m].argv[3].value.intVal;	
				x2 = materials[m].argv[4].value.intVal;	
				y2 = materials[m].argv[5].value.intVal;	
				x3 = materials[m].argv[6].value.intVal;	
				y3 = materials[m].argv[7].value.intVal;	
				
				float d1, d2, d3;
				bool has_neg, has_pos;
				for (int y = 0; y < simulation->height; y++) {
					for (int x = 0; x < simulation->width; x++) {
						index = y * simulation->width + x;
						
						// Calculate barycentric coordinates for the point
						// (x, y) in the simulation grid
						d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2);
						d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3);
						d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1);

						has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
						has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
						
						// Check if (x, y) is inside the triangular region
						if (!(has_neg && has_pos)) {
							field->Epsilon[index] *= rel_eps;
							field->Mu[index] *= rel_mu;
						}
					}
				}
				break;
			case MG_CIRCLE:
				rel_eps = materials[m].argv[0].value.floatVal;
				rel_mu = materials[m].argv[1].value.floatVal;

				int cx, cy, R;
				cx = materials[m].argv[2].value.intVal;
				cy = materials[m].argv[3].value.intVal;
				R = materials[m].argv[4].value.intVal;
				
				float d;
				for (int y = 0; y < simulation->height; y++) {
					for (int x = 0; x < simulation->width; x++) {
						index = y * simulation->width + x;
						d = (x - cx) * (x - cx);
						d += (y - cy) * (y - cy);
						if (d < R * R) {
							field->Epsilon[index] *= rel_eps;
							field->Mu[index] *= rel_mu;
						}
					}
				}
				break;
			default:
				break;
		}
	}
	printf("\rApplying material characteristics... done.\n");
}

float randNormalFloat(void) {
	return (float)rand() / (float)RAND_MAX;
}

void glfw_error_callback(int error, const char* desc) {
	// Log GLFW errors to stderr
	fprintf(stderr, "GLFW error %d: %s\n", error, desc);
}

float gaussianPulse(float t, float t0, float spread) {
	return exp(-pow(t / (t0 * spread), 2));
}

void iterateFieldsOnCPU(Field* field, Simulation* simulation) { 
	int index;
	float eps, mu;

   	
	// Update E field
    for (int j = 1; j < simulation->height - 1; j++) {
       	for (int i = 1; i < simulation->width - 1; i++) {
			index = j * simulation->width + i;
			eps = field->Epsilon[index];
       	    field->Ez[index] += (simulation->dt / eps) 
					* ((field->Hy[index] 
					- field->Hy[index - 1]) 
					/ simulation->dx - (field->Hx[index] 
					- field->Hx[(j - 1) * simulation->width + i]) 
					/ simulation->dy)
					- (simulation->dt * field->Sigma[index] * field->Ez[index]
					/ field->Epsilon[index]);
       	}
   	}


	// Update H field
	for (int j = 0; j < simulation->height - 1; j++) {
		for (int i = 0; i < simulation->width; i++) {
			index = j * simulation->width + i;
			mu = field->Mu[index];
			if (i < simulation->width - 1) {
				field->Hx[index] -= simulation->dt 
						/ (mu * simulation->dy) 
						* (field->Ez[index + simulation->width] 
						- field->Ez[index]);
			}
			if (j < simulation->height - 1) {
				field->Hy[index] += simulation->dt 
						/ (mu * simulation->dx) 
						* (field->Ez[index + 1] - field->Ez[index]);
			}
		}
	}
}

void iterateFieldsOnGPU(Field* field, Simulation* simulation) {
	size_t global_size[2] = {simulation->width, simulation->height};

	clEnqueueWriteBuffer(simulation->queue, simulation->Epsilon_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Epsilon, 0, NULL, NULL);	
	clEnqueueWriteBuffer(simulation->queue, simulation->Mu_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Mu, 0, NULL, NULL);	
	clEnqueueWriteBuffer(simulation->queue, simulation->Ez_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Ez, 0, NULL, NULL);	
	clEnqueueWriteBuffer(simulation->queue, simulation->Hx_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Hx, 0, NULL, NULL);	
	clEnqueueWriteBuffer(simulation->queue, simulation->Hy_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Hy, 0, NULL, NULL);	
	clEnqueueWriteBuffer(simulation->queue, simulation->Sigma_kbuf, CL_TRUE,
			0, sizeof(float) * simulation->width * simulation->height,
			field->Sigma, 0, NULL, NULL);

	clSetKernelArg(simulation->E_kernel, 0, sizeof(cl_mem), 
			&simulation->Hx_kbuf);
	clSetKernelArg(simulation->E_kernel, 1, sizeof(cl_mem), 
			&simulation->Hy_kbuf);
	clSetKernelArg(simulation->E_kernel, 2, sizeof(cl_mem), 
			&simulation->Ez_kbuf);
	clSetKernelArg(simulation->E_kernel, 3, sizeof(cl_mem), 
			&simulation->Epsilon_kbuf);
	clSetKernelArg(simulation->E_kernel, 4, sizeof(cl_mem), 
			&simulation->Mu_kbuf);
	clSetKernelArg(simulation->E_kernel, 5, sizeof(float), &simulation->dt);
	clSetKernelArg(simulation->E_kernel, 6, sizeof(float), &simulation->dy);
	clSetKernelArg(simulation->E_kernel, 7, sizeof(float), &simulation->dx);
	clSetKernelArg(simulation->E_kernel, 8, sizeof(int), &simulation->width);
	clSetKernelArg(simulation->E_kernel, 9, sizeof(int), &simulation->height);
	clSetKernelArg(simulation->E_kernel, 10, sizeof(cl_mem), 
			&simulation->Sigma_kbuf);

	clEnqueueNDRangeKernel(simulation->queue, simulation->E_kernel, 2, NULL, 
			global_size, NULL, 0, NULL, NULL);
	clFinish(simulation->queue);
	
	clSetKernelArg(simulation->H_kernel, 0, sizeof(cl_mem), 
			&simulation->Hx_kbuf);
	clSetKernelArg(simulation->H_kernel, 1, sizeof(cl_mem), 
			&simulation->Hy_kbuf);
	clSetKernelArg(simulation->H_kernel, 2, sizeof(cl_mem), 
			&simulation->Ez_kbuf);
	clSetKernelArg(simulation->H_kernel, 3, sizeof(cl_mem), 
			&simulation->Epsilon_kbuf);
	clSetKernelArg(simulation->H_kernel, 4, sizeof(cl_mem), 
			&simulation->Mu_kbuf);
	clSetKernelArg(simulation->H_kernel, 5, sizeof(float), &simulation->dt);
	clSetKernelArg(simulation->H_kernel, 6, sizeof(float), &simulation->dy);
	clSetKernelArg(simulation->H_kernel, 7, sizeof(float), &simulation->dx);
	clSetKernelArg(simulation->H_kernel, 8, sizeof(int), &simulation->width);
	clSetKernelArg(simulation->H_kernel, 9, sizeof(int), &simulation->height);

	clEnqueueNDRangeKernel(simulation->queue, simulation->H_kernel, 2, NULL, 
			global_size, NULL, 0, NULL, NULL);
	clFinish(simulation->queue);

	clEnqueueReadBuffer(simulation->queue, simulation->Ez_kbuf, CL_TRUE, 
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Ez, 0, NULL, NULL);
	clEnqueueReadBuffer(simulation->queue, simulation->Hx_kbuf, CL_TRUE, 
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Hx, 0, NULL, NULL);
	clEnqueueReadBuffer(simulation->queue, simulation->Hy_kbuf, CL_TRUE, 
			0, sizeof(float) * simulation->width * simulation->height, 
			field->Hy, 0, NULL, NULL);
}

void updateFields(Field* field, Simulation* simulation, Source* sources) {
	// Increment simulation time
	simulation->time += simulation->dt;
	simulation->frame++;

	// Add contributions from user-specified sources
	for (int i = 0; i < simulation->sourcec; i++) {
		float sourceVal = 0.0f;
		switch (sources[i].fxn) {
			case SINELINFREQ:
				// Calculate source value for a linear-frequency sinusoid
				sourceVal = sin(2 * M_PI * sources[i].argv[2].value.floatVal 
						* simulation->time 
						+ sources[i].argv[3].value.floatVal);
				int index = simulation->height 
						* sources[i].argv[1].value.intVal
						+ sources[i].argv[0].value.intVal;
				
				// Add source value to the specified field component
				switch (sources[i].fc) {
					case FC_EZ:
						field->Ez[index] += sourceVal;
						break;
					case FC_HX:
						field->Hx[index] += sourceVal;
						break;
					case FC_HY:
						field->Hy[index] += sourceVal;
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}
	}

	if (gpu_support) {
		iterateFieldsOnGPU(field, simulation);
	} else {
		iterateFieldsOnCPU(field, simulation);
	}
	
	if (simulation->boundary_condition == BC_PEC) {
		int index;
		for (int i = 0; i < simulation->height; i++) {
			for (int j = 0; j < simulation->width; j++) {
				index = i * simulation->width + j;
				if (i == 0 || i == simulation->height - 1 || j == 0 
						|| j == simulation->width - 1) {
					field->Ez[index] = 0; 
					field->Hx[index] = 0;
					field->Hy[index] = 0;
				}
			}
		}
	}
}

void visualizeOnCPU(Field* field, Simulation* simulation) { 
	int index;

	// Initialize min and max field values for normalization
	float ezMin;
	float hxMin;
	__attribute__((unused)) float hyMin;
	__attribute__((unused)) float ezMax;
	float hxMax;
	float hyMax;
	ezMin = MAX_FIELD;
	hxMin = MAX_FIELD;
	hyMin = MAX_FIELD;
	ezMax = MIN_FIELD;
	hxMax = MIN_FIELD;
	hyMax = MIN_FIELD;

	// Find the min and max field values
	for (int i = 0; i < simulation->width 
			* simulation->height; ++i) {
		if (field->Ez[i] < ezMin) ezMin = field->Ez[i];
		if (field->Ez[i] < ezMin) ezMin = field->Ez[i];
		if (field->Hx[i] < hxMin) hxMin = field->Hx[i];
		if (field->Hx[i] > hxMax) hxMax = field->Hx[i];
		if (field->Hy[i] > hyMax) hyMax = field->Hy[i];
		if (field->Hy[i] > hyMax) hyMax = field->Hy[i];
	}

	// Normalize field values and store in the image buffer
	for (int y = 0; y < simulation->height; ++y) {
		for (int x = 0; x < simulation->width; ++x) {
			index = y * simulation->width + x;
			float ezVal = field->Ez[index];
			float hxVal = field->Hx[index];
			float hyVal = field->Hy[index];
	
			// Apply user-selected visualization function
			switch (simulation->vis_fxn) {		
				case VIS_TE_1:
					float normVal = (ezVal - -1e1) / (1e2 - -1e1);
					simulation->image[3 * index + 2] = normVal < 0.5 
							? 2 * normVal : 1.0;
					simulation->image[3 * index + 0] = normVal < 0.5 
							? 2 * normVal : 2 * (1 - normVal);
					simulation->image[3 * index + 1] = normVal > 0.5 
							? 2 * (normVal - 0.5) : 0.0;
					break;
				case VIS_TE_2:
					simulation->image[3 * index] = (ezVal*ezVal - MIN_FIELD) 
							/ (MAX_FIELD - MIN_FIELD);
					simulation->image[3 * index + 1] = (hxVal*hxVal 
							- MIN_FIELD) / (MAX_FIELD - MIN_FIELD);
					simulation->image[3 * index + 2] = (hyVal*hyVal 
							- MIN_FIELD) / (MAX_FIELD - MIN_FIELD); 
					break;
				default:
					break;
			}
		}
	}
	
	// If material boundary rendering is enabled, draw them over the image
	if (draw_material_boundaries) {
		for (int i = 0; i < simulation->width * simulation->height; i++) {
			if (simulation->matBoundMask[i] == 1) {
				simulation->image[3 * i] = 0;
				simulation->image[3 * i + 1] = 0;
				simulation->image[3 * i + 2] = 0;
			}
		}
	}
}

void visualizeOnGPU(Field* field, Simulation* simulation) { 
	size_t global_size[2] = {simulation->width, simulation->height};

	cl_int err;
	float minField, maxField;
	switch (simulation->vis_fxn) {
		case VIS_TE_1:
			switch (err = clEnqueueWriteBuffer(simulation->queue, 
					simulation->image_kbuf, CL_TRUE, 0, sizeof(float) 
					* simulation->width * simulation->height * 3, 
					simulation->image, 0 , NULL, NULL)) {
				case CL_SUCCESS:
					break;
				default:
					fprintf(stderr, "Error writing image_kbuf: %d\n", err);
		
			}
			clEnqueueWriteBuffer(simulation->queue, simulation->Hx_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Hx, 0 , NULL, NULL);
			clEnqueueWriteBuffer(simulation->queue, simulation->Hy_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Hy, 0 , NULL, NULL);
			clEnqueueWriteBuffer(simulation->queue, simulation->Ez_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Ez, 0 , NULL, NULL);

			minField = -1e1;
			maxField = 1e2;

			clSetKernelArg(simulation->VIS_TE_1_kernel, 0, sizeof(cl_mem), 
					&simulation->image_kbuf);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 1, sizeof(cl_mem),
					&simulation->Hx_kbuf);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 2, sizeof(cl_mem),
					&simulation->Hy_kbuf);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 3, sizeof(cl_mem),
					&simulation->Ez_kbuf);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 4, sizeof(float),
					&minField);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 5, sizeof(float),
					&maxField);
			clSetKernelArg(simulation->VIS_TE_1_kernel, 6, sizeof(int),
					&simulation->width);

			clEnqueueNDRangeKernel(simulation->queue, 
					simulation->VIS_TE_1_kernel, 2, NULL, global_size, NULL, 
					0, NULL, NULL);
			clFinish(simulation->queue);
			
			clEnqueueReadBuffer(simulation->queue, simulation->image_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height * 3, simulation->image, 0, NULL, 
					NULL);
			break;
		case VIS_TE_2:
			switch (err = clEnqueueWriteBuffer(simulation->queue, 
					simulation->image_kbuf, CL_TRUE, 0, sizeof(float) 
					* simulation->width * simulation->height * 3, 
					simulation->image, 0 , NULL, NULL)) {
				case CL_SUCCESS:
					break;
				default:
					fprintf(stderr, "Error writing image_kbuf: %d\n", err);
		
			}
			clEnqueueWriteBuffer(simulation->queue, simulation->Hx_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Hx, 0 , NULL, NULL);
			clEnqueueWriteBuffer(simulation->queue, simulation->Hy_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Hy, 0 , NULL, NULL);
			clEnqueueWriteBuffer(simulation->queue, simulation->Ez_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height, field->Ez, 0 , NULL, NULL);

			minField = (float)MIN_FIELD;
			maxField = (float)MAX_FIELD;

			clSetKernelArg(simulation->VIS_TE_2_kernel, 0, sizeof(cl_mem), 
					&simulation->image_kbuf);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 1, sizeof(cl_mem),
					&simulation->Hx_kbuf);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 2, sizeof(cl_mem),
					&simulation->Hy_kbuf);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 3, sizeof(cl_mem),
					&simulation->Ez_kbuf);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 4, sizeof(float),
					&minField);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 5, sizeof(float),
					&maxField);
			clSetKernelArg(simulation->VIS_TE_2_kernel, 6, sizeof(int),
					&simulation->width);

			clEnqueueNDRangeKernel(simulation->queue, 
					simulation->VIS_TE_2_kernel, 2, NULL, global_size, NULL, 
					0, NULL, NULL);
			clFinish(simulation->queue);
			
			clEnqueueReadBuffer(simulation->queue, simulation->image_kbuf, 
					CL_TRUE, 0, sizeof(float) * simulation->width 
					* simulation->height * 3, simulation->image, 0, NULL, 
					NULL);
			break;
		default:
			break;
	}

	if (draw_material_boundaries) {
		clEnqueueWriteBuffer(simulation->queue, simulation->image_kbuf, 
				CL_TRUE, 0, sizeof(float) * simulation->width
				* simulation->height * 3, simulation->image, 0, NULL, NULL);
		clEnqueueWriteBuffer(simulation->queue, simulation->matBoundMask_kbuf,
				CL_TRUE, 0, sizeof(float) * simulation->width 
				* simulation->height, simulation->matBoundMask, 0, NULL,
				NULL);	

		clSetKernelArg(simulation->drawMatBounds_kernel, 0, sizeof(cl_mem),
				&simulation->image_kbuf);
		clSetKernelArg(simulation->drawMatBounds_kernel, 1, sizeof(cl_mem),
				&simulation->matBoundMask_kbuf);
		clSetKernelArg(simulation->drawMatBounds_kernel, 2, sizeof(int),
				&simulation->width);
		
		clEnqueueNDRangeKernel(simulation->queue, 
				simulation->drawMatBounds_kernel, 2, NULL, global_size, NULL,
				0, NULL, NULL);
		clFinish(simulation->queue);
		
		clEnqueueReadBuffer(simulation->queue, simulation->image_kbuf, 
				CL_TRUE, 0, sizeof(float) * simulation->width 
				* simulation->height * 3, simulation->image, 0, NULL, NULL);
	}
}

void updateImage(Field* field, Simulation* simulation, Source* sources) { 
	updateFields(field, simulation, sources);	
	
	if (gpu_support) {
		visualizeOnGPU(field, simulation);
	} else {
		visualizeOnCPU(field, simulation);
	}

	// Update OpenGL texture with the new image data
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, simulation->width, 
			simulation->height, 0, GL_RGB, GL_FLOAT, simulation->image);
}

void computeMaterialBoundary(Simulation* simulation, Material* material) {
	int index;
	switch (material->geom) {
		case MG_UNKNOWN:
 			break;
		case MG_TRIANGLE:
			// Extract vertex coordinates
			int x1, y1, x2, y2, x3, y3;
			x1 = material->argv[2].value.intVal;
			y1 = material->argv[3].value.intVal;
			x2 = material->argv[4].value.intVal;
			y2 = material->argv[5].value.intVal;
			x3 = material->argv[6].value.intVal;
			y3 = material->argv[7].value.intVal;
			
			// Compute coefficients for line equations corresponding to 
			// the edges of the triangle
			int A1, B1, C1, A2, B2, C2, A3, B3, C3;
			A1 = y1 - y2;
			B1 = x2 - x1;
			C1 = (x1 * y2) - (x2 * y1);
			A2 = y2 - y3;
			B2 = x3 - x2;
			C2 = (x2 * y3) - (x3 * y2);
			A3 = y3 - y1;
			B3 = x1 - x3;
			C3 = (x3 * y1) - (x1 * y3);
			
			// Store the min and max extent of the triangle edges
			int L1x1, L1x2, L1y1, L1y2;
			int L2x1, L2x2, L2y1, L2y2;
			int L3x1, L3x2, L3y1, L3y2;
			L1x1 = min(x1, x2);
			L1x2 = max(x1, x2);
			L1y1 = min(y1, y2);
			L1y2 = max(y1, y2);
			L2x1 = min(x2, x3);
			L2x2 = max(x2, x3);
			L2y1 = min(y2, y3);
			L2y2 = max(y2, y3);
			L3x1 = min(x1, x3);
			L3x2 = max(x1, x3);
			L3y1 = min(y1, y3);
			L3y2 = max(y1, y3);

			float d1, d2, d3;
			for (int y = 0; y < simulation->height; y++) {
				for (int x = 0; x < simulation->width; x++) {
					index = y * simulation->width + x;
					
					d1 = abs(A1 * x + B1 * y + C1) / sqrt(A1*A1 
							+ B1*B1);
					d2 = abs(A2 * x + B2 * y + C2) / sqrt(A2*A2 
							+ B2*B2);
					d3 = abs(A3 * x + B3 * y + C3) / sqrt(A3*A3 
							+ B3*B3);
					
					if ((d1 < MX_MAT_BOUNDARY_PX 
							&& !(x < L1x1 || x > L1x2 || y < L1y1 
							|| y > L1y2))
							|| (d2 < MX_MAT_BOUNDARY_PX 
							&& !(x < L2x1 || x > L2x2 || y < L2y1 
							|| y > L2y2))
							|| (d3 < MX_MAT_BOUNDARY_PX
							&& !(x < L3x1 || x > L3x2 || y < L3y1 
							|| y > L3y2))) {
						material->boundary[index] = 1;
					} 
				}
			} 
			break;
		case MG_CIRCLE:
			int xc, yc, R;
			xc = material->argv[2].value.intVal;
			yc = material->argv[3].value.intVal;
			R = material->argv[4].value.intVal;
			float d;
			for (int y = 0; y < simulation->height; y++) {
				for (int x = 0; x < simulation->width; x++) {
					index = y * simulation->width + x;
					d = (x - xc) * (x - xc);
					d += (y - yc) * (y - yc);
					d = sqrt(d);
					if (fabs(d - R) < MX_MAT_BOUNDARY_PX) {
						material->boundary[index] = 1;
					}
				}
			}
		default:
			break;	
	}
}

int main(int argc, char** argv) {
	// Ensure a simulation description file has been provided
	if (argc < 2) {
		fprintf(stderr, "Invalid number of arguments.\n");
		fprintf(stderr, "Usage: %s sim_file\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// Initialize the simulation parameter data structure
	Simulation simulation;
	simulation.width = -1;
	simulation.height = -1;
	simulation.time = 0.0f;
	simulation.dx = 1e1;
	simulation.dy = 1e1;
	// Calculate timestep based on grid spacing and c - must satisfy stability
	// condition
	simulation.dt = MX_DT_SCALE / (SPEED_OF_LIGHT * sqrt(1 / (simulation.dx 
			* simulation.dx) + 1 / (simulation.dy * simulation.dy)));
	simulation.sourcec = 0;	
	simulation.vis_fxn = VIS_TE_1;
	simulation.frame = 0;
	simulation.boundary_condition = BC_UNK;
	simulation.pml_layers = -1;
	simulation.pml_conductivity = -1;
	simulation.pml_sigma_polyorder = -1;

	// Open the file for parsing
	const int nsections = MX_SIMDEF_NSEC;
	const char* sections[] = {"[Simulation]", "[Sources]", "[Materials]"};

	Source sources[MX_MAX_SOURCES];  
	Material materials[MX_MAX_MATERIALS];

	for (int s = 0; s < nsections; s++) {
		FILE* sim_file = fopen(argv[1], "r");
		if (sim_file == NULL) {
			fprintf(stderr, "Error opening file at %s.\n", argv[1]);
			exit(EXIT_FAILURE);
		}
		
		char line[MX_SIMFILE_MAX_LINEL]; 
		enum { SEARCH, READ, DONE } state = SEARCH;

		while (fgets(line, sizeof(line), sim_file) && state != DONE) {
			line[strcspn(line, "\n")] = 0;
			if (state == SEARCH) {
				if (strcmp(line, sections[s]) == 0) state = READ;
			} else if (state == READ) {
				if (line[0] == '[' && line[strlen(line) - 1] == ']') {
					state = DONE;
				} else {
					if (strcmp(line, "") == 0) continue;
					if (strcmp(sections[s], "[Simulation]") == 0) {
						char key[MX_SIMFILE_MAX_LINEL];
						char ROL[MX_SIMFILE_MAX_LINEL];
						if (sscanf(line, "%255s %[^\n]", key, ROL) < 1) {
							fprintf(stderr, "Error: Invalid simulation\n");
							fclose(sim_file);
							exit(EXIT_FAILURE);
						}
						if (strcmp(key, "Width") == 0) {
							if (sscanf(ROL, "%d", &simulation.width) != 1) {
								fprintf(stderr, "Error: Invalid format for "
										"Simulation.Width\n");
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
						} else if (strcmp(key, "Height") == 0) {
							if (sscanf(ROL, "%d", &simulation.height) != 1) {
								fprintf(stderr, "Error: Invalid format for "
										"Simulation.Height\n");
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
						} else if (strcmp(key, "ComputeOn") == 0) {
							if (strcmp(ROL, "CPU") == 0) {
								trying_gpu = false;
							}
						} else if (strcmp(key, "Boundary") == 0) {
							if (sscanf(ROL, "%255s %[^\n]", key, ROL) < 1) {
								fprintf(stderr, "Error: Invalid boudnary "
										"conditions.\n");
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							if (strcmp(key, "Natural") == 0) {
								printf("Using natural boundaries.\n");
								simulation.boundary_condition = BC_NAT;
							} else if (strcmp(key, "PEC") == 0) {
								printf("Using PEC boundaries.\n");
								simulation.boundary_condition = BC_PEC;
							} else if (strcmp(key, "PML") == 0) {
								printf("Using PML boundaries.\n");
								simulation.boundary_condition = BC_PML;
								if (strcmp(ROL, "PML") == 0) {
									printf("No arguments specified for PML "
											"boundary - using defaults.\n");
								} else if (sscanf(ROL, "%d %f %d", 
										&simulation.pml_layers,
										&simulation.pml_conductivity,
										&simulation.pml_sigma_polyorder) 
										!= MX_BC_PML_ARGC) {
									fprintf(stderr, "Warning: Improper number"
											" of arguments specified for PML "
											"boundary - using defaults.\n");
									simulation.pml_layers = -1;
									simulation.pml_conductivity = -1;
									simulation.pml_sigma_polyorder = -1;
								}
							}
						} else {
							fprintf(stderr, "Warning: Unknown key: "
									"Simulation.%s - ignoring\n", key);
						}
					} else if (strcmp(sections[s], "[Sources]") == 0) {
						char key[MX_SIMFILE_MAX_LINEL];
						char ROL[MX_SIMFILE_MAX_LINEL];
						if (sscanf(line, "%255s %[^\n]", key, ROL) < 1) {
							fprintf(stderr, "Error: Invalid simulation\n");
							fclose(sim_file);
							exit(EXIT_FAILURE);
						}	
						if (strcmp(key, "SineLinFreq") == 0) {
							char fc_str[MX_FC_STRL]; 
							Source source;
							source.fxn = SINELINFREQ;
							source.argc = MX_SRC_ARGC_SINELINFREQ;
							source.argv[0].type = TYPE_INT;
							source.argv[1].type = TYPE_INT;
							source.argv[2].type = TYPE_FLOAT; 
							source.argv[3].type = TYPE_FLOAT;
							int x, y;
							float f, phi;
							if (sscanf(ROL, "%s %d %d %f %f", fc_str, 
									&x, &y, &f, &phi) != 1 
									+ MX_SRC_ARGC_SINELINFREQ) {
								fprintf(stderr, "Error: Invalid format for "
										"Source #%d: SineLinFreq\n", 
										simulation.sourcec);
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							if (strcmp(fc_str, "Ez") == 0) {
								source.fc = FC_EZ;
							} else if (strcmp(fc_str, "Hx") == 0) {
								source.fc = FC_HX;
							} else if (strcmp(fc_str, "Hy") == 0) {
								source.fc = FC_HY;
							} else {
								fprintf(stderr, "Warning: Unknown field "
										"component for Source #%d - "
										"defaulting to Ez\n", 
										simulation.sourcec);
								source.fc = FC_EZ;
								
							}
							source.argv[0].value.intVal = x;
							source.argv[1].value.intVal = y;
							source.argv[2].value.floatVal = f;
							source.argv[3].value.floatVal = phi;
							
							sources[simulation.sourcec] = source;
							simulation.sourcec++;
						}
					} else if (strcmp(sections[s], "[Materials]") == 0) {
						char key[MX_SIMFILE_MAX_LINEL];
						char ROL[MX_SIMFILE_MAX_LINEL];
						if (sscanf(line, "%255s %[^\n]", key, ROL) < 1) {
							fprintf(stderr, "Error: Invalid simulation\n");
							fclose(sim_file);
							exit(EXIT_FAILURE);
						}
						if (strcmp(key, "Triangle") == 0) {
							Material material;
							material.geom = MG_TRIANGLE;
							material.argc = MX_MAT_ARGC_TRIANGLE;
							material.argv[0].type = TYPE_FLOAT;
							material.argv[1].type = TYPE_FLOAT;
							material.argv[2].type = TYPE_INT;
							material.argv[3].type = TYPE_INT;
							material.argv[4].type = TYPE_INT;
							material.argv[5].type = TYPE_INT;
							material.argv[6].type = TYPE_INT;
							material.argv[7].type = TYPE_INT;
							float rel_eps, rel_mu;
							int x1, x2, x3, y1, y2, y3;
							if (sscanf(ROL, "%f %f %d %d %d %d %d %d", 
									&rel_eps, &rel_mu, &x1, &y1, &x2, &y2, 
									&x3, &y3) != MX_MAT_ARGC_TRIANGLE) {
								fprintf(stderr, "Error: Invalid format for "
										"Material #%d: Triangle\n", 
										simulation.materialc);
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							material.argv[0].value.floatVal = rel_eps;
							material.argv[1].value.floatVal = rel_mu;
							material.argv[2].value.intVal = x1;
							material.argv[3].value.intVal = y1;
							material.argv[4].value.intVal = x2;
							material.argv[5].value.intVal = y2;
							material.argv[6].value.intVal = x3;
							material.argv[7].value.intVal = y3;
						
							material.boundary = (int*)
									calloc(simulation.width 
									* simulation.height, sizeof(int));
							if (material.boundary == NULL) {
								fprintf(stderr, "Failed to allocate memory "
										"for Material #%d's boundary mask.\n",
										simulation.materialc);
								for (int m = 0; m < simulation.materialc; 
										m++) {
									free(materials[m].boundary);
								}
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							computeMaterialBoundary(&simulation, &material);
							materials[simulation.materialc] = material;
							simulation.materialc++;
						}
						if (strcmp(key, "Circle") == 0) {
							Material material;
							material.geom = MG_CIRCLE;
							material.argc = MX_MAT_ARGC_CIRCLE;
							material.argv[0].type = TYPE_FLOAT;
							material.argv[1].type = TYPE_FLOAT;
							material.argv[2].type = TYPE_INT;
							material.argv[3].type = TYPE_INT;
							material.argv[4].type = TYPE_INT;
							float rel_eps, rel_mu;
							int x, y, R;
							if (sscanf(ROL, "%f %f %d %d %d", &rel_eps, 
									&rel_mu, &x, &y, &R) 
									!= MX_MAT_ARGC_CIRCLE) {
								fprintf(stderr, "Error: Invalid format for "
										"Material #%d: Circle\n", 
										simulation.materialc);
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							material.argv[0].value.floatVal = rel_eps;
							material.argv[1].value.floatVal = rel_mu;
							material.argv[2].value.intVal = x;
							material.argv[3].value.intVal = y;
							material.argv[4].value.intVal = R;

							material.boundary = (int*)
									calloc(simulation.width
									* simulation.height, sizeof(int));
							if (material.boundary == NULL) {
								fprintf(stderr, "Failed to allocate memory "
										"for Material #%d's boundary mask.\n",
										simulation.materialc);
								for (int m = 0; m < simulation.materialc;
										m++) {
									free(materials[m].boundary);
								}
								fclose(sim_file);
								exit(EXIT_FAILURE);
							}
							computeMaterialBoundary(&simulation, &material);
							materials[simulation.materialc] = material;
							simulation.materialc++;
						}
					} else {
						fprintf(stderr, "Unknown configuation section\n"); 
					}
				}
			}
		}		
		fclose(sim_file);
	}

	if (simulation.boundary_condition == BC_UNK) {
		fprintf(stderr, "Warning: No boundary conditions specified - "
				"defaulting to natural.\n");
		simulation.boundary_condition = BC_NAT;
	} else if (simulation.boundary_condition == BC_PML) {
		if (simulation.pml_layers == -1) 
				simulation.pml_layers = MX_BC_PML_DEF_LAYERS;
		if (simulation.pml_conductivity == -1) 
				simulation.pml_conductivity = MX_BC_PML_DEF_SIGMA;
		if (simulation.pml_sigma_polyorder == -1) 
				simulation.pml_sigma_polyorder = MX_BC_PML_DEF_SIGPOLYORDER;
	}

	// Initialize GLFW
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize glfw\n");
		exit(EXIT_FAILURE);
	}

	// Create a GLFW window for displaying simulation
	GLFWwindow* window;
	window = glfwCreateWindow(simulation.width, simulation.height, "Maxwell", 
			NULL, NULL);
	if (!window) {
		fprintf(stderr, "Failed to create glfw window\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);

	// Allocate memory for field components
	Field field;
	float* Epsilon = (float*)malloc(simulation.width * simulation.height
			* sizeof(float));
	float* Mu = (float*)malloc(simulation.width * simulation.height
			* sizeof(float));
	float* Ex = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Ey = (float*)malloc(simulation.width * simulation.height
			* sizeof(float));
	float* Ez = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Hx = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Hy = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Hz = (float*)malloc(simulation.width * simulation.height
			* sizeof(float));
	float* Sigma = (float*)malloc(simulation.width * simulation.height
			* sizeof(float));

	if (Epsilon == NULL || Mu == NULL || Ex == NULL || Ey == NULL 
			|| Ez == NULL || Hx == NULL || Hy == NULL || Hz == NULL
			|| Sigma == NULL) {
		fprintf(stderr, "Failed to allocate memory for field object.\n");
		if (Epsilon != NULL) free(Epsilon);
		if (Mu != NULL) free(Mu);
		if (Ex != NULL) free(Ex);
		if (Ey != NULL) free(Ey);
		if (Ez != NULL) free(Ez);
		if (Hx != NULL) free(Hx);
		if (Hy != NULL) free(Hy);
		if (Hz != NULL) free(Hz);
		if (Sigma != NULL) free(Sigma);
		for (int m = 0; m < simulation.materialc; m++) {
			free(materials[m].boundary);
		}
		glfwDestroyWindow(window);
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	field.Epsilon = Epsilon;
	field.Mu = Mu;
	field.Ex = Ex;
	field.Ey = Ey;
	field.Ez = Ez;
	field.Hx = Hx;
	field.Hy = Hy;
	field.Hz = Hz;
	field.Sigma = Sigma;

	// Initialize the field components and add any user-defined materials
	initFields(&field, &simulation);
	addMaterials(&field, &simulation, materials);

	// Allocate memory for simulation image buffer
	simulation.image = (float*)malloc(3 * simulation.width 
			* simulation.height * sizeof(float));
	if (simulation.image == NULL) {
		fprintf(stderr, "Failed to allocate memory for simulation image "
				"buffer.\n");
		free(Epsilon);
		free(Mu);
		free(Ex);
		free(Ey);
		free(Ez);
		free(Hx);
		free(Hy);
		free(Hz);
		for (int m = 0; m < simulation.materialc; m++) {
			free(materials[m].boundary);
		}
		glfwDestroyWindow(window);
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	// Initialize OpenCL
	cl_platform_id platform;
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel E_kernel;
	cl_kernel H_kernel;
	cl_kernel VIS_TE_1_kernel;
	cl_kernel VIS_TE_2_kernel;
	cl_kernel drawMatBounds_kernel;
	cl_int err;

	if (trying_gpu) {
		printf("Attempting to set up GPU... ");
	} else {
		gpu_support = false;
		printf("Running on CPU.\n");
	}

	if (gpu_support) {
		switch (clGetPlatformIDs(1, &platform, NULL)) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error getting platform IDs.\n");
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
		switch(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, 
				NULL)) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error getting device IDs.\n");
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
		context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating OpenCL context.\n");
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
		cl_command_queue_properties properties[] = {
			CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE,
			0
		};
		queue = clCreateCommandQueueWithProperties(context, device, 
				properties, &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating OpenCL command queue.\n");
				gpu_support = false;
		}
	}

	FILE* kernel_file;
	if (gpu_support) {
		const char* kernel_filename = "kernel.cl";
		kernel_file = fopen(kernel_filename, "r");
		if (!kernel_file) {
			fprintf(stderr, "Failed to load kernel source file.\n");
			gpu_support = false;
		}
	}

	size_t kernelSize;
	char* kernelSource;
	if (gpu_support) {
		fseek(kernel_file, 0, SEEK_END);
		kernelSize = ftell(kernel_file);
		rewind(kernel_file);
		kernelSource = (char*)malloc(kernelSize + 1);
	
		if (!kernelSource) {
			fprintf(stderr, "Error allocating memory for kernel source.\n");
			fclose(kernel_file);
			free(Epsilon);
			free(Mu);
			free(Ex);
			free(Ey);
			free(Ez);
			free(Hx);
			free(Hy);
			free(Hz);
			for (int m = 0; m < simulation.materialc; m++) {
				free(materials[m].boundary);
			}
			glfwDestroyWindow(window);
			glfwTerminate();
			exit(EXIT_FAILURE);
		}
		fread(kernelSource, 1, kernelSize, kernel_file);
		kernelSource[kernelSize] = '\0';
		fclose(kernel_file);
	}	

	if (gpu_support) {
		size_t kernelSourceSize = strlen(kernelSource);
		const char* kernelSourceArr[] = {kernelSource};
		program = clCreateProgramWithSource(context, 1, kernelSourceArr, 
				&kernelSourceSize, &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating program from OpenCL kernel."
						"\n");
				free(kernelSource);
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
	    switch (err = clBuildProgram(program, 1, &device, NULL, NULL, NULL)) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error building OpenCL program: %d\n", err);
				size_t log_size;
				clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
						0, NULL, &log_size);
				
				char* log = (char*)malloc(log_size);
				if (log == NULL) {
					fprintf(stderr, "Failed to allocate memory for OpenCL "
							"kernel build log.\n");
					free(kernelSource);
					free(Epsilon);
					free(Mu);
					free(Ex);
					free(Ey);
					free(Ez);
					free(Hx);
					free(Hy);
					free(Hz);
					for (int m = 0; m < simulation.materialc; m++) {
						free(materials[m].boundary);
					}
					glfwDestroyWindow(window);
					glfwTerminate();
					exit(EXIT_FAILURE);
				}
				clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
						log_size, log, NULL);
				fprintf(stderr, "Build log:\n%s\n", log);
				free(kernelSource);
				gpu_support = false;
		}
	}

	if (gpu_support) {
		E_kernel = clCreateKernel(program, "updateEFields", &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating OpenCL kernel: %d\n", err);
				free(kernelSource);
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
		H_kernel = clCreateKernel(program, "updateHFields", &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating OpenCL kernel: %d\n", err);
				free(kernelSource);
				gpu_support = false;
		}
	}

	if (gpu_support) {
		VIS_TE_1_kernel = clCreateKernel(program, "visualizeTE1", &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating visualization kernel: TE1\n");
				free(kernelSource);
				gpu_support = false;
		}
	}

	if (gpu_support) {
		VIS_TE_2_kernel = clCreateKernel(program, "visualizeTE2", &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating visualization kernel: TE2\n");
				free(kernelSource);
				gpu_support = false;
		}
	}

	if (gpu_support) {
		drawMatBounds_kernel = clCreateKernel(program, 
				"drawMaterialBoundaries", &err);
		switch (err) {
			case CL_SUCCESS:
				break;
			default:
				fprintf(stderr, "Error creating material boundary rendering "
						"kernel.\n");
				free(kernelSource);
				gpu_support = false;
		}
	}
	
	if (gpu_support) {
		cl_mem Epsilon_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE, 
				sizeof(float) * simulation.width * simulation.height, NULL, 
				&err);
		cl_mem Mu_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE, 
				sizeof(float) * simulation.width * simulation.height, NULL, 
				&err);
		cl_mem Ez_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE, 
				sizeof(float) * simulation.width * simulation.height, NULL, 
				&err);
		cl_mem Hx_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE, 
				sizeof(float) * simulation.width * simulation.height, NULL, 
				&err);
		cl_mem Hy_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE, 
				sizeof(float) * simulation.width * simulation.height, NULL, 
				&err);
		cl_mem image_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE,
				sizeof(float) * simulation.width * simulation.height * 3, 
				NULL, &err);
		cl_mem matBoundMask_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE,
				sizeof(float) * simulation.width * simulation.height, NULL,
				&err);
		cl_mem Sigma_kbuf = clCreateBuffer(context, CL_MEM_READ_WRITE,
				sizeof(float) * simulation.width * simulation.height, NULL,
				&err);		

		simulation.Epsilon_kbuf = Epsilon_kbuf;
		simulation.Mu_kbuf = Mu_kbuf;
		simulation.Ez_kbuf = Ez_kbuf;
		simulation.Hx_kbuf = Hx_kbuf;
		simulation.Hy_kbuf = Hy_kbuf;
		simulation.image_kbuf = image_kbuf;	
		simulation.matBoundMask_kbuf = matBoundMask_kbuf;
		simulation.Sigma_kbuf = Sigma_kbuf;
		simulation.context = context;
		simulation.queue = queue;
		simulation.program = program;
		simulation.E_kernel = E_kernel;
		simulation.H_kernel = H_kernel;
		simulation.VIS_TE_1_kernel = VIS_TE_1_kernel;
		simulation.VIS_TE_2_kernel = VIS_TE_2_kernel;
		simulation.drawMatBounds_kernel = drawMatBounds_kernel;
	}
	
	if (trying_gpu) {
		if (gpu_support) {
			printf("Good.\n");
		} else {
			printf("Failed. Falling back to CPU.\n");
		}
	} 

	// Aggregate together all material boundaries
	float* matBoundMask = (float*)calloc(simulation.width * simulation.height,
			sizeof(float));
	if (matBoundMask == NULL) {
		fprintf(stderr, "Failed to allocate memory for aggregated material "
				"boundary mask.\n");
		if (gpu_support) free(kernelSource);
		free(Epsilon);
		free(Mu);
		free(Ex);
		free(Ey);
		free(Ez);
		free(Hx);
		free(Hy);
		free(Hz);
		for (int m = 0; m < simulation.materialc; m++) {
			free(materials[m].boundary);
		}
		glfwDestroyWindow(window);
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	for (int m = 0; m < simulation.materialc; m++) {
		for (int i = 0; i < simulation.width * simulation.height; i++) {
			matBoundMask[i] = matBoundMask[i] || materials[m].boundary[i];
		}
	}
	simulation.matBoundMask = matBoundMask;

	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	simulation.start_time = clock();
	
	// Begin main simulation loop
	while (!glfwWindowShouldClose(window)) {
		if (report_framerate) {
			double framerate = ((double) (clock() - simulation.start_time)) 
					/ CLOCKS_PER_SEC;
			framerate = simulation.frame / framerate;
			printf("Simulation averaging %d FPS since last interrupt.\n", 
					(int)framerate);
			report_framerate = false;
		}
		if (cycle_vis) {
			simulation.vis_fxn++;
			if (simulation.vis_fxn == VIS_MAX) simulation.vis_fxn = 0;
			cycle_vis = false;
		}
		if (reset_sim) {
			initFields(&field, &simulation);
			simulation.time = 0.0f;
			updateImage(&field, &simulation, sources);
			reset_sim = false;
		}
		if (just_resumed) {
			simulation.start_time = clock();
			simulation.frame = 0;
			just_resumed = false;
		}
		if (sim_running) updateImage(&field, &simulation, sources);
		
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

	// Clean up, release allocated resources
	glfwDestroyWindow(window);
	glfwTerminate();

	free(Epsilon);
	free(Mu);
	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

	for (int m = 0; m < simulation.materialc; m++) {
		free(materials[m].boundary);
	}

	if (gpu_support) free(kernelSource);

	free(simulation.image);
	free(matBoundMask);

	printf("Goodbye!\n");
		
	exit(EXIT_SUCCESS);
}
