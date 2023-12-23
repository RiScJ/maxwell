#include "maxwell.h"

bool sim_running = true;
bool reset_sim = false;
bool cycle_vis = false;
bool draw_material_boundaries = true;

int min(int a, int b) {
	return b ^ ((a ^ b) & -(a < b));
}

int max (int a, int b) {
	return a ^ ((a ^ b) & -(a < b));
}

void key_callback(GLFWwindow* window, int key, int __attribute__((unused)) 
		scancode, int action, int mods) {
	if (key == GLFW_KEY_C && mods == GLFW_MOD_CONTROL && (action == GLFW_PRESS 
			|| action == GLFW_REPEAT)) {
		printf("Caught interrupt - exiting...\n");
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		if (sim_running) {
			printf("Pausing simulation.\n");
		} else {
			printf("Resuming simulation.\n");
		}
		sim_running = !sim_running;
	}
	if (key == GLFW_KEY_B && action == GLFW_PRESS) {
		if (draw_material_boundaries) {
			printf("Disabling material boundary rendering.\n");
		} else {
			printf("Enabling material boundary rendering.\n");
		}
		draw_material_boundaries = !draw_material_boundaries;
	}
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		printf("Resetting simulation.\n");
		sim_running = false;
		reset_sim = true;
	}
	if (key == GLFW_KEY_V && action == GLFW_PRESS) {
		printf("Advancing to next visualization function.\n");
		cycle_vis = true;
	}
}

void initFields(Field* field, Simulation* simulation) {
	int index;
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
		}
	}
}

void addMaterials(Field* field, Simulation* simulation, Material* materials) {
	for (int m = 0; m < simulation->materialc; m++) {
		printf("\rApplying material characteristics... (%d/%d)", m, 
			simulation->materialc);
		switch (materials[m].geom) {
			case MG_UNKNOWN:
				break;
			case MG_TRIANGLE:
				float rel_eps, rel_mu;
				rel_eps = materials[m].argv[0].value.floatVal;
				rel_mu = materials[m].argv[1].value.floatVal;

				int x1, y1, x2, y2, x3, y3;
				x1 = materials[m].argv[2].value.intVal;	
				y1 = materials[m].argv[3].value.intVal;	
				x2 = materials[m].argv[4].value.intVal;	
				y2 = materials[m].argv[5].value.intVal;	
				x3 = materials[m].argv[6].value.intVal;	
				y3 = materials[m].argv[7].value.intVal;	
				
				float d1, d2, d3;
				bool has_neg, has_pos;
				int index;
				for (int y = 0; y < simulation->height; y++) {
					for (int x = 0; x < simulation->width; x++) {
						index = y * simulation->width + x;
						
						d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2);
						d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3);
						d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1);

						has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
						has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
						
						if (!(has_neg && has_pos)) {
							field->Epsilon[index] *= rel_eps;
							field->Mu[index] *= rel_mu;
						}
					}
				}
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
	fprintf(stderr, "GLFW error %d: %s\n", error, desc);
}

float gaussianPulse(float t, float t0, float spread) {
	return exp(-pow(t / (t0 * spread), 2));
}

void updateFields(Field* field, Simulation* simulation, Source* sources) {
	simulation->time += simulation->dt;

	// Add sources
	for (int i = 0; i < simulation->sourcec; i++) {
		float sourceVal = 0.0f;
		switch (sources[i].fxn) {
			case SINELINFREQ:
				sourceVal = sin(2 * M_PI * sources[i].argv[2].value.floatVal 
						* simulation->time 
						+ sources[i].argv[3].value.floatVal);
				int index = simulation->height 
						* sources[i].argv[1].value.intVal
						+ sources[i].argv[0].value.intVal;
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

	int index;
	float eps, mu;

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
					/ simulation->dy);
       	}
   	}

    // Apply PEC boundary conditions: E field is zero at all boundaries
    // Top and bottom boundaries
    for (int i = 0; i < simulation->width; i++) {
        field->Ez[i] = 0; // Top boundary
        field->Ez[(simulation->height - 1) * simulation->width + i] = 0;
    }

    // Left and right boundaries
    for (int j = 0; j < simulation->height; j++) {
        field->Ez[j * simulation->width] = 0; // Left boundary
        field->Ez[j * simulation->width + (simulation->width - 1)] = 0; 
    }
}

void updateImage(Field* field, Simulation* simulation, Source* sources, 
			Material* materials) {
	updateFields(field, simulation, sources);	
	
	int index;

	float logMax = log10(MAX_FIELD) - 6;
	float logMin = log10(MIN_FIELD);

	float ezMin, hxMin, hyMin;
	float ezMax, hxMax, hyMax;
	ezMin = MAX_FIELD;
	hxMin = MAX_FIELD;
	hyMin = MAX_FIELD;
	ezMax = MIN_FIELD;
	hxMax = MIN_FIELD;
	hyMax = MIN_FIELD;

	for (int i = 0; i < simulation->width * simulation->height; ++i) {
		if (field->Ez[i] < ezMin) ezMin = field->Ez[i];
		if (field->Ez[i] < ezMin) ezMin = field->Ez[i];
		if (field->Hx[i] < hxMin) hxMin = field->Hx[i];
		if (field->Hx[i] > hxMax) hxMax = field->Hx[i];
		if (field->Hy[i] > hyMax) hyMax = field->Hy[i];
		if (field->Hy[i] > hyMax) hyMax = field->Hy[i];
	}

	for (int y = 0; y < simulation->height; ++y) {
		for (int x = 0; x < simulation->width; ++x) {
			index = y * simulation->width + x;
			float ezVal = field->Ez[index];
			float hxVal = field->Hx[index];
			float hyVal = field->Hy[index];
	
			switch (simulation->vis_fxn) {		
				case VIS_TE_LIN_RGB:
					simulation->image[3 * index + 2] = (ezVal - ezMin) 
							/ (ezMax - ezMin);
					simulation->image[3 * index + 1] = (hyVal - hyMin) 
							/ (hyMax - hyMin);
					simulation->image[3 * index] = (hxVal - hxMin) 
							/ (hxMax - hxMin);
					break;
				case VIS_TE_LIN_EZ_RGB:
					float normVal = (ezVal - ezMin) / (ezMax - ezMin);
					simulation->image[3 * index + 2] = normVal < 0.5 
							? 2 * normVal : 1.0;
					simulation->image[3 * index + 0] = normVal < 0.5 
							? 2 * normVal : 2 * (1 - normVal);
					simulation->image[3 * index + 1] = normVal > 0.5 
							? 2 * (normVal - 0.5) : 0.0;
					break;
				case VIS_TE_SQR_RGB:
					simulation->image[3 * index + 1] = (hyVal*hyVal 
							- hyMin*hyMin) / (hyMax*hyMax - hyMin*hyMin); 
					simulation->image[3 * index] = (hxVal*hxVal 
							- hxMin*hxMin) / (hxMax*hxMax - hxMin*hxMin);
					simulation->image[3 * index + 2] = (ezVal*ezVal 
							- ezMin*ezMin) / (ezMax*ezMax - ezMin*ezMin); 
					break;				
				case VIS_TE_SQR2_RGB:
					simulation->image[3 * index] = (ezVal*ezVal - MIN_FIELD) 
							/ (MAX_FIELD - MIN_FIELD);
					simulation->image[3 * index + 1] = (hxVal*hxVal 
							- MIN_FIELD) / (MAX_FIELD - MIN_FIELD);
					simulation->image[3 * index + 2] = (hyVal*hyVal 
							- MIN_FIELD) / (MAX_FIELD - MIN_FIELD); 
					break;
				case VIS_TE_LOG_RGB:
					simulation->image[3 * index] = (log10(ezVal) - logMin) 
							/ (logMax - logMin);
					simulation->image[3 * index + 1] = (log10(hxVal) - logMin) 
							/ (logMax - logMin);
					simulation->image[3 * index + 2] = (log10(hyVal) - logMin) 
							/ (logMax - logMin);
					break;
				default:
					break;
			}
		}
	}
	if (draw_material_boundaries) {
		for (int m = 0; m < simulation->materialc; m++) {
			for (int y = 0; y < simulation->height; y++) {
				for (int x = 0; x < simulation->width; x++) {
					index = y * simulation->width + x;
					if (materials[m].boundary[index] == 1) {
						simulation->image[3 * index] = 0;
						simulation->image[3 * index + 1] = 0;
						simulation->image[3 * index + 2] = 0;
					}
				}
			}	
		}
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, simulation->width, 
			simulation->height, 0, GL_RGB, GL_FLOAT, simulation->image);
}

void computeMaterialBoundary(Simulation* simulation, Material* material) {
	for (int m = 0; m <= simulation->materialc; m++) {
		switch (material->geom) {
			case MG_UNKNOWN:
	 			break;
			case MG_TRIANGLE:
				int x1, y1, x2, y2, x3, y3;
				x1 = material->argv[2].value.intVal;
				y1 = material->argv[3].value.intVal;
				x2 = material->argv[4].value.intVal;
				y2 = material->argv[5].value.intVal;
				x3 = material->argv[6].value.intVal;
				y3 = material->argv[7].value.intVal;
				
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
	
				int index;
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
			default:
				break;	
		}
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
	simulation.dt = MX_DT_SCALE / (SPEED_OF_LIGHT * sqrt(1 / (simulation.dx 
			* simulation.dx) + 1 / (simulation.dy * simulation.dy)));
	simulation.sourcec = 0;	
	simulation.vis_fxn = VIS_TE_SQR2_RGB;

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
					} else {
						fprintf(stderr, "Unknown configuation section\n"); 
					}
				}
			}
		}		
		fclose(sim_file);
	}

	// Initialize GLFW
	glfwSetErrorCallback(glfw_error_callback);

	GLFWwindow* window;
	
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize glfw\n");
		return -1;
	}

	window = glfwCreateWindow(simulation.width, simulation.height, "Maxwell", 
			NULL, NULL);
	if (!window) {
		fprintf(stderr, "Failed to create glfw window\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);


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

	if (Epsilon == NULL || Mu == NULL || Ex == NULL || Ey == NULL 
			|| Ez == NULL || Hx == NULL || Hy == NULL || Hz == NULL) {
		fprintf(stderr, "Failed to allocate memory for field object.\n");
		if (Epsilon != NULL) free(Epsilon);
		if (Mu != NULL) free(Mu);
		if (Ex != NULL) free(Ex);
		if (Ey != NULL) free(Ey);
		if (Ez != NULL) free(Ez);
		if (Hx != NULL) free(Hx);
		if (Hy != NULL) free(Hy);
		if (Hz != NULL) free(Hz);
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

	initFields(&field, &simulation);
	addMaterials(&field, &simulation, materials);

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

	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	while (!glfwWindowShouldClose(window)) {
		if (cycle_vis) {
			simulation.vis_fxn++;
			if (simulation.vis_fxn == VIS_MAX) simulation.vis_fxn = 0;
			cycle_vis = false;
		}
		if (reset_sim) {
			initFields(&field, &simulation);
			simulation.time = 0.0f;
			updateImage(&field, &simulation, sources, materials);
			reset_sim = false;
		}
		if (sim_running) updateImage(&field, &simulation, sources, materials);
		
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

	free(simulation.image);

	printf("Goodbye!\n");
		
	exit(EXIT_SUCCESS);
}
