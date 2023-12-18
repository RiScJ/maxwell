#include "maxwell.h"

bool sim_running = true;
bool reset_sim = false;
float time_dir = 1;

void key_callback(GLFWwindow* window, int key, int scancode, int action, 
		int mods) {
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
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		printf("Resetting simulation.\n");
		sim_running = false;
		reset_sim = true;
	}
}

void initFields(Field* field, Simulation* simulation) {
	int index;
	for (int y = 0; y < simulation->height; ++y) {
		for (int x = 0; x < simulation->width; ++x) {
			index = y * simulation->width + x;
			field->Ez[index] = 0;
			field->Hx[index] = 0;
			field->Hy[index] = 0;
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

void updateFields(Field* field, Simulation* simulation, Source* sources) {
	const float eps = 8.854e-12;
	const float mu = 1.2566e-6;

	simulation->time += time_dir * simulation->dt;

	float t0 = 1e-9;
	float spread = 1e-2;
	int sourceX = simulation->width / 2;
	int sourceY = simulation->height / 2;

	float frequency = 1.5e6;
	float omega = 2 * M_PI * frequency;

	// Add sources
	for (int i = 0; i < simulation->sourcec; i++) {
		float sourceVal = 0.0f;
		switch (sources[i].fxn) {
			case SINELINFREQ:
				sourceVal = sin(2 * M_PI * sources[i].argv[2].value.floatVal 
						* simulation->time 
						+ sources[i].argv[3].value.floatVal);
				int index = simulation->height * sources[i].argv[1].value.intVal
						+ sources[i].argv[0].value.intVal;
				switch (sources[i].fc) {
					case FC_EZ:
						field->Ez[index] += time_dir * sourceVal;
						break;
					case FC_HX:
						field->Hx[index] += time_dir * sourceVal;
						break;
					case FC_HY:
						field->Hy[index] += time_dir * sourceVal;
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}
	}

	// Update H field
	for (int j = 0; j < simulation->height - 1; j++) {
		for (int i = 0; i < simulation->width; i++) {
			int index = j * simulation->width + i;
			if (i < simulation->width - 1) {
				field->Hx[index] -= simulation->dt / (mu * simulation->dy) 
						* (field->Ez[index + simulation->width] 
						- field->Ez[index]);
			}
			if (j < simulation->height - 1) {
				field->Hy[index] += simulation->dt / (mu * simulation->dx) 
						* (field->Ez[index + 1] - field->Ez[index]);
			}
		}
	}
   	
	// Update E field
    for (int j = 1; j < simulation->height - 1; j++) {
       	for (int i = 1; i < simulation->width - 1; i++) {
       	    field->Ez[j * simulation->width + i] += (simulation->dt / eps) 
					* ((field->Hy[j * simulation->width + i] 
					- field->Hy[j * simulation->width + i - 1]) 
					/ simulation->dx - (field->Hx[j * simulation->width + i] 
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
        field->Ez[j * simulation->width + (simulation->width - 1)] = 0; // Right boundary
    }
}

void updateImage(Field* field, Simulation* simulation, Source* sources) {
	updateFields(field, simulation, sources);	
	
	printf("%f\n", field->Ez[512 * simulation->width + 512]);
	float* visualData = (float*)malloc(3 * simulation->width * simulation->height * sizeof(float));
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
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, simulation->width, 
			simulation->height, 0, GL_RGB, GL_FLOAT, visualData);
	free(visualData);
}

int main(int argc, char** argv) {
	// Ensure a simulation description file has been provided
	if (argc != 2) {
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

	// Open the file for parsing
	const int nsections = MX_SIMDEF_NSEC;
	const char* sections[] = {"[Simulation]", "[Sources]"};

	Source sources[MX_MAX_SOURCES];  

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
					if (sections[s] == "[Simulation]") {
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
					} else if (sections[s] == "[Sources]") {
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
							source.argc = MX_SRC_ARGC_SINELINFREQ ;
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

	window = glfwCreateWindow(simulation.width, simulation.height, "Example", 
			NULL, NULL);
	if (!window) {
		fprintf(stderr, "Failed to create glfw window\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);


	Field field;
	float* Ez = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Hx = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	float* Hy = (float*)malloc(simulation.width * simulation.height 
			* sizeof(float));
	field.Ez = Ez;
	field.Hx = Hx;
	field.Hy = Hy;

	initFields(&field, &simulation);

	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	while (!glfwWindowShouldClose(window)) {
		if (reset_sim) {
			initFields(&field, &simulation);
			simulation.time = 0.0f;
			updateImage(&field, &simulation, sources);
			reset_sim = false;
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

	glfwDestroyWindow(window);
	glfwTerminate();

	free(Ez);
	free(Hx);
	free(Hy);

	printf("Goodbye!\n");
		
	exit(EXIT_SUCCESS);
}
