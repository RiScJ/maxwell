#ifndef MAXWELL_H
#define MAXWELL_H

#include <GLFW/glfw3.h>
#include <CL/cl.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SPEED_OF_LIGHT 299792458.0
#define VACUUM_PERMITTIVITY 8.854e-12
#define VACUUM_PERMEABILITY 1.2566e-6

#define DT 1e-2
#define MAX_FIELD 1e2
#define MIN_FIELD 0
#define MX_DT_SCALE 0.9

#define MX_MAX_MATERIALS 1000
#define MX_MAX_MAT_ARGS 10
#define MX_MAT_ARGC_TRIANGLE 8

#define MX_MAT_BOUNDARY_PX 1

#define MX_STRING_ARGL 255
#define MX_MAX_SRC_ARGS 10
#define MX_SIMFILE_MAX_LINEL 256
#define MX_SIMDEF_NSEC 3
#define MX_MAX_SOURCES 1000
#define MX_FC_STRL 10

#define MX_SRC_ARGC_SINELINFREQ 4

typedef enum {
	VIS_TE_LIN_RGB = 0,
	VIS_TE_LIN_EZ_RGB,
	VIS_TE_SQR_RGB,
	VIS_TE_SQR2_RGB,
	VIS_TE_LOG_RGB,
	VIS_MAX
} VisualizationFunction;

typedef struct {
	float* Epsilon;
	float* Mu;
	float* Ex;
	float* Ey;
	float* Ez;
	float* Hx;
	float* Hy;
	float* Hz;
} Field;

typedef struct {
	int width;
	int height;
	float time;
	float dt;
	float dx;
	float dy;
	int sourcec;
	int materialc;
	VisualizationFunction vis_fxn;
	float* image;
	int frame;
	clock_t start_time;
	cl_mem Ez_kbuf;
	cl_mem Hx_kbuf;
	cl_mem Hy_kbuf;
	cl_mem Epsilon_kbuf;
	cl_mem Mu_kbuf;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel E_kernel;
	cl_kernel H_kernel;
} Simulation;

typedef enum {
	TYPE_INT,
	TYPE_DOUBLE,
	TYPE_FLOAT,
	TYPE_STRING
} ArgType;

typedef union {
	int intVal;
	double doubleVal;
	float floatVal;
	char stringVal[MX_STRING_ARGL];
} ArgValue;

typedef struct {
	ArgType type;
	ArgValue value;
} Argument;

typedef enum {
	SINELINFREQ
} SourceFunction;

typedef enum {
	FC_EPS,
	FC_MU,
	FC_EX,
	FC_EY,
	FC_EZ,
	FC_HX,
	FC_HY,
	FC_HZ
} FieldComponent;

typedef struct {
	SourceFunction fxn;
	FieldComponent fc;
	int argc;
	Argument argv[MX_MAX_SRC_ARGS];
} Source;

typedef enum {
	MG_UNKNOWN,
	MG_TRIANGLE
} MaterialGeometry;

typedef struct {
	MaterialGeometry geom;
	float rel_eps;
	float rel_mu;
	int argc;
	Argument argv[MX_MAX_MAT_ARGS];
	int* boundary;
} Material;

#endif
