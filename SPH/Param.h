#pragma once

#define NUM_PHASE 2

#define SMOOTHING_LENGTH 0.01f    
#define SEARCH_RADIUS (1*SMOOTHING_LENGTH)    
#define N_LOOPS 2   
#define CUBE_LEN_X 10    
#define CUBE_LEN_Y 10    
#define LEN_Z 20
#define N_PARTICLES (CUBE_LEN_X*CUBE_LEN_Y*LEN_Z)    
#define EPSILON 0.0001f    
#define VISCOSITY 0.5f//0.2f    
#define STIFF 1.5f    
#define MASS 0.00020543f    
#define RANDF() ((rand()/(GLfloat)RAND_MAX - 0.5f))
#define ISO_RADUIS 0.0115f
#define MC_GRID_LEN (0.5*SMOOTHING_LENGTH)
#define TIMESTEP 0.005f
#define STATICDENSITY 1000.0f//¾²Ë®ÃÜ¶È
#define N_STEPS 2