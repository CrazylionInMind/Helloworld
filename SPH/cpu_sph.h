#ifndef _SPH_CPU_H_ 
#define _SPH_CPU_H_ 
 
#include "sph_common.h" 
#include "SPHSystem.h"
typedef struct 
{ 
	int start; 
	int end; 
	matrix4 rot; // Rotation matrix 
	vector3 rg;  // Center of gravity 
} sph_rigid_body; 
 
typedef struct 
{ 
	// Global parametesr of fluid 
	float viscosity; // Viscosity of fluid 
	float smoothlen; // Smoothing length 
	float timestep;  // Timestep 
	float stiff;     // Stiffness 
	int n_particles; // Number of particles; 
	float r_search;  // Search radius of neighbour list 
 
	// Region of particles attribures 
	int fluid_start;   // Primary index of fluid particle 
	int n_fluidp;      // Number of fluid particles 
 
	// Attributes of particles 
	vector3* normal;   // Normal 
	vector3* pos;      // Position 
	vector3* vel;      // Velocity 
	vector3* vel_half; // Velocity half steps forwarded 
	vector3* acc;      // Acceleration 
	float* curvature;  // Curvature at particle location���� 
	float* density;    // Density and its reciprocal at particle location 
	float* pressure;   // Pressure at particle location 
	float* mass;       // Mass of particle 
	int  * mark;//	ÿ����ı��
	sph_neighbour_list n_list; // Neighbour list 
 
	sph_grid grid; // Grid for neighbour search 
 
	int n_loops; /* n_loops equation-solving for 1 neighbour search performed */ 
 
	// Coefficients for kernel 
	float poly6_coef; 
	float lap_poly6_coef; 
	float grad_poly6_coef; 
	float grad_spiky_coef; 
	float lap_vis_coef; 
 
	// Collision handling 
	matrix4 mat_col; 
	matrix4 mat_inv_col; 

	//��Ӷ���ģ��Ĳ���
	float num_phase;//�����Ŀ
	float *A_k;//��������������洢�ռ�num_phase*���ӵ�����
	float *C_k;//������������������ͬ�����������
	vector3 *U_mk;//Ư���ٶ�
	float *P_k;//ÿ���ѹ��
	float *change_for_volue;
	float *sa_A_k;
	float *sa_P_k;
	vector3 *xiangliang_for_Pk;
	vector3 *xiangliang_for_Ak;
	int frame;
} cpu_sph; 
 
/** 
 * Create sph simulation context 
 */ 
void cpu_sph_create(cpu_sph* sph,//���� 
					int size, const vector3* pos, const vector3* vel, 
					float smoothlen, float viscosity, float mass, float stiff, float search_radius, int n_loops); 
 
void cpu_sph_elapse(cpu_sph* sph, float t,int Frame);//ģ�� 
void cpu_sph_transform_obstacles(cpu_sph* sph, const matrix4* m);//�ƶ������� 
void cpu_sph_transform_obstacles(SPHNewSystem * sph, const matrix4* m);
void cpu_sph_get_pos(cpu_sph* sph, vector3* pos);//��ȡ����λ�� 
void mult_create(cpu_sph* sph, int n,float *a_k);
#endif