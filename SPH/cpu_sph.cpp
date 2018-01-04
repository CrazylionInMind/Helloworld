#include <MEMORY.H>    
#include <MATH.H>    
#include <WINDOWS.H>    
#include <ASSERT.H>    

#include "cpu_sph.h"    


//#define RADIUS 0.03f     
const float PI = 3.1415926535f;

#define POS(i) sph->pos[i]    
#define VEL(i) sph->vel[i]    
#define ACC(i) sph->acc[i]    
#define VELH(i) sph->vel_half[i]    

void mult_create(cpu_sph* sph, int n, float *a_k)
{
	//��n��
	sph->num_phase = n;
	int num = sph->n_particles;//���ӵ���Ŀ
	sph->A_k = (float *)malloc(n*num * sizeof(float));//�����n�����������Ŀռ�
	sph->C_k = (float *)malloc(n*num * sizeof(float));//��������
	sph->U_mk = (vector3 *)malloc(n*num * sizeof(vector3));//ƫ���ٶ�

	sph->P_k = (float *)malloc(n*num * sizeof(float));//����k��ѹ��
	sph->change_for_volue = (float *)malloc(n*num * sizeof(float));//��������������仯

	for (int i = 0; i < n*num; i++)
	{
		sph->A_k[i] = a_k[i];
		sph->U_mk[i].x = sph->U_mk[i].y = sph->U_mk[i].z = 0.0f;
		sph->P_k[i] = 0.0f;
	}

	sph->sa_A_k = (float *)malloc(sph->num_phase *sph->n_particles * sizeof(float));//���������ɢ��
	sph->sa_P_k = (float *)malloc(sph->num_phase *sph->n_particles * sizeof(float));//ÿ��ѹ����ɢ��
	sph->xiangliang_for_Ak = (vector3 *)malloc(n*sph->n_particles * sizeof(vector3));
	sph->xiangliang_for_Pk = (vector3 *)malloc(n*sph->n_particles * sizeof(vector3));
	for (int i = 0; i < num; i++)
	{
		sph->xiangliang_for_Ak[i].x = sph->xiangliang_for_Ak[i].y = sph->xiangliang_for_Ak[i].z = 0.0f;
		sph->xiangliang_for_Pk[i].x = sph->xiangliang_for_Pk[i].y = sph->xiangliang_for_Pk[i].z = 0.0f;
	}
}

void cpu_sph_create(cpu_sph* sph, //����SPH����    
	int size, const vector3* pos, const vector3* vel,
	float smoothlen, float viscosity, float mass, float stiff, float r_search, int n_loops)
{
	int i;
	float h = smoothlen;
	float* mem;
	sph->viscosity = viscosity;
	sph->n_particles = size;
	sph->smoothlen = smoothlen;
	sph->n_loops = n_loops;
	sph->stiff = stiff;
	
	sph_neighbour_list_create(&sph->n_list, size);//ΪSPH��list����ռ䣬��sph��sizeͬ��С    

	sph_grid_create(&sph->grid, sph->n_particles, h);//Ϊ�������ռ�    
	sph->grid.grid_len = r_search;
	sph->r_search = r_search;

	// Allocate memory for attributes of particles    
	mem = (float*)malloc(49 * size * sizeof(float));
	sph->mass = (float*)&mem[0];
	sph->density = (float*)&mem[size];
	sph->pressure = (float*)&mem[2 * size];
	sph->pos = (vector3*)&mem[3 * size];
	sph->vel = (vector3*)&mem[6 * size];
	sph->vel_half = (vector3*)&mem[9 * size];
	sph->acc = (vector3*)&mem[12 * size];
	sph->normal = (vector3*)&mem[15 * size];
	sph->curvature = (float*)&mem[18 * size];
	sph->mark = (int*)&mem[19 * size];
	memset(sph->mark, 0, size);
	sph->fluid_start = 0;
	sph->n_fluidp = size - sph->fluid_start;

	for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)
	{
		sph->mass[i] = mass;
	}

	memcpy(sph->pos, pos, sph->n_particles * sizeof(vector3));
	memcpy(sph->vel, vel, sph->n_particles * sizeof(vector3));
	memcpy(sph->vel_half, vel, sph->n_particles * sizeof(vector3));

	// Precompute kernel coefficients    
	sph->poly6_coef = 315.0f / (64.0f*PI*(float)pow(h, 9));
	sph->grad_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	sph->lap_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	sph->grad_spiky_coef = -45.0f / (PI*h*h*h*h*h*h);
	sph->lap_vis_coef = 45.0f / (PI*h*h*h*h*h*h);

	/**** FOR TEST ****/
	mat4_set_identity(&sph->mat_col);//�Ѳ���������0    
	sph_grid_clear(&sph->grid, sph->pos, 0, sph->n_particles - 1);
}

void cpu_sph_transform_obstacles(cpu_sph* sph, const matrix4* m)//�����ϰ�����Ϣ���������ӵ�λ��    
{
	mat4_mul(&sph->mat_col, m, &sph->mat_col);
	mat4_invert(&sph->mat_inv_col, &sph->mat_col);
}

void cpu_sph_transform_obstacles(SPHNewSystem* sph, const matrix4* m)//�����ϰ�����Ϣ���������ӵ�λ��    
{
	mat4_mul(&sph->mat_col, m, &sph->mat_col);
	mat4_invert(&sph->mat_inv_col, &sph->mat_col);
}

float SPH_function(cpu_sph *sph, int i,int j,float jfunc,float ifunc=0.0f)
{
	//�ܶȵļ�����ʱ����
	float result;
	vector3 dis;
	vec3_set(&dis, 0, 0, 0);

	vec3_sub(&dis, &POS(i), &POS(j));
	float r2 = vec3_dot(&dis, &dis);

	float h = sph->smoothlen;
	float h2 = h*h;


	float temp = 0.0f;
	temp = sph->mass[j] / sph->density[j] * sph->poly6_coef*(h2-r2)*(h2 - r2)*(h2 - r2);
	result =temp* (ifunc + jfunc);

	return result;
}

vector3 SPH_function_gradient(cpu_sph *sph, int i, int j, float jfunc,  float ifunc=0)
{
	vector3 result;
	vec3_set(&result,0,0,0);

	vector3 dis;
	vec3_set(&dis, 0, 0, 0);

	vec3_sub(&dis, &POS(i), &POS(j));
	float r = sqrtf(vec3_dot(&dis, &dis));//���ӵľ���
	float h = sph->smoothlen;
	
	float temp = sph->mass[j]*(1.0/sph->density[j])* (jfunc+ifunc);
	temp *= sph->grad_spiky_coef*(h - r)*(h - r) / r;
	vec3_scale(&result, temp, &dis);	

	return result;
}

float SPH_function_gradient(cpu_sph *sph, int i, int j, vector3 jfunc, vector3 ifunc = zeros())//ɢ��
{
	float result = 0.0f;

	vector3 dis;
	vec3_set(&dis, 0, 0, 0);

	vec3_sub(&dis, &POS(i), &POS(j));
	float r = sqrtf(vec3_dot(&dis, &dis));
	float h = sph->r_search;

	float temp = sph->mass[j] * (1.0 / sph->density[j])*sph->grad_spiky_coef*(h - r)*(h - r) / r;

	vector3 res;
	vec3_add(&res, &jfunc, &ifunc);

	result = vec3_dot(&res, &vec3_scale(temp, &dis));
	return result;
}

void cpu_sph_compute_density(cpu_sph* sph, bool first)//�����ܶȣ���ģ�⺯������    
{
	int i;
	int j;

	int nindex;
	sph_neighbour_list* nlist;

	memset(sph->density, 0, sph->n_particles * sizeof(float));

	nlist = &sph->n_list;

	if (first)
	{
		for (i = 0; i < sph->n_particles; i++)
		{
			for (j = 0; j < nlist->sizes[i]; j++)
			{
				float density;
				float distsq;
				float h2 = sph->smoothlen*sph->smoothlen;

				nindex = nlist->p[i][j].index;
				distsq = nlist->p[i][j].distsq;

				if (h2 > distsq)
				{
					float h2_r2;
					h2_r2 = h2 - distsq;
					density = h2_r2*h2_r2*h2_r2;
					sph->density[i] += sph->mass[nindex] * density;
					if (i != nindex)
					{
						sph->density[nindex] += sph->mass[i] * density;
					}
				}
			}
		}
	}
	else
	{
		for (i = 0; i < sph->n_particles; i++)
		{
			for (j = 0; j < nlist->sizes[i]; j++)
			{
				float density;
				float distsq;
				float h2 = sph->smoothlen*sph->smoothlen;

				nindex = nlist->p[i][j].index;
				distsq = vec3_distsq(&POS(i), &POS(nindex));
				nlist->p[i][j].distsq = distsq;

				if (h2 > distsq)
				{
					float h2_r2;
					h2_r2 = h2 - distsq;
					density = h2_r2*h2_r2*h2_r2;
					sph->density[i] += sph->mass[nindex] * density;
					if (i != nindex)
						sph->density[nindex] += sph->mass[i] * density;
				}
			}
		}
	}


	for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)
	{
		sph->density[i] *= sph->poly6_coef;
		sph->pressure[i] = sph->stiff*(sph->density[i] - 1000.0f);

		//C_k ��Ϊ���������ܵ�ʱ����Ҫ����ɫ�ı仯
		//�����赥����ܶ����ϵ��ܶ�һ��,����Ϊ����ˮ���ӣ�ֻ����ɫֵ��һ��
		//��������ֵ�����������������������ͬ��
		//��k�����������
		for (int j = 0; j < sph->num_phase; j++)
		{
			sph->C_k[i + j*sph->n_particles] = sph->A_k[i + j*sph->n_particles];
		}

		//�ָ��һ�� ����P_k
		for (int j = 0; j < sph->num_phase; j++)
		{
			//��ǰ����i�е�ÿһ���ѹ��p_k
			sph->P_k[i + j*sph->n_particles] = sph->pressure[i] * sph->A_k[i + j*sph->n_particles];
		}

		//ע�⣡���������Ѿ�ȡ����
		//sph->density[i] = 1.0f / sph->density[i];
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//ע��ѵ���ȥ����
	}
}


void cpu_sph_compute_force(cpu_sph* sph)//������������ģ�⺯������    
{
	int i;
	int j;
	float h;

	vector3 vdiff;

	sph_neighbour_list* nlist;

	nlist = &sph->n_list;

	h = sph->smoothlen;
	memset(sph->acc, 0, sph->n_particles * sizeof(vector3));

	for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)
	{
		for (j = 1; j < nlist->sizes[i]; j++)
		{
			vector3 force;
			float r;
			int nindex;

			nindex = nlist->p[i][j].index;

			r = glg_sqrt(nlist->p[i][j].distsq);


			if (r < h)
			{
				float h_r;

				vector3 diff;


				h_r = h - r;
				vec3_sub(&diff, &POS(i), &POS(nindex));


				vec3_scale(&force, -0.5f*(sph->pressure[i] + sph->pressure[nindex])*sph->grad_spiky_coef*h_r / r, &diff);

				vec3_sub(&vdiff, &VEL(nindex), &VEL(i));

				vec3_scale(&vdiff, sph->viscosity*sph->lap_vis_coef, &vdiff);

				vec3_add(&force, &force, &vdiff);

				vec3_scale(&force, h_r/sph->density[i] / sph->density[nindex], &force);

				// Apply force    
				vec3_scaleadd(&ACC(i), &ACC(i), sph->mass[nindex], &force);
				vec3_scaleadd(&ACC(nindex), &ACC(nindex), -sph->mass[i], &force);
			}
		}
	}
}


_inline void compute_col(vector3* col,
	const vector3* vel, const vector3* n,
	float diff, float stiff, float damp)
{
	float v0;
	float reverse;

	v0 = vec3_dot(n, vel);//���ڻ�    
	reverse = stiff*diff - damp*v0;
	vec3_scaleadd(col, col, reverse, n);
}


void cpu_sph_glass_collision(vector3* p, vector3* col, const vector3* vel,//ģ��Ͳ�������ײ    
	const matrix4* mat, const matrix4* mat_inv,
	float radius, float stiff, float damp)
{
	float GLASS_R = 0.05f;
	float GLASS_BOTTOM = 0.0f;
	float GLASS_TOP = 10.22f;

	vector3 p_col;
	vector3 n;
	float diff; /** Distance between object and fluid particle **/

	vec3_set(&n, 0.0f, 0.0f, 0.0f);
	mat4_mulvec3(&p_col, mat_inv, p);// �������ӵ�ǰλ�ü���p_col    

	diff = 2.0f*radius - (GLASS_R - (float)sqrt(p_col.x*p_col.x + p_col.y*p_col.y));//����diff    

	if (((diff < 8.0f*radius) && (diff > EPSILON)) && (p_col.z < GLASS_TOP))
	{
		vec3_set(&n, -p_col.x, -p_col.y, 0.0f);
		vec3_normalize(&n, &n);
		mat4_mulvec3_as_mat3(&n, mat, &n);
		compute_col(col, vel, &n, diff, stiff, damp);//����diff����col    
	}

	diff = 2.0f*radius - (p_col.z - GLASS_BOTTOM);

	if (diff > EPSILON)
	{
		vec3_set(&n, 0.0f, 0.0f, 1.0f);
		mat4_mulvec3_as_mat3(&n, mat, &n);
		compute_col(col, vel, &n, diff, stiff, damp);
	}
}

void cpu_sph_process_collision(cpu_sph* sph)//��ײ    
{
	int i;
	//  float e= 1.0f;    
	float sphere_radius = 0.004f;
	float stiff = 30000.0f;
	float damp = 128.0f;

	_ALIGNED vector3 pre_p; /** Predicted pos **/
	vector3 col;
	vec3_set(&pre_p, 0, 0, 0);

	for (i = 0; i < sph->n_particles; i++)
	{
		vec3_set(&col, 0.0f, 0.0f, 0.0f);//������colȫ����0    
		vec3_scaleadd(&pre_p, &POS(i), sph->timestep, &VELH(i));//���㵱ǰʱ�����ӵ�λ��    
		cpu_sph_glass_collision(&pre_p, &col, &VEL(i), &sph->mat_col, &sph->mat_inv_col, sphere_radius, stiff, damp);//�����벣����ײ    
		vec3_add(&ACC(i), &ACC(i), &col);//������ӣ��������Ӽ��ٶ�    
	}

}
//�����ָ幫ʽ��һ����벿��
void cpu_sph_compute_drift(cpu_sph *sph)
{
	sph_neighbour_list* nlist;

	float h = sph->smoothlen;

	nlist = &sph->n_list;
	//�洢p��ѹ�����ݶ�
	vector3 *sa_p_k = (vector3 *)malloc(sph->n_particles*sph->num_phase * sizeof(vector3));
	//�洢ÿ�����ӵ���������仯�ݶ�
	vector3 *sa_a_k = (vector3 *)malloc(sph->n_particles*sph->num_phase * sizeof(vector3));
	for (int m = 0; m < sph->n_particles*sph->num_phase; m++)
	{
		vec3_set(&sa_p_k[m], 0, 0, 0);
		vec3_set(&sa_a_k[m], 0, 0, 0);
		vec3_set(&sph->xiangliang_for_Ak[m], 0, 0, 0);
		vec3_set(&sph->xiangliang_for_Pk[m], 0, 0, 0);//ÿ�����¼���
	}
	int n = sph->n_particles;//��֮����������

	//����ÿ������
	for (int i = sph->fluid_start; i < sph->fluid_start + n; i++)
	{
		//����ڽӱ�
		for (int j = 1; j < nlist->sizes[i]; j++)
		{
			float r;
			int nindex;
			float h_r;
			nindex = nlist->p[i][j].index;

			r = glg_sqrt(nlist->p[i][j].distsq);
			vector3 temp_v1, temp_v2;
			vector3 temp_a1, temp_a2;
			if (r < h)
			{
				//����д��ĺ��������⻬�Ĺ�ʽ
				vector3 diff;
				h_r = h - r;

				vec3_sub(&diff, &POS(i), &POS(nindex));

				for (int k = 0; k < sph->num_phase; k++)
				{
					//û��ʹ�õ�����ʽ��������������������������������������������
					vector3 temp_v1 = SPH_function_gradient(sph, i, nindex, sph->P_k[nindex + k*n]);//��ѹ��ʹ�õ�����Ĺ�ʽ�������������
					vec3_scaleadd(&sph->xiangliang_for_Pk[i + k*n], &sph->xiangliang_for_Pk[i + k*n], 1, &temp_v1);

					//������������ݶȣ���������Ǳ�����
					vec3_set(&temp_a1, 0.0f, 0.0f, 0.0f);

					
					/*if (sph->frame == 17&&i==591)
					{
						printf_s("This is 17 frame");
					}*/
					//temp_a1 = SPH_function_gradient(sph, i, nindex, sph->A_k[nindex + k*n],-sph->A_k[i+k*n]);//ƽ���أ����Ч����Щ����
					temp_a1 = SPH_function_gradient(sph, i, nindex, sph->A_k[nindex + k*n]);
					vec3_scaleadd(&sph->xiangliang_for_Ak[i + k*n], &sph->xiangliang_for_Ak[i + k*n], 1, &temp_a1);
					//if (i==591&&k==0)
					//{
					//	printf_s("591����������ݶ��ڵ�һ�%f,%f,%f......��������˸ı��ǵ�%d�����Ӳ�����Ӱ��\n\n",sph->xiangliang_for_Ak[i + k*n].x, sph->xiangliang_for_Ak[i + k*n].y, sph->xiangliang_for_Ak[i + k*n].z,nindex);
					//}

					//�Ϳ������費��Ҫƽ����,��Ҫ�Ļ�����������
					//�����ƽ�����µ����������A_k�����Ӷ�һЩ

				}

			}
		}
	}
	//ÿ��i���Ӽ������ÿ�������ݶȺ�ѹ���ݶȺ�
	for (int i = 0; i < n; i++)
	{
		float tao = 0.0000001f;//tao��ȡֵ������������������ 10e-8��10e-6
		float o = 0.0001f;//�ĸ���ȡֵΪ10e-4��10e-3

		float temp1 = 0.0f;
		//float temp2 = 0.0f;
		//float temp3 = 0.0f;
		vector3 temp2;
		vector3 temp3;
		temp2.x = 0.0f;temp2.y = 0.0f;temp2.z = 0.0f;//������0


		//��ʽ����Ư���ٶ�u_mk,�ȼ����ۼ�ֵ
		for (int k = 0; k < sph->num_phase; k++)
		{
			temp1 += sph->C_k[i + k*n] * sph->density[i];//����ǰ��ļ��裬sph->density[i]�ͱ�ʾΪk����ܶ�
			//temp2 += sph->C_k[i + j*n] * sph->sa_P_k[i + j*n];
			//temp3 += sph->sa_A_k[i + j*n] * sph->C_k[i + j*n] / sph->A_k[i + j*n];
			vec3_scaleadd(&temp2, &temp2, sph->C_k[i + k*n], &sph->xiangliang_for_Pk[i + k*n]);

			//��ʱ��������ɢ�����������û�е�����

			//if (sph->A_k[i + k*n]<=0)
			//{
			//	vector3 zero;
			//	vec3_set(&zero, 0, 0, 0);
			//	vec3_scaleadd(&temp3, &temp3, 1.0f, &zero);
			//}
			//else
			//{
			//	vec3_scaleadd(&temp3, &temp3, sph->C_k[i + k*n] / sph->A_k[i + k*n], &sph->xiangliang_for_Ak[i + k*n]);
			//}
			////vec3_scaleadd(&temp3, &temp3, sph->C_k[i + k*n] / sph->A_k[i + k*n], &sph->xiangliang_for_Ak[i + k*n]);
		}


		for (int k = 0; k < sph->num_phase; k++)//��һ��������������forѭ��
		{
			//��u_mk��һ�����
			//�ұ߹�ʽ��һ��
			temp3.x = 0.0f;temp3.y = 0.0f;temp3.z = 0.0f;//������0
			vector3 acc, gravity;
			vec3_set(&gravity, -0.0f, -0.0f, -0.0f);
			//vec3_set(&gravity, -0.0f, -0.0f, -9.8f);
			vec3_add(&acc, &ACC(i), &gravity);//��Ϊsph�ݴ�ļ��ٶ���û������
			vec3_scale(&sph->U_mk[i + k*n], tao*(sph->density[i] - temp1), &acc);//��Ϊǰ�����������ÿһ��k���ܶȶ����ڻ�����ӵ��ܶ�

			//�ұ߹�ʽ�ڶ���
			vector3 temp4;
			vec3_set(&temp4, 0.0f, 0.0f, 0.0f);
			vec3_sub(&temp4, &sph->xiangliang_for_Pk[i + k*n], &temp2);
			vec3_scaleadd(&sph->U_mk[i + k*n], &sph->U_mk[i + k*n], -tao, &temp4);

			vec3_sub(&sph->U_mk[i + k*n], &sph->U_mk[i + k*n], &sph->vel[i]);
		//�ұ߹�ʽ������(������ʱ������)
			//temp2�Ѿ��ù������ˣ���temp2���洢�µļ����ݴ�
			//        1/����������ܵ������������
			//if (sph->A_k[i + k*n]==0)
			//{
			//	vec3_set(&temp2, 0.0f, 0.0f, 0.0f);
			//}
			//else
			//{
			//	vec3_scale(&temp2, 1.0 / sph->A_k[i + k*n], &sph->xiangliang_for_Ak[i + k*n]);
			//}
			//vec3_sub(&temp3, &temp2, &temp3);
			//vec3_scaleadd(&sph->U_mk[i + k*n], &sph->U_mk[i + k*n], -o, &temp3);//�õ��������Ư���ٶ�

		}
	}
	//�洢p��ѹ�����ݶ�
	free(sa_p_k);
	free(sa_a_k);
}

void cpu_sph_pingliu_fracion(cpu_sph *sph)
{
	int i, j;
	float h;
	int n = sph->n_particles;
	vector3 vdiff;

	sph_neighbour_list* nlist;

	nlist = &sph->n_list;

	h = sph->smoothlen;
	int num = sph->num_phase*sph->n_particles;

	float *zancun = new float[num];
	float *pingliu_a_k = new float[num];
	for (int m = 0; m < num; m++)
	{
		pingliu_a_k[m] = 0.0f;
		zancun[m] = 0.0f;
	}

	for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)
	{
		vector3 temp1, temp2, temp3;//����ÿһ������
		temp1.x = temp1.y = temp1.z = 0.0f;
		temp2.x = temp2.y = temp2.z = 0.0f;
		temp3.x = temp3.y = temp3.z = 0.0f;
		float t1 = 0.0f, t2[2] = { 0.0f ,0.0f }, t3 = 0.0f;

		for (j = 1; j < nlist->sizes[i]; j++)
		{
			float r;
			int nindex;
			vector3 ttttemp;
			nindex = nlist->p[i][j].index;

			r = glg_sqrt(nlist->p[i][j].distsq);

			if (r < h)
			{
				float h_r;
				vector3 diff;
				h_r = h - r;
				vec3_sub(&diff, &POS(i), &POS(nindex));//�ݴ��������				

				//�ָ�ڶ����й�ʽ1��2��3
				for (int k = 0; k < sph->num_phase; k++)
				{
					//��ʽ��һ���һ�μ���
					vec3_set(&ttttemp, 0.0f, 0.0f, 0.0f);
					vec3_sub(&ttttemp, &sph->vel[nindex], &sph->vel[i]);//ע��ǰ��д���ˣ������ǻ���ٶȵĲ�ֵ
					temp1 = vec3_scale(0.5*(sph->A_k[i + k*n] + sph->A_k[nindex + k*n]), &ttttemp);
					t1 = 0.0f;
					t1 = SPH_function_gradient(sph, i, nindex, temp1);//���Ų�ʹ��t1�ĵ�����ʽ
					//t1 = SPH_function_gradient(sph,i,nindex,vec3_scale(0.5*(sph->A_k[i + k*n] + sph->A_k[nindex + k*n]),&sph->vel[nindex]));

					//�ұ�ʽ�ӵ�2��
					t2[k] = 0.0f;
					//��ʽ�ڶ��һ�μ�����
					t2[k] = SPH_function_gradient(sph, i, nindex, vec3_scale(sph->A_k[nindex + k*n], &sph->U_mk[nindex + k*n]), vec3_scale(sph->A_k[i + k*n], &sph->U_mk[i + k*n]));
					zancun[i + k*n] = 0.0f;


					//��߹�ʽ�ڶ���
					zancun[i + k*n] = vec3_dot(&sph->vel[i], &sph->xiangliang_for_Ak[i + k*n]);

					//vec3_add(&zancun_uki[i + k*n], &zancun_uki[i + k*n], &zancun_akumk[i + k*n]);
					float tempadd = t1 + t2[k] + zancun[i + k*n];
					pingliu_a_k[i + k*n] += tempadd;
					
					if (i == 591)
					{
						//printf_s("��%d֡����������ֱ�Ϊ%f,%f\n", sph->frame, sph->A_k[591], sph->A_k[4591]);
						if (k == 0 && (tempadd > 0.001|| tempadd<-0.001))
						{
							printf_s("����Ӱ�췧ֵ:%f,�����ӵ����Ϊ%d\n", tempadd,nindex);
							sph->mark[nindex] = 1;
						}
						//printf_s("pingliu_a_k[ %d +n]:%f\n",i,pingliu_a_k[i + 1 * n]);
						//printf_s("\n\n");
					}
					//Ϊʲôƽ�ⶼ��û��Ч���ģ�ע���������������Ч��һ��
					/*if (i!=nindex)
					{
						pingliu_a_k[nindex + k*n] += -t1 - t2[k] - zancun[i + k*n];
					}*/
					
					/*if ((nindex==1767||i==1767)&&k==1&&t2[k]!=0)
					{
						printf_s("t2 is change\n");
					}*/

				}
			}
		}//����j����

	}
	for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)
	{
		for (int k = 0; k < sph->num_phase; k++)
		{
			sph->A_k[i + k*n] += pingliu_a_k[i + k*n] * sph->timestep;
			//if (i==675)
			//{
			//	printf_s("����675��%d�%f",k, sph->A_k[i + k*n]);
			//	if (k == 1)
			//		printf_s("\n\n");
			//}
		}
	}

	//�������У��
	float *temp_change = new float[num];
	for (int i = 0; i < sph->n_particles; i++)
	{
		float a_i = 0.0f;
		for (int j = 0; j < sph->num_phase; j++)
		{
			if (sph->A_k[i + j*n] < 0.0000000000001)
			{
				sph->A_k[i + j*n] = 0.0;
			}//���С��0����Ϊ0
			else
				a_i += sph->A_k[i + j*n];
		}//����������������ۼ�

		//���������ԭ����һ�������� ÿ����������������0
		//�о�����ΪС��0����� ����ֱ�ӽ����������Ϊ0 ����ʱ�䲽�ֱȽϴ󣬵����˿������඼�����������Ϊ0�����

		//��������
		for (int j = 0; j < sph->num_phase; j++)
		{
			if (a_i < 0.00000001)
			{
				temp_change[i + j*n] = sph->A_k[i + j*n] - 0.0;
				sph->A_k[i + j*n] = 0.0f;
				continue;
			}
			temp_change[i + j*n] = sph->A_k[i + j*n] - sph->A_k[i + j*n] / a_i;
			sph->A_k[i + j*n] = sph->A_k[i + j*n] / a_i;
		}

		//if (i <= 2100&&i>=1900)
		//{
		//	printf_s("����%d�������������Ϊ%f and %f\n", i, sph->A_k[i], sph->A_k[i + sph->n_fluidp]);
		//}

		//ѹ��У��

		//��ʱ������ѹ��

	}
	delete[] temp_change;
	delete[] pingliu_a_k;
	delete[] zancun;
}

//float tx,ty=0;    
//matrix4 m;    
void cpu_sph_elapse(cpu_sph* sph, float t,int Frame)
{
	int i;
	int j;
	_ALIGNED vector3 gravity;

	sph->frame = Frame;
	sph->timestep = t;
	vec3_set(&gravity, -0.0f, -0.0f, -9.8f);

	for (j = 0; j < sph->n_loops; j++)
	{
		if (j == 0)
		{
			sph_grid_clear(&sph->grid, sph->pos, 0, sph->n_particles - 1);
			sph_grid_alloc(&sph->grid, sph->pos, 0, sph->fluid_start - 1);
			sph_grid_get_neighbours(&sph->grid, sph->pos, sph->fluid_start, sph->n_particles - 1, &sph->n_list, sph->r_search);
		}

		cpu_sph_compute_density(sph, (j == 0));//�����ܶȣ�����֮��ľ���   

		cpu_sph_compute_drift(sph);//����Ư���ٶ�
		cpu_sph_pingliu_fracion(sph);//ƽ���������

		cpu_sph_compute_force(sph);//��������    
		cpu_sph_process_collision(sph);//������ײ�������&������    

		//tx += 0.000001f;    
		//ty += 0.00001f;    
		//mat4_set_translate(&m, -tx, 0, -ty);    
		//cpu_sph_transform_obstacles(sph, &m);    
		for (i = sph->fluid_start; i < sph->fluid_start + sph->n_fluidp; i++)// ����ģ��    
		{
			_ALIGNED vector3 v_half;
			_ALIGNED vector3 final_acc;

			vec3_add(&final_acc, &sph->acc[i], &gravity);//�������ٶȸı���ٶ�    
			vec3_scaleadd(&v_half, &VELH(i), t, &final_acc);//���ٶȸı��ٶ�    
			vec3_scaleadd(&POS(i), &POS(i), t, &v_half);//�ٶȸı�λ��    
			vec3_add(&VEL(i), &VELH(i), &v_half);//�ٶȱ仯    
			vec3_scale(&VEL(i), 0.5f, &VEL(i));//�ٶȼ���    
			vec3_copy(&VELH(i), &v_half);//�ٶ�    

		}
	}
}

void cpu_sph_get_pos(cpu_sph* sph, vector3* pos)
{
	memcpy(pos, sph->pos, sph->n_particles * sizeof(vector3));
}