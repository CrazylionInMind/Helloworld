#pragma once
#include"vector.h"
#include"sph_common.h"
#include"cpu_sph.h"
#include<GL\glut.h>

typedef struct
{
	Vector p;//λ��
	Vector n;//����
	GLfloat val;//����f��x��y,z��=a�����е�a
} cellvertex;

typedef struct
{
	cellvertex* v[8];//vһ���vmem��ΪʲôҪ���һ����Ҳ���Ǻ�����
	cellvertex vmem[8];//vmem[i]�൱��һ����Ԫ���㻺����
} gridcell;

struct Net_grid
{
	GLuint width;//����������ж��ٸ�����
	GLuint height;//ͬ�Ϸ���Ϊ��
	GLuint depth;//ͬwidth����Ϊ��
	Vector	grid_position;//�������λ��
	GLfloat voxelsize;//һ������Ŀ��
	GLfloat* density;//�����ܶ�ָ�� ���ڵ�ֵ�����ȡʹ�� 
	gridcell g;//ӵ��8�������һ��С����
	GLfloat ISO_radius;
};

void init_mc_paramter(Net_grid* N, GLfloat iso_radius, GLfloat grid_Lens, float k, cpu_sph* N_cpu);//׼��������������񳤶ȣ�iso�İ뾶������ռ�
void sph_render_draw_surface(cpu_sph* N_cpu, float iso_surface,Net_grid* N);
void render_marching_cubes(const Net_grid * N, int x0, int x1, int y0, int y1, int z0, int z1, float iso_value);
void set_cellvertex(cellvertex*& v ,GLfloat* density,GLfloat voxelsize, GLuint width,GLuint height,GLuint xn,GLuint yn,GLuint zn);
void lerp_vertex(Vector* p, Vector* n, cellvertex * v0, cellvertex * v1, GLfloat isolevel);
void polygonize(gridcell* grid,GLfloat threshold);


GLfloat ABS(GLfloat x);