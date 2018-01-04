#pragma once
#include"vector.h"
#include"sph_common.h"
#include"cpu_sph.h"
#include<GL\glut.h>

typedef struct
{
	Vector p;//位置
	Vector n;//法线
	GLfloat val;//建立f（x，y,z）=a方程中的a
} cellvertex;

typedef struct
{
	cellvertex* v[8];//v一般绑定vmem，为什么要多绑定一下我也不是很明白
	cellvertex vmem[8];//vmem[i]相当于一个体元顶点缓冲区
} gridcell;

struct Net_grid
{
	GLuint width;//网格宽方向上有多少个网格
	GLuint height;//同上方向为高
	GLuint depth;//同width方向为深
	Vector	grid_position;//网格基点位置
	GLfloat voxelsize;//一个网格的宽度
	GLfloat* density;//连续密度指针 用于等值面的提取使用 
	gridcell g;//拥有8个顶点的一个小网格
	GLfloat ISO_radius;
};

void init_mc_paramter(Net_grid* N, GLfloat iso_radius, GLfloat grid_Lens, float k, cpu_sph* N_cpu);//准备在这里存入网格长度，iso的半径，申请空间
void sph_render_draw_surface(cpu_sph* N_cpu, float iso_surface,Net_grid* N);
void render_marching_cubes(const Net_grid * N, int x0, int x1, int y0, int y1, int z0, int z1, float iso_value);
void set_cellvertex(cellvertex*& v ,GLfloat* density,GLfloat voxelsize, GLuint width,GLuint height,GLuint xn,GLuint yn,GLuint zn);
void lerp_vertex(Vector* p, Vector* n, cellvertex * v0, cellvertex * v1, GLfloat isolevel);
void polygonize(gridcell* grid,GLfloat threshold);


GLfloat ABS(GLfloat x);