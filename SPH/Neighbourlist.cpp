#include "Neighbourlist.h"
#include<limits>
#include <iostream>

Grid::Grid()
{
	capacity = 40 * 40 * 40;//最多容纳
}

Grid::~Grid()
{
}

void Grid::CreateGrid(Fluid input)
{
	float fmin_x;
	float fmax_x;
	float fmin_y;
	float fmax_y;
	float fmin_z;
	float fmax_z;
	fmin_x = fmin_y = fmin_z = numeric_limits<float>::max();//1000;    
	fmax_x = fmax_y = fmax_z = numeric_limits<float>::min();//-1000;   
	int n = input.GetParticleNum();
	for (int i = 0; i < n; i++)
	{
		const Particle* temp = input.GetParticleByIndex(i);
		if (!temp)
			continue;
		if (fmin_x > temp->pos.x)
			fmin_x = temp->pos.x;
		if (fmax_x < temp->pos.x)
			fmax_x = temp->pos.x;
		if (fmin_y > temp->pos.y)
			fmin_y = temp->pos.y;
		if (fmax_y < temp->pos.y)
			fmax_y = temp->pos.y;
		if (fmin_z > temp->pos.z)
			fmin_z = temp->pos.z;
		if (fmax_z < temp->pos.z)
			fmax_z = temp->pos.z;
	}

	//构建网格属性
	width = (int)((fmax_x - fmin_x + grid_len) / grid_len);
	height = (int)((fmax_y - fmin_y + grid_len) / grid_len);
	depth = (int)((fmax_z - fmin_z + grid_len) / grid_len);

	minx = fmin_x;
	miny = fmin_y;
	minz = fmin_z;
	if (width*height*depth >= 10000)
	{
		cout << "构造的网格太多，现在需要构造的网格数为:"<< (width*height*depth) << endl;
		system("pause");
		exit(0);
	}
	particles.clear();
	particles.resize(width*height*depth);

	//把粒子添加到网格中
	for (int i = 0; i < input.GetParticleNum(); i++)
	{
		int gx;
		int gy;
		int gz;
		int gindex;

		const Particle* temp = input.GetParticleByIndex(i);
		gx = (int)((temp->pos.x - fmin_x) / grid_len);
		gy = (int)((temp->pos.y - fmin_y) / grid_len);
		gz = (int)((temp->pos.z - fmin_z) / grid_len);

		gindex = gx + gy * width + gz * width * height;

		if (particles[gindex].size() > 100)
			std::cout << "Warm:Grid No." << gindex << " 's particles is over 100." << std::endl;

		particles[gindex].push_back(i);//particles[i][j]表示第i个网格中第j个粒子的索引
	}
}

void Neighbour_list::SetSize(int m)
{
	if(p.size() < m)
	{
		p.clear();//清除等下重新创建
		m_sizes = m;
		p.resize(m_sizes);
	}

}
