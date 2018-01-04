#pragma once

#include <vector>
#include "sph_common.h"
#include "Fluid.h"
using std::vector;

class Neighbour_list
{
public:
	void SetSize(int m);

	int m_sizes;
	vector<vector<sph_neighbour>> p;
};

class Grid
{
public:
	Grid();
	~Grid();
	void SetGrid_len(float len) { grid_len = len; };//�趨���񳤶���Ҫ��CreateGrid֮ǰ
	void CreateGrid(const Fluid input);
	

	vector<vector<int>> particles; /** particles indexed by grid **/
	//vector<int> sizes;      /** number of particles grid has **/
	vector<int> caps;       /** capacity of each grid **/
	int capacity;
	int width;
	int height;
	int depth;
	float grid_len;  /** length of grid edge */

	float minx;
	float miny;
	float minz;
};


