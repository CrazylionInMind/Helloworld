#ifndef _SPH_COMMON_H_ 
#define _SPH_COMMON_H_ 
 
#include "glg.h" 
 
#define INITIAL_GRID_CAPACITY 60 
 
typedef int particle_index; 
 
typedef struct 
{ 
	void* prealloc; 
	int used_size; 
	int capacity; 
} sph_mem_pool; 
 
typedef struct 
{ 
	int index; 
	float distsq; 
} sph_neighbour; 
 
typedef struct 
{ 
	int* sizes; 
	sph_neighbour** p; /** **/ 
	void* pool; 
	int used_by_rigid; // Neighbour list memory used by rigid body particles 
	int n_poolused; // Used size of memory pool 
} sph_neighbour_list; 
 
 
typedef struct 
{ 
	particle_index** particles; /** particles indexed by grid **/ 
	int* sizes;      /** number of particles grid has **/ 
	int* caps;       /** capacity of each grid **/ 
	int capacity; 
	int width; 
	int height; 
	int depth; 
	float grid_len;  /** length of grid edge */ 
 
	float minx; 
	float miny; 
	float minz; 
	sph_mem_pool mempool; 
} sph_grid; 
 
typedef struct 
{ 
	int gx; 
	int gy; 
	int gz; 
} grid_vector3i; 
 
void sph_mem_pool_create(sph_mem_pool* pool, int size); 
void* sph_mem_pool_alloc(sph_mem_pool* pool, int size); 
void sph_mem_pool_free(sph_mem_pool* pool); 
 
void sph_neighbour_list_create(sph_neighbour_list* l, int n_particles); 
 
void sph_grid_create(sph_grid* g, int n_particles, float grid_len); 
void sph_grid_clear(sph_grid* g, const vector3* pos, int start, int end); 
void sph_grid_alloc(sph_grid* g, const vector3* pos, int start, int end); 
void sph_grid_get_neighbours(sph_grid* g, const vector3* pos, int start, int end, 
							 sph_neighbour_list* n_list, float search_radius); 
#endif
