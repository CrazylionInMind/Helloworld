#pragma once

#include <iostream>
#include<list>
#include "Fluid.h"
#include "sph_common.h"
#include "Neighbourlist.h"
#include "Single.h"
using std::list;
class Calculate;
class SingleCal;

class SPHSystem :public Single<SPHSystem>
{
public:
	SPHSystem();
	~SPHSystem();

	void SetFluid(Fluid *temp) { m_fluid = temp; }
	void SetParam(float b, float c, float d, float f);
	void CreateNighbourList();
	void GridToNeighbour();
	void CalSphOneStep();
	void CalSphforsinglemodel();
	vector3* GetParticlePosByIndex(int num);
	matrix4 mat_col;
	matrix4 mat_inv_col;
	Neighbour_list n_list; //邻接表
	Grid grid; //Grid for neighbour search 

	Fluid* GetFluid() { return m_fluid; }
	float GetStiff() { return m_Stiff; }
	float GetTimeStep() { return m_Timestep; }
private:
	Fluid *m_fluid;
	Calculate *m_Calculate;
	SingleCal *m_SingleCal;
	//系统参数

	float m_Smoothlen; // 光滑长度
	float m_Timestep;  // 时间步
	float m_Stiff;     // 刚性
	float m_R_Search;  // 邻接表搜索半径
};

class Calculate
{
public:
	Calculate();
	~Calculate() = default;

	virtual void CalPressure(SPHSystem &sph);//计算密度
	virtual void CalDriftVel(SPHSystem &sph);//计算漂移速度
	virtual void CalFraction(SPHSystem &sph);//计算体积分数
	virtual void CalForce(SPHSystem &sph);//计算受力
	virtual void compute_col(vector3* col, const vector3* vel, const vector3* n, float diff, float stiff, float damp);
	virtual void glass_collision(vector3* p, vector3* col, const vector3* vel, const matrix4* mat, const matrix4* mat_inv, float radius, float stiff, float damp);/*模拟和玻璃杯碰撞 */
	virtual void CalCollision(SPHSystem &sph);//计算碰撞
	virtual void CalPosition(SPHSystem &sph);//最后计算最后的位置

	void AddFrame() { _Frame++; }
	int GetFrame() {
		return _Frame;
	}
	void setH(float temph);
protected:
	vector3 SPH_function(Particle i, Particle j, vector3 jfunc, vector3 ifunc = zeros());
	float SPH_function(Particle i, Particle j, float jfunc, float ifunc = 0.0f);
	vector3 SPH_function_gradient(Particle i, Particle j, float jfunc, float ifunc = 0.0f);
	float SPH_function_gradient(Particle i, Particle j, vector3 jfunc, vector3 ifunc = zeros(), float Adjustment = 1.0f);
	vector3 SPH_function_gradient_2(Particle i, Particle j, vector3 jfunc, vector3 ifunc = zeros(), float Adjustment = 1.0f);
	const float PI = 3.1415926535f;
	float h;
	float poly6_coef;
	float grad_poly6_coef;
	float lap_poly6_coef;
	float grad_spiky_coef;
	float lap_vis_coef;
	float _tao;
	float sigma;
	int _Frame;//记录第几次计算的参数
};

class SingleCal :public Calculate
{
public:
	virtual void CalPressure(SPHSystem &sph);
	virtual void CalForce(SPHSystem &sph);
};