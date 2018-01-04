#pragma once
#include <vector>
#include <list>
#include "glg.h"
using namespace std;

class Particle
{
public:
	Particle(int k);
	~Particle() {};

	void SetK_phase(int k);
	void SetPosition(vector3 position);
	void SetVel(vector3 vel_temp);
	void SetA_k(vector<float> ak);
	void SetDensity(float d) { density = d; }
	void SetPressure(float p) { pressure = p; }
	void SetMass(float m) { mass = m; }
	void SetDriftVel(vector3 umk, int k);
	void SetDriftVel(vector<vector3> umk) { Um_k = umk; }
	void SetViscosity(float Viscosity) { m_Viscosity = Viscosity; }
	void SetstaticDensity(float a) { staticDensity = a; }
	vector3 GetPosition() { return pos; };
	vector3 GetVel() { return vel; };
	vector<float> GetA_k() { return A_k; };
	float GetDensity() { return density; }
	float GetPressure() { return pressure; }
	float GetMass() { return mass; }
	vector<vector3> GetAllDriftVel() { return Um_k; }
	vector3 GetKDriftVel(int k) { if (k < Um_k.size()) return Um_k[k]; }

	//先这样吧，以后优化才改为private

public:
	int k_phase;//相的数目

	//单相参数
	vector3 normal;//法向
	vector3 pos;//粒子位置
	vector3 vel;//速度
	vector3 vel_half;//半速度
	vector3 acc;//加速度
	float density;//密度
	float pressure;//压力
	float mass;//质量
	float m_Viscosity; // 粘度
	float staticDensity;//静态密度
	vector3 P_m_Gradient;

	//多相参数
	vector<float> Density_k;//每相的密度
	vector<float> A_k;//储存体积分数，存储空间num_phase*粒子的索引
	vector<float> C_k;//储存质量分数，储存同理与体积分数
	vector<vector3> Um_k;//漂移速度
	vector<float> P_k;//每相的压力

	vector<float> A_k_Divergence;//散度
	vector<float> P_k_Divergence;
	vector<vector3> P_k_Gradient;//梯度
	vector<vector3> A_k_Gradient;
	vector<float> Vel_Divergence;//速度的散度

	vector3 temp_A;//论文公式（10）的那个参数
};

class Fluid
{
public:
	Fluid(int num, int k);
	~Fluid();
	Particle* GetParticleByIndex(int index);
	int GetParticleNum() { return m_ParticleVector.size(); };
	void DeleteParticleRoomspace();

public:
	int m_num;//粒子数量
	vector<Particle*> m_ParticleVector;
};



