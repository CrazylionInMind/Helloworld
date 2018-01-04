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

	//�������ɣ��Ժ��Ż��Ÿ�Ϊprivate

public:
	int k_phase;//�����Ŀ

	//�������
	vector3 normal;//����
	vector3 pos;//����λ��
	vector3 vel;//�ٶ�
	vector3 vel_half;//���ٶ�
	vector3 acc;//���ٶ�
	float density;//�ܶ�
	float pressure;//ѹ��
	float mass;//����
	float m_Viscosity; // ճ��
	float staticDensity;//��̬�ܶ�
	vector3 P_m_Gradient;

	//�������
	vector<float> Density_k;//ÿ����ܶ�
	vector<float> A_k;//��������������洢�ռ�num_phase*���ӵ�����
	vector<float> C_k;//������������������ͬ�����������
	vector<vector3> Um_k;//Ư���ٶ�
	vector<float> P_k;//ÿ���ѹ��

	vector<float> A_k_Divergence;//ɢ��
	vector<float> P_k_Divergence;
	vector<vector3> P_k_Gradient;//�ݶ�
	vector<vector3> A_k_Gradient;
	vector<float> Vel_Divergence;//�ٶȵ�ɢ��

	vector3 temp_A;//���Ĺ�ʽ��10�����Ǹ�����
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
	int m_num;//��������
	vector<Particle*> m_ParticleVector;
};



