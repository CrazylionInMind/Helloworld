#pragma once

#include <iostream>
#include<list>
#include "Fluid.h"
#include "sph_common.h"
#include "Neighbourlist.h"
#include "Single.h"
#include "Param.h"
using std::list;
class Calculate_Reverce;
class SingleCal;

class SPHNewSystem :public Single<SPHNewSystem>
{
public:
	void init();
	SPHNewSystem();
	~SPHNewSystem();

	void SetFluid(Fluid *temp) { m_fluid = temp; }
	void SetParam(float b, float c, float d, float f);
	void CreateNighbourList();
	void GridToNeighbour();
	void CalSphOneStep();
	vector3* GetParticlePosByIndex(int num);
	matrix4 mat_col;
	matrix4 mat_inv_col;
	Neighbour_list n_list; //�ڽӱ�
	Grid grid; //Grid for neighbour search 

	Fluid* GetFluid() { return m_fluid; }
	float GetStiff() { return m_Stiff; }
	float GetTimeStep() { return m_Timestep; }
	float GetH2() { return m_Smoothlen*m_Smoothlen; }
	float GetH() { return m_Smoothlen; }
private:
	Fluid *m_fluid;
	Calculate_Reverce *m_Calculate;
	SingleCal *m_SingleCal;
	//ϵͳ����

	float m_Smoothlen; // �⻬����
	float m_Timestep;  // ʱ�䲽
	float m_Stiff;     // ����
	float m_R_Search;  // �ڽӱ������뾶
};

//ע���������еļ��㶼��Ҫ����
class Calculate_Reverce
{
public:
	Calculate_Reverce();
	~Calculate_Reverce() = default;

	virtual void CalPressure(SPHNewSystem &sph,int whattime);//�����ܶ�
	virtual void CalDriftVel(SPHNewSystem &sph);//����Ư���ٶ�
	virtual void CalFraction(SPHNewSystem &sph);//�����������
	virtual void CalForce(SPHNewSystem &sph);//��������
	virtual void compute_col(vector3* col, const vector3* vel, const vector3* n, float diff, float stiff, float damp);
	virtual void glass_collision(vector3* p, vector3* col, const vector3* vel, const matrix4* mat, const matrix4* mat_inv, float radius, float stiff, float damp);/*ģ��Ͳ�������ײ */
	virtual void CalCollision(SPHNewSystem &sph);//������ײ
	virtual void CalPosition(SPHNewSystem &sph);//����������λ��

	void AddFrame() { _Frame++; }
	int GetFrame() {
		return _Frame;
	}
	void setH(float temph);
protected:
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
	int _Frame;//��¼�ڼ��μ���Ĳ���
};
