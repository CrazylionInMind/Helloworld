#include "Fluid.h"

Particle::Particle(int k)
{
	SetK_phase(k);
}

void Particle::SetK_phase(int k)
{
	k_phase = k;
	A_k.resize(k_phase);
	C_k.resize(k_phase);
	Um_k.resize(k_phase);
	P_k.resize(k_phase);

	A_k_Divergence.resize(k_phase);
	P_k_Divergence.resize(k_phase);
	P_k_Gradient.resize(k_phase);
	A_k_Gradient.resize(k_phase);
	Density_k.resize(k_phase);
	Vel_Divergence.resize(k_phase);
}

void Particle::SetDriftVel(vector3 umk, int k)
{
	if (k > Um_k.size())
		return;
	else
		Um_k[k] = umk;
}

void Particle::SetPosition(vector3 position)
{
	pos = position;
}

void Particle::SetVel(vector3 vel_temp)
{
	vel = vel_temp;
}

void Particle::SetA_k(vector<float> ak)
{
	A_k = ak;//vector中两相体积分数都传进来
}

Fluid::Fluid(int num, int k = 1)
{
	for (int i = 0; i < num; i++)
	{
		Particle *temp = new Particle(k);
		m_ParticleVector.push_back(temp);
	}
	m_num = num;
}

Fluid::~Fluid()
{

}

Particle* Fluid::GetParticleByIndex(int index)
{
	if (m_ParticleVector.size() <= index)
		return NULL;

	return m_ParticleVector[index];
}

void Fluid::DeleteParticleRoomspace()
{
	if (m_ParticleVector.size() > 0)
	{
		std::vector<Particle* >::iterator itor = m_ParticleVector.begin();
		while (itor != m_ParticleVector.end())
		{
			if ((*itor) != NULL)
			{
				delete (*itor);
			}
			itor++;
		}
		m_ParticleVector.clear();
	}
}
