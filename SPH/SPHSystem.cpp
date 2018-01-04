#include "SPHSystem.h"
#include "cpu_sph.h"
#include "MyfileRW.h"
#include <string>
using std::string;
#define PARTICLE(i) SPHSystem::getinstance()->GetFluid()->GetParticleByIndex(i) 
#define GETNUM() SPHSystem::getinstance()->GetFluid()->GetParticleNum()

SPHSystem::SPHSystem()
{
	m_fluid = NULL;
	mat4_set_identity(&mat_col);//�Ѳ���������0    
	m_Calculate = new Calculate;
	m_SingleCal = new SingleCal;
}

void SPHSystem::SetParam(float b, float c, float d, float f)
{
	m_Smoothlen = b;
	m_Timestep = c;
	m_Stiff = d;
	m_R_Search = f;

	if (m_Calculate)
		m_Calculate->setH(m_Smoothlen);

	if (m_SingleCal)
		m_SingleCal->setH(m_Smoothlen);

	grid.SetGrid_len(m_R_Search);
}

void SPHSystem::CreateNighbourList()
{
	n_list.SetSize(m_fluid->GetParticleNum());//�൱�ڳ�ʼ���ڽӱ��������������
	grid.CreateGrid(*m_fluid);//�������д�������
	GridToNeighbour();//�������д����ڽӱ�
}

void SPHSystem::GridToNeighbour()
{
	int gx;
	int gy;
	int gz;
	int gindex;
	int neighbour_grid;
	for (int i = 0; i < n_list.m_sizes; i++)
	{
		//sph_neighbour temp;
		//temp.index = i;
		//temp.distsq = 0.0f;
		//n_list.p[i].push_back(temp);//��ʼ����һ������Ϊ�Լ�����

		gx = (int)((m_fluid->GetParticleByIndex(i)->pos.x - grid.minx) / grid.grid_len);
		gy = (int)((m_fluid->GetParticleByIndex(i)->pos.y - grid.miny) / grid.grid_len);
		gz = (int)((m_fluid->GetParticleByIndex(i)->pos.z - grid.minz) / grid.grid_len);

		gindex = gx + gy * grid.width + gz * grid.width * grid.height;//��ǰ�����������е�����

		for (gz = -1; gz <= 1; gz++)
			for (gy = -1; gy <= 1; gy++)
				for (gx = -1; gx <= 1; gx++)
				{
					neighbour_grid = gindex + gx + gy * grid.width + gz * grid.width * grid.height;//�ڽ����������
					if ((neighbour_grid < 0) || (neighbour_grid >= grid.width*grid.depth*grid.height))
						continue;//������Χ����

					for (int j = 0; j < grid.particles[neighbour_grid].size(); j++)
					{
						float dis;
						int pindex = grid.particles[neighbour_grid][j];
						dis = vec3_distsq(GetParticlePosByIndex(i), GetParticlePosByIndex(pindex));

						if (dis < m_R_Search*m_R_Search)
						{
							sph_neighbour temp1;
							temp1.index = pindex;
							temp1.distsq = dis;
							n_list.p[i].push_back(temp1);//�������Ҫ����ӵ���ǰ�ڽӱ��vector��
						}
					}
				}

		//grid.particles[gindex].push_back(i);//����е����,��ע�͵���ǰ���Ѿ���ӹ���
	}
}

void SPHSystem::CalSphOneStep()
{
	CreateNighbourList();//����ǰ�ȴ����ڽӱ�
	if (m_Calculate)
	{
		m_Calculate->AddFrame();
		m_Calculate->CalPressure(*SPHSystem::getinstance());
		m_Calculate->CalDriftVel(*SPHSystem::getinstance());
		m_Calculate->CalFraction(*SPHSystem::getinstance());
		m_Calculate->CalForce(*SPHSystem::getinstance());
		m_Calculate->CalCollision(*SPHSystem::getinstance());
		m_Calculate->CalPosition(*SPHSystem::getinstance());
	}
}

void SPHSystem::CalSphforsinglemodel()
{
	CreateNighbourList();//����ǰ�ȴ����ڽӱ�
	if (m_SingleCal)
	{
		m_SingleCal->AddFrame();
		m_SingleCal->CalPressure(*SPHSystem::getinstance());
		m_SingleCal->CalForce(*SPHSystem::getinstance());
		m_SingleCal->CalCollision(*SPHSystem::getinstance());
		m_SingleCal->CalPosition(*SPHSystem::getinstance());
	}
}

vector3* SPHSystem::GetParticlePosByIndex(int num)
{
	if (!m_fluid)
		return NULL;
	return &(m_fluid->GetParticleByIndex(num)->pos);
}

SPHSystem::~SPHSystem()
{
	if (m_fluid)
	{
		delete m_fluid;
	}
}

Calculate::Calculate()
{
	poly6_coef = 0.0f;
	grad_poly6_coef = 0.0f;
	lap_poly6_coef = 0.0f;
	grad_spiky_coef = 0.0f;
	lap_vis_coef = 0.0f;
	_tao = (float)pow(0.1, 7);
	sigma = 0.0001f;
	_Frame = 0;
}

void Calculate::CalPressure(SPHSystem &sph)
{
	Fluid *temp = sph.GetFluid();
	if (!temp)
	{
		cout << "CalPressure Fluid is NULL." << endl;
		return;
	}
	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		float density = 0.0f;
		for (int j = 0; j < sph.n_list.p[i].size(); j++)
		{
			int j_Index = sph.n_list.p[i][j].index;
			//�����ܶȹ�ʽ(22)
			density += SPH_function(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index), 1.0);
		}
		PARTICLE(i)->density = density;//���ÿ�����ӵ��ܶ�ƽ��

		//��ԭ�Ĳ����ĵط������ܶȼ��㵥�����ӵĻ��ѹ����24��,û���õ�����ʽ
		temp->GetParticleByIndex(i)->pressure = sph.GetStiff()*(temp->GetParticleByIndex(i)->density - PARTICLE(i)->staticDensity);

		for (int k = 0; k < temp->GetParticleByIndex(i)->k_phase; k++)
		{
			//������������Ϊ�������
			PARTICLE(i)->C_k[k] = PARTICLE(i)->A_k[k];
			//��ԭ�Ĳ����ĵط���ÿһ���ܶ�
			temp->GetParticleByIndex(i)->Density_k[k] = temp->GetParticleByIndex(i)->density;//������趼Ϊˮ����ֻ����ɫ��ͬ����
			//���ѹ������ÿһ���ѹ����ʽ��12��
			if (temp->GetParticleByIndex(i)->P_k.size() > k)
				temp->GetParticleByIndex(i)->P_k[k] = temp->GetParticleByIndex(i)->pressure * temp->GetParticleByIndex(i)->A_k[k];
		}
	}
	//string name = "density" + to_string(_Frame) + ".txt";
	//MyFileRW::getinstance()->SetWriteFileName(name);
	//for (int i = 0; i < temp->GetParticleNum(); i++)
	//{
	//	MyFileRW::getinstance()->WriteSomething("����");
	//	MyFileRW::getinstance()->WriteSomething(i);
	//	MyFileRW::getinstance()->WriteSomething(":", 1);
	//	MyFileRW::getinstance()->WriteSomething(PARTICLE(i)->density, 2);
	//}
	//MyFileRW::getinstance()->CloseWriteFile();
	//�����Ĺ�ʽ��10����ǰ�ڹ���
	vector<vector3> gradient_P_m;
	vector<vector3> gradient_T_m;
	vector<vector3> gradient_T_Dm;
	gradient_P_m.resize(temp->GetParticleNum());
	gradient_T_m.resize(temp->GetParticleNum());
	gradient_T_Dm.resize(temp->GetParticleNum());

	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		for (int j = 0; j < sph.n_list.p[i].size(); j++)
		{
			int j_Index = sph.n_list.p[i][j].index;
			vector3 dist;
			vec3_sub(&dist, &temp->GetParticleByIndex(i)->pos, &temp->GetParticleByIndex(j_Index)->pos);

			//gradient_P_m ���㣬��ʽ��20�����ٸ��ݹ�ʽ��8���������һ����ϵ��ܶȣ���i���ܶȣ�
			gradient_P_m[i] += SPH_function_gradient(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index),
				((temp->GetParticleByIndex(i)->pressure + temp->GetParticleByIndex(j_Index)->pressure)
					/
					(2.0f * temp->GetParticleByIndex(j_Index)->density * PARTICLE(i)->staticDensity))
			);

			//gradient_T_m ����Ƚ�����,����д,��ʽ��21��
			float r = sqrtf(vec3_dot(&dist, &dist));
			vector3 vel_dis;
			vec3_sub(&vel_dis, &temp->GetParticleByIndex(j_Index)->vel, &temp->GetParticleByIndex(i)->vel);//ע��Ϊj_Index-i
			//��ԭ�Ĳ����ĵط�����Ϊ�⻬�˺����ݶ���rj-ri��Լȥ������rj-ri/(rj-ri)*(rj-ri)
			if (i != j_Index)
				vec3_scaleadd(&gradient_T_m[i], &gradient_T_m[i],
					-temp->GetParticleByIndex(j_Index)->mass / temp->GetParticleByIndex(j_Index)->density
					*(temp->GetParticleByIndex(j_Index)->m_Viscosity + temp->GetParticleByIndex(j_Index)->m_Viscosity)*
					grad_spiky_coef*(h - r)*(h - r) / r,/*��һ���Ǻ˺������ݶȵ�������*/
					&vel_dis);

			vector3 gradient_W_ij;
			gradient_W_ij = vec3_scale(grad_spiky_coef*(h - r)*(h - r), &dist);

			for (int k = 0; k < temp->GetParticleByIndex(i)->k_phase; k++)
			{
				//ע������k-���ܶ��Ҽ���ͻ���ܶ�һ���ˣ�����ˮ���ܶȣ��Ժ���������Ҫ�������ʽ
				//��ʽ��19��
				vector3 temp1(0.0f, 0.0f, 0.0f), temp2(0.0f, 0.0f, 0.0f);
				vec3_scaleadd(&temp1, &zeros(), temp->GetParticleByIndex(i)->A_k[k] * vec3_dot(&temp->GetParticleByIndex(i)->Um_k[k], &gradient_W_ij), &temp->GetParticleByIndex(i)->Um_k[k]);
				vec3_scaleadd(&temp2, &zeros(), temp->GetParticleByIndex(j_Index)->A_k[k] * vec3_dot(&temp->GetParticleByIndex(j_Index)->Um_k[k], &gradient_W_ij), &temp->GetParticleByIndex(j_Index)->Um_k[k]);
				vector3 tempsum;
				vec3_add(&tempsum, &temp1, &temp2);
				tempsum *= PARTICLE(i)->staticDensity;//��Ϊ����Ϊ��ˮ�ܶ��ˣ��������Ҫ����Ҫ��Particle����϶���ľ�ˮ�ܶ��������
				vec3_scaleadd(&gradient_T_Dm[i], &gradient_T_Dm[i], -PARTICLE(j_Index)->mass / (PARTICLE(j_Index)->density), &tempsum);//��ԭ�Ĳ����ĵط���k-���ܶ���ȥ

				//��ʽ(13)
				PARTICLE(i)->P_k_Gradient[k] += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index), (PARTICLE(j_Index)->P_k[k] - PARTICLE(i)->P_k[k]) / PARTICLE(j_Index)->density);
				//��ʽ(14)
				PARTICLE(i)->A_k_Gradient[k] += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index), (PARTICLE(j_Index)->A_k[k] - PARTICLE(i)->A_k[k]) / PARTICLE(j_Index)->density);
			}
		}

		PARTICLE(i)->P_m_Gradient = gradient_P_m[i];

		PARTICLE(i)->temp_A = (zeros() - gradient_T_m[i] /*- gradient_T_Dm[i]*/)*(1.0f / PARTICLE(i)->staticDensity);
	}
	//string name = "P_m_Gradient" + to_string(_Frame) + ".txt";
	//MyFileRW::getinstance()->SetWriteFileName(name);//�����P_m_Gradient��ԭ��������һ����
	//for (size_t i = 0; i < GETNUM(); i++)
	//{
	//	MyFileRW::getinstance()->WriteSomething("����");
	//	MyFileRW::getinstance()->WriteSomething(i);
	//	MyFileRW::getinstance()->WriteSomething(" :", 2);
	//	MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->density, 2);
	//	MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->pressure, 2);
	//	MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.x, 1);
	//	MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.y, 1);
	//	MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.z, 2);
	//}
	//MyFileRW::getinstance()->CloseWriteFile();
}

void Calculate::CalDriftVel(SPHSystem &sph)
{
	for (int i = 0; i < GETNUM(); i++)
	{
		float density_m = 0.0f;
		vector3 sum_C_k_P_k;
		//��ʽ��9��
		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			vec3_scaleadd(&sum_C_k_P_k, &sum_C_k_P_k, PARTICLE(i)->C_k[k], &PARTICLE(i)->P_k_Gradient[k]);
			density_m += PARTICLE(i)->C_k[k] * PARTICLE(i)->staticDensity;
		}

		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			//��һ����ʦ�Ӻͦ�ȡ�˶���
			//��ʽ��һ��
			vec3_scaleadd(&PARTICLE(i)->Um_k[k], &PARTICLE(i)->Um_k[k], _tao * (PARTICLE(i)->staticDensity - density_m), &(PARTICLE(i)->temp_A + PARTICLE(i)->P_m_Gradient));
			//��ʽ�ڶ���
			vector3 tempsub;
			vec3_sub(&tempsub, &PARTICLE(i)->P_k_Gradient[k], &sum_C_k_P_k);
			vec3_scaleadd(&PARTICLE(i)->Um_k[k], &PARTICLE(i)->Um_k[k], sigma, &tempsub);
			//������ʽ������
		}
	}
}

void Calculate::CalFraction(SPHSystem &sph)
{
	//��������Ϊ��ʽ��7��
	for (int i = 0; i < GETNUM(); i++)
	{
		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			float tempresult1 = 0.0f;
			float tempresult2 = 0.0f;
			float tempresult3 = 0.0f;
			for (int j = 0; j < sph.n_list.p[i].size(); j++)
			{
				int j_Index = sph.n_list.p[i][j].index;
				//���㹫ʽ��17�����˹�ʽΪ������ʽ���ǵ�����ʽΪ��15��
				vector3 temp = zeros() - PARTICLE(i)->vel;
				tempresult1 += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
					PARTICLE(j_Index)->vel, temp,
					(PARTICLE(i)->A_k[k] + PARTICLE(j_Index)->A_k[k]) / (PARTICLE(j_Index)->density * 2.0));
				//��ʽ��18�����˹�ʽΪ������ʽ���ǵ�����ʽΪ��16��
				tempresult2 += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
					PARTICLE(j_Index)->Um_k[k] * PARTICLE(j_Index)->A_k[k], PARTICLE(i)->Um_k[k] * PARTICLE(i)->A_k[k],
					1.0 / PARTICLE(j_Index)->density);

			}
			//��ʽ��7��
			PARTICLE(i)->A_k[k] += (tempresult2 - tempresult1)*sph.GetTimeStep();//�����������

		}//for k end

		float sum_A_k = 0.0f;
		//ÿһ��������,У���������
		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			if (PARTICLE(i)->A_k[k] < 0)
				PARTICLE(i)->A_k[k] = 0;
			sum_A_k += PARTICLE(i)->A_k[k];
		}

		//���������ѹ��		
		//float tempPressureChange;
		//for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		//{
		//	float tempak = PARTICLE(i)->A_k[k];
		//	PARTICLE(i)->A_k[k] = PARTICLE(i)->A_k[k] / sum_A_k;
		//	tempPressureChange += -sph.GetStiff()*PARTICLE(i)->staticDensity * (PARTICLE(i)->A_k[k] - tempak);
		//}
		////У��ѹ��
		//PARTICLE(i)->pressure += tempPressureChange;
	}

}

void Calculate::CalForce(SPHSystem &sph)
{
	//��������Ϊ��ʽ��8��,��ҪΪ����u_m��ʱ��ĵ���
	vector3 G(0.0f, 0.0f, -9.8f);

	//vector<vector3> tempumum_gredient;
	//tempumum_gredient.resize(GETNUM());
	//������ʽ�ڶ���
	for (size_t i = 0; i < GETNUM(); i++)
	{
		//PARTICLE(i)->P_m_Gradient = zeros();
		//float temp_um_gredient = 0.0f;
		//for (int j = 0; j < sph.n_list.p[i].size(); j++)
		//{
			//int j_Index = sph.n_list.p[i][j].index;
			//temp_um_gredient += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index), PARTICLE(j_Index)->vel, PARTICLE(i)->vel, 0.5);
			//PARTICLE(i)->P_m_Gradient += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
			//	((PARTICLE(i)->pressure + PARTICLE(j_Index)->pressure)
			//		/
			//	(2.0f * PARTICLE(j_Index)->density))
			//);//��Ϊ������ѹ�����¼��㣬��ʽ��20��
		//}
		//tempumum_gredient[i] = PARTICLE(i)->vel * temp_um_gredient;

		//Ϊ�˲�����ֵд����ʱ����
		vector3 test_temp_A(PARTICLE(i)->temp_A);
		vector3 test_P_m_Gradient(PARTICLE(i)->P_m_Gradient);
		//��ʽ��8��
		//PARTICLE(i)->acc = /*G*/zeros() - tempumum_gredient[i] - test_temp_A - PARTICLE(i)->P_m_Gradient*(1.0f/PARTICLE(i)->staticDensity)/*PARTICLE(i)->temp_A*/;
		PARTICLE(i)->acc = G - PARTICLE(i)->temp_A + PARTICLE(i)->P_m_Gradient;
	}
	string name = "FrameForce" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
	for (size_t i = 0; i < GETNUM(); i++)
	{
		MyFileRW::getinstance()->WriteSomething("����");
		MyFileRW::getinstance()->WriteSomething(i);
		MyFileRW::getinstance()->WriteSomething(" :", 2);

		MyFileRW::getinstance()->WriteSomething("�ܶȣ�", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->density, 2);
		MyFileRW::getinstance()->WriteSomething("ѹ����", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->pressure, 2);

		MyFileRW::getinstance()->WriteSomething("���ٶȣ�", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);

		MyFileRW::getinstance()->WriteSomething("ѹ��Ӱ�죺", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.z, 2);

		MyFileRW::getinstance()->WriteSomething("ճ��Ӱ�죺", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->temp_A.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->temp_A.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->temp_A.z, 2);
	}
	MyFileRW::getinstance()->CloseWriteFile();
}
void Calculate::compute_col(vector3* col,
	const vector3* vel, const vector3* n,
	float diff, float stiff, float damp)
{
	float v0;
	float reverse;

	v0 = vec3_dot(n, vel);//���ڻ�    
	reverse = stiff*diff - damp*v0;
	vec3_scaleadd(col, col, reverse, n);
}

void Calculate::glass_collision(vector3* p, vector3* col, const vector3* vel,//ģ��Ͳ�������ײ    
	const matrix4* mat, const matrix4* mat_inv,
	float radius, float stiff, float damp)
{
	float GLASS_R = 0.05f;
	float GLASS_BOTTOM = -0.08f;
	float GLASS_TOP = 10.22f;

	vector3 p_col;
	vector3 n;
	float diff; /** Distance between object and fluid particle **/

	vec3_set(&n, 0.0f, 0.0f, 0.0f);
	mat4_mulvec3(&p_col, mat_inv, p);// �������ӵ�ǰλ�ü���p_col    

	diff = 2.0f*radius - (GLASS_R - (float)sqrt(p_col.x*p_col.x + p_col.y*p_col.y));//����diff    

	if (((diff < 8.0f*radius) && (diff > 0.00001f)) && (p_col.z < GLASS_TOP))
	{
		vec3_set(&n, -p_col.x, -p_col.y, 0.0f);
		vec3_normalize(&n, &n);
		mat4_mulvec3_as_mat3(&n, mat, &n);
		compute_col(col, vel, &n, diff, stiff, damp);//����diff����col    
	}

	diff = 2.0f*radius - (p_col.z - GLASS_BOTTOM);

	if (diff > 0.00001f)
	{
		vec3_set(&n, 0.0f, 0.0f, 1.0f);
		mat4_mulvec3_as_mat3(&n, mat, &n);
		compute_col(col, vel, &n, diff, stiff, damp);
	}
}


void Calculate::CalCollision(SPHSystem &sph)
{
	float sphere_radius = 0.004f;
	float stiff = 30000.0f;
	float damp = 128.0f;
	//MyFileRW::getinstance()->SetWriteFileName("col.txt");
	//string w = "This is NO."+ to_string(_Frame);	
	//MyFileRW::getinstance()->WriteSomething(w,2);
	for (int i = 0; i < GETNUM(); i++)
	{
		vector3 pre_p; /** Ԥ��λ�� **/
		vector3 col;
		vec3_set(&pre_p, 0, 0, 0);
		vec3_set(&col, 0, 0, 0);
		pre_p += PARTICLE(i)->pos + PARTICLE(i)->vel_half * sph.GetTimeStep();
		glass_collision(&pre_p, &col, &PARTICLE(i)->vel, &sph.mat_col, &sph.mat_inv_col, sphere_radius, stiff, damp);
		PARTICLE(i)->acc += col;

		//MyFileRW::getinstance()->WriteSomething("����", 1);
		//MyFileRW::getinstance()->WriteSomething(i, 1);
		//MyFileRW::getinstance()->WriteSomething(":", 2);
		//MyFileRW::getinstance()->WriteFloat(col.x, 1);
		//MyFileRW::getinstance()->WriteFloat(col.y, 1);
		//MyFileRW::getinstance()->WriteFloat(col.z, 2);
	}
	//MyFileRW::getinstance()->CloseWriteFile();
}



void Calculate::CalPosition(SPHSystem &sph)
{
#define PRINTSOME
#ifdef PRINTSOME
	string name = "finalacc" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
#endif
	for (int i = 0; i < GETNUM(); i++)
	{
		vector3 v_half;

		v_half = PARTICLE(i)->vel_half + PARTICLE(i)->acc * sph.GetTimeStep();
		vector3 ShowOneStepDistance = v_half*sph.GetTimeStep();
		PARTICLE(i)->pos += ShowOneStepDistance;
#ifdef PRINTSOME
		MyFileRW::getinstance()->WriteSomething("����");
		MyFileRW::getinstance()->WriteSomething(i, 1);
		MyFileRW::getinstance()->WriteSomething(":",2);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);
		string positioninformation = "λ��:" + to_string((float)PARTICLE(i)->pos.x) + "  " + to_string((float)PARTICLE(i)->pos.y) + "  " + to_string((float)PARTICLE(i)->pos.z);
		MyFileRW::getinstance()->WriteSomething(positioninformation,2);
#endif 
		PARTICLE(i)->vel = (PARTICLE(i)->vel + PARTICLE(i)->vel_half)*0.5;
		PARTICLE(i)->vel_half = v_half;
	}
#ifdef PRINTSOME
	MyFileRW::getinstance()->CloseWriteFile();
#endif 
}

void Calculate::setH(float temph)
{
	h = temph;
	poly6_coef = 315.0f / (64.0f*PI*(float)pow(h, 9));
	grad_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	lap_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	grad_spiky_coef = -45.0f / (PI*h*h*h*h*h*h);
	lap_vis_coef = 45.0f / (PI*h*h*h*h*h*h);
}

vector3 Calculate::SPH_function(Particle i, Particle j, vector3 jfunc, vector3 ifunc /*= zeros()*/)
{
	return zeros();//������ʱ�ò���
}

float Calculate::SPH_function(Particle i, Particle j, float jfunc, float ifunc /*= 0.0f*/)
{
	float h_2 = h * h;
	vector3 temp;
	vec3_set(&temp, 0.0f, 0.0f, 0.0f);
	vec3_sub(&temp, &i.pos, &j.pos);
	//if (temp == zeros())
	//{
	//	return 0.0f;
	//}
	float r_2 = vec3_dot(&temp, &temp);

	float _result = j.mass * poly6_coef * (h_2 - r_2) * (h_2 - r_2) * (h_2 - r_2);
	_result *= (jfunc + ifunc);//ע������û�г���j���ӵ��ܶȣ���jfunc�Լ�ȥ�������û��jfunc������Ϊ1.0

	return _result;
}

vector3 Calculate::SPH_function_gradient(Particle i, Particle j, float jfunc, float ifunc)
{
	vector3 temp;
	vec3_set(&temp, 0.0f, 0.0f, 0.0f);
	vec3_sub(&temp, &i.pos, &j.pos);
	if (temp == zeros())
	{
		return zeros();
	}
	float r = sqrtf(vec3_dot(&temp, &temp));

	vector3 _result;
	vec3_scaleadd(&_result, &_result, -j.mass * grad_spiky_coef * (h - r)*(h - r) / r * (jfunc + ifunc), &temp);

	return _result;
}

float Calculate::SPH_function_gradient(Particle i, Particle j, vector3 jfunc, vector3 ifunc /*= zeros()*/, float Adjustment /*= 1.0f*/)
{
	vector3 dist;
	vec3_set(&dist, 0.0f, 0.0f, 0.0f);
	vec3_sub(&dist, &i.pos, &j.pos);
	if (dist == zeros())
	{
		return 0.0f;
	}
	float _result;
	_result = -j.mass * Adjustment * grad_spiky_coef * vec3_dot(&(jfunc + ifunc), &dist);//Adjustment���ڳ���j���ӵ��ܶȻ������ڵ�����ʽ�õ�
	return _result;
}

vector3 Calculate::SPH_function_gradient_2(Particle i, Particle j, vector3 jfunc, vector3 ifunc /*= zeros()*/, float Adjustment /*= 1.0f*/)
{
	vector3 dist;
	vec3_set(&dist, 0.0f, 0.0f, 0.0f);
	vec3_sub(&dist, &i.pos, &j.pos);
	if (dist == zeros())
	{
		return zeros();
	}
	vector3 _result;
	float r = sqrtf(vec3_dot(&dist, &dist));
	_result = (jfunc + ifunc) * j.mass * lap_vis_coef * (h - r) * Adjustment;
	return _result;
}

void SingleCal::CalPressure(SPHSystem &sph)
{
	Fluid *temp = sph.GetFluid();
	if (!temp)
	{
		cout << "CalPressure Fluid is NULL." << endl;
		return;
	}
	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		float density = 0.0f;
		for (int j = 0; j < sph.n_list.p[i].size(); j++)
		{
			int j_Index = sph.n_list.p[i][j].index;
			//�����ܶȹ�ʽ(22)
			density += SPH_function(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index), 1.0);
		}
		PARTICLE(i)->density = density;//���ÿ�����ӵ��ܶ�ƽ��

		temp->GetParticleByIndex(i)->pressure = sph.GetStiff()*(temp->GetParticleByIndex(i)->density - PARTICLE(i)->staticDensity);
	}
#define PIRNTTEST
#ifdef PIRNTTEST
	string name = "density" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		MyFileRW::getinstance()->WriteSomething("����");
		MyFileRW::getinstance()->WriteSomething(i);
		MyFileRW::getinstance()->WriteSomething(":", 1);
		MyFileRW::getinstance()->WriteSomething(PARTICLE(i)->density, 2);
	}
	MyFileRW::getinstance()->CloseWriteFile();
#endif
}

void SingleCal::CalForce(SPHSystem &sph)
{
	vector3 G(0.0f, 0.0f, -9.8f);
	Fluid *temp = sph.GetFluid();
#ifdef PIRNTTEST
	string name = "singleacc" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
#endif
	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		vector3 pressureforce;
		vector3 viscosityforce;
		for (int j = 0; j < sph.n_list.p[i].size(); j++)
		{
			int j_Index = sph.n_list.p[i][j].index;
			//����ѹ��
			pressureforce += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
				(PARTICLE(i)->pressure + PARTICLE(j_Index)->pressure) / (2.0f*PARTICLE(j_Index)->density * PARTICLE(i)->density)
			);

			vector3 vel_diff = PARTICLE(j_Index)->vel - PARTICLE(i)->vel;
			viscosityforce += SPH_function_gradient_2(*PARTICLE(i), *PARTICLE(j_Index), vel_diff, zeros(),
				PARTICLE(j_Index)->m_Viscosity / (PARTICLE(j_Index)->density*PARTICLE(i)->density));
		}
		PARTICLE(i)->acc = G + pressureforce + viscosityforce;
#ifdef PIRNTTEST
		string something = "����" + to_string(i) + ":";
		MyFileRW::getinstance()->WriteSomething(something, 2);
		something = "ѹ�����ٶȣ�" + to_string((float)pressureforce.x)+"  "+to_string((float)pressureforce.y)+"  "+to_string((float)pressureforce.z);
		MyFileRW::getinstance()->WriteSomething(something, 2);
		something = "ճ�ȼ��ٶȣ�" + to_string((float)viscosityforce.x) + "  " + to_string((float)viscosityforce.y) + "  " + to_string((float)viscosityforce.z);;
		MyFileRW::getinstance()->WriteSomething(something, 2);
#endif
	}
#ifdef PIRNTTEST
	MyFileRW::getinstance()->CloseWriteFile();
#endif
}
