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
	mat4_set_identity(&mat_col);//把参数内容置0    
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
	n_list.SetSize(m_fluid->GetParticleNum());//相当于初始化邻接表，里面有清除操作
	grid.CreateGrid(*m_fluid);//从流体中创建网格
	GridToNeighbour();//从网格中创建邻接表
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
		//n_list.p[i].push_back(temp);//初始化第一个粒子为自己本身

		gx = (int)((m_fluid->GetParticleByIndex(i)->pos.x - grid.minx) / grid.grid_len);
		gy = (int)((m_fluid->GetParticleByIndex(i)->pos.y - grid.miny) / grid.grid_len);
		gz = (int)((m_fluid->GetParticleByIndex(i)->pos.z - grid.minz) / grid.grid_len);

		gindex = gx + gy * grid.width + gz * grid.width * grid.height;//当前粒子在网格中的索引

		for (gz = -1; gz <= 1; gz++)
			for (gy = -1; gy <= 1; gy++)
				for (gx = -1; gx <= 1; gx++)
				{
					neighbour_grid = gindex + gx + gy * grid.width + gz * grid.width * grid.height;//邻接网格的索引
					if ((neighbour_grid < 0) || (neighbour_grid >= grid.width*grid.depth*grid.height))
						continue;//超出范围跳过

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
							n_list.p[i].push_back(temp1);//满足距离要求添加到当前邻接表的vector中
						}
					}
				}

		//grid.particles[gindex].push_back(i);//这句有点奇怪,先注释掉，前面已经添加过了
	}
}

void SPHSystem::CalSphOneStep()
{
	CreateNighbourList();//计算前先创建邻接表
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
	CreateNighbourList();//计算前先创建邻接表
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
			//计算密度公式(22)
			density += SPH_function(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index), 1.0);
		}
		PARTICLE(i)->density = density;//获得每个粒子的密度平均

		//与原文不符的地方：从密度计算单个粒子的混合压力（24）,没有用调整公式
		temp->GetParticleByIndex(i)->pressure = sph.GetStiff()*(temp->GetParticleByIndex(i)->density - PARTICLE(i)->staticDensity);

		for (int k = 0; k < temp->GetParticleByIndex(i)->k_phase; k++)
		{
			//质量分数更新为体积分数
			PARTICLE(i)->C_k[k] = PARTICLE(i)->A_k[k];
			//与原文不符的地方：每一项密度
			temp->GetParticleByIndex(i)->Density_k[k] = temp->GetParticleByIndex(i)->density;//这里假设都为水分子只是颜色不同而已
			//混合压力计算每一相的压力公式（12）
			if (temp->GetParticleByIndex(i)->P_k.size() > k)
				temp->GetParticleByIndex(i)->P_k[k] = temp->GetParticleByIndex(i)->pressure * temp->GetParticleByIndex(i)->A_k[k];
		}
	}
	//string name = "density" + to_string(_Frame) + ".txt";
	//MyFileRW::getinstance()->SetWriteFileName(name);
	//for (int i = 0; i < temp->GetParticleNum(); i++)
	//{
	//	MyFileRW::getinstance()->WriteSomething("粒子");
	//	MyFileRW::getinstance()->WriteSomething(i);
	//	MyFileRW::getinstance()->WriteSomething(":", 1);
	//	MyFileRW::getinstance()->WriteSomething(PARTICLE(i)->density, 2);
	//}
	//MyFileRW::getinstance()->CloseWriteFile();
	//求论文公式（10）的前期工作
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

			//gradient_P_m 计算，公式（20），再根据公式（8），多除了一个混合的密度（即i的密度）
			gradient_P_m[i] += SPH_function_gradient(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index),
				((temp->GetParticleByIndex(i)->pressure + temp->GetParticleByIndex(j_Index)->pressure)
					/
					(2.0f * temp->GetParticleByIndex(j_Index)->density * PARTICLE(i)->staticDensity))
			);

			//gradient_T_m 计算比较特殊,单独写,公式（21）
			float r = sqrtf(vec3_dot(&dist, &dist));
			vector3 vel_dis;
			vec3_sub(&vel_dis, &temp->GetParticleByIndex(j_Index)->vel, &temp->GetParticleByIndex(i)->vel);//注意为j_Index-i
			//与原文不符的地方：因为光滑核函数梯度中rj-ri，约去了最后的rj-ri/(rj-ri)*(rj-ri)
			if (i != j_Index)
				vec3_scaleadd(&gradient_T_m[i], &gradient_T_m[i],
					-temp->GetParticleByIndex(j_Index)->mass / temp->GetParticleByIndex(j_Index)->density
					*(temp->GetParticleByIndex(j_Index)->m_Viscosity + temp->GetParticleByIndex(j_Index)->m_Viscosity)*
					grad_spiky_coef*(h - r)*(h - r) / r,/*这一行是核函数的梯度的所带的*/
					&vel_dis);

			vector3 gradient_W_ij;
			gradient_W_ij = vec3_scale(grad_spiky_coef*(h - r)*(h - r), &dist);

			for (int k = 0; k < temp->GetParticleByIndex(i)->k_phase; k++)
			{
				//注意这里k-相密度我假设和混合密度一样了，都是水的密度，以后有问题需要改这个公式
				//公式（19）
				vector3 temp1(0.0f, 0.0f, 0.0f), temp2(0.0f, 0.0f, 0.0f);
				vec3_scaleadd(&temp1, &zeros(), temp->GetParticleByIndex(i)->A_k[k] * vec3_dot(&temp->GetParticleByIndex(i)->Um_k[k], &gradient_W_ij), &temp->GetParticleByIndex(i)->Um_k[k]);
				vec3_scaleadd(&temp2, &zeros(), temp->GetParticleByIndex(j_Index)->A_k[k] * vec3_dot(&temp->GetParticleByIndex(j_Index)->Um_k[k], &gradient_W_ij), &temp->GetParticleByIndex(j_Index)->Um_k[k]);
				vector3 tempsum;
				vec3_add(&tempsum, &temp1, &temp2);
				tempsum *= PARTICLE(i)->staticDensity;//因为都设为静水密度了，如果这样要用又要在Particle里加上多相的静水密度这个参数
				vec3_scaleadd(&gradient_T_Dm[i], &gradient_T_Dm[i], -PARTICLE(j_Index)->mass / (PARTICLE(j_Index)->density), &tempsum);//与原文不符的地方：k-相密度消去

				//公式(13)
				PARTICLE(i)->P_k_Gradient[k] += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index), (PARTICLE(j_Index)->P_k[k] - PARTICLE(i)->P_k[k]) / PARTICLE(j_Index)->density);
				//公式(14)
				PARTICLE(i)->A_k_Gradient[k] += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index), (PARTICLE(j_Index)->A_k[k] - PARTICLE(i)->A_k[k]) / PARTICLE(j_Index)->density);
			}
		}

		PARTICLE(i)->P_m_Gradient = gradient_P_m[i];

		PARTICLE(i)->temp_A = (zeros() - gradient_T_m[i] /*- gradient_T_Dm[i]*/)*(1.0f / PARTICLE(i)->staticDensity);
	}
	//string name = "P_m_Gradient" + to_string(_Frame) + ".txt";
	//MyFileRW::getinstance()->SetWriteFileName(name);//这里的P_m_Gradient与原来给的是一样的
	//for (size_t i = 0; i < GETNUM(); i++)
	//{
	//	MyFileRW::getinstance()->WriteSomething("粒子");
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
		//公式（9）
		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			vec3_scaleadd(&sum_C_k_P_k, &sum_C_k_P_k, PARTICLE(i)->C_k[k], &PARTICLE(i)->P_k_Gradient[k]);
			density_m += PARTICLE(i)->C_k[k] * PARTICLE(i)->staticDensity;
		}

		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			//问一下老师τ和σ取了多少
			//右式第一项
			vec3_scaleadd(&PARTICLE(i)->Um_k[k], &PARTICLE(i)->Um_k[k], _tao * (PARTICLE(i)->staticDensity - density_m), &(PARTICLE(i)->temp_A + PARTICLE(i)->P_m_Gradient));
			//右式第二项
			vector3 tempsub;
			vec3_sub(&tempsub, &PARTICLE(i)->P_k_Gradient[k], &sum_C_k_P_k);
			vec3_scaleadd(&PARTICLE(i)->Um_k[k], &PARTICLE(i)->Um_k[k], sigma, &tempsub);
			//暂无右式第三项
		}
	}
}

void Calculate::CalFraction(SPHSystem &sph)
{
	//整个方法为公式（7）
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
				//先算公式（17），此公式为调整公式，非调整公式为（15）
				vector3 temp = zeros() - PARTICLE(i)->vel;
				tempresult1 += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
					PARTICLE(j_Index)->vel, temp,
					(PARTICLE(i)->A_k[k] + PARTICLE(j_Index)->A_k[k]) / (PARTICLE(j_Index)->density * 2.0));
				//公式（18），此公式为调整公式，非调整公式为（16）
				tempresult2 += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
					PARTICLE(j_Index)->Um_k[k] * PARTICLE(j_Index)->A_k[k], PARTICLE(i)->Um_k[k] * PARTICLE(i)->A_k[k],
					1.0 / PARTICLE(j_Index)->density);

			}
			//公式（7）
			PARTICLE(i)->A_k[k] += (tempresult2 - tempresult1)*sph.GetTimeStep();//更新体积分数

		}//for k end

		float sum_A_k = 0.0f;
		//每一相计算完后,校正体积分数
		for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		{
			if (PARTICLE(i)->A_k[k] < 0)
				PARTICLE(i)->A_k[k] = 0;
			sum_A_k += PARTICLE(i)->A_k[k];
		}

		//如果不矫正压力		
		//float tempPressureChange;
		//for (int k = 0; k < PARTICLE(i)->k_phase; k++)
		//{
		//	float tempak = PARTICLE(i)->A_k[k];
		//	PARTICLE(i)->A_k[k] = PARTICLE(i)->A_k[k] / sum_A_k;
		//	tempPressureChange += -sph.GetStiff()*PARTICLE(i)->staticDensity * (PARTICLE(i)->A_k[k] - tempak);
		//}
		////校正压力
		//PARTICLE(i)->pressure += tempPressureChange;
	}

}

void Calculate::CalForce(SPHSystem &sph)
{
	//整个方法为公式（8）,主要为计算u_m随时间的导数
	vector3 G(0.0f, 0.0f, -9.8f);

	//vector<vector3> tempumum_gredient;
	//tempumum_gredient.resize(GETNUM());
	//计算左式第二项
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
			//);//因为更新了压力重新计算，公式（20）
		//}
		//tempumum_gredient[i] = PARTICLE(i)->vel * temp_um_gredient;

		//为了测试数值写的临时变量
		vector3 test_temp_A(PARTICLE(i)->temp_A);
		vector3 test_P_m_Gradient(PARTICLE(i)->P_m_Gradient);
		//公式（8）
		//PARTICLE(i)->acc = /*G*/zeros() - tempumum_gredient[i] - test_temp_A - PARTICLE(i)->P_m_Gradient*(1.0f/PARTICLE(i)->staticDensity)/*PARTICLE(i)->temp_A*/;
		PARTICLE(i)->acc = G - PARTICLE(i)->temp_A + PARTICLE(i)->P_m_Gradient;
	}
	string name = "FrameForce" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
	for (size_t i = 0; i < GETNUM(); i++)
	{
		MyFileRW::getinstance()->WriteSomething("粒子");
		MyFileRW::getinstance()->WriteSomething(i);
		MyFileRW::getinstance()->WriteSomething(" :", 2);

		MyFileRW::getinstance()->WriteSomething("密度：", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->density, 2);
		MyFileRW::getinstance()->WriteSomething("压力：", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->pressure, 2);

		MyFileRW::getinstance()->WriteSomething("加速度：", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);

		MyFileRW::getinstance()->WriteSomething("压力影响：", 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->P_m_Gradient.z, 2);

		MyFileRW::getinstance()->WriteSomething("粘度影响：", 1);
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

	v0 = vec3_dot(n, vel);//求内积    
	reverse = stiff*diff - damp*v0;
	vec3_scaleadd(col, col, reverse, n);
}

void Calculate::glass_collision(vector3* p, vector3* col, const vector3* vel,//模拟和玻璃杯碰撞    
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
	mat4_mulvec3(&p_col, mat_inv, p);// 根据粒子当前位置计算p_col    

	diff = 2.0f*radius - (GLASS_R - (float)sqrt(p_col.x*p_col.x + p_col.y*p_col.y));//计算diff    

	if (((diff < 8.0f*radius) && (diff > 0.00001f)) && (p_col.z < GLASS_TOP))
	{
		vec3_set(&n, -p_col.x, -p_col.y, 0.0f);
		vec3_normalize(&n, &n);
		mat4_mulvec3_as_mat3(&n, mat, &n);
		compute_col(col, vel, &n, diff, stiff, damp);//根据diff计算col    
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
		vector3 pre_p; /** 预测位置 **/
		vector3 col;
		vec3_set(&pre_p, 0, 0, 0);
		vec3_set(&col, 0, 0, 0);
		pre_p += PARTICLE(i)->pos + PARTICLE(i)->vel_half * sph.GetTimeStep();
		glass_collision(&pre_p, &col, &PARTICLE(i)->vel, &sph.mat_col, &sph.mat_inv_col, sphere_radius, stiff, damp);
		PARTICLE(i)->acc += col;

		//MyFileRW::getinstance()->WriteSomething("粒子", 1);
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
		MyFileRW::getinstance()->WriteSomething("粒子");
		MyFileRW::getinstance()->WriteSomething(i, 1);
		MyFileRW::getinstance()->WriteSomething(":",2);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);
		string positioninformation = "位置:" + to_string((float)PARTICLE(i)->pos.x) + "  " + to_string((float)PARTICLE(i)->pos.y) + "  " + to_string((float)PARTICLE(i)->pos.z);
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
	return zeros();//好像暂时用不到
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
	_result *= (jfunc + ifunc);//注意这里没有除以j粒子的密度，在jfunc自己去除，如果没有jfunc可以设为1.0

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
	_result = -j.mass * Adjustment * grad_spiky_coef * vec3_dot(&(jfunc + ifunc), &dist);//Adjustment用于除以j粒子的密度或者用于调整公式用的
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
			//计算密度公式(22)
			density += SPH_function(*temp->GetParticleByIndex(i), *temp->GetParticleByIndex(j_Index), 1.0);
		}
		PARTICLE(i)->density = density;//获得每个粒子的密度平均

		temp->GetParticleByIndex(i)->pressure = sph.GetStiff()*(temp->GetParticleByIndex(i)->density - PARTICLE(i)->staticDensity);
	}
#define PIRNTTEST
#ifdef PIRNTTEST
	string name = "density" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
	for (int i = 0; i < temp->GetParticleNum(); i++)
	{
		MyFileRW::getinstance()->WriteSomething("粒子");
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
			//计算压力
			pressureforce += SPH_function_gradient(*PARTICLE(i), *PARTICLE(j_Index),
				(PARTICLE(i)->pressure + PARTICLE(j_Index)->pressure) / (2.0f*PARTICLE(j_Index)->density * PARTICLE(i)->density)
			);

			vector3 vel_diff = PARTICLE(j_Index)->vel - PARTICLE(i)->vel;
			viscosityforce += SPH_function_gradient_2(*PARTICLE(i), *PARTICLE(j_Index), vel_diff, zeros(),
				PARTICLE(j_Index)->m_Viscosity / (PARTICLE(j_Index)->density*PARTICLE(i)->density));
		}
		PARTICLE(i)->acc = G + pressureforce + viscosityforce;
#ifdef PIRNTTEST
		string something = "粒子" + to_string(i) + ":";
		MyFileRW::getinstance()->WriteSomething(something, 2);
		something = "压力加速度：" + to_string((float)pressureforce.x)+"  "+to_string((float)pressureforce.y)+"  "+to_string((float)pressureforce.z);
		MyFileRW::getinstance()->WriteSomething(something, 2);
		something = "粘度加速度：" + to_string((float)viscosityforce.x) + "  " + to_string((float)viscosityforce.y) + "  " + to_string((float)viscosityforce.z);;
		MyFileRW::getinstance()->WriteSomething(something, 2);
#endif
	}
#ifdef PIRNTTEST
	MyFileRW::getinstance()->CloseWriteFile();
#endif
}
