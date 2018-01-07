#include "SPHSystem.h"
#include "cpu_sph.h"
#include "MyfileRW.h"
#include <string>
using std::string;
#define PARTICLE(i) SPHNewSystem::getinstance()->GetFluid()->GetParticleByIndex(i) 
#define GETNUM() SPHNewSystem::getinstance()->GetFluid()->GetParticleNum()

void SPHNewSystem::init()
{
	int i;
	int x;
	int y;
	int z;
	float s;
	float s2;
	float cx;
	float cy;
	float cz;

	s = 0.006f;
	cx = cy = 0.0f;
	cz = 0.05f;

	s2 = 0.001f;
	//初始化流体,流体中N_PARTICLES个粒子，粒子有NUM_PHASE相

	//保留功能的时候，这里记得释放空间，暂时还没有释放空间
	Fluid *input = new Fluid(N_PARTICLES, NUM_PHASE);

#ifdef FILEINPUT
	//test io
	MyFileRW::getinstance()->SetReadNmae("./Data/input.txt");//SetName和close要成对出现
	MyFileRW::getinstance()->ReadFile(input);
	MyFileRW::getinstance()->CloseReadFile();
#else
	for (x = 0; x < CUBE_LEN_X; x++)
	{
		for (y = 0; y < CUBE_LEN_Y; y++)
		{
			for (z = 0; z < LEN_Z; z++)
			{
				i = x + y*CUBE_LEN_X + z*CUBE_LEN_X*CUBE_LEN_Y;
				if (input->m_ParticleVector.size() < i)
					continue;
				//设置位置
				//vec3_set(&(input->GetParticleByIndex(i)->pos), s*(x - CUBE_LEN_X / 2) - cx, s*(y - CUBE_LEN_Y / 2) - cy, 0.8*s*z - cz + 0.1f);
				vec3_set(&(input->GetParticleByIndex(i)->pos), s*(x - CUBE_LEN_X / 2) - cx, s*(y - CUBE_LEN_Y / 2) - cy, s*z - cz);
				//设置速度
				vec3_set(&(input->GetParticleByIndex(i)->vel), 0.0f, 0.0f, 0.0f);
				input->GetParticleByIndex(i)->vel_half = input->GetParticleByIndex(i)->vel;
				input->GetParticleByIndex(i)->SetstaticDensity(STATICDENSITY);
				//以x轴为分界，x取正为红色，x取负为蓝色
				if (x >= CUBE_LEN_X / 2)
				{
					//设置体积分数
					input->GetParticleByIndex(i)->A_k[0] = (1.0f);//设立体积分数，第一项为1
					input->GetParticleByIndex(i)->A_k[1] = (0.0f);//第二项为0
																  //设置质量分数
					input->GetParticleByIndex(i)->C_k[0] = (1.0f);//设立体积分数，第一项为1
					input->GetParticleByIndex(i)->C_k[1] = (0.0f);//第二项为0
				}
				else
				{
					input->GetParticleByIndex(i)->A_k[0] = (0.0f);//设立体积分数，第一项为0
					input->GetParticleByIndex(i)->A_k[1] = (1.0f);//第二项为1

					input->GetParticleByIndex(i)->C_k[0] = (0.0f);//设立体积分数，第一项为0
					input->GetParticleByIndex(i)->C_k[1] = (1.0f);//第二项为1
				}
				//设置每个粒子的质量
				input->GetParticleByIndex(i)->SetMass(MASS);
				//设置粘度
				input->GetParticleByIndex(i)->SetViscosity(VISCOSITY);
			}
		}
	}
#endif
	SPHNewSystem::getinstance()->SetFluid(input);
	SPHNewSystem::getinstance()->SetParam(SMOOTHING_LENGTH, TIMESTEP, STIFF, SEARCH_RADIUS);
}

SPHNewSystem::SPHNewSystem()
{
	m_fluid = NULL;
	mat4_set_identity(&mat_col);//把参数内容置0    
	m_Calculate = new Calculate_Reverce;
}

void SPHNewSystem::SetParam(float b, float c, float d, float f)
{
	m_Smoothlen = b;
	m_Timestep = c;
	m_Stiff = d;
	m_R_Search = f;

	if (m_Calculate)
		m_Calculate->setH(m_Smoothlen);

	grid.SetGrid_len(m_R_Search);
}

void SPHNewSystem::CreateNighbourList()
{
	n_list.SetSize(m_fluid->GetParticleNum());//相当于初始化邻接表，里面有清除操作
	grid.CreateGrid(*m_fluid);//从流体中创建网格
	GridToNeighbour();//从网格中创建邻接表
}

void SPHNewSystem::GridToNeighbour()
{
	int gx;
	int gy;
	int gz;
	int gindex;
	int neighbour_grid;
	for (int i = 0; i < n_list.m_sizes; i++)
	{
		sph_neighbour temp;
		temp.index = i;
		temp.distsq = 0.0f;
		n_list.p[i].clear();//使用前先清空
		n_list.p[i].push_back(temp);//初始化第一个粒子为自己本身

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
						if(pindex==i)
							continue;
						if (dis < m_R_Search*m_R_Search)
						{
							sph_neighbour temp1;
							temp1.index = pindex;
							temp1.distsq = dis;
							n_list.p[i].push_back(temp1);//满足距离要求添加到当前邻接表的vector中
						}
					}
				}
		grid.particles[gindex].push_back(i);
	}
}

void SPHNewSystem::CalSphOneStep()
{
	//重新按照原来的循环方式和网格构造模式
	for (int i = 0;i < N_STEPS;i++)
	{
		if (i == 0)
			CreateNighbourList();//计算前先创建邻接表

		if (m_Calculate)
		{
			m_Calculate->AddFrame();

			m_Calculate->CalPressure(*SPHNewSystem::getinstance(), i);
			m_Calculate->CalDriftVel(*SPHNewSystem::getinstance());
			m_Calculate->CalFraction(*SPHNewSystem::getinstance());
			m_Calculate->CalForce(*SPHNewSystem::getinstance());
			m_Calculate->CalCollision(*SPHNewSystem::getinstance());
			m_Calculate->CalPosition(*SPHNewSystem::getinstance());
		}
	}
}

vector3* SPHNewSystem::GetParticlePosByIndex(int num)
{
	if (!m_fluid)
		return NULL;
	return &(m_fluid->GetParticleByIndex(num)->pos);
}

SPHNewSystem::~SPHNewSystem()
{
	if (m_fluid)
	{
		delete m_fluid;
	}
}

Calculate_Reverce::Calculate_Reverce()
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
//今天的目标单相模型重写一遍
void Calculate_Reverce::CalPressure(SPHNewSystem &sph, int whattime)
{
	//清空密度重新计算
	for (int i = 0;i < GETNUM();i++)
	{
		PARTICLE(i)->density = 0;
	}
	float h2 = sph.GetH2();
	if (whattime == 0)//表示第一次
	{
		for (int i = 0;i < GETNUM();i++)
		{
			for (int j = 0;j < sph.n_list.p[i].size();j++)
			{
				float dist = 0.0f;
				int Jindex = sph.n_list.p[i][j].index;//邻接表的粒子索引
				dist = sph.n_list.p[i][j].distsq;

				if (h2 > dist)
				{
					//为了更好的看出公式的错误，尽量公式完整表示
					//公式(22)
					PARTICLE(i)->density += SPH_function(*PARTICLE(i), *PARTICLE(Jindex), 1.0f);

					if (Jindex != i)
						PARTICLE(Jindex)->density += SPH_function(*PARTICLE(Jindex), *PARTICLE(i), 1.0f);
				}
			}
		}
	}
	else
	{
		for (int i = 0;i < GETNUM();i++)
		{
			for (int j = 0;j < sph.n_list.p[i].size();j++)
			{
				float dist = 0.0f;
				int Jindex = sph.n_list.p[i][j].index;//邻接表的粒子索引
				//因为进行了一次计算，我们假设邻居还是原来的邻居，但是位置发生了一些改变，所以更新粒子的距离
				dist = vec3_distsq(&PARTICLE(i)->pos, &PARTICLE(Jindex)->pos);
				sph.n_list.p[i][j].distsq = dist;

				if (h2 > dist)
				{
					//为了更好的看出公式的错误，尽量公式完整表示
					//公式(22)
					PARTICLE(i)->density += SPH_function(*PARTICLE(i), *PARTICLE(Jindex), 1.0f);

					if (Jindex != i)
						PARTICLE(Jindex)->density += SPH_function(*PARTICLE(Jindex), *PARTICLE(i), 1.0f);
				}
			}
		}
	}

//#define TESTPRINT
#ifdef TESTPRINT
	string name = "./Data/yali" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
#endif	
	//计算压力
	for (size_t i = 0; i < GETNUM(); i++)
	{
		PARTICLE(i)->pressure = sph.GetStiff()*(PARTICLE(i)->density - PARTICLE(i)->staticDensity);//减去静态密度（1000）

#ifdef TESTPRINT
		string some;
		some = "粒子" + to_string(i) + ":";
		MyFileRW::getinstance()->WriteSomething(some, 2);
		some = "密度:" + to_string(PARTICLE(i)->density);
		MyFileRW::getinstance()->WriteSomething(some, 2);
		some = "压力:" + to_string(PARTICLE(i)->pressure);
		MyFileRW::getinstance()->WriteSomething(some, 2);
#endif	
	}
#ifdef TESTPRINT
	MyFileRW::getinstance()->CloseWriteFile();
#endif

}

void Calculate_Reverce::CalDriftVel(SPHNewSystem &sph)
{
}

void Calculate_Reverce::CalFraction(SPHNewSystem &sph)
{

}

void Calculate_Reverce::CalForce(SPHNewSystem &sph)
{
	//清空加速度，重新计算
	for (int i = 0;i < GETNUM();i++)
	{
		PARTICLE(i)->acc = zeros();
	}
	float h = sph.GetH();
	for (int i = 0;i < GETNUM();i++)
	{
		for (int j = 1;j < sph.n_list.p[i].size();j++)
		{
			int Jindex = sph.n_list.p[i][j].index;//邻接表的粒子索引
			vector3 posdiff = PARTICLE(i)->pos - PARTICLE(Jindex)->pos;
			float r = sqrtf(vec3_dot(&posdiff, &posdiff));
			//压力影响

			//公式（20），并根据公式（8）多除以混合密度（对于单相，混合密度为 粒子i的密度）
			vector3 tempPacc;
			tempPacc = posdiff / r * (PARTICLE(i)->pressure + PARTICLE(Jindex)->pressure) / (2.0f*PARTICLE(i)->density*PARTICLE(Jindex)->density) * -grad_spiky_coef* (h - r)*(h - r);

			//一次计算
			PARTICLE(i)->acc += tempPacc * PARTICLE(Jindex)->mass;
			PARTICLE(Jindex)->acc -= tempPacc * PARTICLE(i)->mass;


			//粘度影响
			vector3 veldiff = PARTICLE(Jindex)->vel - PARTICLE(i)->vel;//速度差
			vector3 tempVisacc;
#define SINGCAL
#ifdef SINGCAL
			//使用2阶梯度
			tempVisacc = veldiff / (PARTICLE(i)->density * PARTICLE(Jindex)->density)* PARTICLE(Jindex)->m_Viscosity  * lap_vis_coef*(h - r);
#else
			//公式（21），同理除以混合密度
			tempVisacc = veldiff * ((PARTICLE(Jindex)->m_Viscosity + PARTICLE(i)->m_Viscosity) / (PARTICLE(i)->density*PARTICLE(Jindex)->density)) *  -grad_spiky_coef* (h - r)*(h - r) / r;
#endif
			//一次计算
			PARTICLE(i)->acc += tempVisacc * PARTICLE(Jindex)->mass;
			PARTICLE(Jindex)->acc -= tempVisacc * PARTICLE(i)->mass;
		}
	}
//#define TESTPRINT
#ifdef TESTPRINT
	string name = "./Data/Force" + to_string(_Frame) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
	for (int i = 0;i < GETNUM();i++)
	{
		string some;
		some = "粒子" + to_string(i) + ":";
		MyFileRW::getinstance()->WriteSomething(some, 2);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 2);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);
	}
	MyFileRW::getinstance()->CloseWriteFile();
#endif
}
void Calculate_Reverce::compute_col(vector3* col,
	const vector3* vel, const vector3* n,
	float diff, float stiff, float damp)
{
	float v0;
	float reverse;

	v0 = vec3_dot(n, vel);//求内积    
	reverse = stiff*diff - damp*v0;
	vec3_scaleadd(col, col, reverse, n);
}

void Calculate_Reverce::glass_collision(vector3* p, vector3* col, const vector3* vel,//模拟和玻璃杯碰撞    
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


void Calculate_Reverce::CalCollision(SPHNewSystem &sph)
{
	float sphere_radius = 0.004f;
	float stiff = 30000.0f;
	float damp = 128.0f;
//#define PRINTCOLL
#ifdef PRINTCOLL
	string name = "./Data/collision" + to_string(_Frame - 1) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
#endif	
	for (int i = 0; i < GETNUM(); i++)
	{
		vector3 pre_p; /** 预测位置 **/
		vector3 col;

		pre_p = PARTICLE(i)->pos + PARTICLE(i)->vel_half * sph.GetTimeStep();
		glass_collision(&pre_p, &col, &PARTICLE(i)->vel, &sph.mat_col, &sph.mat_inv_col, sphere_radius, stiff, damp);
		PARTICLE(i)->acc += col;
#ifdef PRINTCOLL
		string some = "粒子" + to_string(i)+":";
		MyFileRW::getinstance()->WriteSomething(some, 2);
		some = "碰撞引起的加速度:"+to_string(col.x)+" " + to_string(col.y)+" " + to_string(col.z);
		MyFileRW::getinstance()->WriteSomething(some, 2);
#endif
	}

#ifdef PRINTCOLL
	MyFileRW::getinstance()->CloseWriteFile();
#endif	
}



void Calculate_Reverce::CalPosition(SPHNewSystem &sph)
{
//#define PRINTSOME
#ifdef PRINTSOME
	string name = "./Data/finalacc" + to_string(_Frame-1) + ".txt";
	MyFileRW::getinstance()->SetWriteFileName(name);
#endif
	vector3 G(0.0f, 0.0f, -9.8f);
	for (int i = 0; i < GETNUM(); i++)
	{
		vector3 v_half;
		PARTICLE(i)->acc += G;
		v_half = PARTICLE(i)->vel_half + PARTICLE(i)->acc * sph.GetTimeStep();
		vector3 ShowOneStepDistance = v_half*sph.GetTimeStep();
		PARTICLE(i)->pos += ShowOneStepDistance;
#ifdef PRINTSOME
		MyFileRW::getinstance()->WriteSomething("粒子");
		MyFileRW::getinstance()->WriteSomething(i, 1);
		MyFileRW::getinstance()->WriteSomething(":", 2);
		MyFileRW::getinstance()->WriteSomething("加速度 ： ");
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.x, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.y, 1);
		MyFileRW::getinstance()->WriteFloat(PARTICLE(i)->acc.z, 2);
		string positioninformation = "位置:" + to_string((float)PARTICLE(i)->pos.x) + "  " + to_string((float)PARTICLE(i)->pos.y) + "  " + to_string((float)PARTICLE(i)->pos.z);
		MyFileRW::getinstance()->WriteSomething(positioninformation, 2);
#endif 
		PARTICLE(i)->vel = (PARTICLE(i)->vel + v_half)*0.5;//有点问题
		PARTICLE(i)->vel_half = v_half;
	}
#ifdef PRINTSOME
	MyFileRW::getinstance()->CloseWriteFile();
#endif 
}

void Calculate_Reverce::setH(float temph)
{
	h = temph;
	poly6_coef = 315.0f / (64.0f*PI*(float)pow(h, 9));
	grad_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	lap_poly6_coef = 945.0f / (32.0f*PI*(float)pow(h, 9));
	grad_spiky_coef = -45.0f / (PI*h*h*h*h*h*h);
	lap_vis_coef = 45.0f / (PI*h*h*h*h*h*h);
}

float Calculate_Reverce::SPH_function(Particle i, Particle j, float jfunc, float ifunc /*= 0.0f*/)
{
	float h_2 = h * h;
	vector3 temp;
	vec3_set(&temp, 0.0f, 0.0f, 0.0f);
	vec3_sub(&temp, &i.pos, &j.pos);

	float r_2 = vec3_dot(&temp, &temp);

	float _result = j.mass * poly6_coef * (h_2 - r_2) * (h_2 - r_2) * (h_2 - r_2);
	_result *= (jfunc + ifunc);//注意这里没有除以j粒子的密度，在jfunc自己去除，如果没有jfunc可以设为1.0

	return _result;
}

vector3 Calculate_Reverce::SPH_function_gradient(Particle i, Particle j, float jfunc, float ifunc)
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

float Calculate_Reverce::SPH_function_gradient(Particle i, Particle j, vector3 jfunc, vector3 ifunc /*= zeros()*/, float Adjustment /*= 1.0f*/)
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

vector3 Calculate_Reverce::SPH_function_gradient_2(Particle i, Particle j, vector3 jfunc, vector3 ifunc /*= zeros()*/, float Adjustment /*= 1.0f*/)
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
