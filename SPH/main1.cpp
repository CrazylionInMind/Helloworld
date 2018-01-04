//#include <MEMORY.H>    
//#include <WINDOWS.H>    
//#include <MATH.H>    
//#include <STDIO.H>    
//#include <ASSERT.H>    
//#define GLUT_DISABLE_ATEXIT_HACK
//#include <GL/glut.h>     
//#include<GL\GLAUX.H>//new add   
//
//#include<vector>
//#include<string>
//
//#include "cpu_sph.h"  
//#include"Marching_Cube.h"
//#include "SPHSystem.h"
//
//using namespace std;
//GLUquadricObj *quadratic;
//#define NUM_PHASE 2
//
//#define SMOOTHING_LENGTH 0.01f    
//#define SEARCH_RADIUS (1*SMOOTHING_LENGTH)    
//#define N_LOOPS 2   
//#define CUBE_LEN_X 10    
//#define CUBE_LEN_Y 10    
//#define LEN_Z 40   
//#define N_PARTICLES (CUBE_LEN_X*CUBE_LEN_Y*LEN_Z)    
//#define EPSILON 0.0001f    
//#define VISCOSITY 0.2f    
//#define STIFF 1.5f    
//#define MASS 0.00020543f    
//#define RANDF() ((rand()/(GLfloat)RAND_MAX - 0.5f))
//#define ISO_RADUIS 0.0115f
//#define MC_GRID_LEN (0.5*SMOOTHING_LENGTH)
//
////光照
//bool light;
//GLfloat LightAmbient[] = { 0.5f,0.5f,0.5f,1.0f };
//GLfloat LightDiffuse[] = { 1.0f,0.0f,0.0f,1.0f };
//GLfloat Lightposition[] = { 0.0f,0.0f,1.0f,1.0f };
//
//float anglex = 0.0f;
//float anglez = 0.0f;
//float tranx = 0.0f;
//float trany = 0.0f;
//float tranz = -0.07f;
//
//BOOL cal = TRUE;
//BOOL run = FALSE;
//BOOL bStepOne = FALSE;
//const float radius = 0.002f;
//vector<vector3> vdata;//暂存粒子位置供渲染    
//int n_steps;
////cpu_sph cpu;
//int iCurFrame = 0;
//Net_grid mcgrid;
//double viewmatrix[16], modelviewmatrix[16];
//
//void set_projview(GLfloat fov)
//{
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//
//	gluPerspective(fov, (GLfloat)800 / (GLfloat)600, 0.005, 30.0);
//
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//
//	gluLookAt(
//		0.0, 0.15, 0.0,
//		0.0, 0.0, 0.0,
//		0.0, 0.0, 1.0);
//}
//
//
//void init_obstacles()
//{
//	vector3 tmp;
//	matrix4 m;
//	vec3_set(&tmp, 0.0f, 0.0f, -0.0f);//tmp3个分量全部赋值0    
//	mat4_set_translate(&m, tmp.x, tmp.y, tmp.z);
//	cpu_sph_transform_obstacles(SPHSystem::getinstance(), &m);
//}
//
//void init_particles()
//{
//	int i;
//	int x;
//	int y;
//	int z;
//	GLfloat s;
//	GLfloat s2;
//	GLfloat s3;
//	GLfloat cx;
//	GLfloat cy;
//	GLfloat cz;
//
//	s = 0.006f;
//	cx = cy = 0.0f;
//	cz = 0.035f;
//
//	s2 = 0.001f;
//	s3 = 0.0001f;
//
//	//初始化流体,流体中N_PARTICLES个粒子，粒子有NUM_PHASE相
//	
//	//保留功能的时候，这里记得释放空间，正常结束释放在SPHSystem的析构函数中
//	Fluid *input = new Fluid(N_PARTICLES, NUM_PHASE);
//
//	for (x = 0; x < CUBE_LEN_X; x++)
//	{
//		for (y = 0; y < CUBE_LEN_Y; y++)
//		{
//			for (z = 0; z < LEN_Z; z++)
//			{
//				i = x + y*CUBE_LEN_X + z*CUBE_LEN_X*CUBE_LEN_Y;
//				if(input->m_ParticleVector.size() < i)
//					continue;
//				//设置位置
//				vec3_set(&(input->GetParticleByIndex(i)->pos), s*(x - CUBE_LEN_X / 2) - cx, s*(y - CUBE_LEN_Y / 2) - cy, 0.8*s*z - cz + 0.1f);
//				//设置速度
//				vec3_set(&(input->GetParticleByIndex(i)->vel), 0.0f, 0.0f, 0.0f);
//				input->GetParticleByIndex(i)->vel_half = input->GetParticleByIndex(i)->vel;
//				//以x轴为分界，x取正为红色，x取负为蓝色
//				if (x >= CUBE_LEN_X / 2)
//				{
//					//设置体积分数
//					input->GetParticleByIndex(i)->A_k[0]=(1.0f);//设立体积分数，第一项为1
//					input->GetParticleByIndex(i)->A_k[1]=(0.0f);//第二项为0
//					//设置质量分数
//					input->GetParticleByIndex(i)->C_k[0] = (1.0f);//设立体积分数，第一项为1
//					input->GetParticleByIndex(i)->C_k[1] = (0.0f);//第二项为0
//				}
//				else
//				{
//					input->GetParticleByIndex(i)->A_k[0]=(0.0f);//设立体积分数，第一项为1
//					input->GetParticleByIndex(i)->A_k[1]=(1.0f);//第二项为0
//
//					input->GetParticleByIndex(i)->C_k[0] = (0.0f);//设立体积分数，第一项为1
//					input->GetParticleByIndex(i)->C_k[1] = (1.0f);//第二项为0
//				}	
//				//设置每个粒子的质量
//				input->GetParticleByIndex(i)->SetMass(MASS);
//				//设置粘度
//				input->GetParticleByIndex(i)->SetViscosity(VISCOSITY);
//			}
//		}
//	}
//	SPHSystem::getinstance()->SetFluid(input);
//	SPHSystem::getinstance()->SetParam(SMOOTHING_LENGTH, 0.0008f,STIFF,SEARCH_RADIUS);
//	iCurFrame = 0;
//}
//
//int init()
//{
//	n_steps = 0;
//	init_particles();
//	init_obstacles();
//
//	vdata.resize(N_PARTICLES);
//
//	glClearColor(0.0, 0.0, 0.0, 0.0);
//
//	glShadeModel(GL_SMOOTH);
//	glClearDepth(1.0f);
//	glEnable(GL_DEPTH_TEST);
//	glDepthFunc(GL_LEQUAL);
//	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
//
//	//开下光照
//	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
//	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
//	glLightfv(GL_LIGHT1, GL_POSITION, Lightposition);
//	//glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT1);
//
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
//
//	return TRUE;
//}
//
//
//void reshape(int w, int h)
//{
//	glViewport(0, 0, w, h);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//}
//
//void color_by_power(float p)
//{
//	if (p < 1100.0f)
//	{
//		glColor3f(0.0f, 0.0f, p / 1100.0f);
//	}
//	else if (p < 1400.0f)
//	{
//		glColor3f(0.0f, p / 1400.0f, p / 1100.0f);
//	}
//	else if (p < 1900.0f)
//	{
//		glColor3f(p / 2100.0f, p / 1500.0f, 0.0f);
//	}
//	else glColor3f(p / 2100.0f, 0.1f, 0.0f);
//}
////mc密度使用的
////void color_by_power(float p)
////{
////	if (p < 700.0f)
////	{
////		glColor3f(0.0f, 0.0f, 1.0f);
////	}
////	else if (p < 900.0f)
////	{
////		glColor3f(0.0f, 1.0f, 0.0f);
////	}
////	else if (p < 1200.0f)
////	{
////		glColor3f(1.0f, 0.0f, 0.0f);
////	}
////	else glColor3f(1.0f, 1.0f, 0.5f);
////}
//
//void render_fluid()
//{
//	Fluid *temp = SPHSystem::getinstance()->GetFluid();
//	if (!temp)
//		return;
//	glPointSize(6.0);//以像素点为单位   
//	glBegin(GL_POINTS);
//	for (size_t i = 0; i < temp->GetParticleNum(); i++)
//	{
//		Particle *nowpar = temp->GetParticleByIndex(i);
//		if(!nowpar)
//			continue;
//		vector3 nowpos = nowpar->GetPosition();
//		glVertex3f(nowpos.x, nowpos.y, nowpos.z);
//	}
//	glEnd();
//}
//void drawxyzline(float length)
//{
//	glBegin(GL_LINES);
//	glColor3f(1.0f, 0.0f, 0.0f);
//	glVertex3f(0.0f, 0.0f, 0.0f);
//	glVertex3f(length, 0.0f, 0.0f);
//
//	glColor3f(0.0f, 1.0f, 0.0f);
//	glVertex3f(0.0f, 0.0f, 0.0f);
//	glVertex3f(0.0f, length, 0.0f);
//
//	glColor3f(0.0f, 0.0f, 1.0f);
//	glVertex3f(0.0f, 0.0f, 0.0f);
//	glVertex3f(0.0f, 0.0f, length);
//	glEnd();
//}
//void display(void)
//{
//	if (cal || bStepOne)
//	{
//		//cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);//模拟SPH液体和固体 0.003f
//		SPHSystem::getinstance()->CalSphOneStep();
//		iCurFrame++;
//		if (!cal) bStepOne = FALSE;
//	}
//	glEnable(GL_DEPTH_TEST);
//	/*   glEnable(GL_CULL_FACE); */
//	glPushMatrix();
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	set_projview(60);
//
//	glTranslatef(tranx, trany, tranz);
//	glRotatef(anglex, 1, 0, 0);
//	glRotatef(anglez, 0, 0, 1);
//	glGetDoublev(GL_MODELVIEW_MATRIX, viewmatrix);
//
//	glEnable(GL_BLEND);
//	render_fluid();//绘制流体  
//	glDisable(GL_BLEND);
//	
//	glPopMatrix();
//	glutSwapBuffers();
//}
//
//void idle(void)
//{
//	display();
//}
//
//void key(unsigned char key, int x, int y)
//{
//	switch (key)
//	{
//	case 'a':
//		tranx -= 0.01f;
//		break;
//
//	case 'd':
//		tranx += 0.01f;
//		break;
//
//	case 'q':
//		trany += 0.01f;
//		break;
//
//	case 'e':
//		trany -= 0.01f;
//		break;
//
//	case 's':
//		tranz += 0.01f;
//		break;
//
//	case 'w':
//		tranz -= 0.01f;
//		break;
//
//	case 'i':
//		anglex += 1.0f;
//		break;
//
//	case 'k':
//		anglex -= 1.0f;
//		break;
//
//	case 'j':
//		anglez += 1.0f;
//		break;
//
//	case 'l':
//		anglez -= 1.0f;
//		break;
//	case '0':
//		if (cal)
//		{
//			cal = FALSE;
//		}
//		else
//		{
//			cal = TRUE;
//		}
//		break;
//	case '1':
//		bStepOne = TRUE;
//		break;
//	case '2':
//	{
//		//int *mark = new int[cpu.n_particles];
//		//for (int i = 0;i < cpu.n_particles;i++)
//		//{
//		//	mark[i] = cpu.mark[i];
//		//}
//		//int n = max(iCurFrame - 1, 0);
//		//iCurFrame = 0;
//		//init_particles();
//		//for (int i = 0;i < n;i++)
//		//{
//		//	cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);//模拟SPH液体和固体 0.003f
//		//	iCurFrame++;
//		//}
//		//for (int i = 0;i < cpu.n_particles;i++)
//		//{
//		//	cpu.mark[i] = mark[i];
//		//}
//		//delete[]mark;
//		//printf("curframe:%d\n", iCurFrame);
//		//glutPostRedisplay();
//	}
//	case ' ':
//		if (run)
//		{
//			glutIdleFunc(NULL);
//			run = FALSE;
//		}
//		else
//		{
//			run = TRUE;
//			glutIdleFunc(idle);
//		}
//		glutPostRedisplay();
//		break;
//
//	case 'r':
//		init_particles();
//		glutPostRedisplay();
//		break;
//
//
//	case 'x':
//		exit(0);
//	case 'f':
//		if (light != true)
//		{
//			glEnable(GL_LIGHTING);
//			light = true;
//		}
//		else
//		{
//			glDisable(GL_LIGHTING);
//			light = false;
//		}
//	}
//}
//void findpoint(float x, float y, float z, cpu_sph* cpu, int mousekey)
//{
//	int *mark = new int[cpu->n_particles];
//	for (int i = 0;i < cpu->n_particles;i++)
//	{
//		mark[i] = cpu->mark[i];
//	}
//	float min = 999.0f;
//	int intomark = 0;
//	for (int i = 0; i < cpu->n_particles; i++)
//	{
//		cpu->mark[i] = 0;//这里重置了
//		float dis = (cpu->pos[i].x - x)*(cpu->pos[i].x - x) + (cpu->pos[i].y - y)*(cpu->pos[i].y - y) + (cpu->pos[i].z - z)*(cpu->pos[i].z - z);
//		if (min > dis)
//		{
//			intomark = i;
//			min = dis;
//		}
//	}
//	if (mousekey == 1)
//	{
//		delete[]mark;
//		return;
//	}
//	if (min>10)
//	{
//		for (int i = 0;i < cpu->n_particles;i++)
//		{
//			cpu->mark[i] = mark[i];
//		}
//		delete[]mark;
//		return;
//	}
//	if (mousekey == 0)
//	{
//		mark[intomark] = 2;
//		printf_s("This is the number of No. %d\n", intomark);
//	}
//	else if (mousekey == 2)
//	{
//		printf_s("第%d号粒子第一项的体积分数为：%f;第二项的体积分数为:%f\n\n", intomark, cpu->A_k[intomark], cpu->A_k[intomark + 4000]);
//	}
//
//	for (int i = 0;i < cpu->n_particles;i++)
//	{
//		cpu->mark[i] = mark[i];
//	}
//	delete[]mark;
//
//
//}
//
//void mytimefun(int sm)
//{
//	display();
//	string title = "This is the NO.";
//	title += to_string(iCurFrame);
//	title += " Frame";
//	glutSetWindowTitle(title.c_str());
//	glutTimerFunc(10, mytimefun, 0);
//}
//int main(int argc, char** argv)
//{
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
//	glutInitWindowSize(800, 600);
//	glutCreateWindow("Liquid");//argv[0]);    
//
//	glutReshapeFunc(reshape);
//	glutKeyboardFunc(key);//使用空格控制粒子的模拟    
//
//	init();//初始化粒子，刚体    
//
//	//glutIdleFunc(idle);//已经按时间开始计算了
//	glutTimerFunc(10, mytimefun, 0);
//	glutDisplayFunc(display);//绘制粒子，刚体  
//	//glutMouseFunc(&myMouseFunc);
//
//	glutMainLoop();
//	return 0;
//}
