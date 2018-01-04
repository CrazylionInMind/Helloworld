/*
*/

#include <MEMORY.H>    
#include <WINDOWS.H>    
#include <MATH.H>    
#include <STDIO.H>    
#include <ASSERT.H>    
#define GLUT_DISABLE_ATEXIT_HACK
#include <GL/glut.h>     
#include<GL\GLAUX.H>//new add   

#include<vector>
#include<string>

#include "cpu_sph.h"  
#include"Marching_Cube.h"
#include "SPHSystem.h"
#include "MyfileRW.h"
using namespace std;
GLUquadricObj *quadratic;

//#define FILEINPUT
//光照
bool light;
GLfloat LightAmbient[] = { 0.5f,0.5f,0.5f,1.0f };
GLfloat LightDiffuse[] = { 1.0f,0.0f,0.0f,1.0f };
GLfloat Lightposition[] = { 0.0f,0.0f,1.0f,1.0f };

float anglex = 0.0f;
float anglez = 0.0f;
float tranx = 0.0f;
float trany = 0.0f;
float tranz = -0.07f;

BOOL cal = TRUE;
BOOL run = FALSE;
BOOL bStepOne = FALSE;
const float radius = 0.002f;
int iCurFrame = 0;
double viewmatrix[16], modelviewmatrix[16];

void set_projview(GLfloat fov)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(fov, (GLfloat)800 / (GLfloat)600, 0.005, 30.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(
		0.0, 0.15, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0);
}

void init_obstacles()
{
	vector3 tmp;
	matrix4 m;
	vec3_set(&tmp, 0.0f, 0.0f, -0.0f);//tmp3个分量全部赋值0    
	mat4_set_translate(&m, tmp.x, tmp.y, tmp.z);
	cpu_sph_transform_obstacles(SPHNewSystem::getinstance(), &m);
}

void init_particles()
{
	SPHNewSystem::getinstance()->init();
	iCurFrame = 0;
}

int init()
{
	init_particles();
	init_obstacles();

	glClearColor(0.0, 0.0, 0.0, 0.0);

	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	//开下光照
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, Lightposition);
	//glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT1);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	return TRUE;
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void color_by_power(float p)
{
	if (p < 1100.0f)
	{
		glColor3f(0.0f, 0.0f, p / 1100.0f);
	}
	else if (p < 1400.0f)
	{
		glColor3f(0.0f, p / 1400.0f, p / 1100.0f);
	}
	else if (p < 1900.0f)
	{
		glColor3f(p / 2100.0f, p / 1500.0f, 0.0f);
	}
	else glColor3f(p / 2100.0f, 0.1f, 0.0f);
}

void render_fluid()
{
	Fluid *temp = SPHNewSystem::getinstance()->GetFluid();
	if (!temp)
		return;
	glPointSize(6.0);//以像素点为单位   
	glBegin(GL_POINTS);
	for (size_t i = 0; i < temp->GetParticleNum(); i++)
	{
		Particle *nowpar = temp->GetParticleByIndex(i);
		if (!nowpar)
			continue;
		vector3 nowpos = nowpar->GetPosition();
		glVertex3f(nowpos.x, nowpos.y, nowpos.z);
	}
	glEnd();
}

void drawxyzline(float length)
{
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(length, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, length, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, length);
	glEnd();
}

void display(void)
{
	if (cal || bStepOne)
	{
		//cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);
		SPHNewSystem::getinstance()->CalSphOneStep();
		iCurFrame++;
		if (!cal) bStepOne = FALSE;

#ifdef FILEINPUT
		MyFileRW::getinstance()->SetWriteFileName("output.txt");

		Fluid *tempfluid = SPHNewSystem::getinstance()->GetFluid();
		if(tempfluid)
			for (int i = 0; i < tempfluid->GetParticleNum(); i++)
			{
				Particle *temp = tempfluid->GetParticleByIndex(i);
				MyFileRW::getinstance()->WriteSomething(i);
				MyFileRW::getinstance()->WriteSomething(temp->pos.x);
				MyFileRW::getinstance()->WriteSomething(temp->pos.y);
				MyFileRW::getinstance()->WriteSomething(temp->pos.z,1);
			}
		MyFileRW::getinstance()->CloseWriteFile();
#endif
	}
	glEnable(GL_DEPTH_TEST);
	/*   glEnable(GL_CULL_FACE); */
	glPushMatrix();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	set_projview(60);

	glTranslatef(tranx, trany, tranz);
	glRotatef(anglex, 1, 0, 0);
	glRotatef(anglez, 0, 0, 1);
	glGetDoublev(GL_MODELVIEW_MATRIX, viewmatrix);

	glEnable(GL_BLEND);
	render_fluid();//绘制流体  
	glDisable(GL_BLEND);

	glPopMatrix();
	glutSwapBuffers();
}

void idle(void)
{
	display();
}

void key(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'a':
		tranx -= 0.01f;
		break;

	case 'd':
		tranx += 0.01f;
		break;

	case 'q':
		trany += 0.01f;
		break;

	case 'e':
		trany -= 0.01f;
		break;

	case 's':
		tranz += 0.01f;
		break;

	case 'w':
		tranz -= 0.01f;
		break;

	case 'i':
		anglex += 1.0f;
		break;

	case 'k':
		anglex -= 1.0f;
		break;

	case 'j':
		anglez += 1.0f;
		break;

	case 'l':
		anglez -= 1.0f;
		break;
	case '0':
		if (cal)
		{
			cal = FALSE;
		}
		else
		{
			cal = TRUE;
		}
		break;
	case '1':
		bStepOne = TRUE;
		break;
	case '2':
	{
		//int *mark = new int[cpu.n_particles];
		//for (int i = 0;i < cpu.n_particles;i++)
		//{
		//	mark[i] = cpu.mark[i];
		//}
		//int n = max(iCurFrame - 1, 0);
		//iCurFrame = 0;
		//init_particles();
		//for (int i = 0;i < n;i++)
		//{
		//	cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);//模拟SPH液体和固体 0.003f
		//	iCurFrame++;
		//}
		//for (int i = 0;i < cpu.n_particles;i++)
		//{
		//	cpu.mark[i] = mark[i];
		//}
		//delete[]mark;
		//printf("curframe:%d\n", iCurFrame);
		//glutPostRedisplay();
	}
	case ' ':
		if (run)
		{
			glutIdleFunc(NULL);
			run = FALSE;
		}
		else
		{
			run = TRUE;
			glutIdleFunc(idle);
		}
		glutPostRedisplay();
		break;

	case 'r':
		init_particles();
		glutPostRedisplay();
		break;


	case 'x':
		exit(0);
	case 'f':
		if (light != true)
		{
			glEnable(GL_LIGHTING);
			light = true;
		}
		else
		{
			glDisable(GL_LIGHTING);
			light = false;
		}
	}
}

void findpoint(float x, float y, float z, cpu_sph* cpu, int mousekey)
{
	int *mark = new int[cpu->n_particles];
	for (int i = 0;i < cpu->n_particles;i++)
	{
		mark[i] = cpu->mark[i];
	}
	float min = 999.0f;
	int intomark = 0;
	for (int i = 0; i < cpu->n_particles; i++)
	{
		cpu->mark[i] = 0;//这里重置了
		float dis = (cpu->pos[i].x - x)*(cpu->pos[i].x - x) + (cpu->pos[i].y - y)*(cpu->pos[i].y - y) + (cpu->pos[i].z - z)*(cpu->pos[i].z - z);
		if (min > dis)
		{
			intomark = i;
			min = dis;
		}
	}
	if (mousekey == 1)
	{
		delete[]mark;
		return;
	}
	if (min > 10)
	{
		for (int i = 0;i < cpu->n_particles;i++)
		{
			cpu->mark[i] = mark[i];
		}
		delete[]mark;
		return;
	}
	if (mousekey == 0)
	{
		mark[intomark] = 2;
		printf_s("This is the number of No. %d\n", intomark);
	}
	else if (mousekey == 2)
	{
		printf_s("第%d号粒子第一项的体积分数为：%f;第二项的体积分数为:%f\n\n", intomark, cpu->A_k[intomark], cpu->A_k[intomark + 4000]);
	}

	for (int i = 0;i < cpu->n_particles;i++)
	{
		cpu->mark[i] = mark[i];
	}
	delete[]mark;
}

void mytimefun(int sm)
{
	display();
	string title = "This is the NO.";
	title += to_string(iCurFrame);
	title += " Frame";
	glutSetWindowTitle(title.c_str());
	glutTimerFunc(10, mytimefun, 0);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Liquid");//argv[0]);    

	glutReshapeFunc(reshape);
	glutKeyboardFunc(key);//使用空格控制粒子的模拟    

	init();//初始化粒子，刚体    

	//glutIdleFunc(idle);//已经按时间开始计算了
	glutTimerFunc(10, mytimefun, 0);
	glutDisplayFunc(display);//绘制粒子，刚体  
	//glutMouseFunc(&myMouseFunc);

	glutMainLoop();
	return 0;
}
