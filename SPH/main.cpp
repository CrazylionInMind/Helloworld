//#include <MEMORY.H>    
//#include <WINDOWS.H>    
//#include <MATH.H>    
//#include <STDIO.H>    
//#include <ASSERT.H>    
//#define GLUT_DISABLE_ATEXIT_HACK
//#include <GL/glut.h>     
//#include<GL\GLAUX.H>//new add   
//#include "cpu_sph.h"  
//#include"Marching_Cube.h"
//#include<vector>
//#include<string>
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
//vector3* vdata;//暂存粒子位置供渲染    
//int n_steps;   
//cpu_sph cpu;   
//int iCurFrame = 0;
//Net_grid mcgrid;
//double viewmatrix[16], modelviewmatrix[16];
//
////GLuint texture[1];//暂存纹理
////				  //Loads A Bitmap mage
////AUX_RGBImageRec *LoadBMP(char *Filename)				// Loads A Bitmap Image
////{
////	FILE *File = NULL;									// File Handle
////	if (!Filename)										// Make Sure A Filename Was Given
////	{
////		return NULL;									// If Not Return NULL
////	}
////	File = fopen(Filename, "r");							// Check To See If The File Exists
////	if (File)											// Does The File Exist?
////	{
////		fclose(File);									// Close The Handle
////		return auxDIBImageLoad(Filename);				// Load The Bitmap And Return A Pointer
////	}
////	return NULL;										// If Load Failed Return NULL
////}
////先试试简单的绑定纹理
////int LoadGLTextures()									// Load Bitmaps And Convert To Textures
////{
////	int Status = FALSE;									// Status Indicator
////
////	AUX_RGBImageRec *TextureImage[1];					// Create Storage Space For The Texture
////
////	memset(TextureImage, 0, sizeof(void *) * 1);           	// Set The Pointer To NULL
////
////															// Load The Bitmap, Check For Errors, If Bitmap's Not Found Quit
////	if (TextureImage[0] = LoadBMP("Data/metal3.bmp"))//注意在这里改名字
////	{
////		Status = TRUE;									// Set The Status To TRUE
////
////		glGenTextures(1, &texture[0]);					// Create The Texture
////
////														// Typical Texture Generation Using Data From The Bitmap
////		glBindTexture(GL_TEXTURE_2D, texture[0]);
////		glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage[0]->sizeX, TextureImage[0]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[0]->data);
////		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
////		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
////	}
////
////	if (TextureImage[0])									// If Texture Exists
////	{
////		if (TextureImage[0]->data)							// If Texture Image Exists
////		{
////			free(TextureImage[0]->data);					// Free The Texture Image Memory
////		}
////
////		free(TextureImage[0]);								// Free The Image Structure
////	}
////
////	return Status;										// Return The Status
////}
//
//
//
//
//void set_projview(GLfloat fov)   
//{   
//    glMatrixMode(GL_PROJECTION);   
//    glLoadIdentity();   
//   
//    gluPerspective(fov, (GLfloat)800/(GLfloat)600,0.005,30.0);   
//   
//    glMatrixMode(GL_MODELVIEW);   
//    glLoadIdentity();   
//       
//    gluLookAt(   
//        0.0, 0.15, 0.0,   
//        0.0, 0.0, 0.0,   
//        0.0, 0.0, 1.0);   
//}   
//   
//  
//void init_obstacles()   
//{   
//    vector3 tmp;   
//    matrix4 m;   
//    vec3_set(&tmp, 0.0f, 0.0f, -0.0f);//tmp3个分量全部赋值0    
//    mat4_set_translate(&m, tmp.x, tmp.y, tmp.z);   
//    cpu_sph_transform_obstacles(&cpu, &m);   
//}   
//   
///**  
// * Initialize the states of particles  
// */   
//void init_particles()   
//{   
//    int i;   
//    int x;   
//    int y;   
//    int z;   
//    GLfloat s;   
//    GLfloat s2;   
//    GLfloat s3;   
//    GLfloat cx;   
//    GLfloat cy;   
//    GLfloat cz;   
//   
//    vector3 pos[N_PARTICLES];   
//    vector3 vel[N_PARTICLES];   
//	float a_k[NUM_PHASE*N_PARTICLES];
//
//    s = 0.006f;   
//    cx = cy = 0.0f;   
//    cz = 0.035f;   
//   
//    s2 = 0.001f;   
//    s3 = 0.0001f;   
//   
//    for (x = 0; x < CUBE_LEN_X; x++)   
//    {   
//        for (y = 0; y < CUBE_LEN_Y; y++)   
//        {          
//            for (z = 0; z < LEN_Z; z++)   
//            {   
//                i = x + y*CUBE_LEN_X + z*CUBE_LEN_X*CUBE_LEN_Y;   
//               // vec3_set(&pos[i], s*(x -CUBE_LEN_X/2)- cx + s2*RANDF(), s*(y - CUBE_LEN_Y/2) - cy + s2*RANDF(),0.8*s*z - cz + 0.1f);  
//				vec3_set(&pos[i], s*(x - CUBE_LEN_X / 2) - cx, s*(y - CUBE_LEN_Y / 2) - cy, 0.8*s*z - cz + 0.1f);
//                vec3_set(&vel[i], 0.0f,  0.0f, 0.0f);   
//				//以x轴为分界，x取正为红色，x取负为蓝色
//				if (x>=CUBE_LEN_X/2)
//				{
//					a_k[i] = 1.0f;//第一项为1
//					a_k[i + N_PARTICLES] = 0.0f;//第二项为0
//				}
//				else
//				{
//					a_k[i] = 0.0f;//第一项为0
//					a_k[i + N_PARTICLES] = 1.0f;//第二项为1
//				}
//            }   
//        }   
//    }   
//   
//    cpu_sph_create(&cpu, N_PARTICLES, pos, vel, SMOOTHING_LENGTH, VISCOSITY, MASS, STIFF, SEARCH_RADIUS, N_LOOPS);   
//	
//	//对于多相申请空间
//	mult_create(&cpu , NUM_PHASE,a_k);
//    
//	cpu.n_particles = N_PARTICLES;   
//	iCurFrame = 0;
//}   
//   
//int init()   
//{   
//    n_steps = 0;   
//    init_particles();   
//    init_obstacles();   
//   
//    vdata = vec3_alloc(N_PARTICLES);   
//	//if (!LoadGLTextures())								// Jump To Texture Loading Routine ( NEW )
//	//{
//	//	return FALSE;									// If Texture Didn't Load Return FALSE
//	//}
//	glClearColor(0.0, 0.0, 0.0, 0.0);
//
//	/*glEnable(GL_TEXTURE_2D);*/
//	glShadeModel(GL_SMOOTH);
//	glClearDepth(1.0f);
//	glEnable(GL_DEPTH_TEST);
//	glDepthFunc(GL_LEQUAL);
//	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
//	
//
//	////创建球形多边形的申明，并绑定自动生成纹理
//	//quadratic = gluNewQuadric();							// Create A Pointer To The Quadric Object (Return 0 If No Memory)
//	//gluQuadricNormals(quadratic, GLU_SMOOTH);			// Create Smooth Normals 
//	//gluQuadricTexture(quadratic, GL_TRUE);
//	//glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP); // Set The Texture Generation Mode For S To Sphere Mapping (NEW)
//	//glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP); // Set The Texture Generation Mode For T To Sphere Mapping (NEW)
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
//    glViewport(0, 0, w, h);   
//    glMatrixMode(GL_PROJECTION);   
//    glLoadIdentity();   
//    glMatrixMode(GL_MODELVIEW);   
//    glLoadIdentity();   
//}   
//   
//void color_by_power(float p)   
//{   
//   if (p < 1100.0f)   
//    {  
//        glColor3f(0.0f, 0.0f, p/1100.0f);   
//    }   
//    else if (p < 1400.0f)   
//    {   
//        glColor3f(0.0f, p/1400.0f, p/1100.0f);   
//    }   
//    else if (p < 1900.0f)   
//    {   
//        glColor3f(p / 2100.0f, p/1500.0f, 0.0f);   
//    }   
//    else glColor3f(p / 2100.0f, 0.1f, 0.0f);  
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
//    //int i;   
//    cpu_sph_get_pos(&cpu, vdata);//获取粒子位置供渲染
//	
//	float *ak = cpu.A_k;
//	//从这里开始编写Marching Cube 算法
//	//已知条件仅为：粒子坐标      等值面所需要的值
//
//	
//	//init_mc_paramter(&mcgrid, ISO_RADUIS, MC_GRID_LEN, MASS*315.0f / ((64.0f*M_PI)*pow(ISO_RADUIS, 9)), &cpu);
//	//sph_render_draw_surface(&cpu,700.0f,&mcgrid);
//	
//	
//	//glPointSize(4.0);//以网格左下角为粒子中心，查看颜色划分
//	//glBegin(GL_POINTS);   
//	//for (int x = 0; x < mcgrid.width; x++)
//	//{
//	//	for (int y = 0; y < mcgrid.height; y++)
//	//	{
//	//		for (int z = 0; z < mcgrid.depth; z++)
//	//		{
//	//			color_by_power( mcgrid.density[x+y*mcgrid.width+z*mcgrid.width*mcgrid.height]);
//	//			glVertex3f(mcgrid.grid_position.getX()+x*mcgrid.voxelsize, mcgrid.grid_position.getY()+y*mcgrid.voxelsize, mcgrid.grid_position.getZ()*mcgrid.voxelsize);
//	//		}
//	//	}
//	//}	    	  
//	//glEnd();   
//
//
//	//glMatrixMode(GL_MODELVIEW);
//	//glColor3ub(149, 195, 203);
//	////glBegin(GL_LINE_STRIP);
//	//glBegin(GL_QUAD_STRIP);
//	//for ( i = 0; i < cpu.n_particles; i++)
//	//{
//	//	glVertex3f(vdata[i].x, vdata[i].y, vdata[i].z);
//	//}
//	//glVertex3f(vdata[0].x, vdata[0].y, vdata[0].z);
//	//
//	////glEnable(GL_TEXTURE_GEN_S);
//	////glEnable(GL_TEXTURE_GEN_T);
//	////glBindTexture(GL_TEXTURE_2D, texture[0]);
//	////glTranslatef(0.0f, -0.01f, .0f);
//	////for (i = 0; i < cpu.n_particles; i++)
//	////{
//	////	
//	////	glPushMatrix();
//	////	glTranslatef(vdata[i].x, vdata[i].y, vdata[i].z);
//	////	gluSphere(quadratic, 0.005f, 16, 16);
//	////	glPopMatrix();
//	////}
//	////glDisable(GL_TEXTURE_GEN_S);
//	////glDisable(GL_TEXTURE_GEN_T);
//	//glEnable(GL_TEXTURE_GEN_S);
//	//glEnable(GL_TEXTURE_GEN_T);
//	//glBindTexture(GL_TEXTURE_2D, texture[0]);
//	//glPushMatrix();
//	//glTranslatef(0.0f, -30.0f, 0.0f);//第二参数是前后，第三参数是上下，第一参数是左右
//	//gluSphere(quadratic, .40f, 32, 32);
//	//glPopMatrix();
//	//glDisable(GL_TEXTURE_GEN_S);
//	//glDisable(GL_TEXTURE_GEN_T);
//
//    glPointSize(6.0);//以像素点为单位   
//    glBegin(GL_POINTS);   
//	std::vector<int> num_for_sur;
//    for (int i = 0; i< cpu.n_particles; i++)   
//    {   
//		//vector3 colors_every;
//		//vec3_set(&colors_every, 0, 0, 0);
//		float red=0.0f;
//		float blue=0.0f;
//		float green = 0.0f;
//		for (int j = 0; j < cpu.num_phase; j++)
//		{
//			if (j==0)
//			{
//				red = 1.0f*ak[i + j*cpu.n_particles];//体积分数
//			}
//			else
//			{
//				blue = 1.0f*ak[i + j*cpu.n_particles];
//			}
//		}
//	
//		glColor3f(red, 0, blue);
//		//if (blue - red>0.5)
//		//{
//		//	red = 1.0f;
//		//	green = 1.0f;
//		//	glColor3f(red, green, 0.0);
//		//}
//        //color_by_power(1.0f/cpu.density[i]);   
//		glVertex3fv((GLfloat*)&vdata[i]);   
//    }   
//    glEnd();   
//	for (int i = 0; i < cpu.n_particles; i++)
//	{
//		if (cpu.mark[i] == 1)
//		{
//			glColor3f(0, 1, 0);
//			glPushMatrix();
//			glTranslatef(cpu.pos[i].x, cpu.pos[i].y, cpu.pos[i].z);
//			glutWireSphere(0.002, 10, 10);
//			glPopMatrix();
//		}
//		if (cpu.mark[i]==2)
//		{
//			glColor3f(1, 0, 1);
//			glPushMatrix();
//			glTranslatef(cpu.pos[i].x, cpu.pos[i].y, cpu.pos[i].z);
//			glutWireSphere(0.002, 10, 10);
//			glPopMatrix();
//		}
//	}
//}   
//void drawxyzline(float length)
//{
//	glBegin(GL_LINES);
//		glColor3f(1.0f, 0.0f, 0.0f);	
//		glVertex3f(0.0f, 0.0f, 0.0f);
//		glVertex3f(length, 0.0f, 0.0f);
//
//		glColor3f(0.0f, 1.0f, 0.0f);
//		glVertex3f(0.0f, 0.0f, 0.0f);
//		glVertex3f(0.0f, length, 0.0f);
//
//		glColor3f(0.0f, 0.0f, 1.0f);
//		glVertex3f(0.0f, 0.0f, 0.0f);
//		glVertex3f(0.0f, 0.0f, length);
//	glEnd();
//}
//void display(void)   
//{  
//	if (cal||bStepOne)
//	{
//		cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);//模拟SPH液体和固体 0.003f
//		iCurFrame++;
//		if (!cal) bStepOne = FALSE;
//	}
//    glEnable(GL_DEPTH_TEST);   
// /*   glEnable(GL_CULL_FACE); */  
//	glPushMatrix();
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);   
//    set_projview(60);   
//	 
//	glTranslatef(tranx,trany,tranz);
//    glRotatef(anglex,1,0,0);   
//    glRotatef(anglez,0,0,1);
//	glGetDoublev(GL_MODELVIEW_MATRIX, viewmatrix);
//
//	//drawxyzline(10.0f);
//	//设置一个物体体现透明
//	//glPushMatrix();
//	//glTranslatef(0.0f, -10.0f, 0.0f);
//	//glColor3f(0.0f, 0.0f, 0.0f);
//	//glutSolidSphere(5.0f, 10, 10);
//	//glPopMatrix();
//		
//    /*glTranslatef(tranx,trany - 0.2f,tranz); */  
//	
//	//glColor4f(0.8f, 1.0f, 1.0f,0.5f);
//
////	glDisable(GL_DEPTH_TEST);
//	glEnable(GL_BLEND);
//    render_fluid();//绘制流体  
//	glDisable(GL_BLEND);
////	glEnable(GL_DEPTH_TEST);
//	glPopMatrix();
//    glutSwapBuffers();   
//}   
//   
//void idle(void)   
//{   
//   // cpu_sph_elapse(&cpu, 0.003f);//模拟SPH液体和固体    
//   // glutPostRedisplay();  
//	display();
//}   
//   
//void key(unsigned char key, int x, int y)   
//{   
//    switch(key)   
//    {   
//    case 'a':   
//        tranx -=0.01f;   
//        break;   
//   
//    case 'd':   
//        tranx +=0.01f;   
//        break;   
//   
//    case 'q':   
//    trany +=0.01f;   
//        break;   
//   
//    case 'e':   
//    trany -=0.01f;   
//        break;   
//   
//    case 's':   
//        tranz +=0.01f;   
//        break;   
//   
//    case 'w':   
//        tranz -=0.01f;   
//        break;   
//   
//    case 'i':   
//        anglex +=1.0f;   
//        break;   
//   
//    case 'k':   
//        anglex -=1.0f;   
//        break;   
//   
//    case 'j':   
//        anglez +=1.0f;   
//        break;   
//   
//    case 'l':   
//        anglez -=1.0f;   
//        break;  
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
//		int *mark = new int[cpu.n_particles];
//		for (int i = 0;i < cpu.n_particles;i++)
//		{
//			mark[i] = cpu.mark[i];
//		}
//		int n = max(iCurFrame - 1, 0);
//		iCurFrame = 0;
//		init_particles();
//		for (int i = 0;i < n;i++)
//		{
//			cpu_sph_elapse(&cpu, 0.0008f, iCurFrame);//模拟SPH液体和固体 0.003f
//			iCurFrame++;
//		}
//		for (int i = 0;i < cpu.n_particles;i++)
//		{
//			cpu.mark[i] = mark[i];
//		}
//		delete[]mark;
//		printf("curframe:%d\n", iCurFrame);
//		glutPostRedisplay();
//	}
//    case ' ':   
//        if (run)   
//        {   
//            glutIdleFunc(NULL);   
//            run = FALSE;   
//        }   
//        else   
//        {   
//            run = TRUE;   
//            glutIdleFunc(idle);   
//        }   
//        glutPostRedisplay();   
//        break;   
//   
//    case 'r':   
//    init_particles();   
//    glutPostRedisplay();   
//        break;   
//   
//       
//    case 'x':   
//        exit(0);  
//	case 'f':
//		if (light!=true)
//		{
//			glEnable(GL_LIGHTING);
//			light = true;
//		}
//		else
//		{
//			glDisable(GL_LIGHTING);
//			light = false;
//		}	
//    }   
//}   
//void findpoint(float x,float y,float z, cpu_sph* cpu,int mousekey)
//{
//	int *mark = new int[cpu->n_particles];
//	for (int i = 0;i < cpu->n_particles;i++)
//	{
//		mark[i] = cpu->mark[i];
//	}
//	float min = 999.0f;
//	int intomark=0;
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
//	if (mousekey==1)
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
//	if (mousekey==0)
//	{
//		mark[intomark] = 2;
//		printf_s("This is the number of No. %d\n", intomark);
//	}
//	else if (mousekey==2)
//	{
//		printf_s("第%d号粒子第一项的体积分数为：%f;第二项的体积分数为:%f\n\n",intomark,cpu->A_k[intomark], cpu->A_k[intomark + 4000]);
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
//void myMouseFunc(int button, int state, int x, int y)
//{
//	if (state)
//	{
//		return;
//	}
//	float val;
//	double modelview[16], project[16], pos[3];
//	int viewport[4];
//	//glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
//	glGetDoublev(GL_PROJECTION_MATRIX, project);
//	glGetIntegerv(GL_VIEWPORT, viewport);
//	for (int i = 0; i < 4; i++)
//	{
//		printf_s("%d  ", viewport[i]);
//		if (i==3)
//		{
//			printf_s("\n");
//		}
//	}
//	y = viewport[3] - y;
//	glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &val);
//	gluUnProject(x, y, val, viewmatrix, project, viewport, &pos[0], &pos[1], &pos[2]);
//
//	printf("%d:%d\t(%d:%d)\t%f\t(%.2f,%.2f,%.2f)\n", button, state, x, y, val, pos[0], pos[1], pos[2]);
//	findpoint(pos[0],pos[1],pos[2],&cpu,button);
//
//}
//
//void mytimefun(int sm)
//{
//	display();
//	string title="This is the NO.";
//	title += to_string(iCurFrame);
//	title += " Frame";
//	glutSetWindowTitle(title.c_str());
//	glutTimerFunc(10, mytimefun, 0);
//}
//int main(int argc, char** argv)   
//{   
//    glutInit(&argc, argv);   
//    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);   
//    glutInitWindowSize(800, 600);   
//    glutCreateWindow("Liquid");//argv[0]);    
//   
//    glutReshapeFunc(reshape);   
//    glutKeyboardFunc(key);//使用空格控制粒子的模拟    
//   
//    init();//初始化粒子，刚体    
//   
//    //glutIdleFunc(idle);//已经按时间开始计算了
//	glutTimerFunc(10, mytimefun, 0);
//    glutDisplayFunc(display);//绘制粒子，刚体  
//	glutMouseFunc(&myMouseFunc);
//   
//    glutMainLoop();   
//    return 0;   
//}
