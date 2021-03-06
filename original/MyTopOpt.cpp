// MyTopOpt.cpp: 定义控制台应用程序的入口点。
//
#include<windows.h>
#include <GL/GLU.H>
#include <GL/freeglut.h>
#include<GL/GL.H>
#include"TopOpt.h"

const static int WINDOW_WIDTH = 700;//显示区域
const static int WINDOW_HEIGHT = 700;

TopOpt* mytop;

void render()
{
	//mytop->draw();
	Eigen::MatrixXd x;
	x = mytop->X;

	int nely = mytop->nely; int nelx = mytop->nelx;
	//glBegin();
	for (int i = 0; i < nely; i++)
	{
		for (int j = 0; j < nelx; j++)
		{
			float color = 1 - x(i, j);
			glColor3f(color, color, color);
			glRectf(float(-nelx + j) / 10, float(nely - i) / 10, float(-nelx + j + 1) / 10, float(nely - i - 1) / 10);
		}
	}
	glFlush();
	//glEnd();
}

void update()
{

}

int main(int argc, char** argv)
{
	long t1 = GetTickCount();
	
	
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(0, 0);
	glutInit(&argc, argv);  //对GLUT进行初始化，这个函数必须在其它的GLUT使用之前调用一次
	glutCreateWindow("MPMElastic");
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);//设置显示方式，其中GLUT_RGB表示使用RGB颜色，与之对应的还有GLUT_INDEX（表示使用索引颜色）。GLUT_SINGLE表示使用单缓冲，与之对应的还有GLUT_DOUBLE（使用双缓冲）


	//mytop = new TopOpt();
	mytop = new TopOpt(3, 2, 0.4, 3.0, 2.0);
	mytop->compute();
	long t2 = GetTickCount();
	cout << "运行时间：" << (t2 - t1) << endl;


	glutDisplayFunc(render);
	glutMainLoop();


	system("pause");
    return 0;
}

