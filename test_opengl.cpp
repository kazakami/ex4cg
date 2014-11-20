#include <iostream>
#include <memory>
#include <random>

#include <stdlib.h>
#include <GL/glut.h>

#include "image.hpp"
#include "vrml.hpp"

constexpr int sx = 256;
constexpr int sy = 256;

typedef unsigned char uchar;


kazakami::Image3d im(sx, sy);
kazakami::Polygon test;

int GLframe = 0; //フレーム数
int GLtimenow = 0;//経過時間
int GLtimebase=0;//計測開始時間



/*------ render_viewer functions ------------------*/
static void redraw_image(void)
{
  static double l = 0;
  im.ResetZBuffer();
  for (int i = 0; i < sy; i++)
    for (int j = 0; j < sx; j++)
    {
      //*
      im.at(i, j, kazakami::colorRGB::B) = i*0;
      im.at(i, j, kazakami::colorRGB::G) = j*0;
      im.at(i, j, kazakami::colorRGB::R) = 0;
      //*/
    }

  l += 0.01;

  im.LoadIdentity();
  im.Perspective(256);
  //im.Scale(1 + 0.5*sin(l), 1 + 0.5*cos(l), 1);

  im.Rotate(l, 0, 0, 1);
  im.Scale(2, 2, 1);
  if ((int)(l*10) % 3 == 0)
  {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_int_distribution<int> dice(-8,8);
    im.dif = dice(mt);
  }
  //im.SetLightPos({250*cos(l), 0, 250+400*sin(l)});
  //im.LookAt({0, 0, 0}, {0, 0, 1}, {cos(l), sin(l), 0});
  //im.LookAt({500*sin(l), 0, 500*cos(l)}, {0, 0, 0}, {0, 1, 0});
  im.LookAt({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  test.SuperDraw(&im);
  



  /* Just for Safe */
  glPushMatrix();{

    /* Just for Safe Again */
    glLoadIdentity();

    /* With Identity Projection Matrix,                                                                                                       
          Window Screen Range is (-1,-1)[bottomleft] to (1,1)[upright].                                                                       
          Write the Buffer from (-1,-1), bottom left. */
    glRasterPos3d(-1.0, -1.0, 0);
    
    glDrawPixels(sx, sy, GL_BGR, GL_UNSIGNED_BYTE,
                 (const GLvoid *)(im.ptr()));
    
    glPopMatrix();
  }


  GLframe++; //フレーム数を＋１
  GLtimenow = glutGet(GLUT_ELAPSED_TIME);//経過時間を取得
        
  if (GLtimenow - GLtimebase > 1000)      //１秒以上たったらfpsを出力
  {
    std::cout << "\rfps:" << GLframe*1000.0/(GLtimenow-GLtimebase) << std::flush;
    GLtimebase = GLtimenow;//基準時間を設定                
    GLframe = 0;//フレーム数をリセット
  }


  /* Double Buffering */
  glutSwapBuffers();
}


int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "please filenameeeeeee!!" << std::endl;
    return 1;
  }
  if (argc >= 3)
  {
    auto tex = std::make_shared<kazakami::Texture>();
    tex->LoadPPM(argv[2]);
    if (tex->LoadPPM(argv[2]))
    {
      std::cerr << "Enable environment mapping" << std::endl;
      im.SetEnvMapTex(tex);
    }
    else
    {
      std::cerr << "invalid environment texture." << std::endl;
    }
  }
  im.SetOffset(sx/2,sy/2);
  im.LoadIdentity();
  im.SetLightPos({-500, -500, 500});
  im.SetLightCol(1, 1, 1);
  //Zバッファを有効にする
  im.Enable(kazakami::Image3d::Function::DepthTest);
  //隠面処理をenableにする
  im.Disable(kazakami::Image3d::Function::CullFace);
  /*
  if (argc >= 3)
    im.Translate(0, std::atoi(argv[2]), 0);
  */
#ifdef _OPENMP
  puts("MP is enable");
#else
  puts("MP is disable");
#endif
  
  test.ReadVRML(argv[1]);
  test.AddVertex(&im);

  /* Prepareing the GL and GLUT */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
  glutInitWindowSize(sx, sy);
  glutCreateWindow("LE4 CG");
  glutDisplayFunc(redraw_image);
  glutIdleFunc(redraw_image);

  /* MainLoop */
  glutMainLoop();


  return 0;
}
