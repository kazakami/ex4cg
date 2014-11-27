#include <iostream>
#include <memory>
#include <random>

#include <stdlib.h>
#include <GL/glut.h>

#include "image.hpp"
#include "vrml.hpp"

double ratio = 1;

int sx = 256;
int sy = 256;

int tmp_sx = sx;
int tmp_sy = sy;

typedef unsigned char uchar;

kazakami::Image3d image(sx, sy);
kazakami::Polygon test;

int GLframe = 0; //フレーム数
int GLtimenow = 0;//経過時間
int GLtimebase=0;//計測開始時間



/*------ render_viewer functions ------------------*/
static void redraw_image(void)
{
  //オブジェクトの回転角度
  static double l = 0;
  //Zバッファの初期化
  image.ResetZBuffer();
  //スクリーンの初期化　黒で塗りつぶす。
  for (int i = 0; i < sy; i++)
    for (int j = 0; j < sx; j++)
    {
      image.at(i, j, kazakami::colorRGB::B) = 0;
      image.at(i, j, kazakami::colorRGB::G) = 0;
      image.at(i, j, kazakami::colorRGB::R) = 0;
    }

  //変換行列を初期化
  image.LoadIdentity();
  //焦点距離256の透視投影
  image.Perspective(256);

  //(0, 0, 1)を軸としてlラジアンだけ回転
  image.Rotate(l, 0, 0, 1);
  image.Scale(2*ratio, 2*ratio, 1);
  //カメラの位置の設定
  //カメラは座標(0, 0, 0)にあり、座標(0, 0, 1)に向いていて、(0, 1, 0)方向を上向きとしている。
  image.LookAt({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
  test.Draw(&image);
  
  //オブジェクトの回転角度
  l += 0.01;


  /* Just for Safe */
  glPushMatrix();{

    /* Just for Safe Again */
    glLoadIdentity();

    glRasterPos3d(-1.0, -1.0, 0);
    
    glDrawPixels(sx, sy, GL_RGB, GL_UNSIGNED_BYTE,
                 (const GLvoid *)(image.ptr()));
    
    glPopMatrix();
  }


  GLframe++;
  GLtimenow = glutGet(GLUT_ELAPSED_TIME);//経過時間を取得
  
  //１秒以上たったらfpsを出力
  if (GLtimenow - GLtimebase > 1000)
  {
    std::cout << "\rfps:" 
	      << GLframe*1000.0/(GLtimenow-GLtimebase)
	      << std::flush;
    GLtimebase = GLtimenow;//基準時間を設定
    GLframe = 0;//フレーム数をリセット
  }

  //sx = tmp_sx;
  //sy = tmp_sy;
  //image.Resize(sx, sy);
  //image.SetOffset(sx/2,sy/2);

  /* Double Buffering */
  glutSwapBuffers();
}

void errExit(const std::string & str, int e = 1)
{
  std::cerr << str << std::endl;
  exit(e);
}

void Resize(int w, int h)
{
  tmp_sx = w;
  tmp_sy = h;
  glutReshapeWindow(sx, sy);
}

int main(int argc, char* argv[])
{
  //読み込むVRMLファイル名
  std::string filename;
  bool filenameIs = false;
  //環境マッピング用のテクスチャファイル名
  std::string envFilename;
  bool envFilenameIs = false;
  //出力するファイル名
  std::string outputFilename = "hoge";
  bool outputFilenameIs = false;
  //どのように出力するか
  enum class OutputType
  {
    WINDOW,
    PPM,
    BMP,
  }outputType = OutputType::WINDOW;
  int i = 1;
  while (i < argc)
  {
    //std::cout << i << " : " << argv[i] << std::endl;
    //扱い易くするためにstd::stringにする。
    std::string str(argv[i]);
    if (str == "-e")
    {
      if (i+1 < argc)
	envFilename = argv[i+1];
      else
	errExit("missing filename after '-e'");
      envFilenameIs = true;
      i+=2;
    }
    else if (str == "-o")
    {
      if (i+1 < argc)
	outputFilename = argv[i+1];
      else
	errExit("missing filename after '-o'");
      outputFilenameIs = true;
      i+=2;
    }
    else if (str == "-s")
    {
      if (i+1 < argc)
	ratio = std::stod(argv[i+1]);
      else
	errExit("");
      i+=2;
    }
    else if (str == "--bmp")
    {
      outputType = OutputType::BMP;
      i++;
    }
    else if (str == "--ppm")
    {
      outputType = OutputType::PPM;
      i++;
    }
    else
    {
      filename = argv[i];
      filenameIs = true;
      i++;
    }
  }
  tmp_sx = sx = 256*ratio;
  tmp_sy = sy = 256*ratio;
  image.Resize(sx, sy);
  /*
  if (filenameIs)
    std::cout << "filename: " << filename << std::endl;
  if (envFilenameIs)
    std::cout << "env: " << envFilename << std::endl;
  if (outputFilenameIs)
    std::cout << "out: " << outputFilename << std::endl;
  std::cout << "screen size: " << sx << "x" << sy << std::endl;
  */

  if (!filenameIs)
  {
    errExit("no input file");
  }
  if (envFilenameIs)
  {
    auto tex = std::make_shared<kazakami::Texture>();
    if (tex->LoadPPM(envFilename))
    {
      std::cerr << "Enable environment mapping" << std::endl;
      image.SetEnvMapTex(tex);
    }
    else
    {
      std::cerr << "invalid environment texture." << std::endl;
    }
  }
  image.SetOffset(sx/2,sy/2);
  image.LoadIdentity();
  image.SetLightPos({-500, -500, 500});
  image.SetLightCol(1, 1, 1);
  //Zバッファを有効にする
  image.Enable(kazakami::Image3d::Function::DepthTest);
  //隠面処理を無効にする
  image.Disable(kazakami::Image3d::Function::CullFace);
#ifdef _OPENMP
  puts("MP is enable");
#else
  puts("MP is disable");
#endif
  
  test.ReadVRML(filename);
  test.AddVertex(&image);


  if (outputType ==  OutputType::BMP
      || outputType ==  OutputType::PPM)
  {
    image.LoadIdentity();
    //焦点距離256の透視投影
    image.Perspective(256);
    
    image.Scale(2*ratio, 2*ratio, 1);
    //カメラの位置の設定
    //カメラは座標(0, 0, 0)にあり、座標(0, 0, 1)に向いていて、(0, 1, 0)方向を上向きとしている。
    image.LookAt({0, 0, 0}, {0, 0, 1}, {0, 1, 0});
    test.Draw(&image);
  }

  if (outputType ==  OutputType::BMP)
  {
    image.writeOut(outputFilename);
    return 0;
  }

  if (outputType ==  OutputType::PPM)
  {
    image.writeOutWithPPM(outputFilename);
    return 0;
  }


  /* Prepareing the GL and GLUT */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
  glutInitWindowSize(sx, sy);
  glutCreateWindow("LE4 CG");
  glutReshapeFunc(Resize);
  glutDisplayFunc(redraw_image);
  glutIdleFunc(redraw_image);

  /* MainLoop */
  glutMainLoop();


  return 0;
}
