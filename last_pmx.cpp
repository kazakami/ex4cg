#include <GL/glut.h>
#include <float.h>


#include "vmdLoader.hpp"
#include "pmxLoader.hpp"
#include "image.hpp"
#include "imageLoader.hpp"


double ratio = 1;


int sx = 256 * ratio;
int sy = 256 * ratio;

bool depth = false;

typedef unsigned char uchar;
typedef unsigned int uint;

kazakami::Image3d im(sx, sy);
kazakami::PMXLoader pmx;
std::vector<std::shared_ptr<kazakami::Texture>> textures;
std::vector<kazakami::Material> materials;

int GLframe = 0; //フレーム数
int GLtimenow = 0;//経過時間
int GLtimebase = 0;//計測開始時間

double l = 0, ll = 0, size = 10.0; //視点とかに使う。
double f = 20; //焦点距離?
double fd = 20; //絞り?


/*------ render_viewer functions ------------------*/
static void redraw_image(void)
{
  im.ResetZBuffer();
  for (int i = 0; i < sy; i++)
    for (int j = 0; j < sx; j++)
    {
      im.at(i, j, kazakami::colorRGB::B) = i / ratio;
      im.at(i, j, kazakami::colorRGB::G) = j / ratio;
      im.at(i, j, kazakami::colorRGB::R) = 0;
    }

  im.LoadIdentity();

  //static double hogehoge = 0;
  //hogehoge += 0.01;
  im.Rotate(3.141592653589793238, 0, 1, 0);

  //im.Parallel(30);
  im.Perspective(400*ratio*(size/10));

  //im.Scale(size * ratio, size * ratio, size * ratio);
  im.Translate(0, -10, 0);
  /*
  if ((int)(l*10) % 3 == 0)
  {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_int_distribution<int> dice(-8,8);
    im.dif = dice(mt);
  }
  */

  im.LookAt(kazakami::Vector3d(40*sin(l)*cos(ll), 40*sin(ll), 40*cos(l)*cos(ll)),
	    kazakami::Vector3d(0, 0, 0),
	    kazakami::Vector3d(0, 1, 0));
  
  //床を描画
  //適当なマテリアル
  /*
  static const kazakami::Material hoge =
  {
    {1, 1, 1},  //ambient RGB
    {0, 0, 0}, //diffuse RGB
    {0, 0, 0},        //specular RGB
    50                      //shininess
  };
  //マテリアルを指定
  im.SetMaterial(hoge);
  static const kazakami::Vector3d groundPos[] =
    {
      {20, 0, 20},
      {20, 0, -20},
      {-20, 0, -20},
      {-20, 0, 20},
    };
  */
  /*
  im.Disable(kazakami::Image3d::Function::DepthTest);
  im.drawTriangle_constant(groundPos[0], groundPos[1], groundPos[2],
			   {0, 1, 0});
  im.drawTriangle_constant(groundPos[2], groundPos[3], groundPos[0],
			   {0, 1, 0});
  im.Enable(kazakami::Image3d::Function::DepthTest);
  */
  


  int drawnVertexNum = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  im.CulcVertex();


/*
#ifdef _OPENMP
#pragma omp for
#endif
  for (int i = 0; i < pmx.faceNum/3; i++)
  {
    im.drawTriangle_Phong(pmx.faceIndex[i*3],
			  pmx.faceIndex[i*3+1],
			  pmx.faceIndex[i*3+2]);
 }*/
 
  for (int mtr_i = 0; mtr_i < pmx.mtrNum; mtr_i++)
  {
#ifdef _OPENMP
#pragma omp single
#endif
    {
      im.SetMaterial(materials[mtr_i]);
      const kazakami::PMXMaterial & pm = pmx.materials[mtr_i];
      const int textureNumber = pm.normalTexture;
      if (textureNumber == -1)
      {
	im.SetTexture(nullptr);
      }
      else
      {
	im.SetTexture(textures[textureNumber]);
      }
      if (pm.drawFlag & 0x01)//両面描画か
	im.Disable(kazakami::Image3d::Function::CullFace);
      else
	im.Enable(kazakami::Image3d::Function::CullFace);
    }// end single
    //この材質に対応する面を描画
    const int loopCount = pmx.materials[mtr_i].vertexNum/3;
#ifdef _OPENMP
#pragma omp for
#endif
    for (int face_i = 0; face_i < loopCount; face_i++)
    {
      im.drawTriangle_Phong(pmx.faceIndex[drawnVertexNum+face_i*3],
			    pmx.faceIndex[drawnVertexNum+face_i*3+1],
			    pmx.faceIndex[drawnVertexNum+face_i*3+2]);
    }
#ifdef _OPENMP
#pragma omp single
#endif
    {
      drawnVertexNum += pmx.materials[mtr_i].vertexNum;
    }
  }
}//parallel

 if (depth)
 {
   static const uint kernel[5][5] =
   {
      {1, 4, 6, 4, 1},
      {4, 16, 24, 16, 4},
      {6, 24, 36, 24, 6},
      {4, 16, 24, 16, 4},
      {1, 4, 6, 4, 1}
    };
   std::vector<uchar> vignette; //BGR
   vignette.resize(sx * sy * 3);

#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < sx; i++)
       for (int j = 0; j < sy; j++)
       {
	 uint r = 0, g = 0, b = 0;
	 int div = 0;
	 for (int k = -2; k <= 2; k++)
	   for (int l = -2; l <= 2; l++)
	   {
	     if (i + k >= 0 && i + k < sx && j + l >= 0 && j + l < sy)
	     {
	       r += im.at(i+k, j+l, kazakami::colorRGB::R)
		 * kernel[k+2][l+2];
	       g += im.at(i+k, j+l, kazakami::colorRGB::G)
		 * kernel[k+2][l+2];
	       b += im.at(i+k, j+l, kazakami::colorRGB::B)
		 * kernel[k+2][l+2];
	       div += abs(kernel[k+2][l+2]);
	     }
	   }
	 vignette[3 * (i + j * sx) + 2] = r / div;
	 vignette[3 * (i + j * sx) + 1] = g / div;
	 vignette[3 * (i + j * sx) + 0] = b / div;
       }

  
     //深度で遊ぼう
#pragma omp for
     for (int i = 0; i < sx; i++)
       for (int j = 0; j < sy; j++)
       {
	 double z_val = im.z_get(i, j);
	 double val; //重み?
	 if (z_val == DBL_MAX)
	 {
	   val = 0;
	 }
	 else
	 {
	   val = 1.0 / (abs(z_val+f)/fd + 1);
	 }
	 im.at(i, j, kazakami::colorRGB::R) = 
	   im.at(i, j, kazakami::colorRGB::R) * val
	   + vignette[3 * (i + j * sx) + 2] * (1 - val);
	 im.at(i, j, kazakami::colorRGB::G) = 
	   im.at(i, j, kazakami::colorRGB::G) * val
	   + vignette[3 * (i + j * sx) + 1] * (1 - val);
	 im.at(i, j, kazakami::colorRGB::B) = 
	   im.at(i, j, kazakami::colorRGB::B) * val
	   + vignette[3 * (i + j * sx) + 0] * (1 - val);
       }
   }
 }

  //*/

  /* Just for Safe */
  glPushMatrix();{

    /* Just for Safe Again */
    glLoadIdentity();

    /* With Identity Projection Matrix,                                                                                                       
          Window Screen Range is (-1,-1)[bottomleft] to (1,1)[upright].                                                                       
          Write the Buffer from (-1,-1), bottom left. */
    glRasterPos3d(-1.0, -1.0, 0);
    
    glDrawPixels(sx, sy, GL_RGB, GL_UNSIGNED_BYTE,
		 (const GLvoid *)(im.ptr()));
                 //(const GLvoid *)(vignette.data()));
    
    glPopMatrix();
  }


  GLframe++; //フレーム数を＋１
  GLtimenow = glutGet(GLUT_ELAPSED_TIME);//経過時間を取得
        
  if (GLtimenow - GLtimebase > 1000)      //１秒以上たったらfpsを出力
  {
    std::cout << "\rfps:" << GLframe*1000.0/(GLtimenow-GLtimebase) 
	      << "   f = " << f << std::flush;
    GLtimebase = GLtimenow;//基準時間を設定                
    GLframe = 0;//フレーム数をリセット
  }


  /* Double Buffering */
  glutSwapBuffers();
}

bool right = false, left = false;


void mouse(int button, int state, int x, int y)
{
  if (button == GLUT_RIGHT_BUTTON)
  {
    if (state == GLUT_DOWN)
      right = true;
    else if (state == GLUT_UP)
      right = false;
  }
  else if (button == GLUT_LEFT_BUTTON)
  {
    if (state == GLUT_DOWN)
      left = true;
    else if (state == GLUT_UP)
      left = false;
  }
  //std::cout << "mouse" << std::endl;
}

void motion(int x, int y)
{
  //std::cout << x << ":" << y << std::endl;
  if (left)
  {
    l = 2 * 3.141592653 * (1.0 * x / sx);
    ll = 3.141592653 * (1.0 * y / sy - 0.5);
  }
  if (right)
  {
    f = y*2-256;
    fd = abs(x) / 2.0;
  }
}

void keyBoard(unsigned char key, int x, int y)
{
  if (key == 'w')
    size = 20.0;
  else if (key == 's')
    size = 10.0;
  if (key == 'q')
    f += 1;
  if (key == 'e')
    f -= 1;
}

void errExit(const std::string & str, int e = 1)
{
  std::cerr << str << std::endl;
  exit(e);
}


int main(int argc, char* argv[])
{  //読み込むVRMLファイル名
  std::string filename;
  bool filenameIs = false;
  //モーションデータ
  std::string motionFilename;
  bool motionFilenameIs = false;
  int i = 1;
  while (i < argc)
  {
    //std::cout << i << " : " << argv[i] << std::endl;
    //扱い易くするためにstd::stringにする。
    std::string str(argv[i]);
    if (str == "-m")
    {
      if (i+1 < argc)
	motionFilename = argv[i+1];
      else
	errExit("missing filename after '-m'");
      motionFilenameIs = true;
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
    else if (str == "--depth")
    {
      depth = true;
      i++;
    }
    else
    {
      filename = argv[i];
      filenameIs = true;
      i++;
    }
  }
  sx = 256*ratio;
  sy = 256*ratio;
  im.Resize(sx, sy);
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

#ifdef _WIN32
  const std::string directory
    = filename.substr(0, filename.find_last_of("\\") + 1);
#else
  const std::string directory
    = filename.substr(0, filename.find_last_of("/") + 1);
#endif

  
  //モーション読み込み。
  if (motionFilenameIs)
  {
    auto vmd = std::make_shared<kazakami::vmdLoader>();
    if (!vmd->Load(motionFilename))
      std::cerr << "no such a file : " << motionFilename << std::endl;
    else
      im.SetMotion(vmd);
  }


  pmx.readPMX(filename);
  
  //頂点を登録
  for (int i = 0; i < pmx.vertexNum; i++)
  {
    const kazakami::PMXVertex & v = pmx.vertexes[i];
	/*
    im.AddVertex({v.pos[0], v.pos[1], v.pos[2]},
		 {v.normal[0], v.normal[1], v.normal[2]},
		 {v.UV[0], v.UV[1]});
		 */
    im.AddVertex(kazakami::Vector3d(v.pos[0], v.pos[1], v.pos[2]),
		 kazakami::Vector3d(v.normal[0], v.normal[1], v.normal[2]),
		 kazakami::Vector2d(v.UV[0], v.UV[1]),
		 v.weightType,
		 v.boneIndex1, v.boneIndex2, v.boneIndex3, v.boneIndex4,
		 v.weight1, v.weight2, v.weight3, v.weight4);
  }

  //ボーンの登録
  for (int i = 0; i < pmx.boneNum; i++)
  {
    const kazakami::PMXBone & b = pmx.bones[i];
    im.AddBone(b.name, b.nameEng,
	       b.pos[0], b.pos[1], b.pos[2],
	       b.parentBoneIndex);
  }
  //ボーンの初期化
  im.InitBone();
  
  //テクスチャの読み込み
  for (int i = 0; i < pmx.texNum; i++)
  {
    const std::string & str = pmx.textureNames[i];
    kazakami::imageLoader pic;
    pic.load(directory + str);
    auto tex = std::make_shared<kazakami::Texture>();
    tex->LoadFromVector(pic.BGR_data(), pic.Width(), pic.Height());
    textures.push_back(tex);
  }

  //材質の読み込み
  for (int i = 0; i < pmx.mtrNum; i++)
  {
    const kazakami::PMXMaterial & mtr = pmx.materials[i];
    kazakami::Material m =
    {
      {mtr.ambient[0], mtr.ambient[1], mtr.ambient[2]},
      {mtr.diffuse[0], mtr.diffuse[1], mtr.diffuse[2]},
      {mtr.specular[0], mtr.specular[1], mtr.specular[2]},
      mtr.shininess
    };
    materials.push_back(std::move(m));
  }

  im.SetOffset(sx/2,sy/2);
  im.LoadIdentity();
  //im.SetLightPos({-500, -500, 500});
  im.SetLightPos(kazakami::Vector3d(-500, -500, -500));
  im.SetLightCol(1, 1, 1);
  //Zバッファを有効にする
  im.Enable(kazakami::Image3d::Function::DepthTest);
  //隠面処理をdisableにする
  //im.Disable(kazakami::Image3d::Function::CullFace);


#ifdef _OPENMP
  puts("MP is enable");
#else
  puts("MP is disable");
#endif
  

  /* Prepareing the GL and GLUT */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
  glutInitWindowSize(sx, sy);
  glutCreateWindow("LE4 CG");
  glutDisplayFunc(redraw_image);
  glutIdleFunc(redraw_image);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyBoard);
  glutMouseFunc(mouse);

  /* MainLoop */
  glutMainLoop();


  return 0;
}
