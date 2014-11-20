#include <iostream>
#include <stdlib.h>
#include <random>
#include <memory>

#include "image.hpp"
#include "vrml.hpp"

//画像の大きさの指定
constexpr int sx = 256;
constexpr int sy = 256;

typedef unsigned char uchar;

int main(int argc, char* argv[])
{
  //キャンパスの生成。
  kazakami::Image3d image(sx, sy);
  //原点を中心に描画するためにずらす。
  image.SetOffset(sx/2,sy/2);
  //変換行列を初期化
  image.LoadIdentity();
  //焦点距離256の透視射影
  image.Perspective(256);
  //光源位置を(-1, -1, 2)に指定
  image.SetLightPos({-1, -1, 2});
  //光源色を(1, 1, 1)に指定
  image.SetLightCol(1, 1, 1);
  //カメラの状態を指定
  //eye, lookat, up
  image.LookAt({0, 0, 0}, {0, 0, 1}, {0, 1, 0});

  image.Scale(2, 2, 1);

  //Zバッファを有効にする
  image.Enable(kazakami::Image3d::Function::DepthTest);
  //隠面処理を無効にする
  image.Enable(kazakami::Image3d::Function::CullFace);


  if (argc >= 3)
  {
    auto tex = std::make_shared<kazakami::Texture>();
    if (tex->LoadPPM(argv[2]))
    {
      std::cerr << "Enable environment mapping" << std::endl;
      image.SetEnvMapTex(tex);
    }
    else
    {
      std::cerr << "invalid environment texture" << std::endl;
    }
  }

  if (argc >= 2)
  {
    kazakami::Polygon polygon;
    polygon.ReadVRML(argv[1]);
    polygon.Draw(&image);
  }

  //画像ファイルを出力する
  image.writeOutWithPPM("hoge.ppm");

  //.bmp 形式でも画像を出力可能
  image.writeOut("piyo.bmp");
  return 0;
}
