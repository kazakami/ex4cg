#include <iostream>
#include <stdlib.h>
#include <random>

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

  //Zバッファを有効にする
  image.Enable(kazakami::Image3d::Function::DepthTest);
  //隠面処理を無効にする
  image.Disable(kazakami::Image3d::Function::CullFace);

  //三角形の描画数の指定。デフォルトは10
  int count = 10;
  //コマンドライン引数のファイル名の後に数字を指定することで
  //描画する三角形の個数を指定可能
  if (argc >= 2)
    count = std::stoi(argv[1]);

  //乱数の初期シード用
  std::random_device rd;
  //メルセンヌツイスターによる乱数生成
  std::mt19937 mt(rd());
  //x, y座標用の一様分布乱数生成器 範囲は [-200, 200]
  std::uniform_real_distribution<double> position_XY(-200.0, 200.0);
  //z座標用の一様分布乱数生成器 範囲は [200, 400]
  std::uniform_real_distribution<double> position_Z(200.0, 400.0);

  //適当なマテリアル
  kazakami::Material hoge =
  {
    {0.810, 0.893, 0.931},  //ambient RGB
    {0.114, 0.514, 0.1919}, //diffuse RGB
    {0.3, 0.3, 0.4},        //specular RGB
    24                      //shininess
  };

  //マテリアルを指定
  image.SetMaterial(hoge);


  for (int i = 0; i < count; i++)
  {
    //3点の座標をランダムに生成
    kazakami::Vector3d
      a(position_XY(mt),
	position_XY(mt),
	position_Z(mt)),
      b(position_XY(mt),
	position_XY(mt),
	position_Z(mt)),
      c(position_XY(mt),
	position_XY(mt),
	position_Z(mt));
    std::cerr << "Triangle" << i << "\n"
	      << "\t" << a << "\n"
	      << "\t" << b << "\n"
	      << "\t" << c << "\n";
    //三角形の描画
    image.drawTriangle(a, b, c);
  }
  
  //画像の出力
  image.writeOutWithPPM("hoge.ppm");

  //Bit Map 形式でも画像を出力可能
  //image.writeOut("hoge.bmp");
  return 0;
}
