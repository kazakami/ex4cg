#include <iostream>

#include <stdlib.h>

#include "image.hpp"
#include "vrml.hpp"

constexpr int sx = 1024;
constexpr int sy = 1024;

typedef unsigned char uchar;

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "please filenameeeeeee!!" << std::endl;
    return 1;
  }
  kazakami::Image3d im(sx, sy);
  im.SetOffset(sx/2,sy/2);
  im.LoadIdentity();
  im.Perspective(1024);
  im.SetLightPos({-1, -1, 2});
  im.SetLightCol(1, 1, 1);
  im.Scale(1, 1, 1);
  im.Translate(0, 0, 0);
  for (int i = 0; i < sy; i++)
    for (int j = 0; j < sx; j++)
    {
      im.at(i, j, kazakami::colorRGB::B) = i/4;
      im.at(i, j, kazakami::colorRGB::G) = j/4;
    }
  
#ifdef _OPENMP
  puts("MP is enable");
#else
  puts("MP is disable");
#endif
  
  kazakami::Polygon test;
  test.ReadVRML(argv[1]);

  int repeat = 1;
  if (argc >= 3)
    repeat = std::stoi(argv[2]);
  for (int i = 0; i < repeat; i++)
    test.Draw(&im);
  
  /*
  kazakami::Material ruby =
  {
    {0.1745,   0.01175,  0.01175},
    {0.61424,  0.04136,  0.04136},
    {0.727811, 0.626959, 0.626959},
    76.8
  };
  kazakami::Material white =
  {
    {1, 1, 1},
    {1, 1, 1},
    {1, 1, 1},
    10
  };
  */
  std::cout << im.writeOut("fuga.bmp") << std::endl;
  //im.writeOutWithPPM("hoge.ppm");
  return 0;
}
