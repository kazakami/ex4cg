#include <fstream>
#include <iostream>
#include <vector>
#include <memory>

#include "Vector.hpp"
#include "image.hpp"

namespace kazakami
{
	typedef unsigned int uint;
  /*
class Triangle
{
  Vector3d pos_a, pos_b, pos_c;
  Vector3d normal_a, normal_b, normal_c;
  Material mat;
public:
  Triangle(const Vector3d & pa, const Vector3d & na,
	   const Vector3d & pb, const Vector3d & nb,
	   const Vector3d & pc, const Vector3d & nc);
  void Draw(Image3d * im, const Material & mtr);
};*/

class Polygon
{
  std::vector<Vector3d> pos;
  std::vector<Vector2d> posTex;
  std::vector<Vector3d> pos_normal;
  //face number which contains point
  std::vector<std::vector<int>> pos_belong_to;
  bool isTexed;
  struct face
  {
    uint a, b, c;
    face(uint _a, uint _b, uint _c)
    {
      a = _a;
      b = _b;
      c = _c;
    }
  };
  std::vector<face> faces;
  std::vector<face> facesTex;
  std::vector<Vector3d> face_normal;
  Material mtr;
  void AddPos(double x, double y, double z);
  void AddFace(uint a, uint b, uint c);
  void AddPosTex(double x, double y);
  void AddFaceTex(uint a, uint b, uint c);
  bool CheckIndex() const;
public:
  std::shared_ptr<Texture> tex;
  enum class ShadingType
  {
    Constant,
    Gouraud,
    Phong,
  };
  ShadingType shadingType;
  Polygon();
  //Multi Thread
  //void ReadVRML_MT(const char * filename);
  void ReadVRML(const char * filename);
  void AddVertex(Image3d * im);
  void SuperDraw(Image3d * im);
  void Draw(Image3d * im);
};





}
