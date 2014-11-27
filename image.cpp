#include <string>
#include <sstream>
#include <float.h>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "image.hpp"



namespace kazakami
{


void Texture::LoadFromVector(const std::vector<uchar> & d, size_t x, size_t y)
{
  //std::cerr << "it`s omoi" << std::endl;
  std::copy(d.begin(), d.end(), std::back_inserter(data));
  sx = x;
  sy = y;
}

void Texture::LoadFromVector(std::vector<uchar> && d, size_t x, size_t y)
{
  std::swap(data, d);
  sx = x;
  sy = y;
}

bool Texture::LoadPPM(const std::string & filename)
{
  return LoadPPM(filename.c_str());
}

bool Texture::LoadPPM(const char * filename)
{
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);

  if (!ifs)
    return false;

  std::string str;
  //種類の取得。今はP6のみ対応
  if (getline(ifs, str))
  {
    //std::cout << str << std::endl;
    if (str != "P6")
      return false;
  }
  //サイズの取得
  if (getline(ifs, str))
  {
    //std::cout << str << std::endl;
    std::size_t pos = str.find(' ');
    if (pos !=std::string::npos)
    {
      sx = std::stoi(str.substr(0, pos));
      sy = std::stoi(str.substr(pos+1));
      data.reserve(sx * sy * 3);
    }
    else
      return false;
  }
  //最大値の取得
  if (getline(ifs, str))
  {
    //割とどうでもいい
    //std::cout << str << std::endl;
  }

  uchar u;
  //ifs.read((char*)&u, sizeof(uchar));
  while (!ifs.eof())
  {
    ifs.read((char*)&u, sizeof(uchar));
    data.push_back(u);
  }
  data.pop_back();
  std::reverse(data.begin(), data.end());
  //std::cout << data.size() << std::endl;

  return true;
}



uchar * Texture::ptr()
{
  return data.data();
}

const uchar * Texture::ptr() const
{
  return data.data();
}

size_t Texture::height() const
{
  return sy;
}

size_t Texture::width() const
{
  return sx;
}

Vector3d Texture::GetColor(double x, double y) const
{
  x = (x + 1) / 2;
  y = (y + 1) / 2;
  double px_dbl = x * sx;
  double py_dbl = y * sy;
  int px = static_cast<int>(px_dbl);
  int py = static_cast<int>(py_dbl);
  double dx = px_dbl - px;
  double dy = py_dbl - py;
  if (px < 0 || py < 0 || px >= static_cast<int>(sx) || py >= static_cast<int>(sy))
      return Vector3d(0, 0, 0);
  Vector3d h00(data.at(3 * (px + sx * py) + 2),
	       data.at(3 * (px + sx * py) + 1),
	       data.at(3 * (px + sx * py) + 0));
  if (px == static_cast<int>(sx) - 1 || py == static_cast<int>(sy) - 1)
    return h00 * (1.0 / 255);
  Vector3d h01(data.at(3 * (px + sx * (py + 1)) + 2),
	       data.at(3 * (px + sx * (py + 1)) + 1),
	       data.at(3 * (px + sx * (py + 1)) + 0));
  Vector3d h10(data.at(3 * (px + 1 + sx * py) + 2),
	       data.at(3 * (px + 1 + sx * py) + 1),
	       data.at(3 * (px + 1 + sx * py) + 0));
  Vector3d h11(data.at(3 * (px + 1 + sx * (py + 1)) + 2),
	       data.at(3 * (px + 1 + sx * (py + 1)) + 1),
	       data.at(3 * (px + 1 + sx * (py + 1)) + 0));
  return blend((1 - dx) * (1 - dy) * (1.0 / 255), h00,
	             dx * (1 - dy) * (1.0 / 255), h10,
	       (1 - dx) *       dy * (1.0 / 255), h01,
	             dx *       dy * (1.0 / 255), h11);
}



Image::Image(size_t x, size_t y)
  :sx(x), sy(y)
{
  data.resize(sx * sy * 3);
}

void Image::resize(size_t x, size_t y)
{
  data.resize(x * y * 3);
  sx = x;
  sy = y;
}

uchar & Image::at(size_t x, size_t y, colorRGB col)
{
  return data.at(3 * (x + sx * y) + static_cast<int>(col));
}

const uchar & Image::at(size_t x, size_t y, colorRGB col) const
{
  return data.at(3 * (x + sx * y) + static_cast<int>(col));
}

uchar * Image::ptr()
{
  return data.data();
}

const uchar * Image::ptr() const
{
  return data.data();
}

size_t Image::height() const
{
  return sy;
}

size_t Image::width() const
{
  return sx;
}

//bmp形式で画像を出力する。
int Image::writeOut(const char * filename) const
{
  std::ofstream fout;
  fout.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);

  if (!fout)
  {
    fout.close();
    return -1;
  }

  char s[2];
  s[0] = 'B';
  s[1] = 'M';
  short reserved = 0;
  long offset = 14 + 40;
  long size = sx * sy + offset;

  fout.write(s, 2);
  fout.write((char*)&size, 4);
  fout.write((char*)&reserved, 2);
  fout.write((char*)&reserved, 2);
  fout.write((char*)&offset, 4);

  long var_long;
  short var_short;

  var_long = 40;
  fout.write((char*)&var_long, 4);
  var_long = (long)sx;
  fout.write((char*)&var_long, 4);
  var_long = -(long)sy;
  fout.write((char*)&var_long, 4);

  var_short = 1;
  fout.write((char*)&var_short, 2);
  //biBitCount
  var_short = 32;
  fout.write((char*)&var_short, 2);
  //biCOmpression
  var_long = 0;
  fout.write((char*)&var_long, 4);
  //biSizeImage
  var_long = 3780;
  fout.write((char*)&var_long, 4);
  //biXPixPerMeter
  fout.write((char*)&var_long, 4);
  //biYPixPerMeter
  fout.write((char*)&var_long, 4);
  //biClrUsed
  var_long = 0;
  fout.write((char*)&var_long, 4);
  //biCirImportant
  fout.write((char*)&var_long, 4);

  //char var_byte = 255;
  char var_byte_zero = 0;
  int end = static_cast<int>(sx * sy);
  for (int i = 0; i < end; i++)
  {
    fout.write((const char *)&data.at(3 * i + static_cast<int>(colorRGB::B)), 1);
    fout.write((const char *)&data.at(3 * i + static_cast<int>(colorRGB::G)), 1);
    fout.write((const char *)&data.at(3 * i + static_cast<int>(colorRGB::R)), 1);
    fout.write((char*)&var_byte_zero, 1);
  }

  fout.close();
  return size;
}

int Image::writeOutWithPPM(const char * filename) const
{
  std::ofstream ofs(filename);

  ofs << "P3" << std::endl;
  ofs << sx << "\t" << sy << std::endl;
  ofs << 255 << std::endl;
  for (int j = 0; j < (int)sy; j++)
  {
    for ( int i = 0; i < (int)sx ; i++)
    {
      ofs << " " << std::to_string(at(i, j, colorRGB::R))
	  << " " << std::to_string(at(i, j, colorRGB::G))
	  << " " << std::to_string(at(i, j, colorRGB::B));
    }
    ofs << std::endl;
  }
  return 0;
}




Image3d::Image3d(size_t x, size_t y)
  :image(x, y),
   enableZBuffer(false),
   enableCulling(false),
   cam_pos(0, 0, 0),
   isParallelLight(true),
   offsetX(0), offsetY(0),
   worldMatrix(IdentityMat()),
   viewMatrix(IdentityMat()),
   projMatrix(IdentityMat()),
   vMatrix(IdentityMat()),
   covMatrix(IdentityMat()),
   envTex(nullptr),
   texture(nullptr),
   vmdData(nullptr),
   motionFrame(0)
{
  zBuffer.resize(x * y);
  int s = static_cast<int>(x * y);
  for (int i = 0; i < s; i++)
    zBuffer.at(i) = DBL_MAX;
}

void Image3d::Resize(size_t x, size_t y)
{
  image.resize(x, y);
  zBuffer.resize(x * y);
}

size_t Image3d::height() const
{
  return image.height();
}

size_t Image3d::width() const
{
  return image.width();
}

void Image3d::PushMatrix()
{
  worldMatrixStack.push(worldMatrix);
  viewMatrixStack.push(viewMatrix);
  projMatrixStack.push(projMatrix);
  vMatrixStack.push(vMatrix);
  covMatrixStack.push(covMatrix);
}
void Image3d::PopMatrix()
{
  worldMatrix = worldMatrixStack.top();
  viewMatrix = viewMatrixStack.top();
  projMatrix = projMatrixStack.top();
  vMatrix = vMatrixStack.top();
  covMatrix = covMatrixStack.top();
  
  worldMatrixStack.pop();
  viewMatrixStack.pop();
  projMatrixStack.pop();
  vMatrixStack.pop();
  covMatrixStack.pop();
}

void Image3d::AddVertex(const Vector3d & pos,
			const Vector3d & normal,
			const Vector2d & texPos)
{
  Vertex v;
  v.pos = pos;
  v.normal = normal;
  v.texPos = texPos;
  v.normalPlusPos = pos + normal;
  vertexes.push_back(v);
}

void Image3d::AddVertex(const Vector3d & pos,
			const Vector3d & normal,
			const Vector2d & texPos,
			unsigned char weightType,
			int boneIndex1,
			int boneIndex2,
			int boneIndex3,
			int boneIndex4,
			float weight1,
			float weight2,
			float weight3,
			float weight4)
{
  Vertex v;
  v.pos = pos;
  v.normal = normal;
  v.texPos = texPos;
  v.normalPlusPos = pos + normal;
  v.weightType = weightType;
  v.boneIndex1 = boneIndex1;
  v.boneIndex2 = boneIndex2;
  v.boneIndex3 = boneIndex3;
  v.boneIndex4 = boneIndex4;
  v.weight1 = weight1;
  v.weight2 = weight2;
  v.weight3 = weight3;
  v.weight4 = weight4;
  //emplace_backにした方がいい?
  vertexes.push_back(v);
}

void Image3d::AddBone(std::string name,
		      std::string nameEng,
		      float X,
		      float Y,
		      float Z,
		      int parentIndex)
{
  /*
  Bone b = 
  {
    name,
    nameEng,
    {X, Y, Z},
    {X, Y, Z},
    parentIndex,
    IdentityMat(),
    IdentityMat(),
    IdentityMat(),
    IdentityQuaternion(),
    std::vector<int>()
  };
  //*/
  //*
  Bone b(name,
	 nameEng,
	 X, Y, Z,
	 parentIndex);
  //*/
  //std::stringあるしmoveにする意味あるよね
  //そもそもemplace_backにすべき?
  //頻繁に呼ばれる処理じゃないし、まぁいけるやろ
  bones.push_back(std::move(b));
}

void Image3d::InitBone()
{
  const uint boneNum = bones.size();
  //std::cout << boneNum << std::endl;
  for (uint i = 0; i < boneNum; i++)
  {
    //初期姿勢行列の計算
    bones[i].initMat = TranslateMat(bones[i].pos[0],
				    bones[i].pos[1],
				    bones[i].pos[2]);
    bones[i].offsetMat = TranslateMat(-bones[i].pos[0],
				      -bones[i].pos[1],
				      -bones[i].pos[2]);
    //親子関係の構築　
    const int parentIndex = bones[i].parentBoneIndex;
    if (parentIndex == -1)
      rootBoneIndexes.push_back(i);
    else
      bones[parentIndex].children.push_back(i);
    //マップへ登録
    boneNameMap.insert(make_pair(bones[i].name, i));
  }
}


void Image3d::SetMotion(std::shared_ptr<vmdLoader> vmd)
{
  vmdData = vmd;
}

void Image3d::MoveBone()
{
  /*
  {//test
    static double i = 0;
    i += 0.1;
    bones[9].rotation.t = sin(i) / 4;
    bones[48].rotation.t = - sin(i) / 4;
    bones[17].rotation.t = -(1 - sin(i)) / 4;
    bones[55].rotation.t = (1 - sin(i)) / 4;
    bones[84].rotation.v.Set(1, 0, 0);
    bones[84].rotation.t = sin(i) / 2;
    bones[85].rotation.v.Set(1, 0, 0);
    bones[85].rotation.t=  -(1 - sin(i)) / 3;
    bones[45].rotation.v.Set(1, 0, 0);
    bones[45].rotation.t = sin(-i) / 2;
    bones[46].rotation.v.Set(1, 0, 0);
    bones[46].rotation.t = -(1 - sin(-i)) / 3;
    return;
  }
  //*/
  //*
  if (vmdData)
  {
    for (const auto & b : vmdData->boneMotion)
    {
      const auto it = boneNameMap.find(b.GetName());
      //モーションデータが指定するボーンを全て保持しているとは限らない
      if (it != boneNameMap.end()) //実在ボーン
      {
	const int index = it->second;
	auto rotation = b.GetRotation(motionFrame);
	//std::cout << rotation.toMat() << std::endl;
	bones[index].rotation = rotation;
	Vector3d v = b.GetPos(motionFrame);
        bones[index].pos[0] = bones[index].initPos[0] + v.X();
        bones[index].pos[1] = bones[index].initPos[1] + v.Y();
        bones[index].pos[2] = bones[index].initPos[2] + v.Z();
      }
    }
    motionFrame++;
    //motionFrame=0;
    if (vmdData->GetMaxFrame() + 10 < motionFrame)
    {
      std::cout << "Note: motion loop restart at frame "
		<< motionFrame
		<< std::endl;
      motionFrame = 0;
    }
    //std::cout << motionFrame++ << std::endl;
  }
  //*/
}

//ボーン動かしたい!!
void Image3d::CulcBone()
{
  const uint rootNum = rootBoneIndexes.size();
#ifdef _OPENMP
#pragma omp for
#endif
  for (uint index = 0; index < rootNum; index++)
  {
    Bone & r = bones[rootBoneIndexes[index]];
    r.boneMat = TranslateMat(r.pos[0], r.pos[1], r.pos[2])
      * r.offsetMat;
    CulcBone(rootBoneIndexes[index]);
  }
}

void Image3d::CulcBone(uint index)
{
  Bone & b = bones[index];
  for (int i : b.children)
  {
    Bone & c = bones[i]; //children
    //ボーンが保持すべき変換の順番
    //行列の変換は右から行われる事に注意。
    //----------------------------
    //自身のオフセット
    //自身の回転
    //自身の現在位置への移動
    //親のオフセット       |
    //親の回転             +---親の保持する変換行列
    //親の現在位置への移動 |
    //親の親の・・・       |
    //　　・               |
    //　　・               |
    //　　・               |
    //----------------------------
    c.boneMat = b.boneMat                           //親の変行列
      * TranslateMat(c.pos[0], c.pos[1], c.pos[2])  //自身の現在位置
      * c.rotation.toMat()                          //自身の回転
      * c.offsetMat;                                //自身のオフセット
    CulcBone(i);
  }
}

void Image3d::CulcVertex()
{
#ifdef _OPENMP
#pragma omp single
#endif
  {
    MoveBone();
  }//end single
  CulcBone();
  /*
  //ボーン上のローカル座標から世界座標への変換
  std::vector<Matrix4x4d> worldBones;
  //ボーン上のローカル座標からスクリーン座標への変換
  //別途スクリーンオフセットを足す必要がある。
  std::vector<Matrix4x4d> convBones;
  int boneNum = bones.size();
  worldBones.resize(boneNum);
  convBones.resize(boneNum);
#ifdef _OPENMP
#pragma omp for
#endif
  for (int i = 0; i < boneNum; i++)
  {
    worldBones[i] = worldMatrix * bones[i].boneMat;
    convBones[i] = covMatrix * bones[i].boneMat;
  }
  */
  int num = vertexes.size();
#ifdef _OPENMP
#pragma omp for
#endif
  for (int i = 0; i < num; i++)
  {
    Vertex & v = vertexes[i];
    Matrix4x4d boneMat;
    if (v.weightType == 0)//BDEF1
    {
      boneMat = bones[v.boneIndex1].boneMat;
    }
    else if (v.weightType == 1 || v.weightType == 3)//BDEF2 or SDEF
    {
      //SDEFがよく分からないのでBDEF2として変形する。
      boneMat = v.weight1 * bones[v.boneIndex1].boneMat
	+ (1 - v.weight2) * bones[v.boneIndex2].boneMat;
    }
    else if (v.weightType == 2)//BDEF4
    {
      boneMat = v.weight1 * bones[v.boneIndex1].boneMat
	+ v.weight2 * bones[v.boneIndex2].boneMat
	+ v.weight3 * bones[v.boneIndex3].boneMat
	+ v.weight4 * bones[v.boneIndex4].boneMat;
    }
    //*
    Vector3d lpos = matConv(boneMat, v.pos);
    Vector3d lnormalPlusPos = matConv(boneMat, v.normalPlusPos);
    v.worldPos = worldConv(lpos);
    v.convedPos = conversion(lpos);
    v.worldNormal = worldConv(lnormalPlusPos) - v.worldPos;
    v.worldNormal.Normalise();
    //*/
    /*
    v.worldPos = worldConv(v.pos);
    //v.viewPos = viewConv(v.pos);
    v.convedPos = conversion(v.pos);
    v.worldNormal = worldConv(v.normalPlusPos) - v.worldPos;
    v.worldNormal.Normalise();
    //*/
  }
}

void Image3d::SetEnvMapTex(std::shared_ptr<Texture> tex)
{
  envTex = tex;
}

void Image3d::SetTexture(std::shared_ptr<Texture> tex)
{
  texture = tex;
}

void Image3d::Enable(Function f)
{
  switch (f)
  {
  case Function::DepthTest:
    enableZBuffer = true;
    break;
  case Function::CullFace:
    enableCulling = true;
    break;
  }
}

void Image3d::Disable(Function f)
{
  switch (f)
  {
  case Function::DepthTest:
    enableZBuffer = false;
    break;
  case Function::CullFace:
    enableCulling = false;
    break;
  }
}

uchar & Image3d::at(size_t x, size_t y, colorRGB col)
{
  return image.at(x, y, col);
}
const uchar & Image3d::at(size_t x, size_t y, colorRGB col) const
{
  return image.at(x, y, col);
}

uchar * Image3d::ptr()
{
  return image.ptr();
}

const uchar * Image3d::ptr() const
{
  return image.ptr();
}

void Image3d::ResetZBuffer()
{
  std::fill(zBuffer.begin(), zBuffer.end(), DBL_MAX);
  /*
  int range = image.height() * image.width();
  for (int i = 0; i < range; i++)
  {
    zBuffer.at(i) = DBL_MAX;
  }
  */
}

double Image3d::z_get(size_t x, size_t y)
{
  return zBuffer.at(x + width() * y);
}

double & Image3d::z_at(size_t x, size_t y)
{
  return zBuffer.at(x + width() * y);
}

const double & Image3d::z_at(size_t x, size_t y) const
{
  return zBuffer.at(x + width() * y);
}

  int Image3d::writeOut(const std::string & filename) const
{
  return image.writeOut(filename.c_str());
}

  int Image3d::writeOutWithPPM(const std::string & filename) const
{
  return image.writeOutWithPPM(filename.c_str());
}

int Image3d::writeOut(const char * filename) const
{
  return image.writeOut(filename);
}

int Image3d::writeOutWithPPM(const char * filename) const
{
  return image.writeOutWithPPM(filename);
}

void Image3d::LoadIdentity()
{
  worldMatrix = IdentityMat();
  viewMatrix = IdentityMat();
  projMatrix = IdentityMat();
  vMatrix = IdentityMat();
  covMatrix = IdentityMat();
}

void Image3d::Parallel(double z)
{
  projMatrix = projMatrix *  ParallelProj(z);
  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

void Image3d::Perspective(double z)
{
  projMatrix = projMatrix * PerspectiveProj(z);
  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

void Image3d::Translate(double x, double y, double z)
{
  worldMatrix = worldMatrix * TranslateMat(x, y, z);
  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

void Image3d::Scale(double x, double y, double z)
{
  worldMatrix = worldMatrix * ScaleMat(x, y, z);
  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

void Image3d::LookAt(const Vector3d & eye,
		     const Vector3d & lookat,
		     const Vector3d & up)
{
  Vector3d ez = eye - lookat;
  ez.Normalise();
  Vector3d ex = outerProduct(up, ez);
  ex.Normalise();
  Vector3d ey = outerProduct(ex, ez);
  ey.Normalise();

  viewMatrix = Matrix4x4d(ex.X(), ex.Y(), ex.Z(), 0,
			  ey.X(), ey.Y(), ey.Z(), 0,
			  ez.X(), ez.Y(), ez.Z(), 0,
			  0     , 0     , 0     , 1);
  viewMatrix = viewMatrix * Matrix4x4d(1, 0, 0, -eye.X(),
				       0, 1, 0, -eye.Y(),
				       0, 0, 1, -eye.Z(),
				       0, 0, 0, 1);
  //std::cout << viewMatrix << std::endl;
  cam_pos = eye;

  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

void Image3d::Rotate(double angle, double x, double y, double z)
{
  worldMatrix = worldMatrix * RotateMat(angle, x, y, z);
  vMatrix = viewMatrix * worldMatrix;
  covMatrix = projMatrix * vMatrix;
}

Vector3d Image3d::homo(const Vector4d & vec) const
{
  return Vector3d(vec.X / vec.W,
		  vec.Y / vec.W,
		  vec.Z / vec.W);
}

void Image3d::SetOffset(int x, int y)
{
  offsetX = x;
  offsetY = y;
}

Vector3d Image3d::matConv(const Matrix4x4d & mat, const Vector3d & vec) const
{
  double X = mat.At(0, 0) * vec.X() + mat.At(0, 1) * vec.Y()
    + mat.At(0, 2) * vec.Z() + mat.At(0, 3);
  double Y = mat.At(1, 0) * vec.X() + mat.At(1, 1) * vec.Y()
    + mat.At(1, 2) * vec.Z() + mat.At(1, 3);
  double Z = mat.At(2, 0)*vec.X() + mat.At(2, 1) * vec.Y()
    + mat.At(2, 2) * vec.Z() + mat.At(2, 3);
  double W = mat.At(3, 0) * vec.X() + mat.At(3, 1) * vec.Y()
    + mat.At(3, 2) * vec.Z() + mat.At(3, 3);
  return Vector3d(X/W, Y/W, Z/W);
}

double Image3d::matConvZ(const Matrix4x4d & mat, const Vector3d & vec) const
{
  double Z = mat.At(2, 0)*vec.X() + mat.At(2, 1) * vec.Y()
    + mat.At(2, 2) * vec.Z() + mat.At(2, 3);
  double W = mat.At(3, 0) * vec.X() + mat.At(3, 1) * vec.Y()
    + mat.At(3, 2) * vec.Z() + mat.At(3, 3);
  return Z/W;
}

Vector3d Image3d::conversion(const Vector3d& vec) const
{/*
  return homo(covMatrix * Vector4d(vec, 1))
  + Vector3d(offsetX, offsetY, 0);*/
  auto v = matConv(covMatrix, vec);
  v.Add(offsetX, offsetY, 0);
  return v;
}

Vector3d Image3d::lightConversion(const Vector3d & vec) const
{
  throw "Not Implement";
}

Vector3d Image3d::worldConv(const Vector3d & vec) const
{
  return matConv(worldMatrix, vec);
}

double Image3d::worldConvZ(const Vector3d & vec) const
{
  return matConvZ(worldMatrix, vec);
}

Vector3d Image3d::viewConv(const Vector3d & vec) const
{
  return matConv(viewMatrix, vec);
}

double Image3d::viewConvZ(const Vector3d & vec) const
{
  return matConvZ(viewMatrix, vec);
}

void Image3d::Line(const Vector3d & a, const Vector3d & b)
{
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  drawLine(va.X(), va.Y(), vb.X(), vb.Y(), 255, 255, 255, 255, 255, 255);
}

void Image3d::Triangle(const Vector3d & a, const Vector3d & b, const Vector3d & c)
{
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  Vector3d vc = conversion(c);
  drawTriangle(va.X(), va.Y(), 255, 0, 0,
	       vb.X(), vb.Y(), 0, 255, 0,
	       vc.X(), vc.Y(), 0, 0, 255);
}

void Image3d::Triangle(const Vector3d & a, const Vector3d & b, const Vector3d & c,
		       uchar col_r, uchar col_g, uchar col_b)
{
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  Vector3d vc = conversion(c);
  drawTriangle(va.X(), va.Y(), col_r, col_g, col_b,
	       vb.X(), vb.Y(), col_r, col_g, col_b,
	       vc.X(), vc.Y(), col_r, col_g, col_b);
}

void Image3d::drawTriangle(const Vector3d & a,
			   const Vector3d & b,
			   const Vector3d & c)
{
  Vector3d normal = GetNormal(a, b, c);
  drawTriangle_constant(a, b, c, normal);
}

//シャドウマップに使うあれ
void Image3d::drawTriangle_shadow(const Vector3d & a,
				  const Vector3d & b,
				  const Vector3d & c)
{
  Vector3d va = lightConversion(a);
  Vector3d vb = lightConversion(b);
  Vector3d vc = lightConversion(c);
  const Vector3d *pos_a = &a, *pos_b = &b, *pos_c = &c;
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(pos_a, pos_b);
  }
  if (vb.Y() > vc.Y())
  {
    std::swap(vb, vc);
    std::swap(pos_b, pos_c);
  }
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(pos_a, pos_b);
  }
  double b_aX = vb.X() - va.X();
  double c_aX = vc.X() - va.X();
  double b_aY = vb.Y() - va.Y();
  double c_aY = vc.Y() - va.Y();
  int upperStart = va.Y() < 0 ? 0 : va.Y();
  int upperEnd = std::min((int)(vb.Y())+1, (int)(image.height()));
  int lowerStart = vb.Y() < 0 ? 0 : vb.Y();
  int lowerEnd = std::min((int)(vc.Y())+1, (int)(image.height()));
  for (int i = upperStart; i < upperEnd; i++)
  {
    double s = (double)(i - va.Y()) / b_aY;// a ~ b
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = va.X() + b_aX * (i - va.Y()) / b_aY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_a)) + s*(-viewConvZ(*pos_b));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    drawLine_shadow(px, pz, qx, qz, i);
  }
  double b_cX = vb.X() - vc.X();
  double b_cY = vb.Y() - vc.Y();
  double c_bY = -b_cY;
  for (int i = lowerStart; i < lowerEnd; i++)
  {
    double s = (double)(i - vb.Y()) / c_bY;// b ~ c
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = vc.X() + b_cX * (vc.Y() - i) / c_bY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_b)) + s*(-viewConvZ(*pos_c));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    drawLine_shadow(px, pz, qx, qz, i);
  }

}

void Image3d::drawLine_shadow(int px, double pz,
			      int qx, double qz,
			      int h)
{
}


//材質・法線を指定できる三角形の描画。Constant Shading
void Image3d::drawTriangle_constant(const Vector3d & a,
				    const Vector3d & b,
				    const Vector3d & c,
				    const Vector3d & normal)
{
  /*
  if (enableCulling && false &&
      outerProduct(viewConv(b) - viewConv(a),
		   viewConv(c) - viewConv(a)).Z() > 0)
    return;
  */
  Vector3d center = (a + b + c)*(1.0/3);
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  Vector3d vc = conversion(c);
  if (enableCulling &&
      outerProduct(vb - va,
		   vc - va).Z() > 0)
    return;
  Vector3d n = Normalise(normal);
  const Vector3d *pos_a = &a, *pos_b = &b, *pos_c = &c;
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(pos_a, pos_b);
  }
  if (vb.Y() > vc.Y())
  {
    std::swap(vb, vc);
    std::swap(pos_b, pos_c);
  }
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(pos_a, pos_b);
  }
  //色
  Vector3d col = simLight(center, n, nullptr);
  double b_aX = vb.X() - va.X();
  double c_aX = vc.X() - va.X();
  double b_aY = vb.Y() - va.Y();
  double c_aY = vc.Y() - va.Y();
  int upperStart = va.Y() < 0 ? 0 : va.Y();
  int upperEnd = std::min((int)(vb.Y())+1, (int)(image.height()));
  int lowerStart = vb.Y() < 0 ? 0 : vb.Y();
  int lowerEnd = std::min((int)(vc.Y())+1, (int)(image.height()));
  for (int i = upperStart; i < upperEnd; i++)
  {
    double s = (double)(i - va.Y()) / b_aY;// a ~ b
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = va.X() + b_aX * (i - va.Y()) / b_aY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_a)) + s*(-viewConvZ(*pos_b));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    drawLine(px, pz, qx, qz, col, i);
  }
  double b_cX = vb.X() - vc.X();
  double b_cY = vb.Y() - vc.Y();
  double c_bY = -b_cY;
  for (int i = lowerStart; i < lowerEnd; i++)
  {
    double s = (double)(i - vb.Y()) / c_bY;// b ~ c
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = vc.X() + b_cX * (vc.Y() - i) / c_bY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_b)) + s*(-viewConvZ(*pos_c));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    drawLine(px, pz, qx, qz, col, i);
  }
}

// in Constant Shading
void Image3d::drawLine(int px, double pz,
		       int qx, double qz,
		       const Vector3d & col,
		       int h)
{
  if (px == qx || h < 0 || h >= (int)image.height())
    return;
  int ax, bx;
  double az, bz;
  if (px < qx)
  {
    ax = px;
    az = pz;
    bx = qx;
    bz = qz;
  }
  else
  {
    ax = qx;
    az = qz;
    bx = px;
    bz = pz;
  }
  int start = ax < 0 ? 0 : ax;
  int end = bx < (int)image.width() ? bx : image.width();
  for (int i = start; i < end; i++)
  {
    double s = (double)(i - ax) / (bx - ax);
    //    Vector3d vr = (1-s)*(*va) + s*(*vb);
    double z_val = (1-s)*az + s*bz;
    if (enableZBuffer)
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
	if (z_val < z_at(i, h))
	{
	  z_at(i, h) = z_val;
	  drawPoint(i, h, 255*col.X(),
		          255*col.Y(),
		          255*col.Z());
	}
      }//critical end
    }
    else
    {
      drawPoint(i, h, 255*col.X(),
		      255*col.Y(),
		      255*col.Z());      
    }
  }
}



//材質・法線を指定できる三角形の描画。Gouraud Shading
void Image3d::drawTriangle_Gouraud(const Vector3d & a, const Vector3d & normal_a,
				   const Vector3d & b, const Vector3d & normal_b,
				   const Vector3d & c, const Vector3d & normal_c)
{
  /*
  if (enableCulling &&
      outerProduct(viewConv(b) - viewConv(a),
		   viewConv(c) - viewConv(a)).Z() > 0)
    return;
  */
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  Vector3d vc = conversion(c);
  if (enableCulling &&
      outerProduct(vb - va,
		   vc - va).Z() > 0)
    return;
  Vector3d na = Normalise(normal_a);
  Vector3d nb = Normalise(normal_b);
  Vector3d nc = Normalise(normal_c);
  const Vector3d *pos_a = &a, *pos_b = &b, *pos_c = &c;
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(na, nb);
    std::swap(pos_a, pos_b);
  }
  if (vb.Y() > vc.Y())
  {
    std::swap(vb, vc);
    std::swap(nb, nc);
    std::swap(pos_b, pos_c);
  }
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(na, nb);
    std::swap(pos_a, pos_b);
  }
  //点a, b, cの色
  Vector3d col_a = simLight(*pos_a, na, nullptr);
  Vector3d col_b = simLight(*pos_b, nb, nullptr);
  Vector3d col_c = simLight(*pos_c, nc, nullptr);
  double b_aX = vb.X() - va.X();
  double c_aX = vc.X() - va.X();
  double b_aY = vb.Y() - va.Y();
  double c_aY = vc.Y() - va.Y();
  int upperStart = va.Y() < 0 ? 0 : va.Y();
  int upperEnd = std::min((int)(vb.Y())+1, (int)(image.height()));
  int lowerStart = vb.Y() < 0 ? 0 : vb.Y();
  int lowerEnd = std::min((int)(vc.Y())+1, (int)(image.height()));
  for (int i = upperStart; i < upperEnd; i++)
  {
    double s = (double)(i - va.Y()) / b_aY;// a ~ b
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = va.X() + b_aX * (i - va.Y()) / b_aY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_a)) + s*(-viewConvZ(*pos_b));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    //点p, qの線形補完により求められた色。
    Vector3d col_p = (1-s)*col_a + s*col_b;
    Vector3d col_q = (1-t)*col_a + t*col_c;
    drawLine(px, pz, col_p, qx, qz, col_q, i);
  }
  double b_cX = vb.X() - vc.X();
  double b_cY = vb.Y() - vc.Y();
  double c_bY = -b_cY;
  for (int i = lowerStart; i < lowerEnd; i++)
  {
    double s = (double)(i - vb.Y()) / c_bY;// b ~ c
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = vc.X() + b_cX * (vc.Y() - i) / c_bY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //
    double pz = (1-s)*(-viewConvZ(*pos_b)) + s*(-viewConvZ(*pos_c));
    double qz = (1-t)*(-viewConvZ(*pos_a)) + t*(-viewConvZ(*pos_c));
    //点p, qの線形補完により求められた色。
    Vector3d col_p = (1-s)*col_b + s*col_c;
    Vector3d col_q = (1-t)*col_a + t*col_c;
    drawLine(px, pz, col_p, qx, qz, col_q, i);
  }
}
//in Gouraud Shading
void Image3d::drawLine(int px, double pz, const Vector3d & col_p,
		       int qx, double qz, const Vector3d & col_q,
		       int h)
{
  if (px == qx || h < 0 || h >= (int)image.height())
    return;
  int ax, bx;
  double az, bz;
  const Vector3d * col_a;
  const Vector3d * col_b;
  if (px < qx)
  {
    ax = px;
    az = pz;
    col_a = &col_p;
    bx = qx;
    bz = qz;
    col_b = &col_q;
  }
  else
  {
    ax = qx;
    az = qz;
    col_a = &col_q;
    bx = px;
    bz = pz;
    col_b = &col_p;
  }
  int start = ax < 0 ? 0 : ax;
  int end = bx < (int)image.width() ? bx : image.width();
  for (int i = start; i < end; i++)
  {
    double s = (double)(i - ax) / (bx - ax);
    //    Vector3d vr = (1-s)*(*va) + s*(*vb);
    double z_val = (1-s)*az + s*bz;
    if (enableZBuffer)
    {
      if (z_val < z_at(i, h))
      {
	Vector3d col = (1-s)*(*col_a) + s*(*col_b);
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	  if (z_val < z_at(i, h))
	  {
	    z_at(i, h) = z_val;
	    drawPoint(i, h, 255*col.X(),
	  	            255*col.Y(),
		            255*col.Z());
	  }
	}//critical end
      }
    }
    else
    {
      Vector3d col = (1-s)*(*col_a) + s*(*col_b);
      drawPoint(i, h, 255*col.X(),
		      255*col.Y(),
		      255*col.Z());
    } 
  }
}


void Image3d::drawTriangle_Phong(uint a, uint b, uint c)
{
  if (enableCulling &&
      outerProduct(vertexes[b].convedPos - vertexes[a].convedPos,
		   vertexes[c].convedPos - vertexes[a].convedPos).Z() > 0)
    return;
  
  if (vertexes[a].convedPos.Y() > vertexes[b].convedPos.Y())
  {
    std::swap(a, b);
  }
  if (vertexes[b].convedPos.Y() > vertexes[c].convedPos.Y())
  {
    std::swap(b, c);
  }
  if (vertexes[a].convedPos.Y() > vertexes[b].convedPos.Y())
  {
    std::swap(a, b);
  }
  const Vector3d & va = vertexes[a].convedPos;
  const Vector3d & vb = vertexes[b].convedPos;
  const Vector3d & vc = vertexes[c].convedPos;
  double b_aX = vb.X() - va.X();
  double c_aX = vc.X() - va.X();
  double b_aY = vb.Y() - va.Y();
  double c_aY = vc.Y() - va.Y();
  int upperStart = va.Y() < 0 ? 0 : va.Y();
  int upperEnd = std::min((int)(vb.Y())+1, (int)(image.height()));
  for (int i = upperStart; i < upperEnd; i++)
  {
    double s = (double)(i - va.Y()) / b_aY;// a ~ b
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = va.X() + b_aX * (i - va.Y()) / b_aY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //線形補完により求められた世界座標での点p, qの座標。
    //光源とか深度とかの計算に使う。
    Vector3d pwp = blend(1-s, vertexes[a].worldPos,
			  s, vertexes[b].worldPos);
    Vector3d qwp = blend(1-t, vertexes[a].worldPos,
			  t, vertexes[c].worldPos);
    //点p, qの線形補完により求められた世界座標での法線。光源とか(ry
    Vector3d pwn = blend(1-s, vertexes[a].worldNormal,
			   s, vertexes[b].worldNormal);
    Vector3d qwn = blend(1-t, vertexes[a].worldNormal,
			   t, vertexes[c].worldNormal);
    //TexCoord
    Vector2d ptp = blend(1-s, vertexes[a].texPos,
			   s, vertexes[b].texPos);
    Vector2d qtp = blend(1-t, vertexes[a].texPos,
			   t, vertexes[c].texPos);
    if (px < qx)
      drawLine_Phong(px, pwp, pwn, ptp,
		     qx, qwp, qwn, qtp,
		     i);
    else if (qx < px)
      drawLine_Phong(qx, qwp, qwn, qtp,
		     px, pwp, pwn, ptp,
		     i);
  }
  double b_cX = vb.X() - vc.X();
  double b_cY = vb.Y() - vc.Y();
  double c_bY = -b_cY;
  int lowerStart = vb.Y() < 0 ? 0 : vb.Y();
  int lowerEnd = std::min((int)(vc.Y())+1, (int)(image.height()));
  for (int i = lowerStart; i < lowerEnd; i++)
  {
    double s = (double)(i - vb.Y()) / c_bY;// b ~ c
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
      continue;
    //画像平面上での点p, qのx座標
    int px = vc.X() + b_cX * (vc.Y() - i) / c_bY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //世界座標での点p, qの座標。光源とかの計算に使う。
    Vector3d pwp = blend(1-s, vertexes[b].worldPos,
			   s, vertexes[c].worldPos);
    Vector3d qwp = blend(1-t, vertexes[a].worldPos,
			   t, vertexes[c].worldPos);
    //点p, qの線形補完により求められた世界座標での法線。光源とか(ry
    Vector3d pwn = blend(1-s, vertexes[b].worldNormal,
			   s, vertexes[c].worldNormal);
    Vector3d qwn = blend(1-t, vertexes[a].worldNormal,
			   t, vertexes[c].worldNormal);
    //TexCoord
    Vector2d ptp = blend(1-s, vertexes[b].texPos,
			  s, vertexes[c].texPos);
    Vector2d qtp = blend(1-t, vertexes[a].texPos,
			  t, vertexes[c].texPos);
    if (px < qx)
      drawLine_Phong(px, pwp, pwn, ptp,
		     qx, qwp, qwn, qtp,
		     i);
    else if (qx < px)
      drawLine_Phong(qx, qwp, qwn, qtp,
		     px, pwp, pwn, ptp,
		     i);
  }
}

void Image3d::drawLine_Phong(int ax, const Vector3d & awp,
			     const Vector3d & awn, const Vector2d & atp,
			     int bx, const Vector3d & bwp,
			     const Vector3d & bwn, const Vector2d & btp,
			     int h)
{
  int start = ax < 0 ? 0 : ax;
  int end = bx < (int)image.width() ? bx : image.width();
  int b_a = bx - ax;
  for (int i = start; i < end; i++)
  {
    double s = (double)(i - ax) / b_a;
    double ss = 1-s;
    Vector3d rwp = blend(ss, awp, s, bwp);
    double z_val = -viewConvZ(rwp);
    Vector2d rtp = blend(ss, atp, s, btp);
    if (enableZBuffer)
    {
      if (z_val < z_at(i, h))
      {
	Vector3d rwn = blend(ss, awn, s, bwn);
	Vector3d col = simColor(rwp, rwn, &rtp);
#ifdef _OPENMP
#pragma omp critical
#endif
	{//start critical
	  if (z_val < z_at(i, h))
	  {
	    z_at(i, h) = z_val;
	    drawPoint(i,
		      h,
		      255*col.X(),
		      255*col.Y(),
		      255*col.Z());
	  }
	}//end critical
      }
    }
    else
    {
      Vector3d rwn = blend(ss, awn, s, bwn);
      Vector3d col = simColor(rwp, rwn, &rtp);
      z_at(i, h) = z_val;
      drawPoint(i,
		h,
		255*col.X(),
		255*col.Y(),
		255*col.Z());
    }
  }
}

Vector3d Image3d::simColor(const Vector3d & worldPos,
			   const Vector3d & worldNormal,
			   const Vector2d * texPos)
{
  Vector3d light(mtrAmbX,
		 mtrAmbY,
		 mtrAmbZ);
  //Vector3d i(Light0.pos.X, Light0.pos.Y, Light0.pos.Z);
  //i.Normalise();

  const Vector3d & n = worldNormal;

  Vector3d s = cam_pos - worldPos;
  s.Normalise();
  s.Sub(normalised_i0);
  s.Normalise();

  //diffusion
  double dot = -innerProduct(normalised_i0, n);
  if (dot > 0)
    light.Add(dot * mtrDifX,
	      dot * mtrDifY,
	      dot * mtrDifZ);


  //specular
  dot = innerProduct(s, n);
  if (dot > 0)
  {
    double d = pow(dot, mtrShininess);
    light.Add(d * mtrSpcX,
	      d * mtrSpcY,
	      d * mtrSpcZ);
  }



  //環境マッピング
  if (envTex)
  {
    Vector3d u = worldPos - cam_pos;
    Vector3d f = u - 2 * n * innerProduct(n, u);
    f.Normalise();
    light = envTex->GetColor(f.X(), f.Y());
  }


  if (light.X() > 1)
    light.X() = 1;
  if (light.Y() > 1)
    light.Y() = 1;
  if (light.Z() > 1)
    light.Z() = 1;
  
  //テクスチャマッピング
  if (texture && texPos != nullptr)
  {
    Vector3d tex = texture->GetColor((texPos->X * 2) - 1,
				     (texPos->Y * 2) - 1);
    light.Set(tex.Z() * light.X(), tex.Y() * light.Y(), tex.X() * light.Z());
  }

  
  
  return light;

}


//材質・法線を指定できる三角形の描画。Phong Shading
void Image3d::drawTriangle_Phong(const Vector3d & a, const Vector3d & normal_a, const Vector2d & ta,
				 const Vector3d & b, const Vector3d & normal_b, const Vector2d & tb,
				 const Vector3d & c, const Vector3d & normal_c, const Vector2d & tc)
{
  /*
  if (enableCulling &&
      outerProduct(viewConv(b) - viewConv(a),
		   viewConv(c) - viewConv(a)).Z() > 0)
		   return;*/
  Vector3d va = conversion(a);
  Vector3d vb = conversion(b);
  Vector3d vc = conversion(c);
  //たぶん、これでも大丈夫だよね
  if (enableCulling &&
      outerProduct(vb - va,
		   vc - va).Z() > 0)
    return;
  Vector3d na = Normalise(normal_a);
  Vector3d nb = Normalise(normal_b);
  Vector3d nc = Normalise(normal_c);
  const Vector3d *pos_a = &a, *pos_b = &b, *pos_c = &c;
  const Vector2d *tex_a = &ta, *tex_b = &tb, *tex_c = &tc;
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(na, nb);
    std::swap(pos_a, pos_b);
    std::swap(tex_a, tex_b);
  }
  if (vb.Y() > vc.Y())
  {
    std::swap(vb, vc);
    std::swap(nb, nc);
    std::swap(pos_b, pos_c);
    std::swap(tex_b, tex_c);
  }
  if (va.Y() > vb.Y())
  {
    std::swap(va, vb);
    std::swap(na, nb);
    std::swap(pos_a, pos_b);
    std::swap(tex_a, tex_b);
  }
  double b_aX = vb.X() - va.X();
  double c_aX = vc.X() - va.X();
  double b_aY = vb.Y() - va.Y();
  double c_aY = vc.Y() - va.Y();
  int upperStart = va.Y() < 0 ? 0 : va.Y();
  int upperEnd = std::min((int)(vb.Y())+1, (int)(image.height()));
  for (int i = upperStart; i < upperEnd; i++)
  {
    double s = (double)(i - va.Y()) / b_aY;// a ~ b
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
    {
      //std::cerr << va.Y << ":" << vb.Y << ":" << vc.Y << std::endl;
      continue;
    }
    //画像平面上での点p, qのx座標
    int px = va.X() + b_aX * (i - va.Y()) / b_aY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //元の空間での点p, qの座標。光源とかの計算に使う。
    Vector3d vp = blend(1-s, *pos_a, s, *pos_b);
    Vector3d vq = blend(1-t, *pos_a, t, *pos_c);
    //点p, qの線形補完により求められた法線。光源とか(ry
    Vector3d np = blend(1-s, na, s, nb);
    Vector3d nq = blend(1-t, na, t, nc);
    //TexCoord
    Vector2d tp = blend(1-s, ta, s, tb);
    Vector2d tq = blend(1-t, ta, t, tc);
    drawLine(px, vp, np, tp,
	     qx, vq, nq, tq, i);
  }
  double b_cX = vb.X() - vc.X();
  double b_cY = vb.Y() - vc.Y();
  double c_bY = -b_cY;
  int lowerStart = vb.Y() < 0 ? 0 : vb.Y();
  int lowerEnd = std::min((int)(vc.Y())+1, (int)(image.height()));
  for (int i = lowerStart; i < lowerEnd; i++)
  {
    double s = (double)(i - vb.Y()) / c_bY;// b ~ c
    double t = (double)(i - va.Y()) / c_aY;// a ~ c
    if (s < 0 || s > 1 || t < 0 || t > 1)
    {
      //std::cerr << va.Y << ":" << vb.Y << ":" << vc.Y << std::endl;
      continue;
    }
    //画像平面上での点p, qのx座標
    int px = vc.X() + b_cX * (vc.Y() - i) / c_bY;
    int qx = va.X() + c_aX * (i - va.Y()) / c_aY;
    //元の空間での点p, qの座標。光源とかの計算に使う。
    Vector3d vp = blend(1-s, *pos_b, s, *pos_c);
    Vector3d vq = blend(1-t, *pos_a, t, *pos_c);
    //点p, qの線形補完により求められた法線。光源とか(ry
    Vector3d np = blend(1-s, nb, s, nc);
    Vector3d nq = blend(1-t, na, t, nc);
    //std::cout << s << std::endl;
    Vector2d tp = blend(1-s, tb, s, tc);
    Vector2d tq = blend(1-t, ta, t, tc);
    drawLine(px, vp, np, tp,
	     qx, vq, nq, tq, i);
  }
}


//法線・材質指定できる。 in Phong
void Image3d::drawLine(int px, const Vector3d & vp, const Vector3d & np, const Vector2d & tp,
		       int qx, const Vector3d & vq, const Vector3d & nq, const Vector2d & tq,
		       int h)
{
  //std::cout << px << ":" << qx << ":" << h << std::endl;
  if (px == qx || h < 0 || h >= (int)image.height())
    return;
  int ax, bx;
  const Vector3d * va;
  const Vector3d * na;
  const Vector3d * vb;
  const Vector3d * nb;
  const Vector2d * ta;
  const Vector2d * tb;
  if (px < qx)
  {
    ax = px;
    va = &vp;
    na = &np;
    ta = &tp;
    bx = qx;
    vb = &vq;
    nb = &nq;
    tb = &tq;
  }
  else
  {
    ax = qx;
    va = &vq;
    na = &nq;
    ta = &tq;
    bx = px;
    vb = &vp;
    nb = &np;
    tb = &tp;
  }
  int start = ax < 0 ? 0 : ax;
  int end = bx < (int)image.width() ? bx : image.width();
  int b_a = bx - ax;
  for (int i = start; i < end; i++)
  {
    double s = (double)(i - ax) / b_a;
    double ss = 1-s;
    Vector3d vr = blend(ss, *va, s, *vb);
    double z_val = -viewConvZ(vr);
    Vector2d tr = blend(ss, *ta, s, *tb);
    if (enableZBuffer)
    {
      if (z_val < z_at(i, h))
      {
	Vector3d nr = blend(ss, *na, s, *nb);
	Vector3d col = simLight(vr, nr, &tr);
#ifdef _OPENMP
#pragma omp critical
#endif
	{//start critical
	  if (z_val < z_at(i, h))
	  {
	    z_at(i, h) = z_val;
	    drawPoint(i, h, 255*col.X(),
		            255*col.Y(),
		            255*col.Z());
	  }
	}//critical end
      }
    }
    else
    {
      Vector3d nr = blend(ss, *na, s, *nb);
      Vector3d col = simLight(vr, nr, &tr);
      z_at(i, h) = z_val;
      drawPoint(i, h, 255*col.X(),
		      255*col.Y(),
		      255*col.Z());      
    }
  }
}


void Image3d::drawPoint(int x, int y, uchar r, uchar g, uchar b)
{
  image.at(x, y, kazakami::colorRGB::R) = r;
  image.at(x, y, kazakami::colorRGB::G) = g;
  image.at(x, y, kazakami::colorRGB::B) = b;
/*
  dif = 0;
  int rx = x - dif;
  int ry = y - dif;
  int bx = x + dif;
  int by = y - dif;
  //dif = 6;
  int gx = x * (1 + dif/20.0);//x + dif;
  int gy = y * (1 + dif/20.0);//y - dif;
  if (rx > 0 && rx < (int)width() && ry > 0 && ry < (int)height())
    image.at(rx, ry, kazakami::colorRGB::R) = r;
  if (gx > 0 && gx < (int)width() && gy > 0 && gy < (int)height())
    image.at(gx, gy, kazakami::colorRGB::G) = g;
  if (bx > 0 && bx < (int)width() && by > 0 && by < (int)height())
    image.at(bx, by, kazakami::colorRGB::B) = b;
 */
}

void Image3d::drawLine(int ax, int ay, int bx, int by,
		       uchar ar, uchar ag, uchar ab,
		       uchar br, uchar bg, uchar bb)
{
  int deltax = abs(bx - ax), deltay = abs(by - ay);
  if (deltax > deltay)
  {
    if (ax > bx)
    {
      std::swap(ax, bx);
      std::swap(ay, by);
      std::swap(ar, br);
      std::swap(ag, bg);
      std::swap(ab, bb);
    }
    for (int i = 0; i < deltax; i++)
    {
      if (ax + i < 0)
	i = -ax;
      if (ax + i >= (int)image.width())
	break;
      int x = ax + i;
      int y = ay + ((by - ay) * i) / deltax;
      if (y < 0 || y >= (int)image.height())
	continue;
      drawPoint(x, y,
		ar + ((br - ar) * i) / deltax,
		ag + ((bg - ag) * i) / deltax,
		ab + ((bb - ab) * i) / deltax);
    }
  }
  else
  {
    if (ay > by)
    {
      std::swap(ay, by);
      std::swap(ay, by);
      std::swap(ar, br);
      std::swap(ag, bg);
      std::swap(ab, bb);
    }
    for (int i = 0; i < deltay; i++)
    {
      if (ay + i < 0)
	i = -ay;
      if (ay + i >= (int)image.height())
	break;
      int x = ax + ((bx - ax) * i) / deltay;
      int y = ay + i;
      if (x < 0 || x >= (int)image.width())
	continue;
      drawPoint(x, y,
		ar + ((br - ar) * i) / deltay,
		ag + ((bg - ag) * i) / deltay,
		ab + ((bb - ab) * i) / deltay);
    }
  }
}

void Image3d::drawTriangle(int ax, int ay, uchar ar, uchar ag, uchar ab,
			   int bx, int by, uchar br, uchar bg, uchar bb,
			   int cx, int cy, uchar cr, uchar cg, uchar cb)
{
  if (ay > by)
  {
    std::swap(ax, bx);
    std::swap(ay, by);
    std::swap(ar, br);
    std::swap(ag, bg);
    std::swap(ab, bb);
  }
  if (by > cy)
  {
    std::swap(cx, bx);
    std::swap(cy, by);
    std::swap(cr, br);
    std::swap(cg, bg);
    std::swap(cb, bb);
  }
  if (ay > by)
  {
    std::swap(ax, bx);
    std::swap(ay, by);
    std::swap(ar, br);
    std::swap(ag, bg);
    std::swap(ab, bb);
  }
  uint start = ay < 0 ? 0 : ay;
  uint end = by < (int)image.height() ? by : image.height();
  for (uint i = start; i < end; i++)
  {
    int px = ax + ((bx - ax) * (i - ay)) / (by - ay);
    int qx = ax + ((cx - ax) * (i - ay)) / (cy - ay);
    uchar pr = ar + ((br - ar) * (i - ay)) / (by - ay);
    uchar pg = ag + ((bg - ag) * (i - ay)) / (by - ay);
    uchar pb = ab + ((bb - ab) * (i - ay)) / (by - ay);
    uchar qr = ar + ((cr - ar) * (i - ay)) / (cy - ay);
    uchar qg = ag + ((cg - ag) * (i - ay)) / (cy - ay);
    uchar qb = ab + ((cb - ab) * (i - ay)) / (cy - ay);
    drawLine(px, i, qx, i, pr, pg, pb, qr, qg, qb);
  }
  start = by < 0 ? 0 : by;
  end = cy < (int)image.height() ? cy : image.height();
  for (uint i = start; i < end; i++)
  {
    int px = cx + ((bx - cx) * (cy - i)) / (cy - by);
    int qx = ax + ((cx - ax) * (i - ay)) / (cy - ay);
    uchar pr = cr + ((br - cr) * (cy - i)) / (cy - by);
    uchar pg = cg + ((bg - cg) * (cy - i)) / (cy - by);
    uchar pb = cb + ((bb - cb) * (cy - i)) / (cy - by);
    uchar qr = ar + ((cr - ar) * (i - ay)) / (cy - ay);
    uchar qg = ag + ((cg - ag) * (i - ay)) / (cy - ay);
    uchar qb = ab + ((cb - ab) * (i - ay)) / (cy - ay);
    drawLine(px, i, qx, i, pr, pg, pb, qr, qg, qb);
  }
}

double Image3d::innerProduct(const Vector3d & a, const Vector3d & b) const
{
  return a.X() * b.X() + a.Y() * b.Y() + a.Z() * b.Z();
}

Vector3d Image3d::outerProduct(const Vector3d & a, const Vector3d & b) const
{
  return Vector3d(a.Y() * b.Z() - a.Z() * b.Y(),
		  a.Z() * b.X() - a.X() * b.Z(),
		  a.X() * b.Y() - a.Y() * b.X());
}

void Image3d::SetLightPos(const Vector3d & pos)
{
  Light0.pos = worldMatrix * Vector4d(pos, 1);
  normalised_i0.Set(Light0.pos.X, Light0.pos.Y, Light0.pos.Z);
  normalised_i0.Normalise();
}

void Image3d::SetLightCol(double r, double g, double b)
{
  Light0.col[0] = r;
  Light0.col[1] = g;
  Light0.col[2] = b;
}

void Image3d::SetMaterial(const Material & mtr)
{
  mtrAmbX = mtr.ambient[0] * Light0.col[0];
  mtrAmbY = mtr.ambient[1] * Light0.col[1];
  mtrAmbZ = mtr.ambient[2] * Light0.col[2];
  mtrDifX = mtr.diffuse[0] * Light0.col[0];
  mtrDifY = mtr.diffuse[1] * Light0.col[1];
  mtrDifZ = mtr.diffuse[2] * Light0.col[2];
  mtrSpcX = mtr.specular[0] * Light0.col[0];
  mtrSpcY = mtr.specular[1] * Light0.col[1];
  mtrSpcZ = mtr.specular[2] * Light0.col[2];
  mtrShininess = mtr.shininess;
}

Vector3d Image3d::simLight(const Vector3d & pos,
			   const Vector3d & normal,
			   const Vector2d * texPos) const
{
  //(-(i, n)*kd + (s, n)^m*ks + ka) * Il

  //environment
  Vector3d light(mtrAmbX,
		 mtrAmbY,
		 mtrAmbZ);

  Vector3d convedPos = worldConv(pos);
  Vector3d i(Light0.pos.X, Light0.pos.Y, Light0.pos.Z);

  double a = 1.0;

  if (!isParallelLight)
  {
    i.Sub(convedPos);
    a = 1 / (i.Norm() * i.Norm() / 500000);
  }
  

  i.Normalise();

  Vector3d n = worldConv(pos + normal) - convedPos;
  n.Normalise();

  Vector3d s = cam_pos - convedPos;
  s.Normalise();
  s.Sub(i);
  s.Normalise();

  //diffusion
  double dot = -innerProduct(i, n) * a;
  if (dot > 0)
    light.Add(dot * mtrDifX,
	      dot * mtrDifY,
	      dot * mtrDifZ);


  //specular
  dot = innerProduct(s, n) * a;
  if (dot > 0)
  {
    double d = pow(dot, mtrShininess);
    light.Add(d * mtrSpcX,
	      d * mtrSpcY,
	      d * mtrSpcZ);
  }



  //環境マッピング
  if (envTex)
  {
    Vector3d u = convedPos - cam_pos;
    Vector3d f = u - 2 * n * innerProduct(n, u);
    f.Normalise();
    light = envTex->GetColor(f.X(), f.Y());
    //light = blend(0.5, envTex->GetColor(f.X(), f.Y()), 0.5, light);
  }

  //
  if (texture && texPos != nullptr)
  {
    light = texture->GetColor((texPos->X * 2) - 1,
			      (texPos->Y * 2) - 1);
    light.Set(light.Z(), light.Y(), light.X());
  }

  
  if (light.X() > 1)
    light.X() = 1;
  if (light.Y() > 1)
    light.Y() = 1;
  if (light.Z() > 1)
    light.Z() = 1;
  
  return light;
}


}
