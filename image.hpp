#pragma once

#include <vector>
//#include <cstddef>
#include <fstream>
#include <math.h>
#include <memory>
#include <map>
#include <stack>


#include "Matrix4x4d.hpp"
#include "Quaternion.hpp"
#include "vmdLoader.hpp"


namespace kazakami
{

typedef unsigned char uchar;
typedef unsigned int uint;



enum class colorRGB
{
  B = 2, G = 1, R = 0
};

struct Material
{
  double ambient[3];
  double diffuse[3];
  double specular[3];
  double shininess;
};


class Texture
{
  std::vector<uchar> data;
  size_t sx;
  size_t sy;
public:
  uchar * ptr();
  const uchar * ptr() const;
  size_t height() const;
  size_t width() const;
  bool LoadPPM(const char * filename);
  void LoadFromVector(const std::vector<uchar> & d, size_t x, size_t y);
  void LoadFromVector(std::vector<uchar> && d, size_t x, size_t y);
  //画像サイズをx, yともに[-1, 1]とし、線形補間された色を返す。
  Vector3d GetColor(double x, double y) const;
};

class Image
{
  size_t sx, sy;
  std::vector<uchar> data;
public:
  size_t height() const;
  size_t width() const;
  Image(size_t x, size_t y);
  uchar & at(size_t x, size_t y, colorRGB col);
  const uchar & at(size_t x, size_t y, colorRGB col) const;
  uchar * ptr();
  const uchar * ptr() const;
  //画像データをbmp形式で出力する。
  //返り値はファイルサイズ(バイト単位)で、失敗した場合負値となる。
  int writeOut(const char * filename) const;
  //ppm形式で画像を表示する。
  //返り値は0で成功、負値で失敗。
  int writeOutWithPPM(const char * filename) const;
};

//3D図形の描画機能を持つImage
class Image3d
{
  Image image;
  std::vector<double> shadowMap;
  std::vector<double> zBuffer;
  //Zバッファを有効にするか
  bool enableZBuffer;
  //隠面処理を有効にするか
  bool enableCulling;
  //カメラ位置
  Vector3d cam_pos;
  //平行光源か
  bool isParallelLight;
  //画面左上の撮像領域上での座標
  int offsetX, offsetY;
  //ワールド座標変換行列
  Matrix4x4d worldMatrix;
  std::stack<Matrix4x4d> worldMatrixStack;
  //ビュー座標変換行列
  Matrix4x4d viewMatrix;
  std::stack<Matrix4x4d> viewMatrixStack;
  //プロジェクション座標変換行列
  Matrix4x4d projMatrix;
  std::stack<Matrix4x4d> projMatrixStack;
  //viewMatrix * worldMatrix
  Matrix4x4d vMatrix;
  std::stack<Matrix4x4d> vMatrixStack;
  //projMatrix * viewMatrix * worldMatrix
  Matrix4x4d covMatrix;
  std::stack<Matrix4x4d> covMatrixStack;
  //光源視点座標への変換行列
  //普通の描画のビュー変換に対応する行列
  Matrix4x4d lightMatrix;
  //環境マッピング用のテクスチャ
  std::shared_ptr<Texture> envTex;
  //テクスチャマッピング用のテクスチャ
  std::shared_ptr<Texture> texture;
  //設定されている材質に対するambient
  double mtrAmbX, mtrAmbY, mtrAmbZ;
  //設定されている材質に対するdidffuse
  double mtrDifX, mtrDifY, mtrDifZ;
  //設定されている材質に対するspecular
  double mtrSpcX, mtrSpcY, mtrSpcZ;
  //設定されている材質のshininess
  double mtrShininess;
  //光源
  struct Light
  {
    Vector4d pos;
    double col[3]; 
  } Light0;
  //光源方向を単位化したもの
  Vector3d normalised_i0;
  //深度値を書き換え・取得する。
  double & z_at(size_t x, size_t y);
  const double & z_at(size_t x, size_t y) const;
  //シャドウマップの値を書き換え・取得する。
  double & shadow_at(size_t x, size_t y);
  const double & shadow_at(size_t x, size_t y) const;
  //ただの内積
  double innerProduct(const Vector3d & a, const Vector3d & b) const;
  //多田野外積
  Vector3d outerProduct(const Vector3d & a, const Vector3d & b) const;
  //三次元座標を変換する。
  Vector3d matConv(const Matrix4x4d & mat, const Vector3d & vec) const;
  //三次元座標を変換した時のZ座標を求める。
  double matConvZ(const Matrix4x4d & mat, const Vector3d & vec) const;
  //同次座標を3次元座標に変換する。homogenious
  Vector3d homo(const Vector4d & vec) const;
  //空間内の座標をキャンパス上の座標に変換する。
  Vector3d conversion(const Vector3d& vec) const;
  //世界座標をview座標に変換する
  Vector3d viewConv(const Vector3d & vec) const;
  //世界座標をview座標に変換したZ座標を求める。
  double viewConvZ(const Vector3d & vec) const;
  //ローカル座標を世界座標に変換する
  Vector3d worldConv(const Vector3d & vec) const;
  //ローカル座標を世界座標に変換したZ座標を求める。
  double worldConvZ(const Vector3d & vec) const;
  //ローカル座標を光源スクリーンに投影した座標を求める。
  Vector3d lightConversion(const Vector3d & vec) const;
  //ローカル座標を光源座標に変換する。
  Vector3d lightConv(const Vector3d & vec) const;
  //ローカル座標を光源座標に変換したZ座標を求める。
  double lightConvZ(const Vector3d & vec) const;
  //光源の影響を計算する。
  Vector3d simLight(const Vector3d & pos,
		    const Vector3d & norm,
		    const Vector2d * texPos) const;
  //モーション
  std::shared_ptr<vmdLoader> vmdData;
  //モーション用のフレーム
  uint motionFrame;
  //ボーン!!
  struct Bone
  {
    std::string name;
    std::string nameEng;
    float pos[3];
    float initPos[3];
    int parentBoneIndex;
    Matrix4x4d initMat;
    Matrix4x4d boneMat;
    Matrix4x4d offsetMat;
    Quaternion rotation;
    std::vector<int> children;
  };
  //全てのボーン
  std::vector<Bone> bones;
  //ルートボーンのインデックスを保持しておく
  //ルートはひとつとは限らない。
  std::vector<int> rootBoneIndexes;
  //ボーンの名前からインデックスを取るためのマップ
  std::map<std::string, int> boneNameMap;
  //ボーンを動かす。
  void MoveBone();
  //ボーンの姿勢行列の更新
  //CulcBone()は各ルートボーンについてCulcBone(const Bone&)を呼び出す。
  void CulcBone();
  //子ボーンを再帰的に呼び出して更新する。
  void CulcBone(uint index);
  //内部で頂点データを持っておく
  struct Vertex
  {
    //ローカル座標
    Vector3d pos;
    //世界座標
    Vector3d worldPos;
    //カメラ座標 <- 使ってない
    Vector3d viewPos;
    //画面上の座標
    Vector3d convedPos;
    //ローカル法線
    Vector3d normal;
    //ローカルでの座標＋法線
    Vector3d normalPlusPos;
    //世界座標での法線
    Vector3d worldNormal;
    Vector2d texPos;
    unsigned char weightType;
    int boneIndex1, boneIndex2, boneIndex3, boneIndex4;
    float weight1, weight2, weight3, weight4;
  };
  std::vector<Vertex> vertexes;
public:
  //なんかあれ
  void SetMotion(std::shared_ptr<vmdLoader> vmd);
  void PushMatrix();
  void PopMatrix();
  double z_get(size_t x, size_t y);
  int dif;
  Image3d(size_t x, size_t y);
  enum class Function
  {
    DepthTest,
    CullFace,
  };
  void AddVertex(const Vector3d & pos,
		 const Vector3d & normal,
		 const Vector2d & texPos);
  void AddVertex(const Vector3d & pos,
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
		 float weight4);
  void AddBone(std::string name,
	       std::string nameEng,
	       float X,
	       float Y,
	       float Z,
	       int parentIndex);
  //子から親に張られているインデックスを元に、
  //親から子に張られるインデックスを作る。
  //初期姿勢ボーンを計算したりする。
  void InitBone();
  void CulcVertex();
  void SetMaterial(const Material & mtr);
  void SetEnvMapTex(std::shared_ptr<Texture> tex);
  void SetTexture(std::shared_ptr<Texture> tex);
  void Enable(Function f);
  void Disable(Function f);
  size_t height() const;
  size_t width() const;
  void SetOffset(int x, int y);
  uchar & at(size_t x, size_t y, colorRGB col);
  const uchar & at(size_t x, size_t y, colorRGB col) const;
  uchar * ptr();
  const uchar * ptr() const;
  //void clear(uchar r = 0, uchar g = 0, uchar b = 0);
  void ResetZBuffer();
  int writeOut(const char * filename) const;
  int writeOutWithPPM(const char * filename) const;
  void LoadIdentity();
  void Parallel(double z);
  void Perspective(double z);
  void LookAt(const Vector3d & eye,
	      const Vector3d & lookat,
	      const Vector3d & up);
  void Translate(double x, double y, double z);
  void Scale(double x, double y, double z);
  void Rotate(double angle, double x, double y, double z);
  void SetLightPos(const Vector3d & pos);
  void SetLightCol(double r, double g, double b);

  void Line(const Vector3d & a, const Vector3d & b);
  void Triangle(const Vector3d & a, const Vector3d & b, const Vector3d & c);
  void Triangle(const Vector3d & a, const Vector3d & b, const Vector3d & c,
		uchar col_r, uchar col_b, uchar col_g);

  void drawTriangle_shadow(const Vector3d & a, const Vector3d & b, const Vector3d & c);
  void drawLine_shadow(int px, double pz,
		       int qx, double qz,
		       int h);

  void drawTriangle_constant(const Vector3d & a,
			     const Vector3d & b,
			     const Vector3d & c,
			     const Vector3d & normal);
  void drawTriangle_Gouraud(const Vector3d & a, const Vector3d & normal_a,
			    const Vector3d & b, const Vector3d & normal_b,
			    const Vector3d & c, const Vector3d & normal_c);
  void drawTriangle_Phong(const Vector3d & a, const Vector3d & normal_a, const Vector2d & ta,
			  const Vector3d & b, const Vector3d & normal_b, const Vector2d & tb,
			  const Vector3d & c, const Vector3d & normal_c, const Vector2d & tc);
  //登録された頂点番号で三角形を描画する。
  void drawTriangle_Phong(uint a, uint b, uint c);
  void drawLine_Phong(int ax, const Vector3d & awp,
		      const Vector3d & awn, const Vector2d & atp,
		      int bx, const Vector3d & bwp,
		      const Vector3d & bwn, const Vector2d & btp,
		      int h);
  Vector3d simColor(const Vector3d & worldPos,
		    const Vector3d & worldNormal,
		    const Vector2d * texPos);

  void drawTriangle(const Vector3d & a, const Vector3d & b, const Vector3d & c);

  //in Phong
  void drawLine(int px, const Vector3d & vp, const Vector3d & np, const Vector2d & tp,
		int qx, const Vector3d & vq, const Vector3d & nq, const Vector2d & tq,
		int h);
  //in Gouraud
  void drawLine(int px, double pz, const Vector3d & col_p,
		int qx, double qz, const Vector3d & col_q,
		int h);
  //in constant
  void drawLine(int px, double pz,
		int qx, double qz,
		const Vector3d & col,
		int h);
  void drawPoint(int x, int y, uchar r, uchar g, uchar b);
  void drawTriangle(int ax, int ay, uchar ar, uchar ag, uchar ab,
		    int bx, int by, uchar br, uchar bg, uchar bb,
		    int cx, int cy, uchar cr, uchar cg, uchar cb);
  void drawLine(int ax, int ay, int bx, int by,
		uchar ar, uchar ag, uchar ab, uchar br, uchar bg, uchar bb);
  
};



}
