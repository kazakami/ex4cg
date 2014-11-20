#ifdef _OPENMP
#include <omp.h>
#endif

#include "vrml.hpp"
#include "image.hpp"

namespace kazakami
{

Polygon::Polygon()
  :isTexed(false),
   tex(nullptr)
{
  mtr.diffuse[0] = 1.0;
  mtr.diffuse[1] = 1.0;
  mtr.diffuse[2] = 1.0;
  mtr.specular[0] = 0.0;
  mtr.specular[1] = 0.0;
  mtr.specular[2] = 0.0;
  mtr.ambient[0] = 0.5;
  mtr.ambient[1] = 0.5;
  mtr.ambient[2] = 0.5;
  mtr.shininess = 0.2;
}

bool Polygon::CheckIndex() const
{
  size_t max = pos.size();
  for (const face & f : faces)
  {
    if (f.a >= max)
	return false;
    if (f.b >= max)
	return false;
    if (f.c >= max)
	return false;
  }
  return true;
}

Vector2d Zero(0, 0);


void Polygon::AddVertex(Image3d * im)
{
  uint size = pos.size();
  for (uint i = 0; i < size; i++)
  {
    if (isTexed)
      im->AddVertex(pos[i], pos_normal[i], posTex[i]);
    else
      im->AddVertex(pos[i], pos_normal[i], Zero);
  }
}



void Polygon::SuperDraw(Image3d * im)
{
  im->SetMaterial(mtr);
  im->SetTexture(tex);
  uint size = faces.size();
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    im->CulcVertex();
#ifdef _OPENMP
#pragma omp for
#endif
    for (uint i = 0; i < size; i++)
    {
      const face & f = faces[i];
      im->drawTriangle_Phong(f.a, f.b, f.c);
    }
  }
}

void Polygon::Draw(Image3d * im)
{
  /*
  for (const face & f : faces)
    im->drawTriangle(pos[f.a], pos[f.b], pos[f.c], mtr);
  */
  im->SetMaterial(mtr);
  im->SetTexture(tex);
  uint size = faces.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (uint i = 0; i < size; i++)
    {
    const face & f = faces[i];
    //*
    if (isTexed)
    {
      const face & ft = facesTex[i];
      im->drawTriangle_Phong(pos[f.a], pos_normal[f.a], posTex[ft.a],
			     pos[f.b], pos_normal[f.b], posTex[ft.b],
			     pos[f.c], pos_normal[f.c], posTex[ft.c]);
    }
    else
    {
      im->drawTriangle_Phong(pos[f.a], pos_normal[f.a], Zero,
			     pos[f.b], pos_normal[f.b], Zero,
			     pos[f.c], pos_normal[f.c], Zero);
    }
    //*/
    /*
    im->drawTriangle_constant(pos[f.a],// pos_normal[f.a],
			      pos[f.b],// pos_normal[f.b],
			      pos[f.c],// pos_normal[f.c],
			      face_normal[i]);
    //*/
  }
}


void Polygon::AddPos(double x, double y, double z)
{
  pos.emplace_back(x, y, z);
}

void Polygon::AddFace(uint a, uint b, uint c)
{
  faces.emplace_back(a, b, c);
  face_normal.push_back(GetNormal(pos[a], pos[b], pos[c]));
}

void Polygon::AddPosTex(double x, double y)
{
  posTex.emplace_back(x, y);
}


void Polygon::AddFaceTex(uint a, uint b, uint c)
{
  facesTex.emplace_back(a, b, c);
}

std::vector<std::string> splitString(const std::string & str,
				     char delim, size_t start)
{
  std::vector<std::string> splits;
  std::string s;
  for (size_t i = start; i < str.size(); i++)
  {
    char c = str.at(i);
    if (c == delim)
    {
      if (!s.empty())
      {
	splits.push_back(s);
	s.clear();
      }
    }
    else
    {
      s += c;
    }
  }
  if (!s.empty())
    splits.push_back(s);
  return splits;
}

/*
void Polygon::ReadVRML_MT(const char * filename)
{
  std::vector<std::string> file_data;
  std::ifstream ifs(filename);
  std::string str;
  std::vector<std::string>::iterator 
  if (ifs.fail())
  {
    std::cerr << "fail to open " << filename << std::endl;
    exit(1);
  }
}
*/

void Polygon::ReadVRML(const char * filename) 
{
  std::vector<std::string> file_data;
  std::ifstream ifs(filename);
  std::string str;
  if (ifs.fail())
  {
    std::cerr << "failed to open " << filename << std::endl;
    exit(1);
  }
  while (getline(ifs, str))
  {
    file_data.push_back(str);
  }
  
  if (file_data.empty() || file_data[0] != "#VRML V2.0 utf8")
    std::cerr << "not a vrml v2.0 utf8 file" << std::endl;

  auto line_iterator = file_data.begin();
  while (line_iterator != file_data.end())
  {
    //shapeの読み取り開始
    if ((*line_iterator).find("Shape") != std::string::npos)
      while (line_iterator != file_data.end())
      {
	//appearance Appearanceの読み取り開始 in shape
	if ((*line_iterator).find("appearance Appearance") != std::string::npos)
	  while (line_iterator != file_data.end())
	  {
	    //material Materialの読み取り開始 in app in shape
	    if ((*line_iterator).find("material") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		//diffuseColorやで in material in app in shape
		size_t pos;
		if ((pos = (*line_iterator).find("diffuseColor"))
		    != std::string::npos)
		{
		  auto s = splitString((*line_iterator), ' ', pos + 12);
		  if (s.size() != 3)
		  {
		    std::cerr << "diffuse value is not correct" << std::endl;
		    exit(1);
		  }
		  for (int i = 0; i < 3; i++)
		    mtr.diffuse[i] = std::stod(s[i]);
		}//diffuseColor読み取り終了 in material in app in shape
		//specularColorの読み取り in material in app in shape
		if ((pos = (*line_iterator).find("specularColor"))
		    != std::string::npos)
		{
		  auto s = splitString((*line_iterator), ' ', pos + 13);
		  if (s.size() != 3)
		  {
		    std::cerr << "diffuse value is not correct" << std::endl;
		    exit(1);
		  }
		  for (int i = 0; i < 3; i++)
		    mtr.specular[i] = std::stod(s[i]);
		}//specularColor読み取り終了 in material in app in shape
		//shinessの読み取り in material in app in shape
		if ((pos = (*line_iterator).find("shininess"))
		    != std::string::npos)
		{
		  mtr.shininess = std::stod((*line_iterator).substr(pos + 9)) * 128;
		}//shiness読み取り終了 in material in app in shape
		//ambientの読み取り in material in app in shape
		if ((pos = (*line_iterator).find("ambientIntensity"))
		    != std::string::npos)
		{
		  auto d = std::stod((*line_iterator).substr(pos + 16));
		  for (int i = 0; i < 3; i++)
		    mtr.ambient[i] = d;
		}//ambient読み取り終了 in material in app in shape
		line_iterator++;
		if ((*line_iterator).find("}") != std::string::npos)
		  break;
	      }//material Materialの読み取り終了 in app in shape
	    //texture PixelTexture start in app in shape
	    if ((*line_iterator).find("texture PixelTexture") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		isTexed = true;
		size_t pos;
		if ((pos = (*line_iterator).find("image")) != std::string::npos)
		{
		  auto s = splitString((*line_iterator), ' ', pos + 5);
		  std::vector<uchar> d;
		  d.reserve(std::stoi(s[0]) * std::stoi(s[1]) * std::stoi(s[2]));
		  line_iterator++;
		  //start to load texture data
		  while (line_iterator != file_data.end())
		  {
		    if ((*line_iterator).find("}") != std::string::npos)
		    {
		      line_iterator--;
		      break;
		    }
		    auto ss = splitString((*line_iterator), ' ', 0);
		    for (const auto & str: ss)
		    {
		      uint ui = std::stoi(str, nullptr, 16);
		      uchar r = (ui & 0xFF0000) >> 16;
		      uchar g = (ui & 0x00FF00) >> 8;
		      uchar b = (ui & 0x0000FF);
		      //std::cout << "(" << ((int)r) << ", " << ((int)g) << ", " << ((int)b) << ")";
		      d.push_back(r);
		      d.push_back(g);
		      d.push_back(b);
		    }
		    //std::cout << std::endl;
		    line_iterator++;
		  }//end loading tex
		  tex = std::make_shared<Texture>();
		  tex->LoadFromVector(std::move(d), std::stoi(s[0]), std::stoi(s[1]));
		}
		line_iterator++;
		if ((*line_iterator).find("}") != std::string::npos)
		  break;
	      }//texture PixelTexture end in app in shape
	    line_iterator++;
	    if ((*line_iterator).find("}") != std::string::npos)
	      break;
	  }//appearance Appearanceの読み取り終了 in shape
	//geometry IndexedFaceSetの読み取り開始 in shape
	if ((*line_iterator).find("geometry IndexedFaceSet") != std::string::npos)
	  while (line_iterator != file_data.end())
	  {
	    //std::cout << *line_iterator << std::endl;
	    //texCoord TextureCoordinate start to load
	    if ((*line_iterator).find("texCoord TextureCoordinate") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		//pointの読み取り開始 in Tex coord in geo in shape
		if ((*line_iterator).find("point") != std::string::npos)
		  while (line_iterator != file_data.end())
		  {
		    line_iterator++;
		    if ((*line_iterator).find("]") != std::string::npos)
		      break;
		    else
		    {
		      auto s = splitString((*line_iterator), ' ', 0);
		      if (s.size() != 2)
		      {
			std::cerr << "tex coord point is not 2d" << std::endl;
			exit(1);
		      }
		      AddPosTex(std::stod(s[0]), std::stod(s[1]));
		      //std::cout << std::stod(s[0]) << std::stod(s[1]) << std::endl;
		    }
		  }//pointの読み取り終了 in Tex coord in geo in shape
		line_iterator++;
		if ((*line_iterator).find("}") != std::string::npos)
		  break;
	      }//texCoord TextureCoordinate end loading
	    //coord Coordinateの読み取り開始 in geo in shape
	    if ((*line_iterator).find("coord Coordinate") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		//pointの読み取り開始 in coord in geo in shape
		if ((*line_iterator).find("point") != std::string::npos)
		  while (line_iterator != file_data.end())
		  {
		    line_iterator++;
		    if ((*line_iterator).find("]") != std::string::npos)
		      break;
		    else
		    {
		      auto s = splitString((*line_iterator), ' ', 0);
		      if (s.size() != 3)
		      {
			std::cerr << "point is not 3d" << std::endl;
			exit(1);
		      }
		      AddPos(std::stod(s[0]),
			     std::stod(s[1]),
			     std::stod(s[2]));
		    }
		  }//pointの読み取り終了 in coord in geo in shape
		line_iterator++;
		if ((*line_iterator).find("}") != std::string::npos)
		  break;
	      }//coord Coordinateの読み取り終了 in geo in shape
	    //coordIndexの読み取り開始 in geo in shape
	    if ((*line_iterator).find("coordIndex") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		line_iterator++;
		if ((*line_iterator).find("]") != std::string::npos)
		  break;
		else
		{
		  auto s = splitString((*line_iterator), ',', 0);
		  if (s.size() != 4 && std::stoi(s[3]) != -1)
		  {
		    std::cerr << "incorect coordIndex" << std::endl;
		    exit(1);
		  }
		  AddFace(std::stoi(s[0]),
			  std::stoi(s[1]),
			  std::stoi(s[2]));
		}
	      }//coordIndexの読み取り終了 in geo in shape
	    //texCoordIndexの読み取り開始 in geo in shape
	    if ((*line_iterator).find("texCoordIndex") != std::string::npos)
	      while (line_iterator != file_data.end())
	      {
		line_iterator++;
		if ((*line_iterator).find("]") != std::string::npos)
		  break;
		else
		{
		  auto s = splitString((*line_iterator), ',', 0);
		  if (s.size() != 4 && std::stoi(s[3]) != -1)
		  {
		    std::cerr << "incorect texCoordIndex" << std::endl;
		    exit(1);
		  }
		  AddFaceTex(std::stoi(s[0]),
			     std::stoi(s[1]),
			     std::stoi(s[2]));
		}
	      }//texCoordIndexの読み取り終了 in geo in shape
	    line_iterator++;
	    if ((*line_iterator).find("}") != std::string::npos)
	      break;
	  }//geometry IndexedFaceSetの読み取り終了 in shape
	line_iterator++;
	if ((*line_iterator).find("}") != std::string::npos)
	  break;
      }//Shapeの読み取り終了 at root
    line_iterator++;
  }
  if (!CheckIndex())
  {
    std::cerr << "incorrect vertex index" << std::endl;
    exit(1);
  }
  std::cerr << "VRML loaded: face_num = " << faces.size() << std::endl;
  std::cerr << "ft_num = " << facesTex.size() << std::endl;
  std::cerr << "pos_num = " << pos.size() << std::endl;
  std::cerr << "pt_num = " << posTex.size() << std::endl;
  
  int max = faces.size();
  pos_normal.resize(pos.size());
  pos_belong_to.resize(pos.size());
  for (int i = 0; i < max; i++)
  {
    const face & f = faces[i];
    face_normal.push_back(GetNormal(pos[f.a],
				    pos[f.b],
				    pos[f.c]));
    pos_belong_to[f.a].push_back(i);
    pos_belong_to[f.b].push_back(i);
    pos_belong_to[f.c].push_back(i);
  }
  max = pos.size();
  for (int i = 0; i < max; i++)
  {
    Vector3d normal;
    for (int j : pos_belong_to[i])
      normal += face_normal[j];
    pos_normal[i] = Normalise(normal);
  }
}


}
