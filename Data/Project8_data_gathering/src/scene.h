//-------------------------------------------------------------------------------
///
/// \file       scene.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    8.0
/// \date       October 16, 2017
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _SCENE_H_INCLUDED_
#define _SCENE_H_INCLUDED_

//-------------------------------------------------------------------------------

#define TEXTURE_SAMPLE_COUNT 16

//-------------------------------------------------------------------------------

//RAY DIFFERENTIALS
//comment this out if you don't want ray differentials
#define RAY_DIFFERENTIALS
//we do, in fact, need this #define in this file too
//https://stackoverflow.com/a/16445543/752878

//comment this out if you don't want to collect ray data
#define COLLECT_DATA

#include <string.h>
#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>
#include <atomic>

#include "lodepng.h"

#include "cyPoint.h"
typedef cyPoint2f Point2;
typedef cyPoint3f Point3;
typedef cyPoint4f Point4;

#include "cyMatrix.h"
typedef cyMatrix3f Matrix3;

#include "cyColor.h"
typedef cyColor Color;
typedef cyColorA ColorA;
typedef cyColor24 Color24;

typedef unsigned char uchar;

//-------------------------------------------------------------------------------

#ifndef min
# define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
# define max(a,b) ((a)>(b)?(a):(b))
#endif

#define BIGFLOAT 1.0e30f

//-------------------------------------------------------------------------------

/*
For exporting trace data to JSON
*/
struct RayInfo {
	Point3 Origin; //Origin
	Point3 Direction; //Direction
	Point3 HitPoint;	//Hit point
	//TODO				  //Object hit
	unsigned int numRayBoxIntersections; //Number of ray-box intersection tests
	unsigned int numRayPrimitiveIntersections;	//Number of ray-primitive intersection tests
	int rayType;	//Ray type (0: camera, 1: shadow, 2: reflection, 3: refraction)
	cy::Color colorContribution;	//bounce color contribution

	std::vector<RayInfo*> shadowRays;	//vector of shadow ray info(s)
	RayInfo *refractionRay;	//pointer to refraction ray info
	RayInfo	*reflectionRay;	//pointer to reflection ray info

	RayInfo() { Init(); }
	void Init(){
		numRayBoxIntersections = 0;
		numRayPrimitiveIntersections = 0;
		rayType = -1;
		//TODO: initialize the vector of shadowRays?
		refractionRay = NULL;
		reflectionRay = NULL;

		Origin = Point3();
		Direction = Point3();
		HitPoint = Point3();
	}

	int numSecondaryRays() {
		int result = 0;
		//reflection ray
		if (reflectionRay != NULL) {
			result++;
			//its secondary rays
			result += reflectionRay->numSecondaryRays();
		}
		//refraction ray
		if (refractionRay != NULL) {
			result++;
			//its secondary rays
			result += refractionRay->numSecondaryRays();
		}
		//shadow rays
		result += shadowRays.size();

		return result;
	}

	std::string toJSON() {
		//o, d, p, ob, rb, rp, t, c, ch
		std::string json = "{";
		//json += "\"o\":" + pointToJSON(Origin);
		//json += ",\"d\":" + pointToJSON(Direction);
		//json += ",\"p\":" + pointToJSON(HitPoint);
		//json += ",\"ob\":0";	//TODO: object index
		json += "\"rb\":" + std::to_string(numRayBoxIntersections);
		json += ",\"rp\":" + std::to_string(numRayPrimitiveIntersections);
		json += ",\"t\":" + std::to_string(rayType);
		json += ",\"c\":" + colorToJSON(colorContribution);

		json += ",\"ch\":[";
		int numShadows = shadowRays.size();
		for (int i = 0; i < numShadows; i++) {
			if (i != 0) {
				json += ",";
			}
			json += shadowRays[i]->toJSON();
		}
		if (refractionRay != NULL) {
			json += "," + refractionRay->toJSON();
		}
		if (reflectionRay != NULL) {
			json += "," + reflectionRay->toJSON();
		}
		json += "]";

		json += "}";
		return json;
	}
	//for serializing the Point3 to JSON
	/*std::string toJSON() {
		std::string result("[");
		result += std::to_string(Element(0)) + ",";
		result += std::to_string(Element(1)) + ",";
		result += std::to_string(Element(2)) + "]";
		return result;

		return std::string("["
			+ std::to_string(x) + ","
			+ std::to_string(y) + ","
			+ std::to_string(z) + "]");
	}*/
	std::string pointToJSON(Point3 p) {
		return std::string("["
			+ std::to_string(p.x) + ","
			+ std::to_string(p.y) + ","
			+ std::to_string(p.z) + "]");
	}
	std::string colorToJSON(Color c) {
		return std::string("["
			+ std::to_string(c.r) + ","
			+ std::to_string(c.g) + ","
			+ std::to_string(c.b) + "]");
	}
};
//BRAIN DUMP:
//add this struct as a member of the Ray class
//populate this struct as the ray traverses the scene
//Write a recursive function which builds a JSON string for a RayInfo struct and its children
//When loading the scene, give each object an integer index so we can store that in the struct, instead of storing a name or something
//We can't build the RayInfos for shadow rays directly, but we can assume enough data to build a reasonable guess:
//		Without the actual ray, we can get the start point, the direction (from the light's position), and the target endpoint, and whether or not it reached the light

struct PixelInfo {
	//Final pixel color
	//Time taken to render
	double renderTime;
	//Sample count
	unsigned int sampleCount;
	//Depth buffer?
	//Variance of the samples
	float sampleVariance;
	//Number of secondary rays
	unsigned int numSecondaryRays;

	unsigned int numRayBoxIntersections;

	unsigned int numRayPrimitiveIntersections;

	//Vector of RayInfos
	std::vector<RayInfo*> samples;

	PixelInfo() { Init(); }
	void Init() {
		renderTime = 0.f;
		sampleCount = 0;
		numSecondaryRays = 0;
		sampleVariance = 0;
		//TODO: Initialize the vector of samples?
	}
};
//BRAIN DUMP:
//We need to add a map of these to RenderImage maybe, one for each pixel


class Ray
{
public:
	Point3 p, dir;
#ifdef RAY_DIFFERENTIALS
	Point3 p_u, dir_u, p_v, dir_v;	//the ray origins and directions for the U and V differential rays
	Ray(const Point3 &_p, const Point3 &_dir, const Point3 &_p_u, const Point3 &_dir_u, const Point3 &_p_v, const Point3 &_dir_v) : p(_p), dir(_dir), p_u(_p_u), dir_u(_dir_u), p_v(_p_v), dir_v(_dir_v) {}
	Ray(const Point3 &_p, const Point3 &_dir) : p(_p), dir(_dir) {}
	Ray(const Ray &r) : p(r.p), dir(r.dir), p_u(r.p_u), dir_u(r.dir_u), p_v(r.p_v), dir_v(r.dir_v) {}
	void Normalize() { dir.Normalize(); dir_u.Normalize(); dir_v.Normalize(); }
#else
	Ray(const Point3 &_p, const Point3 &_dir) : p(_p), dir(_dir) {}
	Ray(const Ray &r) : p(r.p), dir(r.dir) {}
	void Normalize() { dir.Normalize(); }
#endif

	Ray() {}
};

//-------------------------------------------------------------------------------

class Box
{
public:
	Point3 pmin, pmax;

	// Constructors
	Box() { Init(); }
	Box(const Point3 &_pmin, const Point3 &_pmax) : pmin(_pmin), pmax(_pmax) {}
	Box(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax ) : pmin(xmin,ymin,zmin), pmax(xmax,ymax,zmax) {}
	Box(const float *dim) : pmin(dim), pmax(&dim[3]) {}

	// Initializes the box, such that there exists no point inside the box (i.e. it is empty).
	void Init() { pmin.Set(BIGFLOAT,BIGFLOAT,BIGFLOAT); pmax.Set(-BIGFLOAT,-BIGFLOAT,-BIGFLOAT); }

	// Returns true if the box is empty; otherwise, returns false.
	bool IsEmpty() const { return pmin.x>pmax.x || pmin.y>pmax.y || pmin.z>pmax.z; }

	// Returns one of the 8 corner point of the box in the following order:
	// 0:(x_min,y_min,z_min), 1:(x_max,y_min,z_min)
	// 2:(x_min,y_max,z_min), 3:(x_max,y_max,z_min)
	// 4:(x_min,y_min,z_max), 5:(x_max,y_min,z_max)
	// 6:(x_min,y_max,z_max), 7:(x_max,y_max,z_max)
	Point3 Corner( int i ) const	// 8 corners of the box
	{
		Point3 p;
		p.x = (i & 1) ? pmax.x : pmin.x;
		p.y = (i & 2) ? pmax.y : pmin.y;
		p.z = (i & 4) ? pmax.z : pmin.z;
		return p;
	}

	// Enlarges the box such that it includes the given point p.
	void operator += (const Point3 &p)
	{
		for ( int i=0; i<3; i++ ) {
			if ( pmin[i] > p[i] ) pmin[i] = p[i];
			if ( pmax[i] < p[i] ) pmax[i] = p[i];
		}
	}

	// Enlarges the box such that it includes the given box b.
	void operator += (const Box &b)
	{
		for ( int i=0; i<3; i++ ) {
			if ( pmin[i] > b.pmin[i] ) pmin[i] = b.pmin[i];
			if ( pmax[i] < b.pmax[i] ) pmax[i] = b.pmax[i];
		}
	}

	// Returns true if the point is inside the box; otherwise, returns false.
	bool IsInside(const Point3 &p) const { for ( int i=0; i<3; i++ ) if ( pmin[i] > p[i] || pmax[i] < p[i] ) return false; return true; }

	// Returns true if the ray intersects with the box for any parameter that is smaller than t_max; otherwise, returns false.
	bool IntersectRay(const Ray &r, float t_max, float &tval) const;
};

//-------------------------------------------------------------------------------

inline float Halton(int index, int base)
{
	float r = 0;
	float f = 1.0f / (float)base;
	for ( int i=index; i>0; i/=base ) {
		r += f * (i%base);
		f /= (float) base;
	}
	return r;
}

//-------------------------------------------------------------------------------

class Node;

#define HIT_NONE			0
#define HIT_FRONT			1
#define HIT_BACK			2
#define HIT_FRONT_AND_BACK	(HIT_FRONT|HIT_BACK)

struct HitInfo
{
	float z;			// the distance from the ray center to the hit point
	Point3 p;			// position of the hit point
	Point3 N;			// surface normal at the hit point
	Point3 uvw;			// texture coordinate at the hit point
	Point3 duvw[2];		// derivatives of the texture coordinate
	const Node *node;	// the object node that was hit
	bool front;			// true if the ray hits the front side, false if the ray hits the back side
	int mtlID;			// sub-material index

#ifdef RAY_DIFFERENTIALS
	Point3 p_u;			// positions of the differential hit points
	Point3 p_v;
#endif
	HitInfo() { Init(); }
	void Init() { z=BIGFLOAT; node=NULL; front=true; uvw.Set(0.5f,0.5f,0.5f); duvw[0].Zero(); duvw[1].Zero(); mtlID=0; }
};

//-------------------------------------------------------------------------------

class ItemBase
{
private:
	char *name;					// The name of the item

public:
	ItemBase() : name(NULL) {}
	virtual ~ItemBase() { if ( name ) delete [] name; }

	const char* GetName() const { return name ? name : ""; }
	void SetName(const char *newName)
	{
		if ( name ) delete [] name;
		if ( newName ) {
			int n = strlen(newName);
			name = new char[n+1];
			for ( int i=0; i<n; i++ ) name[i] = newName[i];
			name[n] = '\0';
		} else { name = NULL; }
	}
};

template <class T> class ItemList : public std::vector<T*>
{
public:
	virtual ~ItemList() { DeleteAll(); }
	void DeleteAll() { int n=(int)this->size(); for ( int i=0; i<n; i++ ) if ( this->at(i) ) delete this->at(i); }
};


template <class T> class ItemFileList
{
public:
	void Clear() { list.DeleteAll(); }
	void Append( T* item, const char *name ) { list.push_back( new FileInfo(item,name) ); }
	T* Find( const char *name ) const { int n=list.size(); for ( int i=0; i<n; i++ ) if ( list[i] && strcmp(name,list[i]->GetName())==0 ) return list[i]->GetObj(); return NULL; }

private:
	class FileInfo : public ItemBase
	{
	private:
		T *item;
	public:
		FileInfo() : item(NULL) {}
		FileInfo(T *_item, const char *name) : item(_item) { SetName(name); }
		~FileInfo() { Delete(); }
		void Delete() { if (item) delete item; item=NULL; }
		void SetObj(T *_item) { Delete(); item=_item; }
		T* GetObj() { return item; }
	};

	ItemList<FileInfo> list;
};

//-------------------------------------------------------------------------------

class Transformation
{
private:
	Matrix3 tm;						// Transformation matrix to the local space
	Point3 pos;						// Translation part of the transformation matrix
	mutable Matrix3 itm;			// Inverse of the transformation matrix (cached)
public:
	Transformation() : pos(0,0,0) { tm.SetIdentity(); itm.SetIdentity(); }
	const Matrix3& GetTransform() const { return tm; }
	const Point3& GetPosition() const { return pos; }
	const Matrix3&	GetInverseTransform() const { return itm; }

	Point3 TransformTo( const Point3 &p ) const { return itm * (p - pos); }	// Transform to the local coordinate system
	Point3 TransformFrom( const Point3 &p ) const { return tm*p + pos; }	// Transform from the local coordinate system

	// Transforms a vector to the local coordinate system (same as multiplication with the inverse transpose of the transformation)
	Point3 VectorTransformTo( const Point3 &dir ) const { return TransposeMult(tm,dir); }

	// Transforms a vector from the local coordinate system (same as multiplication with the inverse transpose of the transformation)
	Point3 VectorTransformFrom( const Point3 &dir ) const { return TransposeMult(itm,dir); }

	void Translate(Point3 p) { pos+=p; }
	void Rotate(Point3 axis, float degree) { Matrix3 m; m.SetRotation(axis,degree*(float)M_PI/180.0f); Transform(m); }
	void Scale(float sx, float sy, float sz) { Matrix3 m; m.Zero(); m[0]=sx; m[4]=sy; m[8]=sz; Transform(m); }
	void Transform(const Matrix3 &m) { tm=m*tm; pos=m*pos; tm.GetInverse(itm); }

	void InitTransform() { pos.Zero(); tm.SetIdentity(); itm.SetIdentity(); }

private:
	// Multiplies the given vector with the transpose of the given matrix
	static Point3 TransposeMult( const Matrix3 &m, const Point3 &dir )
	{
		Point3 d;
		d.x = m.GetColumn(0) % dir;
		d.y = m.GetColumn(1) % dir;
		d.z = m.GetColumn(2) % dir;
		return d;
	}
};

//-------------------------------------------------------------------------------

class Material;

// Base class for all object types
class Object
{
public:
	virtual bool IntersectRay( const Ray &ray, HitInfo &hInfo, RayInfo *_rInfo, int hitSide=HIT_FRONT ) const=0;
	virtual Box  GetBoundBox() const=0;
	virtual void ViewportDisplay(const Material *mtl) const {}	// used for OpenGL display
};

typedef ItemFileList<Object> ObjFileList;

//-------------------------------------------------------------------------------

class Light : public ItemBase
{
public:
	virtual Color	Illuminate(const Point3 &p, const Point3 &N) const=0;
	virtual Point3	Direction (const Point3 &p) const=0;
	virtual bool	IsAmbient () const { return false; }
	virtual void	SetViewportLight(int lightID) const {}	// used for OpenGL display
};

class LightList : public ItemList<Light> {};

//-------------------------------------------------------------------------------

class Material : public ItemBase
{
public:
	// The main method that handles the shading by calling all the lights in the list.
	// ray: incoming ray,
	// hInfo: hit information for the point that is being shaded, lights: the light list,
	// bounceCount: permitted number of additional bounces for reflection and refraction.
	virtual Color Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, RayInfo *_rInfo) const=0;


	virtual void SetViewportMaterial(int subMtlID=0) const {}	// used for OpenGL display
};

class MaterialList : public ItemList<Material>
{
public:
	Material* Find( const char *name ) { int n=size(); for ( int i=0; i<n; i++ ) if ( at(i) && strcmp(name,at(i)->GetName())==0 ) return at(i); return NULL; }
};

//-------------------------------------------------------------------------------

class Texture : public ItemBase
{
public:
	// Evaluates the color at the given uvw location.
	virtual Color Sample(const Point3 &uvw) const=0;

	// Evaluates the color around the given uvw location using the derivatives duvw
	// by calling the Sample function multiple times.
	virtual Color Sample(const Point3 &uvw, const Point3 duvw[2], bool elliptic=true) const
	{
		Color c = Sample(uvw);
		if ( duvw[0].LengthSquared() + duvw[1].LengthSquared() == 0 ) return c;
		for ( int i=1; i<TEXTURE_SAMPLE_COUNT; i++ ) {
			float x = Halton(i,2);
			float y = Halton(i,3);
			if ( elliptic ) {
				float r = sqrtf(x)*0.5f;
				x = r*sinf(y*(float)M_PI*2);
				y = r*cosf(y*(float)M_PI*2);
			} else {
				if ( x > 0.5f ) x-=1;
				if ( y > 0.5f ) y-=1;
			}
			c += Sample( uvw + x*duvw[0] + y*duvw[1] );
		}
		return c / float(TEXTURE_SAMPLE_COUNT);
	}

	virtual bool SetViewportTexture() const { return false; }	// used for OpenGL display

protected:

	// Clamps the uvw values for tiling textures, such that all values fall between 0 and 1.
	static Point3 TileClamp(const Point3 &uvw)
	{
		Point3 u;
		u.x = uvw.x - (int) uvw.x;
		u.y = uvw.y - (int) uvw.y;
		u.z = uvw.z - (int) uvw.z;
		if ( u.x < 0 ) u.x += 1;
		if ( u.y < 0 ) u.y += 1;
		if ( u.z < 0 ) u.z += 1;
		return u;
	}
};

typedef ItemFileList<Texture> TextureList;

//-------------------------------------------------------------------------------

// This class handles textures with texture transformations.
// The uvw values passed to the Sample methods are transformed
// using the texture transformation.
class TextureMap : public Transformation
{
public:
	TextureMap() : texture(NULL) {}
	TextureMap(Texture *tex) : texture(tex) {}
	void SetTexture(Texture *tex) { texture = tex; }

	virtual Color Sample(const Point3 &uvw) const { return texture ? texture->Sample(TransformTo(uvw)) : Color(0,0,0); }
	virtual Color Sample(const Point3 &uvw, const Point3 duvw[2], bool elliptic=true) const
	{
		if ( texture == NULL ) return Color(0,0,0);
		Point3 u = TransformTo(uvw);
		Point3 d[2];
		d[0] = TransformTo(duvw[0]+uvw)-u;
		d[1] = TransformTo(duvw[1]+uvw)-u;
		return texture->Sample(u,d,elliptic);
	}

	bool SetViewportTexture() const { if ( texture ) return texture->SetViewportTexture(); return false; }	// used for OpenGL display

private:
	Texture *texture;
};

//-------------------------------------------------------------------------------

// This class keeps a TextureMap and a color. This is useful for keeping material
// color parameters that can also be textures. If no texture is specified, it
// automatically uses the color value. Otherwise, the texture value is multiplied
// by the color value.
class TexturedColor
{
private:
	Color color;
	TextureMap *map;
public:
	TexturedColor() : color(0,0,0), map(NULL) {}
	TexturedColor(float r, float g, float b) : color(r,g,b), map(NULL) {}
	virtual ~TexturedColor() { if ( map ) delete map; }

	void SetColor(const Color &c) { color=c; }
	void SetTexture(TextureMap *m) { if ( map ) delete map; map=m; }

	Color GetColor() const { return color; }
	const TextureMap* GetTexture() const { return map; }

	Color Sample(const Point3 &uvw) const { return ( map ) ? color*map->Sample(uvw) : color; }
	Color Sample(const Point3 &uvw, const Point3 duvw[2], bool elliptic=true) const { return ( map ) ? color*map->Sample(uvw,duvw,elliptic) : color; }

	// Returns the color value at the given direction for environment mapping.
	Color SampleEnvironment(const Point3 &dir) const
	{ 
		float z = asinf(-dir.z)/float(M_PI)+0.5f;
		float x = dir.x / (fabs(dir.x)+fabs(dir.y));
		float y = dir.y / (fabs(dir.x)+fabs(dir.y));
		return Sample( Point3(0.5f,0.5f,0.0f) + z*(x*Point3(0.5f,0.5f,0) + y*Point3(-0.5f,0.5f,0)) );
	}

};

//-------------------------------------------------------------------------------

class Node : public ItemBase, public Transformation
{
private:
	Node **child;				// Child nodes
	int numChild;				// The number of child nodes
	Object *obj;				// Object reference (merely points to the object, but does not own the object, so it doesn't get deleted automatically)
	Material *mtl;				// Material used for shading the object
	Box childBoundBox;			// Bounding box of the child nodes, which does not include the object of this node, but includes the objects of the child nodes
public:
	Node() : child(NULL), numChild(0), obj(NULL), mtl(NULL) {}
	virtual ~Node() { DeleteAllChildNodes(); }

	void Init() { DeleteAllChildNodes(); obj=NULL; mtl=NULL; childBoundBox.Init(); SetName(NULL); InitTransform(); } // Initialize the node deleting all child nodes

	// Hierarchy management
	int	 GetNumChild() const { return numChild; }
	void SetNumChild(int n, int keepOld=false)
	{
		if ( n < 0 ) n=0;	// just to be sure
		Node **nc = NULL;	// new child pointer
		if ( n > 0 ) nc = new Node*[n];
		for ( int i=0; i<n; i++ ) nc[i] = NULL;
		if ( keepOld ) {
			int sn = min(n,numChild);
			for ( int i=0; i<sn; i++ ) nc[i] = child[i];
		}
		if ( child ) delete [] child;
		child = nc;
		numChild = n;
	}
	const Node*	GetChild(int i) const		{ return child[i]; }
	Node*		GetChild(int i)				{ return child[i]; }
	void		SetChild(int i, Node *node)	{ child[i]=node; }
	void		AppendChild(Node *node)		{ SetNumChild(numChild+1,true); SetChild(numChild-1,node); }
	void		RemoveChild(int i)			{ for ( int j=i; j<numChild-1; j++) child[j]=child[j+1]; SetNumChild(numChild-1); }
	void		DeleteAllChildNodes()		{ for ( int i=0; i<numChild; i++ ) { child[i]->DeleteAllChildNodes(); delete child[i]; } SetNumChild(0); }

	// Bounding Box
	const Box& ComputeChildBoundBox()
	{
		childBoundBox.Init();
		for ( int i=0; i<numChild; i++ ) {
			Box childBox = child[i]->ComputeChildBoundBox();
			Object *cobj = child[i]->GetNodeObj();
			if ( cobj ) childBox += cobj->GetBoundBox();
			if ( ! childBox.IsEmpty() ) {
				// transform the box from child coordinates
				for ( int j=0; j<8; j++ ) childBoundBox += child[i]->TransformFrom( childBox.Corner(j) );
			}
		}
		return childBoundBox;
	}
	const Box& GetChildBoundBox() const { return childBoundBox; }

	// Object management
	const Object*	GetNodeObj() const { return obj; }
	Object*			GetNodeObj() { return obj; }
	void			SetNodeObj(Object *object) { obj=object; }

	// Material management
	const Material* GetMaterial() const { return mtl; }
	void			SetMaterial(Material *material) { mtl=material; }

	// Transformations
	Ray ToNodeCoords( const Ray &ray ) const
	{
		Ray r;
		r.p   = TransformTo(ray.p);
		r.dir = TransformTo(ray.p + ray.dir) - r.p;
#ifdef RAY_DIFFERENTIALS
		r.p_u = TransformTo(ray.p_u);
		r.dir_u = TransformTo(ray.p_u + ray.dir_u) - r.p_u;
		r.p_v = TransformTo(ray.p_v);
		r.dir_v = TransformTo(ray.p_v + ray.dir_v) - r.p_v;
#endif
		return r;
	}
	void FromNodeCoords( HitInfo &hInfo ) const
	{
		hInfo.p = TransformFrom(hInfo.p);
		hInfo.N = VectorTransformFrom(hInfo.N).GetNormalized();
#ifdef RAY_DIFFERENTIALS
		hInfo.p_u = TransformFrom(hInfo.p_u);
		hInfo.p_v = TransformFrom(hInfo.p_v);
#endif
	}
};

//-------------------------------------------------------------------------------

class Camera
{
public:
	Point3 pos, dir, up;
	float fov;
	int imgWidth, imgHeight;

	void Init()
	{
		pos.Set(0,0,0);
		dir.Set(0,0,-1);
		up.Set(0,1,0);
		fov = 40;
		imgWidth = 200;
		imgHeight = 150;
	}
};

//-------------------------------------------------------------------------------

class RenderImage
{
private:
	Color24	*img;
	float	*zbuffer;
	uchar	*zbufferImg;
	uchar	*sampleCount;
	uchar	*sampleCountImg;
	int		width, height;
	std::atomic<int> numRenderedPixels;
public:
//TODO: ray data here
	PixelInfo **pixInfo;
	RenderImage() : img(NULL), zbuffer(NULL), zbufferImg(NULL), width(0), height(0), numRenderedPixels(0), pixInfo(0) {}
	void Init(int w, int h)
	{
		width=w;
		height=h;
		if (img) delete [] img;
		img = new Color24[width*height];
		if (zbuffer) delete [] zbuffer;
		zbuffer = new float[width*height];
		if (zbufferImg) delete [] zbufferImg;
		zbufferImg = NULL;
		if ( sampleCount ) delete [] sampleCount;
		sampleCount = new uchar[width*height];;
		if ( sampleCountImg ) delete [] sampleCountImg;
		sampleCountImg = NULL;
		if (pixInfo) delete[] pixInfo;
		pixInfo = new PixelInfo*[width*height];
		ResetNumRenderedPixels();
	}

	int			GetWidth() const	{ return width; }
	int			GetHeight() const	{ return height; }
	Color24*	GetPixels()			{ return img; }
	float*		GetZBuffer()		{ return zbuffer; }
	uchar*		GetZBufferImage()	{ return zbufferImg; }
	uchar*		GetSampleCount()	{ return sampleCount; }
	uchar*		GetSampleCountImage(){ return sampleCountImg; }
	PixelInfo**	GetPixelInfoMap() { return pixInfo; }

	void	ResetNumRenderedPixels()		{ numRenderedPixels=0; }
	int		GetNumRenderedPixels() const	{ return numRenderedPixels; }
	void	IncrementNumRenderPixel(int n)	{ numRenderedPixels+=n; }
	bool	IsRenderDone() const			{ return numRenderedPixels >= width*height; }

	void	ComputeZBufferImage()
	{
		int size = width * height;
		if (zbufferImg) delete [] zbufferImg;
		zbufferImg = new uchar[size];

		float zmin=BIGFLOAT, zmax=0;
		for ( int i=0; i<size; i++ ) {
			if ( zbuffer[i] == BIGFLOAT ) continue;
			if ( zmin > zbuffer[i] ) zmin = zbuffer[i];
			if ( zmax < zbuffer[i] ) zmax = zbuffer[i];
		}
		for ( int i=0; i<size; i++ ) {
			if ( zbuffer[i] == BIGFLOAT ) zbufferImg[i] = 0;
			else {
				float f = (zmax-zbuffer[i])/(zmax-zmin);
				int c = int(f * 255);
				if ( c < 0 ) c = 0;
				if ( c > 255 ) c = 255;
				zbufferImg[i] = c;
			}
		}
	}

	int ComputeSampleCountImage()
	{
		int size = width * height;
		if (sampleCountImg) delete [] sampleCountImg;
		sampleCountImg = new uchar[size];

		uchar smin=255, smax=0;
		for ( int i=0; i<size; i++ ) {
			if ( smin > sampleCount[i] ) smin = sampleCount[i];
			if ( smax < sampleCount[i] ) smax = sampleCount[i];
		}
		if ( smax == smin ) {
			for ( int i=0; i<size; i++ ) sampleCountImg[i] = 0;
		} else {
			for ( int i=0; i<size; i++ ) {
				int c = (255*(sampleCount[i]-smin))/(smax-smin);
				if ( c < 0 ) c = 0;
				if ( c > 255 ) c = 255;
				sampleCountImg[i] = c;
			}
		}
		return smax;
	}

	bool SaveImage (const char *filename) const { return SavePNG(filename,&img[0].r,3); }
	bool SaveZImage(const char *filename) const { return SavePNG(filename,zbufferImg,1); }
	bool SaveSampleCountImage(const char *filename) const { return SavePNG(filename,sampleCountImg,1); }

	bool SavePixInfo(const char *dirname) const {
		std::string SceneJSON, PixelDataJSON, RayDataJSON;

		//Scene TODO
		//PixelDataJSON += "\"Scene\":{}";
			//camera
			//objects
		
		
		PixelDataJSON = "{";
		//PixelData
		PixelDataJSON += "\"PixelData\":{";
		//dimensions [width,height]
		PixelDataJSON += "\"dimensions\":["+ std::to_string(width) +","+ std::to_string(height) +"],";
		//Color [r,g,b, r,g,b, r,g,b,...]
		PixelDataJSON += "\"Color\":[";
		unsigned int idx = 0;

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				//these are 8-bit chars
				PixelDataJSON += std::to_string(img[idx].r);
				PixelDataJSON += "," + std::to_string(img[idx].g);
				PixelDataJSON += "," + std::to_string(img[idx].b);
				idx++;
			}
		}
		PixelDataJSON += "],";
		//render time
		idx = 0;
		PixelDataJSON += "\"Render_time\":[";
		double minRenderTime = BIGFLOAT;
		double maxRenderTime = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(pixInfo[idx]->renderTime);
				minRenderTime = min(minRenderTime, pixInfo[idx]->renderTime);
				maxRenderTime = max(maxRenderTime, pixInfo[idx]->renderTime);
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max render time
		PixelDataJSON += "\"min_render_time\":" + std::to_string(minRenderTime) + ",";
		PixelDataJSON += "\"max_render_time\":" + std::to_string(maxRenderTime) + ",";

		/* render time image */
		std::vector<uchar> renderTimeImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				renderTimeImage[x + width * y] = ((pixInfo[x + width * y]->renderTime - minRenderTime) / (maxRenderTime - minRenderTime)) * 255;
			}
		}
		std::cout << "renderTime.png " << SavePNG("renderTime.png", renderTimeImage.data(), 1) << std::endl;;


		//secondary rays
		idx = 0;
		PixelDataJSON += "\"Secondary_rays\":[";
		int minSecondaryRays = INT32_MAX;
		int maxSecondaryRays = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(pixInfo[idx]->numSecondaryRays);
				minSecondaryRays = min(minSecondaryRays, pixInfo[idx]->numSecondaryRays);
				maxSecondaryRays = max(maxSecondaryRays, pixInfo[idx]->numSecondaryRays);
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max secondary rays
		PixelDataJSON += "\"min_secondary_rays\":" + std::to_string(minSecondaryRays) + ",";
		PixelDataJSON += "\"max_secondary_rays\":" + std::to_string(maxSecondaryRays) + ",";

		/* Secondary rays image */
		std::vector<uchar> secondaryRaysImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				secondaryRaysImage[x + width * y] = ((pixInfo[x + width * y]->numSecondaryRays - minSecondaryRays) / (float)(maxSecondaryRays - minSecondaryRays)) * 255;
			}
		}
		std::cout << "secondaryRays.png " << SavePNG("secondaryRays.png", secondaryRaysImage.data(), 1) << std::endl;

		//sample count
		idx = 0;
		PixelDataJSON += "\"Sample_count\":[";
		int minSampleCount = INT32_MAX;
		int maxSampleCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(pixInfo[idx]->samples.size());
				minSampleCount = min(minSampleCount, pixInfo[idx]->samples.size());
				maxSampleCount = max(maxSampleCount, pixInfo[idx]->samples.size());
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max sample count
		PixelDataJSON += "\"min_sample_count\":" + std::to_string(minSampleCount) + ",";
		PixelDataJSON += "\"max_sample_count\":" + std::to_string(maxSampleCount) + ",";

		/* Sample count image */
		std::vector<uchar> sampleCountImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				sampleCountImage[x + width * y] = ((pixInfo[x + width * y]->samples.size() - minSampleCount) / (float)(maxSampleCount - minSampleCount)) * 255;
			}
		}
		std::cout << "sampleCount.png " << SavePNG("sampleCount.png", sampleCountImage.data(), 1) << std::endl;

		//depth buffer
		idx = 0;
		PixelDataJSON += "\"Depth_buffer\":[";
		float minDepth = BIGFLOAT;
		float maxDepth = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (zbuffer[x + width * y] >= BIGFLOAT) continue;

				minDepth = min(minDepth, zbuffer[x + width * y]);
				maxDepth = max(maxDepth, zbuffer[x + width * y]);
			}
		}


		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(min(maxDepth, zbuffer[idx]));
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max depth
		PixelDataJSON += "\"min_depth\":" + std::to_string(minDepth) + ",";
		PixelDataJSON += "\"max_depth\":" + std::to_string(maxDepth) + ",";

		/* depth image */
		std::vector<uchar> depthImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				depthImage[x + width * y] = ((min(maxDepth, zbuffer[x + y * width]) - minDepth) / (maxDepth - minDepth)) * 255;
			}
		}
		std::cout << "depth.png " << SavePNG("depth.png", depthImage.data(), 1) << std::endl;


		//Ray Box intersections
		idx = 0;
		PixelDataJSON += "\"Box_intersections\":[";
		int minNumRayBoxIntersections = INT32_MAX;
		int maxNumRayBoxIntersections = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int numSamples = pixInfo[idx]->samples.size();
				int numRayBoxIntersections = 0;
				for (int i = 0; i < numSamples; ++i) {
					numRayBoxIntersections += pixInfo[idx]->samples[i]->numRayBoxIntersections;
				}


				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(numRayBoxIntersections);
				minNumRayBoxIntersections = min(minNumRayBoxIntersections, numRayBoxIntersections);
				maxNumRayBoxIntersections = max(maxNumRayBoxIntersections, numRayBoxIntersections);
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max numRayBoxIntersections
		PixelDataJSON += "\"min_box_intersections\":" + std::to_string(minNumRayBoxIntersections) + ",";
		PixelDataJSON += "\"max_box_intersections\":" + std::to_string(maxNumRayBoxIntersections) + ",";

		/* Ray Box Intersections image */
		std::vector<uchar> boxIntersectionsImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int numSamples = pixInfo[x + width * y]->samples.size();
				int numRayBoxIntersections = 0;
				for (int i = 0; i < numSamples; ++i) {
					numRayBoxIntersections += pixInfo[x + width * y]->samples[i]->numRayBoxIntersections;
				}

				boxIntersectionsImage[x + width * y] = ((numRayBoxIntersections - minNumRayBoxIntersections) / (float)(maxNumRayBoxIntersections - minNumRayBoxIntersections)) * 255;
			}
		}
		std::cout << "boxIntersections.png " << SavePNG("boxIntersections.png", boxIntersectionsImage.data(), 1) << std::endl;

		//Ray primitive intersections
		idx = 0;
		PixelDataJSON += "\"Primitive_intersections\":[";
		int minNumRayPrimitiveIntersections = INT32_MAX;
		int maxNumRayPrimitiveIntersections = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int numSamples = pixInfo[idx]->samples.size();
				int numRayPrimitiveIntersections = 0;
				for (int i = 0; i < numSamples; ++i) {
					numRayPrimitiveIntersections += pixInfo[idx]->samples[i]->numRayPrimitiveIntersections;
				}


				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(numRayPrimitiveIntersections);
				minNumRayPrimitiveIntersections = min(minNumRayPrimitiveIntersections, numRayPrimitiveIntersections);
				maxNumRayPrimitiveIntersections = max(maxNumRayPrimitiveIntersections, numRayPrimitiveIntersections);
				idx++;
			}
		}
		PixelDataJSON += "],";

		//min/max numRayBoxIntersections
		PixelDataJSON += "\"min_primitive_intersections\":" + std::to_string(minNumRayPrimitiveIntersections) + ",";
		PixelDataJSON += "\"max_primitive_intersections\":" + std::to_string(maxNumRayPrimitiveIntersections) + ",";

		/* Ray Obj Intersections image */
		std::vector<uchar> objIntersectionsImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int numSamples = pixInfo[x + width * y]->samples.size();
				int numRayPrimitiveIntersections = 0;
				for (int i = 0; i < numSamples; ++i) {
					numRayPrimitiveIntersections += pixInfo[x + width * y]->samples[i]->numRayPrimitiveIntersections;
				}

				objIntersectionsImage[x + width * y] = ((numRayPrimitiveIntersections - minNumRayPrimitiveIntersections) / (float)(maxNumRayPrimitiveIntersections - minNumRayPrimitiveIntersections)) * 255;
			}
		}
		std::cout << "objIntersections.png " << SavePNG("objIntersections.png", objIntersectionsImage.data(), 1) << std::endl;

		//variance
		idx = 0;
		PixelDataJSON += "\"Variance\":[";
		float minVariance = BIGFLOAT;
		float maxVariance = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (idx != 0)
					PixelDataJSON += ",";
				PixelDataJSON += std::to_string(pixInfo[idx]->sampleVariance);
				minVariance = min(minVariance, pixInfo[idx]->sampleVariance);
				maxVariance = max(maxVariance, pixInfo[idx]->sampleVariance);
				idx++;
			}
		}
		PixelDataJSON += "],";

		PixelDataJSON += "\"min_variance\":" + std::to_string(minVariance) + ",";
		PixelDataJSON += "\"max_variance\":" + std::to_string(maxVariance);

		/* Variance image */
		std::vector<uchar> varianceImage(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				varianceImage[x + width * y] = ((pixInfo[x + width * y]->sampleVariance - minVariance) / (float)(maxVariance - minVariance)) * 255;
			}
		}
		std::cout << "variance.png " << SavePNG("variance.png", varianceImage.data(), 1) << std::endl;

		PixelDataJSON += "}";
		PixelDataJSON += "}";

		std::ofstream out(dirname + std::string("pixeldata.json"));
		out << PixelDataJSON;
		out.close();

		////RayData
		int numPixels = width * height;
		int pixelsPerFile = (width * height) / 32;
		int numFiles = numPixels / pixelsPerFile;

		/* ASSUMING SCREEN SIZE IS DIVISABLE BY 64! */
		for (int fileidx = 0; fileidx < numFiles; fileidx++) {
			RayDataJSON = "{";
			RayDataJSON += "\"RayData\":[";
			
			//for each pixel
			for (int i = fileidx * pixelsPerFile; i < (fileidx * pixelsPerFile) + pixelsPerFile; i++) {
				int numSamples = pixInfo[i]->samples.size();
				if (i != fileidx * pixelsPerFile) RayDataJSON += ",";
				
				RayDataJSON += "{\"c\":" + std::string("["
					+ std::to_string(img[i].r/255.0) + ","
					+ std::to_string(img[i].g/255.0) + ","
					+ std::to_string(img[i].b/255.0) + "],");

				RayDataJSON += "\"ch\":[";

				

				//for each sample
				for (int s = 0; s < numSamples; s++) {
					RayDataJSON += pixInfo[i]->samples[s]->toJSON();
					if (s != numSamples - 1) RayDataJSON += ",";
				}
				RayDataJSON += "]}";
			}

			RayDataJSON += "]}";
			std::ofstream out2(dirname + std::string("raydata" + std::to_string(fileidx) + ".json"));
			out2 << RayDataJSON;
			out2.close();

		}


		return true;
	};

private:
	bool SavePNG(const char *filename, uchar *data, int compCount) const
	{
		LodePNGColorType colortype;
		switch( compCount ) {
			case 1: colortype = LCT_GREY; break;
			case 3: colortype = LCT_RGB;  break;
			default: return false;
		}
		unsigned int error = lodepng::encode(filename,data,width,height,colortype,8);
		return error == 0;
	}
};

//-------------------------------------------------------------------------------

#endif
