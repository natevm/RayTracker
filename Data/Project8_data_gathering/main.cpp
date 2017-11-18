#include <stdio.h>	//printf
#include <ctime>	//clock, clock_t
#include <thread>	//thread

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

//ADAPTIVE ANTIALIASING
#define MIN_SAMPLES 4
#define MAX_SAMPLES 8
#define VARIANCE_THRESHOLD 0.001f

#define MAX_BOUNCES 4

//RAY DIFFERENTIALS
//comment this out if you don't want ray differentials
#define RAY_DIFFERENTIALS
//there's also one in scene.h

//comment this out if you don't want to collect ray data
//#define COLLECT_DATA
//there's also one in scene.h

#include "src/scene.h"
#include "src/objects.h"
#include "include/lodepng/lodepng.h"
#include "src/materials.h"

#include "PixelIterator.h"
#include "src/lights.h"
#include <chrono>
//#include <atomic>	//std::atomic

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
TriObj theMesh;
Plane thePlane;
LightList lights;
MaterialList materials;
ObjFileList objList;
ItemFileList<Texture> textureList;
TexturedColor background;
TexturedColor environment;

PixelIterator pixelIterator;
//PixelIterator& pixIter;

int RenderID;	//a flag to tell threads to stop rendering

#define T_MIN 0.001f	//the minimum intersect distance

void Render();
bool Trace(Ray, HitInfo &, RayInfo *);
//void RenderPixels(PixelIterator &it, int rID);

void getHeatMapColor(float value, float &red, float &green, float &blue){
	const int NUM_COLORS = 4;
	static float color[NUM_COLORS][3] = { { 0,0,1 },{ 0,1,0 },{ 1,1,0 },{ 1,0,0 } };
	// A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.

	int idx1;        // |-- Our desired color will be between these two indexes in "color".
	int idx2;        // |
	float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

	if (value <= 0) { idx1 = idx2 = 0; }    // accounts for an input <=0
	else if (value >= 1) { idx1 = idx2 = NUM_COLORS - 1; }    // accounts for an input >=0
	else
	{
		value = value * (NUM_COLORS - 1);        // Will multiply value by 3.
		idx1 = floor(value);                  // Our desired color will be after this index.
		idx2 = idx1 + 1;                        // ... and before this index (inclusive).
		fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
	}

	red = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
	green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
	blue = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
}

void BeginRender() {
	RenderID = 0;
	//numBoxTests = 0;
	//numIntersectionTests = 0;
	//maxNum = 6696;	//Elaborate scene
	//maxNum = 2578;	//Elaborate scene without Happy Bhudda
	//maxNum = 4277;	//Teapot scene
	//maxNum = 0;
	//initialize the PixelIterator
	//pixIter.setSize(camera.imgWidth, camera.imgHeight);
	//start a timer
	const clock_t begin_time = clock();
	Render();
	//stop a timer
	const clock_t end_time = clock();
	float duration = float(end_time - begin_time) / (float)CLOCKS_PER_SEC;
	printf("Render time: %f\n", duration);

	//printf("Max: %u\n", maxNum);
	//printf("Box tests: %u\n", numBoxTests);
	//printf("Intersection tests: %u\n", numIntersectionTests);

	renderImage.ComputeZBufferImage();

	//compute number of samples image
	renderImage.ComputeSampleCountImage();

	//save the image to disk
	renderImage.SaveImage("render.png");
	//save the z-buffer to disk
	renderImage.SaveZImage("render_z.png");
	//save number of samples image to disk
	renderImage.SaveSampleCountImage("render_samples.png");


	//TODO: save ray info to a JSON file
	//TODO: save numIntersectionTests maps to a JSON file. (Or a PNG?)
	renderImage.SavePixInfo("");

}

void StopRender() {
	RenderID++;	//tell the threads to stop rendering
}

bool TriObj::TraceBVHNode(const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID, RayInfo *_rInfo) const {
	/*
	//RECURSIVE SOLUTION
	//get the node in question
	//get its bounding box
	Box curBBox = Box(this->bvh.GetNodeBounds(nodeID));
	float tval;
	//if the ray doesn't hit its bounding box
	if (!curBBox.IntersectRay(ray, hInfo.z, tval)) {
		//return false
		return false;
	}

	bool hitPrimitive = false;
	//if the node is a leaf node
	if (this->bvh.IsLeafNode(nodeID)) {
		//trace every primitive in the node
		int numElements = this->bvh.GetNodeElementCount(nodeID);
		const unsigned int* elements = this->bvh.GetNodeElements(nodeID);
		for (int i = 0; i < numElements; i++) {
			hitPrimitive |= this->IntersectTriangle(ray, hInfo, hitSide, elements[i]);
		}
	}
	//else
	else {
		//get the nodeIDs for each child
		//call TraceBVHNode on each child
		unsigned int ch1, ch2;
		this->bvh.GetChildNodes(nodeID, ch1, ch2);
		
		//TODO: check which child is closer and trace that one first

		hitPrimitive |= this->TraceBVHNode(ray, hInfo, hitSide, ch1);
		hitPrimitive |= this->TraceBVHNode(ray, hInfo, hitSide, ch2);
	}

	return hitPrimitive;*/


	//RayInfo rinfo = ray.rInfo;


	bool hitPrimitive = false;
	
	//create an array of size 50 (or the maximum depth of any bvh tree). This will be our stack.
	unsigned int stack[50];	//50 array elements will handle A LOT of triangles
	unsigned int next = 0;	//the next available stack position
	float tval;	//the distance from the ray to the box we just intersected (in the future)
	//do ray-box intersection on the root node's bounding box
	Box curBBox = Box(this->bvh.GetNodeBounds(nodeID));


	//increment the ray-box intersection test count for this ray
	_rInfo->numRayBoxIntersections++;
	//if the ray doesn't hit the root node
	if (!curBBox.IntersectRay(ray, hInfo.z, tval)) {
		//return false
		return false;
	}
	//add the root node to the stack
	stack[next++] = nodeID;

	//while the stack is not empty
	while (next > 0) {
		//get an item off the top of the stack
		unsigned int curNode = stack[--next];
		curBBox = this->bvh.GetNodeBounds(curNode);
		//do ray-box intersection on the node's bounding box

		//increment the ray-box intersection test count for this ray
		_rInfo->numRayBoxIntersections++;
		//if the t-value for the bbox is further than the current closest hit
		if (!curBBox.IntersectRay(ray, hInfo.z, tval)) {
			//continue
			continue;
		}

		//if it's a leaf node
		if (this->bvh.IsLeafNode(curNode)) {
			//trace every primitive in the node
			int numElements = this->bvh.GetNodeElementCount(curNode);
			const unsigned int* elements = this->bvh.GetNodeElements(curNode);
			for (int i = 0; i < numElements; i++) {
				//increment the ray-primitive intersection test count for this ray
				_rInfo->numRayPrimitiveIntersections++;

				hitPrimitive |= this->IntersectTriangle(ray, hInfo, hitSide, elements[i]);
			}
		}
		//else
		else {
			//get the item's child nodes
			unsigned int ch1, ch2;
			this->bvh.GetChildNodes(curNode, ch1, ch2);
			//t values for each child node
			float t1, t2;
			bool hit1 = false, hit2 = false;
			//do ray-box intersection on each child
			curBBox = this->bvh.GetNodeBounds(ch1);

			//increment the ray-box intersection test count for this ray
			_rInfo->numRayBoxIntersections++;
			//if the ray hits the first one
			if (curBBox.IntersectRay(ray, hInfo.z, t1)) {
				//push it on the stack
				stack[next++] = ch1;
				hit1 = true;
			}
			curBBox = this->bvh.GetNodeBounds(ch2);

			//increment the ray-box intersection test count for this ray
			_rInfo->numRayBoxIntersections++;
			//if the ray hits the second one
			if (curBBox.IntersectRay(ray, hInfo.z, t2)) {
				//push it on the stack
				stack[next++] = ch2;
				hit2 = true;
			}
			//if both children were hit and if the first hit was closer than the second hit
			if (hit1 && hit2 && t1 < t2) {
				//swap the top and second elements in the stack
				unsigned int tmp = stack[next - 1];
				stack[next - 1] = stack[next - 2];
				stack[next - 2] = tmp;
			}
		}
	}
	//return whether or not we hit a primitive
	return hitPrimitive;
}

bool Sphere::IntersectRay(const Ray &ray, HitInfo &hInfo, RayInfo *_rInfo, int hitSide) const {
	//numIntersectionTests++;	//for profiling BVH efficiency
	
	//The ray is in model space
	//Compute the coefficients for the quadratic equation
	float a = ray.dir % ray.dir;	//dot product, right?
	float b = 2.0f * (float)(ray.dir % ray.p);
	float c = (ray.p % ray.p) - 1.0f;	//all spheres are unit spheres in object space, so radius^2 == 1

										//discriminant
	float disc = b * b - 4.0f * a * c;

	//if there are no real roots
	if (disc < 0.0f) {
		//missed the sphere
		return false;
	}

	//compute q
	float discSqrt = std::sqrt(disc);
	float q;
	if (b < 0.0f) {
		q = -0.5f * (b - discSqrt);
	}
	else {
		q = -0.5f * (b + discSqrt);
	}

	//compute intersection distances
	float t0 = q / a;
	float t1;
	if (q != 0.0f) {
		t1 = c / q;
	}
	else {
		t1 = BIGFLOAT;
	}

	//get the closer point
	if (t1 < t0) {
		float temp = t0;
		t0 = t1;
		t1 = temp;
	}
	//make sure it's in front of the ray
	float dist;
	bool hitFront = false;
	//if (t1 < 0.0f) {
	if (t1 < T_MIN) {
		//sphere is entirely behind the ray
		return false;
	}
	//if (t0 >= 0.0f) {
	if (t0 >= T_MIN) {
		//t0 is in front of the ray
		dist = t0;
		hitFront = true;
	}
	//else if (t1 >= 0.0f) {
	else if (t1 >= T_MIN) {
		//ray originated inside the sphere
		dist = t1;
		hitFront = false;
	}
	else {
		//Not sure we need this
		return false;
	}

	//if it's closer than the hit in _hit
	if (dist < hInfo.z) {
		//set the new closest hit distance
		hInfo.z = dist;
		//set the hit side
		hInfo.front = hitFront;
		//fill in the new values of hInfo
		hInfo.p = ray.p + ray.dir * dist;	//the point of intersection in model space
		hInfo.N = hInfo.p;	//the surface normal at the point of intersection in model space
		//since all spheres are unit spheres in model space, the surface normal is equal to the point of intersection

		//calculate the texture space coordinate: get spherical coordinates from the hit point
		//theta = inverse tan(y/x)
		float theta = atan(hInfo.p.y / hInfo.p.x);
		//phi = inverse tan( sqrt(x^2 + y^2)/z )
		float phi = atan2(sqrt(hInfo.p.x * hInfo.p.x + hInfo.p.y * hInfo.p.y), hInfo.p.z);

		hInfo.uvw = Point3(theta / (M_PI * 2.f), phi / M_PI, 1.f);

#ifdef RAY_DIFFERENTIALS
		//Since spheres are rarely (never?) used in production renderers,
		//we'll take a shortcut here and just make up some numbers
		hInfo.duvw[0] = Point3(0.001f, 0.001f, 0);
		hInfo.duvw[1] = Point3(0.001f, 0.001f, 0);
		//this is all a dirty hack
		hInfo.p_u = hInfo.p;
		hInfo.p_v = hInfo.p;
#endif

		return true;
	}
	//it wasn't the closest hit
	return false;
}

bool Plane::IntersectRay(const Ray &ray, HitInfo &hInfo, RayInfo *_rInfo, int hitSide) const {
	//numIntersectionTests++;	//for profiling BVH efficiency
	
	//this is a unit plane in model space where the origin is (0,0,0) and the normal is (0,0,1)
	Point3 N = Point3(0, 0, 1);
	Point3 S = Point3(0, 0, 0);
	//get surface normal DOT ray direction, NdotD
	float NdotD = N % ray.dir;
	//if |NdotD| < 0.00001
	const float epsilon = 0.00001f;	//0.0001 is too big
	if (abs(NdotD) < epsilon) {
		//the ray is parallel, so return false
		return false;
	}

	//get the t-distance from the ray origin to the plane:
	//t = -(P_N / d_N) = -( (ray origin - plane origin) DOT N ) / ( ray direction DOT N )
	float P_N = (ray.p - S) % N;
	float d_N = NdotD;
	float t = -P_N / d_N;

	//if the surface is too close or behind the ray or isn't the closest hit
	if (t < T_MIN || t > hInfo.z) {
		return false;
	}

	//get the hit point, q
	Point3 q = ray.p + ray.dir * t;
	//see if it's within the bounds of the plane
	if (q.x > 1.0 || q.x < -1.0 || q.y > 1.0 || q.y < -1.0) {
		return false;
	}

	//if t is less than hitInfo's t value
	if (t < hInfo.z) {
		//record a hit
		hInfo.z = t;
		hInfo.front = NdotD < 0.0f;
		//hInfo.front = true;
		hInfo.p = q;
		hInfo.N = N;
		//hInfo.N = NdotD < 0.0 ? N : -N;

		//get texture coordinates by transforming the hit point from [-1,1] to [0,1]
		hInfo.uvw = Point3((q.x + 1.f) / 2.f, (q.y + 1.f) / 2.f, 1.f);

#ifdef RAY_DIFFERENTIALS
		//ray differentials

		//calculate NdotD_u
		float NdotD_u = N % ray.dir_u;
		//if dir_u is parallel to the plane
		if (NdotD_u < epsilon) {
			//make it not so parallel
			//TODO: how?
		}
		//calculate t_u using P_N_u and d_N_u
		float P_N_u = (ray.p_u - S) % N;
		float d_N_u = NdotD_u;
		float t_u = -P_N_u / d_N_u;
		//calculate q_u using t_u
		Point3 q_u = ray.p_u + ray.dir_u * t_u;
		//record this point in hInfo
		hInfo.p_u = q_u;
		//get the texture coordinates for q_u
		Point3 uvw_u = Point3((q_u.x + 1.f) / 2.f, (q_u.y + 1.f) / 2.f, 1.f);
		//calculate the ray differential for u
		hInfo.duvw[0] = (uvw_u - hInfo.uvw) * 2.f;

		//calculate NdotD_v
		float NdotD_v = N % ray.dir_v;
		//if dir_v is parallel to the plane
		if (NdotD_v < epsilon) {
			//make it not so parallel
			//TODO: how?
		}
		//calculate t_v using P_N_v and d_N_v
		float P_N_v = (ray.p_v - S) % N;
		float d_N_v = NdotD_v;
		float t_v = -P_N_v / d_N_v;
		//calculate q_v using t_v
		Point3 q_v = ray.p_v + ray.dir_v * t_v;
		//record this point in hInfo
		hInfo.p_v = q_v;
		//get the texture coordinates for q_v
		Point3 uvw_v = Point3((q_v.x + 1.f) / 2.f, (q_v.y + 1.f) / 2.f, 1.f);
		//calculate the ray differential for v
		hInfo.duvw[1] = (uvw_v - hInfo.uvw) * 2.f;
#endif

		return true;
	}
	return false;
}

bool TriObj::IntersectTriangle(const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const {
	//numIntersectionTests++;	//for profiling BVH efficiency

	//get the triangle in question
	TriFace curTri = this->F(faceID);
	//get triangle vertices, A, B, and C
	Point3 A = this->V(curTri.v[0]);
	Point3 B = this->V(curTri.v[1]);
	Point3 C = this->V(curTri.v[2]);

	//get the triangle face normal, N
	Point3 N = ((B - A) ^ (C - A)).GetNormalized();
	//check for nan Normal values. This is probably a degenerate triangle
	if (isnan(N.x)) {
		return false;
	}/**/
	
	float NdotD = ray.dir % N;
	//if the ray direction is parallel to the plane of the triangle
	float epsilon = 0.00001f;	//0.0001 is too big for planes, so probably for triangles too?
	if (abs(NdotD) < epsilon) {
		return false;
	}
	//get the t-distance from the ray origin to the plane:
	//t = -(P_N / d_N) = -( (ray origin - plane origin) DOT N ) / ( ray direction DOT N )
	float P_N = (ray.p - A) % N;
	float d_N = NdotD;
	float t = -P_N / d_N;
	//do a ray-plane intersection test to get the t value

	//if t isn't the closest hit
	if (t < T_MIN || t > hInfo.z) {
		//we don't want it
		return false;
	}

	//calculate the hit point, q
	Point3 q = ray.p + ray.dir * t;
	
	//get |N|, absNorm
	Point3 absNorm = Point3(N);
	absNorm.Abs();
	//get absNorm's max value index using Point3::MaxID(), Nmax
	int Nmax = absNorm.MaxID();
	//also see Point3::[]

	//create 2D projected vertices, pA, pB, and pC
	Point2 pA, pB, pC;
	//create 2D projected hit point, pq
	Point2 pq;
	switch (Nmax) {
	//if Nmax is 0
	case 0:
		//project A, B, C, and q to the YZ plane
		pA = Point2(A.y, A.z);
		pB = Point2(B.y, B.z);
		pC = Point2(C.y, C.z);
		pq = Point2(q.y, q.z);
		break;
	//else if Nmax is 1
	case 1:
		//project A, B, C, and q to the XZ plane
		pA = Point2(A.x, A.z);
		pB = Point2(B.x, B.z);
		pC = Point2(C.x, C.z);
		pq = Point2(q.x, q.z);
		break;
	//else
	case 2:
		//project A, B, C, and q to the XY plane
		pA = Point2(A.x, A.y);
		pB = Point2(B.x, B.y);
		pC = Point2(C.x, C.y);
		pq = Point2(q.x, q.y);
		break;
	}

	//calculate the triangle's area, triArea
	float triArea = ((pB - pA) ^ (pC - pA)) / 2.0;

	//CALCULATE ALPHA
	//calculate the area of triangle (pB, pq, pC), alphaArea
	float alphaArea = ((pB - pq) ^ (pC - pq)) / 2.0;	
	//alpha = alphaArea / triArea
	float alpha = alphaArea / triArea;
	//if alpha is outside of [0, 1]
	if (alpha < 0.0 || alpha > 1.0) {
		return false;
	}

	//CALCULATE BETA
	//calculate the area of triangle (pA, pq, pC), betaArea
	float betaArea = ((pC - pq) ^ (pA - pq)) / 2.0;
	//beta = betaArea / triArea
	float beta = betaArea / triArea;
	//if beta is outside of [0, 1]
	if (beta < 0.0 || beta > 1.0) {
		return false;
	}

	//Calculate Gamma = 1 - (alpha + beta)
	float gamma = 1.0 - (alpha + beta);
	//if gamma < 0
	if (gamma < 0.0) {
		//return false
		return false;
	}

	//if you've gotten to this point, then it's the closest hit!
	//store hit info
	//t value
	hInfo.z = t;
	//front hit
	hInfo.front = NdotD < 0.0;
	//hit point
	hInfo.p = q;
	//surface normal (see Point3f::GetNormal() and normalize it)
	Point3 bc = Point3(alpha, beta, gamma);
	hInfo.N = this->GetNormal(faceID, bc).GetNormalized();
	//hInfo.N = N;	//uncomment this if you want flat shading

	//get texture coordinates
	hInfo.uvw = this->GetTexCoord(faceID, bc);

#ifdef RAY_DIFFERENTIALS
	//ray differentials
	//assuming this is an infinite plane, calculate the barycentric coordinates for the U and V rays
	//calculate NdotD_u
	float NdotD_u = N % ray.dir_u;
	float P_N_u, d_N_u, t_u;
	//if dir_u is parallel to the triangle
	if (abs(NdotD_u) < epsilon) {
		//set t_u to something big, like t * 10 or something
		//just pick an arbitrary large t-value
		t_u = 10000.f;
	}
	else {
		//calculate q_u using P_N_u, d_N_u, to get t_u
		P_N_u = (ray.p_u - A) % N;
		d_N_u = NdotD_u;
		t_u = -P_N_u / d_N_u;
	}
	Point3 q_u = ray.p_u + ray.dir_u * t_u;
	//record this point in hInfo
	hInfo.p_u = q_u;

	//calculate NdotD_v
	float NdotD_v = N % ray.dir_v;
	float P_N_v, d_N_v, t_v;
	//if dir_v is parallel to the triangle
	if (abs(NdotD_v) < epsilon) {
		//set t_v to something big, like t * 10 or something
		//just pick an arbitrary large t-value
		t_v = 10000.f;
	}
	else {
		//calculate q_v using P_N_v, d_N_v, to get t_v
		P_N_v = (ray.p_v - A) % N;
		d_N_v = NdotD_v;
		t_v = -P_N_v / d_N_v;
	}
	Point3 q_v = ray.p_v + ray.dir_v * t_v;
	//record this point in hInfo
	hInfo.p_v = q_v;

	//use Nmax to project q_u and q_v to Pq_u and Pq_v
	Point2 pq_u, pq_v;
	switch (Nmax) {
	case 0:
		//Project q_u and q_v to the YZ plane
		pq_u = Point2(q_u.y, q_u.z);
		pq_v = Point2(q_v.y, q_v.z);
		break;
	case 1:
		//Project q_u and q_v to the XZ plane
		pq_u = Point2(q_u.x, q_u.z);
		pq_v = Point2(q_v.x, q_v.z);
		break;
	case 2:
		//Project q_u and q_v to the XY plane
		pq_u = Point2(q_u.x, q_u.y);
		pq_v = Point2(q_v.x, q_v.y);
		break;
	}
	float invTriArea = 1.0f / triArea;	//TODO: is this actually faster?
	//calculate Alpha_u, Beta_u, and Gamma_u
	float alpha_u = (((pB - pq_u) ^ (pC - pq_u)) / 2.0f) * invTriArea;
	float beta_u = (((pC - pq_u) ^ (pA - pq_u)) / 2.0f)  * invTriArea;
	float gamma_u = (((pA - pq_u) ^ (pB - pq_u)) / 2.0f) * invTriArea;
	//calculate Alpha_v, Beta_v, and Gamma_v
	float alpha_v = (((pB - pq_v) ^ (pC - pq_v)) / 2.0f) * invTriArea;
	float beta_v = (((pC - pq_v) ^ (pA - pq_v)) / 2.0f)  * invTriArea;
	float gamma_v = (((pA - pq_v) ^ (pB - pq_v)) / 2.0f) * invTriArea;
	//get the texture coordinates for these barycentric coordinates, bc_u
	Point3 bc_u = Point3(alpha_u, beta_u, gamma_u);
	Point3 uvw_u = this->GetTexCoord(faceID, bc_u);
	//get the texture coordinates for these barycentric coordinates, bc_v
	Point3 bc_v = Point3(alpha_v, beta_v, gamma_v);
	Point3 uvw_v = this->GetTexCoord(faceID, bc_v);
	//calculate du/ds: (uvw_u - uvw) * 2.0
	hInfo.duvw[0] = (uvw_u - hInfo.uvw) * 2.0f;
	//calculate dv/dt: (uvw_v - uvw) * 2.0
	hInfo.duvw[1] = (uvw_v - hInfo.uvw) * 2.0f;
#endif

	return true;
}

bool TriObj::IntersectRay(const Ray &ray, HitInfo &hInfo, RayInfo *_rInfo, int hitSide) const {
	
	/*bool hitObj = false;
	
	//WITHOUT BVH
	unsigned int numTris = this->NF();
	//for each triangle
	for (unsigned int i = 0; i < numTris; i++) {
		//if the ray hit the triangle
			//set hitObj to true
		hitObj |= this->IntersectTriangle(ray, hInfo, hitSide, i);
	}
	return hitObj;/**/
	
	
	//WITH BVH
	//Call TraceBVHNode on our BVH root node
	return this->TraceBVHNode(ray, hInfo, hitSide, this->bvh.GetRootNodeID(), _rInfo);
	
}

bool Box::IntersectRay(const Ray &ray, float t_max, float &tval) const {
	//numBoxTests++;	//for profiling BVH efficiency
	Point3 absDir = Point3(ray.dir);
	absDir.Abs();
	float tmin = -BIGFLOAT, tmax = BIGFLOAT;

	//X, Y, and Z planes in a loop
	for (int i = 0; i < 3; i++) {
		//if the ray direction[i] is not 0
		if (absDir[i] > 0.00001f) {	//0.0001 is too big
			//get t0 and t1
			float t0 = (pmin[i] - ray.p[i]) / ray.dir[i];
			float t1 = (pmax[i] - ray.p[i]) / ray.dir[i];
			if (t0 > t1) {
				float tmp = t0;
				t0 = t1;
				t1 = tmp;
			}
			tmin = max(tmin, t0);
			tmax = min(tmax, t1);
		}
		//else
		else {
			//The ray is parallel to the planes
			//if the ray origin i-component falls outside the ith plane's bounds, return false
			if (ray.p[i] > pmax[i] || ray.p[i] < pmin[i]) {
				return false;
			}
		}
		//do an early rejection test
		if (tmin > tmax) {
			return false;
		}
	}

	//if you've already found a closer hit
	if (tmin > t_max) {
		return false;
	}

	//if the whole box is behind you
	if (tmax < T_MIN) {
		return false;
	}

	//set the tval so the caller knows how far the intersection was
	tval = tmin;

	//if the ray started inside the box
	/*if (tmin < T_MIN) {
		return true;
	}*/	//this test is redundant since you'll return true whether or not you enter this code block

	return true;
}

cy::Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, RayInfo *_rInfo) const {
	//get the hit point and surface normal (they should already be in world space)
	Point3 w_p = hInfo.p;
	Point3 w_N = hInfo.N;
	//Get the view vector, V
	//Point3 V = (ray.p - w_p).GetNormalized();	//this is too complicated!
	//isn't this the same thing?: I think it is
	Point3 V = -ray.dir;
	float NdotV = w_N % V;

#ifdef RAY_DIFFERENTIALS
	Point3 V_u = -ray.dir_u;
	Point3 V_v = -ray.dir_v;

	float NdotV_u = w_N % V_u;
	float NdotV_v = w_N % V_v;
#endif

	//create a Color return value, shCol
	cy::Color shCol = cy::Color();
	shCol.SetBlack();
	//get this material's diffuse and specular colors, kd and ks
	cy::Color kd;
#ifdef RAY_DIFFERENTIALS
	//use ray differentials to get a texture sample
	kd = diffuse.Sample(hInfo.uvw, hInfo.duvw, true);	//true for elliptical sampling
#else
	kd = diffuse.Sample(hInfo.uvw);
#endif

	cy::Color ks = this->specular.GetColor();	//TODO: sample instead of GetColor

	int lightCount = lights.size();
	//for each light
	for (int i = 0; i < lightCount; i++) {
		Light *curLight = lights[i];
		//Get the intensity of the light, Li
		cy::Color Li = curLight->Illuminate(w_p, w_N);
		//if it's an ambient light
		if (curLight->IsAmbient()) {
			//add the light's ambient term to shCol
			shCol += kd * Li;
			//shCol += Li;
		}
		else {
			//Get the light vector, L
			Point3 L = -(curLight->Direction(w_p));

			//FOR REGULAR PHONG:
			//Get the reflection vector, R
			//Point3 R = (-L) - 2 * ((-L) % w_N) * w_N;
			//END PHONG STUFF

			//Calculate the half vector, H
			Point3 H = (L + V).GetNormalized();
			//Get the shininess of the light, alpha
			float alpha = this->glossiness;
			//Calculate the Blinn-Phong BRDF:
			//Li * kd * (N dot L) + Li * ks * (N dot H)^alpha
			//add it to shCol
			//shCol += Li * kd * max(0.0, (w_N % L)) + Li * ks * pow(max(0.0, (w_N % H)), alpha);	//use N dot L only on the diffuse term
			shCol += (Li * kd + Li * ks * pow(max(0.0, (w_N % H)), alpha)) * max(0.0, (w_N % L));	//use N dot L on everything

			//REGULAR PHONG:
			//shCol += Li * kd * max(0.0, (w_N % L)) + Li * ks * pow(max(0.0, (V % R)), alpha);

			//create a rayInfo for the shadow ray
			RayInfo *shadowInfo = new RayInfo();
			shadowInfo->Init();
			//we can get origin, direction, colorContribution, and rayType
			shadowInfo->Origin = w_p;
			shadowInfo->Direction = L;
			shadowInfo->colorContribution = Li;
			shadowInfo->rayType = 1;	//for shadow ray
			//we can't get hitPoint, numRayBoxIntersectionTests, numRayPrimitiveIntersectionTests
			
			_rInfo->shadowRays.push_back(shadowInfo);

		}
	}

	//if we're at the end of our bounces
	if (bounceCount <= 0) {
		//return what we have
		return shCol;
	}

	//get the reflective component, kr
	cy::Color kr = this->reflection.GetColor();	//TODO: sample instead of GetColor

	//get the reflection direction ready
	Point3 R_v = Point3();
#ifdef RAY_DIFFERENTIALS
	//get the differential reflection directions ready
	Point3 R_v_u = Point3();
	Point3 R_v_v = Point3();
#endif
	bool REFLECTION_DIRECTION_DONE = false;

	//Refractive component
	cy::Color kt = this->refraction.GetColor();
	//if there's a refraction component
	if (kt.r + kt.g + kt.b > 0) {
		//CALCULATE THE REFRACTION DIRECTION
		
		//if N dot V is positive, we're on the outside going in
		if (NdotV > 0) {
			//n1 is 1, n2 is the material's IOR
			float n1 = 1.f;
			float n2 = this->ior;
			//calculate sin(theta1) = sqrt( 1 - (N dot V)^2 )
			//make sure the stuff inside the sqrt is positive (if it's negative, it should probably be 0)
			float sinTheta1 = sqrt(max(0.f, 1.f - (NdotV * NdotV)));
			//calculate sin(theta2) = (n1 / n2) * sinTheta1
			float sinTheta2 = (n1 / n2) * sinTheta1;
			//calculate cos(theta1)
			float cosTheta1 = NdotV;
			//calculate cos(theta2) = sqrt( 1 - (sinTheta2)^2 )	//make sure it's positive
			float cosTheta2 = sqrt(max(0.f, 1.f - (sinTheta2 * sinTheta2)));
			//calculate S: N x (N x V) (normalized)
			Point3 S = (w_N ^ (w_N ^ V)).GetNormalized();
			//float len = S.Length();
			//calculate T: -N * cosTheta2 + S * sinTheta2
			Point3 T = -w_N * cosTheta2 + S * sinTheta2;
			//len = T.Length();
			//Point3 T = (-w_N * cosTheta2 + S * sinTheta2).GetNormalized();

#ifdef RAY_DIFFERENTIALS
			/*float sinTheta1_u = sqrt(max(0.f, 1.f - (NdotV_u * NdotV_u)));
			float sinTheta2_u = (n1 / n2) * sinTheta1_u;
			float cosTheta1_u = NdotV_u;
			float cosTheta2_u = sqrt(max(0.f, 1.f - (sinTheta2_u * sinTheta2_u)));*/
			Point3 S_u = (w_N ^ (w_N ^ V_u)).GetNormalized();
			//Point3 T_u = -w_N * cosTheta2_u + S_u * sinTheta2_u;
			Point3 T_u = -w_N * cosTheta2 + S_u * sinTheta2;

			/*float sinTheta1_v = sqrt(max(0.f, 1.f - (NdotV_v * NdotV_v)));
			float sinTheta2_v = (n1 / n2) * sinTheta1_v;
			float cosTheta1_v = NdotV_v;
			float cosTheta2_v = sqrt(max(0.f, 1.f - (sinTheta2_v * sinTheta2_v)));*/
			Point3 S_v = (w_N ^ (w_N ^ V_v)).GetNormalized();
			//Point3 T_v = -w_N * cosTheta2_v + S_v * sinTheta2_v;
			Point3 T_v = -w_N * cosTheta2 + S_v * sinTheta2;
#endif
			//calculate R0 = ( (n1 - n2) / (n1 + n2) )^2
			float R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
			//calculate the Fresnel reflection term, Fr_r = R0 + (1 - R0)(1 - cosTheta1)^5
			float Fr_r = R0 + (1.f - R0) * pow(1.f - cosTheta1, 5);
			//calculate the Fresnel refraction term, Fr_t = 1 - Fr_r
			float Fr_t = 1.f - Fr_r;

			//if Fr_r > 0
			if (Fr_r > 0) {
				//calculate the reflection direction around N
				R_v = 2.0f * w_N * (NdotV) - V;

#ifdef RAY_DIFFERENTIALS
				//calculate reflection directions for differential rays, R_v_u and R_v_v
				R_v_u = 2.0f * w_N * (NdotV_u) - V_u;
				R_v_v = 2.0f * w_N * (NdotV_v) - V_v;
#endif
				//set REFLECTION_DIRECTION_DONE to true
				REFLECTION_DIRECTION_DONE = true;
				//? add Fr_r to kr?
				kr += cy::Color(Fr_r);	//TODO: is this right?
			}

			//trace a ray in the direction of T
#ifdef RAY_DIFFERENTIALS
			Ray transRay = Ray(w_p, T, hInfo.p_u, T_u, hInfo.p_v, T_v);
#else
			Ray transRay = Ray(w_p, T);
#endif

			RayInfo *transInfo = new RayInfo();
			transInfo->Origin = transRay.p;
			transInfo->Direction = transRay.dir;
			_rInfo->refractionRay = transInfo;
			transInfo->rayType = 3;	//for refraction

			HitInfo transHit = HitInfo();
			Trace(transRay, transHit, transInfo);
			//if it hit something
			if (transHit.node != NULL) {
				//transRay.rInfo	//TODO: store the object index that it hit
				transInfo->HitPoint = transHit.p;

				//shade that point
				const Material* transMat = transHit.node->GetMaterial();
				cy::Color transCol = transMat->Shade(transRay, transHit, lights, bounceCount - 1, transInfo);
				//multiply it by the refraction color
				//calculate the length of the transmission ray, h
				float transLen = (w_p - transHit.p).Length();
				cy::Color absorption = this->absorption;
				//calculate the absorption factor: e^(-h * absorption coefficient)
				float abs_fac_r = exp(-transLen * absorption.r);
				float abs_fac_g = exp(-transLen * absorption.g);
				float abs_fac_b = exp(-transLen * absorption.b);
				cy::Color abs_fac = cy::Color(abs_fac_r, abs_fac_g, abs_fac_b);
				//multiply the shaded point by the absorption factor
				//add that into shCol
				shCol += transCol * kt * abs_fac;

				transInfo->colorContribution = transCol * kt * abs_fac;

			}
			//otherwise, why not? You should have hit something
			else {
				
			}
		}
		//else if N dot V is negative, we're on the inside going out
		else {
			//n1 is the material's IOR, n2 is 1
			float n1 = this->ior;
			float n2 = 1.f;
			//calculate sinTheta1 = sqrt( 1 - (-N dot V)^2 )	//make sure it's positive
			float sinTheta1 = sqrt(max(0.f, 1.f - (-NdotV * -NdotV)));
			//calculate sinTheta2 = (n1 / n2) * sinTheta1
			float sinTheta2 = (n1 / n2) * sinTheta1;
			//if sinTheta2 > 1
			if (sinTheta2 > 1.f) {
				//we have TOTAL INTERNAL REFLECTION
				//calculate the reflection direction around -N
				Point3 ttl_R_v = 2.0f * -w_N * (-NdotV) - V;
#ifdef RAY_DIFFERENTIALS
				//calculate reflection directions around -N for differential rays
				Point3 ttl_R_v_u = 2.0f * -w_N * (-NdotV_u) - V_u;
				Point3 ttl_R_v_v = 2.0f * -w_N * (-NdotV_v) - V_v;
				//create the ray with differential rays
				Ray ttl_reflRay = Ray(w_p, ttl_R_v, hInfo.p_u, ttl_R_v_u, hInfo.p_v, ttl_R_v_v);
#else
				Ray ttl_reflRay = Ray(w_p, ttl_R_v);
#endif
				RayInfo *ttl_reflInfo = new RayInfo();
				ttl_reflInfo->Init();
				ttl_reflInfo->Origin = ttl_reflRay.p;
				ttl_reflInfo->Direction = ttl_reflRay.dir;
				_rInfo->refractionRay = ttl_reflInfo;
				ttl_reflInfo->rayType = 3;	//for refraction

				HitInfo ttl_reflHit = HitInfo();
				//trace a ray in the total internal reflection direction
				Trace(ttl_reflRay, ttl_reflHit, ttl_reflInfo);
				//if it hit something
				if (ttl_reflHit.node != NULL) {

					//ttl_reflRay.rInfo	//TODO: store the object index that it hit
					ttl_reflInfo->HitPoint = ttl_reflHit.p;

					//shade that point
					const Material* ttlReflMat = ttl_reflHit.node->GetMaterial();
					cy::Color ttl_reflCol = ttlReflMat->Shade(ttl_reflRay, ttl_reflHit, lights, bounceCount - 1, ttl_reflInfo);
					//calculate absorption for this edge of the path
					float ttlReflLen = (w_p - ttl_reflHit.p).Length();
					cy::Color absorption = this->absorption;
					cy::Color abs_fac = cy::Color();
					abs_fac.r = exp(-ttlReflLen * absorption.r);
					abs_fac.g = exp(-ttlReflLen * absorption.g);
					abs_fac.b = exp(-ttlReflLen * absorption.b);
					//DEBUG
					//abs_fac.SetWhite();
					//add that into shCol
					shCol += ttl_reflCol * abs_fac;

					ttl_reflInfo->colorContribution = ttl_reflCol * abs_fac;

				}
			}
			//else
			else {
				//calculate cosTheta2
				float cosTheta2 = sqrt(max(0.f, 1.f - (sinTheta2 * sinTheta2)));
				//calculate S: -N x (-N x V) (normalized)
				Point3 S = (w_N ^ (w_N ^ V)).GetNormalized();
				//calculate T: N * cosTheta2 + S * sinTheta2
				Point3 T = w_N * cosTheta2 + S * sinTheta2;
				//T = -V;	//to simulate the bug where back hits are ignored
#ifdef RAY_DIFFERENTIALS
				//calculate refraction directions for differential rays
				//float sinTheta1_u = sqrt(max(0.f, 1.f - (-NdotV_u * -NdotV_u)));
				//float sinTheta2_u = (n1 / n2) * sinTheta1_u;
				////cosTheta1_u == NdotV_u
				//float cosTheta2_u = sqrt(max(0.f, 1.f - (sinTheta2_u * sinTheta2_u)));
				Point3 S_u = (w_N ^ (w_N ^ V_u)).GetNormalized();
				//Point3 T_u = w_N * cosTheta2_u + S_u * sinTheta2_u;
				Point3 T_u = w_N * cosTheta2 + S_u * sinTheta2;

				//float sinTheta1_v = sqrt(max(0.f, 1.f - (-NdotV_v * -NdotV_v)));
				//float sinTheta2_v = (n1 / n2) * sinTheta1_v;
				////cosTheta1_v == NdotV_v
				//float cosTheta2_v = sqrt(max(0.f, 1.f - (sinTheta2_v * sinTheta2_v)));
				Point3 S_v = (w_N ^ (w_N ^ V_v)).GetNormalized();
				//Point3 T_v = w_N * cosTheta2_v + S_v * sinTheta2_v;
				Point3 T_v = w_N * cosTheta2 + S_v * sinTheta2;
				//create a ray with differential rays
				Ray transRay = Ray(w_p, T, hInfo.p_u, T_u, hInfo.p_v, T_v);
#else
				//trace a ray in the direction of T
				Ray transRay = Ray(w_p, T);
#endif
				RayInfo *transInfo = new RayInfo();
				transInfo->Init();
				transInfo->Origin = transRay.p;
				transInfo->Direction = transRay.dir;
				_rInfo->refractionRay = transInfo;
				transInfo->rayType = 3;	//for refraction

				HitInfo transHit = HitInfo();
				Trace(transRay, transHit, transInfo);
				//if it hit something
				if (transHit.node != NULL) {

					//TODO:	get the object index
					transInfo->HitPoint = transHit.p;

					//shade that point
					const Material* transMat = transHit.node->GetMaterial();
					cy::Color transCol = transMat->Shade(transRay, transHit, lights, bounceCount - 1, transInfo);
					//add that into shCol
					shCol += transCol;

					transInfo->colorContribution = transCol;

				}
				//otherwise, see if there's an environment map
				else {
					shCol += environment.SampleEnvironment(transRay.dir);

					transInfo->colorContribution = environment.SampleEnvironment(transRay.dir);

				}
			}
		}
	}
		
	
	//if the reflection color isn't black
	if (kr.r + kr.g + kr.b > 0.0f) {
		//if REFLECTION_DIRECTION_DONE is false
		if (REFLECTION_DIRECTION_DONE == false) {
			//calculate the reflection direction around N
			R_v = 2.0f * w_N * (NdotV) - V;
#ifdef RAY_DIFFERENTIALS
			//calculate reflection direction for differential rays
			R_v_u = 2.0f * w_N * (NdotV_u) - V_u;
			R_v_v = 2.0f * w_N * (NdotV_v) - V_v;
		}
		//create a ray with differential ray directions
		Ray reflRay = Ray(w_p, R_v, hInfo.p_u, R_v_u, hInfo.p_v, R_v_v);
#else
		}
		Ray reflRay = Ray(w_p, R_v);
#endif
		RayInfo *reflInfo = new RayInfo();
		reflInfo->Init();
		reflInfo->Origin = reflRay.p;
		reflInfo->Direction = reflRay.dir;
		reflInfo->rayType = 2;	//for reflection
		_rInfo->reflectionRay = reflInfo;

		HitInfo reflHit = HitInfo();
		//fire a ray in the direction of reflection
		Trace(reflRay, reflHit, reflInfo);
		//if it hit something
		if (reflHit.node != NULL) {

			//TODO: get object index
			reflInfo->HitPoint = reflHit.p;

			const Material* pixMat = reflHit.node->GetMaterial();
			cy::Color reflCol = pixMat->Shade(reflRay, reflHit, lights, bounceCount - 1, reflInfo);
			shCol += kr * reflCol;

			reflInfo->colorContribution = kr * reflCol;

		}
		//otherwise, see if there's an environment map
		else {
			shCol += kr * environment.SampleEnvironment(reflRay.dir);

			reflInfo->colorContribution = kr * environment.SampleEnvironment(reflRay.dir);

		}
	}

	return shCol;
}

void ShowViewport();
int LoadScene(const char *filename);

int main() {

	//pixIter = pixelIterator;

	//load the scene
	//LoadScene("scene_5.xml");
	//LoadScene("simple_scene_5.xml");
	//LoadScene("scene_5_hypergrid.xml");
	//LoadScene("scene_6_ajax.xml");
	LoadScene("scene_7.xml");
	//LoadScene("scene_5.xml");

	//Show the viewport
	ShowViewport();

	return 0;
}

bool TraceNode(Node *_n, Ray _r, HitInfo &_hit, RayInfo *_rInfo) {

	//if _n has an object
	Object *nodeObj = _n->GetNodeObj();
	//transform the ray into node space
	Ray trans_ray = _n->ToNodeCoords(_r);
	bool hitThisNode = false;
	if (nodeObj != NULL) {
		//do a ray-bounding box intersection test first
		Box nodeBox = nodeObj->GetBoundBox();
		//we need to pass a variable to get the t value of a possible intersection
		float tval;

		_rInfo->numRayBoxIntersections++;

		//if it hit the bounding box

		if (nodeBox.IntersectRay(trans_ray, _hit.z, tval)) {

			//do a ray-node intersection test
			//if there was a new closest intersection, set _hit's node

			_rInfo->numRayPrimitiveIntersections++;
			if (nodeObj->IntersectRay(trans_ray, _hit, _rInfo)) {

				_hit.node = _n;
				hitThisNode = true;
			}
		}
	}
	bool childNodeHit = false;
	//for each child node
	size_t numChildren = _n->GetNumChild();
	for (size_t i = 0; i < numChildren; i++) {
		//call TraceNode on the child
		childNodeHit |= TraceNode(_n->GetChild(i), trans_ray, _hit, _rInfo);
	}
	//if the ray intersected this node
	if (hitThisNode || childNodeHit) {
		//transform the hitInfo from this node's space to its parent space
		_n->FromNodeCoords(_hit);
	}
	//return whether or not an object was hit by this node or one of its descendants
	return hitThisNode | childNodeHit;
}


bool Trace(Ray _r, HitInfo &_hit, RayInfo *_rInfo) {

	//get the root node
	Node *root = &rootNode;
	//trace the ray against that node

	return TraceNode(root, _r, _hit, _rInfo);

}

void RenderPixels(PixelIterator &it, int rID) {
	//Grab pointers to the pixel array and the zbuffer array, img and zbuff
	Color24 *img = renderImage.GetPixels();
	float *zbuff = renderImage.GetZBuffer();
	uchar *sampCount = renderImage.GetSampleCount();
	PixelInfo **pixInfos = renderImage.GetPixelInfoMap();
	
	int imgWidth = camera.imgWidth;
	int imgHeight = camera.imgHeight;
	//distance to the viewing plane, d
	float d = 1.0f;
	//calculate the height of the viewing plane, h
	float h = 2.0f * d * tan((camera.fov * M_PI / 180.0f) / 2.0f);
	//calculate the width of the viewing plane, w
	float w = (h * (float)imgWidth) / (float)imgHeight;
	//calculate camera viewing vectors:
	//y = camera.up
	Point3 y = camera.up;
	//z = -camera.dir
	Point3 z = -camera.dir;
	//x = camera.dir cross camera.up
	Point3 x = (camera.dir ^ camera.up).GetNormalized();
	//calculate top-middle point of the viewing plane, B
	Point3 B = camera.pos + (d * -z) + (h / 2.0f) * y;
	//calculate top-left corner of the viewing plane, A
	Point3 A = B + (w / 2.0f) * -x;
	//calculate horizontal and vertical pixel vectors, u and v
	Point3 u = x * (w / (float)imgWidth);
	Point3 v = -y * (h / (float)imgHeight);
	//int numRendered = 0;
	int i_b_x, i_b_y,	//the x and y indices of the pixel blocks
		nx, ny,			//number of pixels in x and y of the current pixel block
		pixIdx;
	int i_x, i_y;	//the x and y coordinates of the pixel
	while (it.GetPixel(i_b_x, i_b_y, nx, ny)) {
		//for each row of pixels in the pixel block
		for (int dy = 0; dy < ny; dy++) {
			i_y = i_b_y + dy;	//get the y pixel location
			//for each column of pixels in the pixel block
			for (int dx = 0; dx < nx; dx++) {
				i_x = i_b_x + dx;	//get the x pixel location

				pixIdx = i_y * imgWidth + i_x;

				//476 395
				/*if (i_x == 476 && i_y == 395) {
				int one = 2;
				}*/

				//16 16
				/*if (i_x == 16 && i_y == 16) {
					int one = 2;
				}*/
				if (pixIdx == 3) {
					int one = 2;
				}

				//create a pixelInfo for this pixel
				PixelInfo *pInfo = new PixelInfo();

				//start the clock on this pixel
				const clock_t pixStart = clock();
				auto start = std::chrono::high_resolution_clock::now();

				//adaptive antialiasing
				//keep track of the number of samples so far, numSamples
				unsigned int numSamples = 0;
				//keep track of the sample variance for this pixel, pixVariance
				float pixVariance = 0.f;
				cy::Color pixColor = cy::Color();
				pixColor.SetBlack();
				cy::Color subPixColor = cy::Color();
				subPixColor.SetBlack();

				float sampleTotal = 0.f;		//the total accumulated sample value (for finding the mean)
				float sampleSquaredTotal = 0.f;	//the total accumulated square sample value

				//while numSamples < MIN_SAMPLES OR (numSamples < MAX_SAMPLES AND pixVariance > VARIANCE_THRESHOLD)
				while(numSamples < MIN_SAMPLES || (numSamples < MAX_SAMPLES && pixVariance > VARIANCE_THRESHOLD)){
					//use the Halton Sequence to generate a pixel sample location
					float haltonX = Halton(numSamples, 2);
					float haltonY = Halton(numSamples, 3);	//TODO: offset these by a random amount per pixel, make sure to wrap around
					
					//trace a camera ray and get the color, add it to pixColor
					//create a camera ray
					Point3 q = A + ((float)(i_x)+haltonX) * u + ((float)(i_y)+haltonY) * v;

#ifdef RAY_DIFFERENTIALS
					//TODO: HOW DO WE DO THIS WITH MULTISAMPLING???
					//calculate differential ray directions
					/*Point3 q_u = A + ((float)(i_x)+1.0f) * u + ((float)(i_y)+0.5f) * v;
					Point3 q_v = A + ((float)(i_x)+0.5f) * u + ((float)(i_y)+1.0f) * v;*/
					Point3 q_u = A + ((float)(i_x)+haltonX+0.5f) * u + ((float)(i_y)+haltonY) * v;
					Point3 q_v = A + ((float)(i_x)+haltonX) * u + ((float)(i_y)+haltonY+0.5f) * v;

					Ray qRay = Ray(camera.pos, q - camera.pos, camera.pos, q_u - camera.pos, camera.pos, q_v - camera.pos);
#else
					Ray qRay = Ray(camera.pos, q - camera.pos);
#endif

					qRay.Normalize();

					//create a rayInfo to pass through
					RayInfo *sampleRayInfo = new RayInfo();
					sampleRayInfo->Origin = qRay.p;
					sampleRayInfo->Direction = qRay.dir;
					sampleRayInfo->rayType = 0;	//for camera
					
					
					/*qRay.rInfo.Init();
					qRay.rInfo.Origin = qRay.p;
					qRay.rInfo.Direction = qRay.dir;
					qRay.rInfo.rayType = 0;	//for camera*/

					//create a hitInfo object
					HitInfo	qHit = HitInfo();

					//trace the ray
					Trace(qRay, qHit, sampleRayInfo);
					//if it hit something
					if (qHit.node != NULL) {

						//TODO: store object index
						sampleRayInfo->HitPoint = qHit.p;

						//shade the pixel
						const Material* pixMat = qHit.node->GetMaterial();
						subPixColor = pixMat->Shade(qRay, qHit, lights, MAX_BOUNCES, sampleRayInfo);
						//for BVH profiling, map the intersection test counts to colors
						//maxNum = max(maxNum, numBoxTests + numIntersectionTests);


						sampleRayInfo->colorContribution = subPixColor;

					}
					else {
						//set the pixel color to the background image color
						Point3 texSpace = Point3(float(i_b_x + dx + haltonX) / float(imgWidth), float(i_b_y + dy + haltonY) / float(imgHeight), 0.f);

						subPixColor = background.Sample(texSpace);
					}
					zbuff[pixIdx] = qHit.z;

					pixColor += subPixColor;

					sampleRayInfo->colorContribution = subPixColor;


					//increment numSamples
					numSamples++;
					
					//combine the color values into one value, sample
					float subPixSample = subPixColor.Sum() / 3.f;
					sampleTotal += subPixSample;
					sampleSquaredTotal += subPixSample * subPixSample;
					//calculate the sample variance of this pixel, pixVariance
					//get the average of the squared samples
					float A = sampleSquaredTotal / float(numSamples);
					//get the average of the samples
					float C_bar = sampleTotal / float(numSamples);
					//get ((2 * average of the samples) / numSamples) * the sum of the samples, B
					float B = ((2.f * C_bar) / float(numSamples)) * sampleTotal;
					//combine them: A + C_bar^2 - B
					pixVariance = A + (C_bar * C_bar) - B;

					//store this pixel variance on the pixelInfo
					pInfo->sampleVariance = pixVariance;

					//store the number of secondary rays
					pInfo->numSecondaryRays += sampleRayInfo->numSecondaryRays();

					//push the sample info onto the PixelInfo's vector of samples
					pInfo->samples.push_back(sampleRayInfo);
				}
				
				//divide the accumulated pixel color by the number of samples
				pixColor = pixColor / float(numSamples);

				//end the clock on this pixel
				const clock_t pixEnd = clock();
				auto finish = std::chrono::high_resolution_clock::now();
				//std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() << "ns\n";

				//record the time taken for this pixel
				pInfo->renderTime = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();// (pixEnd - pixStart) / (float)CLOCKS_PER_SEC;
				//store it in the map
				pixInfos[pixIdx] = pInfo;

				//store the pixel color
				img[pixIdx].r = (unsigned char)(pixColor.r * 255.f);
				img[pixIdx].g = (unsigned char)(pixColor.g * 255.f);
				img[pixIdx].b = (unsigned char)(pixColor.b * 255.f);

				//record the number of samples in its map
				sampCount[pixIdx] = numSamples;

				//numRendered++;
				//if the user kills the render early
				if (rID != RenderID) {
					return;
				}
			}
		}
	}
	//printf("%d\n", numRendered);
}



void Render() {
	//get the number of available threads
	int nThreads = std::thread::hardware_concurrency();
	if (nThreads == 0) {
		nThreads = 1;
	}
	nThreads = 1;
	printf("Rendering with %d threads\n", nThreads);
	PixelIterator& pixIter = pixelIterator;
	pixIter.Init();
	pixIter.setSize(camera.imgWidth, camera.imgHeight, 16);

	std::vector<std::thread> threads;
	
	
	int& renderid = RenderID;
	//spawn a handful of threads
	for (int i = 0; i < nThreads; i++) {
		threads.push_back(std::thread(RenderPixels, std::ref(pixIter), RenderID));	//we call detach for some reason
		//threads.push_back(th);
		//th.detach();
		//std::thread(threadedFunction, i, std::ref(pixIter)).detach();
	}
	//RenderPixels(pixIter, RenderID);
	//wait for everyone to finish
	for (int i = 0; i < nThreads; i++) {
		threads[i].join();
	}
}


class GenLight;
float GenLight::Shadow(Ray ray, float t_max) {
	//TODO:
		//implement a shadow ray intersect function

	RayInfo *shadowInfo = new RayInfo();	//TODO: how do we get this one's info?
	HitInfo shadowHit = HitInfo();
	shadowHit.z = t_max;
	//push the ray origin a tiny amount in the ray direction
	//ray.p += ray.dir * 0.001f;	//THIS IS BAD
	//fire ray into the scene
	Trace(ray, shadowHit, shadowInfo);
	//if it hits something before t_max
	if (shadowHit.z < t_max) {
		//return 0
		return 0.0f;
	}
	//return 1
	return 1.0f;
}

