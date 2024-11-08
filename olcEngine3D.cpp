#include "olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <windows.h> // For finding the pwd
#include <algorithm>

struct vec3d
{
	float x = 0, y = 0, z = 0, w = 1; // originally only was x,y,z
};

struct triangle
{
	vec3d p[3];

	wchar_t sym;
	short col;
};

struct mesh
{
	std::vector<triangle> tris;
	 
	bool LoadFromObjectFile(std::string sFilename)   // Requires revisiting
	{
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		std::vector<vec3d> verts;

		while (f.eof() == false)
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;   // Sample : "v 1.000000 -1.000000 -1.000000"

			char identifier;

			if (line[0] == 'v')
			{
				// strstream automatically converts data separated by " ".
				vec3d v;
				s >> identifier >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> identifier >> f[0] >> f[1] >> f[2];	
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}	
		}
		return true;
	}
};

// Template Matrix
struct mat4x4
{
	float m[4][4]{ 0 };
};

////////////////////////////////////////////////////////////////////////////

class olcEngine3D : public olcConsoleGameEngine
{

public:
	olcEngine3D()
	{
		m_sAppName = L"3D Demo";
	}

private:
	mesh Cube;
	mat4x4 matProj; // Projecion Matrix

	vec3d vCamera;
	vec3d vLookDir; // Travels along the direction we want the camera to point
	
	float fYaw; // rotation in Y axis



	float fTheta;

	// passing 2 vectors to not mess with the input one
	vec3d  Matrix_MultiplyVector( mat4x4& m, vec3d& in) // usual matrix multiplication
	{
		vec3d v;
		v.x = in.x * m.m[0][0] + in.y * m.m[1][0] + in.z * m.m[2][0] + in.w * m.m[3][0];
		v.y = in.x * m.m[0][1] + in.y * m.m[1][1] + in.z * m.m[2][1] + in.w * m.m[3][1];
		v.z = in.x * m.m[0][2] + in.y * m.m[1][2] + in.z * m.m[2][2] + in.w * m.m[3][2];
		// din cauza ca am pus out.x in loc de out.z m-am futut 2 ore
		v.w = in.x * m.m[0][3] + in.y * m.m[1][3] + in.z * m.m[2][3] + in.w * m.m[3][3];

		// diding out parameters by w , because now we work with 4D, but we work with 3D Cartesian space

		return v;
	}
	mat4x4 Matrix_MakeIdentity()
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationX(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[1][2] = sinf(fAngleRad);
		matrix.m[2][1] = -sinf(fAngleRad);
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationY(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][2] = sinf(fAngleRad);
		matrix.m[2][0] = -sinf(fAngleRad);
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationZ(float fAngleRad)
	{
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][1] = sinf(fAngleRad);
		matrix.m[1][0] = -sinf(fAngleRad);
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeTranslation(float x, float y, float z)
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		matrix.m[3][0] = x;
		matrix.m[3][1] = y;
		matrix.m[3][2] = z;
		return matrix;
	}
	mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspecRatio, float f_ZNear, float f_ZFar)
	{
		mat4x4 matrix;
		//float fFovRad = 1.0f / tan( fFov / 2);  // Theory
		float pi = 3.14159f;
		// Practice - converting from degrees to rads
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5 / 180.0f * 3.14159f);
		// Setting Projection Matrix
		matrix.m[0][0] = fAspecRatio * fFovRad;
		matrix.m[1][1] = fFovRad;
		matrix.m[2][2] = f_ZFar / (f_ZFar - f_ZNear);
		matrix.m[2][3] = 1.0f;
		matrix.m[3][2] = (-f_ZFar * f_ZNear) / (f_ZFar - f_ZNear);
		matrix.m[3][3] = 0.0f;
		return matrix;
	}
	mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
	{
		mat4x4 matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
		return matrix;
	}
	mat4x4 Matrix_View(mat4x4 &m) // Just invertes the PoinAt matrix
	{ 
		mat4x4 matrix;
		matrix.m[0][0] = m.m[0][0] ; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = m.m[3][0];
		matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = m.m[3][1];
		matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = m.m[3][2];
		matrix.m[3][0] = -(m.m[3][0] * m.m[0][0] + m.m[3][1] * m.m[0][1] + m.m[3][2] * m.m[0][2]);
		matrix.m[3][1] = -(m.m[3][0] * m.m[1][0] + m.m[3][1] * m.m[1][1] + m.m[3][2] * m.m[1][2]);
		matrix.m[3][2] = -(m.m[3][0] * m.m[2][0] + m.m[3][1] * m.m[2][1] + m.m[3][2] * m.m[2][2]);
		matrix.m[3][3] = m.m[3][3];
		return matrix;
	}
	
	// To revisit the newUp computing, didn't understand it
	mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up)
	{
		vec3d newForward = Vector_Sub(target, pos);
		newForward = Vector_Normalise(newForward);

		// Pitch = X axis

		// Verifying how similar newForward and up are, and adjusting it to the newUp
		// The Up vector Could change, due to a change in pitch
		// If changed, adjusting this to be again orthogonal to the newForward
		vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
		vec3d newUp = Vector_Sub(up, a);  // !!! Didn't understand the computing of newUp
		newUp = Vector_Normalise(newUp);

		vec3d newRight = Vector_CrossProduct(newUp, newForward);

		// Construct Dimensioning and Translation Matrix ! Must be inverted !
		mat4x4 matrix;
		// It's not like in theory : Forward = z , Right = x, Up = y -- Our convention
		matrix.m[0][0] = newRight.x;	 matrix.m[0][1] = newRight.y;	  matrix.m[0][2] = newRight.z;	  matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = newUp.x;		 matrix.m[1][1] = newUp.y;		  matrix.m[1][2] = newUp.z;		  matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = newForward.x;	 matrix.m[2][1] = newForward.y;  matrix.m[2][2] = newForward.z;  matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = pos.x;			 matrix.m[3][1] = pos.y;		  matrix.m[3][2] = pos.z;		  matrix.m[3][3] = 1.0f;
		return matrix;

	}


	vec3d Vector_Add(vec3d& v1, vec3d& v2)
	{
		return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	}
	vec3d Vector_Sub(vec3d& v1, vec3d& v2)
	{
		return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
	}
	vec3d Vector_Mul(vec3d& v1, float k)
	{
		return { v1.x * k, v1.y * k, v1.z * k };
	}
	vec3d Vector_Div(vec3d& v1, float k)
	{
		return { v1.x / k, v1.y / k, v1.z / k };
	}
	float Vector_DotProduct(vec3d& v1, vec3d& v2)
	{
		return  v1.x*v2.x +  v1.y*v2.y + v1.z*v2.z;
	}
	float Vector_Length(vec3d& v)
	{
		return sqrtf(Vector_DotProduct(v, v));
	}
	vec3d Vector_Normalise(vec3d& v)
	{
		float l = Vector_Length(v);
		return { v.x / l, v.y / l, v.z / l };

	}
	vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2)
	{
		vec3d v;
		v.x = v1.y * v2.z - v1.z * v2.y;
		v.y = v1.z * v2.x - v1.x * v2.z;
		v.z = v1.x * v2.y - v1.y * v2.x;
		return v;
	}

	// To revisit - How do lines interescts with planes
	vec3d Vector_IntersectPlane(vec3d &plane_p, vec3d &plane_n, vec3d &lineStart, vec3d &lineEnd)
	{
		// plane_p = point on the plane
		// plane_n = the normal of that plane

		plane_n = Vector_Normalise(plane_n);

		float plane_d = -Vector_DotProduct(plane_n, plane_p);
		float ad = Vector_DotProduct(lineStart, plane_n);
		float bd = Vector_DotProduct(lineEnd, plane_n);

		float t = (-plane_d - ad) / (bd - ad);
		vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
		vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
		
		return Vector_Add(lineStart, lineToIntersect);
	}
	
	// To revisit - What the hell is happening ???
	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n,triangle &in_tri, triangle &out_tri1, triangle &out_tri2)
	{
		plane_n = Vector_Normalise(plane_n);

		auto dist = [&](vec3d& p)  // lambda function
			{
				vec3d n = Vector_Normalise(p);
				// dot product plane_n ** p
				return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
			};

		vec3d* inside_points[3]; int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;

		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
		if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
		if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
		else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

		if (nInsidePointCount == 0)
		{
			return 0;
		}

		if (nInsidePointCount == 3)
		{
			out_tri1 = in_tri;

			return 1;
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{

			// For Debugging
			out_tri1.col = FG_BLUE;


			// copy the appeareance
			//out_tri1.col = in_tri.col; 
			out_tri1.sym = in_tri.sym;

			out_tri1.p[0] = *inside_points[0];

			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

			return 1;

		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)// quad
		{
			// For Debugging
			out_tri1.col = FG_RED;

			//out_tri1.col = in_tri.col; 
			out_tri1.sym = in_tri.sym;

			// For Debugging
			out_tri2.col = FG_RED; 

			//out_tri2.col = in_tri.col;
			out_tri2.sym = in_tri.sym;

			// just splitting the triangle into 2 new one
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);


			out_tri2.p[0] = *inside_points[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

			return 2;
		}


	}
	CHAR_INFO GetColour(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(lum * 13.0f);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;   // The blackest one

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;   // The most white
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}


public:
	bool OnUserCreate() override
	{
		Cube.LoadFromObjectFile("Mountains.txt");

		// Projection Matrix changes only once because the screen dimension and FOV will not change in the process
		float fNear = 0.1f;
		float fFar = 1000.0f;
		float fFov_Degrees = 120.0f;   // theta 
		float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
		// Setting Projection Matrix

		matProj = Matrix_MakeProjection(fFov_Degrees, fAspectRatio, fNear, fFar);
		return true;
	}



	bool OnUserUpdate(float fElapsedTime) override
	{
		// User Input

		if (GetKey(VK_UP).bHeld)
		{
			vCamera.y += 8.0f * fElapsedTime;
		}

		if (GetKey(VK_DOWN).bHeld)
		{
			vCamera.y -= 8.0f * fElapsedTime;
		}

		// Why are they still inverted?
		if (GetKey(VK_LEFT).bHeld)
		{
			//vCamera.x -= 8.0f * fElapsedTime;
			vCamera.x += 8.0f * fElapsedTime;
		}
		if (GetKey(VK_RIGHT).bHeld)
		{
			//vCamera.x += 8.0f * fElapsedTime;
			vCamera.x -= 8.0f * fElapsedTime;
		}

		vec3d vForward = Vector_Mul(vLookDir, 8.0f * fElapsedTime);

		if (GetKey(L'W').bHeld)
		{
			vCamera = Vector_Add(vCamera, vForward);
		}
		if (GetKey(L'S').bHeld)
		{
			vCamera = Vector_Sub(vCamera, vForward);
		}

		
		if (GetKey(L'A').bHeld)  // IDK why it's inverted, in the video it's different
		{
			fYaw -= 2.0f * fElapsedTime;
		}
		if (GetKey(L'D').bHeld)
		{
			fYaw += 2.0f * fElapsedTime;
		}
									// our target vector , velocity(speed)




		// CLEAR THE SCREEN

		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		// DRAW TRIANGLES 

		mat4x4 matRotZ, matRotX;
		//fTheta += 1.0f * fElapsedTime;

		matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);
		matRotX = Matrix_MakeRotationX(fTheta);

		mat4x4 matTrans;
		matTrans = Matrix_MakeTranslation(0.0f, -10.0f, 5.0f);

		mat4x4 matWorld;
		matWorld = Matrix_MakeIdentity();
		matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
		matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);


	
		vec3d vUp = { 0,1,0 };  // Are they normalised by default ??
		vec3d vTarget = { 0,0,1 };  // Are they normalised by default ??
		mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
		vLookDir = Matrix_MultiplyVector(matCameraRot, vTarget);
		vTarget = Vector_Add(vCamera, vLookDir);  // offsetting to the current camera location


		mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);

		// Making view matrix from camera
		mat4x4 matView = Matrix_View(matCamera); // Not LookAt  ?


		std::vector<triangle> vecTrianglesToRaster;


		// Rotation =  POINT AT SYSTEM


		// Draw Triangles
		for (auto tri : Cube.tris)
		{
			triangle triProjected, triTransformed, triViewed;


			// Applying all transformation,  Rotation, Translation
			triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
			triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
			triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);

			// Getting the normal of the triangle
			
				//- Triangle Normal
				vec3d normal, line1, line2;

				//- Getting the lines
				line1 = Vector_Sub(triTransformed.p[1], triTransformed.p[0]);   // Achtung - The order of line is very important, The direction of normal could be inverted
				line2 = Vector_Sub(triTransformed.p[2], triTransformed.p[0]);

				//- Cross product
				normal = Vector_CrossProduct(line1, line2);

				//- Normalize the normal
				normal = Vector_Normalise(normal);

			// Getting the distance between Camera and Triangle  
			vec3d vCameraRay = Vector_Sub(triTransformed.p[0], vCamera);


			// if (normal.z < 0)  // didn't understand fully this, as camera could be in different locations, and it wouldn't work more
			if ( Vector_DotProduct(normal, vCameraRay) < 0.0f ) // comparing Camera projection vector and triangle's normal. Could be any point of tri
			{

				// Illumination
				vec3d light_direction = { 0.0f, 1.0f, -1.0f }; // shining from the camera
				// Normalising the vector
				light_direction = Vector_Normalise(light_direction);
				float dp = max(0.1f, Vector_DotProduct(light_direction, normal)); // could be done without max() function

				
				// if it faces the light direction, it should be 1, so black. Why it is the opposite?

				CHAR_INFO c = GetColour(dp); // The value is normalized , so dp = (0,1)   
				// std::cout << "DP = " << dp << std::endl;
				triTransformed.col = c.Attributes;
				triTransformed.sym = c.Char.UnicodeChar;


				// Convert		World Space --> View Space
				triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
				triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
				triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);
				triViewed.col = triTransformed.col;
				triViewed.sym = triTransformed.sym;

				// Clip Viewed Triangle against near plane - could form 2 additional triangles
				int nClippedTriangles = 0;
				triangle clipped[2];
				//nClippedTriangles = Triangle_ClipAgainstPlane( { 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1] ) ;   // Default
				// For testing purpose
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

				for (int n = 0; n < nClippedTriangles; n++)
				{
					//// Project triangles from   3D -> 2D
					triProjected.p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);  // x Coord of triangle
					triProjected.p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
					triProjected.p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);
					triProjected.col = clipped[n].col;
					triProjected.sym = clipped[n].sym;


					// Normalising the Projected vector

					triProjected.p[0] = Vector_Div(triProjected.p[0], triProjected.p[0].w);
					triProjected.p[1] = Vector_Div(triProjected.p[1], triProjected.p[1].w);
					triProjected.p[2] = Vector_Div(triProjected.p[2], triProjected.p[2].w);


					// X,Y = inverted . Fix it
					triProjected.p[0].x *= -1.0f;
					triProjected.p[1].x *= -1.0f;
					triProjected.p[2].x *= -1.0f;
					triProjected.p[0].y *= -1.0f;
					triProjected.p[1].y *= -1.0f;
					triProjected.p[2].y *= -1.0f;





					//// before that moment, our coordinates are between 1 and -1, normalized. We should scale them

					//// Scale into view


					// Offset verts into visible normalised space (0,2)
					vec3d vOffsetView = { 1,1,0 };
					triProjected.p[0] = Vector_Add(triProjected.p[0], vOffsetView);
					triProjected.p[1] = Vector_Add(triProjected.p[1], vOffsetView);
					triProjected.p[2] = Vector_Add(triProjected.p[2], vOffsetView);


					triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[2].y *= 0.5f * (float)ScreenHeight();




					vecTrianglesToRaster.push_back(triProjected);
				}


				//Can't use those functions anymore because the triangles are not drawn in a specific order. It happens that some triangles that are behind are drawn in top of the front 
				/*FillTriangle(
					triProjected.p[0].x, triProjected.p[0].y,
					triProjected.p[1].x, triProjected.p[1].y,
					triProjected.p[2].x, triProjected.p[2].y,
					triProjected.sym,  triProjected.col
				);
				DrawTriangle(   // just for the outline
					triProjected.p[0].x, triProjected.p[0].y,
					triProjected.p[1].x, triProjected.p[1].y,
					triProjected.p[2].x, triProjected.p[2].y,
					PIXEL_SOLID, FG_BLACK
				);*/
			}
		}


		// Sorting triangles from back to front
		std::sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(),
			[](triangle& t1, triangle& t2) // lambda function
			{
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});



		for (auto& triToRaster : vecTrianglesToRaster)
		{
			// Doing Clipping against our screen, all 4 directions

			triangle clipped[2];
			std::list<triangle> listTriangles;
			listTriangles.push_back(triToRaster);
			int nNewTriangles = 1;

			for(int p = 0; p < 4; p++)
			{ 
			
				int nTrisToAdd = 0;

				while (nNewTriangles > 0)
				{
					triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;


					switch (p)
					{
						// To Revisit
					case 0: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;  // Y  up  plane
					case 1: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float) ScreenHeight() - 1, 0.0f}, {0.0f, -1.0f, 0.0f}, test, clipped[0], clipped[1]); break;  // Y down plane
					case 2: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;  // X left plane
					case 3: nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, test, clipped[0], clipped[1]); break;  // X right plane


					}

					for (int w = 0; w < nTrisToAdd; w++)
					{
						listTriangles.push_back(clipped[w]);
					}
				}
				nNewTriangles = listTriangles.size();
			}


			for (auto& t : listTriangles)
			{
				FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y,  t.sym , t.col);


				// For Debugging
				DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y,  PIXEL_SOLID, FG_BLACK);

			}


			//FillTriangle(
			//		triProjected.p[0].x, triProjected.p[0].y,
			//		triProjected.p[1].x, triProjected.p[1].y,
			//		triProjected.p[2].x, triProjected.p[2].y,
			//		triProjected.sym,  triProjected.col
			//	);

			//DrawTriangle(   // just for the outline
			//	triProjected.p[0].x, triProjected.p[0].y,
			//	triProjected.p[1].x, triProjected.p[1].y,
			//	triProjected.p[2].x, triProjected.p[2].y,
			//	PIXEL_SOLID, FG_BLACK
			//); 
		
		}

		return true;
	}


};

int main()
{
	olcEngine3D demo;
	if (demo.ConstructConsole(300, 200, 4, 4))
		demo.Start();
	return 0;
}



// Current(Defalut) path: O:\OpenGL\OverView\Project1
// 
//wchar_t buffer[MAX_PATH];
//GetCurrentDirectory(MAX_PATH, buffer);
//std::wcout << L"Current path: " << buffer << std::endl;




// !!!      Problems

// When we approach objects, our Z value gets smaller and smaller, approaching to 0
// In result, our Triangles size get bigger and bigger, hence the artifacts, memory corruption

// When camera moves forward too much, on the corners of x axis appears a black triangles,lines,rectangles , and it consumes too much detail of the scene ???


// !!!		QUESTIONS

// Are we doing the clipping because the triangles are too damn big?