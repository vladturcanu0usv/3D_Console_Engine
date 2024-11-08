//#include "olcConsoleGameEngine.h"
//
//struct vec3d
//{
//	float x, y, z;
//};
//
//struct triangle
//{
//	vec3d p[3];
//};
//
//struct mesh
//{
//	std::vector<triangle> tris;
//};
//
//
//// Template Matrix
//struct mat4x4
//{
//	float m[4][4]{ 0 };
//};
//
//
//
//class olcEngine3D : public olcConsoleGameEngine
//{
//
//
//
//public:
//	olcEngine3D()
//	{
//		m_sAppName = L"3D Demo";
//	}
//
//
//
//private:
//	mesh Cube;
//	mat4x4 matProj; // Projecion Matrix
//	float fTheta;
//
//	// passing 2 vectors to not mess with the input one
//	void MultiplyMatrixVector(vec3d& in, vec3d& out, mat4x4& m) // usual matrix multiplication
//	{
//		out.x = in.x * m.m[0][0] + in.y * m.m[1][0] + in.z * m.m[2][0] + m.m[3][0];
//		out.y = in.x * m.m[0][1] + in.y * m.m[1][1] + in.z * m.m[2][1] + m.m[3][1];
//		out.x = in.x * m.m[0][2] + in.y * m.m[1][2] + in.z * m.m[2][2] + m.m[3][2];
//		float w = in.x * m.m[0][3] + in.y * m.m[1][3] + in.z * m.m[2][3] + m.m[3][3];
//
//		// diding out parameters by w , because now we work with 4D, but we work with 3D Cartesian space
//
//		if (w != 0.0f)
//		{
//			out.x /= w; out.y /= w; out.z /= w;
//		}
//
//	}
//
//
//public:
//	bool OnUserCreate() override
//	{
//		Cube.tris =
//		{
//			// SOUTH
//
//			{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,	1.0f, 1.0f, 0.0f },
//			{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,	1.0f, 0.0f, 0.0f },
//
//			// EAST
//
//			{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,	1.0f, 1.0f, 1.0f },
//			{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,	1.0f, 0.0f, 1.0f },
//
//			// WEST
//
//			{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,	0.0f, 1.0f, 0.0f },
//			{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,	0.0f, 0.0f, 0.0f },
//
//			// NORTH
//
//			{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,	0.0f, 1.0f, 1.0f },
//			{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,	0.0f, 0.0f, 1.0f },
//
//			// TOP
//
//			{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,	1.0f, 1.0f, 1.0f },
//			{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,	1.0f, 1.0f, 0.0f },
//
//			// BOTTOM
//
//			{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,	0.0f, 0.0f, 0.0f },
//			{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,	1.0f, 0.0f, 0.0f },
//
//		};
//
//
//		// Projection Matrix changes only once because the screen dimension and FOV will not change in the process
//		float fNear = 0.1f;
//		float fFar = 1000.0f;
//		float fFov = 90.0f;   // theta 
//		float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
//
//		//float fFovRad = 1.0f / tan( fFov / 2);  // Theory
//		float pi = 3.14159f;
//		// Practice - converting from degrees to rads
//		float fFovRad = 1.0f / tanf(fFov * 0.5 / 180.0f * 3.14159f);
//
//		// Setting Projection Matrix
//
//		matProj.m[0][0] = fAspectRatio * fFovRad;
//		matProj.m[1][1] = fFovRad;
//		matProj.m[2][2] = fFar / (fFar - fNear);
//		matProj.m[2][3] = 1.0f;
//		matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
//		matProj.m[3][3] = 0.0f;
//
//		return true;
//	}
//
//
//	float count = 0.01f;
//
//	bool OnUserUpdate(float fElapsedTime) override
//	{
//		// CLEAR THE SCREEN
//
//		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);
//
//		// DRAW TRIANGLES 
//
//
//
//		mat4x4 matRotZ, matRotX;
//		fTheta += 1.0f * fElapsedTime;
//
//
//		// Rotation Z
//		matRotZ.m[0][0] = cosf(fTheta);
//		matRotZ.m[0][1] = sinf(fTheta);
//		matRotZ.m[1][0] = -sinf(fTheta);
//		matRotZ.m[1][1] = cosf(fTheta);
//		matRotZ.m[2][2] = 1;
//		matRotZ.m[3][3] = 1;
//
//		// Rotation X
//		matRotX.m[0][0] = 1;
//		matRotX.m[1][1] = cosf(fTheta * 0.5f);
//		matRotX.m[1][2] = sinf(fTheta * 0.5f);
//		matRotX.m[2][1] = -sinf(fTheta * 0.5f);
//		matRotX.m[2][2] = cosf(fTheta * 0.5f);
//		matRotX.m[3][3] = 1;
//
//
//
//
//		// Draw Triangles
//		for (auto tri : Cube.tris)
//		{
//			triangle triProjected, triTranslate, triRotate_Z, triRotate_X;
//
//
//
//			MultiplyMatrixVector(tri.p[0], triRotate_Z.p[0], matRotZ);  // x Coord of triangle
//			MultiplyMatrixVector(tri.p[1], triRotate_Z.p[1], matRotZ);  // y Coord of triangle
//			MultiplyMatrixVector(tri.p[2], triRotate_Z.p[2], matRotZ);  // z Coord of triangle
//
//			MultiplyMatrixVector(triRotate_Z.p[0], triRotate_X.p[0], matRotX);  // x Coord of triangle
//			MultiplyMatrixVector(triRotate_Z.p[1], triRotate_X.p[1], matRotX);  // y Coord of triangle
//			MultiplyMatrixVector(triRotate_Z.p[2], triRotate_X.p[2], matRotX);  // z Coord of triangle
//
//
//
//			// REVISITED CODE
//			triTranslate = triRotate_X;
//			triTranslate.p[0].z = triRotate_X.p[0].z + 3.0f;
//			triTranslate.p[1].z = triRotate_X.p[1].z + 3.0f;
//			triTranslate.p[2].z = triRotate_X.p[2].z + 3.0f;
//
//
//			//// Last thing to do is translate
//			MultiplyMatrixVector(triTranslate.p[0], triProjected.p[0], matProj);  // x Coord of triangle
//			MultiplyMatrixVector(triTranslate.p[1], triProjected.p[1], matProj);  // y Coord of triangle
//			MultiplyMatrixVector(triTranslate.p[2], triProjected.p[2], matProj);
//
//
//			//// before that moment, our coordinates are between 1 and -1, normalized. We should scale them
//
//			//// Scale into view
//
//			triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
//			triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
//			triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;  // now it's betwe
//
//
//			//// Again, we want to normalize it, ad then scale it properly to the screen size
//			triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
//			triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
//			triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
//			triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
//			triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
//			triProjected.p[2].y *= 0.5f * (float)ScreenHeight();  // now it's betwenn 0 and 2
//
//
//
//			DrawTriangle(
//				triProjected.p[0].x, triProjected.p[0].y,
//				triProjected.p[1].x, triProjected.p[1].y,
//				triProjected.p[2].x, triProjected.p[2].y,
//				PIXEL_SOLID, FG_WHITE
//			);
//
//
//		}
//
//
//		return true;
//	}
//
//
//};
//
//int main()
//{
//	olcEngine3D demo;
//	if (demo.ConstructConsole(200, 200, 4, 4))
//		demo.Start();
//	return 0;
//}








//// Second one
//#include "olcConsoleGameEngine.h"
//#include <fstream>
//#include <strstream>
//#include <windows.h> // For finding the pwd
//#include <algorithm>
//
//struct vec3d
//{
//	float x, y, z, w; // originally only was x,y,z
//};
//
//struct triangle
//{
//	vec3d p[3];
//
//	wchar_t sym;
//	short col;
//};
//
//struct mesh
//{
//	std::vector<triangle> tris;
//
//	bool LoadFromObjectFile(std::string sFilename)   // Requires revisiting
//	{
//		std::ifstream f(sFilename);
//		if (!f.is_open())
//			return false;
//
//		std::vector<vec3d> verts;
//
//		while (f.eof() == false)
//		{
//			char line[128];
//			f.getline(line, 128);
//
//			std::strstream s;
//			s << line;   // Sample : "v 1.000000 -1.000000 -1.000000"
//
//			char identifier;
//
//			if (line[0] == 'v')
//			{
//				// strstream automatically converts data separated by " ".
//				vec3d v;
//				s >> identifier >> v.x >> v.y >> v.z;
//				verts.push_back(v);
//			}
//
//			if (line[0] == 'f')
//			{
//				int f[3];
//				s >> identifier >> f[0] >> f[1] >> f[2];
//				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
//			}
//		}
//		return true;
//	}
//};
//
//// Template Matrix
//struct mat4x4
//{
//	float m[4][4]{ 0 };
//};
//
//////////////////////////////////////////////////////////////////////////////
//
//class olcEngine3D : public olcConsoleGameEngine
//{
//
//public:
//	olcEngine3D()
//	{
//		m_sAppName = L"3D Demo";
//	}
//
//private:
//	mesh Cube;
//	mat4x4 matProj; // Projecion Matrix
//
//	vec3d vCamera;
//
//	float fTheta;
//
//	// passing 2 vectors to not mess with the input one
//	vec3d MultiplyMatrixVector(mat4x4& m, vec3d& in) // usual matrix multiplication
//	{
//		vec3d v;
//		v.x = in.x * m.m[0][0] + in.y * m.m[1][0] + in.z * m.m[2][0] + in.w * m.m[3][0];
//		v.y = in.x * m.m[0][1] + in.y * m.m[1][1] + in.z * m.m[2][1] + in.w * m.m[3][1];
//		v.z = in.x * m.m[0][2] + in.y * m.m[1][2] + in.z * m.m[2][2] + in.w * m.m[3][2];
//		// din cauza ca am pus out.x in loc de out.z m-am futut 2 ore
//		v.w = in.x * m.m[0][3] + in.y * m.m[1][3] + in.z * m.m[2][3] + in.w * m.m[3][3];
//
//		// diding out parameters by w , because now we work with 4D, but we work with 3D Cartesian space
//
//		return v;
//	}
//
//
//	mat4x4 Matrix_MakeIdentity()
//	{
//		mat4x4 matrix;
//		matrix.m[0][0] = 1.0f;
//		matrix.m[1][1] = 1.0f;
//		matrix.m[2][2] = 1.0f;
//		matrix.m[3][3] = 1.0f;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MakeRotationX(float fAngleRad)
//	{
//		mat4x4 matrix;
//		matrix.m[0][0] = 1.0f;
//		matrix.m[1][1] = cosf(fAngleRad);
//		matrix.m[1][2] = sinf(fAngleRad);
//		matrix.m[2][1] = -sinf(fAngleRad);
//		matrix.m[2][2] = cosf(fAngleRad);
//		matrix.m[3][3] = 1.0f;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MakeRotationY(float fAngleRad)
//	{
//		mat4x4 matrix;
//		matrix.m[0][0] = cosf(fAngleRad);
//		matrix.m[0][2] = sinf(fAngleRad);
//		matrix.m[2][0] = -sinf(fAngleRad);
//		matrix.m[1][1] = 1.0f;
//		matrix.m[2][2] = cosf(fAngleRad);
//		matrix.m[3][3] = 1.0f;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MakeRotationZ(float fAngleRad)
//	{
//		mat4x4 matrix;
//		matrix.m[0][0] = cosf(fAngleRad);
//		matrix.m[0][1] = sinf(fAngleRad);
//		matrix.m[1][0] = -sinf(fAngleRad);
//		matrix.m[1][1] = cosf(fAngleRad);
//		matrix.m[2][2] = 1.0f;
//		matrix.m[3][3] = 1.0f;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MakeTranslation(float x, float y, float z)
//	{
//		mat4x4 matrix;
//		matrix.m[0][0] = 1.0f;
//		matrix.m[1][1] = 1.0f;
//		matrix.m[2][2] = 1.0f;
//		matrix.m[3][3] = 1.0f;
//		matrix.m[3][0] = x;
//		matrix.m[3][1] = y;
//		matrix.m[3][2] = z;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspecRatio, float f_ZNear, float f_ZFar)
//	{
//		mat4x4 matrix;
//		//float fFovRad = 1.0f / tan( fFov / 2);  // Theory
//		float pi = 3.14159f;
//		// Practice - converting from degrees to rads
//		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5 / 180.0f * 3.14159f);
//		// Setting Projection Matrix
//		matrix.m[0][0] = fAspecRatio * fFovRad;
//		matrix.m[1][1] = fFovRad;
//		matrix.m[2][2] = f_ZFar / (f_ZFar - f_ZNear);
//		matrix.m[2][3] = 1.0f;
//		matrix.m[3][2] = (-f_ZFar * f_ZNear) / (f_ZFar - f_ZNear);
//		matrix.m[3][3] = 0.0f;
//		return matrix;
//	}
//
//	mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
//	{
//		mat4x4 matrix;
//		for (int c = 0; c < 4; c++)
//			for (int r = 0; r < 4; r++)
//				matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
//		return matrix;
//	}
//
//
//
//	vec3d Vector_Add(vec3d& v1, vec3d& v2)
//	{
//		return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
//	}
//
//	vec3d Vector_Sub(vec3d& v1, vec3d& v2)
//	{
//		return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
//	}
//	vec3d Vector_Mul(vec3d& v1, float k)
//	{
//		return { v1.x * k, v1.y * k, v1.z * k };
//	}
//	vec3d Vector_Div(vec3d& v1, float k)
//	{
//		return { v1.x / k, v1.y / k, v1.z / k };
//	}
//	float Vector_DotProduct(vec3d& v1, vec3d& v2)
//	{
//		return  v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
//	}
//
//	float Vector_Length(vec3d& v)
//	{
//		return sqrtf(Vector_DotProduct(v, v));
//	}
//
//	vec3d Vector_Normalise(vec3d& v)
//	{
//		float l = Vector_Length(v);
//		return { v.x / l, v.y / l, v.z / l };
//
//	}
//
//	vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2)
//	{
//		vec3d v;
//		v.x = v1.y * v2.z - v1.z * v2.y;
//		v.y = v1.z * v2.x - v1.x * v2.z;
//		v.z = v1.x * v2.y - v1.y * v2.x;
//		return v;
//	}
//
//	CHAR_INFO GetColour(float lum)
//	{
//		short bg_col, fg_col;
//		wchar_t sym;
//		int pixel_bw = (int)(lum * 13.0f);
//		switch (pixel_bw)
//		{
//		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;   // The blackest one
//
//		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
//		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
//		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
//		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;
//
//		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
//		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
//		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
//		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;
//
//		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
//		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
//		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
//		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;   // The most white
//		default:
//			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
//		}
//
//		CHAR_INFO c;
//		c.Attributes = bg_col | fg_col;
//		c.Char.UnicodeChar = sym;
//		return c;
//	}
//
//
//public:
//	bool OnUserCreate() override
//	{
//
//		Cube.LoadFromObjectFile("VideoShip.txt");
//
//		// Projection Matrix changes only once because the screen dimension and FOV will not change in the process
//		float fNear = 0.1f;
//		float fFar = 1000.0f;
//		float fFov = 90.0f;   // theta 
//		float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
//
//		//float fFovRad = 1.0f / tan( fFov / 2);  // Theory
//		float pi = 3.14159f;
//		// Practice - converting from degrees to rads
//		float fFovRad = 1.0f / tanf(fFov * 0.5 / 180.0f * 3.14159f);
//
//		// Setting Projection Matrix
//
//		matProj = Matrix_MakeProjection(fFovRad, fAspectRatio, fNear, fFar);
//
//		return true;
//	}
//
//
//	float count = 0.01f;
//
//	bool OnUserUpdate(float fElapsedTime) override
//	{
//		// CLEAR THE SCREEN
//
//		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);
//
//		// DRAW TRIANGLES 
//
//		mat4x4 matRotZ, matRotX;
//		fTheta += 1.0f * fElapsedTime;
//
//		matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);
//		matRotX = Matrix_MakeRotationX(fTheta);
//
//
//		mat4x4 matTrans;
//		matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 16.0f);
//
//		mat4x4 matWorld;
//		matWorld = Matrix_MakeIdentity();
//		matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
//		matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);
//
//
//
//		std::vector<triangle> vecTrianglesToRaster;
//
//
//		// Draw Triangles
//		for (auto tri : Cube.tris)
//		{
//			triangle triProjected, triTransformed;
//
//
//			MultiplyMatrixVector(tri.p[0], triRotate_Z.p[0], matRotZ);  // x Coord of triangle
//			MultiplyMatrixVector(tri.p[1], triRotate_Z.p[1], matRotZ);  // y Coord of triangle
//			MultiplyMatrixVector(tri.p[2], triRotate_Z.p[2], matRotZ);  // z Coord of triangle
//
//			MultiplyMatrixVector(triRotate_Z.p[0], triRotate_X.p[0], matRotX);  // x Coord of triangle
//			MultiplyMatrixVector(triRotate_Z.p[1], triRotate_X.p[1], matRotX);  // y Coord of triangle
//			MultiplyMatrixVector(triRotate_Z.p[2], triRotate_X.p[2], matRotX);  // z Coord of triangle
//
//
//			// REVISITED CODE
//			triTranslate = triRotate_X;
//			triTranslate.p[0].z = triRotate_X.p[0].z + 6.0f;
//			triTranslate.p[1].z = triRotate_X.p[1].z + 6.0f;
//			triTranslate.p[2].z = triRotate_X.p[2].z + 6.0f;
//
//
//			vec3d normal, line1, line2;
//
//			line1.x = triTranslate.p[1].x - triTranslate.p[0].x;
//			line1.y = triTranslate.p[1].y - triTranslate.p[0].y;
//			line1.z = triTranslate.p[1].z - triTranslate.p[0].z;
//
//
//			line2.x = triTranslate.p[2].x - triTranslate.p[0].x;
//			line2.y = triTranslate.p[2].y - triTranslate.p[0].y;
//			line2.z = triTranslate.p[2].z - triTranslate.p[0].z;
//
//
//			normal.x = line1.y * line2.z - line1.z * line2.y;
//			normal.y = line1.z * line2.x - line1.x * line2.z;
//			normal.z = line1.x * line2.y - line1.y * line2.x;
//
//
//			// normalizing the normal
//
//			float l = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
//			normal.x /= l; normal.y /= l; normal.z /= l;  // works fine without normalising, just an optimisation
//
//
//			// if (normal.z < 0)  // didn't understand fully this, as camera could be in different locations, and it wouldn't work more
//			if (normal.x * (triTranslate.p[0].x - vCamera.x) +
//				normal.y * (triTranslate.p[0].y - vCamera.y) +
//				normal.z * (triTranslate.p[0].z - vCamera.z) < 0.0f) // comparing Camera projection vector and triangle's normal. Could be any point of tri
//			{
//
//				// Illumination
//				vec3d light_direction = { 0.0f, 0.0f, -1.0f }; // shining from the camera
//				// Normalising the vector
//				float l = sqrt(light_direction.x * light_direction.x + light_direction.y * light_direction.y + light_direction.z * light_direction.z);
//				light_direction.x /= l;
//				light_direction.y /= l;
//				light_direction.z /= l;
//
//				float dp = normal.x * (light_direction.x) +
//					normal.y * light_direction.y +
//					normal.z * light_direction.z;
//
//
//				// if it faces the light direction, it should be 1, so black. Why it is the opposite?
//
//
//				CHAR_INFO c = GetColour(dp); // The value is normalized , so dp = (0,1)   
//				// std::cout << "DP = " << dp << std::endl;
//				triTranslate.col = c.Attributes;
//				triTranslate.sym = c.Char.UnicodeChar;
//
//
//				//// Last thing to do is translate
//				MultiplyMatrixVector(triTranslate.p[0], triProjected.p[0], matProj);  // x Coord of triangle
//				MultiplyMatrixVector(triTranslate.p[1], triProjected.p[1], matProj);  // y Coord of triangle
//				MultiplyMatrixVector(triTranslate.p[2], triProjected.p[2], matProj);
//
//
//				triProjected.col = triTranslate.col;
//				triProjected.sym = triTranslate.sym;
//
//
//				//// before that moment, our coordinates are between 1 and -1, normalized. We should scale them
//
//				//// Scale into view
//
//
//				triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
//				triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
//				triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;  // now it's betwen 0 and 2
//
//
//				triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
//				triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
//				triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
//				triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
//				triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
//				triProjected.p[2].y *= 0.5f * (float)ScreenHeight();
//
//
//
//
//				vecTrianglesToRaster.push_back(triProjected);
//
//
//
//				//Can't use those functions anymore because the triangles are not drawn in a specific order. It happens that some triangles that are behind are drawn in top of the front 
//				/*FillTriangle(
//					triProjected.p[0].x, triProjected.p[0].y,
//					triProjected.p[1].x, triProjected.p[1].y,
//					triProjected.p[2].x, triProjected.p[2].y,
//					triProjected.sym,  triProjected.col
//				);
//				DrawTriangle(   // just for the outline
//					triProjected.p[0].x, triProjected.p[0].y,
//					triProjected.p[1].x, triProjected.p[1].y,
//					triProjected.p[2].x, triProjected.p[2].y,
//					PIXEL_SOLID, FG_BLACK
//				);*/
//			}
//		}
//
//
//		// Sorting triangles from back to front
//		std::sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(),
//			[](triangle& t1, triangle& t2) // lambda function
//			{
//				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
//				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
//				return z1 > z2;
//			});
//
//
//
//		for (auto& triProjected : vecTrianglesToRaster)
//		{
//			FillTriangle(
//				triProjected.p[0].x, triProjected.p[0].y,
//				triProjected.p[1].x, triProjected.p[1].y,
//				triProjected.p[2].x, triProjected.p[2].y,
//				triProjected.sym, triProjected.col
//			);
//
//		}
//
//		return true;
//	}
//
//
//};
//
//int main()
//{
//	olcEngine3D demo;
//	if (demo.ConstructConsole(300, 200, 4, 4))
//		demo.Start();
//	return 0;
//}
//
//
//
//// Current(Defalut) path: O:\OpenGL\OverView\Project1
//// 
////wchar_t buffer[MAX_PATH];
////GetCurrentDirectory(MAX_PATH, buffer);
//std::wcout << L"Current path: " << buffer << std::endl;
