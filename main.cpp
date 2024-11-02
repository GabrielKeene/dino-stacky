#include <vector>
#include <algorithm>
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/Text.hpp>
#include <fstream>
#include <strstream>
#include <iostream>
#include <time.h>

using namespace std;

float gnd_height=0.5f;
float rotateZ=0;
struct vec3D
{
	float x = 0;
	float y = 0;
	float z = 0;
};
float Vec3Dlength(vec3D in) {
	return sqrtf(in.x * in.x + in.y * in.y + in.z * in.z);
}
struct color {
	int r = 0;
	int g=0;
	int b = 0;
};
struct triangle
{
	vec3D p[3];
	float lum = -1;
	color col ;
};
struct mesh
{
	vector<triangle> tris;
};
bool LoadObjFile(vector<triangle> &tris,std::string sFilename) {

	ifstream f(sFilename);
	if (!f.is_open()) {
		return false;
	}
	vector<vec3D> vertex;
	while (!f.eof()) {
		char line[128];
		f.getline(line, 128);
		std::strstream s;
		s << line;

		char trash;

		if (line[0] == 'v') {
			vec3D v;
			s >> trash >> v.x >> v.y >> v.z;
			vertex.push_back(v);
		}if (line[0] == 'f') {
			triangle tri;
			int a, b, c;
			s >> trash >> a >> b >> c;
			tri.p[0] = vertex[a - 1];
			tri.p[1] = vertex[b - 1];
			tri.p[2] = vertex[c - 1];
			tris.push_back(tri);
		}
	}

	return true;
}
void MoveObj(vector<triangle>& tris,double x, double y, double z) {
	for (int i = 0; i < tris.size(); i++) {
		for (int j = 0; j < 3; j++) {
			tris[i].p[j].x += x;
			tris[i].p[j].y += y;
			tris[i].p[j].z += z;
		}
	}
}
void ResizeObj(vector<triangle>& tris,double s) {

	for (int i = 0; i < tris.size(); i++) {
		for (int j = 0; j < 3; j++) {
			tris[i].p[j].x /= s;
			tris[i].p[j].y /= s;
			tris[i].p[j].z /= s;
		}
	}
}
void SetTerrainCol(vector<triangle>& tris) {
	for (int i = 0; i < tris.size(); i++) {
		color mix;
		float p = max(tris[i].p[2].y, max(tris[i].p[0].y, tris[i].p[1].y));
		p -= min(tris[i].p[2].y, min(tris[i].p[0].y, tris[i].p[1].y));
		p = min(p / 20, (float)1);
		p = p * p;
		mix.r = p * 153;
		mix.g = (1 - p) * 153 + p * 102;
		mix.b = (1 - p) * 51 + p * 51;
		tris[i].col = mix;
	}
}
void SetCol(vector<triangle>& tris,color col) {
	for (int i = 0; i < tris.size(); i++) {
		tris[i].col = col;
	}
}
void SetColStacky(vector<triangle>& tris) {
	for (int i = 0; i < tris.size(); i++) {
		if (tris[i].p[0].y < -22) {
			tris[i].col = { 200,20,10 };
			if (tris[i].p[0].y < -68) {
				tris[i].col = { 250,250,250 };
			}
		}
		if (tris[i].p[0].z < -22 && tris[i].p[0].x < 8 * 2) {
			tris[i].col = { 200,20,10 };
		}

	}
}
void RotateObjY(vector<triangle>& tris,double angle) {
	for (int i = 0; i < tris.size(); i++) {
		for (int j = 0; j < 3; j++) {
			double x = tris[i].p[j].x;
			double z = tris[i].p[j].z;
			tris[i].p[j].x = x * cos(angle) + z * sin(angle);
			tris[i].p[j].z = z * cos(angle) - x * sin(angle);
		}
	}
}
void RotateObjZ(vector<triangle>& tris,double angle) {
	for (int i = 0; i < tris.size(); i++) {
		for (int j = 0; j < 3; j++) {
			double y = tris[i].p[j].y;
			double x = tris[i].p[j].x;
			tris[i].p[j].y = y * cos(angle) + x * sin(angle);
			tris[i].p[j].x = x * cos(angle) - y * sin(angle);
		}
	}
}
struct mat4x4
{
	float m[4][4] = { 0 };
};
vec3D Vec3DAddfloat(const vec3D& a, const float& b)
{
	vec3D c;
	c.x = a.x + b;
	c.y = a.y + b;
	c.z = a.z + b;
	return c;
}
vec3D Vec3DOverFloat(const vec3D& a, const float& b)
{
	vec3D c;
	c.x = a.x / b;
	c.y = a.y / b;
	c.z = a.z / b;
	return c;
}
vec3D Vec3DTimeFloat(const vec3D& a, const float& b)
{
	vec3D c;
	c.x = a.x * b;
	c.y = a.y * b;
	c.z = a.z * b;
	return c;
}
vec3D Vec3DAddVec3D(const vec3D& a, const vec3D& b)
{
	vec3D c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
}
vec3D Vec3DSubVec3D(const vec3D& a, const vec3D& b)
{
	vec3D c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}
vec3D Vec3DTimeVec3D(const vec3D& a, const vec3D& b)
{
	vec3D c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	return c;
}
float Vec_Dot_Product(vec3D& a, vec3D& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
void MultiplyMatrixVector(vec3D& i, vec3D& o, mat4x4 m)
{
	o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
	o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
	o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
	float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

	if (w != 0.0f)
	{
		o.x /= w; o.y /= w; o.z /= w;
	}
}
mat4x4 PointAtMatrix(vec3D& pos, vec3D& target, vec3D& up) {

	vec3D newForward = Vec3DSubVec3D( target , pos);
	newForward = Vec3DOverFloat(newForward , Vec3Dlength(newForward));

	vec3D a = Vec3DTimeFloat(newForward , Vec_Dot_Product(up, newForward));
	vec3D newUp = Vec3DSubVec3D(up , a);
	newUp = Vec3DOverFloat(newUp , Vec3Dlength(newUp));

	vec3D newRight = Vec3DTimeVec3D(newUp , newForward);

	mat4x4 mat;
	mat.m[0][0] = newRight.x;	mat.m[0][1] = newRight.y;	mat.m[0][2] = newRight.z;	mat.m[0][3] = 0.0f;
	mat.m[1][0] = newUp.x;		mat.m[1][1] = newUp.y;		mat.m[1][2] = newUp.z;		mat.m[1][3] = 0.0f;
	mat.m[2][0] = newForward.x;	mat.m[2][1] = newForward.y;	mat.m[2][2] = newForward.z;	mat.m[2][3] = 0.0f;
	mat.m[3][0] = pos.x;		mat.m[3][1] = pos.y;		mat.m[3][2] = pos.z;		mat.m[3][3] = 1.0f;
	return mat;
}
mat4x4 inverseGaussJordan4x4(const mat4x4& matrice) {
	int taille = 4;

	vector<std::vector<float>> augmentee(taille, std::vector<float>(2 * taille, 0.0f));
	for (int i = 0; i < taille; ++i) {
		for (int j = 0; j < taille; ++j) {
			augmentee[i][j] = matrice.m[i][j];
		}
		augmentee[i][i + taille] = 1.0f;
	}

	for (int i = 0; i < taille; ++i) {
		int pivot = i;
		for (int j = i + 1; j < taille; ++j) {
			if (abs(augmentee[j][i]) > abs(augmentee[pivot][i])) {
				pivot = j;
			}
		}

		if (pivot != i) {
			std::swap(augmentee[i], augmentee[pivot]);
		}

		float valeurPivot = augmentee[i][i];
		if (valeurPivot != 0.0f) {
			for (int j = 0; j < 2 * taille; ++j) {
				augmentee[i][j] /= valeurPivot;
			}
		}

		for (int j = 0; j < taille; ++j) {
			if (j != i) {
				float facteur = augmentee[j][i];
				for (int k = 0; k < 2 * taille; ++k) {
					augmentee[j][k] -= facteur * augmentee[i][k];
				}
			}
		}
	}

	mat4x4 inverse;
	for (int i = 0; i < taille; ++i) {
		for (int j = 0; j < taille; ++j) {
			inverse.m[i][j] = augmentee[i][j + taille];
		}
	}

	return inverse;
}
mat4x4 createRotationMatrixZ(float angle) {
	mat4x4 matRotZ;
	matRotZ.m[0][0] = cosf(angle);
	matRotZ.m[0][1] = sinf(angle);
	matRotZ.m[1][0] = -sinf(angle);
	matRotZ.m[1][1] = cosf(angle);
	matRotZ.m[2][2] = 1;
	matRotZ.m[3][3] = 1;
	return matRotZ;
}
mat4x4 createRotationMatrixX(float angle) {
	mat4x4 matRotX;
	matRotX.m[0][0] = 1;
	matRotX.m[1][1] = cosf(angle);
	matRotX.m[1][2] = sinf(angle);
	matRotX.m[2][1] = -sinf(angle);
	matRotX.m[2][2] = cosf(angle);
	matRotX.m[3][3] = 1;
	return matRotX;
}
mat4x4 createRotationMatrixY(float angle) {
	mat4x4 matRotY;
	matRotY.m[0][0] = cosf(angle);
	matRotY.m[0][2] = sinf(angle);
	matRotY.m[1][1] = 1;
	matRotY.m[2][0] = -sinf(angle);
	matRotY.m[2][2] = cosf(angle);
	matRotY.m[3][3] = 1;
	return matRotY;
} 
mat4x4 createTranslationMatrix(float x, float y, float z)
{
	mat4x4 translationMatrix;
	translationMatrix.m[0][0] = 1.0f;
	translationMatrix.m[1][1] = 1.0f;
	translationMatrix.m[2][2] = 1.0f;
	translationMatrix.m[3][3] = 1.0f;
	translationMatrix.m[3][0] = x;
	translationMatrix.m[3][1] = y;
	translationMatrix.m[3][2] = z;
	return translationMatrix;
}
mat4x4 createProjectionMatrix(float fNear, float fFar, float fFov, float fAspectRatio) {
	mat4x4 matProj;

	float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);

	matProj.m[0][0] = fAspectRatio * fFovRad;
	matProj.m[1][1] = fFovRad;
	matProj.m[2][2] = fFar / (fFar - fNear);
	matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matProj.m[2][3] = 1.0f;
	matProj.m[3][3] = 0.0f;
	return matProj;
}
float max(float& a, float& b) {
	if (a > b)
		return a;
	return b;
}
mat4x4 MultiplyMatrixMatrix(mat4x4 m1, mat4x4 m2)
{
	mat4x4 matrix; 
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrix.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] + m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
	return matrix;
}
void NormalizeVector(vec3D& i) {
	float l = Vec3Dlength(i);
	i.x /= l; i.y /= l; i.z /= l;
}
vec3D IntersectionPlaneVector(vec3D& plane_point, vec3D& plane_normal, vec3D& lineA, vec3D& lineB) {
	NormalizeVector(plane_normal);
	float plane_d = -Vec_Dot_Product(plane_normal, plane_point);
	float ad = Vec_Dot_Product(lineA, plane_normal);
	float bd = Vec_Dot_Product(lineB, plane_normal);
	float t = (-plane_d - ad) / (bd - ad);
	vec3D StartToEnd = Vec3DSubVec3D(lineB , lineA);
	vec3D lineToIntersect = Vec3DTimeFloat(StartToEnd , t);
	return Vec3DAddVec3D(lineA , lineToIntersect);
}
int clipping(vec3D point1, vec3D point2, vec3D point3) {
	int pointsOutScreen = 0;
	if (point1.x < -1 || point1.x>1 || point1.y < -1 || point1.y>1)
		pointsOutScreen++;
	if (point2.x < -1 || point2.x>1 || point2.y < -1 || point2.y>1)
		pointsOutScreen++;
	if (point3.x < -1 || point3.x>1 || point3.y < -1 || point3.y>1)
		pointsOutScreen++;
	return pointsOutScreen;

}
float dist_clip(vec3D point,vec3D normal,vec3D p)
{
	vec3D n = p;
	n.x = n.x / Vec3Dlength(n);
	n.y = n.y / Vec3Dlength(n);
	n.z = n.z / Vec3Dlength(n);

	return (normal.x * p.x + normal.y * p.y + normal.z * p.z - Vec_Dot_Product(normal, point));
}
int clipTriangle(vec3D point, vec3D normal, triangle& tri, triangle& tri_out1, triangle& tri_out2) {

	NormalizeVector(point);
	NormalizeVector(normal);

	vec3D* point_inside[3];
	vec3D* point_outside[3];
	int nb_point_inside = 0;
	int nb_point_outside = 0;
	
	float dis1 = dist_clip(point,normal, tri.p[0]);
	float dis2 = dist_clip(point,normal, tri.p[1]);
	float dis3 = dist_clip(point,normal, tri.p[2]);

	if (dis1 > 0) {
		point_inside[nb_point_inside] = &tri.p[0];
		nb_point_inside++;
	}else{
		point_outside[nb_point_outside] = &tri.p[0];
		nb_point_outside++;
	}
	
	if (dis2 > 0) {
		point_inside[nb_point_inside] = &tri.p[1];
		nb_point_inside++;
	}	else {
		point_outside[nb_point_outside] = &tri.p[1];
		nb_point_outside++;
	}

	if (dis3 > 0) {
		point_inside[nb_point_inside] = &tri.p[2];
		nb_point_inside++;
	}	else {
		point_outside[nb_point_outside] = &tri.p[2];
		nb_point_outside++;
	}


	if (nb_point_inside == 0 || nb_point_inside == 3) {
		tri_out1 = tri;


		return (int)(nb_point_inside == 3);
	}
	
	if (nb_point_inside == 1) {
		tri_out1.col = tri.col;
		tri_out1.p[0] = *point_inside[0];

		tri_out1.p[1] = IntersectionPlaneVector(point, normal, *point_inside[0], *point_outside[0]);
		tri_out1.p[2] = IntersectionPlaneVector(point, normal, *point_inside[0], *point_outside[1]);
		
		return 1;
	}
	if (nb_point_inside == 2) {
		tri_out1.col = tri.col;
		tri_out2.col = tri.col;

		tri_out1.p[0] = *point_inside[0];
		tri_out1.p[1] = *point_inside[1];
		tri_out1.p[2] = IntersectionPlaneVector(point, normal, *point_inside[0], *point_outside[0]);
		

		tri_out2.p[0] = *point_inside[1];
		tri_out2.p[1] = tri_out1.p[2];
		tri_out2.p[2] = IntersectionPlaneVector(point, normal, *point_inside[1], *point_outside[0]);
		
		return 2;
	}
	return -1;
}
void AddCactus(vector<mesh> &gnd) {

	mesh cactus;
	LoadObjFile(cactus.tris,"cactus1.obj");
	gnd.push_back(cactus);
	ResizeObj(gnd[gnd.size() - 1].tris,(rand() % 2 +1)*0.05);
	SetCol(gnd[gnd.size() - 1].tris,{ 20,200,20 });
	double x = ((rand() % 3200*4) - 1600*4);
	double z = ((rand() % 3200*4) - 1600*4);
	double shortest = 1000000000;
	double height = 0;
	MoveObj(gnd[gnd.size() - 1].tris,x, 0, z);

	int a = 0;
		for (int j = 0; j < gnd[4].tris.size(); j++) {
			triangle tri = gnd[4].tris[j];
			double dist_tri = ((tri.p[0].x - x) * (tri.p[0].x - x) + (tri.p[0].z - z) * (tri.p[0].z - z));
			if (dist_tri < shortest) {
				shortest = dist_tri;

				height = tri.p[0].y;
				a = j;
			}
		}
		MoveObj(gnd[gnd.size() - 1].tris, 0, height, 0);
}

int main()
{
	FILE* save_read;
	fopen_s(&save_read, "save.txt", "r");
	int max_score = 0;
	fscanf_s(save_read, "%d", &max_score);
	cout << max_score << endl;
	double score = 0;
	fclose(save_read);


start:
	if (1) {
		bool night = false;
		FILE* save_w;
		fopen_s(&save_w, "save.txt", "w+");
		srand(time(NULL));
		sf::Text Score;
		sf::Text HiScore;

		sf::Font font;
		font.loadFromFile("Roboto-Bold.ttf");
		Score.setFont(font);
		Score.setString("Hello world");
		Score.setCharacterSize(24);
		Score.setFillColor(sf::Color::Red);
		Score.setStyle(sf::Text::Bold | sf::Text::Underlined);
		Score.setPosition({ 0,0 });

		HiScore.setFont(font);
		HiScore.setString("HI: "+to_string(max_score));
		HiScore.setCharacterSize(24);
		HiScore.setFillColor(sf::Color::Red);
		HiScore.setStyle(sf::Text::Bold | sf::Text::Underlined);
		HiScore.setPosition({0,50 });

		mesh Gnd_in;
		mesh Gnd_out;
		mat4x4 matProj = createProjectionMatrix(0.1f, 1000.0f, 90.0f, 1080.f / 1920.f);
		vec3D vCam;
		vec3D VCam;
		vec3D aCam;
		vCam = { 0,-200,0 };
		vec3D g = { 0.0f,9.81f,0.0f };
		aCam = g;

		vec3D vCamPointDir;
		float fCamYaw = 0;

		float fCamPitch = 0;
		string filename;
		int render_distance = 50;
		printf("Choisisez la render distance (entre 5 et 100 unites subjectives, 50 par defaut): ");
		scanf_s("%i", &render_distance);
		if (render_distance > 100 || render_distance < 5) {
			printf("He ho, un entier entre 5 et 100 j'avais dit \n");
			render_distance = 50;
		}
		int diff = 1;
		printf("Choisisez la difficulte entre 1 et 5 (1 par defaut): ");
		scanf_s("%i", &diff);
		printf("Attention, un stacki peut cacher un cactus !\n");

		printf("Generating...\n");
		if (diff < 1 || diff>5)
			diff = 1;
		LoadObjFile(Gnd_in.tris,"gnd.obj");
		LoadObjFile(Gnd_out.tris,"gnd_out.obj");

		color green;
		green = { 0,153,51 };
		vector<mesh> gnd;
		int inv_size = 2;

		for (int i = -1; i < 2; i++) {

			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) {
					gnd.push_back(Gnd_in);

				}
				else {
					gnd.push_back(Gnd_out);
				}
				ResizeObj(gnd[i * 3 + j + 4].tris,1);

				MoveObj(gnd[i * 3 + j + 4].tris,i * 800 * 4 * 4, 0, j * 4 * 800 * 4);
				SetTerrainCol(gnd[i * 3 + j + 4].tris);


			}
		}

		mesh stackiBody;

		LoadObjFile(stackiBody.tris, "stackiBody.obj");

		ResizeObj(stackiBody.tris,0.0625 * 0.25);
		SetCol(stackiBody.tris,{ 148,67,37 });

		SetColStacky(stackiBody.tris);

		gnd.push_back(stackiBody);


		for (int i = 0; i < 100*diff; i++) {
			AddCactus(gnd);
		}

		printf("Loading...\n");
		float fTheta = 3.1415926539 * 2;
		sf::ContextSettings settings;
		settings.antialiasingLevel = 0;
		sf::RenderWindow window(sf::VideoMode(1920, 1080), "SFML works!", sf::Style::Default, settings);

		sf::Clock Clock;

		vec3D light = { 0.5f,-1.0f,0.2f };
		mat4x4 matTranslation = createTranslationMatrix(0.0f, 0.0f, 5.0f);
		vCam.x = 0;
		vCam.y = 0;
		float speed = 20;
		double t = -0.0001;
		double stacki_height_long = 1000;
		bool dead = false;

		while (!dead && window.isOpen())
		{
			Score.setString("Score: "+ to_string((int)score));
			if ((int)score % 1000 < 500) {
				night = false;
			}
			else {
				night = true;
			}
			float fElapsedTime = Clock.restart().asSeconds();
			score += speed /diff * fElapsedTime;

			if (t > 6) {
				t = -0.0001;
			}
			if (t > 0)
				t += fElapsedTime;

			if (rand() % 100 == 12) {
				speed+=diff/5;
				if (rand() % 100 == 12) {
					AddCactus(gnd);
				}
			}
			vCam.y = max(gnd_height - 150 + ((gnd_height - 150) - vCam.y) * fElapsedTime, min(gnd_height, vCam.y + VCam.y * fElapsedTime));
			VCam.y = (VCam.y + (aCam.y * 20) * fElapsedTime);

			sf::Event event;
			vec3D vForwardCam = Vec3DTimeFloat(vCamPointDir , (8.0f * fElapsedTime));
			vForwardCam.y = 0;
			vec3D vRightCam;

			MultiplyMatrixVector(vForwardCam, vRightCam, createRotationMatrixY(3.141592f / 2.f));
			vRightCam = Vec3DTimeFloat( vRightCam , (fElapsedTime * 8.0f));
			fCamYaw = ((sf::Mouse::getPosition().x - 960) * 2) / 960.0f * 3.141592f;
			vec3D vUp = { 0,1,0 };
			vec3D vTarget = { 0,0,1 };
			if (vCam.x > 1600 * 4) {
				vCam.x -= 3200 * 4;
			}
			else if (vCam.x < -1600 * 4) {
				vCam.x += 3200 * 4;
			}
			if (vCam.z > 1600 * 4) {
				vCam.z -= 3200 * 4;
			}
			else if (vCam.z < -1600 * 4) {
				vCam.z += 3200 * 4;
			}

			vCam = Vec3DAddVec3D(vCam , Vec3DTimeFloat( vForwardCam , speed));
			for (int i = 10; i < gnd.size(); i++) {
				if (Vec3Dlength(Vec3DSubVec3D(gnd[i].tris[0].p[0] , gnd[9].tris[0].p[0])) < 50) {
					printf("He non, ce cactus n'est pas comestible");
					dead = true;
				}
			}
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::KeyPressed) {
					switch (event.key.code) {
					case sf::Keyboard::Z:
						vCam = Vec3DAddVec3D(vCam , Vec3DTimeFloat(vForwardCam , 20));
						break;
					case sf::Keyboard::Q:
						vCam = Vec3DAddVec3D(vCam , Vec3DTimeFloat (vRightCam , 30));
						break;
					case sf::Keyboard::S:
						vCam = Vec3DSubVec3D(vCam , Vec3DTimeFloat(vForwardCam , 20));
						break;
					case sf::Keyboard::D:
						vCam = Vec3DSubVec3D(vCam , Vec3DTimeFloat(vRightCam , 10));
						break;
					case sf::Keyboard::Space:
						if (t < 0)
							t = fElapsedTime;
						break;
					case sf::Keyboard::LShift:
						vCam = Vec3DAddVec3D(vCam , Vec3DTimeFloat(vUp , 0.1f));
						break;
					}
				}

				if (event.type == sf::Event::Closed)
					window.close();
			}
			if (!night) {
				window.clear({ 0,153,255,0 });
			}
			else {
				window.clear({ 0,60,80,0 });

			}
			window.setTitle(to_string(1.0f / fElapsedTime));


			mat4x4 matRotZ, matRotX, matWorld;

			matRotZ = createRotationMatrixZ(fTheta);
			matRotX = createRotationMatrixX(fTheta);
			matWorld = MultiplyMatrixMatrix(matRotZ, matRotX);
			matWorld = MultiplyMatrixMatrix(matWorld, matTranslation);


			mat4x4 matCamRot = MultiplyMatrixMatrix(createRotationMatrixY(fCamYaw), createRotationMatrixX(fCamPitch));

			MultiplyMatrixVector(vTarget, vCamPointDir, matCamRot);
			vTarget = Vec3DAddVec3D(vCam , vCamPointDir);

			mat4x4 matCamera = PointAtMatrix(vCam, vTarget, vUp);
			mat4x4 matView = inverseGaussJordan4x4(matCamera);
			double closest = 100000000000;
			vector<triangle> vecTriangleToDraw;
			int nb_tri = 1;
			gnd[9] = stackiBody;
			RotateObjZ(gnd[9].tris,rotateZ / 700);

			RotateObjY(gnd[9].tris,-fCamYaw + 3.14159265 / 2);

			double etpaf = 5.0;
			double etpif = 0;


			MoveObj(gnd[9].tris,vCam.x - gnd[9].tris[0].p[0].x + etpif, vCam.y - gnd[9].tris[0].p[0].y, vCam.z - gnd[9].tris[0].p[0].z - etpaf);
			float shortest_sist_staky_gnd = 100000;
			double stacki_height = 314;
			double stacki_height_front = 314;

			MoveObj(gnd[9].tris, vCamPointDir.x * 100, 0, vCamPointDir.z * 100);


			for (int i = 0; i < gnd[4].tris.size(); i++) {
				float cur_dist = ((gnd[4].tris[i].p[0].x - gnd[9].tris[0].p[0].x) * (gnd[4].tris[i].p[0].x - gnd[9].tris[0].p[0].x) + (gnd[4].tris[i].p[0].z - gnd[9].tris[0].p[0].z) * (gnd[4].tris[i].p[0].z - gnd[9].tris[0].p[0].z));
				if (shortest_sist_staky_gnd > cur_dist) {
					shortest_sist_staky_gnd = cur_dist;
					double avg = gnd[9].tris[0].p[0].y + gnd[9].tris[0].p[1].y + gnd[9].tris[0].p[2].y;
					stacki_height = avg / 3 - gnd[4].tris[i].p[0].y;

				}
			}
			shortest_sist_staky_gnd = 100000;


			double dist_front = 150;
			for (int i = 0; i < gnd[4].tris.size(); i++) {
				float cur_dist = ((gnd[4].tris[i].p[0].x - gnd[9].tris[0].p[0].x - vCamPointDir.x * dist_front) * (gnd[4].tris[i].p[0].x - gnd[9].tris[0].p[0].x - vCamPointDir.x * dist_front) + (gnd[4].tris[i].p[0].z - gnd[9].tris[0].p[0].z - vCamPointDir.z * dist_front) * (gnd[4].tris[i].p[0].z - gnd[9].tris[0].p[0].z - vCamPointDir.z * dist_front));
				if (shortest_sist_staky_gnd > cur_dist) {
					shortest_sist_staky_gnd = cur_dist;
					double avg = gnd[9].tris[0].p[0].y + gnd[9].tris[0].p[1].y + gnd[9].tris[0].p[2].y;
					stacki_height_front = avg / 3 - gnd[4].tris[i].p[0].y;

				}
			}

			stacki_height_long += fElapsedTime * 20 * (stacki_height - stacki_height_long);

			MoveObj(gnd[9].tris,0, -stacki_height_long - 10, 0);
			rotateZ -= (rotateZ - (stacki_height_front - stacki_height)) * 10 * fElapsedTime;
			stacki_height += stacki_height_front;
			stacki_height /= 2;
			gnd_height = 0;
			MoveObj(gnd[9].tris, 0, min((float)0, (float)(-480 * (-(t - 0.5) * (t - 0.5) + 0.25))), 0);


			vCam.y = gnd[9].tris[0].p[0].y - 50;
			for (int a = 0; a < gnd.size(); a++) {
				mesh obj = gnd[a];
				for (int b = 0; b < obj.tris.size(); b++)
				{
					triangle tri = obj.tris[b];
					double dist_cam_tri = ((tri.p[0].x - vCam.x) * (tri.p[0].x - vCam.x) + (tri.p[0].z - vCam.z) * (tri.p[0].z - vCam.z));

					if (dist_cam_tri < 1200 * 1200 * render_distance / 40) {
						if (a < 9 && closest > dist_cam_tri) {
							closest = dist_cam_tri;
							gnd_height = (tri.p[0].y + tri.p[1].y + tri.p[2].y) / 3;
							if (dist_cam_tri < 10) {
								gnd_height += (tri.p[0].y + tri.p[1].y + tri.p[2].y) / 3;
								nb_tri++;
							}
						}
						triangle triProjected, triTranformed, triViewed;

						MultiplyMatrixVector(tri.p[0], triTranformed.p[0], matWorld);	
						MultiplyMatrixVector(tri.p[1], triTranformed.p[1], matWorld);
						MultiplyMatrixVector(tri.p[2], triTranformed.p[2], matWorld);


						vec3D line1, line2, normal;

						line1 = Vec3DSubVec3D( triTranformed.p[1] , triTranformed.p[0]);
						line2 = Vec3DSubVec3D( triTranformed.p[2] , triTranformed.p[0]);
						normal = Vec3DTimeVec3D(line1 , line2);

						NormalizeVector(normal);

						vec3D Ray = Vec3DSubVec3D(triTranformed.p[0] , vCam);

						if (/*Vec_Dot_Product(Ray, normal) < 0.0f && Vec_Dot_Product(vForwardCam, Ray) > 0.0f*/true) {
							if (tri.lum < 0) {
								vec3D lightNormalized = light;
								float l = Vec3Dlength(lightNormalized);
								lightNormalized = Vec3DOverFloat(lightNormalized , l);

								float dp = normal.x * lightNormalized.x + normal.y * lightNormalized.y + normal.z * lightNormalized.z;
								triTranformed.lum = max(0.05f, dp);
							}
							MultiplyMatrixVector(triTranformed.p[0], triViewed.p[0], matView);
							MultiplyMatrixVector(triTranformed.p[1], triViewed.p[1], matView);
							MultiplyMatrixVector(triTranformed.p[2], triViewed.p[2], matView);
							int nb_clipped_tri = 0;
							triangle clipped[2];
							nb_clipped_tri = clipTriangle({ 0.0f,0.0f,0.005f }, { 0.0f,0.0f,1.0f }, triViewed, clipped[0], clipped[1]);

							for (int i = 0; i < nb_clipped_tri; i++) {
								triViewed = clipped[i];
								MultiplyMatrixVector(triViewed.p[0], triProjected.p[0], matProj);
								MultiplyMatrixVector(triViewed.p[1], triProjected.p[1], matProj);
								MultiplyMatrixVector(triViewed.p[2], triProjected.p[2], matProj);


								triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
								triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
								triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;
								triProjected.p[0].x *= 0.5f * (float)window.getSize().x;
								triProjected.p[0].y *= 0.5f * (float)window.getSize().y;
								triProjected.p[1].x *= 0.5f * (float)window.getSize().x;
								triProjected.p[1].y *= 0.5f * (float)window.getSize().y;
								triProjected.p[2].x *= 0.5f * (float)window.getSize().x;
								triProjected.p[2].y *= 0.5f * (float)window.getSize().y;

								triProjected.lum = triTranformed.lum;
								triProjected.col = tri.col;

								vecTriangleToDraw.push_back(triProjected);
							}
						}
					}
				}
			}
			gnd_height = gnd_height / nb_tri;
			gnd_height -= 130;
			std::sort(vecTriangleToDraw.begin(), vecTriangleToDraw.end(), [](triangle& t1, triangle& t2) {
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z);
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z);
				return z2 < z1;
				});
			for (int i = 0; i < vecTriangleToDraw.size();i++) {
				triangle triProjected = vecTriangleToDraw[i];
				sf::ConvexShape triangle;
				triangle.setPointCount(3);

				triangle.setPoint(0, sf::Vector2f(triProjected.p[0].x, triProjected.p[0].y));
				triangle.setPoint(1, sf::Vector2f(triProjected.p[1].x, triProjected.p[1].y));
				triangle.setPoint(2, sf::Vector2f(triProjected.p[2].x, triProjected.p[2].y));

				triProjected.lum = max(triProjected.lum, (float)0.2);

				if (night) {
					triProjected.lum /= 3;
				}
				triangle.setFillColor(sf::Color(triProjected.col.r * triProjected.lum, triProjected.col.g * triProjected.lum, triProjected.col.b * triProjected.lum));

				window.draw(triangle);

			}
			window.draw(Score);
			window.draw(HiScore);

			window.display();
		}
		max_score = (int)max((float)score, (float)max_score);

		fprintf(save_w, "%d", (int)max_score);
		fclose(save_w);

	}



	score = 0;
	printf("\n On repart pour un tour ? (1: Oui )\n");
	int answ = 0;
	scanf_s("%i", &answ);

	if (answ == 1) {
		goto start;
	}

	return 0;
}
