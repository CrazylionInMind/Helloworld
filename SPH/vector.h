#pragma once
#include<math.h>
class Vector//3ά����
{
public:
	Vector();//���캯��
	Vector(float a, float b, float c);//���캯��

	Vector V_normal();//��һ��
	float getX() { return x; }
	float getY() { return y; }
	float getZ() { return z; }
	void SetX(float r_x) { x = r_x; }
	void SetY(float r_y) { y = r_y; }
	void SetZ(float r_z) { z = r_z; }
	~Vector();
	//�������
	Vector Vector::operator+(const Vector &b) const;
	Vector Vector::operator-(const Vector &b) const;
	Vector Vector::operator*(float b) const;
	float dot(const Vector &b) const;
	Vector Vector::operator*(const Vector &b) const;
private:
	float x, y, z;

};


