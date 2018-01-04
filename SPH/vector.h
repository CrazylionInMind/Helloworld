#pragma once
#include<math.h>
class Vector//3维向量
{
public:
	Vector();//构造函数
	Vector(float a, float b, float c);//构造函数

	Vector V_normal();//归一化
	float getX() { return x; }
	float getY() { return y; }
	float getZ() { return z; }
	void SetX(float r_x) { x = r_x; }
	void SetY(float r_y) { y = r_y; }
	void SetZ(float r_z) { z = r_z; }
	~Vector();
	//向量相加
	Vector Vector::operator+(const Vector &b) const;
	Vector Vector::operator-(const Vector &b) const;
	Vector Vector::operator*(float b) const;
	float dot(const Vector &b) const;
	Vector Vector::operator*(const Vector &b) const;
private:
	float x, y, z;

};


