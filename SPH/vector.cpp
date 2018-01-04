#include"vector.h"

Vector Vector::V_normal()//¹éÒ»»¯
{
	Vector temp1;
	float temp;
	temp = this->x*this->x + this->y*this->y + this->z*this->z;
	temp = (float)sqrt(temp);
	temp1.x = this->x / temp;
	temp1.y = this->y / temp;
	temp1.z = this->z / temp;
	return temp1;
}
Vector::Vector()
{
	this->x = 0.0f;
	this->y = 0.0f;
	this->z = 0.0f;
}

Vector::Vector(float a, float b, float c)
{
	this->x = a;
	this->y = b;
	this->z = c;
}

Vector Vector::operator+(const Vector &b) const
{
	Vector sum;
	sum.x = x + b.x;
	sum.y = y + b.y;
	sum.z = z + b.z;
	return sum;
}

Vector Vector::operator-(const Vector &b) const
{
	Vector sum;
	sum.x = x - b.x;
	sum.y = y - b.y;
	sum.z = z - b.z;
	return sum;
}

Vector Vector::operator*(float b) const
{
	Vector temp;
	temp.x = b*x;
	temp.y = b*y;
	temp.z = b*z;
	return temp;
}

float Vector::dot(const Vector & b) const
{
	float temp;
	temp = this->x*b.x + this->y*b.y + this->z*b.z;
	return temp;
}

Vector Vector::operator*(const Vector & b) const
{
	Vector temp;
	temp.x = y*b.z - z*b.y;
	temp.y = z*b.x - x*b.z;
	temp.z = x*b.y - y*b.x;
	return temp;
}



Vector::~Vector()
{
}

