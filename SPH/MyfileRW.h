#pragma once
#include <string>
#include <iosfwd>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "Single.h"
#include "Fluid.h"
class MyFileRW :public Single<MyFileRW>
{
public:
	MyFileRW();
	~MyFileRW();
	void SetReadNmae(const std::string rname);
	void SetWriteFileName(const std::string wname);
	virtual bool ReadFile(Fluid *input);
	void CloseReadFile();
	void CloseWriteFile();
	template<typename M>
	bool WriteSomething(M input, int isChangeLine = 0, std::string filename = _WriteFileName);
	void WriteFloat(float input, int isChangeLine = 0, std::string filename = _WriteFileName);
private:
	std::string _ReadFilename;
	static std::string _WriteFileName;
	std::ofstream *out;
	std::ifstream *in;
};

template<typename M>
bool MyFileRW::WriteSomething(M input, int isChangeLine /*= 0*/, std::string filename /*= _WriteFileName*/)
{
	if (_WriteFileName == "" && filename == "")
	{
		cout << "没有为文件名赋值" << endl;
	}
	else if (_WriteFileName != filename)
	{
		if (_WriteFileName=="")
		{
			out->close();
			delete out;
		}
		_WriteFileName = filename;
		out = new std::ofstream(_WriteFileName);
	}

	if (_WriteFileName == filename && out->is_open())
	{
		if (isChangeLine ==0)
		{
			*out << input ;
		}
		else if (isChangeLine ==1)
		{
			*out << input<<" ";
		}
		else if(isChangeLine == 2)
		{
			*out << input <<"\n";
		}
		return true;
	}
	return false;
}

