#include "MyfileRW.h"
std::string MyFileRW::_WriteFileName = "";

MyFileRW::MyFileRW()
{
	_ReadFilename.clear();
	_WriteFileName.clear();
	in = NULL;
	out = NULL;
}

MyFileRW::~MyFileRW()
{
}

void MyFileRW::SetReadNmae(const std::string rname)
{
	if (in)
	{
		if (_ReadFilename != rname)
		{
			(*in).close();
			delete in;
			_ReadFilename = rname;
			in = new std::ifstream(_ReadFilename);
		}
	}
	else
	{
		_ReadFilename = rname;
		in = new std::ifstream(_ReadFilename);
	}
}

void MyFileRW::SetWriteFileName(const std::string wname)
{
	_WriteFileName = wname;
	if (_WriteFileName == "")
	{
		out->close();
		delete out;
	}
	_WriteFileName = wname;
	out = new std::ofstream(_WriteFileName);
}

bool MyFileRW::ReadFile(Fluid * input)
{
	if (!input)
		return false;
	if (!in || (*in).is_open() == NULL)
		return false;
	
	while (!(*in).eof())
	{
		float pos[3];
		for (int i=0;i<3;i++)
		{
			(*in) >> pos[i];
		}//读取位置

		float vel[3];
		for (int i = 0;i < 3;i++)
		{
			(*in) >> vel[i];
		}//读取速度

		float mass;
		(*in) >> mass;		
	}
	return true;
}

void MyFileRW::CloseReadFile()
{
	if (in)
	{
		(*in).close();
		delete in;
		_ReadFilename = "";
	}
}

void MyFileRW::CloseWriteFile()
{
	out->close();
	delete out;
	_WriteFileName = "";
}

void MyFileRW::WriteFloat(float input, int isChangeLine /*= 0*/, std::string filename /*= _WriteFileName*/)
{
	if (_WriteFileName == "" && filename == "")
	{
		cout << "没有为文件名赋值" << endl;
		return;
	}
	else if (_WriteFileName != filename)
	{
		if (_WriteFileName == "")
		{
			out->close();
			delete out;
		}
		_WriteFileName = filename;
		out = new std::ofstream(_WriteFileName);
	}

	if (_WriteFileName == filename && out->is_open())
	{
		out->setf(ios::fixed, ios::floatfield);
		out->precision(8);
		if (isChangeLine == 0)
		{
			*out << input;
		}
		else if (isChangeLine == 1)
		{
			*out << input << " ";
		}
		else if (isChangeLine == 2)
		{
			*out << input << ".\n";
		}
		return;
	}
	return;
}

