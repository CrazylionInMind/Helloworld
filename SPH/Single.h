#pragma once
template<class T>
class Single
{
public:
	static T* getinstance()
	{
		if (_single == NULL)
			_single = new T;
		return _single;
	}
	~Single() = default;
	Single() = default;

private:
	static T* _single;
};
template<class T>
T* Single<T>::_single = NULL;