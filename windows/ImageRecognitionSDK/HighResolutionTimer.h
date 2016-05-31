#include <exception>
#include <windows.h>
#include <cassert>

using namespace std;

class HighResolutionTimer
{
public:
	HighResolutionTimer()
	{
		BOOL ret = QueryPerformanceFrequency(&_frequency);
		if (ret == FALSE)
		{
			// be careful here!!
			throw exception("Do not support QueryPerformanceFrequency!");
		}
	}

	bool Reset()
	{
		BOOL ret = QueryPerformanceCounter(&_begin);
		return (ret == TRUE);
	}

	// unit in ms
	double Duration()
	{
		BOOL ret = QueryPerformanceCounter(&_end);
		assert(ret == TRUE);
		if(ret != TRUE)
		{
			throw exception("Do not support QueryPerformanceCounter!");
		}

		return 1000.0 * (_end.QuadPart - _begin.QuadPart) / _frequency.QuadPart;
	}

private:
	LARGE_INTEGER _frequency;
	LARGE_INTEGER _begin;
	LARGE_INTEGER _end;
};

