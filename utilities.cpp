#include "utilities.h"
#include <math.h>

float std_dev(float arr[], float mean, int start, int end)
{
	if (start >= end)
	{
		return 0;
	}

	float _stddev = 0.0;

	for (int i = start; i <= end; i++)
	{
		_stddev += (arr[i] - mean)*(arr[i] - mean);
	}
	_stddev /= (end - start + 1);

	return sqrt(_stddev);
};
