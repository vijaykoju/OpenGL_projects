#include "Plot.h"

//--------------------------constructor----------------------------------//
// input parameters : none 
// set default size to 20
Plot::Plot()
{
	size = 20;
}

//--------------------------constructor----------------------------------//
// input parameters : int obj_num         --> object number
// set default size to 20
Plot::Plot(int obj_number)
{
	size = 20;
	object_num = obj_number;
}	
//--------------------------constructor----------------------------------//
// input parameters : double data_input[] --> array of the whole data
//                  : int size_input      --> size of output data array
//                  : int obj_num         --> object number
Plot::Plot(double data_input[],int size_input, int obj_num)
{
	size = size_input; // set size 
	object_num = obj_num; // set object_num
	data = new double[size]; // dynamically allocate the size of data[]
	for (int i=0; i<size; ++i)
		data[i] = data_input[(size*object_num)+i]; // set the values of data
}


//--------------------------destructor-----------------------------------//
Plot::~Plot()
{
}

//-----------------------------------------------------------------------//
// Method name     : SetData()           --> set data
// input parameter : double data_input[] --> full data array 
void Plot::SetData(double data_input[])
{
	data = new double[size];
	for (int i=0; i<size; ++i)
		data[i] = data_input[(size*object_num)+i];
}

//-----------------------------------------------------------------------//
// Method name     : SetObjectNumber      --> set object number
// input parameter : int obj_num          --> object number
void Plot::SetObjectNumber(int obj_num)
{
	object_num = obj_num;
}

//-----------------------------------------------------------------------//
// Method name : GetMax() --> return the maximum value of the data array
double Plot::GetMax()
{
	double max_val = data[0]; // set the first element to be the maximum
	for (int i=1; i<size; ++i)
	{
		if (data[i]>max_val) 
			max_val = data[i]; // update the max_val if the next element is bigger
	}
	return max_val;
}

//-----------------------------------------------------------------------//
// Method name : GetMax() --> return the maximum value of the data array
double Plot::GetMin()
{
	double min_val = data[0]; // set the first element to be the minimum
	for (int i=1; i<size; ++i)
	{
		if (data[i]<min_val)
			min_val = data[i]; // update the min_val if the next elemenet is small
	}
	return min_val;
}

//-----------------------------------------------------------------------//
// Method name : GetRange() --> return the range of the data
double Plot::GetRange()
{
	double max_val = GetMax(); // get the maximum value
	double min_val = GetMin(); // get the minimum valu
	return max_val-min_val;    // return the difference
}
//-----------------------------------------------------------------------//
// Method name : GetFirstVal() --> return the first value of the object
double Plot::GetFirstVal()
{
	return data[0]; // return the first element of the array
}


//-----------------------------------------------------------------------//
// Method name : GetMean() --> compute and return the mean of data in the 
//               object
double Plot::GetMean()
{
	double totSum = 0; // intialize the total sum to zero
	for (int i=0; i<size; ++i) 
		totSum += data[i]; // summation of all the data
	return totSum/size; // mean
}

//-----------------------------------------------------------------------//
// Method name : GetStdDev() --> compute and return the unbiased standard
//              deviation of data in the object
double Plot::GetStdDev()
{
	double sqDisFrMean = 0; 
	double sMean = GetMean(); // compute mean
	for (int i=0; i<size; ++i)
		sqDisFrMean += (data[i]-sMean)*(data[i]-sMean); // deviation from mean
	return sqrt(sqDisFrMean/(size-1)); // standard deviation
}

//-----------------------------------------------------------------------//
// Method name : GetData()  --> get data array associated with the object
// input param : size_t s   --> size of the data
// output      : pointer to the data array
double* Plot::GetData(size_t s)
{
	double *dat = new double[s];
	for (size_t i=0; i<s; ++i)
		dat[i] = data[i];
	return dat;
}

//-----------------------------------------------------------------------//
// Method name : printData() --> print out the data to standard output
void Plot::printData()
{
	cout << "The data is: " << size << endl;
	for (int i=0; i<size; ++i)
		cout << data[i] << " ";
	cout << endl;
}

//-----------------------------------------------------------------------//
// Method name : dealloc() --> free dynamically allocated arrays
void Plot::dealloc()
{
	delete[] data;
}

