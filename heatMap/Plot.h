#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std;
class Plot
{
	public:
		Plot();                                                 // constructor
		Plot(int object_number);                                // constructor
		Plot(double data_input[], int size_input, int obj_num); // constructor
		~Plot();                                                // deconstructor
		
		void SetData(double data_input[]); // set data
		void SetObjectNumber(int obj_num); // set object number
		double GetMax();                   // get Maximum value
		double GetMin();                   // get Minimum value
		double GetRange();                 // get Range
		double GetFirstVal();              // get the first value of the data
		double GetMean();                  // get Mean
		double GetStdDev();                // get Standard deviation
		double *GetData(size_t s);         // get Data
		void printData();                  // print data
		void dealloc();                    // deallocate arrays

	private:
		int size;                          // size of the data array
		double* data;                      // pointer to the data array
		int object_num;                    // object number
};
