#ifndef ST_ID_h
#define ST_ID_h 1

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <math.h>
#include <iterator>



class ST_ID 
{

private:
	
	int id = 0;
	
	
public:

	ST_ID()
	{};
	
	~ST_ID()
	{};
		
	

	int createST_BarID(int LayerNum, int STNumber)
{	
	
	int ST_Number = 0;
	
	int ST_LayerNum = 0;
	
	
	ST_Number = STNumber;
		
	ST_LayerNum = LayerNum<<8;
	
	//creation ID
	
	id = ST_Number + ST_LayerNum;
	
	
return id;
}

int get_LayerNum(int id)
{
	int ST_LayerNum = 0;
	
	ST_LayerNum = (id>>8)&0x7;
	
	return ST_LayerNum;
}
	
int get_TubeNum(int id)
{
	int ST_Number = 0;
	
	ST_Number = id&0xff;
	
	return ST_Number;
}
	

	
};

#endif
