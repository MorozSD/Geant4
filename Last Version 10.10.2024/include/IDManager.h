/*

Straw Tracker description

ID_Manager.h

Author: Aytadzh Allakhverdieva
Created on: September, 2023

@class IDManager


 */
 
#ifndef IDManager_h
#define IDManager_h 1

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <math.h>
#include <iterator>
#include <string>
#include <sqlite3.h>

// GeoModel includes
#include "GeoModelKernel/GeoIdentifierTag.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoAlignableTransform.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/Query.h"
#include "GeoModelDBManager/GMDBManager.h"


class IDManager 
{

private:

	int BarrelID = 0;
	int ForECID = 1;
	int BackECID = 2;
	
	int SubSystemID = 0;
	
	int ST_subsystemID = 12;
	
	std::string m_dbpath;

	
	
	
public:
	
	IDManager(std::string str_path);
	
	~IDManager();
	
	
	
	void printVolume_info(int VolumeID);
	
	void printST_info(int VolumeID);
	
	int getST_bar_layer_num(int VolumeID);
	
protected:
  class Imp;
  Imp * m_d;
	
	
	
};

#endif

class IDManager::Imp {
public:
    // constructor
  Imp (IDManager* dbm) : idmanager(dbm), db(nullptr)
  {}

  // The class
  IDManager* idmanager;

  // Pointer to SQLite connection
  sqlite3* db;


};

IDManager::IDManager(std::string str_path) : m_dbpath(str_path), m_d(new Imp(this))
{
  	const char* path = str_path.c_str();
	
	int opendb = sqlite3_open(path, &m_d->db);
	
	if (opendb == SQLITE_OK) 
	{
    		std::cout << "Database has been opened successfully!" << std::endl;
    	}
    	
    	else 
    	{
    		std::cout << "Error! Database can't open! Check your DB!" << std::endl;
    	}
  


}


IDManager::~IDManager()
{
  sqlite3_close(m_d->db);
  m_d->db = nullptr;
  delete m_d;
  m_d = nullptr;
}






int my_callback(void* counter_ptr,int,char**,char**);

void IDManager::printVolume_info(int VolumeID)
{	
	
	std::string sql1 = "SELECT * FROM IdentifierTags WHERE identifier = ";
	std::string sql2 = std::to_string(VolumeID);	
	std::string sql3 = ";";
	
	std::string str_sql = sql1 + sql2 + sql3;

	const char* sql = str_sql.c_str();
	
	
	int counter = 0;
	
	int statement = sqlite3_exec(m_d->db, sql, my_callback, &counter, NULL);
	
	
	if (counter != SQLITE_OK )
	{
		std::cout << "ID is founded!"<< std::endl;
		std::cout << "Printing info from: "<< VolumeID  << "."<< std::endl;
		
		SubSystemID = VolumeID & 0xf;
		
		if (SubSystemID == ST_subsystemID)
		{
			printST_info(VolumeID);
		}
		
			
	}
	
	else
	{
	
		std::string err1 = "There is no information about volume with ID = ";
		std::string error =  err1 + sql2;
		
		std::string error_msg =  err1 + sql2;
		
		std::cout << error_msg << std::endl;
		std::cout << "Check your db!" << std::endl;
	}
	
	
	
}


void IDManager::printST_info(int VolumeID)
{
	std::string str_StrawTracker = "straw tracker";
	
	int barrel_endcap = 0;
	
	barrel_endcap = (VolumeID>>4)&0x3;
	
	if( barrel_endcap == BarrelID )
	{
		std::string str_barrel = "barrel";
		
		int octant = (VolumeID>>6)&0x7;
			
		int layer = (VolumeID>>9)&0x3f;
		
		int angle = (VolumeID>>15)&0x3;
			
		int stangle = 0;
			
		if (angle == 0)
		{
			stangle = 0;
		}
			
		if (angle == 1)
		{
		
			stangle = + 5;
		}
			
		if (angle == 2)
		{
		
			stangle = - 5;
		}
		
		int tubenumber = (VolumeID>>17)&0x7f;
		
	
		std::cout << "Volume with ID: " << VolumeID << " is on " << str_StrawTracker << " in " << str_barrel  << " in octant # "<< octant <<  " in layer # "<< layer << " with angle = " << stangle << " with tube number = " << tubenumber << std::endl;
	}
	
	if( barrel_endcap == ForECID )
	{
		std::string str_forendcap = "forward end-cap"; 	
	}
	
	if( barrel_endcap == BackECID )
	{
		std::string str_backendcap = "bacward end-cap";
	}
		

}


int IDManager::getST_bar_layer_num(int VolumeID)
{

	std::string sql1 = "SELECT * FROM IdentifierTags WHERE identifier = ";
	std::string sql2 = std::to_string(VolumeID);	
	std::string sql3 = ";";
	
	std::string str_sql = sql1 + sql2 + sql3;

	const char* sql = str_sql.c_str();
	
	
	int counter = 0;
	
	int statement = sqlite3_exec(m_d->db, sql, my_callback, &counter, NULL);
	
	
	if (counter != SQLITE_OK )
	{
		std::cout << "ID is founded!"<< std::endl;
		std::cout << "Printing info from: "<< VolumeID  << "."<< std::endl;
		
		int layer = (VolumeID>>9)&0x3f;
		
		return layer;
		
			
	}
	
	else
	{
	
		std::string err1 = "There is no information about volume with ID = ";
		std::string error =  err1 + sql2;
		
		std::string error_msg =  err1 + sql2;
		
		std::cout << error_msg << std::endl;
		std::cout << "Check your db!" << std::endl;
	}

return 0;	
}








int my_callback(void* counter_ptr,int,char**,char**)
{

	(*((int*)counter_ptr))++;
	//std::cout << "ID is find!" << std::endl;
		
	return SQLITE_OK;
	
}
