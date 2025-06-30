#include <iostream>
#include <fstream>
#include <string>

int main()
{
	
	std::string line;
	double BAD_Tracks = 0, ALL_FIT_Tracks = 0, ALL_Tracks = 0, RSS_z_Tracks = 0, RSS_Tracks = 0, RSS_Tracks_Bad = 0, RSS_z_Tracks_Bad = 0, All_Fit_Vert = 0;
	double Wrong_vert = 0, Wrong = 0, All_Vert = 0;
	double x = 0, y = 0, z = 0, k = 0, n = 0, m = 0, p = 0;
	
	//std::ifstream in("Last Version 29.04.25/build/Output.txt");
	std::ifstream in("LV_19.12/build/Clust.txt");
	
	if(in.is_open())
	{
		//while(std::getline(in, line))
		while(!in.eof())
		{
			//std::cout << line << std::endl;
			in >> x;
			in >> y;
			in >> z;
			in >> k;
			in >> n;
			in >> m;
			in >> p;
			
			//for Clust
			
			/*in >> x;
			in >> y;
			in >> z;*/
			
			//ALL_Tracks +=x; Wrong += y; Wrong_vert += z;
			//BAD_Tracks +=x; ALL_FIT_Tracks += m; ALL_Tracks += p; RSS_z_Tracks +=k; RSS_Tracks +=y; RSS_Tracks_Bad +=z; RSS_z_Tracks_Bad +=n;
			ALL_Tracks +=x; ALL_FIT_Tracks += y; RSS_z_Tracks +=z; RSS_z_Tracks_Bad +=k; All_Fit_Vert+=n; Wrong += m; Wrong_vert += p;
			//std::cout << " After RSS " << RSS << " BAD " << BAD << " All_Fit " << All_Fit << " All " << All << std::endl;
			
		
			
		}
		ALL_Tracks -=x; ALL_FIT_Tracks -= y; RSS_z_Tracks -=z; RSS_z_Tracks_Bad -=k; All_Fit_Vert-=n; Wrong -= m; Wrong_vert -= p;
		std::cout << " All_Tracks " << ALL_Tracks << " Dphi " << ALL_FIT_Tracks/ALL_Tracks*100 << " RSS_Z " << (ALL_FIT_Tracks-RSS_z_Tracks)/ALL_Tracks*100 << " BAD " << RSS_z_Tracks_Bad/ALL_Tracks*100 << " ALL Fit Vert " << All_Fit_Vert << " Wrong " << Wrong << " Vert_wrong " << Wrong_vert <<  " Eff " << 100 - (Wrong+Wrong_vert)/All_Fit_Vert*100<< std::endl;
		
		//BAD_Tracks -=x; ALL_FIT_Tracks -= m; ALL_Tracks -= p; RSS_z_Tracks -=k; RSS_Tracks -=y; RSS_Tracks_Bad -=z; RSS_z_Tracks_Bad -=n;
		//std::cout << " BAD " << BAD_Tracks<<  " RSS " << RSS_Tracks << " RSS_BAD " << RSS_Tracks_Bad << " RSS_z " << RSS_z_Tracks << " RSS_z_BAD " << RSS_z_Tracks_Bad  <<" ALL_FIT " << ALL_FIT_Tracks << " ALL " << ALL_Tracks << std::endl; 
	//	ALL_Tracks -=x; Wrong -= y; Wrong_vert -= z;
		//std::cout << " ALL " << ALL_Tracks << " Wrong " << Wrong << " Vert_wrong " << Wrong_vert << std::endl; 
		//std::cout << " ALL " << ALL_Tracks << " ALL_GOOD " << Wrong << " GOOD_But_>3*mm " << Wrong_vert << std::endl; 
	}
	// BAD 78 RSS 319 RSS_BAD 15 RSS_z 19 RSS_z_BAD 14 ALL_FIT 8677 ALL 10858

	in.close();
	
	/*std::ifstream Vert("LV_19.12/build/Vertex.txt");
	
	if(Vert.is_open())
	{
		//while(std::getline(in, line))
		while(!Vert.eof())
		{
			//std::cout << line << std::endl;
			//Vert >> p;
			
			
			All_Vert++;			//std::cout << " After RSS " << RSS << " BAD " << BAD << " All_Fit " << All_Fit << " All " << All << std::endl;
			
		
			
		}
		
		std::cout << " All_Vert " << All_Vert << std::endl;
		
		
	}
	// BAD 78 RSS 319 RSS_BAD 15 RSS_z 19 RSS_z_BAD 14 ALL_FIT 8677 ALL 10858

	Vert.close();*/
	
	
	return 0;
}
