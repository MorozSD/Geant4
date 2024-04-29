#include <omp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Point
{
private:
    int pointId, clusterId;
    int dimensions;
    vector<double> values;

    vector<double> lineToVec(string &line)
    {
        vector<double> values;
        string tmp = "";

        for (int i = 0; i < (int)line.length(); i++)
        {
            if ((48 <= int(line[i]) && int(line[i])  <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
            {
                tmp += line[i];
            }
            else if (tmp.length() > 0)
            {

                values.push_back(stod(tmp));
                tmp = "";
            }
        }
        if (tmp.length() > 0)
        {
            values.push_back(stod(tmp));
            tmp = "";
        }

        return values;
    }

public:
    Point(int id, string line)
    {
        pointId = id;
        values = lineToVec(line);
        dimensions = values.size();
        clusterId = 0; // Initially not assigned to any cluster
    }

    int getDimensions() { return dimensions; }

    int getCluster() { return clusterId; }

    int getID() { return pointId; }

    void setCluster(int val) { clusterId = val; }

    double getVal(int pos) { return values[pos]; }
};

class Cluster
{
private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster(int clusterId, Point centroid)
    {
        this->clusterId = clusterId;
        for (int i = 0; i < centroid.getDimensions(); i++)
        {
            this->centroid.push_back(centroid.getVal(i));
        }
        this->addPoint(centroid);
    }

    void addPoint(Point p)
    {
        p.setCluster(this->clusterId);
        points.push_back(p);
    }

    bool removePoint(int pointId)
    {
        int size = points.size();

        for (int i = 0; i < size; i++)
        {
            if (points[i].getID() == pointId)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    void removeAllPoints() { points.clear(); }

    int getId() { return clusterId; }

    Point getPoint(int pos) { return points[pos]; }

    int getSize() { return points.size(); }

    double getCentroidByPos(int pos) { return centroid[pos]; }

    void setCentroidByPos(int pos, double val) { this->centroid[pos] = val; }
};

class KMeans
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;
    string output_dir;

    void clearClusters()
    {
        for (int i = 0; i < K; i++)
        {
            clusters[i].removeAllPoints();
        }
    }

    int getNearestClusterId(Point point)
    {
        double sum = 0.0, min_dist;
        int NearestClusterId;
        if(dimensions==1) {
            min_dist = abs(clusters[0].getCentroidByPos(0) - point.getVal(0));
        }	
        else 
        {
          for (int i = 0; i < dimensions; i++)
          {
             sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
             // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
          }
          min_dist = sqrt(sum);
        }
        NearestClusterId = clusters[0].getId();

        for (int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;
            
            if(dimensions==1) {
                  dist = abs(clusters[i].getCentroidByPos(0) - point.getVal(0));
            } 
            else {
              for (int j = 0; j < dimensions; j++)
              {
                  sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
                  // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
              }

              dist = sqrt(sum);
              // dist = sum;
            }
            if (dist < min_dist)
            {
                min_dist = dist;
                NearestClusterId = clusters[i].getId();
            }
        }

        return NearestClusterId;
    }

public:
    KMeans(int K, int iterations, string output_dir)
    {
        this->K = K;
        this->iters = iterations;
        this->output_dir = output_dir;
    }

    void run(vector<Point> &all_points, vector<Cluster> &my_clusters)
    {
        total_points = all_points.size();
        dimensions = all_points[0].getDimensions();

        // Initializing Clusters
        vector<int> used_pointIds;

        for (int i = 1; i <= K; i++)
        {
            while (true)
            {
                int index = rand() % total_points;

                if (find(used_pointIds.begin(), used_pointIds.end(), index) ==
                    used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    all_points[index].setCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        cout << "Clusters initialized = " << clusters.size() << endl
             << endl;

        cout << "Running K-Means Clustering.." << endl;

        int iter = 1;
        while (true)
        {
            cout << "Iter - " << iter << "/" << iters << endl;
            bool done = true;

            // Add all points to their nearest cluster
            #pragma omp parallel for reduction(&&: done) num_threads(16)
            for (int i = 0; i < total_points; i++)
            {
                int currentClusterId = all_points[i].getCluster();
                int nearestClusterId = getNearestClusterId(all_points[i]);

                if (currentClusterId != nearestClusterId)
                {
                    all_points[i].setCluster(nearestClusterId);
                    done = false;
                }
            }

            // clear all existing clusters
            clearClusters();

            // reassign points to their new clusters
            for (int i = 0; i < total_points; i++)
            {
                // cluster index is ID-1
                clusters[all_points[i].getCluster() - 1].addPoint(all_points[i]);
            }

            // Recalculating the center of each cluster
            for (int i = 0; i < K; i++)
            {
                int ClusterSize = clusters[i].getSize();

                for (int j = 0; j < dimensions; j++)
                {
                    double sum = 0.0;
                    if (ClusterSize > 0)
                    {
                        #pragma omp parallel for reduction(+: sum) num_threads(16)
                        for (int p = 0; p < ClusterSize; p++)
                        {
                            sum += clusters[i].getPoint(p).getVal(j);
                        }
                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
                    }
                }
            }

            if (done || iter >= iters)
            {
                cout << "Clustering completed in iteration : " << iter << endl
                     << endl;
                break;
            }
            iter++;
        }

        ofstream pointsFile;
        pointsFile.open(output_dir + "/" + to_string(K) + "-points.txt", ios::out);

        for (int i = 0; i < total_points; i++)
        {
            pointsFile << all_points[i].getCluster() << endl;
        }

        pointsFile.close();

        // Write cluster centers to file
        ofstream outfile;
        outfile.open(output_dir + "/" + to_string(K) + "-clusters.txt");
        if (outfile.is_open())
        {
            for (int i = 0; i < K; i++)
            {
                cout << "Cluster " << clusters[i].getId() << " centroid : ";
                my_clusters.push_back(clusters[i]);
                for (int j = 0; j < dimensions; j++)
                {
                    cout << clusters[i].getCentroidByPos(j) << " ";    // Output to console
                    outfile << clusters[i].getCentroidByPos(j) << " "; // Output to file
                }
                cout << endl;
                outfile << endl;
            }
            outfile.close();
        }
        else
        {
            cout << "Error: Unable to write to clusters.txt";
        }
        
    }
};

int main(int argc, char **argv)
{
    // Need 3 arguments (except filename) to run, else exit
    if (argc != 4)
    {
        cout << "Error: command-line argument count mismatch. \n ./kmeans <INPUT> <K> <OUT-DIR>" << endl;
        return 1;
    }

    string output_dir = argv[3];

    // Fetching number of clusters
    int K = atoi(argv[2]);

    // Open file for fetching points
    string filename = argv[1];
    ifstream infile(filename.c_str());

    if (!infile.is_open())
    {
        cout << "Error: Failed to open file." << endl;
        return 1;
    }

    // Fetching points from file
    int pointId = 1;
    vector<Point> all_points;
    vector<Cluster> my_clusters;
    string line;

    while (getline(infile, line))
    {
        Point point(pointId, line);
        all_points.push_back(point);
        pointId++;
    }
    
    infile.close();
    cout << "\nData fetched successfully!" << endl
         << endl;

    // Return if number of clusters > number of points
    if ((int)all_points.size() < K)
    {
        cout << "Error: Number of clusters greater than number of points." << endl;
        return 1;
    }

    // Running K-Means Clustering
    int iters = 100;

    KMeans kmeans(K, iters, output_dir);
    kmeans.run(all_points,my_clusters);
    
    //seluette my attempt
    //all_points - vector of points
    double a = 0., A = 0., r, R = 0., X, Y, b = 0., s, sumS = 0., TotalS = 0.;
    int N = all_points.size(), I = 0;
    int ID_c_prev, ID_c_next;
   
        
    if(K>1){ 
    // Nearest clasters
    double** Near = new double *[K];
    	for (int i = 0; i < K; ++i)
	Near[i] = new double [K];
		for (int i = 0; i < K; i++)
		{
			X = my_clusters[i].getCentroidByPos(0);
			Y = my_clusters[i].getCentroidByPos(1);
			cout << " i " << i << X << " " << Y << endl; 
			Near[i][i] = 0;
			for (int j = I+1; j < K; j++)
			{
			        r = sqrt(pow((X-my_clusters[j].getCentroidByPos(0)),2)+pow((Y-my_clusters[j].getCentroidByPos(1)),2));
			        Near[i][j] = r;
			        Near[j][i] = r;
			        cout << " i " << i << " j " << j << " X " << X << " Y " << Y << " X1 " <<my_clusters[j].getCentroidByPos(0) << " Y1 " << my_clusters[j].getCentroidByPos(1) << " r " << r << endl; 
			}
		I++;	
		}
		
	/*	for (int i = 0; i < K; i++)
		{
			for (int j = 0; j < K; j++)
			{
			       cout << " i " << i << " j "<< j << " Near " << Near[i][j] << endl; 
			}
			
		}*/
		// ai
		cout << " Near is ready! " << endl;
		
		double Xi, Yi, Xj, Yj, Xk, Yk, amed, bmed;
		
		for(int i = 0; i<K; i++){
			//int i = 1;
			int j = 0;
			double  r_min=Near[i][0], r;
			if (!i) r_min = Near[i][1];
			cout << " i cluster " << i << endl;
			for (int jj = 1; jj<K; jj++){
				if(i!=jj){
					 r = Near[i][jj];
					// cout << " r " << r << " r_min " << r_min << endl;
					 if(r < r_min){ 
					 	r_min = r; 
					 	j = jj;
					 //	cout << " r_min " << r_min << " r " << r << " j " << j << endl;
					 }
				 }
			
			}
			cout << " Nearest cluster " << j << endl;
			int ClusterSizei = my_clusters[i].getSize(), ClusterSizej = my_clusters[j].getSize();;
			cout << " ClusterSizei " << ClusterSizei << " ClusterSizej " << ClusterSizej << endl;
			
			for (int k = 0; k<ClusterSizei; k++)
			{
				Xi = my_clusters[i].getPoint(k).getVal(0);
				Yi = my_clusters[i].getPoint(k).getVal(1);
				for (int p1 = 0; p1 < ClusterSizei; p1++) // can be optimized to not count dis between equel points
		                	{
		                	    Xj = my_clusters[i].getPoint(p1).getVal(0);
		                	    Yj = my_clusters[i].getPoint(p1).getVal(1);
		                	    a+=sqrt(pow((Xi-Xj),2)+pow((Yi-Yj),2));
		                //	    cout << " ai " << a <<" bi " << b << " xi " << Xi << " xj " << Xj << " p1 " << p1 << endl;
		                	}
		                for (int p2 = 0; p2 < ClusterSizej; p2++)
		                	{
		                	    Xj = my_clusters[j].getPoint(p2).getVal(0);
		                	    Yj = my_clusters[j].getPoint(p2).getVal(1);
		                	    b+=sqrt(pow((Xi-Xj),2)+pow((Yi-Yj),2));
		                //	    cout << " bi " << b << " xi " << Xi << " xj " << Xj<<" p2 " << p2 << endl;
		                	}
		                	///cout << " sum ai " << a << " sum bi " << b << endl;
		                	amed = a/ClusterSizei; bmed = b/ClusterSizej;
		                	cout << " amed " << amed << " bmed " << bmed << endl;
		                	if (amed>bmed) s = (bmed-amed)/amed;
		                	else s = (bmed-amed)/bmed;
		                	cout << " s " << s << endl;  
		                	sumS +=s;
		                	a = 0.;
		                	b = 0.;
				 
			}
			sumS = sumS/ClusterSizei;
			cout << " SumS " << sumS <<" ClusterSizei " << ClusterSizei << endl;
			TotalS += sumS;
			sumS = 0.;
			
		}
		TotalS = TotalS/K;
		cout << " TotalS " << TotalS << endl;
		/*for (int i = 0; i<N; i++){
    			ID_c_prev = all_points[i].getCluster();
    			for (int j = 1; j<N; j++){
    				cout << "point" << all_points[j].getID() << " id_c " << all_points[j].getCluster() << " x " << all_points[j].getVal(0) << " y " <<all_points[j].getVal(1) << endl;
    				ID_c_next = all_points[j].getCluster();
    		
    				if (ID_c_prev==ID_c_next){
    					r = pow((all_points[j-1].getVal(0)-all_points[j].getVal(0)),2)+pow((all_points[j-1].getVal(1)-all_points[j].getVal(1)),2);
    				a+= sqrt(r);

    				}
    		
    		
    			}	
    		cout << " ai " << a << endl;
    		a = 0;
    		}*/
		
		
		
		
for (int i = 0; i < K; i++)
delete[] Near[i];
delete [] Near;




     }	
     else{
     //K==1
     b = 0;
     }
     
	


   
    cout << " Length " << N << endl;
  /*  for (int i = 0; i < K; i++)
    {
        point.Get 
    }
    a =  */
    
    

    return 0;
}
