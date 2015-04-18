#include<stdio.h>
#include<iostream>
#include<string>
#include<vector>
#include<time.h>
#include<stdlib.h>
#include<sstream>
#include<random>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<ios>
#include<cmath>
#include<array>
#include <iterator> 
#include <cstdlib>
#include <time.h>
#include <new>

using namespace std;

double roundoff(double n) {
	return (floor(n*100))/100;
}

int signums(int k) {
	if (k>0)
		return 1;
	else if (k<0)
		return -1;
	else
		return 0;
}

double magnitude(vector<int> a) {
	double k=sqrt(a[0]*a[0] + a[1]*a[1]);
	return k;
}

double Distance(int x1, int y1, int x2, int y2){
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

int simulation_count;
int simulation_type;
int simulation_time;
int mode;
int bird_types;

//ENVT_SIZE
int grid_height;
//ENVT_SIZE
int grid_width;

vector<int> bird_num;
vector<int> max_vel;
int bird_size;

int detectionRange;
int separationRange;

double SEP_CONST;
double ALIGN_CONST;
double COH_CONST;

int iteration_size;
int groups;
int DOIS[3];
int history_size;

double avg_Lemergent_dist;
double avg_Lemergent_ic;
double avg_Lemergent_otf;

double std_dev_Lemergent_dist;
double std_dev_Lemergent_ic;
double std_dev_Lemergent_otf;

vector<double> sum_arr1;
vector<double> sum_arr2;
vector<double> sum_arr3;

vector<double> avg_dist;
vector<double> avg_ic;
vector<double> avg_otf;

vector<double> std_dev_dist;
vector<double> std_dev_ic;
vector<double> std_dev_otf;

vector<double> run_data_Lemergent_dist;
vector<double> run_data_Lemergent_ic;
vector<double> run_data_Lemergent_otf; 

vector<int> sims;
vector<vector<int>> run_data_dist;
vector<vector<int>> run_data_ic;
vector<vector<int>> run_data_otf;

//---
int Lwhole_states;
int Lzero_dist;
int Lzero_ic;
int Lzero_otf;

int min_doi_ic;
int max_doi_ic;
double min_doi_dist;
double max_doi_dist;

vector<int> dist;
vector<int> ic;
vector<int> otf;

//ENVT
vector<int> width;
vector<vector<int>> grid;
vector<int> history;
//---

void resize_mainrun_vectors() {

	avg_Lemergent_dist=0;
	avg_Lemergent_ic=0;
	avg_Lemergent_otf=0;

	std_dev_Lemergent_dist=0;
	std_dev_Lemergent_ic=0;
	std_dev_Lemergent_otf=0;

	sum_arr1.resize(groups);
	sum_arr2.resize(groups);
	sum_arr3.resize(groups);

	avg_dist.resize(groups);
	avg_ic.resize(groups);
	avg_otf.resize(groups);

	std_dev_dist.resize(groups);
	std_dev_ic.resize(groups);
	std_dev_otf.resize(groups);

	run_data_Lemergent_dist.resize(simulation_count);
	run_data_Lemergent_ic.resize(simulation_count);
	run_data_Lemergent_otf.resize(simulation_count);

	sims.resize(simulation_count);
	run_data_dist.resize(groups,sims);
	run_data_ic.resize(groups,sims);
	run_data_otf.resize(groups,sims);

	dist.resize(groups);
	ic.resize(groups);
	otf.resize(groups);
	
	width.resize(grid_width);
	grid.resize(grid_height,width);

}

void set_bird_size() {
	for (int i=0;i<bird_types;i++)
		bird_size += bird_num[i];
}

void initialize_DOI_stats() {

	for (int i=0;i<groups;i++) {
		avg_dist[i]=0;
		avg_ic[i]=0;
		avg_otf[i]=0;
		std_dev_dist[i]=0;
		std_dev_ic[i]=0;
		std_dev_otf[i]=0;
		sum_arr1[i]=0;
		sum_arr2[i]=0;
		sum_arr3[i]=0;
	}
}

void read_input_parameter_file() {

	int k;
	ifstream input_parameter_file;
	input_parameter_file.open("input_parameters.txt", ios::in);

	if (input_parameter_file.is_open()) {
		while (!input_parameter_file.eof()) {

			input_parameter_file >> simulation_count;
			input_parameter_file >> simulation_type;
			input_parameter_file >> mode;
			input_parameter_file >> bird_types;
			for (int i=0;i<bird_types;i++) {
				input_parameter_file >> k;
				bird_num.push_back(k);
			}
			input_parameter_file >> grid_height;
			input_parameter_file >> grid_width;
			input_parameter_file >> detectionRange;
			input_parameter_file >> separationRange;
			input_parameter_file >> SEP_CONST;
			input_parameter_file >> ALIGN_CONST;
			input_parameter_file >> COH_CONST;
			for (int i=0;i<bird_types;i++) {
				input_parameter_file >> k;
				max_vel.push_back(k);
			}
			if (simulation_type == 1)
				input_parameter_file >> iteration_size;
			else if (simulation_type == 2)
				input_parameter_file >> simulation_time;
			else
				cout << "ERROR - Wrong Simulation Type\n";
			input_parameter_file >> groups;
			input_parameter_file >> DOIS[0];
			input_parameter_file >> DOIS[1];
			input_parameter_file >> DOIS[2];
			input_parameter_file >> history_size;
		}
		input_parameter_file.close();
		set_bird_size();
	}
	else
		cout << "Cannot open Input Paramter File \n";
}

void generate_final_statistics_file(clock_t end) {

	double sum1=0, sum2=0, sum3=0;

	initialize_DOI_stats();

	for (int i=0;i<simulation_count;i++) {
		avg_Lemergent_dist += run_data_Lemergent_dist[i];
		avg_Lemergent_ic += run_data_Lemergent_ic[i];
		avg_Lemergent_otf += run_data_Lemergent_otf[i];

		for (int j=0;j<groups;j++) {
			avg_dist[j] += run_data_dist[j][i];
			avg_ic[j] += run_data_ic[j][i];
			avg_otf[j] += run_data_otf[j][i];
		}
	}

	avg_Lemergent_dist /= simulation_count;
	avg_Lemergent_ic /= simulation_count;
	avg_Lemergent_otf /= simulation_count;
	for (int j=0;j<groups;j++) {
		avg_dist[j] /= simulation_count;
		avg_ic[j] /= simulation_count;
		avg_otf[j] /= simulation_count;
	}

	for (int i=0;i<simulation_count;i++) {
		sum1 += pow(avg_Lemergent_dist - run_data_Lemergent_dist[i],2);
		sum2 += pow(avg_Lemergent_ic - run_data_Lemergent_ic[i],2);
		sum3 += pow(avg_Lemergent_otf - run_data_Lemergent_otf[i],2);
	}
	std_dev_Lemergent_dist = sqrt(sum1/simulation_count);
	std_dev_Lemergent_ic = sqrt(sum2/simulation_count);
	std_dev_Lemergent_otf = sqrt(sum3/simulation_count);

	for (int j=0;j<groups;j++) {
		for (int i=0;i<simulation_count;i++) {
			sum_arr1[j] += pow(avg_dist[j] - run_data_dist[j][i],2);
			sum_arr2[j] += pow(avg_ic[j] - run_data_ic[j][i],2);
			sum_arr3[j] += pow(avg_otf[j] - run_data_otf[j][i],2);
		}
		std_dev_dist[j] = sqrt(sum_arr1[j]/simulation_count);
		std_dev_ic[j] = sqrt(sum_arr2[j]/simulation_count);
		std_dev_otf[j] = sqrt(sum_arr3[j]/simulation_count);
	}

	ofstream output_file;
	output_file.open("final_statistics.txt", ios::out);

	if (DOIS[0]==1) {
		output_file << "\nAverage DOI Distribution with Min-Distance Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << avg_dist[i] << endl;
		output_file << "Average no of states in L-emergent (Min Distance Method) : " << avg_Lemergent_dist << endl;

		output_file << "\nStandard Deviation of DOI Distribution with Min-Distance Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << std_dev_dist[i] << endl;
		output_file << "Standard Deviation of no of states in L-emergent (Min Distance Method) : " << std_dev_Lemergent_dist << endl;
	}
	if (DOIS[1]==1) {
		output_file << "\nAverage DOI Distribution with Interaction Count Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << avg_ic[i] << endl;
		output_file << "Average no of states in L-emergent (Interaction Count Method) : " << avg_Lemergent_ic << endl;

		output_file << "\nStandard Deviation of DOI Distribution with Interaction Count Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << std_dev_ic[i] << endl;
		output_file << "Standard Deviation of no of states in L-emergent (Interaction Count Method) : " << std_dev_Lemergent_ic << endl;
	}
	if (DOIS[2]==1) {
		output_file << "\nAverage DOI Distribution with on-the-fly Interaction Count Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << avg_otf[i] << endl;
		output_file << "Average no of states in L-emergent (on-the-fly Interaction Count Method) : " << avg_Lemergent_otf << endl;

		output_file << "\nStandard Deviation of DOI Distribution with on-the-fly Interaction Count Method ... " << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/groups << " - " << (double)(i+1)/groups << ") : " << std_dev_otf[i] << endl;
		output_file << "Standard Deviation of no of states in L-emergent (on-the-fly Interaction Count Method) : " << std_dev_Lemergent_otf << endl;
	}
	output_file << "\nTime Taken for completion: " << (double)end/(double)CLOCKS_PER_SEC << endl;
}

class Bird{
public:

	//AGENT_TYPE
	int type;
	//AGENT_INSTANCE
	int No;
	//MOBILE_ATTRIBUTE
	int Velocity[2];
	//MOBILE_ATTRIBUTE
	int Coordinates[2];

	int maxVel;
	vector<int> separationVec;
	vector<int> cohesionVec;
	vector<int> alignmentVec;
	vector<int> MainVector;

	void resize_vectors() {
		separationVec.resize(2);
		cohesionVec.resize(2);
		alignmentVec.resize(2);
		MainVector.resize(2);
	}

	void update_attributes(int n) {

		Velocity[0] = MainVector[0];
		Velocity[1] = MainVector[1];

		int newx = Coordinates[0] + Velocity[0];
		int newy = Coordinates[1] + Velocity[1];

		if(newx > grid_width-1)
			newx = newx - grid_width;					
		else if(newx < 0)
			newx = newx + grid_width;

		if(newy > grid_height-1)
			newy = newy - grid_height;					
		else if(newy < 0)
			newy = newy + grid_height;

		Coordinates[0] = newx;
		Coordinates[1] = newy;

	}

	double BirdDistance(Bird b){
		double k=sqrt((b.Coordinates[0]-Coordinates[0])*(b.Coordinates[0]-Coordinates[0]) + (b.Coordinates[1]-Coordinates[1])*(b.Coordinates[1]-Coordinates[1]));
		return k;
	}

	vector<int> closestVect(Bird b){

		int ax = Coordinates[0], ay = Coordinates[1];
		int bx = b.Coordinates[0], by = b.Coordinates[1];
		vector<int> c(2);

		if (bx>ax) {
			if (abs(bx-ax) <= abs(bx-ax-grid_width))
				c[0] = bx-ax;
			else
				c[0] = bx-ax-grid_width;
		}
		else {
			if (abs(bx-ax) <= abs(bx-ax+grid_width))
				c[0] = bx-ax;
			else
				c[0] = bx-ax+grid_width;
		}
		if (by>ay) {
			if (abs(by-ay) <= abs(by-ay-grid_height))
				c[1] = by-ay;
			else
				c[1] = by-ay-grid_height;
		}
		else {
			if (abs(by-ay) <= abs(by-ay+grid_height))
				c[1] = by-ay;
			else
				c[1] = by-ay+grid_height;
		}	

		return c;
	}

	bool equate_boids(Bird b) {
		if (type == b.type && No == b.No) {
			if (Velocity[0] == b.Velocity[0] && Velocity[1] == b.Velocity[1] && Coordinates[0] == b.Coordinates[0] && Coordinates[1] == b.Coordinates[1])
				return true;
			else
				return false;
		}
		else
			return false;
	}

};

//----------------------------------
typedef vector<Bird> State;
State Flock;

State initial_state;
vector<State> states;
//----------------------------------

//MOBILE_RULE
vector<int> separation(int n,Bird& b){

	vector<int> separationVector(2);
	int x=0,y=0;
	int count=0;

	for(int i=0; i<Flock.size(); i++) {
		Bird& bird2 = Flock[i];
		if (b.type == bird2.type) {
			if (n==0) {
				if (b.BirdDistance(bird2) <= separationRange) {
					x = -b.closestVect(bird2)[0];
					y = -b.closestVect(bird2)[1];
					separationVector[0] += x;
					separationVector[1] += y;
					//count++;
				}
			}
			else {
				if (magnitude(b.closestVect(bird2)) <= separationRange) {
					x = -b.closestVect(bird2)[0];
					y = -b.closestVect(bird2)[1];
					separationVector[0] += x;
					separationVector[1] += y;
					//count++;
				}
			}
		}
		else {
			if (n==0) {
				if (b.BirdDistance(bird2) <= separationRange) {
					x = -b.closestVect(bird2)[0];
					y = -b.closestVect(bird2)[1];
					separationVector[0] += x;
					separationVector[1] += y;
					//count++;
				}
			}
			else {
				if (magnitude(b.closestVect(bird2)) <= separationRange) {
					x = -b.closestVect(bird2)[0];
					y = -b.closestVect(bird2)[1];
					separationVector[0] += x;
					separationVector[1] += y;
					//count++;
				}
			}
		}
	}

	separationVector[0] /= SEP_CONST;
	separationVector[1] /= SEP_CONST;

	return separationVector;
}

//MOBILE_RULE
vector<int> alignment(int n,Bird& b){

	vector<int> alignmentVector(2);
	int number = 0;

	for(int i=0; i<Flock.size(); i++){		                         
		Bird& bird2 = Flock[i];
		/*if (n==0) {
		if (BirdDistance(bird2) <= r.detectionRange){

		if(bird2.type == type && bird2.No != No){

		number++;
		alignmentVector[0] += bird2.Velocity[0];
		alignmentVector[1] += bird2.Velocity[1];
		}
		}
		}
		else {*/
		if (magnitude(b.closestVect(bird2)) <= detectionRange){

			if(bird2.type == b.type && bird2.No != b.No){

				number++;
				alignmentVector[0] += bird2.Velocity[0];
				alignmentVector[1] += bird2.Velocity[1];
			}
		}		
		/*}*/
	}

	if (number>0) {
		alignmentVector[0] = (alignmentVector[0]/number - b.Velocity[0])/ALIGN_CONST;
		alignmentVector[1] = (alignmentVector[1]/number - b.Velocity[1])/ALIGN_CONST;
	}

	return alignmentVector;
}

//MOBILE_RULE
vector<int> cohesion(int n,Bird& b){

	vector<int> cohesionVector(2);
	int number = 0; 

	for(int i=0; i<Flock.size(); i++){
		Bird& bird2 = Flock[i];
		/*if (n==0) {
		if (BirdDistance(bird2) <= r.detectionRange){
		if(bird2.type == type && bird2.No != No){
		number++;
		cohesionVector[0] += closestVect(bird2)[0];
		cohesionVector[1] += closestVect(bird2)[1];
		}
		}
		}
		else {*/
		if (magnitude(b.closestVect(bird2)) <= detectionRange){
			if(bird2.type == b.type && bird2.No != b.No){
				number++;
				cohesionVector[0] += b.closestVect(bird2)[0]/* + Coordinates[0]*/;
				cohesionVector[1] += b.closestVect(bird2)[1]/* + Coordinates[1]*/;
			}
		}		
		/*}*/
	}

	if (number>0) {
		cohesionVector[0] = (cohesionVector[0]/number)/* - Coordinates[0]*//COH_CONST;
		cohesionVector[1] = (cohesionVector[1]/number)/* - Coordinates[1]*//COH_CONST;
	}

	return cohesionVector;
}

void update_main_vector(int n,Bird& b) {

	b.separationVec = separation(n,b);
	b.cohesionVec = cohesion(n,b);
	b.alignmentVec = alignment(n,b);

	b.MainVector[0] = b.Velocity[0] + b.separationVec[0] + b.cohesionVec[0] + b.alignmentVec[0];
	b.MainVector[1] = b.Velocity[1] + b.separationVec[1] + b.cohesionVec[1] + b.alignmentVec[1];
	
	std::ofstream fil("Vectors.txt", std::ios_base::app | std::ios_base::out);
	fil << "Bird (Type , Instance) : ( " << b.type << " , " << b.No << " )\n";
	fil << "Position : ( " << b.Coordinates[0] << " , " << b.Coordinates[1] << " )\n";
	fil << "Separation : ( " << b.separationVec[0] << " , " << b.separationVec[1] << " )\n";
	fil << "Alignment : ( " << b.alignmentVec[0] << " , " << b.alignmentVec[1] << " )\n";
	fil << "Cohesion : ( " << b.cohesionVec[0] << " , " << b.cohesionVec[1] << " )\n";
	fil << "Old Velocity : ( " << b.Velocity[0] << " , " << b.Velocity[1] << " )\n";

	for(int d =0; d<2; d++) {
		if(b.MainVector[d] > b.maxVel)
			b.MainVector[d] = b.maxVel;
		else if(b.MainVector[d] < -b.maxVel)
			b.MainVector[d] = -b.maxVel;
		else;
	}

	if (b.MainVector[0] == 0 && b.MainVector[1] == 0) {
		if (b.Velocity[0]<0)
			b.MainVector[0] = -1;
		if (b.Velocity[0]>0)
			b.MainVector[0] = 1;
		if (b.Velocity[1]<0)
			b.MainVector[1] = -1;
		if (b.Velocity[1]>0)
			b.MainVector[1] = 1;
	}
	fil << "Main Vector (New Velocity) : ( " << b.MainVector[0] << " , " << b.MainVector[1] << " )\n";
	fil << "\n";

}

bool compare_states(State s1,State s2) {
	for (int i=0;i<bird_size;i++) {
		if (!(s1[i].equate_boids(s2[i])))
			return false;
	}
	return true;
}

void save_state() {
	states.push_back(Flock);
}

void read_initial_condition_file(int t) {

	ifstream initial_condition_file;
	initial_condition_file.open("initial_conditions_" + to_string(t) + ".txt", ios::in);

	if (initial_condition_file.is_open()) {

		for (int i=0;i<bird_types;i++) {
			for (int j=0;j<bird_num[i];j++) {

				Bird b;
				b.resize_vectors();
				initial_condition_file >> b.Coordinates[0];
				initial_condition_file >> b.Coordinates[1];
				initial_condition_file >> b.Velocity[0];
				initial_condition_file >> b.Velocity[1];

				if (grid[b.Coordinates[1]][b.Coordinates[0]] == 0)
					grid[b.Coordinates[1]][b.Coordinates[0]] = i+1;
				else
					cout << "Problem! Trying to put 2 birds in same cell \n";

				b.maxVel = max_vel[i];
				b.type = (i+1);
				b.No = (j+1);
				Flock.push_back(b);
			}
		}
		initial_condition_file.close();
	}
	else
		cout << "Cannot open Initial Condition File \n";
}

void generate_initial_condition(int t) {

	std::ofstream initial_file("initial_conditions_" + to_string(t) + ".txt", /*std::ios_base::app | */std::ios_base::out);
	for(int i=0;i<bird_types;i++) {
		for (int j=0;j<bird_num[i];j++) {

			Bird b;
			b.resize_vectors();
			b.Velocity[0] = (rand() % (2*max_vel[i]+1)) - max_vel[i]; 
			b.Velocity[1] = (rand() % (2*max_vel[i]+1)) - max_vel[i]; 

			while (1) {
				b.Coordinates[0] = rand() % grid_height; 
				b.Coordinates[1] = rand() % grid_width; 
				if (grid[b.Coordinates[1]][b.Coordinates[0]] == 0) {
					grid[b.Coordinates[1]][b.Coordinates[0]] = i+1;
					break;
				}
			}
			b.maxVel = max_vel[i];
			b.type = (i+1);
			b.No = (j+1);
			Flock.push_back(b);
			initial_file << b.Coordinates[0] << " " << b.Coordinates[1] << " " << b.Velocity[0] << " " << b.Velocity[1] << endl;
		}
	}			
	initial_file.close();
}

void create_log(int t) {

	std::ofstream vis("Visual_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);
	vis << bird_size << " " << grid_width << " " << grid_height << endl;			
	vis.close();

}

void generate_output_file(int t) {

	ofstream output_file;
	output_file.open("output_file_" + to_string(t) + ".txt", ios::out);

	if (simulation_type==1)
		output_file << "Simulation Type - Fix number of Iterations\n";
	else
		output_file << "Simualtion Type - Fix execution time\n";

	output_file << "Types of Birds : " << bird_types << endl;
	for (int i=0;i<bird_types;i++)
		output_file << "Instance of Bird " << i+1 << " : " << bird_num[i] << endl;
	output_file << "Total Number of Birds : " << bird_size << endl; 
	output_file << "Grid Dimensions : " << grid_height << "x" << grid_width << endl;
	output_file << "Detection Range : " << detectionRange << endl;
	output_file << "Separation Range : " << separationRange << endl;
	output_file << "Separation Constant : " << SEP_CONST << endl;
	output_file << "Alignment Constant : " << ALIGN_CONST << endl;
	output_file << "Cohesion Constant : " << COH_CONST << endl;
	for (int i=0;i<bird_types;i++)
		output_file << "Maximum Velocity of Bird " << i+1 << " : " <<  max_vel[i] << endl;

	if (simulation_type==1)
		output_file << "Number of Iterations : " << iteration_size << endl;
	else
		output_file << "Execution Time : " << simulation_time << endl;

	output_file << "Number of groups in DOI Distribution : " << groups << endl;
	if (DOIS[2]==1)
		output_file << "History Size : " << history_size << endl;

	output_file << "\nInitial Locations ...\n" << endl;
	for(int z=0; z<grid_height; z++) {
		for(int x=0; x<grid_width; x++)
			output_file << " " << grid[z][x] << " ";	
		output_file << endl;
	}

	output_file << "\nNo of states in L-whole : " << Lwhole_states << endl;

	if (DOIS[0]==1) {
		output_file << "\nDOI Distribution with Min-Distance Method ..." << endl;
		output_file << "Max DOI value : " << max_doi_dist << endl;
		output_file << "Min DOI value : " << min_doi_dist << endl;
		output_file << "States in L-zero : " << Lzero_dist << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/(double)groups << " - " << (double)(i+1)/(double)groups << ") : " << dist[i] << endl;
		output_file << "No of states in L-emergent : " << Lwhole_states-Lzero_dist << endl;

		run_data_Lemergent_dist[t-1] = (Lwhole_states-Lzero_dist);
	}
	if (DOIS[1]==1) {
		output_file << "\nDOI Distribution with Interaction Count Method ..." << endl;
		output_file << "Max DOI value : " << max_doi_ic << endl;
		output_file << "Min DOI value : " << min_doi_ic << endl;
		output_file << "States in L-zero : " << Lzero_ic << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/(double)groups << " - " << (double)(i+1)/(double)groups << ") : " << ic[i] << endl;
		output_file << "No of states in L-emergent : " << Lwhole_states-Lzero_ic << endl;

		run_data_Lemergent_ic[t-1] = (Lwhole_states-Lzero_ic);
	}
	if (DOIS[2]==1) {
		output_file << "\nDOI Distribution with on-the-fly Interaction Count Method ..." << endl;
		output_file << "States in L-zero : " << Lzero_otf << endl;
		for (int i=0;i<groups;i++)
			output_file << "\t D in (" << (double)i/(double)groups << " - " << (double)(i+1)/(double)groups << ") : " << otf[i] << endl;
		output_file << "No of states in L-emergent : " << Lwhole_states-Lzero_otf << endl;

		run_data_Lemergent_otf[t-1] = (Lwhole_states-Lzero_otf);
	}
}

void DOI_dist(int t) {

	for (int i=0;i<groups;i++)
		dist[i]=0;

	double k,ans;
	vector<double> arr(Lwhole_states);
	ifstream min_dist_file;
	min_dist_file.open("MinDist_" + to_string(t) + ".txt", ios::in);

	if (min_dist_file.is_open()) {

		for(int i=0;i<Lwhole_states;i++) {
			min_dist_file >> k;
			arr[i]=k; 
		}
		min_doi_dist=arr[0];
		for(int i=0; i<Lwhole_states; i++) {
			if(max_doi_dist < arr[i])
				max_doi_dist=arr[i];
			if(min_doi_dist > arr[i])
				min_doi_dist=arr[i];
		}

		for (int i=0;i<Lwhole_states;i++) {
			ans = 1-((arr[i]-min_doi_dist)/(max_doi_dist-min_doi_dist));
			if(ans == 0)
				Lzero_dist++;
			for (int j=1;j<=groups;j++) {
				if (ans >= (double)j/(double)groups) {
					if (j==groups) {
						dist[groups-1]++;
						break;
					}
				}
				else {
					dist[j-1]++;
					break;
				}
			}
		}
	}
	for (int i=0;i<groups;i++)
		run_data_dist[i][t-1] = dist[i];
}

void DOI_ic(int t) {

	for (int i=0;i<groups;i++)
		ic[i]=0;

	int k;
	double ans;
	vector<int> arr(Lwhole_states);
	ifstream interaction_count_file;
	interaction_count_file.open("InteractionCount_" + to_string(t) + ".txt", ios::in);

	if (interaction_count_file.is_open()) {

		for(int i=0;i<Lwhole_states;i++) {
			interaction_count_file >> k;
			arr[i]=k; 
		}
		min_doi_ic=arr[0];
		for(int i=0; i<Lwhole_states; i++) {
			if(max_doi_ic < arr[i])
				max_doi_ic=arr[i];
			if(min_doi_ic > arr[i])
				min_doi_ic=arr[i];
		}

		for (int i=0;i<Lwhole_states;i++) {
			ans = (double)(arr[i]-min_doi_ic)/(double)(max_doi_ic-min_doi_ic);
			if(ans == 0)
				Lzero_ic++;
			for (int j=1;j<=groups;j++) {
				if (ans >= (double)j/(double)groups) {
					if (j==groups) {
						ic[groups-1]++;
						break;
					}
				}
				else {
					ic[j-1]++;
					break;
				}
			}
		}
	}
	for (int i=0;i<groups;i++)
		run_data_ic[i][t-1] = ic[i];
}

void initiate_history() {

	history.resize(history_size);
	for (int i=0;i<groups;i++)
		otf[i]=0;
	for (int i=0;i<history_size;i++)
		history[i]=-1;
}

void find_repeated_state() {
	int flag=0;
	for (int i=0;i<iteration_size-1;i++) {
		for (int j=i+1;j<iteration_size;j++) {
			if (compare_states(states[i],states[j])) {
				Lwhole_states=j;
				cout << "Equal : " << i << " and " << j << endl;
				cout << "Size of Lwhole : " << Lwhole_states << endl;
				flag=1;
				break;
			}
		}
		if (flag==1)
			break;
	}
	if (flag==0)
		Lwhole_states=0;
}

void clear_grid() {
	for(int z=0;z<grid_height;z++) {
		for(int x=0; x<grid_width; x++)
			grid[z][x] = 0;		
	}
}

void update_grid() {
	clear_grid();
	for (int i=0;i<Flock.size();i++) {
		Bird& b = Flock[i];
		grid[b.Coordinates[1]][b.Coordinates[0]] = b.type;
	}
}

void add_to_log(int t,int n) {

	std::ofstream vis("Visual_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);

	for (int i=0;i<Flock.size(); i++){
		Bird& a = Flock[i];
		int type;
		if (a.type==1 || a.type==2)
			type=a.type;
		else
			type=0;
		vis << (n+1) << " " << a.Coordinates[0] << " " << a.Coordinates[1] << " " << type << endl;

	}
	if (n==iteration_size-1)
		vis << -1 << endl;
	vis.close();

}

void calculate_next_state(int n) {
		
	std::ofstream fil("Vectors.txt", std::ios_base::app | std::ios_base::out);
	fil << "--------------------------------------------------\n";
	if (n != iteration_size-1)
		fil << "Iteration : " << n+1 << "\n";
	
	for(int i=0; i<Flock.size(); i++) {	
		Bird& a = Flock[i];
		update_main_vector(n,a);
	}
	fil << "--------------------------------------------------\n\n";
	fil.close();
}

void avoidCollision() {

	for (int i=0;i<Flock.size();i++) {
		for (int j=0;j<Flock.size();j++) {
			Bird& a = Flock[i];
			Bird& b = Flock[j];

			if (i!=j && a.Coordinates[0] == b.Coordinates[0] && a.Coordinates[1] == b.Coordinates[1]) {

				int g=signums(b.Velocity[0]);
				int h=signums(b.Velocity[1]);

				/*std::ofstream tmp("Temporary.txt", std::ios_base::app | std::ios_base::out);
				tmp << "From : " << b.Coordinates[0] << "  " << b.Coordinates[1] << "\n";*/

				b.Coordinates[0] += g;
				b.Coordinates[1] += h;

				/*tmp << "To : " << b.Coordinates[0] << "  " << b.Coordinates[1] << "\n\n";*/

				if (b.Coordinates[0]==-1)
					b.Coordinates[0] += grid_width;
				if (b.Coordinates[0]==grid_width)
					b.Coordinates[0] -= grid_width;
				if (b.Coordinates[1]==-1)
					b.Coordinates[1] += grid_height;
				if (b.Coordinates[1]==grid_height)
					b.Coordinates[1] -= grid_height;
			}
		}
	}
}

void update_new_state(int n) {

	for(int i=0; i<Flock.size(); i++) {	
		Bird& a = Flock[i];
		a.update_attributes(n);
	}
	for (int i=0;i<10;i++)
		avoidCollision();
}

void interaction_count(int t) {

	int totalCount=0;

	for(int i=0; i<Flock.size(); i++){
		for(int j=0; j<Flock.size(); j++){
			Bird& b1 = Flock[i];
			Bird& b2 = Flock[j];
			if (b1.type == b2.type && b1.No == b2.No);
			else {
				if(magnitude(b1.closestVect(b2)) <= detectionRange)
					totalCount++;
			}
		}
	}

	std::ofstream counter("InteractionCount_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);
	counter << (totalCount/2) << endl;		//To Avoid Double Counting
	counter.close();

}

void on_the_fly_ic(int t) {

	int totalCount=0,interCount;

	for(int i=0; i<Flock.size(); i++){
		for(int j=0; j<Flock.size(); j++){
			Bird& b1 = Flock[i];
			Bird& b2 = Flock[j];
			if (b1.type == b2.type && b1.No == b2.No);
			else {
				if(magnitude(b1.closestVect(b2)) <= detectionRange)
					totalCount++;
			}
		}
	}
	interCount = totalCount/2;			//To Avoid Double Counting
	std::ofstream counter("OnTheFly_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);

	int flag=0;
	int max_doi_otf = INT_MIN;
	int min_doi_otf = INT_MAX;
	
	int i;
	for(i=0; i<history_size; i++) {
		if (history[i]!=-1) {
			if(max_doi_otf < history[i])
				max_doi_otf=history[i];
			if(min_doi_otf > history[i])
				min_doi_otf=history[i];
		}
		else {
			flag=1;
			break;
		}
	}

	if (flag==0) {
		if (max_doi_otf != min_doi_otf) {
			double ans = (double)(interCount-min_doi_otf)/(double)(max_doi_otf-min_doi_otf);		//DIVIDEBYZERO POSSIBILITY
			counter << ans << endl;
			if(ans == 0)
				Lzero_otf++;
			for (int j=1;j<=groups;j++) {
				if (ans >= (double)j/(double)groups) {
					if (j==groups) {
						otf[groups-1]++;
						break;
					}
				}
				else {
					otf[j-1]++;
					break;
				}
			}
		}
		else {
			otf[groups-1]++;
			counter << 1 << endl;
		}
		for (i=0;i<groups;i++)
			run_data_otf[i][t-1] = otf[i];
	}
	
	if (flag==1 && i==history_size-1) {
		
		for (i=history_size-1;i>0;i--)
			history[i]=history[i-1];
		history[0] = interCount;
		
		max_doi_otf = INT_MIN;
		min_doi_otf = INT_MAX;
		
		for(i=0; i<history_size; i++) {
			if(max_doi_otf < history[i])
				max_doi_otf=history[i];
			if(min_doi_otf > history[i])
				min_doi_otf=history[i];
		}
		
		for (i=history_size-1;i>=0;i--) {
			if (max_doi_otf != min_doi_otf) {
				double ans = (double)(history[i]-min_doi_otf)/(double)(max_doi_otf-min_doi_otf);		//DIVIDEBYZERO POSSIBILITY
				counter << ans << endl;
				if(ans == 0)
					Lzero_otf++;
				for (int j=1;j<=groups;j++) {
					if (ans >= (double)j/(double)groups) {
						if (j==groups) {
							otf[groups-1]++;
							break;
						}
					}
					else {
						otf[j-1]++;
						break;
					}
				}
			}
			else {
				otf[groups-1]++;
				counter << 1 << endl;
			}
		}
		for (i=0;i<groups;i++)
			run_data_otf[i][t-1] = otf[i];
	}
	else {
		for (i=history_size-1;i>0;i--)
			history[i]=history[i-1];
		history[0] = interCount;
	}
	//counter << "History : ";
	//for (i=0;i<history_size;i++)
	//	counter << history[i] << " ";
	//counter << endl;
	counter.close();
}

void min_distance(int t) {

	double min,count=0;
	int flag=0;

	for(int i=0; i<Flock.size(); i++){
		for(int j=0; j<Flock.size(); j++){
			Bird& b1 = Flock[i];
			Bird& b2 = Flock[j];
			if (b1.type == b2.type && b1.No == b2.No);
			else {
				if (flag==0) {
					min = magnitude(b1.closestVect(b2));
					flag=1;
				}
				else {
					if (min > magnitude(b1.closestVect(b2)))
						min = magnitude(b1.closestVect(b2));
				}
			}
		}
		count += min;
	}
	std::ofstream min_dist_file("MinDist_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);
	min_dist_file << count << endl;	
	min_dist_file.close();
}

void record_state(int t, int n) {

	std::ofstream state_file1("Times_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);
	std::ofstream state_file2("Grids_" + to_string(t) + ".txt", std::ios_base::app | std::ios_base::out);
	if (n != iteration_size) {
		state_file1 << "Iteration : " << n << "\n\n";
		state_file2 << "\nIteration : " << n << "\n\n";
	}
	else {
		state_file1 << "\nfinal state\n\n";
		state_file2 << "\nfinal state\n\n";	
	}

	for(int z=0;z<grid_height;z++){
		for(int x=0; x<grid_width; x++) {
			if (grid[z][x] != 0)
				state_file2 << " " << grid[z][x] << " |";
			else
				state_file2 << "   |";
		}
		state_file2 << "\n";
	}

	for(int i=0; i<Flock.size(); i++) {	
		Bird& a = Flock[i];
		state_file1 << "Bird (type , instance) : ( " << a.type << " , " << a.No << " )\n"; 
		state_file1 << "Position : ( " << a.Coordinates[0] << " , " << a.Coordinates[1] << " )\n"; 
		state_file1 << "Velocity : ( " << a.Velocity[0] << " , " << a.Velocity[1] << " )\n\n";
	}
}

int main(int args,char** argv) {

	srand(time(NULL));
	clock_t init, final;
	init=clock();

	read_input_parameter_file();
	resize_mainrun_vectors();

	for (int t=1;t<=simulation_count;t++) {

		if (mode == 1)
			read_initial_condition_file(t);
		else if (mode == 2)
			generate_initial_condition(t);
		else {
			cout << "Invalid Mode";
			exit;
		}

		std::cout << "-- Boids Model Simulator --" << endl;
		std::cout << "Simulation Count : " << t << endl;

		int n=0;
		if (simulation_type == 1) {
			while (n<iteration_size) {

				save_state();
				record_state(t,n);
				if (n==0)
					create_log(t);
				add_to_log(t,n);

				calculate_next_state(n);
				update_new_state(n);
				update_grid();

				if (n==0)
					initiate_history();

				if (DOIS[0] == 1)
					min_distance(t);
				if (DOIS[1] == 1)
					interaction_count(t);
				if (DOIS[2] == 1)
					on_the_fly_ic(t);				
				n++;
			}
			find_repeated_state();
			if (DOIS[0]==1)
				DOI_dist(t);
			if (DOIS[1]==1)
				DOI_ic(t);

			generate_output_file(t);
		}
		else if (simulation_type == 2) {
			while (clock() <= simulation_time){

				clear_grid();
				save_state();
				record_state(t,n);

				calculate_next_state(n);
				update_new_state(n);
				update_grid();

				if (n==0)
					initiate_history();

				if (DOIS[0] == 1)
					min_distance(t);
				if (DOIS[1] == 1)
					interaction_count(t);
				if (DOIS[2] == 1)
					on_the_fly_ic(t);	
				n++;
			}
			find_repeated_state();
			if (DOIS[0]==1)
				DOI_dist(t);
			if (DOIS[1]==1)
				DOI_ic(t);

			generate_output_file(t);
		}			
	}

	final=clock()-init;
	generate_final_statistics_file(final);

	std::cout << "Time Taken for completion: " << (double)final / ((double)CLOCKS_PER_SEC) << endl;
	std::getchar();

}