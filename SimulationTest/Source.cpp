#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <engine.h>

#pragma comment (lib,"libmx.lib")
#pragma comment (lib,"libeng.lib")
#pragma comment (lib,"libmex.lib")
#pragma comment (lib,"libmat.lib")


using namespace std;

const double radius = 1;
const double pi = 3.1415926535897;
const double stepSize = .001;
double xMin = -5, xMax = -xMin, yMin = xMin, yMax = xMax;
const double kb = 1.38e-23, T = 300, eta = 1e-3;
const double D = kb*T / 6 / pi / eta / radius;
const double k = sqrt(2 * D*stepSize);

const double mass = 1.;
const double g=-9.81;


struct Vec3 {

	double x, y, z;
	Vec3() {}
	Vec3(double x, double y, double z) :x(x), y(y), z(z) { }
public:
	Vec3 operator+(const Vec3& other) const {
		Vec3 ans(other.x + x, other.y + y, other.z + z);
		return ans;
	}
	Vec3 operator-(const Vec3& other) const{
		Vec3 ans(-other.x + x, -other.y + y, -other.z + z);
		return ans;
	}
	double operator*(const Vec3& other) const{
		return other.x*x + other.y*y + other.z*z;
	}
	Vec3 operator*(const double scalar) const{
		return Vec3(x*scalar, y*scalar, z*scalar);
	}
	void operator-=(const Vec3& other)  {
		x -= other.x;
		y -= other.y;
	}
	void operator+=(const Vec3& other) {
		x += other.x;
		y += other.y;
	}
	bool operator==(const Vec3 & other) {
		return x == other.x && y == other.y && z == other.z;
	}
	ostream& print(ostream& stream) const {
		return stream << x << " " << y << " " << z << "\n";
	}
	
};

ostream& operator<<(ostream& stream,const Vec3 &vec){
	return vec.print(stream);
}

vector<Vec3> pairDirs;

Vec3 operator*(const double scalar, const Vec3& vector)  {
	return vector*scalar;
}


bool any(const Vec3& vec) {
	return !(vec.x == 0 && vec.y == 0 && vec.z == 0);
}

double normsq(const Vec3 vec) {
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}
double norm(const Vec3& vec) {
	return sqrt(normsq(vec));
}

Vec3 normr(const Vec3& vec) {
	double invModulus = 1/norm(vec);
	return invModulus*vec;
}

Vec3 cross(const Vec3& A, const Vec3 &B) {
	return Vec3(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x);
}


Vec3 periodicBound(const Vec3& position, const bool inclRadius) {
	Vec3 shift(0, 0, 0);
	double eps = inclRadius ? 2 * radius : 0;
	if (position.x > xMax - eps)
		shift.x = -xMax + xMin;
	else if (position.x < xMin + eps)
		shift.x = xMax - xMin;

	if (position.y > yMax - eps)
		shift.y = -yMax + yMin;
	else if (position.y < yMin + eps)
		shift.y = yMax - yMin;

	return shift;
}



struct Particle {
	int number;
	Vec3 position;
	bool deposed = false;
	int numbOfHalt = 0;
	unordered_map <Particle*,Vec3> touching;
	Particle(){}
	Particle(Vec3 pos,int n):position(pos),number(n) {}

	double getTimeToColl(const Vec3& posOfColl, const Vec3& move, const double deltaTime, const Vec3 &shift, Vec3& tShift) {
		Vec3 vel = move*(1 / deltaTime);
		Vec3 deltaR = posOfColl - position;
		tShift = Vec3(0,0,0);
		double t = getTimeToColl(deltaR, vel);
		double finalT = -1.;
		
		if (t > 0) {
			finalT = t;
		}

		if (shift.x != 0 && shift.y != 0) {
			t = getTimeToColl(deltaR - shift, vel);
			if (t > 0 && (finalT == -1 || t < finalT)) {
				finalT = t;
				tShift = shift;
			}
		}
		if (shift.x != 0) {
			t = getTimeToColl(deltaR - Vec3(shift.x, 0, 0), vel);
			if (t > 0 && (finalT == -1 || t < finalT)) {
				tShift = Vec3(shift.x, 0, 0);
				finalT = t;
			}
		}
		if (shift.y != 0) {
			t = getTimeToColl(deltaR - Vec3(0, shift.y, 0), vel);
			if (t > 0 && (finalT == -1 || t < finalT)) {
				tShift = Vec3(0,shift.y,0);
				finalT = t;
			}

		}
		if (finalT < 0){
			if (t > 0)
				cout << "error\n";
			tShift.x = 0;
			tShift.y = 0;
		}

	

		return finalT>0?finalT:t;
	}


	ostream& print(ostream& stream) const {
		return stream << position.x << " " << position.y << " " << position.z << "\n";
	}
private:
	double getTimeToColl(const Vec3& deltaR, const Vec3& vel) {
		if (normsq(deltaR) < 4 * radius*radius)
			cout << "w\n";
		double discriminant = pow(deltaR*vel, 2) - normsq(vel)*(normsq(deltaR) - 4 * radius*radius);
		if (discriminant != discriminant)
			cout << "nan found\n";
		cout << discriminant<<"\n";
		if (discriminant >=0)
			return (deltaR*vel - sqrt(discriminant)) / normsq(vel);
		else return -1.;
	}
};

ostream& operator<<(ostream& stream, const Particle& p) {
	return p.print(stream);
}


double randomDouble(double min, double max) {
	return min+(double)rand() / RAND_MAX*(max-min);
}



typedef pair<Particle*, Vec3> PVecPair;
typedef unordered_map<Particle*, Vec3> PVecMap;

bool fixTouching(const Particle& particle, Vec3& move) {
	int size = particle.touching.size();
	Vec3 newMove;
	
	if (size > 0) {
		Vec3 dir = normr((*(particle.touching.begin())).first->position - particle.position- (*(particle.touching.begin())).second);
		if (dir*move > 0)
			newMove = move - dir*(dir*move);
		else
			newMove = move;
		if (size > 1) {
			bool moveWrong = false,wrongRollMove=false;
			pairDirs.clear();
			Vec3 rollMove(0, 0, 0);
			for (auto iter = particle.touching.begin(); iter != particle.touching.end(); iter++) {
				Vec3 dir1 = normr((*iter).first->position - particle.position-(*iter).second);
				auto iter1 = iter;
				advance(iter1, 1);
				for (iter1; iter1 != particle.touching.end(); iter1++) {
					Vec3 dir2 = (*iter1).first->position - particle.position-(*iter1).second;
					pairDirs.push_back(normr(cross(dir1, dir2)));
				}
				if (dir1*move > 0) {
					if (!any(rollMove))
						rollMove = move - dir1*(dir1*move);
					else {
						rollMove = Vec3(0, 0, 0);
						wrongRollMove = true;
					}
				}
			}
				for (auto iter = particle.touching.begin(); iter != particle.touching.end(); iter++) {
					Vec3 dir1 = (*iter).first->position- (*iter).second - particle.position;
					for (int k = 0; k < pairDirs.size(); k++) {
						if ((pairDirs[k] * (move*pairDirs[k]))*dir1>1e-14) {
							
							pairDirs.erase(pairDirs.begin() + k);
							k--;
							
						}
					}
					if (!wrongRollMove && rollMove*dir1 > 1e-14) {
						wrongRollMove = true;
					}

					if (move*dir1 > 1e-14)
						moveWrong = true;
				
				}
				if (moveWrong) {
					if (!wrongRollMove)
						newMove = rollMove;
					else {
						if (!pairDirs.empty()) {
							newMove = (pairDirs[0] * move)*pairDirs[0];
							if (pairDirs.size() > 1)
								cout << "size>1\n";
						}
						else {
							//newMove = Vec3(0, 0, 0);
							return false;
						}
					}
				}
				else
					newMove = move;

				
			}
		for (auto iter = particle.touching.begin(); iter != particle.touching.end(); iter++) {
			Vec3 dir1 =(*iter).first->position- (*iter).second - particle.position;
			double prod = newMove*dir1;
			if (prod > 1e-14)
				cout << "wrong direction\n";
		}
		move = newMove;
		}
	return true;

	
}

void printToFile(const vector<Particle> &particles,string file) {
	ofstream mFile;
	mFile.open(file);
	if (mFile.is_open()) {
		for each (Particle p in particles)
		{
			mFile << p.position.x << " " << p.position.y << " " << p.position.z << "\n";
		}
		mFile.close();
	}
	else
		cout << "cant open file\n";
}

Vec3 getMove() {
	
	double dx = k*randomDouble(0,1) ;
	double dy = k*randomDouble(0, 1) ;
	double dz = k*randomDouble(0, 1) + D / kb / T*mass*g*stepSize;
	return Vec3(dx, dy, dz);
}

void drawSphereMatlab(Engine* m_pEngine, mxArray* matArray,const Vec3& position, int i) {
	

	double* pToMatArray = mxGetPr(matArray);
	*pToMatArray = position.x;
	*(pToMatArray + 1) = position.y;
	*(pToMatArray + 2) = position.z;
	engPutVariable(m_pEngine, "positions", matArray);
	string del = "delete(f(" + to_string(i+1) + "));";
	engEvalString(m_pEngine,del.c_str());
	string radiusS = to_string(radius);
	string s = "f(" + to_string(i+1) + ")=surf(x*" + radiusS +
		"+positions(1),y*" + radiusS + "+positions(2),z*" + radiusS + "+positions(3));";
	engEvalString(m_pEngine, s.c_str());
	engEvalString(m_pEngine, "hold on");
}

int main() {

	Engine* m_pEngine;
	m_pEngine = engOpen("null");
	engEvalString(m_pEngine, "figure");
	engEvalString(m_pEngine, "axis('equal')");
	engEvalString(m_pEngine, "[x y z]=sphere;");


	engEvalString(m_pEngine, "xlabel('x')");
	engEvalString(m_pEngine, "ylabel('y')");
	engEvalString(m_pEngine, "hold on");

	
	srand(time(NULL));
	vector<Particle> particles;
	
	vector <pair<int,Vec3>> simultColl;
	double height = 10;
	const int time = 20, rate = 2;
	const int limitOfHalt = 3;
	string s = "f=gobjects(" + to_string(time*rate) + ");";
	engEvalString(m_pEngine, s.c_str());
	mxArray* matArray = mxCreateDoubleMatrix(3, 1, mxREAL);
	string limits = "axis([2*" + to_string(xMin) + ",2*" + to_string(xMax) + ",2*" + to_string(yMin) + ",2*" + to_string(yMax) + ",0," + to_string(height) + "]);";
	engEvalString(m_pEngine, limits.c_str());
	//debug
	//double spx[] = { 1.8,0.6,1.2, -0.9, -1.3, -0.5, 1.7, -1.4, 0.8,1.5,1.4,-1.6};
	//double spy[] = { 0.2,-1.7,0.2, 1.4, 0.2, -0.1, 1.8, -0.4, 0.6,-0.1,1.4,1.3 };
	//Vec3 sp[time*rate];
	int n = 0;
	particles.reserve(time*rate);
	vector<pair<int,double>> tested;

	/*ifstream mfile;
	mfile.open("sp.txt");
	if (mfile.is_open()) {
		for (int i = 0; i < time*rate; i++) {
			
			mfile >> sp[i].x>>sp[i].y;
			sp[i].z = height;
		}
		mfile.close();
	}
	else
		cout << "no sp.txt file\n";
	
*/



	for (int t = 0; t < time; t++) {
		for (int p = 0; p < rate; p++) {
			Vec3 pos(randomDouble(xMin, xMax), randomDouble(yMin, yMax), height);
			//Vec3 pos=sp[n]; //debug
			
			particles.push_back(Particle(pos,n++));
			double simulTime = 1 /(double) rate;
			while (simulTime > 0)
			{
				for (int i = 0; i < particles.size(); i++) {
					if (!particles[i].deposed) {
						Vec3 move = getMove();
						double deltaTime = stepSize;
						int numbOfIter = 0;
						
						while (deltaTime > 0 && any(move) && numbOfIter < 5) {
							numbOfIter++;
							int coll = -1;
							double minT = deltaTime;
							double dt;
							
							if (fixTouching(particles[i], move)) {
								Vec3 vel = move*(1 / deltaTime);
								Vec3 shift(0, 0, 0), fShift(0, 0, 0), tShift(0, 0, 0);
								shift = periodicBound(particles[i].position, true);
								Vec3 shift1 = periodicBound(particles[i].position + move, true);
								shift.x = shift1.x != 0 ? shift1.x : shift.x;
								shift.y = shift1.y != 0 ? shift1.y : shift.y;
								tested.clear();
								for (int j = 0; j < particles.size(); j++) {
									if (i != j && particles[i].touching.find(&particles[j]) == particles[i].touching.end()) {
										cout << i << " " << j << "\n";
										//outputs NaN for some reason
											dt = particles[i].getTimeToColl(particles[j].position, move,deltaTime,shift,tShift);
											tested.push_back(make_pair(j, dt));
											if (dt > 0 && dt < minT - 1e-14) {
												minT = dt;
												
												coll = j;
												simultColl.clear();
												fShift = tShift;
											}
											else {
												if (dt > 0 && dt < minT + 1e-14) {
													simultColl.push_back(make_pair(j,tShift));
												}
												else if (dt!=-1. && dt <= 0 && norm(particles[i].position - particles[j].position-tShift) <= 2 * radius)
													cout << "this cant happen\n";
											}

										
									}
								}
								if (particles[i].position.z + move.z < radius) {
									dt = (radius - particles[i].position.z) / vel.z;
									if (dt < minT) {
										minT = dt;
										coll = -2;
										particles[i].deposed = true;
										deltaTime = 0;
									}
								}
								if (coll != -1) {
									deltaTime = deltaTime - minT;
									particles[i].position = particles[i].position + vel*minT;
									if (coll > -1) {
										cout << i << " collided with " << coll << " \n";
										Vec3 dir = normr(particles[coll].position - particles[i].position-fShift);
										move = move - dir*(dir*move);
										particles[i].touching.insert(PVecPair(&particles[coll], fShift));
										for (auto pair : particles[i].touching) {
											if (pair.first->number < 0)
												cout << "negative index\n";
										}
										if (!particles[coll].deposed)
											particles[coll].touching.insert(PVecPair(&particles[i],(-1.)* fShift));

										if (!simultColl.empty()) {
											for (int w = 0; w < simultColl.size(); w++) {
												particles[i].touching.insert(PVecPair(&particles[simultColl[w].first], simultColl[w].second));
												if (!particles[simultColl[w].first].deposed)
													particles[simultColl[w].first].touching.insert(PVecPair(&particles[i], (-1.) * simultColl[w].second));
											}
										}
									}
								}
								else {
									particles[i].position = particles[i].position + move;
									deltaTime = 0;
									

								}
								Vec3 boundShift = periodicBound(particles[i].position, false);
								particles[i].position += boundShift;
								for (PVecMap::iterator iter = particles[i].touching.begin(); iter != particles[i].touching.end();) {
										(*iter).second -= boundShift;
										if(!(*iter).first->deposed && any(boundShift))
											(*iter).first->touching.at(&particles[i]) += boundShift;

									
									double dist = norm((*iter).first->position - (*iter).second - particles[i].position);
									if (dist>(2 + 1e-5)*radius &&!(coll > -1 && (*iter).first == &particles[coll])) {

										if (!(*iter).first->deposed) {
											cout << "erased " << i << " from " << (*iter).first->number <<"dist = "<<dist<< "\n";
											(*iter).first->touching.erase(&particles[i]);

										}
										cout << "erased " << (*iter).first->number << " from " << i << "dist = " << dist << "\n";
										iter = particles[i].touching.erase(iter);

									}
									else {
										iter++;
									}
								}
								for (auto pair : particles[i].touching) {
									if (pair.first->number < 0)
										cout << "negative index\n";
								}
							
								for (int j = 0; j < particles.size(); j++) {
									Vec3 deltaR = particles[i].position - particles[j].position;
									double dist = norm(deltaR), dist1 = norm(deltaR + shift), dist2 = norm(deltaR + Vec3(shift.x, 0, 0)), dist3 = norm(deltaR + Vec3(0, shift.y, 0));

									if (i != j && (dist < 2 * radius || dist1<2*radius || dist2<2*radius||dist3<2*radius)) {
										if (particles[i].touching.find(&particles[j]) == particles[i].touching.end() || dist < (2 - 1e-2)*radius) {

											cout << "error \n";
										}
										//else cout << particles[i]<<particles[j];
									}
								}
								drawSphereMatlab(m_pEngine, matArray, particles[i].position, i);
								////debug
								//printToFile(particles, "results.txt");
								//cout << "writing to file complete\n";
							}
							else {
								if (++particles[i].numbOfHalt > limitOfHalt)
									particles[i].deposed = true;
								break;
							}
						}
					}
				}
		
				simulTime -= stepSize;
			}

			printf("%.2f %% completed\n", ((double)particles.size() / time / rate * 100));
		}
	}
	
	//engEvalString(m_pEngine, "close");
	//printToFile(particles, "results.txt");
	return 0;
}