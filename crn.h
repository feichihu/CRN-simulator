//
// Created by willf on 4/6/2020.
//

#ifndef CRN_SIMULATOR_CRN_H
#define CRN_SIMULATOR_CRN_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <queue>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
#include <chrono>
#include "matplotlibcpp.h"

using namespace std;
using namespace boost;
using namespace std::chrono;
namespace plt = matplotlibcpp;

struct Species {
    Species(){};
    Species(const string& s, const int fa){
        name = s;
        factor = fa;
    }
    string name;
    int factor = 1;
};
struct Slot{
    Slot(int& a, double& b){
       tau = b;
       reaction_idx = a;
    }
    Slot(){};
    int reaction_idx=-1;
    double tau=0.0;
    int seq = 0;
};
struct comparator{
    bool operator()(const Slot& a, const Slot& b){
       return a.tau>b.tau;

    }
};

class PriorityQ{
    //priority queue with update keys
public:
    PriorityQ();

    PriorityQ(vector<Slot>& s);

    void push(Slot a);

    Slot pop();

    Slot top();

    void print();

    void clean();

    bool verifyTop();


private:
    priority_queue<Slot,vector<Slot>,comparator> pq;
    vector<int> cur_seq;
};


class Reaction {
public:
    Reaction();


    Reaction(const string &r, const string &p, double l, int id);

    inline vector <Species> parseSpecies(string a);
    double calculateProp(unordered_map<string,int> &m);
    int adjustCount(unordered_map<string,int> &m);
    static void printSpecies(vector <Species> &v);
    void print();

    vector <Species> reactant;
    vector <Species> product;
    int idx;
    double lambda=1.0;
};

class CRN {
public:
    CRN();

    int addRxn(const string& reactant, const string& product, double lambda);

    int setConc(string species_name, int init_count);

    int simulate(int tmax, bool verbose=false, string mode="DM");//DM or NRM

    int plot(const string name);


    int print();

    void saveResult(const string& filename);

    void generateDependency();
    int clear();



private:
    class Frame{
    public:
        explicit Frame(double t){
            time = t;
        }
        double time;
        vector<int> counts;
    };
    unordered_map<string, int> s_count;
    vector <Reaction> reactions;
    vector<Frame> simulation_result;
    vector<string> species; //All species in the reaction
    vector<vector<int>> dependencyList;
    void saveFrame(double cur_time);



};


#endif //CRN_SIMULATOR_CRN_H
