//
// Created by willf on 4/6/2020.
//

#ifndef CRN_SIMULATOR_CRN_H
#define CRN_SIMULATOR_CRN_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <cstdlib>

using namespace std;
using namespace boost;

struct Species {
    string name;
    int factor = 1;
};




class Reaction {
public:
    Reaction();


    Reaction(const string &r, const string &p, double l);

    inline vector <Species> parseSpecies(string a);
    double calculateProp(unordered_map<string,int> &m);
    int adjustCount(unordered_map<string,int> &m);
    static void printSpecies(vector <Species> &v);
    void print();

    vector <Species> reactant;
    vector <Species> product;
    double lambda=1.0;
};

class CRN {
public:
    CRN();

    int addRxn(const string& reactant, const string& product, double lambda);

    int setConc(string species_name, int init_count);

    int simulate(int tmax);

    int print();

    void saveResult(const string& filename);

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
    void saveFrame(double cur_time);



};


#endif //CRN_SIMULATOR_CRN_H
