//
// Created by willf on 4/6/2020.
//

#include <iostream>
#include <string>
#include "crn.h"

using namespace std;

CRN::CRN() {

}

int CRN::addRxn(const string &reactant, const string &product, double lambda) {
    //x-1+x2, 2x2, 1.5
    reactions.emplace_back(reactant, product, lambda);

}

int CRN::setConc(string species_name, int init_count) {
    s_count[species_name] = init_count;
}

int CRN::simulate(int tmax) {
    cout<<"-------start simulation--------"<<endl;
    //C++11 random number generator (0,1)
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    //dist(mt) will be in (0,1)

    unordered_map<string,int> cur_count(s_count);
    double max_time = (double) tmax;
    double cur_time = 0;
    int num_reactions = reactions.size();
    vector<double> propensities(num_reactions);
    while(cur_time < max_time){
        cout<<"++++++++++at "<<cur_time<<"+++++++++"<<endl;
        double a0{};
        //propensity for each reaction
        for(int i=0;i<num_reactions;i++){
            propensities[i] = reactions[i].calculateProp(s_count);
            reactions[i].print();
            cout<<propensities[i]<<endl;
            a0 += propensities[i];
        }
        //calculate a-1
        double tau = (1/(double)a0)*log(1/dist(mt));
        double r = a0 * (double) dist(mt);
        cur_time += tau;
        //find reaction index
        int mu = 0;//mu is the selected index of reaction
        double sum = 0.0;
        for(;sum<=r;mu++){
            sum += propensities[mu];
        }
        mu--;
        assert(sum>r);
        cout<<"tau="<<tau<<"; r="<<r<<"; a0="<<a0<<"; mu="<<mu<<endl;
        cout<<"selected reaction";
        reactions[mu].print();
        //adjust count
        reactions[mu].adjustCount(s_count);
        SaveFrame(cur_time);
    }
    SaveResult("output.csv");
}

int CRN::clear() {
    s_count.clear();

}

void CRN::SaveFrame(double cur_time) {
    if(species.empty()){
        for(auto i:s_count){
            species.push_back(i.first);
        }
    }
    Frame f(cur_time);
    for(auto i:s_count){
        f.counts.push_back(i.second);
    }
    simulation_result.push_back(f);
}

void CRN::SaveResult(const string &filename) {
    ofstream outfile;
    cout<<"Saving simulation result to "<<filename<<endl;
    outfile.open(filename);
    outfile<<"time";
    for(string i:species){
        outfile<<","<<i;
    }
    outfile<<"\n";
    for(const Frame& i:simulation_result){
        outfile<<i.time<<",";
        for(int c:i.counts){
            outfile<<c<<",";
        }
        outfile<<"\n";
    }
    outfile.close();

}

int CRN::print() {
    cout<<"CRN contains the following reactions:"<<endl;
    for(auto i:reactions){
       i. print();
    }
    cout<<"CRN contains the following reactions:"<<endl;
    for(auto i:s_count){
        cout<<i.first<<":"<<i.second<<endl;
    }
    cout<<"--------------------------------------"<<endl;

    return 0;
}


Reaction::Reaction(const string &r, const string& p, double l) {
    reactant = parseSpecies(r);
    product = parseSpecies(p);
    cout << r << " is converted to " << endl;
    printSpecies(reactant);

    lambda = l;
}

vector<Species> Reaction::parseSpecies(string a) {
    vector <Species> result;
    vector <string> temp;
    Species t;
    split(temp, a, is_any_of("+"));//may need escape
    for (string s:temp) {
        cout<<"parsing:"<<s<<endl;
        int i = 0;
        for (; i < s.size(); i++) {
            if (!isdigit(s[i]))
                break;
        }
        if (i == s.size()) {
            cout << "Parsing Species Error: Species [" << s << "] format is invalid" << endl;
            exit(EXIT_FAILURE);
        }
        if(i>0)
            result.push_back({ s.substr(i),stoi(s.substr(0,i))});
        else
            result.push_back({ s.substr(i), 1});

    }
    return result;
}

double Reaction::calculateProp(unordered_map<string, int> &m) {
    double result = lambda;
    for(const Species& i:reactant){
        int factor = i.factor;
        int s_count = m[i.name];
        if(s_count<factor) return 0.0; //cannot carry reaction
        result *= pow((double)s_count,(double)factor);
    }
    return result;
}

void Reaction::printSpecies(vector<Species> &v) {
   // cout << "Print Species:----------" << endl;
    for (Species i:v) {
        cout << i.factor << '*' << i.name << ' ';
    }
    //cout << "-------------------------" << endl;
}

void Reaction::print() {
   printSpecies(reactant) ;
   cout<<"->";
   printSpecies(product) ;
   cout<< "rate="<<lambda<<endl;

}

int Reaction::adjustCount(unordered_map<string, int> &m) {
    for(Species r:reactant){
       m[r.name] -= r.factor;
       assert(m[r.name]>=0);
    }
    for(Species r:product){
        m[r.name] += r.factor;
        assert(m[r.name]>=0);
    }
    return 0;
}
