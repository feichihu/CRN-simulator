//
// Created by willf on 4/6/2020.
//

#include <iostream>
#include <string>
#include <limits>
#include "crn.h"

using namespace std;
#define DEFAULT_MAX_TIME 50000.0 //max time when tmax not specified

CRN::CRN() {

}

int CRN::addRxn(const string &reactant, const string &product, double lambda) {
    //x-1+x2, 2x2, 1.5
    reactions.emplace_back(reactant, product, lambda, (int)reactions.size());

}

int CRN::setConc(string species_name, int init_count) {
    s_count[species_name] = init_count;
}

int CRN::simulate(int tmax, bool verbose, string mode) {
    if(mode!="DM" && mode!="NRM"){
        cout<<"error: method "<<mode<<" is not implemented, please choose from DM or NRM"<<endl;
        exit(0);
    }
    cout<<"-------start simulation--------"<<endl;
    //C++11 random number generator (0,1)
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    //timing
    auto start = high_resolution_clock::now();
    //dist(mt) will be in (0,1)
    if(mode=="DM"){
        double max_time = tmax>0? (double) tmax : DEFAULT_MAX_TIME;
        double cur_time = 0;
        int num_reactions = reactions.size();
        vector<double> propensities(num_reactions);
        saveFrame(cur_time);
        while(cur_time < max_time){
            if(verbose){
                cout<<"++++++++++at "<<cur_time<<"+++++++++"<<endl;
            }
            double a0{};
            //propensity for each reaction
            for(int i=0;i<num_reactions;i++){
                propensities[i] = reactions[i].calculateProp(s_count);
                if(verbose){
                    reactions[i].print();
                    cout<<propensities[i]<<endl;
                }
                a0 += propensities[i];
            }
            if(tmax<=0 && a0==0.0){
                cout<<"No further reaction occurs, simulation terminates."<<endl;
                cout<<"total time elapsed:"<<cur_time<<endl;
                break;
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
            if(verbose){
                cout<<"tau="<<tau<<"; r="<<r<<"; a0="<<a0<<"; mu="<<mu<<endl;
                cout<<"selected reaction";
                reactions[mu].print();
            }
            if(max_time<cur_time){
                cout<<"Simulation terminates."<<endl;
                cout<<"total time elapsed:"<<cur_time<<endl;
            }
            //adjust count
            reactions[mu].adjustCount(s_count);
            saveFrame(cur_time);
        }

    }
    else if (mode=="NRM"){
        //init
        double inf = numeric_limits<double>::infinity();
        generateDependency();
        vector<Slot> slots;
        vector<double> propensities(reactions.size());
        vector<double> taus(reactions.size());
        for(int i=0;i<reactions.size();i++){
           double alpha = reactions[i].calculateProp(s_count);
           propensities[i] = alpha;
           double tau = 1.0/alpha*log(1/dist(mt));
           taus[i] = tau;
           slots.emplace_back(i,tau);
        }
        PriorityQ pq(slots);
        double cur_time = 0.0;
        saveFrame(cur_time);
        double max_time = tmax>0? (double) tmax : DEFAULT_MAX_TIME;
        while(cur_time<max_time){
           Slot cur_slot =  pq.pop();
           if(cur_slot.tau==inf){
               cout<<"No further reaction occurs, simulation terminates."<<endl;
               cout<<"total time elapsed:"<<cur_time<<endl;
               break;
           }
            if(verbose){
                cout<<"++++++++++at "<<cur_time<<"+++++++++"<<endl;
                cout<<"selected reaction"<<cur_slot.reaction_idx<<" ";
                reactions[cur_slot.reaction_idx].print();
            }
           cur_time = cur_slot.tau;
           reactions[cur_slot.reaction_idx].adjustCount(s_count);
           for(int i:dependencyList[cur_slot.reaction_idx]){
               //for each reaction idx in dependency
               Slot temp;
               if(verbose){
                   cout<<"dependency graph from "<<cur_slot.reaction_idx<<" to "<<i<<endl;
               }
               double alpha = reactions[i].calculateProp(s_count);
               double tau_alpha;
               if(cur_slot.reaction_idx==i){
                  tau_alpha = 1.0/ alpha*log(1.0/dist(mt))+cur_time;
               }
               else{
                   //cout<<"old alpha="<<propensities[i]<<", new alpha="<<alpha<<", tau="<<taus[i]<<endl;
                   if(propensities[i]==0){
                       tau_alpha = 1.0/ alpha*log(1.0/dist(mt))+cur_time;
                   }
                   else{
                       tau_alpha = propensities[i]/alpha *(taus[i]-cur_time)+cur_time;
                   }
                   //cout<<"tau_alpha="<<tau_alpha<<endl;
               }
               propensities[i] = alpha;
               taus[i] = tau_alpha;
               pq.push(Slot(i,tau_alpha));
           }
           saveFrame(cur_time);
        }


    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-start);
    cout<<"simulation duration(computation time)="<<duration.count()<<"ms"<<endl;
}

int CRN::clear() {
    s_count.clear();

}

int CRN::plot(const string name) {
    if(simulation_result.empty()){
        cout<<"No simulation result found, run simulate() first!"<<endl;
        return 0;
    }
    vector<double> time;
    int num_spec = s_count.size();
    vector<vector<double>> counts(num_spec);
    for(Frame &f:simulation_result){
       time.push_back(f.time);
       for(int i=0;i<num_spec;i++){
           counts[i].push_back((double)f.counts[i]);
       }
    }
    //plt::figure_size(1200, 780);
    for(int i=0;i<num_spec;i++){
        // Set the size of output image to 1200x780 pixels
        // Plot a red dashed line from given x and y data.
        //plt::plot(x, counts[i],"r--");
        // Plot a line whose name will show up as "log(x)" in the legend.
        plt::named_plot(species[i], time, counts[i]);
    }
    // Add graph title
    plt::title(name);
    // Enable legend.
    plt::legend();
    // Save the image (file format is determined by the extension)
    plt::save("./"+name+".png");
    plt::show();

}

void CRN::saveFrame(double cur_time) {
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

void CRN::saveResult(const string &filename) {
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

void CRN::generateDependency() {
    unordered_map<string,vector<int>> pool;
    for(auto r:reactions){
        for(auto react:r.reactant){
            pool[react.name].push_back(r.idx);
        }
    }
    dependencyList.clear();
    for(auto r:reactions){
        unordered_set<int> temp;
        for(auto p:r.product){
            if(pool.count(p.name)){
                temp.insert(pool[p.name].begin(),pool[p.name].end());
            }
        }
        vector<int> cur_list;
        for(int i:temp){
            cur_list.push_back(i);
        }
        dependencyList.push_back(cur_list);
    }
    for(int i=0;i<dependencyList.size();i++){
        cout<<i<<":[";
        for(int j:dependencyList[i]){
            cout<<j<<' ';
        }
        cout<<"]"<<endl;
    }
}


Reaction::Reaction(const string &r, const string& p, double l,int id) {
    reactant = parseSpecies(r);
    product = parseSpecies(p);
    //cout << r << " is converted to " << endl;
    //printSpecies(reactant);
    idx = id;
    lambda = l;
}

vector<Species> Reaction::parseSpecies(string a) {
    vector <Species> result;
    if(a.empty()) return result;
    vector <string> temp;
    Species t;
    split(temp, a, is_any_of("+"));//may need escape
    for (string s:temp) {
        //cout<<"parsing:"<<s<<endl;
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
            result.emplace_back( s.substr(i),stoi(s.substr(0,i)));
        else
            result.emplace_back( s.substr(i), 1);

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

void PriorityQ::push(Slot s) {
    cur_seq[s.reaction_idx]++;
    s.seq = cur_seq[s.reaction_idx];
    clean();
    pq.push(s);
    /*
    cout<<"pushing element"<<endl;
    print();
    cout<<"pq size="<<pq.size()<<endl;
     */
    //clean();
}

Slot PriorityQ::pop() {
    clean();
    Slot i = pq.top();
    /*
    cout<<"poping element"<<endl;
    print();
     */
    pq.pop();
    return i;
}

Slot PriorityQ::top() {
    clean();
    return pq.top();
}

void PriorityQ::print() {
    clean();
    Slot i = pq.top();
    cout<<"PQ top element:[idx="<<i.reaction_idx<<", tau="<<i.tau<<", seq="<<i.seq<<endl;
}

PriorityQ::PriorityQ() {

}
PriorityQ::PriorityQ(vector<Slot>& v) {
   pq = priority_queue<Slot,vector<Slot>, comparator>(v.begin(),v.end());
   cur_seq = vector<int>(pq.size(),0);
}

void PriorityQ::clean() {
    while(!verifyTop()){
    }
}

inline bool PriorityQ::verifyTop() {
    const Slot &s = pq.top();
    assert(pq.size()+1>=cur_seq.size());
    if(cur_seq[s.reaction_idx]>s.seq){
        pq.pop();
        return false;
    }
    else{
        return true;
    }
}
