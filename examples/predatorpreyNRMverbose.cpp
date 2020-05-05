//
// Created by willf on 4/21/2020.
//
#include <iostream>
#include "../crn.h"

int main() {
    int n=300;
    int tmax = 50;
    CRN  rsys;
    rsys.addRxn("x1+x2","2x2",1.0/n);
    rsys.addRxn("x1","2x1",1.0);
    rsys.addRxn("x2","",1.0);
    rsys.setConc("x1",n/2);
    rsys.setConc("x2",n/2);
    rsys.print();
    //tmax==0 will continue until no further reaction occur
    rsys.simulate(tmax, true, "NRM");
    rsys.plot("/build/predatorpreyNRM");
    rsys.saveResult("output.csv");
}

