//
// Created by willf on 4/21/2020.
//
#include <iostream>
#include "crn.h"

int main() {
    int n=1000;
    CRN  rsys;
    rsys.addRxn("x+y","x+b",1.0/n);
    rsys.addRxn("x+y","y+b",1.0/n);
    rsys.addRxn("x+b","2x",1.0/n);
    rsys.addRxn("b+y","2y",1.0/n);
    rsys.setConc("x",n/2);
    rsys.setConc("y",n/2);
    rsys.setConc("b",0);
    rsys.print();
    //tmax==0 will continue until no further reaction occur
    rsys.simulate(0, false, "NRM");
    rsys.saveResult("output.csv");
}

