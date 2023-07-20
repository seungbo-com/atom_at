#include <iostream>
#include <string>
#include <vector>

#pragma once

using namespace std;

struct coord_atom{
    string name, x, y, z;
    int position;
};

struct value_hold{
    int one,two;
    double value = 0.0;
    vector <double> vec{};
};


class Atom{
    string file_name; //XYZ file
public:
    Atom(string new_file);
    void atom_position(vector <coord_atom>* atoms, double* lattice_para, int* atom_nu);
    vector <string> atom_types();
    void atom_phys_val(vector <value_hold>* all_dist_atom,  vector <value_hold>* all_vec_atom,vector <string> type);
    void angle(vector <string>& types, string &task, ofstream &file);
};
