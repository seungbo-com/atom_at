#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "gen_tool.h"

using namespace std;


Atom::Atom(string new_file) : file_name(new_file)  {}

void Atom::atom_position(vector <coord_atom>* atoms, double* lattice_para, int* atom_nu) { //Position of Atom

    ifstream inFile(file_name);
    if (!inFile) {
        cerr << "Unable to Open File" << endl;
        exit(1);
    }

    string line;
    int natom = 0;
    double lattice = 0.0;

    if (getline(inFile, line)){
        istringstream iss(line);
        iss >> natom; // # of atoms
    }
    *atom_nu = natom;

    if (getline(inFile, line)){
        istringstream iss(line);
        iss >> lattice;
    }
    *lattice_para = lattice;


    vector <coord_atom> atoms_pre;
    int posit = 1;
    //The atom line
    while (getline(inFile, line)){
        istringstream iss(line);
        coord_atom atom;
        if (!(iss >> atom.name >> atom.x >> atom.y >> atom.z)){
            break;
        }
        atom.position = posit;
        atoms_pre.push_back(atom);
        posit++;
    }
    *atoms = atoms_pre;
}


// Atom kind
vector <string> Atom::atom_types(){
    vector <string> atoms_t;

    vector <coord_atom> atom_vector;
    double lattice;
    int atom_n;

    Atom::atom_position(&atom_vector, &lattice, &atom_n);
    for (const coord_atom&ea_pos: atom_vector){
        auto it = find(atoms_t.begin(), atoms_t.end(),ea_pos.name);
        if (it == atoms_t.end()){
            atoms_t.push_back(ea_pos.name);
        }
    }
    return atoms_t;
}



// minimal image convention
void Atom::atom_phys_val(vector <value_hold>* all_dist_atom,  vector <value_hold>* all_vec_atom, vector <string> type){

    vector <coord_atom> atom_vector; // the vector
    double lattice;
    int atom_n;
    Atom::atom_position(&atom_vector, &lattice, &atom_n);

    //Checking the atom existence
    vector <string> options_atom = Atom::atom_types();
    for (string& ea_atom : type) {
        auto i = find(options_atom.begin(),options_atom.end(), ea_atom);
        if (i == options_atom.end()) {
            cerr << "The atom doesn't exists"<<endl;
            exit(1);
        }
    }

    vector <value_hold> atom_final_dist{}, atom_final_vector{};
    for (coord_atom&atom_one :atom_vector){
        for (coord_atom&atom_two :atom_vector){
            // checking the "Si_BO"
            if (atom_one.name == type[0] && atom_two.name == type[1] && atom_one.position < atom_two.position && atom_one.position + 5 > atom_two.position){
                value_hold dist{}, vector_ea{};
                double adj_dist;
                vector <double> adj_vector;
                pbc_dist_angle(atom_one, atom_two,lattice, &adj_dist, &adj_vector);
                dist.one = vector_ea.one = atom_one.position;
                dist.two = vector_ea.two = atom_two.position;
                dist.value = adj_dist;
                vector_ea.vec = adj_vector;
                atom_final_dist.push_back(dist);
                atom_final_vector.push_back(vector_ea);
            }
        }
    }
    *all_dist_atom = atom_final_dist;
    *all_vec_atom = atom_final_vector;
}

void Atom::angle(vector <string>& types, string &task, ofstream &file){
    vector <double> total_ang{};

    vector <value_hold> dist_vec, vector_vec;
    Atom::atom_phys_val(&dist_vec, &vector_vec,types); //the vector gathering

    vector <coord_atom> atom_vector;
    double lattice;
    int atom_n;
    Atom::atom_position(&atom_vector, &lattice, &atom_n); // the position in the file



    if (task == "Zeolite"){
        for (value_hold& ea_vec: vector_vec){
            for (value_hold& eb_vec: vector_vec){
                if (zeolite_connec(ea_vec, eb_vec, atom_vector,lattice)  && ea_vec.one < eb_vec.one){ // SiO(1) != SiO(2)
                    if (zeolite_comb(ea_vec, atom_vector,types) && zeolite_comb(eb_vec, atom_vector,types)){
                        file << ea_vec.one << " " << ea_vec.two << " " << eb_vec.two << " "<< eb_vec.one <<endl;
                        file << angle_vec(ea_vec,eb_vec)<<"\n\n" << endl;

                        total_ang.push_back(angle_vec(ea_vec,eb_vec));
                    }
                }
            }}}
}