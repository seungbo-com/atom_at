#include <string>
#include <vector>
#include <cmath>

#include "Atom_tool.h"

using namespace std;


void pbc_dist_angle(coord_atom&one, coord_atom&two,double length, double* adj_dist, vector <double>* adj_vector){

    vector <string> one_vec = {one.x, one.y, one.z};
    vector <string> two_vec = {two.x, two.y, two.z};
    vector <double> adj_vv{};
    double tot = 0;
    for (int i  = 0; i < 3; i++){
        double diff = stod(one_vec[i]) - stod(two_vec[i]);
        tot += pow(diff - ( length * round(diff / length)), 2);
        adj_vv.push_back(diff - ( length * round(diff / length)));
    }
    *adj_dist = sqrt(tot);
    *adj_vector = adj_vv;
}


void print_line_Atom(coord_atom& sin_line, int last = 1){
    cout << sin_line.name << " "<< sin_line.x <<" " <<sin_line.y << " "<<sin_line.z << " ";
    if (last) cout << endl;
    else cout << sin_line.position << endl;
}


double angle_vec(value_hold& a, value_hold& b){
    vector <double> a_vec, b_vec;
    a_vec = a.vec;
    b_vec = b.vec;
    double tot_sum, a_sum, b_sum, angle_cos;
    tot_sum= a_sum = b_sum = 0.0;
    for (int i = 0;i < 3; i++) {
        tot_sum += a_vec[i] * b_vec[i];
        a_sum += a_vec[i] * a_vec[i];
        b_sum += b_vec[i] * b_vec[i];
    }
    double cosi = tot_sum / (sqrt(a_sum) * sqrt(b_sum));
    angle_cos = acos(cosi) * 180/M_PI;
    return angle_cos;
}



// zeloite case* : Checking the bond atoms within the vector
bool zeolite_connec(value_hold& a, value_hold& b, vector<coord_atom>& all_coord, double &lattice){
    coord_atom a_second, b_second;
    for (coord_atom& ea_coord: all_coord){
        if (a.two == ea_coord.position) a_second = ea_coord;
        if (b.two == ea_coord.position) b_second = ea_coord;
    }
    double adj_dist;
    vector <double> adj_vector;
    pbc_dist_angle(a_second, b_second,lattice,&adj_dist, &adj_vector);
    //the O or H with Si
    if (a.one < a.two && a.one + 5 > a.two && b.one < b.two && b.one + 5 > b.two && adj_dist < 0.1) return true;
    else return false;
}

// zeloite case* : inputting the type confirm it is specific vector
bool zeolite_comb(value_hold& a, vector<coord_atom>& b, vector <string> &comb){
    string type_a, type_b;
    for (coord_atom &check_b : b){
        if (a.one == check_b.position) type_a = check_b.name;
        else if (a.two == check_b.position) type_b = check_b.name;
    }
    if (type_a == comb[0] && type_b == comb[1]) return true;
    else return false;
}