#include "Atom_tool.h"

#pragma once

void pbc_dist_angle(coord_atom&one, coord_atom&two,double length, double* adj_dist, vector <double>* adj_vector);

void print_line_Atom(coord_atom& sin_line, int last = 1);

double angle_vec(value_hold& a, value_hold& b);

bool zeolite_connec(value_hold& a, value_hold& b, vector<coord_atom>& all_coord, double &lattice);

bool zeolite_comb(value_hold& a, vector<coord_atom>& b, vector <string> &comb);

