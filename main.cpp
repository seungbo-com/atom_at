#include <iostream>
#include <vector>
#include <fstream>

#include "Atom_tool.h"
#include "gen_tool.h"

using namespace std;

int main() {
    Atom tst("file_name");
    vector <coord_atom> atom_vector;
    double lattice;
    int atom_n;

    tst.atom_position(&atom_vector, &lattice, &atom_n);

    ofstream test_one("new_file_name");

    vector <string> atoms_t = tst.atom_types();
    for (const string&types :atoms_t) {
        for (coord_atom&name :atom_vector){
            if (name.name == types) print_line_Atom(name, 0);
        }
    }

    vector <string> types {"Si","O"};
    string dd = "Zeolite";
    tst.angle(types,dd, test_one);

    test_one.close();
}
