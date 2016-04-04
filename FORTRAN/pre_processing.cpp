#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace std;

// A class to store the boundary conditions
class Set {
public:
    string name;
    vector<int> nodes;
    vector<int> elements;
    vector<int> dofs; // 0, 1, 2 represents x, y, z
    vector<double> values; // corresponding to dofs
};

// A class to store the loads
class Surface {
public:
    string name;
    vector<int> elements;
    int face;
    vector<double> traction;
    double pressure;
};

// Helper function to remove the leading and trailing whitespace from a string
void trim(string& s) {
    size_t p = s.find_first_not_of(" \t\r");
    s.erase(0, p);
    p = s.find_last_not_of(" \t\r");
    if (string::npos != p)
        s.erase(p+1);
}

// Helper function to return the number of nodes on a face
int no_facenodes(int nsd, int nen) {
    int n = 0;
    if (nsd == 2) {
        n = 2;
    } else if (nsd == 3) {
        if (nen == 4) {
            n = 3;
        } else if (nen == 8) {
            n = 4;
        }
    }
    return n;
}

// Helper function to return a vector of the numbers of nodes on a face
vector<int> facenodes(int nsd, int nen, int face) {
    vector<int> list(no_facenodes(nsd,nen),0);
    int i3 [] = {2, 3, 1};
    int i4 [] = {2, 3, 4, 1};
    if (nsd == 2) {
        if (nen == 3) {
            int temp [] = {face, i3[face-1]};
            list.assign(temp, temp + 2);
        } else if (nen == 4) {
            int temp [] = {face, i4[face-1]};
            list.assign(temp, temp + 2);
        }
    } else if (nsd == 3) {
        if (nen == 4) {
            if (face == 1) {
                int temp [] = {1, 2, 3};
                list.assign(temp, temp + 3);
            } else if (face == 2) {
                int temp [] = {1, 4, 2};
                list.assign(temp, temp + 3);
            } else if (face == 3) {
                int temp [] = {2, 4, 3};
                list.assign(temp, temp + 3);
            } else if (face == 4) {
                int temp [] = {3, 4, 1};
                list.assign(temp, temp + 3);
            }
        } else if (nen == 8) {
            if (face == 1) {
            } else if (face == 2) {
                int temp [] = {1, 2, 3, 4};
                list.assign(temp, temp + 4);
            } else if (face == 3) {
                int temp [] = {5, 8, 7, 6};
                list.assign(temp, temp + 4);
            } else if (face == 4) {
                int temp [] = {1, 5, 6, 2};
                list.assign(temp, temp + 4);
            } else if (face == 5) {
                int temp [] = {2, 6, 7, 3};
                list.assign(temp, temp + 4);
            } else if (face == 6) {
                int temp [] = {4, 8, 5, 1};
                list.assign(temp, temp + 4);
            }
        }
    }
    return list;
}

int main (int argc, char* argv[]) {
	string load_type("none");
    if (argc != 2 && argc != 3) {
        cerr << "Invalid number of command line arguments!" << endl;
        return 1;
    }
	if (argc == 3) {
		load_type = argv[2];
		if (load_type != "traction" && load_type != "pressure") {
			cerr << "Invalid type of load!" << endl;
			return 1;
		}
	}
    string filename(argv[1]);
    if (filename.compare(filename.size() - 4, 4, ".inp")) {
        cerr << "Invalid file type, .inp is expected." << endl;
        return 1;
    }
    ifstream instr(filename.c_str());
    if (!instr.good()) {
        cerr << "Cannot open input file!" << endl;
        return 1;
    }

    // Find out the element type
    string buffer;
    int nsd = 0, nen = 0;
    while (getline(instr, buffer)) {
		trim(buffer);
        if (!buffer.compare(0, 15, "*Element, type=")) {
            nsd = buffer[16] - '0';
            nen = buffer[18] - '0';
            break;
        }
    }
    instr.close();
    instr.open(filename.c_str());

    // Head line of coordinates
    while (getline(instr, buffer)) {
		trim(buffer);
        if (!buffer.compare(0, 5, "*Node")) {
            break;
        }
    }
    vector<vector<double> > coords;
    int nn = 0;
    while (getline(instr, buffer)) {
		trim(buffer);
        if (!buffer.compare(0, 8, "*Element")) {
            // Tail line of coordinates
            break;
        } else {
            // Read in the coordinates
            istringstream ss(buffer);
            string token;
            vector<double> one_node(nsd, 0.);
            for (int i = 0; i < 1 + nsd; ++i) {
                getline(ss, token, ',');
                if (i == 0) {
                    nn = atoi(token.c_str());
                } else {
                    one_node[i-1] = atof(token.c_str());
                }
            }
            coords.push_back(one_node);
        }
    }
    assert(nn == coords.size());
    for (int i = 0; i < nn; ++i) assert(nsd == coords[i].size());

    // Head line of connectivity
    vector<vector<int> > connect;
    int nel = 0;
    while (getline(instr, buffer)) {
		trim(buffer);
        if (!buffer.compare(0, 1, "*")) {
            // Tail line of connectivity
            break;
        } else {
            // Read in the connectivity
            istringstream ss(buffer);
            string token;
            vector<int> one_element(nen, 0);
            for (int i = 0; i < 1 + nen; ++i) {
                getline(ss, token, ',');
                if (i == 0) {
                    nel = atoi(token.c_str());
                } else {
                    one_element[i-1] = atof(token.c_str());
                }
            }
            connect.push_back(one_element);
        }
    }
    assert(nel == connect.size());
    for (int i = 0; i < nel; ++i) assert(nen == connect[i].size());
	
// ========================================================================================================================
    vector<Set> sets; // to store the boundary conditions
    while (getline(instr, buffer)) {
		trim(buffer);
        // Head line of sets
        if (!buffer.compare(0, 13, "*End Instance")) break;
    }
    while (!instr.eof()) {
		// Sets
		if (!buffer.compare(0, 16, "*Nset, nset=Set-")) {
            // Read in the sets
            Set new_set;
            string temp;
            ostringstream convert;
            convert << sets.size()+1;
            temp = convert.str();
            new_set.name = "Set-" + temp; 
            // Read in the nodes in a set
            if (!buffer.compare(buffer.size()-8, 8, "generate")) {
                // structured mesh
                getline(instr, buffer);
				trim(buffer);
                istringstream ss(buffer);
                string token;
                int node_info[3];
                for (int i = 0; i < 3; ++i) {
                    getline(ss, token, ',');
                    node_info[i] = atoi(token.c_str());
                }
                for (int i = node_info[0]; i <= node_info[1]; i += node_info[2])
                    new_set.nodes.push_back(i);
                getline(instr, buffer);
                trim(buffer);
            } else {
                // unstructured mesh
                while (getline(instr, buffer)) {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_set.nodes.size() != 0) {
                        break;
                    } 
                    istringstream ss(buffer);
                    string token;
                    while (getline(ss, token, ',')) {
                        new_set.nodes.push_back(atoi(token.c_str()));
                    }
                }
            }
            // Read in the elements in a set
            if (!buffer.compare(buffer.size()-8, 8, "generate")) {
                // structured mesh
                getline(instr, buffer);
				trim(buffer);
                istringstream ss(buffer);
                string token;
                int element_info[3];
                for (int i = 0; i < 3; ++i) {
                    getline(ss, token, ',');
                    element_info[i] = atoi(token.c_str());
                }
                for (int i = element_info[0]; i <= element_info[1]; i += element_info[2])
                    new_set.elements.push_back(i);
            } else {
                // unstructured mesh
                while (getline(instr, buffer)) {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_set.elements.size() != 0) {
                        break;
                    }
                    istringstream ss(buffer);
                    string token;
                    while (getline(ss, token, ',')) {
                        new_set.elements.push_back(atoi(token.c_str()));
                    }
                }
            }
			sets.push_back(new_set);
		} else if (!buffer.compare(0, 20, "*Elset, elset=_Surf-")) {
			// Tail line of sets
			break;
		} else {
			getline(instr, buffer);
			trim(buffer);
		}
	}
    instr.close();
    

// ========================================================================================================================
    // Surfaces
    instr.open(filename.c_str());
    while (getline(instr, buffer)) {
        trim(buffer);
        if (!buffer.compare(0, 13, "*End Instance")) {
            break;
        }
    }
    vector<Surface> surfaces;
    while (!instr.eof()) {
        if (!buffer.compare(0, 20, "*Elset, elset=_Surf-")) {
            Surface new_surface;
            string temp;
            ostringstream convert;
            convert << surfaces.size() + 1;
            temp = convert.str();
            new_surface.name = "Surf-" + temp;
            int count = 0;
            for (int i = 0; i < buffer.size(); ++i) {
                if (buffer[i] == 'i') {
                    count = i;
                    break;
                }
            }
            new_surface.face = buffer[count-3] - '0';
            if (!buffer.compare(buffer.size()-8, 8, "generate")) {
                // Structured mesh
                getline(instr, buffer);
                trim(buffer);
                istringstream ss(buffer);
                string token;
                int element_info[3];
                for (int i = 0; i < 3; ++i) {
                    getline(ss, token, ',');
                    element_info[i] = atoi(token.c_str());
                }
                for (int i = element_info[0]; i <= element_info[1]; i += element_info[2]) {
                    new_surface.elements.push_back(i);
                }
            } else {
                // unstructured mesh
                while (getline(instr, buffer)) {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_surface.elements.size() != 0) {
                        break;
                    }
                    istringstream ss(buffer);
                    string token;
                    while (getline(ss, token, ',')) {
                        new_surface.elements.push_back(atoi(token.c_str()));
                    }
                }
            }
			surfaces.push_back(new_surface); 
        } else if (!buffer.compare(0, 13, "*End Assembly")) {
            break;
        } else {
            getline(instr, buffer);
            trim(buffer);
        }
    }
	instr.close();
	
// ========================================================================================================================
	cout << "Specify boundary conditions, to stop enter stop." << endl;
	cout << "Input name of the set, number of dof, and the value, separated by space." << endl;
	cout << "1, 2, 3, 4 represent ux, uy, uz and p" << endl;
	cout << "For example: Set-1 1 0.0" << endl;
	cout << endl;
	
	while (cin >> buffer) {
		if (buffer == "stop") {
			break;
		}
		string user_dof, user_value;
		cin >> user_dof >> user_value;
		int dof = atoi(user_dof.c_str());
		double value = atof(user_value.c_str());
		for (int i = 0; i < sets.size(); ++i) {
			if (sets[i].name == buffer) {
				sets[i].dofs.push_back(dof);
				sets[i].values.push_back(value);
				break;
			}
		}
	}
	cout << endl;
	
	cout << "Specify loads, you may enter traction, pressure or stop" << endl;
	if (load_type == "traction") {
		cout << "Input name of the surface, the three components of traction, separated by space." << endl;
		cout << "For example, Surf-1, 1., 0., 0." << endl;
		while (cin >> buffer) {
			if (buffer == "stop") {
				break;
			}
			bool found = false;
			for (int i = 0; i < surfaces.size(); ++i) {
				if (surfaces[i].name == buffer) {
					found = true;
					for (int j = 0; j < 3; ++j) {
						cin >> buffer;
						double temp = atof(buffer.c_str());
						surfaces[i].traction.push_back(temp);
					}
					break;
				}
			}
			if (!found) {
				cout << "Set not found!" << endl;
			}
		}
	} else if (load_type == "pressure") {
		cout << "Input name of the surface, the value of pressure, separated by space." << endl;
		cout << "For example, Surf-1, 1000." << endl;
		while (cin >> buffer) {
			if (buffer == "stop") {
				break;
			}
			bool found = false;
			for (int i = 0; i < surfaces.size(); ++i) {
				if (surfaces[i].name == buffer) {
					found = true;
					cin >> buffer;
					surfaces[i].pressure = atof(buffer.c_str());
					break;
				}
			}
			if (!found) {
				cout << "Surface not found!" << endl;
			}
		}
	}
	
	cout << endl << "Summary of boundary conditions and loads:" << endl;
	for (int i = 0; i < sets.size(); ++i) {
		assert(sets[i].dofs.size() == sets[i].values.size());
        cout << sets[i].name << " " << sets[i].nodes.size() << " " << sets[i].elements.size() << " ";
		for (int j = 0; j < sets[i].dofs.size(); ++j) {
			cout << "u" << sets[i].dofs[j] << " = " << sets[i].values[j] << " ";
		}
		cout << endl;
	}
	for (int i = 0; i < surfaces.size(); ++i) {
        cout << surfaces[i].name << "-" << surfaces[i].face << " " << surfaces[i].elements.size() << " ";
		if (load_type == "traction") {
			for (int j = 0; j < surfaces[i].traction.size(); ++j) {
				cout << surfaces[i].traction[j] << " ";
			} 
		} else if (load_type == "pressure") {
			cout << surfaces[i].pressure;
		}
		cout << endl;
	}

// ========================================================================================================================
	// Coords
	ofstream outstr("coords.txt");
	cout << setprecision(8) << fixed;
	outstr << setw(13) << nsd << '\t' << setw(13) << nn << endl;
	for (int i = 0; i < coords.size(); ++i) {
		for (int j = 0; j < nsd; ++j) {
			outstr << setw(13) << coords[i][j] << '\t' ;
		}
		outstr << endl;
	}
	outstr.close();
	
	// Connect
	outstr.open("connect.txt");
	outstr << setw(10) << connect.size() << '\t' << setw(10) << nen << endl;
	for (int i = 0; i < connect.size(); ++i) {
		for (int j = 0; j < nen; ++j) {
			outstr << setw(10) << connect[i][j] << '\t';
		}
		outstr << endl;
	}
	outstr.close();
	
	// bc
	outstr.open("bc.txt");
	int count = 0;
	for (int i = 0; i < sets.size(); ++i) {
		count += sets[i].nodes.size()*sets[i].dofs.size();
	}
	outstr << setw(10) << count << endl;
	for (int i = 0; i < sets.size(); ++i) {
		for (int j = 0; j < sets[i].dofs.size(); ++j) {
			for (int k = 0; k < sets[i].nodes.size(); ++k) {
				outstr << setw(10) << sets[i].nodes[k] << '\t' << setw(10) << sets[i].dofs[j] \
					 << '\t' << setw(13) << setprecision(8) << fixed << sets[i].values[j] << endl;
			}
		}
	}
	outstr.close();
	
	// load
	outstr.open("load.txt");
	count = 0;
	for (int i = 0; i < surfaces.size(); ++i) {
		count += surfaces[i].elements.size();
	}
	if (load_type == "traction") {
		outstr << setw(10) << count << '\t' << setw(10) << nsd + 2 << endl;
		for (int i = 0; i < surfaces.size(); ++i) {
			for (int j = 0; j < surfaces[i].elements.size(); ++j) {
				outstr << setw(10) << surfaces[i].elements[j] << '\t' << setw(10) << surfaces[i].face << '\t';
				for (int k = 0; k < nsd; ++k) {
					outstr << setw(13) << setprecision(10) << scientific << surfaces[i].traction[k] << '\t';
				}
				outstr << endl;		
			}
		}
	} else if (load_type == "pressure") {
		outstr << setw(10) << count << '\t' << setw(10) << 3 << endl;
		for (int i = 0; i < surfaces.size(); ++i) {
			for (int j = 0; j < surfaces[i].elements.size(); ++j) {
				outstr << setw(10) << surfaces[i].elements[j] << '\t' << setw(10) << surfaces[i].face << '\t' \
					<< setw(13) << setprecision(10) << scientific << surfaces[i].pressure << endl;
			}
		}
	}
	outstr.close();
	
	
    return 0;
}
