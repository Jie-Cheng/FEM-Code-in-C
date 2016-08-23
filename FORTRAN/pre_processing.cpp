/*
    TODO:
    1. replace assert with exceptions
*/
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <set>
#include <exception>

// A class to store the boundary conditions
class Set 
{
public:
    std::string name;
    std::vector<int> nodes;
    std::vector<int> elements;
    std::vector<int> dofs; // 0, 1, 2 represents x, y, z
    std::vector<double> values; // corresponding to dofs
};

// A class to store the loads
class Surface 
{
public:
    std::string name;
    std::vector<int> elements;
    int face;
    std::vector<double> traction;
    double pressure;
};

// A class to store the nonzero pattern
class Entry 
{
public:
    Entry(int rw, int cl) : row(rw), col(cl) {};
    int row;
    int col;
    bool operator<(const Entry& r) const 
    {
        return ((row < r.row)|| (row == r.row && col < r.col));
    }
    bool operator==(const Entry& r) const 
    {
        return (row == r.row && col == r.col);
    }
    bool operator!=(const Entry& r) const
    {
        return !(*this == r);
    }
};

int nsd = 0;
int nen = 0;
int nn = 0;
int nel = 0;
std::vector<std::vector<double> > coords;
std::vector<std::vector<int> > connect;
std::string load_type("none");
std::vector<Set> sets;
std::vector<Surface> surfaces;

void trim(std::string& s);
int no_facenodes(int nsd, int nen);
std::vector<int> facenodes(int nsd, int nen, int face);
std::string getProcessInfo(int argc, char* argv[]);
void getCoordsAndConnect(const std::string& filename);
void getSetInfo(const std::string& filename);
void getSurfaceInfo(const std::string& filename);
void setBcAndLoads();
void writeFiles();
void writeCRS_v1(); // The first way to write CRS
void writeCRS_v2(); // The second way to write CRS

int main (int argc, char* argv[]) {
    std::string filename;
    std::string buffer;
    try
    {
        filename = getProcessInfo(argc, argv);
    }
    catch (char const* message)
    {
        std::cout << message << std::endl;
        return 1;
    }
    getCoordsAndConnect(filename);
    getSetInfo(filename);
    getSurfaceInfo(filename);
    setBcAndLoads();
    writeFiles();
    writeCRS_v1();
    writeCRS_v2();
    return 0;
}

// Helper function to remove the leading and trailing whitespace from a string
void trim(std::string& s) 
{
    size_t p = s.find_first_not_of(" \t\r");
    s.erase(0, p);
    p = s.find_last_not_of(" \t\r");
    if (std::string::npos != p)
        s.erase(p+1);
}

// Helper function to return the number of nodes on a face
int no_facenodes(int nsd, int nen) 
{
    int n = 0;
    if (nsd == 2) 
    {
        n = 2;
    } 
    else if (nsd == 3) 
    {
        if (nen == 4) 
        {
            n = 3;
        } 
        else if (nen == 8) 
        {
            n = 4;
        }
    }
    return n;
}

// Helper function to return a vector of the numbers of nodes on a face
std::vector<int> facenodes(int nsd, int nen, int face) 
{
    std::vector<int> list(no_facenodes(nsd,nen),0);
    int i3 [] = {2, 3, 1};
    int i4 [] = {2, 3, 4, 1};
    if (nsd == 2) 
    {
        if (nen == 3) 
        {
            int temp [] = {face, i3[face-1]};
            list.assign(temp, temp + 2);
        } 
        else if (nen == 4) 
        {
            int temp [] = {face, i4[face-1]};
            list.assign(temp, temp + 2);
        }
    } 
    else if (nsd == 3) 
    {
        if (nen == 4) 
        {
            if (face == 1) 
            {
                int temp [] = {1, 2, 3};
                list.assign(temp, temp + 3);
            } 
            else if (face == 2) 
            {
                int temp [] = {1, 4, 2};
                list.assign(temp, temp + 3);
            } 
            else if (face == 3) 
            {
                int temp [] = {2, 4, 3};
                list.assign(temp, temp + 3);
            } 
            else if (face == 4) 
            {
                int temp [] = {3, 4, 1};
                list.assign(temp, temp + 3);
            }
        } 
        else if (nen == 8) 
        {
            if (face == 2) 
            {
                int temp [] = {1, 2, 3, 4};
                list.assign(temp, temp + 4);
            } 
            else if (face == 3) 
            {
                int temp [] = {5, 8, 7, 6};
                list.assign(temp, temp + 4);
            } 
            else if (face == 4)
            {
                int temp [] = {1, 5, 6, 2};
                list.assign(temp, temp + 4);
            } 
            else if (face == 5) 
            {
                int temp [] = {2, 6, 7, 3};
                list.assign(temp, temp + 4);
            } 
            else if (face == 6) 
            {
                int temp [] = {4, 8, 5, 1};
                list.assign(temp, temp + 4);
            }
        }
    }
    return list;
}

std::string getProcessInfo(int argc, char* argv[])
{
    if (argc != 2 && argc != 3) 
    {
        throw std::string("Invalid number of command line arguments!");
    }
    if (argc == 3) 
    {
        load_type = argv[2];
        if (load_type != "traction" && load_type != "pressure") 
        {
            throw std::string("Invalid type of load!");
        }
    }
    std::string filename(argv[1]);
    if (filename.compare(filename.size() - 4, 4, ".inp")) 
    {
        throw std::string("Invalid file type, .inp is expected!");
    }
    std::ifstream instr(filename.c_str());
    if (!instr.good()) 
    {
        throw std::string("Cannot open input file!");
    }
    // Find out the element type
    std::string buffer;
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        if (!buffer.compare(0, 15, "*Element, type=")) 
        {
            nsd = buffer[16] - '0';
            nen = buffer[18] - '0';
            break;
        }
    }
    instr.close();
    return filename;
}

void getCoordsAndConnect(const std::string& filename)
{
    std::ifstream instr(filename.c_str());
    std::string buffer;
    // Head line of coordinates
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        if (!buffer.compare(0, 5, "*Node")) 
        {
            break;
        }
    }
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        if (!buffer.compare(0, 8, "*Element")) 
        {
            // Tail line of coordinates
            break;
        } 
        else 
        {
            // Read in the coordinates
            std::istringstream ss(buffer);
            std::string token;
            std::vector<double> one_node(nsd, 0.);
            for (int i = 0; i < 1 + nsd; ++i) 
            {
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
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        if (!buffer.compare(0, 1, "*")) 
        {
            // Tail line of connectivity
            break;
        } 
        else 
        {
            // Read in the connectivity
            std::istringstream ss(buffer);
            std::string token;
            std::vector<int> one_element(nen, 0);
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
    instr.close();
}

void getSetInfo(const std::string& filename)
{
    std::string buffer;
    std::ifstream instr(filename.c_str());
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        // Head line of sets
        if (!buffer.compare(0, 13, "*End Instance")) break;
    }
    while (!instr.eof()) 
    {
        // Sets
        if (!buffer.compare(0, 16, "*Nset, nset=Set-")) 
        {
            // Read in the sets
            Set new_set;
            std::string temp;
            std::ostringstream convert;
            convert << sets.size()+1;
            temp = convert.str();
            new_set.name = "Set-" + temp; 
            // Read in the nodes in a set
            if (!buffer.compare(buffer.size()-8, 8, "generate")) 
            {
                // structured mesh
                getline(instr, buffer);
                trim(buffer);
                std::istringstream ss(buffer);
                std::string token;
                int node_info[3];
                for (int i = 0; i < 3; ++i) 
                {
                    getline(ss, token, ',');
                    node_info[i] = atoi(token.c_str());
                }
                for (int i = node_info[0]; i <= node_info[1]; i += node_info[2])
                    new_set.nodes.push_back(i);
                getline(instr, buffer);
                trim(buffer);
            } 
            else 
            {
                // unstructured mesh
                while (getline(instr, buffer)) 
                {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_set.nodes.size() != 0) 
                    {
                        break;
                    } 
                    std::istringstream ss(buffer);
                    std::string token;
                    while (getline(ss, token, ',')) 
                    {
                        new_set.nodes.push_back(atoi(token.c_str()));
                    }
                }
            }
            // Read in the elements in a set
            if (!buffer.compare(buffer.size()-8, 8, "generate")) 
            {
                // structured mesh
                getline(instr, buffer);
                trim(buffer);
                std::istringstream ss(buffer);
                std::string token;
                int element_info[3];
                for (int i = 0; i < 3; ++i) 
                {
                    getline(ss, token, ',');
                    element_info[i] = atoi(token.c_str());
                }
                for (int i = element_info[0]; i <= element_info[1]; i += element_info[2])
                    new_set.elements.push_back(i);
            } 
            else 
            {
                // unstructured mesh
                while (getline(instr, buffer)) 
                {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_set.elements.size() != 0) 
                    {
                        break;
                    }
                    std::istringstream ss(buffer);
                    std::string token;
                    while (getline(ss, token, ',')) 
                    {
                        new_set.elements.push_back(atoi(token.c_str()));
                    }
                }
            }
            sets.push_back(new_set);
        } 
        else if (!buffer.compare(0, 20, "*Elset, elset=_Surf-")) 
        {
            // Tail line of sets
            break;
        } 
        else 
        {
            getline(instr, buffer);
            trim(buffer);
        }
    }
    instr.close();
}

void getSurfaceInfo(const std::string& filename)
{
    std::string buffer;
    std::ifstream instr(filename.c_str());
    while (getline(instr, buffer)) 
    {
        trim(buffer);
        if (!buffer.compare(0, 13, "*End Instance")) 
        {
            break;
        }
    }
    while (!instr.eof()) 
    {
        if (!buffer.compare(0, 20, "*Elset, elset=_Surf-")) 
        {
            Surface new_surface;
            std::string temp;
            std::ostringstream convert;
            convert << surfaces.size() + 1;
            temp = convert.str();
            new_surface.name = "Surf-" + temp;
            int count = 0;
            for (int i = 0; i < buffer.size(); ++i) 
            {
                if (buffer[i] == 'i') 
                {
                    count = i;
                    break;
                }
            }
            new_surface.face = buffer[count-3] - '0';
            if (!buffer.compare(buffer.size()-8, 8, "generate")) 
            {
                // Structured mesh
                getline(instr, buffer);
                trim(buffer);
                std::istringstream ss(buffer);
                std::string token;
                int element_info[3];
                for (int i = 0; i < 3; ++i) 
                {
                    getline(ss, token, ',');
                    element_info[i] = atoi(token.c_str());
                }
                for (int i = element_info[0]; i <= element_info[1]; i += element_info[2]) 
                {
                    new_surface.elements.push_back(i);
                }
            } 
            else 
            {
                // unstructured mesh
                while (getline(instr, buffer)) 
                {
                    trim(buffer);
                    if (!buffer.compare(0, 1, "*") && new_surface.elements.size() != 0) 
                    {
                        break;
                    }
                    std::istringstream ss(buffer);
                    std::string token;
                    while (getline(ss, token, ',')) 
                    {
                        new_surface.elements.push_back(atoi(token.c_str()));
                    }
                }
            }
            surfaces.push_back(new_surface); 
        } 
        else if (!buffer.compare(0, 13, "*End Assembly")) 
        {
            break;
        } 
        else 
        {
            getline(instr, buffer);
            trim(buffer);
        }
    }
    instr.close();
}

void setBcAndLoads()
{
    std::cout << "Specify boundary conditions, enter 'stop' to stop." << std::endl;
    std::cout << "Input name of the set, number of dof, and the value, separated by space." << std::endl;
    std::cout << "1, 2, 3, 4 represent ux, uy, uz and p" << std::endl;
    std::cout << "For example: Set-1 1 0.0" << std::endl;
    std::cout << std::endl;
    
    std::string buffer;
    while (std::cin >> buffer) 
    {
        if (buffer == "stop") 
        {
            break;
        }
        std::string user_dof, user_value;
        std::cin >> user_dof >> user_value;
        int dof = atoi(user_dof.c_str());
        double value = atof(user_value.c_str());
        for (int i = 0; i < sets.size(); ++i) 
        {
            if (sets[i].name == buffer) 
            {
                sets[i].dofs.push_back(dof);
                sets[i].values.push_back(value);
                break;
            }
        }
    }
    std::cout << std::endl;
    
    std::cout << "Specify loads, you may enter 'traction', 'pressure' or 'stop'" << std::endl;
    if (load_type == "traction") 
    {
        std::cout << "Input name of the surface, the three components of traction, separated by space." << std::endl;
        std::cout << "For example, Surf-1, 1., 0., 0." << std::endl;
        while (std::cin >> buffer) 
        {
            if (buffer == "stop") 
            {
                break;
            }
            bool found = false;
            for (int i = 0; i < surfaces.size(); ++i) 
            {
                if (surfaces[i].name == buffer) 
                {
                    found = true;
                    for (int j = 0; j < 3; ++j) 
                    {
                        std::cin >> buffer;
                        double temp = atof(buffer.c_str());
                        surfaces[i].traction.push_back(temp);
                    }
                    break;
                }
            }
            if (!found) 
            {
                std::cout << "Set not found!" << std::endl;
            }
        }
    } 
    else if (load_type == "pressure") 
    {
        std::cout << "Input name of the surface, the value of pressure, separated by space." << std::endl;
        std::cout << "For example, Surf-1, 1000." << std::endl;
        while (std::cin >> buffer) 
        {
            if (buffer == "stop") 
            {
                break;
            }
            bool found = false;
            for (int i = 0; i < surfaces.size(); ++i) 
            {
                if (surfaces[i].name == buffer) 
                {
                    found = true;
                    std::cin >> buffer;
                    surfaces[i].pressure = atof(buffer.c_str());
                    break;
                }
            }
            if (!found) 
            {
                std::cout << "Surface not found!" << std::endl;
            }
        }
    }
    
    std::cout << std::endl << "Summary of boundary conditions and loads:" << std::endl;
    for (int i = 0; i < sets.size(); ++i) 
    {
        assert(sets[i].dofs.size() == sets[i].values.size());
        std::cout << sets[i].name << " " << sets[i].nodes.size() << " " << sets[i].elements.size() << " ";
        for (int j = 0; j < sets[i].dofs.size(); ++j) 
        {
            std::cout << "u" << sets[i].dofs[j] << " = " << sets[i].values[j] << " ";
        }
        std::cout << std::endl;
    }
    for (int i = 0; i < surfaces.size(); ++i) 
    {
        std::cout << surfaces[i].name << "-" << surfaces[i].face << " " << surfaces[i].elements.size() << " ";
        if (load_type == "traction") 
        {
            for (int j = 0; j < surfaces[i].traction.size(); ++j) 
            {
                std::cout << surfaces[i].traction[j] << " ";
            } 
        } 
        else if (load_type == "pressure") 
        {
            std::cout << surfaces[i].pressure;
        }
        std::cout << std::endl;
    }
}

void writeFiles()
{
    // Coords
    std::ofstream outstr("coords.txt");
    std::cout << std::setprecision(8) << std::fixed;
    outstr << std::setw(13) << nsd << '\t' << std::setw(13) << nn << std::endl;
    for (int i = 0; i < coords.size(); ++i) 
    {
        for (int j = 0; j < nsd; ++j) 
        {
            outstr << std::setw(13) << coords[i][j] << '\t' ;
        }
        outstr << std::endl;
    }
    outstr.close();
    
    // Connect
    outstr.open("connect.txt");
    outstr << std::setw(10) << connect.size() << '\t' << std::setw(10) << nen << std::endl;
    for (int i = 0; i < connect.size(); ++i) 
    {
        for (int j = 0; j < nen; ++j) 
        {
            outstr << std::setw(10) << connect[i][j] << '\t';
        }
        outstr << std::endl;
    }
    outstr.close();
    
    // bc
    outstr.open("bc.txt");
    int count = 0;
    for (int i = 0; i < sets.size(); ++i) 
    {
        count += sets[i].nodes.size()*sets[i].dofs.size();
    }
    outstr << std::setw(10) << count << std::endl;
    for (int i = 0; i < sets.size(); ++i) 
    {
        for (int j = 0; j < sets[i].dofs.size(); ++j) 
        {
            for (int k = 0; k < sets[i].nodes.size(); ++k) 
            {
                outstr << std::setw(10) << sets[i].nodes[k] << '\t' << std::setw(10) << sets[i].dofs[j] \
                     << '\t' << std::setw(13) << std::setprecision(8) << std::fixed << sets[i].values[j] << std::endl;
            }
        }
    }
    outstr.close();
    
    // load
    outstr.open("load.txt");
    count = 0;
    for (int i = 0; i < surfaces.size(); ++i) 
    {
        count += surfaces[i].elements.size();
    }
    if (load_type == "traction") 
    {
        outstr << std::setw(10) << count << '\t' << std::setw(10) << nsd + 2 << std::endl;
        for (int i = 0; i < surfaces.size(); ++i) {
            for (int j = 0; j < surfaces[i].elements.size(); ++j) 
            {
                outstr << std::setw(10) << surfaces[i].elements[j] << '\t' << std::setw(10) << surfaces[i].face << '\t';
                for (int k = 0; k < nsd; ++k) {
                    outstr << std::setw(13) << std::setprecision(10) << std::scientific << surfaces[i].traction[k] << '\t';
                }
                outstr << std::endl;     
            }
        }
    } 
    else if (load_type == "pressure") 
    {
        outstr << std::setw(10) << count << '\t' << std::setw(10) << 3 << std::endl;
        for (int i = 0; i < surfaces.size(); ++i) 
        {
            for (int j = 0; j < surfaces[i].elements.size(); ++j) 
            {
                outstr << std::setw(10) << surfaces[i].elements[j] << '\t' << std::setw(10) << surfaces[i].face << '\t' \
                    << std::setw(13) << std::setprecision(10) << std::scientific << surfaces[i].pressure << std::endl;
            }
        }
    }
    outstr.close();
}

// First version to write CRS
void writeCRS_v1()
{
    int no_nonzeros = 0;
    std::vector<std::set<int> > neighborNodes; // Storing node numbers beginning with 1
    std::vector<std::set<int> > neighborElements; // Storing element numbers beginning with 0
    for (int n = 0; n < nn; ++n) 
    {
        std::set<int> my_neighborNodes;
        std::set<int> my_neighborElements;
        for (int ele = 0; ele < nel; ++ele) 
        {
            for (int i = 0; i < nen; ++i) 
            {
                if (connect[ele][i] == n + 1)
                {
                    my_neighborNodes.insert(connect[ele].begin(), connect[ele].end());
                    my_neighborElements.insert(ele);
                    break;
                }
            }
        }
        neighborNodes.push_back(my_neighborNodes);
        neighborElements.push_back(my_neighborElements);
        no_nonzeros += my_neighborNodes.size()*nsd + my_neighborElements.size();
    }
    no_nonzeros = no_nonzeros*nsd + nel*(nen*nsd + 1); // For each node, every dof is identical. Plus pressure terms.
    
    int count = 0; // the number of nonzeros
    std::vector<int> col_ind; // size = no_nonzeros
    std::vector<int> row_ptr; // size = nsd*nn + nel + 1
    // First nn*nsd rows
    for (int n = 0; n < nn; ++n) 
    {
        // Deal with nsd rows corresponding to one node
        for (int i = 0; i < nsd; ++i) 
        {
            row_ptr.push_back(count);
            // Stress terms
            for (auto j : neighborNodes[n]) 
            {
                for (int k = 0; k < nsd; ++k) 
                {
                    col_ind.push_back(nsd*(j-1) + k);
                    count++;
                }
            }
            // Pressure terms
            for (auto j : neighborElements[n])
            {
                col_ind.push_back(nn*nsd + j);
                count++;
            }
        }
    }
    // Last nel rows
    for (int ele = 0; ele < nel; ++ele) 
    {
        row_ptr.push_back(count);
        // Put the nodes into a set in order to sort them
        std::set<int> myNodes(connect[ele].begin(), connect[ele].end());
        for (auto i : myNodes) 
        {
            for (int j = 0; j < nsd; ++j) 
            {
                col_ind.push_back(nsd*(i-1) + j);
                count++;
            }
        }
        col_ind.push_back(nn*nsd + ele);
        count++;
    }
    row_ptr.push_back(count);

    // Increment col_ind and row_ptr for FORTRAN code to use
    for (auto& i : col_ind) i++;
    for (auto& i : row_ptr) i++;

    assert(count == no_nonzeros);
    assert(row_ptr.size() == nsd*nn + nel + 1);
    assert(col_ind.size() == no_nonzeros);
    assert(row_ptr[row_ptr.size()-1] == no_nonzeros + 1);
    
    std::ofstream ofs("CRS.bin");
    ofs.write(reinterpret_cast<char*>(&no_nonzeros), sizeof(int));
    ofs.write(reinterpret_cast<char*>(&col_ind[0]), no_nonzeros * sizeof(int));
    int size_row_ptr = nsd*nn + nel + 1;
    ofs.write(reinterpret_cast<char*>(&size_row_ptr), sizeof(int));
    ofs.write(reinterpret_cast<char*>(&row_ptr[0]), size_row_ptr * sizeof(int));
    ofs.close();
}

void writeCRS_v2()
{
    std::set<Entry> data;
    int count = 0;
    int row, col;
    for (int ele = 0; ele < nel; ++ele) 
    {
        for (int a = 0; a < nen; ++a) 
        {
            for (int i = 0; i < nsd; ++i) 
            {
                row = nsd*(connect[ele][a] - 1) + i;
                for (int b = 0; b < nen; ++b) 
                {
                    for (int j = 0; j < nsd; ++j) 
                    {
                        col = nsd*(connect[ele][b] - 1) + j;
                        data.insert(Entry(row, col));
                    }
                }
                col = nsd*nn + ele;
                data.insert(Entry(row, col));
                data.insert(Entry(col, row));
            }
        }
        row = nsd*nn + ele;
        col = nsd*nn + ele;
        data.insert(Entry(row, col));
    }
    
    std::vector<int> col_ind;
    std::vector<int> row_ptr;
    std::vector<int> row_ind;
    
    count = 1;
    col_ind.push_back(data.begin()->col+1);
    row_ind.push_back(data.begin()->row+1);
    row_ptr.push_back(1);
    auto itr_pre = data.begin();
    auto itr = data.begin();
    itr++;
    for (; itr != data.end(); ++itr) 
    {
        col_ind.push_back(itr->col + 1); // Fortran convention
        row_ind.push_back(itr->row + 1); // Fortran convention
        if (itr->row != itr_pre->row) 
        {
            row_ptr.push_back(count + 1); // Fortran convention
        }
        itr_pre++;
        count++;
    }
    row_ptr.push_back(data.size()+1);
    assert(col_ind.size() == data.size());
    assert(row_ind.size() == data.size());
    int no_nonzeros = data.size();
    
    int size_row_ptr = row_ptr.size();
    std::ofstream ofs("CRS.bin");
    ofs.write(reinterpret_cast<char*>(&no_nonzeros), sizeof(int));
    ofs.write(reinterpret_cast<char*>(&col_ind[0]), no_nonzeros * sizeof(int));
    ofs.write(reinterpret_cast<char*>(&size_row_ptr), sizeof(int));
    ofs.write(reinterpret_cast<char*>(&row_ptr[0]), size_row_ptr * sizeof(int));
    ofs.close();   
}