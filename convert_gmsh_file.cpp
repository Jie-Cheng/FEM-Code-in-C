#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <iomanip>

class PhysicalObj {
  public:
    int id;
    std::string name;
    int eletype; // based on the format of gmsh
    std::vector< std::vector<int> > connect;
    std::vector<int> nodes;
    friend bool operator< (const PhysicalObj& p1, const PhysicalObj& p2) {
      return (p1.id < p2.id);
    }
    void generate() { // generate nodes from the connect
      std::map<int, int> data;
      std::map<int, int>::iterator itr;
      for (int i = 0; i < connect.size(); ++i) {
        for (int j = 0; j < connect[i].size(); ++j) {
          data[connect[i][j]]++;
        }
      }
      for (itr = data.begin(); itr != data.end(); ++itr) {
        nodes.push_back(itr->first);
      }
    }
};

void trim(std::string& s) {
  size_t p = s.find_first_not_of(" \t\r");
  s.erase(0, p);
  p = s.find_last_not_of(" \t\r");
  if (std::string::npos != p)
    s.erase(p+1);
}

void output (const std::vector< std::vector<double> >& coords, \
    const std::vector<PhysicalObj>& v, const int nsd, const int nen) {

  std::cout << "Begin output...\n";

  // Read the boundary conditions
  // For input: 0 -- essential; 1 -- pressure; 2 -- traction
  int load_type; // For output: 1 -- pressure; 2 or 3 -- traction
  std::ifstream infile("key.txt");
  if (!infile) {
    std::cerr << "Could not open key.txt!" << std::endl;
    exit(0);
  }
  std::string buffer;
  std::vector< std::vector<double> > bc;
  while (std::getline(infile, buffer)) {
    if (buffer[0] == '#') {
      continue; // comment begins with #
    }
    std::vector<double> line;
    std::stringstream ss(buffer);
    double tmp;
    while (ss >> tmp) {
      line.push_back(tmp);
    }
    bc.push_back(line);
    if (line[0] == 1)
      load_type = 1;
    else if (line[0] == 2)
      load_type = nsd;
    line.clear();
  }

  // Write bc.txt and load.txt
  std::ofstream out_bc, out_load, out_coords, out_connect;
  out_bc.open("bc.txt");
  int no_bc = 0;
  out_load.open("load.txt");
  int no_load = 0;

  // Figure out no_bc and no_load
  for (int i = 0; i < bc.size(); ++i) {
    int type = int(bc[i][0]);
    int id = int(bc[i][1]);
    if (type == 0) {
      no_bc += v[id-1].nodes.size();
    } else if (type == 1 || type == 2) {
      no_load += v[id-1].connect.size();
    }
  }
  out_bc << std::setw(10) <<  no_bc << std::endl;
  out_load << std::setw(10) << no_load << '\t' << std::setw(10) <<
    load_type << std::endl;

  for (int i = 0; i < bc.size(); ++i) {
    // type of the boundary condition
    int type = int(bc[i][0]);
    int id = int(bc[i][1]);
    if (type == 0) {
      // essential boundary condition
      int dof = int(bc[i][2]);
      for (int j = 0; j < v[id-1].nodes.size(); ++j) {
        out_bc << std::setw(10) << v[id-1].nodes[j] << '\t' \
          << std::setw(10) << dof << '\t' << std::setw(20) \
          << std::setprecision(8) << std::fixed << bc[i][3] << std::endl;
      }
      no_bc += v[id-1].nodes.size();
    } else if (type == 1) {
      // pressure
      no_load += v[id-1].connect.size();
      for (int j = 0; j < v[id-1].connect.size(); ++j) {
        for (int k = 0; k < v[id-1].connect[j].size(); ++k) {
          out_load << std::setw(10) << v[id-1].connect[j][k] << '\t';
        }
        out_load << std::setw(20) << std::setprecision(8) \
          << std::fixed << bc[i][2] << std::endl;
      }
    } else if (type == 2) {
      // traction
      no_load += v[id-1].connect.size();
      for (int j = 0; j < v[id-1].connect.size(); ++j) {
        for (int k = 0; k < v[id-1].connect[j].size(); ++k) {
          out_load << std::setw(10) << v[id-1].connect[j][k] << '\t';
        }
        for (int k = 0; k < nsd; ++k) {
          out_load << std::setw(20) << std::setprecision(8) \
            << std::fixed << bc[i][2+k] << '\t';
        }
        out_load << std::endl;
      }
    }
  }
  out_bc.close();
  out_load.close();
  std::cout << "Done writing bc and load ...\n";

  // coords.txt
  out_coords.open("coords.txt");
  out_coords << std::setw(13) << nsd << '\t' << std::setw(13) \
    << coords.size() << std::endl; 
  for (int i = 0; i < coords.size(); ++i) {
    for (int j = 0; j < nsd; ++j) {
      out_coords << std::setw(13) << coords[i][j] << '\t';
    }
    out_coords << std::endl;
  }
  out_coords.close();
  std::cout << "Done writing coords ...\n";

  // connect.txt
  out_connect.open("connect.txt");
  out_connect << std::setw(10) << v[0].connect.size() << '\t' \
    << std::setw(10) << nen << std::endl; 
  for (int i = 0; i < v[0].connect.size(); ++i) {
    for (int j = 0; j < nen; ++j) {
      out_connect << std::setw(10) << v[0].connect[i][j] << '\t';
    }
    out_connect << std::endl;
  }
  out_connect.close();
  std::cout << "Done writing connect ...\n";
}

int main (int argc, char* argv[]) {
  int no_physical_objs;
  std::vector< PhysicalObj > physical_objs;
  int nsd, nen, nn;
  int nel = 0;
  int ne; // number of elements in general, as defined in gmsh
  
  std::string target;
  bool target_found = false;

  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == std::string("-target")) {
      target = argv[i+1];
      target_found = true;
    }
  }
  if (!target_found) {
    std::cerr << "Use -target to specify the input mesh file!" << std::endl;
    exit(0);
  } 

  std::vector< std::vector<double> > coords;

  std::ifstream infile(target);
  if (!infile) {
    std::cerr << "Could not open " << target << std::endl;
    exit(0);
  }
  std::string buffer;
  while (std::getline(infile, buffer)) {
    trim(buffer);

    // Physical Objects
    if (buffer == "$PhysicalNames") {
      std::getline(infile, buffer);
      trim(buffer);
      int no_physical_objs = atoi(buffer.c_str());
      for (int i = 0; i < no_physical_objs; ++i) {
        std::getline(infile, buffer);
        trim(buffer);
        PhysicalObj tmp;
        std::istringstream ss(buffer);
        std::string token;
        std::getline(ss, token, ' ');
        std::getline(ss, token, ' ');
        tmp.id = atoi(token.c_str());
        std::getline(ss, token, ' ');
        tmp.name = token;
        physical_objs.push_back(tmp);
      }
    }
    std::sort(physical_objs.begin(), physical_objs.end()); // sort the objects by their ids

    // Nodes
    if (buffer == "$Nodes") {
      std::getline(infile, buffer);
      trim(buffer);
      nn = atoi(buffer.c_str());
      std::cout << "Number of nodes: " << nn << std::endl;
      coords.resize(nn);
      for (int i = 0; i < nn; ++i) {
        std::getline(infile, buffer);
        trim(buffer);
        std::stringstream ss(buffer);
        std::string token;
        coords[i].resize(3);
        double value;
        for (int j = 0; j < 4; ++j) {
          ss >> value;
          if (j == 0) continue;
          coords[i][j-1] = value;
          std::cout << coords[i][j-1] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    // Elements
    if (buffer == "$Elements") {
      std::getline(infile, buffer);
      trim(buffer);
      ne = atoi(buffer.c_str());
      std::cout << "Number of general elements: " << ne << std::endl;
      for (int i = 0; i < ne; ++i) {
        std::vector<int> tmp;
        std::getline(infile, buffer);
        std::stringstream ss(buffer);
        int n;
        while (ss >> n) {
          tmp.push_back(n);
          std::cout << n << " ";
        }
        std::cout << std::endl;
        // Copy the connect info into physical objects
        std::vector<int> connect;
        for (int j = 3 + tmp[2]; j < tmp.size(); ++j) {
          connect.push_back(tmp[j]);
        }
        physical_objs[tmp[3] - 1].connect.push_back(connect);
        physical_objs[tmp[3] - 1].eletype = tmp[1];
        // Check if this gerneral element is the FEM element
        if (tmp[3] == 1) nel++; // The ID for the body should always be 1
        tmp.clear();
      }
      std::cout << std::endl;
    }
  }

  // nsd, nen
  if (physical_objs[0].eletype == 2) {
    nsd = 2;
    nen = 3;
  } else if (physical_objs[0].eletype == 3) { 
    nsd = 2;
    nen = 4;
  } else if (physical_objs[0].eletype == 4) { 
    nsd = 3;
    nen = 4;
  } else if (physical_objs[0].eletype == 5) { 
    nsd = 3;
    nen = 8;
  }

  // Generate the nodes and echo the physical objects
  for (int i = 0; i < physical_objs.size(); ++i) {
    physical_objs[i].generate(); // generate the nodes
    std::cout << "Name: " << physical_objs[i].name << std::endl;
    std::cout << "ID: " << physical_objs[i].id << std::endl;
    std::cout << "Size of connect: " << physical_objs[i].connect.size() << std::endl;
    for (int j = 0; j < physical_objs[i].connect.size(); ++j) {
      for (int k = 0; k < physical_objs[i].connect[j].size(); ++k) {
        std::cout << physical_objs[i].connect[j][k] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "Size of nodes: " << physical_objs[i].nodes.size() << std::endl;
    std::cout << std::endl;
  }

  // Output
  output (coords, physical_objs, nsd, nen);
  return 0; 
}








