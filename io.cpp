#include "mapping.h"
#include "io.h"
#include <iostream>
#include <string>
#include <cstring>

void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
    p.dim=2;
    p.seqs_file="Sets.in";
    p.verb=0;
    p.fix=0;
    while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--input")==0) {
                x++;
                p.seqs_file=argv[x];
        } else if (p_switch.compare("--dim")==0) {
                x++;
                p.dim=atoi(argv[x]);
        } else if (p_switch.compare("--verb")==0) {
                x++;
                p.verb=atoi(argv[x]);
        } else if (p_switch.compare("--fix")==0) {
                x++;
                p.fix=atoi(argv[x]);
        } else {
			cout << "Incorrect usage\n ";
            cout << p_switch << "\n";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ReadFastaAli (run_params p, vector<string>& seqs) {
    ifstream ali_file;
    ali_file.open(p.seqs_file.c_str());
    vector<string> names;
    string seq;
    string str;
    for (int i=0;i<1000000;i++) {
        if (!(ali_file >> str)) break;
        if (str.at(0)=='>') {
            names.push_back(str);
            if (p.verb==1) {
                cout << "Read seqname " << str << "\n";
            }
            if (seq.size()>0) {
                seqs.push_back(seq);
                seq.clear();
            }
        } else {
            seq=seq+str;
        }
    }
    if (seq.size()>0) {
        seqs.push_back(seq);
    }
}

void PrintSeqdists (const vector< vector<double> >& seqdists) {
    for (int i=0;i<seqdists.size();i++) {
        for (int j=0;j<seqdists[i].size();j++) {
            cout << seqdists[i][j] << " ";
        }
        cout << "\n";
    }
}

void PrintSubsets (const vector< vector<int> >& subsets) {
    cout << "Subsets\n";
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            cout << subsets[i][j] << " ";
        }
        cout << "\n";
    }
}

void PrintPointDistanceMatrix (const vector< vector<double> >& pdists) {
    for (int i=0;i<pdists.size();i++) {
        for (int j=0;j<pdists[i].size();j++) {
            cout << pdists[i][j] << " ";
        }
        cout << "\n";
    }
}

void OutputPositions (vector<int>& positions) {
    ofstream p_file;
    p_file.open("Variant_positions.dat");
    for (int i=0;i<positions.size();i++) {
        p_file << positions[i]+1 << "\n";
    }
}

void OutputSeqDistances (vector< vector<double> > seqdists) {
    ofstream sd_file;
    sd_file.open("Sequence_distances.dat");
    for (int i=0;i<seqdists.size();i++) {
        for (int j=0;j<seqdists[i].size();j++) {
            sd_file << seqdists[i][j] << " ";
        }
        sd_file << "\n";
    }
}


void OutputPoints (const vector<string>& seqs, const vector<int>& pset, const vector< vector<double> >& points) {
    ofstream point_file;
    point_file.open("Output_points.dat");
    for (int i=0;i<seqs.size();i++) {
        for (int j=0;j<points[pset[i]].size();j++) {
            point_file << points[pset[i]][j] <<  " ";
        }
        point_file << "\n";
    }
}

void OutputDistance (double dist_best) {
    ofstream dist_file;
    dist_file.open("Distances.dat");
    dist_file << dist_best << "\n";
}

void OutputSubsets (const vector< vector<int> >& subsets) {
    ofstream sub_file;
    sub_file.open("Subsets.out");
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            sub_file << subsets[i][j] << " ";
        }
        sub_file << "\n";
    }
}
