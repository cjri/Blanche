
void GetParameters (run_params& p, int argc, const char **argv);
void ReadFastaAli (run_params p, vector<string>& seqs);
void PrintSeqdists (const vector< vector<double> >& seqdists);
void PrintSubsets (const vector< vector<int> >& subsets);

void PrintPointDistanceMatrix (const vector< vector<double> >& pdists);
void OutputPositions (vector<int>& positions);
void OutputSeqDistances (vector< vector<double> > seqdists);
void OutputPoints (const vector<string>& seqs, const vector<int>& pset, const vector< vector<double> >& points);
void OutputDistance (double dist_best);
void OutputSubsets (const vector< vector<int> >& subsets);
