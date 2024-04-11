
void FindConsensus (string& consensus, vector<string>& seqs);
void FindVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindAmbiguousVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindVariants2 (vector<int>& vpos, vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FixVariants (vector<string>& seqs, vector<sparseseq>& variants);
void RevertToN (vector<int>& xpos, char a1, char a2, char ax, vector<string>& seqs, vector<sparseseq>& variants);
void RevertToN3 (vector<int>& xpos, char a1, char a2, char a3, char ax, vector<string>& seqs, vector<sparseseq>& variants);

void FindNPos (vector<int>& npos, vector<sparseseq>& variants);
void MakeCounts (vector<double>& count, const vector<int>& npos, const vector<sparseseq>& variants);

void FindPairwiseDistances (run_params p, vector<int>& npos, vector<double>& count, vector<char>& count_consensus, vector< vector<double> >& seqdists, vector<sparseseq>& variants, vector<string>& seqs);
void GenerateSubsets(const vector<string>& seqs, const vector< vector<double> >& seqdists, vector< vector<int> >& subsets);
void InitialisePoints(run_params& p, const vector< vector<int> >& subsets, vector< vector<double> >& points);
void GeneratePset (const vector<string>& seqs, const vector< vector<int> >& subsets, vector<int>& pset);

void InitilisePDistMatrix (const vector< vector<double> >& points, vector< vector<double> >& pdists);
void CalculatePointDistanceMatrix (const vector< vector<double> >& points, vector< vector<double> >& pdists);
double PointsDist(const vector<double>& a, const vector<double>& b);
double GetDistanceComparison(const vector< vector<double> >& seqdists, const vector<int>& pset, const vector< vector<double> >& pdists);
