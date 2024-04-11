#include "mapping.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>


void FindConsensus (string& consensus, vector<string>& seqs) {
    consensus=seqs[0];
    if (seqs.size()<2) {
        cout << "Error: Need at least two sequences\n";
    }
    //cout << "Find consensus of all input sequences\n";
    int nA=0;
    int nC=0;
    int nG=0;
    int nT=0;
    for (int pos=0;pos<seqs[0].size();pos++) {
        nA=0;
        nC=0;
        nG=0;
        nT=0;
        for (int seq=0;seq<seqs.size();seq++) {
            if (seqs[seq][pos]=='A') {
                nA++;
            }
            if (seqs[seq][pos]=='C') {
                nC++;
            }
            if (seqs[seq][pos]=='G') {
                nG++;
            }
            if (seqs[seq][pos]=='T') {
                nT++;
            }
        }
        int max=nA;
        consensus[pos]='A';
        if (nC>max) {
            max=nC;
            consensus[pos]='C';
        }
        if (nG>max) {
            max=nG;
            consensus[pos]='G';
        }
        if (nT>max) {
            consensus[pos]='T';
            max=nT;
        }
        if (max==0) {
            consensus[pos]='-';
        }
    }
}

void FindVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //cout << "Find variants\n";
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<seqs[i].size();pos++) {
            if (seqs[i].compare(pos,1,consensus,pos,1)!=0) {
                if (seqs[i].compare(pos,1,"A")==0||seqs[i].compare(pos,1,"C")==0||seqs[i].compare(pos,1,"G")==0||seqs[i].compare(pos,1,"T")==0) {
                    //cout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << seqs[i][pos] << "\n";
                    s.locus.push_back(pos);
                    s.allele.push_back(seqs[i][pos]);
                }
            }
        }
        variants.push_back(s);
    }
}

void FindAmbiguousVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //Variants which are not A, C, G, or T
    vector<int> vpos;
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            vpos.push_back(variants[i].locus[j]);
        }
    }
    sort(vpos.begin(),vpos.end());
    vpos.erase(unique(vpos.begin(),vpos.end()),vpos.end());
    FindVariants2 (vpos,variants,consensus,seqs);
}

void FindVariants2 (vector<int>& vpos, vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    //cout << "Find variants\n";
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<vpos.size();pos++) {
            if (seqs[i].compare(vpos[pos],1,consensus,vpos[pos],1)!=0) {
                int found=0;
                for (int j=0;j<variants[i].locus.size();j++) {
                    if (variants[i].locus[j]==vpos[pos]) {
                        found=1;
                    }
                }
                if (found==0) {
                    variants[i].locus.push_back(vpos[pos]);
                    variants[i].allele.push_back(seqs[i][vpos[pos]]);
                }
            }
        }
    }
}

void FixVariants (vector<string>& seqs, vector<sparseseq>& variants) {
    //Fix for R
    vector<int> rpos;
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='R') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'G','A','R',seqs,variants);
    //Fix for K
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='K') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'T','G','K',seqs,variants);

    //Fix for Y
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='Y') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'C','T','Y',seqs,variants);
    
    //Fix for W
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='W') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'A','T','W',seqs,variants);

    //Fix for S
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='S') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'G','C','S',seqs,variants);

    //Fix for M
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='M') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN (rpos,'A','C','M',seqs,variants);

    //Fix for B
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='B') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN3 (rpos,'C','G','T','B',seqs,variants);

    //Fix for D
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='D') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN3 (rpos,'A','G','T','D',seqs,variants);

    //Fix for H
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='H') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN3 (rpos,'A','C','T','H',seqs,variants);

    //Fix for V
    rpos.clear();
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].locus.size();j++) {
            if (variants[i].allele[j]=='V') {
                rpos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(rpos.begin(),rpos.end());
    rpos.erase(unique(rpos.begin(),rpos.end()),rpos.end());
    RevertToN3 (rpos,'A','C','G','V',seqs,variants);

}

void RevertToN (vector<int>& xpos, char a1, char a2, char ax, vector<string>& seqs, vector<sparseseq>& variants) {
    //cout << "Revert to N\n";
    for (int i=0;i<xpos.size();i++) {
        int n1=0;
        int n2=0;
        for (int j=0;j<seqs.size();j++) {
            //cout << "Here " << j << " " << xpos[i] << " " << seqs[j][xpos[i]] << "\n";
            if (seqs[j][xpos[i]]==a1) {
                n1++;
            } else if (seqs[j][xpos[i]]==a2) {
                n2++;
            }
        }
        //cout << "Count " << n1 << " " << n2 << "\n";
        if (n1>0&&n2>0) {
            //Go through variants and convert uncertain values to N
            for (int j=0;j<variants.size();j++) {
                for (int k=0;k<variants[j].locus.size();k++) {
                    if (variants[j].allele[k]==ax) {
                        variants[j].allele[k]='N';
                    }
                }
            }
        }
    }
}

void RevertToN3 (vector<int>& xpos, char a1, char a2, char a3, char ax, vector<string>& seqs, vector<sparseseq>& variants) {
    //cout << "Revert to N\n";
    for (int i=0;i<xpos.size();i++) {
        int n1=0;
        int n2=0;
        int n3=0;
        for (int j=0;j<seqs.size();j++) {
            //cout << "Here " << j << " " << xpos[i] << " " << seqs[j][xpos[i]] << "\n";
            if (seqs[j][xpos[i]]==a1) {
                n1++;
            } else if (seqs[j][xpos[i]]==a2) {
                n2++;
            } else if (seqs[j][xpos[i]]==a3) {
                n3++;
            }
        }
        int fix=0;
        if (n1>0&&n2>0) {
            fix=1;
        }
        if (n1>0&&n3>0) {
            fix=1;
        }
        if (n2>0&&n3>0) {
            fix=1;
        }
        if (fix==1) {
            //Go through variants and convert uncertain values to N
            for (int j=0;j<variants.size();j++) {
                for (int k=0;k<variants[j].locus.size();k++) {
                    if (variants[j].allele[k]==ax) {
                        variants[j].allele[k]='N';
                    }
                }
            }
        }
    }
}


void FindNPos (vector<int>& npos, vector<sparseseq>& variants) {
    //Find N alleles at variant positions
    for (int i=0;i<variants.size();i++) {
        for (int j=0;j<variants[i].allele.size();j++) {
            if (variants[i].allele[j]=='N') {
                npos.push_back(variants[i].locus[j]);
            }
        }
    }
    sort(npos.begin(),npos.end());
    npos.erase(unique(npos.begin(),npos.end()),npos.end());
}

void MakeCounts (vector<double>& count, const vector<int>& npos, const vector<sparseseq>& variants) {
    for (int i=0;i<npos.size();i++) {
        double c=0;
        for (int j=0;j<variants.size();j++) {
            for (int k=0;k<variants[j].locus.size();k++) {
                if (variants[j].locus[k]==npos[i]) {
                    if (variants[j].allele[k]!='N') {
                        c++;
                    }
                }
            }
        }
        c=c/(variants.size()-1.);
        count.push_back(c);
    }
}

void FindPairwiseDistances (run_params p, vector<int>& npos, vector<double>& count, vector<char>& count_consensus, vector< vector<double> >& seqdists, vector<sparseseq>& variants, vector<string>& seqs) {
    
    //To add: Account for N variables in the sequence.
    //N versus a consensus allele scores a value specified by the count vector
    //N versus a non-consensus allele scores 1-count
    //N versus N scores zero
    //Probably need the consensus allele if we are doing this here.  Or could pre-specify that in count_consensus?
    
    vector<double> zeros(seqs.size(),0.0);
    for (int i=0;i<seqs.size();i++) {
        seqdists.push_back(zeros);
    }
    for (int i=0;i<seqs.size();i++) {
        for (int j=i+1;j<seqs.size();j++) {
            double dist=0;
            //Find unique difference positions between sequences i and j
            vector<int> uniq;
            for (int k=0;k<variants[i].locus.size();k++) {
                uniq.push_back(variants[i].locus[k]);
            }
            for (int k=0;k<variants[j].locus.size();k++) {
                uniq.push_back(variants[j].locus[k]);
            }
            sort(uniq.begin(),uniq.end());
            uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
            for (int k=0;k<uniq.size();k++) {
                if (seqs[i][uniq[k]]=='A'||seqs[i][uniq[k]]=='C'||seqs[i][uniq[k]]=='G'||seqs[i][uniq[k]]=='T') {
                    if (seqs[j][uniq[k]]=='A'||seqs[j][uniq[k]]=='C'||seqs[j][uniq[k]]=='G'||seqs[j][uniq[k]]=='T') {
                        if (seqs[i][uniq[k]]!=seqs[j][uniq[k]]) {
                            dist++;
                        }
                    }
                }
            }
            if (p.fix==1) {//Include an accounting for N nucleotides
                //Look through positions containing an N
                for (int k=0;k<npos.size();k++) {
                    int found=0;
                    if (seqs[i][npos[k]]=='N'&&seqs[j][npos[k]]==count_consensus[k]) {
                     //   cout << "Here1\n";
                        dist=dist+count[k];
                     //   cout << "i " << i << " j " << j << " added " << count[k] << "\n";
                        found=1;
                    }
                    if (seqs[j][npos[k]]=='N'&&seqs[i][npos[k]]==count_consensus[k]) {
                     //   cout << "Here2\n";
                        dist=dist+count[k];
                     //   cout << "i " << i << " j " << j << " added " << count[k] << "\n";
                        found=1;
                    }
                    if (seqs[j][npos[k]]=='N'&&seqs[i][npos[k]]=='N') {
                        found=1;
                    }
                    if (seqs[i][npos[k]]=='N'||seqs[j][npos[k]]=='N') {
                        if (found==0) {
                          //  cout << "Here3\n";
                            dist=dist+(1-count[k]);
                          //  cout << "i " << i << " j " << j << " added " << 1-count[k] << "\n";
                            found=1;
                        }
                    }
                }
            }
            seqdists[i][j]=dist;
            seqdists[j][i]=dist;
        }
    }
    for (int i=0;i<seqs.size();i++) {
        if (seqs[i].size()==0) {
            for (int j=0;j<seqs.size();j++){
                seqdists[i][j]=-1;
                seqdists[j][i]=-1;
            }
            seqdists[i][i]=0;
        }
    }
}

void GenerateSubsets(const vector<string>& seqs, const vector< vector<double> >& seqdists, vector< vector<int> >& subsets) {
    //Identify subsets of identical sequences
    vector<int> found;
    for (int i=0;i<seqs.size();i++) {
        found.push_back(0);
    }
    for (int i=0;i<seqs.size();i++) {
        vector<int> s;
        if (found[i]==0) {
            s.push_back(i);
            for (int j=i+1;j<seqs.size();j++) {
                if (seqdists[i][j]==0) {
                    s.push_back(j);
                    found[j]=1;
                }
            }
            subsets.push_back(s);
        }
    }
}

void InitialisePoints(run_params& p, const vector< vector<int> >& subsets, vector< vector<double> >& points) {
    for (int i=0;i<subsets.size();i++) {
        vector<double> pp;
        for (int j=0;j<p.dim;j++) {
            pp.push_back(0);
        }
        points.push_back(pp);
    }
}

void GeneratePset (const vector<string>& seqs, const vector< vector<int> >& subsets, vector<int>& pset) {
    for (int i=0;i<seqs.size();i++) {
        pset.push_back(0);
    }
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            pset[subsets[i][j]]=i;
        }
    }

}

void InitilisePDistMatrix (const vector< vector<double> >& points, vector< vector<double> >& pdists) {
    for (int i=0;i<points.size();i++) {
        vector<double> p;
        for (int j=0;j<points.size();j++) {
            p.push_back(0);
        }
        pdists.push_back(p);
    }
}

void CalculatePointDistanceMatrix (const vector< vector<double> >& points, vector< vector<double> >& pdists) {
    for (int i=0;i<pdists.size();i++) {
        for (int j=0;j<pdists[i].size();j++) {
            pdists[i][j]=PointsDist(points[i],points[j]);
        }
    }
}

double PointsDist(const vector<double>& a, const vector<double>& b) {
    double d=0;
    for (int i=0;i<a.size();i++) {
        double dx=a[i]-b[i];
        d=d+pow(dx,2);
    }
    d=sqrt(d);
    return d;
}

double GetDistanceComparison(const vector< vector<double> >& seqdists, const vector<int>& pset, const vector< vector<double> >& pdists) {
    double dist=0;
    for (int i=0;i<seqdists.size();i++) {
        for (int j=0;j<seqdists.size();j++) {
            double d=seqdists[i][j]-pdists[pset[i]][pset[j]];
            dist=dist+pow(d,2);
        }
    }
    dist=sqrt(dist);
    return dist;
}
