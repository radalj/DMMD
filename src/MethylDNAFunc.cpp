#include <Rcpp.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List ReadFasta_cpp(List Config) {
    // Reads chromosomes' fasta files and saves the sequences in an R structure
    int NumAutosomes = as<int>(Config["NumAutosomes"]);
    CharacterVector Allosomes = Config["Allosomes"];
    string DirFas = as<string>(Config["DirFas"]);

    List SeqChr;

    // Helper to read a fasta file and return all sequences as a CharacterVector
    auto read_fasta = [](const string& filename) -> CharacterVector {
        CharacterVector seqs;
        ifstream file(filename);
        if (!file.is_open()) {
            Rcpp::Rcout << "Warning: could not open file " << filename << endl;
            return seqs;
        }
        string line, seq;
        while (getline(file, line)) {
            if (!line.empty() && line[0] == '>') {
                if (!seq.empty()) {
                    seqs.push_back(seq);
                    seq.clear();
                }
            } else {
                // change to lowercase
                transform(line.begin(), line.end(), line.begin(), ::tolower);
                seq += line;
            }
        }
        if (!seq.empty()) seqs.push_back(seq);
        file.close();
        return seqs;
    };

    // Read autosomes
    for (int i1 = 0; i1 < NumAutosomes; i1++) {
        string NamFilFas = "chr" + to_string(i1 + 1) + ".fa";
        string FulNamFilFas = DirFas + "/" + NamFilFas;
        CharacterVector seqs = read_fasta(FulNamFilFas);
        SeqChr.push_back(seqs);
    }

    // Read allosomes
    for (int i2 = 0; i2 < Allosomes.size(); i2++) {
        string allosome = as<string>(Allosomes[i2]);
        string NamFilFas = "chr" + allosome + ".fa";
        string FulNamFilFas = DirFas + "/" + NamFilFas;
        CharacterVector seqs = read_fasta(FulNamFilFas);
        SeqChr.push_back(seqs);
    }

    return SeqChr;
}
