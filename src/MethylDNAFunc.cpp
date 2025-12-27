#include <Rcpp.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List CooMov_cpp(List Config, List CooMetFor) {
    // Displaces targets' coordinates (ColCoo) by Config$CooDis
    int numChr = as<int>(Config["NumChr"]);
    int cooDis = as<int>(Config["CooDis"]);

    int nInput = CooMetFor.size();
    if (numChr > nInput) {
        stop("CooMov_cpp: Config$NumChr is greater than length(CooMetFor)");
    }

    List CooMetForUpd(numChr);

    for (int i = 0; i < numChr; ++i) {
        // Each element is expected to be a data.frame with a column named "ColCoo"
        DataFrame df = as<DataFrame>(CooMetFor[i]);

        if (!df.containsElementNamed("ColCoo")) {
            stop("CooMov_cpp: element %d has no 'ColCoo' column", i + 1);
        }

        // Treat ColCoo as numeric; R will upcast integers as needed
        NumericVector col = df["ColCoo"];
        for (int j = 0; j < col.size(); ++j) {
            col[j] += cooDis;
        }
        df["ColCoo"] = col;

        CooMetForUpd[i] = df;
    }

    return CooMetForUpd;
}

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

// Helper: compute reverse complement of a DNA sequence
// Processes only the first 'len' characters
string reverse_complement(const string& seq, int len) {
    int actual_len = min(len, (int)seq.length());
    string result(actual_len, ' ');
    for (int i = 0; i < actual_len; ++i) {
        char c = seq[i];
        char complement;
        switch(c) {
            case 'a': complement = 't'; break;
            case 'c': complement = 'g'; break;
            case 'g': complement = 'c'; break;
            case 't': complement = 'a'; break;
            case 'A': complement = 'T'; break;
            case 'C': complement = 'G'; break;
            case 'G': complement = 'C'; break;
            case 'T': complement = 'A'; break;
            default: complement = c; break;  // preserve other characters (like 'n')
        }
        result[i] = complement;
    }
    // Reverse the string
    std::reverse(result.begin(), result.end());
    return result;
}

// [[Rcpp::export]]
List Rev_cpp(List Config, List SeqMetFreW) {
    // Reverse complement sequences in SeqMetFreW
    int w_min = as<int>(Config["w_min"]);
    int w_max = as<int>(Config["w_max"]);

    List result = clone(SeqMetFreW);

    for (int w = w_min; w <= w_max; ++w) {
        if (result.size() < w || Rf_isNull(result[w - 1])) continue;

        DataFrame df = as<DataFrame>(result[w - 1]);
        if (!df.containsElementNamed("Seq")) continue;

        CharacterVector seqs = df["Seq"];
        CharacterVector revSeqs(seqs.size());

        // Compute Len = 2*w + 2 to match R's Reverse function
        int len = 2 * w + 2;

        for (int i = 0; i < seqs.size(); ++i) {
            string seq_str = as<string>(seqs[i]);
            revSeqs[i] = reverse_complement(seq_str, len);
        }

        // Build new data.frame with reversed sequences
        DataFrame new_df = DataFrame::create(
            Named("Seq") = revSeqs,
            Named("Methyl") = df["Methyl"],
            Named("Freq") = df["Freq"],
            Named("Index") = df["Index"]
        );

        result[w - 1] = new_df;
    }

    return result;
}

// [[Rcpp::export]]
List RevTot_cpp(List Config, List Seqs) {
    // Reverse complement sequences in Seqs list
    int w_min = as<int>(Config["w_min"]);
    int w_max = as<int>(Config["w_max"]);

    List result = clone(Seqs);

    for (int w = w_min; w <= w_max; ++w) {
        if (result.size() < w || Rf_isNull(result[w - 1])) continue;

        CharacterVector seqs = as<CharacterVector>(result[w - 1]);
        CharacterVector revSeqs(seqs.size());

        // Compute Len = 2*w + 2 to match R's Reverse function
        int len = 2 * w + 2;

        for (int i = 0; i < seqs.size(); ++i) {
            string seq_str = as<string>(seqs[i]);
            revSeqs[i] = reverse_complement(seq_str, len);
        }

        result[w - 1] = revSeqs;
    }

    return result;
}

// [[Rcpp::export]]
List DelGaps_cpp(List Config, List SeqMetFreW) {
    int w_min = as<int>(Config["w_min"]);
    int w_max = as<int>(Config["w_max"]);
    List result = clone(SeqMetFreW);

    for (int w = w_min; w <= w_max; ++w) {
        int idx = w - 1;
        if (idx < 0 || idx >= result.size() || Rf_isNull(result[idx])) continue;

        DataFrame df = as<DataFrame>(result[idx]);
        if (!df.containsElementNamed("Seq")) continue;

        CharacterVector Seq = df["Seq"];
        NumericVector Methyl = df["Methyl"];
        IntegerVector Freq = df["Freq"];
        IntegerVector Index = df["Index"];
        int n = Seq.size();

        // Efficiently filter out sequences containing 'n' (case-insensitive)
        vector<int> keep;
        keep.reserve(n);
        for (int i = 0; i < n; ++i) {
            const string& s = as<string>(Seq[i]);
            if (std::none_of(s.begin(), s.end(), [](char c){ return c == 'n' || c == 'N'; })) {
                keep.push_back(i);
            }
        }
        if (!keep.empty() && keep.size() < n) {
            CharacterVector SeqNew(keep.size());
            NumericVector MethylNew(keep.size());
            IntegerVector FreqNew(keep.size());
            IntegerVector IndexNew(keep.size());
            for (size_t j = 0; j < keep.size(); ++j) {
                int i = keep[j];
                SeqNew[j] = Seq[i];
                MethylNew[j] = Methyl[i];
                FreqNew[j] = Freq[i];
                IndexNew[j] = Index[i];
            }
            DataFrame new_df = DataFrame::create(
                Named("Seq") = SeqNew,
                Named("Methyl") = MethylNew,
                Named("Freq") = FreqNew,
                Named("Index") = IndexNew
            );
            result[idx] = new_df;
        }
        else if (keep.empty()) {
            // If no sequences remain, create empty data.frame
            DataFrame new_df = DataFrame::create(
                Named("Seq") = CharacterVector(0),
                Named("Methyl") = NumericVector(0),
                Named("Freq") = IntegerVector(0),
                Named("Index") = IntegerVector(0)
            );
            result[idx] = new_df;
        }
    }
    return result;
}

// [[Rcpp::export]]
List DelGapsTot_cpp(List Config, List Seqs) {
    int w_min = as<int>(Config["w_min"]);
    int w_max = as<int>(Config["w_max"]);
    List result = clone(Seqs);
    for (int w = w_min; w <= w_max; ++w) {
        int idx = w - 1;
        if (idx < 0 || idx >= result.size() || Rf_isNull(result[idx])) continue;
        CharacterVector Seq = as<CharacterVector>(result[idx]);
        int n = Seq.size();
        vector<int> keep;
        keep.reserve(n);
        for (int i = 0; i < n; ++i) {
            const std::string& s = as<std::string>(Seq[i]);
            if (std::none_of(s.begin(), s.end(), [](char c){ return c == 'n' || c == 'N'; })) {
                keep.push_back(i);
            }
        }
        if (!keep.empty() && keep.size() < n) {
            CharacterVector SeqNew(keep.size());
            for (size_t j = 0; j < keep.size(); ++j) {
                SeqNew[j] = Seq[keep[j]];
            }
            result[idx] = SeqNew;
        }
        else if (keep.empty()) {
            result[idx] = CharacterVector(0);
        }
    }
    return result;
}