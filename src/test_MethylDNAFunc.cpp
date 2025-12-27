#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include <Rembedded.h>
#include <Rinterface.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
}

using namespace Rcpp;
using namespace std;

// Declare the C++ function
List ReadFasta_cpp(List Config);
List CooMov_cpp(List Config, List CooMetFor);
List Rev_cpp(List Config, List SeqMetFreW);
List RevTot_cpp(List Config, List Seqs);
List DelGaps_cpp(List Config, List SeqMetFreW);
List DelGapsTot_cpp(List Config, List Seqs);

// Helper: run arbitrary R code in the embedded interpreter
void run_R_code(const char* code) {
    ParseStatus status;
    SEXP cmdSexp = PROTECT(Rf_mkString(code));
    SEXP cmdexpr = PROTECT(R_ParseVector(cmdSexp, -1, &status, R_NilValue));
    if (status != PARSE_OK) {
        Rf_unprotect(2);
        Rprintf("Failed to parse R code: %s\n", code);
        return;
    }
    int exprs = Rf_length(cmdexpr);
    for (int i = 0; i < exprs; ++i) {
        Rf_eval(VECTOR_ELT(cmdexpr, i), R_GlobalEnv);
    }
    Rf_unprotect(2);
}

// Helper: Call the R ReadFasta function
List call_ReadFasta_R(List Config) {
    // Prefer the DMMD namespace if available
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("ReadFasta")) {
            Function ReadFasta = dmmd["ReadFasta"];
            return ReadFasta(Config);
        }
    } catch (...) {
        // Namespace not loaded; fall back to global env
    }

    Environment global = Environment::global_env();
    if (!global.exists("ReadFasta")) {
        stop("ReadFasta not found in DMMD namespace or global environment");
    }
    Function ReadFasta = global["ReadFasta"];
    return ReadFasta(Config);
}

// Helper: Call the R CooMov function
List call_CooMov_R(List Config, List CooMetFor) {
    // Prefer DMMD namespace if loaded
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("CooMov")) {
            Function CooMov = dmmd["CooMov"];
            return CooMov(Config, CooMetFor);
        }
    } catch (...) {
        // Namespace not loaded; fall back to global env
    }

    Environment global = Environment::global_env();
    if (!global.exists("CooMov")) {
        stop("CooMov not found in DMMD namespace or global environment");
    }
    Function CooMov = global["CooMov"];
    return CooMov(Config, CooMetFor);
}

// Helper: Call the R Rev function
List call_Rev_R(List Config, List SeqMetFreW) {
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("Rev")) {
            Function Rev = dmmd["Rev"];
            return Rev(Config, SeqMetFreW);
        }
    } catch (...) {
    }

    Environment global = Environment::global_env();
    if (!global.exists("Rev")) {
        stop("Rev not found in DMMD namespace or global environment");
    }
    Function Rev = global["Rev"];
    return Rev(Config, SeqMetFreW);
}

// Helper: Call the R RevTot function
List call_RevTot_R(List Config, List Seqs) {
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("RevTot")) {
            Function RevTot = dmmd["RevTot"];
            return RevTot(Config, Seqs);
        }
    } catch (...) {
    }

    Environment global = Environment::global_env();
    if (!global.exists("RevTot")) {
        stop("RevTot not found in DMMD namespace or global environment");
    }
    Function RevTot = global["RevTot"];
    return RevTot(Config, Seqs);
}

// Helper: Call the R DelGaps function
List call_DelGaps_R(List Config, List SeqMetFreW) {
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("DelGaps")) {
            Function DelGaps = dmmd["DelGaps"];
            return DelGaps(Config, SeqMetFreW);
        }
    } catch (...) {}
    Environment global = Environment::global_env();
    if (!global.exists("DelGaps")) {
        stop("DelGaps not found in DMMD namespace or global environment");
    }
    Function DelGaps = global["DelGaps"];
    return DelGaps(Config, SeqMetFreW);
}

// Helper: Call the R DelGapsTot function
List call_DelGapsTot_R(List Config, List Seqs) {
    try {
        Environment dmmd = Environment::namespace_env("DMMD");
        if (dmmd.exists("DelGapsTot")) {
            Function DelGapsTot = dmmd["DelGapsTot"];
            return DelGapsTot(Config, Seqs);
        }
    } catch (...) {}
    Environment global = Environment::global_env();
    if (!global.exists("DelGapsTot")) {
        stop("DelGapsTot not found in DMMD namespace or global environment");
    }
    Function DelGapsTot = global["DelGapsTot"];
    return DelGapsTot(Config, Seqs);
}

// Helper: Compare two lists of CharacterVectors
bool compare_seq_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        CharacterVector va = a[i];
        CharacterVector vb = b[i];
        if (va.size() != vb.size()) return false;
        for (int j = 0; j < va.size(); ++j) {
            if (as<string>(va[j]) != as<string>(vb[j])) return false;
        }
    }
    return true;
}

// Helper: Compare two CooMov result lists by ColCoo column
bool compare_coo_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        DataFrame dfA = as<DataFrame>(a[i]);
        DataFrame dfB = as<DataFrame>(b[i]);
        if (!dfA.containsElementNamed("ColCoo") || !dfB.containsElementNamed("ColCoo")) return false;
        NumericVector colA = dfA["ColCoo"];
        NumericVector colB = dfB["ColCoo"];
        if (colA.size() != colB.size()) return false;
        for (int j = 0; j < colA.size(); ++j) {
            if (colA[j] != colB[j]) return false;
        }
    }
    return true;
}

// Helper: Compare two Rev result lists by comparing Seq columns in data.frames
bool compare_rev_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (Rf_isNull(a[i]) && Rf_isNull(b[i])) continue;
        if (Rf_isNull(a[i]) || Rf_isNull(b[i])) return false;
        
        DataFrame dfA = as<DataFrame>(a[i]);
        DataFrame dfB = as<DataFrame>(b[i]);
        
        if (!dfA.containsElementNamed("Seq") || !dfB.containsElementNamed("Seq")) return false;
        
        CharacterVector seqA = dfA["Seq"];
        CharacterVector seqB = dfB["Seq"];
        
        if (seqA.size() != seqB.size()) return false;
        
        for (int j = 0; j < seqA.size(); ++j) {
            if (as<string>(seqA[j]) != as<string>(seqB[j])) return false;
        }
    }
    return true;
}

// Helper: Compare two RevTot result lists by comparing character vectors
bool compare_revtot_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (Rf_isNull(a[i]) && Rf_isNull(b[i])) continue;
        if (Rf_isNull(a[i]) || Rf_isNull(b[i])) return false;
        
        CharacterVector vecA = as<CharacterVector>(a[i]);
        CharacterVector vecB = as<CharacterVector>(b[i]);
        
        if (vecA.size() != vecB.size()) return false;
        
        for (int j = 0; j < vecA.size(); ++j) {
            if (as<string>(vecA[j]) != as<string>(vecB[j])) return false;
        }
    }
    return true;
}

// Helper: Compare two DelGaps result lists by comparing Seq columns in data.frames
bool compare_delgaps_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (Rf_isNull(a[i]) && Rf_isNull(b[i])) continue;
        if (Rf_isNull(a[i]) || Rf_isNull(b[i])) return false;
        DataFrame dfA = as<DataFrame>(a[i]);
        DataFrame dfB = as<DataFrame>(b[i]);
        if (!dfA.containsElementNamed("Seq") || !dfB.containsElementNamed("Seq")) return false;
        CharacterVector seqA = dfA["Seq"];
        CharacterVector seqB = dfB["Seq"];
        if (seqA.size() != seqB.size()) return false;
        for (int j = 0; j < seqA.size(); ++j) {
            if (as<string>(seqA[j]) != as<string>(seqB[j])) return false;
        }
    }
    return true;
}

// Helper: Compare two DelGapsTot result lists by comparing character vectors
bool compare_delgapstot_lists(const List& a, const List& b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (Rf_isNull(a[i]) && Rf_isNull(b[i])) continue;
        if (Rf_isNull(a[i]) || Rf_isNull(b[i])) return false;
        CharacterVector vecA = as<CharacterVector>(a[i]);
        CharacterVector vecB = as<CharacterVector>(b[i]);
        if (vecA.size() != vecB.size()) return false;
        for (int j = 0; j < vecA.size(); ++j) {
            if (as<string>(vecA[j]) != as<string>(vecB[j])) return false;
        }
    }
    return true;
}

// Test function
void test_ReadFasta() {
    Rcout << "Testing ReadFasta_cpp vs ReadFasta (R)... \n";

    // Prepare a dummy Config (adjust paths as needed)
    List Config = List::create(
        Named("NumAutosomes") = 1,
        Named("Allosomes") = CharacterVector::create("X"),
        Named("DirFas") = "test_fasta"
    );
    // Call both functions
    List r_out = call_ReadFasta_R(Config);
    List cpp_out = ReadFasta_cpp(Config);
    // Compare
    bool ok = compare_seq_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}

// Test CooMov_cpp against R CooMov
void test_CooMov() {
    Rcout << "Testing CooMov_cpp vs CooMov (R)... \n";

    // Minimal Config
    List Config = List::create(
        Named("NumChr") = 2,
        Named("CooDis") = 1
    );

    // Build a simple CooMetFor: two chromosomes, each with a ColCoo column
    List CooMetFor(2);

    NumericVector col1 = NumericVector::create(1, 5, 10);
    DataFrame df1 = DataFrame::create(Named("ColCoo") = col1);
    CooMetFor[0] = df1;

    NumericVector col2 = NumericVector::create(2, 7, 20);
    DataFrame df2 = DataFrame::create(Named("ColCoo") = col2);
    CooMetFor[1] = df2;

    // Call both implementations
    List r_out = call_CooMov_R(Config, CooMetFor);
    List cpp_out = CooMov_cpp(Config, CooMetFor);

    bool ok = compare_coo_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}

// Test Rev_cpp against R Rev
void test_Rev() {
    Rcout << "Testing Rev_cpp vs Rev (R)... \n";

    // Minimal Config
    List Config = List::create(
        Named("w_min") = 5,
        Named("w_max") = 6
    );

    // Build SeqMetFreW structure: list with data.frames at positions 5 and 6
    List SeqMetFreW(6);
    
    // Fill positions 0-4 with NULL
    for (int i = 0; i < 5; ++i) {
        SeqMetFreW[i] = R_NilValue;
    }
    
    // Position 5 (w=5): data.frame with Seq, Methyl, Freq, Index
    CharacterVector seqs5 = CharacterVector::create("acgtacgtacgt", "ttttccccaaaa");
    NumericVector methyl5 = NumericVector::create(0.8, 0.3);
    IntegerVector freq5 = IntegerVector::create(10, 5);
    IntegerVector index5 = IntegerVector::create(1, 2);
    DataFrame df5 = DataFrame::create(
        Named("Seq") = seqs5,
        Named("Methyl") = methyl5,
        Named("Freq") = freq5,
        Named("Index") = index5
    );
    SeqMetFreW[4] = df5;
    
    // Position 6 (w=6): data.frame
    CharacterVector seqs6 = CharacterVector::create("ggggggggggtaaaaaaaaa");
    NumericVector methyl6 = NumericVector::create(0.5);
    IntegerVector freq6 = IntegerVector::create(3);
    IntegerVector index6 = IntegerVector::create(3);
    DataFrame df6 = DataFrame::create(
        Named("Seq") = seqs6,
        Named("Methyl") = methyl6,
        Named("Freq") = freq6,
        Named("Index") = index6
    );
    SeqMetFreW[5] = df6;

    // Call both implementations
    List r_out = call_Rev_R(Config, SeqMetFreW);
    List cpp_out = Rev_cpp(Config, SeqMetFreW);


    bool ok = compare_rev_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}

// Test RevTot_cpp against R RevTot
void test_RevTot() {
    Rcout << "Testing RevTot_cpp vs RevTot (R)... \n";

    // Minimal Config
    List Config = List::create(
        Named("w_min") = 3,
        Named("w_max") = 4
    );

    // Build Seqs structure: list with character vectors at positions 3 and 4
    List Seqs(4);
    
    // Fill positions 0-2 with NULL
    for (int i = 0; i < 3; ++i) {
        Seqs[i] = R_NilValue;
    }
    
    // Position 3 (w=3): character vector
    CharacterVector seqs3 = CharacterVector::create("acgtacgtacgt", "ttttcccctttt");
    Seqs[2] = seqs3;
    
    // Position 4 (w=4): character vector
    CharacterVector seqs4 = CharacterVector::create("ggggaaaatttt");
    Seqs[3] = seqs4;

    // Call both implementations
    List r_out = call_RevTot_R(Config, Seqs);
    List cpp_out = RevTot_cpp(Config, Seqs);


    bool ok = compare_revtot_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}

// Test DelGaps_cpp against R DelGaps
void test_DelGaps() {
    Rcout << "Testing DelGaps_cpp vs DelGaps (R)... \n";
    List Config = List::create(
        Named("w_min") = 2,
        Named("w_max") = 4
    );
    List SeqMetFreW(4);
    for (int i = 0; i < 1; ++i) SeqMetFreW[i] = R_NilValue;
    // w=2: all sequences have gaps
    CharacterVector seqs2 = CharacterVector::create("nnnn", "nacg", "acgn");
    NumericVector methyl2 = NumericVector::create(0.1, 0.2, 0.3);
    IntegerVector freq2 = IntegerVector::create(1, 2, 3);
    IntegerVector index2 = IntegerVector::create(1, 2, 3);
    DataFrame df2 = DataFrame::create(
        Named("Seq") = seqs2,
        Named("Methyl") = methyl2,
        Named("Freq") = freq2,
        Named("Index") = index2
    );
    SeqMetFreW[1] = df2;
    // w=3: two sequences, one with gap
    CharacterVector seqs3 = CharacterVector::create("acgtacgt", "ttttnccc");
    NumericVector methyl3 = NumericVector::create(0.8, 0.3);
    IntegerVector freq3 = IntegerVector::create(10, 5);
    IntegerVector index3 = IntegerVector::create(1, 2);
    DataFrame df3 = DataFrame::create(
        Named("Seq") = seqs3,
        Named("Methyl") = methyl3,
        Named("Freq") = freq3,
        Named("Index") = index3
    );
    SeqMetFreW[2] = df3;
    // w=4: one sequence, no gap
    CharacterVector seqs4 = CharacterVector::create("ggggaaaa", "tttt", "cccc");
    NumericVector methyl4 = NumericVector::create(0.5, 0.75, 0.6);
    IntegerVector freq4 = IntegerVector::create(3, 3, 3);
    IntegerVector index4 = IntegerVector::create(3, 4, 5);
    DataFrame df4 = DataFrame::create(
        Named("Seq") = seqs4,
        Named("Methyl") = methyl4,
        Named("Freq") = freq4,
        Named("Index") = index4
    );
    SeqMetFreW[3] = df4;
    List r_out = call_DelGaps_R(Config, SeqMetFreW);
    List cpp_out = DelGaps_cpp(Config, SeqMetFreW);
    bool ok = compare_delgaps_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}

// Test DelGapsTot_cpp against R DelGapsTot
void test_DelGapsTot() {
    Rcout << "Testing DelGapsTot_cpp vs DelGapsTot (R)...\n";
    List Config = List::create(
        Named("w_min") = 2,
        Named("w_max") = 5
    );
    List Seqs(5);
    Seqs[0] = R_NilValue;
    // w=2: all sequences have gaps
    CharacterVector seqs2 = CharacterVector::create("nnnn", "nacg", "acgn");
    Seqs[1] = seqs2;
    // w=3: mixed, some with gaps, some without
    CharacterVector seqs3 = CharacterVector::create("acgtacgt", "ttttnccc", "cccccccc", "nnnnnnnn", "acgtacgt");
    Seqs[2] = seqs3;
    // w=4: all valid
    CharacterVector seqs4 = CharacterVector::create("ggggaaaa", "tttt", "cccc", "acgtacgt");
    Seqs[3] = seqs4;
    // w=5: mix of upper/lower case gaps
    CharacterVector seqs5 = CharacterVector::create("acgtnacgt", "ACGTNACGT", "acgtacgt", "ACGTACGT", "nnnnnnnn", "NNNNNNNN", "acgtacgt", "acgtnnnn");
    Seqs[4] = seqs5;
    List r_out = call_DelGapsTot_R(Config, Seqs);
    List cpp_out = DelGapsTot_cpp(Config, Seqs);
    bool ok = compare_delgapstot_lists(r_out, cpp_out);
    if (ok) {
        Rcout << "\033[32mPASS\033[0m\n";
    } else {
        Rcout << "\033[31mFAIL\033[0m\n";
    }
}


// Main
int main() {
    int argc = 2;
    char *argv[] = {const_cast<char*>("R"), const_cast<char*>("--silent")};
    Rf_initEmbeddedR(argc, argv);

    // Diagnostics: show library paths and Rcpp info inside embedded R
    run_R_code("if (requireNamespace('Rcpp', quietly=TRUE)) {\n  library(Rcpp);\n  cat('Embedded Rcpp version:', as.character(packageVersion('Rcpp')), '\n');\n  print(getLoadedDLLs()[['Rcpp']]);\n} else {\n  cat('Rcpp not found in embedded R .libPaths()\\n');\n}");

    // Ensure the R implementation of ReadFasta is available
    run_R_code("if (requireNamespace('DMMD', quietly=TRUE)) {\n  library(DMMD);\n} else {\n  source('../R/MethylDNAFunc.R');\n}");
    
    test_CooMov();
    test_ReadFasta();
    test_Rev();
    test_RevTot();    
    test_DelGaps();
    test_DelGapsTot();

    Rf_endEmbeddedR(0);
    return 0;
}
