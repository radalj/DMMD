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

// Test function
void test_ReadFasta() {
    Rcout << "Testing ReadFasta_cpp vs ReadFasta (R)... ";

    // Debug: inspect where ReadFasta is defined and its type
    run_R_code("if (exists('ReadFasta', envir=.GlobalEnv)) cat('typeof(ReadFasta) in .GlobalEnv:', typeof(ReadFasta), '\n')");
    run_R_code("if ('DMMD' %in% loadedNamespaces()) {\n  ns <- asNamespace('DMMD');\n  cat('exists(ReadFasta) in DMMD namespace:', exists('ReadFasta', envir=ns, inherits=FALSE), '\n');\n  if (exists('ReadFasta', envir=ns, inherits=FALSE)) cat('typeof(ReadFasta) in DMMD namespace:', typeof(get('ReadFasta', envir=ns)), '\n');\n} else {\n  cat('DMMD namespace not loaded.\\n');\n}");
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

// Main
int main() {
    int argc = 2;
    char *argv[] = {const_cast<char*>("R"), const_cast<char*>("--silent")};
    Rf_initEmbeddedR(argc, argv);

    // Diagnostics: show library paths and Rcpp info inside embedded R
    run_R_code("if (requireNamespace('Rcpp', quietly=TRUE)) {\n  library(Rcpp);\n  cat('Embedded Rcpp version:', as.character(packageVersion('Rcpp')), '\n');\n  print(getLoadedDLLs()[['Rcpp']]);\n} else {\n  cat('Rcpp not found in embedded R .libPaths()\\n');\n}");

    // Ensure the R implementation of ReadFasta is available
    run_R_code("if (requireNamespace('DMMD', quietly=TRUE)) {\n  library(DMMD);\n} else {\n  source('../R/MethylDNAFunc.R');\n}");

    test_ReadFasta();

    Rf_endEmbeddedR(0);
    return 0;
}
