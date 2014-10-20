#ifndef TESTING_H_
#define TESTING_H_

#include <string>
#include <vector>

#include "dlldecl.h"

#define BIT(n) (1 << n)

extern "C" {
	struct sqlite3;
};


namespace scg {

extern const std::string TEST_RESULTS_TABLE_NAME;

// Statistics collected
enum {
	// Stats on test status
	STAT_TESTED,
	STAT_FEASIBLE,
	STAT_RECOGNIZABLE,
	STAT_ATTAINABLE,
	STAT_CORRECT,
	// unhandled internal error occured
	// (typically indicates out-of-memory condition)
	STAT_ABORTED,

	STAT_CORRECT_GROUPS,
	STAT_CORRECT_SYMBOLS,

	// ... on effort units
	STAT_STRUCT_CORRECTIONS,
	STAT_TERM_CORRECTIONS,
	STAT_VIEWS,

	// ... on timings
	STAT_PARSE_TIME,
	STAT_EXTRACT_TIME,
	STAT_SEARCH_TIME,

	// ... on parse complexity
	STAT_EXPR_LEN,
	STAT_EXPR_HEIGHT,
	STAT_EXPR_RELS,
	STAT_NSEGMENTS,
	STAT_NSPANS,
	STAT_NGROUPS,
	STAT_NCELLS,
	STAT_NLINKS,

	STAT__N
};


// Bit flags for the flags field in test_context
enum {
	// In the DEFAULT_MODE case, the tests will use the parser
	// with no special processing.
	TEST_DEFAULT_MODE = 0,

	// In the SKIP_SYMBOL_RECOGNIZER case, the tests will use the
	// ground-truth symbols and bounding boxes	rather than invoking
	// the symbol recognizer.
	TEST_SKIP_SYMBOL_RECOGNIZER = BIT(1)
};


struct DLLDECL test_context {
	// Input: fill these fields in before calling run_tests()
	std::string dbname;
	std::vector<std::string> filenames;
	int flags;

	// Output: the db is created by init_test_context
	sqlite3 *db;

	test_context();
};



// Functions to retrieve textual descriptions of each
// statistic by index.
DLLDECL const char *get_stat_column_name(size_t i);
DLLDECL const char *get_stat_name(size_t i);
DLLDECL const char *get_stat_description(size_t i);

// Functions to control logged output of the testing module.
DLLDECL int enable_test_logging(const std::string &log_filename);
DLLDECL int disable_test_logging();

// Main entry points
DLLDECL int init_test_context(test_context *ctx);
DLLDECL int close_test_context(test_context *ctx);
DLLDECL int run_tests(test_context *ctx);


}



#endif
