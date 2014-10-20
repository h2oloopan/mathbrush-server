#include "mathrecognizer-private.h"
#include "MathRecognizer.h"

#include <fstream>

namespace scg
{


/*
static symbols_db *db = 0;
static symbols_db *autotrain_db = 0;

int
rebuild_global_symbols_db()
{
	delete db;
	db = symbols_db::create();
	if (!db) {
		return get_error().code;
	}
	delete autotrain_db;
	if (IsAutoTrainingEnabled()) {
		std::string path = get_profile_path();
		std::ifstream dbin((path + "/autotrain.symbols").c_str());
		VERBOSE(*verb_out << "autotrain: loading autotrain symbols db\n");
		autotrain_db = symbols_db::create_from_stream(dbin);
		if (!autotrain_db) {
			return get_error().code;
		}
		if (FAILURE(db->augment(*autotrain_db))) {
			return get_error().code;
		}
	}
	return 0;
}


const symbols_db *
global_symbols_db() {
	return global_symbols_db_wr();
}

symbols_db *
global_symbols_db_wr() {
	if (!db) {
		if (rebuild_global_symbols_db() < 0) {
			return 0;
		}
	}
	return db;
}

symbols_db *
global_autotrain_db() {
	if (!db) {
		if (rebuild_global_symbols_db() < 0) {
			return 0;
		}
	}
	return autotrain_db;
}

void
kill_global_symbols_db()
{
	if (db) {
		db->close();
		delete db;
		db = 0;
	}
	if (autotrain_db) {
		std::string path = get_profile_path();
		std::ofstream dbout((path + "/autotrain.symbols").c_str());
		autotrain_db->write_to_stream(dbout);
		autotrain_db->close();
		delete autotrain_db;
		autotrain_db = 0;
	}
}
*/

}

