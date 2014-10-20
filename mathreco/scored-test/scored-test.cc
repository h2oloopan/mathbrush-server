#include "testing.h"
#include "error.h"
#include "annotate.h"
#include "MathRecognizer.h"

#include <iostream>
#include <fstream>
#include <cstring>

int
main(int argc, char *argv[]) {
	scg::RecognizerHandle rh = scg::InitializeRecognizer();

	scg::test_context tctx;
	int dbnamei = 1;
	for (;;) {
		if (!strcmp(argv[dbnamei], "verb")) {
			scg::SetVerbosity(1);
			scg::VerboseOutputToStream(std::cout);
			dbnamei++;
		}
		else if (!strcmp(argv[dbnamei], "auto")) {
			tctx.flags = scg::TEST_SKIP_SYMBOL_RECOGNIZER;
			dbnamei++;
		}
		else if (!strcmp(argv[dbnamei], "dpi")) {
			std::stringstream ss;
			ss << argv[dbnamei + 1];
			size_t dpi;
			ss >> dpi;
			scg::SetTabletResolution(dpi);
			dbnamei += 2;
		}
		else {
			break;
		}
	}

	tctx.dbname = argv[dbnamei];
	//std::ifstream fin(argv[dbnamei + 1]);
	for (int i = dbnamei + 1; i < argc; ++i) {
	//while (fin) {
		/*std::string fname;
		fin >> fname;
		if (!fin) {
			break;
		}*/
		scg::AnnotatedStrokeGroup grp;
		std::ifstream in(argv[i]);//fname.c_str());//argv[i]);
		if (in.is_open()) {
			if (!FAILURE(scg::import_annotated_ink(in, grp)) && grp.symbols_begin() != grp.symbols_end()) {
				tctx.filenames.push_back(argv[i]);//fname);//argv[i]);
			}
		}
	}

	int e;
	e = scg::init_test_context(&tctx);
	if (FAILURE(e)) {
		std::cerr << "stest: init with status " << scg::error_message(e) << std::endl;
		return -1;
	}
	e = scg::run_tests(&tctx);
	if (FAILURE(e)) {
		std::cerr << "stest: tests failed with status " << scg::error_message(e) << std::endl;
	}

	scg::close_test_context(&tctx);

	scg::ShutdownRecognizer(rh);

	return 0;
}
