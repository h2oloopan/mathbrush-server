#define NO_RECO_TYPES
#include "grammar.h"
#include "error.h"
#include <iostream>
#include <fstream>



int
main(int argc, char *argv[])
{
	std::map<std::string, scg::SemanticId> labels;
	scg::Grammar grammar;
	std::cin >> grammar;

	scg::CNFGrammar *cnf = grammar.make_cnf(scg::grammar_conversion_flags::allow_barren_symbols, &labels);


	std::ofstream header_out("grammar-values.h");

	header_out << "// This file was automatically generated.\n";
	header_out << "// It should not be changed directly; instead, change\n";
	header_out << "// the grammar file and run the \"header\" program.\n\n\n";

	header_out << "#ifndef SCG_GRAMMAR_VALUES_H_\n";
	header_out << "#define SCG_GRAMMAR_VALUES_H_\n\n\n";

	header_out << "#include \"dlldecl.h\"\n\n";
	header_out << "#include <string>\n";
	header_out << "#include <map>\n\n\n";


	header_out << "namespace scg {\n\n";

	header_out << "// Semantic type identifiers\n\n";
	header_out << "typedef int SemanticId;\n\n";
	header_out << "DLLDECL const SemanticId InvalidSemanticId = -1;\n\n";
	header_out << "DLLDECL const SemanticId TERMINAL_EXPR = -2;\n";
	header_out << "DLLDECL const SemanticId DEFAULT_EXPR = 0;\n";
	for (std::map<scg::SemanticId, std::string>::const_iterator i = cnf->semantic_ids.begin(); i != cnf->semantic_ids.end(); ++i) {
		header_out << "DLLDECL const SemanticId " << i->second << "_EXPR = " << i->first << ";\n";
	}

	header_out << "\n\n";

	header_out << "// Children labels\n\n" << std::endl;
	for (std::map<std::string, scg::SemanticId>::const_iterator i = labels.begin(); i != labels.end(); ++i) {
		header_out << "DLLDECL extern const unsigned " << i->first << " = " << i->second << ";\n";
	}

	header_out << "\n\n\n// Do not call the functions below this line; they are for internal parser use\n\n";
	header_out << "void initialize_grammar_typemap();\n";
	header_out << "void destroy_grammar_typemap();\n";
	header_out << "std::map<std::string, SemanticId> &grammar_typemap();\n";

	header_out << "\n}\n\n#endif\n";


	header_out.close();


	std::ofstream cc_out("grammar-values.cc");

	cc_out << "// This file was automatically generated.\n";
	cc_out << "// It should not be changed directly; instead, change\n";
	cc_out << "// the grammar file and run the \"header\" program.\n\n\n";

	cc_out << "// map type strings to semantic type values\n\n";
	cc_out << "#include <string>\n";
	cc_out << "#include <map>\n\n";
	cc_out << "#include \"grammar-values.h\"\n\n";

	cc_out << "namespace scg {\n\n";


	cc_out << "std::map<std::string, SemanticId> *grammar_typemap_ = 0;\n\n";

	cc_out << "std::map<std::string, SemanticId> &\ngrammar_typemap()\n{\n";
	cc_out << "\treturn *grammar_typemap_;\n";
	cc_out << "}\n\n";

	cc_out << "void\ninitialize_grammar_typemap()\n{\n";
	cc_out << "\tif (grammar_typemap_) return;\n";
	cc_out << "\tgrammar_typemap_ = new std::map<std::string, SemanticId>;\n";
	for (std::map<scg::SemanticId, std::string>::const_iterator i = cnf->semantic_ids.begin(); i != cnf->semantic_ids.end(); ++i) {
		cc_out << "\t(*grammar_typemap_)[std::string(\"" << i->second << "\")] = " << i->second << "_EXPR;\n";
	}
	cc_out << "}\n\n";
	cc_out << "void\ndestroy_grammar_typemap()\n{\n";
	cc_out << "\tdelete grammar_typemap_;\n";
	cc_out << "\tgrammar_typemap_ = 0;\n";
	cc_out << "}\n\n";
	
	cc_out << "}\n";

	cc_out.close();
	
	delete cnf;

	return 0;
}

