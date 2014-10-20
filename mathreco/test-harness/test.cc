//#include "expat.h"
#include "MathRecognizer.h"
#include "annotate.h"
#include "group.h"
#include "stroke.h"
#include "verb.h"
#include "expr-node.h"
#include "grammar-values.h"
#include "lsp.h"
#include "normalize.h"
#include "segment.h"
#include "symbols.h"
#include "ink-io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <sstream>
#include <stack>
#include <queue>
#include <vector>
#include <cctype>
#include <cstring>


struct ExploreData {
	scg::ExpressionIterator *it;
	const scg::ExpressionTree *expr;
	bool rm;

	ExploreData(scg::ExpressionIterator *it_, const scg::ExpressionTree *expr_, bool rm_ = true)
		: it(it_), expr(expr_), rm(rm_)
		{ }
};


static void
do_pop(std::stack<ExploreData> &xstack)
{
	ExploreData &xd = xstack.top();
	if (xd.it) {
		if (xd.expr && xd.rm) {
			xd.expr->release();
		}
		if (xd.it) {
			xd.it->release();
		}
	}
	xstack.pop();
}

static void
clearstack(std::stack<ExploreData> &xstack) {
	while (!xstack.empty()) do_pop(xstack);
}


static int
build_known_symbols(std::queue<int> &args, const scg::AnnotatedStrokeGroup &strokes, std::vector<scg::MathSymbol> &symbols)
{
	size_t nsymbols = args.size();
	symbols.reserve(nsymbols);
	while (!args.empty()) {
		size_t index = args.front();
		args.pop();
		if (index >= strokes.symbol_annotations.size()) {
			std::cerr << "ts: known symbol index " << index << " is too large; there are only " << strokes.symbol_annotations.size() << " annotated symbols\n";
			return -1;
		}

		scg::AnnotatedStrokeGroup::const_symbol_iterator it = strokes.symbols_begin();
		while (index--) {
			++it;
		}

		const scg::symbol_annotation &A = *it;

		scg::MathSymbol ms;
		scg::symbol *S = scg::symdb_findsymbol_name(A.name);
		if (!S) {
			std::cerr << "ts: annotated symbol " << A.name << " unknown to symbols DB\n";
			return -1;
		}
		ms.symbol = S->unicode;
		ms.bounds = strokes.strokes[A.strokes.front()].bounds();
		std::list<size_t>::const_iterator strokei = A.strokes.begin();
		++strokei;
		for ( ; strokei != A.strokes.end(); ++strokei) {
			ms.bounds = merge(ms.bounds, strokes.strokes[*strokei].bounds());
		}

		symbols.push_back(ms);
	}

	return 0;
}

/*
static std::string inksrc;
static std::vector<scg::RawStroke *> svec;
static const int ST_DEFAULT = 0;
static const int ST_INK = 1;
static const int ST_INKSRC = 2;
static int state = ST_DEFAULT;
static std::string buf;
static bool decimal;


static void
addstroke() {
	size_t n = 0;
	decimal = (buf.find_first_of('.') != std::string::npos);
	std::string::size_type i = buf.find_first_of(',');
	while (i != std::string::npos) {
		n++;
		i = buf.find_first_of(',', i+1);
	}
	n++;

	long *xbuf = new long[n];
	long *ybuf = new long[n];
	i = buf.find_first_of(',');
	std::string::size_type at = 0;
	size_t pti = 0;
	for (size_t j = 0; j < n; ++j) {
		std::stringstream ss;
		if (i != std::string::npos) {
			ss.str(buf.substr(at, i - at));
		}
		else {
			ss.str(buf.substr(at));
		}

		long x, y;
		if (!decimal) {
			ss >> x >> y;
		}
		else {
			double dx, dy;
			ss >> dx >> dy;
			x = (long)(dx*2540);
			y = (long)(dy*2540);
		}
		at = i + 1;
		i = buf.find_first_of(',', at);
		if (pti == 0 || xbuf[pti-1] != x || ybuf[pti-1] != y) {
			xbuf[pti] = x;
			ybuf[pti] = y;
			++pti;
		}
	}

	scg::RawStroke *stk = new scg::RawStroke(xbuf, ybuf, 0, pti);
	svec.push_back(stk);
}

static const char *
getattr(const XML_Char **attrs, const std::string &ptn) {
	while (*attrs && strcasecmp(*attrs, ptn.c_str()) != 0) {
		attrs += 2;
	}
	if (!*attrs) {
		return 0;
	}
	return attrs[1];
}

static void
xmlstart(void *data, const XML_Char *name, const XML_Char **attrs) {
	if (!strcasecmp(name, "trace")) {
		state = ST_INK;
		const char *id = getattr(attrs, "id");
	}
	else if (!strcasecmp(name, "annotation")) {
		const char *type = getattr(attrs, "type");
		if (!strcasecmp(type, "copyright")) {
			state = ST_INKSRC;
		}
	}
}

static void
xmlend(void *data, const XML_Char *name) {
	if (!strcasecmp(name, "trace")) {
		state = ST_DEFAULT;
		addstroke();
		buf.clear();
	}
	else if (state == ST_INKSRC && !strcasecmp(name, "annotation")) {
		state = ST_DEFAULT;
	}
}

static void
xmlcdata(void *data, const XML_Char *cdata, int len) {
	if (state == ST_INK) {
		buf.append(cdata, cdata + len);
	}
	else if (state == ST_INKSRC) {
		inksrc.append(cdata, cdata + len);
	}
}

static int
convert_mathml_to_ink(std::istream &in, scg::AnnotatedStrokeGroup &strokes) {
	char buf[1024];
	XML_Parser xml = XML_ParserCreate(0);
	XML_SetElementHandler(xml, &xmlstart, &xmlend);
	XML_SetCharacterDataHandler(xml, &xmlcdata);
	std::streamsize curr = in.gcount();
	while (in) {
		in.read(buf, sizeof(buf));
		XML_Parse(xml, buf, in.gcount() - curr, XML_FALSE);
	}
	XML_Parse(xml, 0, 0, XML_TRUE);
	XML_ParserFree(xml);

	if (inksrc == "LUNAM/IRCCyN") {
		if (!decimal) {
			scg::SetTabletResolution(400);
		}
		scg::SetWriterPkgid(6);
	}
	else if (inksrc == "CVPRU/ISI") {
		scg::SetTabletResolution(40);
		scg::SetWriterPkgid(5);
	}
	else if (inksrc == "KAIST") {
		scg::SetTabletResolution(2540);
		scg::SetWriterPkgid(7);
	}
	scg::RawStroke *stks = new scg::RawStroke[svec.size()];
	for (size_t i = 0; i < svec.size(); ++i) {
		stks[i] = svec[i]->copy();
	}
	strokes.set_strokes(stks, svec.size());
	return 0;
}*/

static const unsigned NVARS = 8;
struct vartree {
	std::string type;
	std::list<vartree> children;

	void load(const scg::ExpressionTree *expr) {
		const scg::basic_tree *tree = (const scg::basic_tree *)expr;
		if (tree->type() == scg::TERMINAL_EXPR || !tree->nt()) {
			type = tree->str();
		}
		else {
			assert(tree->type() != scg::InvalidSemanticId);
			type = scg::sid_to_string(tree->type());
		}
		for (size_t i = 0; i < expr->nchildren(); ++i) {
			children.push_back(vartree());
			children.back().load(expr->child(i));
		}
	}

	bool operator==(const vartree &rhs) const {
		if (children.size() != rhs.children.size()) return false;
		if (type != rhs.type) return false;
		std::list<vartree>::const_iterator i,j;
		for (i = children.begin(), j = rhs.children.begin(); i != children.end(); ++i, ++j) {
			if (*i != *j) return false;
		}
		return true;
	}
	
	bool operator!=(const vartree &rhs) const {
		return !(*this == rhs);
	}
};

std::ostream &
operator<<(std::ostream &os, const vartree &tree) {
	if (tree.children.empty()) {
		return os << tree.type;
	}
	os << '(';
	os << tree.type;
	for (std::list<vartree>::const_iterator i = tree.children.begin(); i != tree.children.end(); ++i) {
		os << ' ' << *i;
	}
	return os << ')';
}

static bool EXIT_ON_ASSERT = true;

static char *
loadtree(char *p, vartree &tree) {
	static const char *delim = ")( ";
	tree.type.clear();
	if (*p != '(') {
		while (*p && !strchr(delim, *p)) {
			tree.type += *p;
			++p;
		}
		while (*p == ' ') ++p;
		return p;
	}
	else {
		p++;
		while (*p && !strchr(delim, *p)) {
			tree.type += *p;
			++p;
		}
		while (*p == ' ') ++p;
		while (*p != ')') {
			if (!*p) abort();
			tree.children.push_back(vartree());
			p = loadtree(p, tree.children.back());
		}
		++p;
		while (*p == ' ') ++p;
		return p;
	}
}

struct varentry {
	vartree var;
	int aschild;
	varentry() : aschild(-1) { }
	explicit varentry(const vartree &var_) : var(var_), aschild(-1) { }
	varentry(const vartree &var_, int aschild_) : var(var_), aschild(aschild_) { }
};

int
main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "usage: ts [verbose] [nofail] script-file" << std::endl;
		return -1;
	}

	for (int i = 1; i < argc - 1; ++i) {
		if (!strcmp(argv[i], "verbose")) {
			scg::VerboseOutputToStream(std::cout);
			scg::SetVerbosity(1);
		}
		if (!strcmp(argv[i], "ipad")) {
			scg::SetTabletResolution(132);
		}
		else if (!strcmp(argv[i], "nofail")) {
			EXIT_ON_ASSERT = false;
		}
	}

	int e = scg::InitializeRecognizer();
	if (FAILURE(e)) {
		std::cerr << "failed to initialize recognizer: " << scg::error_message(e) << std::endl;
		return -1;
	}
	//scg::SetUserProfilePath("/cygdrive/c/src/mathBrushResearchII/Code/C++/mathreco/TRUNK/src/Documents");

	std::istream *is;
	if (strcmp(argv[argc-1], "-") == 0) {
		is = &std::cin;
	}
	else {
		is = new std::ifstream(argv[argc-1]);
	}

	char buf[256];
	static const char *delim = " \n\r\t";

	std::stack<ExploreData> xstack;
	std::stack<varentry> varstack;

	vartree vars[NVARS];
	std::queue<int> args;
	std::string strarg;
   scg::MathRecognizer *rec = 0;
	scg::AnnotatedStrokeGroup strokes;


new_ink_file:

	clearstack(xstack);

	if (rec) {
		rec->release();
		rec = 0;
	}

	is->getline(buf, sizeof(buf));
	if (is->fail()) {
		std::cerr << "ts: script is empty" << std::endl;
		return -2;
	}

	std::ifstream ink(buf);
	if (!ink.is_open()) {
		std::cerr << "ts: ink file " << buf << " does not exist" << std::endl;
		return -2;
	}
	/*if (strlen(buf) > 5 && !strncmp(buf+strlen(buf)-5, "inkml", 5)) {
		if (convert_mathml_to_ink(ink, strokes) < 0) {
			std::cerr << "ts: could not load mathml file " << buf << std::endl;
			return -2;
		}
	}
	else*/ {
		if (scg::import_annotated_ink(ink, strokes) != 0) {
			std::cerr << "ts: ink file " << buf << " is invalid" << std::endl;
			return -2;
		}
	}

	for (;;) {
		/*
		if (is == &std::cin) {
			std::cout << "> " << std::flush;
		}*/
		is->getline(buf, sizeof(buf));
		if (is->fail()) {
			break;
		}
		
		if (buf[0] == '#') {
			std::cout << "comment: " << buf + 1 << std::endl;
			continue;
		}

		char *p = buf;
		if (buf[0] == '(') {
			varstack.push(varentry());
			p = loadtree(buf, varstack.top().var);
		}
		else if (buf[0] == '\"') {
			char *q = strchr(buf+1, '\"');
			if (!q) {
				std::cerr << "ts: missing closing \" in string\n";
			}
			else {
				strarg.assign(buf+1, q-buf-1);
			}
			p = q+1;
		}
		p = std::strtok(p, delim);
		while (p) {
			VERBOSE(*scg::verb_out << "COMMAND IS " << p << std::endl);
			if (std::isdigit(*p)) {
				std::stringstream ss;
				ss << p;
				unsigned index;
				ss >> index;
				args.push(index);
			}
			else {
				if (!strcmp(p, "new")) {
					while (!args.empty()) {
						args.pop();
					}
					strokes.clear();
					scg::SetVerbosity(1);
					scg::NoVerboseOutput();
					if (rec) {
						rec->release();
						rec = 0;
					}
					goto new_ink_file;
				}
				else if (!strcmp(p, "profile")) {
					if (strarg.empty()) {
						std::cerr << "ts: profile command needs string argument" << std::endl;
					}
					else {
						scg::SetUserProfile(strarg.c_str());
						strarg.clear();
					}
				}
				else if (!strcmp(p, "tree")) {
					if (xstack.empty() || !xstack.top().expr) {
						std::cerr << "ts: no expression on top of stack" << std::endl;
					}
					else {
						varstack.push(varentry());
						varstack.top().var.load(xstack.top().expr);
					}
				}
				else if (!strcmp(p, "tpop")) {
					if (!varstack.empty()) {
						varstack.pop();
					}
				}
				else if (!strcmp(p, "tclear")) {
					while (!varstack.empty()) varstack.pop();
				}
				else if (!strcmp(p, "tchild")) {
					if (varstack.empty()) {
						std::cerr << "ts: no vartree on top of stack" << std::endl;
					}
					else if (args.empty() || args.front() >= varstack.top().var.children.size()) {
						std::cerr << "ts: missing or invalid child index" << std::endl;
						if (!args.empty()) args.pop();
					}
					else {
						std::list<vartree>::const_iterator i = varstack.top().var.children.begin();
						int n = args.front();
						args.pop();
						varentry entry;
						entry.aschild = n;
						while (n--) ++i;
						entry.var = *i;
						varstack.push(entry);
					}
				}
				else if (!strcmp(p, "tshow")) {
					if (varstack.empty()) {
						std::cerr << "ts: no vartree on top of stack" << std::endl;
					}
					else {
						std::cout << varstack.top().var << std::endl;
					}
				}
				else if (!strcmp(p, "trepl")) {
					if (varstack.size() < 2) {
						std::cerr << "ts: no vartree on top of stack" << std::endl;
					}
					else {
						vartree repl = varstack.top().var;
						varstack.pop();
						varentry entry = varstack.top();
						entry.var = repl;
						while (!varstack.empty() && entry.aschild != -1) {
							varstack.pop();
							std::list<vartree>::iterator i = varstack.top().var.children.begin();
							while (entry.aschild--) ++i;
							*i = entry.var;
							entry = varstack.top();
						}
					}
				}
				else if (!strcmp(p, "setvar")) {
					if (args.empty() || args.front() < 0 || args.front() >= NVARS) {
						std::cerr << "ts: must specify variable id (0 through " << NVARS-1 << ")\n";
					}
					else if (varstack.empty()) {
						std::cerr << "ts: no vartree on top of stack" << std::endl;
					}
					else {
						vars[args.front()] = varstack.top().var;
						varstack.pop();
						args.pop();
					}
				}
				else if (!strcmp(p, "getvar")) {
					if (args.empty()) {
						std::cerr << "ts: must specify variable id (0 through " << NVARS-1 << ")\n";
					}
					else {
						int idx = args.front();
						args.pop();
						if (idx < 0 || idx >= NVARS) {
							std::cerr << "ts: must specify variable id (0 through " << NVARS-1 << ")\n";
						}
						else {
							varstack.push(varentry(vars[idx]));
						}
					}
				}
				else if (!strcmp(p, "asserteq")) {
					if (args.size() < 2) {
						std::cerr << "ts: must specify two vartree variables (0 through " << NVARS-1 << ") to compare\n";
						if (!args.empty()) args.pop();
					}
					else {
						int i = args.front();
						args.pop();
						int j = args.front();
						args.pop();
						if (vars[i] == vars[j]) {
							std::cout << "OK.\n";
						}
						else {
							std::cout << "FAIL.\n";
							if (EXIT_ON_ASSERT) {
								exit(0);
							}
						}
					}
				}
				else if (!strcmp(p, "translate")) {
					if (args.size() < 2) {
						std::cerr << "ts: translate require x and y parameters" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					if (!rec) {
						std::cerr << "ts: translate requires an active recognizer" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					int x = args.front();
					args.pop();
					int y = args.front();
					args.pop();
					rec->Translate(x, y);
				}
				else if (!strcmp(p, "setverbosity")) {
					if (args.empty()) {
						scg::SetVerbosity(1);
					}
					else {
						scg::SetVerbosity(args.front());
						while (!args.empty()) {
							args.pop();
						}
					}
					scg::VerboseOutputToStream(std::cout);
				}
				else if (!strcmp(p, "add")) {
					if (!rec) {
						std::cerr << "ts: ignoring add command prior to MathRecognizer creation" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					else if (!args.empty()) {
						unsigned nstrokes = args.size();
						scg::RawStroke *add_strokes = new scg::RawStroke[nstrokes];
						scg::RawStroke *strokep = add_strokes;
						while (!args.empty()) {
							unsigned index = args.front();
							args.pop();
							if (index < strokes.nstrokes) {
								scg::RawStroke S = strokes.strokes[index].copy();
								*(strokep++) = S;
							}
							else {
								std::cerr << "ts: stroke index " << index << " exceeds stroke count " << strokes.nstrokes << std::endl;
							}
						}

						rec->AddStrokes(add_strokes, nstrokes);
						clearstack(xstack);
					}
				}
				else if (!strcmp(p, "dpi")) {
					if (!args.empty()) {
						scg::SetTabletResolution(args.front());
						args.pop();
					}
				}
				else if (!strcmp(p, "addknown")) {
					if (!rec) {
						std::cerr << "ts: ignoring add command prior to MathRecognizer creation" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					else if (!args.empty()) {
						std::vector<scg::MathSymbol> symbols;
						int e = build_known_symbols(args, strokes, symbols);
						if (e < 0) {
							return e;
						}
						rec->AddKnownSymbols(&symbols.front(), symbols.size());
						clearstack(xstack);
					}
				}
				else if (!strcmp(p, "removeknown")) {
					if (!rec) {
						std::cerr << "ts: ignoring add command prior to MathRecognizer creation" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					else if (!args.empty()) {
						std::vector<scg::MathSymbol> symbols;
						int e = build_known_symbols(args, strokes, symbols);
						if (e < 0) {
							return e;
						}
						rec->RemoveKnownSymbols(&symbols.front(), symbols.size());
						clearstack(xstack);
					}
				}
				else if (!strcmp(p, "remove")) {
					if (!rec) {
						std::cerr << "ts: ignoring remove command prior to MathRecognizer creation" << std::endl;
						while (!args.empty()) {
							args.pop();
						}
					}
					else if (args.empty()) {
						if (xstack.empty()) {
							std::cerr << "ts: must either provide removal indices or have an expression on the stack\n";
						}
						else {
							std::vector<size_t> strokes;
							const scg::basic_tree *tree = (const scg::basic_tree *)xstack.top().expr;
							const scg::ordered_segments *segs = tree->span();
							if (!segs) {
								std::cerr << "ts: expression atop stack did not arise from physical strokes\n";
							}
							else {
								for (scg::ordered_segments::ordered_iterator kk = segs->begin(0); kk != segs->end(0); ++kk) {
									const scg::segment *seg = *kk;
									strokes.push_back(seg->pos);
								}
								clearstack(xstack);
								rec->RemoveStrokesByIndex(&strokes[0], strokes.size());
							}
						}
					}
					else {
						size_t nstrokes = args.size();
						size_t *indices = new size_t[nstrokes];
						size_t *indexp = indices;
						while (!args.empty()) {
							*(indexp++) = args.front();
							args.pop();
						}

						rec->RemoveStrokesByIndex(indices, nstrokes);
						clearstack(xstack);
					}
				}
				else if (!strcmp(p, "create")) {
					if (rec) {
						rec->release();
					}

					if (args.empty()) {
						//try {
							rec = scg::CreateMathRecognizer(strokes);
						/*}
						catch (int e) {
							std::cerr << "ts: error creating MathRecognizer: " << scg::error_message(e) << std::endl;
							return -2;
						}*/
					}
					else {
						scg::RawStroke *selected_strokes = new scg::RawStroke[args.size()];
						scg::RawStroke *strokep = selected_strokes;
						unsigned nstrokes = args.size();
						while (!args.empty()) {
							unsigned index = args.front();
							args.pop();

							if (index < strokes.nstrokes) {
								scg::RawStroke S = strokes.strokes[index].copy();
								*(strokep++) = S;
							}
							else {
								std::cerr << "is: ignoring invalid stroke index " << index << std::endl;
								--nstrokes;
							}
						}

						scg::RawStrokeGroup selected_group;
						selected_group.set_strokes(selected_strokes, nstrokes);
						try {
							rec = scg::CreateMathRecognizer(selected_group);
						}
						catch (int e) {
							std::cerr << "ts: error creating MathRecognizer: " << scg::error_message(e) << std::endl;
							return -2;
						}
					}
				}
				else if (!strcmp(p, "report")) {
					unsigned n = 0;
					while (!args.empty()) {
						n = args.front();
						args.pop();
					}

					if (n) {
						++n;
					}

					scg::ExpressionIterator *it;
					if (!xstack.empty()) {
						if (xstack.top().it) {
							it = xstack.top().it;
						}
						else {
							it = rec->CreateSubtreeIterator(xstack.top().expr);
							xstack.top().it = it;
						}
					}
					else {
						it = rec->CreateDefaultIterator();
					}

					if (!it) {
						std::cerr << "ts: error obtaining iterator" << std::endl;
					}
					else {
						for (const scg::ExpressionTree *expr = it->next(); expr; expr = it->next()) {
							std::cout << "report: " << expr->str() << std::endl;
							std::wcout << "wstr: " << expr->wstr() << std::endl;
							std::cout << "long string: " << expr->long_str() << std::endl;
							std::cout << "score: " << expr->score() << std::endl;
							expr->release();

							if (n) {
								if (--n == 1) {
									break;
								}
							}
						}
					}
				}
				else if (!strcmp(p, "train")) {
					if (xstack.empty()) {
						std::cerr << "ts: you must train on an iterator on the exploration stack" << std::endl;
					}
					else {
						std::cerr << "ts: training currently disabled" << std::endl;
						//xstack.top().it->TrainCorrectRelations();
					}
				}
				else if (!strcmp(p, "explore")) {
					if (!rec) {
						std::cerr << "ts: you must recognize the ink before exploring" << std::endl;
					}
					else {
						scg::ExpressionIterator *it = rec->CreateDefaultIterator();
						if (!it) {
							std::cerr << "ts: error creating iterator" << std::endl;
						}
						else {
							xstack.push(ExploreData(it, 0));
						}
					}
				}
				else if (!strcmp(p, "top")) {
					if (!rec) {
						std::cerr << "ts: you must recognize the ink before exploring" << std::endl;
					}
					else {
						const scg::ExpressionTree *toptree = rec->GetTopExpression();
						if (!toptree) {
							std::cerr << "ts: no top tree exists\n";
						}
						else {
							xstack.push(ExploreData(0, toptree));
						}
					}					
				}
				else if (!strcmp(p, "pop")) {
					if (xstack.empty()) {
						std::cerr << "ts: iterator stack empty prior to pop" << std::endl;
					}
					else {
						do_pop(xstack);
					}
				}
				else if (!strcmp(p, "clear")) {
					while (!xstack.empty()) {
						do_pop(xstack);
					}
				}
				else if (!strcmp(p, "rclear")) {
					while (!xstack.empty()) {
						do_pop(xstack);
					}
					if (rec) {
						rec->Clear();
					}
				}
				else if (!strcmp(p, "child")) {
					if (args.empty()) {
						std::cerr << "ts: child index must be specified (from 0)" << std::endl;
					}
					else if (xstack.empty()) {
						std::cerr << "ts: you must explore before extracting children" << std::endl;
					}
					else {
						ExploreData &xd = xstack.top();
						unsigned index = args.front();
						args.pop();
						if (!xd.expr) {
							std::cerr << "ts: no expression exists yet" << std::endl;
						}
						if (index >= xd.expr->nchildren()) {
							std::cerr << "ts: child index " << index << " exceeds child count " << xd.expr->nchildren() << std::endl;
						}
						else {
							const scg::ExpressionTree *child = xd.expr->child(index);
							/*scg::ExpressionIterator *it = rec->CreateSubtreeIterator(child);
							if (!it) {
								std::cerr << "ts: error creating iterator" << std::endl;
							}
							else {
								xstack.push(ExploreData(it, child));
							}*/
							xstack.push(ExploreData(0, child, false));
						}
					}
				}
				else if (!strcmp(p, "next")) {
					if (xstack.empty()) {
						std::cerr << "ts: you must create an iterator before using \'next\'" << std::endl;
					}
					else {
						ExploreData &xd = xstack.top();
						if (!xd.it) {
							//std::cerr << "ts: no iterator on top of stack" << std::endl;
							xd.it = rec->CreateSubtreeIterator(xd.expr);
							if (!xd.it) {
								std::cerr << "ts: error creating iterator" << std::endl;
							}
							else {
								xd.expr = xd.it->next();
							}
						}
						else {
							xd.expr = xd.it->next();
						}
					}
				}
				else if (!strcmp(p, "nth")) {
					if (xstack.empty()) {
						std::cerr << "ts: you must create an iterator before using \'nth\'" << std::endl;
					}
					else {
						ExploreData &xd = xstack.top();
						unsigned n = args.front();
						args.pop();
						if (!xd.it) {
							std::cerr << "ts: no iterator on top of stack" << std::endl;
						}
						else {
							xd.expr = xd.it->nth(n);
						}
					}
				}
				else if (!strcmp(p, "count")) {
					if (xstack.empty()) {
						std::cerr << "ts: you must create an iterator before using \'nth\'" << std::endl;
					}
					else {
						ExploreData &xd = xstack.top();
						if (!xd.it) {
							std::cerr << "ts: no iterator on top of stack" << std::endl;
						}
						else {
							std::cout << xd.it->count() << std::endl;
						}
					}
				}
				else if (!strcmp(p, "save")) {
					if (!rec) {
						std::cerr << "ts: need an active context to save\n";
					}
					else {
						std::ofstream ofs("ctx.out", std::ofstream::binary);
						rec->Save(ofs);
					}
				}
				else if (!strcmp(p, "restore")) {
					if (rec) {
						std::cerr << "ts: need an null context to restore\n";
					}
					else {
						std::ifstream ifs("ctx.out", std::ifstream::binary);
						rec = scg::CreateMathRecognizer(ifs);
						if (!rec) {
							const scg::error &e = scg::get_error();
							if (e.code == E_OUTOFDATE) {
								std::cerr << "ts: file is out of date\n";
							}
							else {
								std::cerr << "ts: restore error\n";
							}
						}
					}
				}
				else if (!strcmp(p, "show")) {
					if (xstack.empty() || !xstack.top().expr) {
						std::cerr << "ts: no expression on top of stack" << std::endl;
					}
					else {
						const scg::ExpressionTree *expr = xstack.top().expr;
						const scg::interpretation *intrp = ((const scg::parsed_tree *)expr)->intrp();
						std::cout << "subset: " << intrp->span->bits << std::endl;
						std::cout << "report: " << expr->str() << std::endl;
						std::cout << "intrp: " << intrp->str() << std::endl;
						//std::wstring ws = expr->wstr();
						//std::wcout << "wstr: " << expr->wstr() << std::endl;
						//std::cout << "long string: " << expr->long_str() << std::endl;
						//std::cout << "latex string: " << expr->latex_str() << std::endl;
						std::cout << "score: " << expr->score() << std::endl;
					}
				}
				else if (!strcmp(p, "nchildren")) {
					if (xstack.empty() || !xstack.top().expr) {
						std::cerr << "ts: no expression on top of stack" << std::endl;
					}
					else {
						const scg::ExpressionTree *expr = xstack.top().expr;
						std::cout << expr->nchildren() << std::endl;
					}
				}
				else if (!strcmp(p, "reset")) {
					if (!rec) {
						std::cerr << "no recognizer to reset\n";
					}
					else {
						std::cerr << "resetting rec " << rec << std::endl;
						std::stringstream ss;
						int e = rec->Save(ss);
						if (FAILURE(e)) {
							std::cerr << "save failed: " << scg::error_message(e) << std::endl;
						}
						else {
							rec->release();
							rec = scg::CreateMathRecognizer(ss);
							if (!rec) {
								std::cerr << "failed to re-create reco from save\n";
							}
							else {
								std::cerr << "reset succeeded: new rec " << rec << std::endl;
							}
						}
					}
				}
				else if (!strcmp(p, "lock")) {
					if (xstack.empty() || !xstack.top().expr || !xstack.top().it) {
						std::cerr << "ts: no expression and/or iterator on top of stack" << std::endl;
					}
					else {
						xstack.top().expr->lock();
					}
				}
				else if (!strcmp(p, "unlock")) {
					if (xstack.empty() || !xstack.top().expr || !xstack.top().it) {
						std::cerr << "ts: no expression and/or iterator on top of stack" << std::endl;
					}
					else {
						xstack.top().expr->unlock();
					}
				}
				else if (!strcmp(p, "islocked")) {
					if (xstack.empty() || !xstack.top().expr || !xstack.top().it) {
						std::cerr << "ts: no expression and/or iterator on top of stack" << std::endl;
					}
					else {
						const scg::ExpressionTree *tree = xstack.top().expr;
						bool locked = tree->is_locked();
						std::cout << (locked ? "yes" : "no") << std::endl;
					}
				}
				else if (!strcmp(p, "haslong")) {
					if (xstack.empty() || !xstack.top().expr || !xstack.top().it) {
						std::cerr << "ts: no expression and/or iterator on top of stack" << std::endl;
					}
					else {
						std::cout << (xstack.top().expr->HasLongForm() ? "yes" : "no") << std::endl;
					}				
				}
				else if (!strcmp(p, "getlong")) {
					if (xstack.empty() || !xstack.top().expr || !xstack.top().it) {
						std::cerr << "ts: no expression and/or iterator on top of stack" << std::endl;
					}
					else {
						scg::ExpressionIterator *it = xstack.top().expr->CreateLongFormIterator();
						if (!it) {
							std::cerr << "ts: error creating long-form iterator\n";
						}
						else {
							xstack.push(ExploreData(it, 0));
						}
					}				
				}
				else if (!strcmp(p, "quit")) {
					goto done_all;
				}
				else {
					std::cerr << "ts: invalid command " << p << std::endl;
				}
			}

			p = std::strtok(0, delim);
		}
	}

done_all:

	clearstack(xstack);

	if (rec) {
		rec->release();
		rec = 0;
	}

	if (is != &std::cin) {
		delete is;
	}

	scg::ShutdownRecognizer();

	return 0;
}
