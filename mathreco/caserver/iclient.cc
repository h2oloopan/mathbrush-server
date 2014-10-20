#include "caclient.h"
#include "MathRecoTypes.h"
#include "grammar-values.h"
#include "annotate.h"
#include "client-priv.h"

#include <stack>

class test_tree : public scg::ExpressionTree {
public:
	test_tree(scg::SemanticId type) : scg::ExpressionTree(), type_(type) { }
	void release() const { 
		for (std::vector<test_tree *>::const_iterator i = children.begin(); i != children.end(); ++i) {
			(*i)->release();
		}
		delete this;
	}
	scg::SemanticId type() const { return type_; }
	const char *str() const { return str_.c_str(); }
	const char *long_str() const { return str_.c_str(); }
	const char *latex_str() const { return str_.c_str(); }
	const char *getstr(int) const { return str_.c_str(); }
	const wchar_t *wstr() const { return wstr_.c_str(); }
	double score() const { return 0.0; }
	const scg::ExpressionBox *box() const { return 0; }
	size_t nchildren() const { return children.size(); }
	const ExpressionTree *child(size_t i) const { return children[i]; }
	int lock() const { return 0; }
	int unlock() const { }
	bool is_locked() const { return false; }
	bool HasLongForm() const { return false; }
	scg::ExpressionIterator *CreateLongFormIterator() const { return 0; }
	
	void add_child(test_tree *ch) { children.push_back(ch); }

public:
	scg::SemanticId type_;
	std::string str_;
	std::wstring wstr_;
	std::vector<test_tree *> children;
};

static test_tree *
pt2exptree(const scg::portable_tree &pt) {
	test_tree *tree = new test_tree(pt.sid);
	if (pt.sid == scg::TERMINAL_EXPR) {
		std::string val;
		std::wstring wval;
		for (std::list<std::string>::const_iterator i = pt.content.begin(); i != pt.content.end(); ++i) {
			scg::symbol *S = scg::symdb_findsymbol_name(*i);
			if (!S) {
				tree->release();
				return 0;
			}
			val += *i;
			wval += S->unicode;
		}
		tree->str_ = val;
		tree->wstr_ = wval;
	}
	else {
		for (std::list<scg::portable_tree>::const_iterator i = pt.children.begin(); i != pt.children.end(); ++i) {
			test_tree *child = pt2exptree(*i);
			if (!child) {
				tree->release();
				return 0;
			}
			tree->add_child(child);
		}	
	}
	return tree;
}

static void
docmd(int cmd) {
	scg::ExpressionTree *expr = exprcmd(cmd);
	if (!expr) {
		std::cerr << "caclient: error: " << casclient_replystring(casclient_geterror()) << std::endl;
	}
	else {
		std::cout << expr->long_str() << std::endl;
		expr->release();
	}
}

int
main(int argc, char *argv[]) {
	char buf[1024];
	size_t n;

	int e;
	e = casclient_connect();
	scg::VerboseOutputToFile("verbout");
	scg::SetVerbosity(1);
	if (e < 0) {
		return e;
	}

	std::stack<test_tree *> trees;

	while (!std::cin.eof()) {
		while (std::isspace(std::cin.peek())) std::cin.get();
		if (std::cin.peek() == '(') {
			scg::portable_tree pt;
			e = scg::read_portable_tree(std::cin, pt);
			test_tree *tree = pt2exptree(pt);
			if (!tree) {
				std::cout << "error creating expression tree from tree input\n";
			}
			else {
				trees.push(tree);
				std::cout << "local tree is " << trees.top()->str() << std::endl;
			}
		}
		else {
			std::string cmd;
			std::cin >> cmd;
			if (cmd == "quit") {
				break;
			}
			else if (cmd == "push") {
				if (trees.empty()) {
					std::cerr << "caclient: tree stack is empty\n";
				}
				else {
					test_tree *tree = trees.top();
					trees.pop();
					if ((e = casclient_pushexpr(tree)) != CASREP_OK) {
						std::cerr << "caclient: error: " << casclient_replystring(e) << std::endl;
					}
					else {
						std::cout << "ok.\n";
					}
					//scg::ExpressionTree *rep = casclient_factor(tree);
					//if (rep) std::cout << "response: " << rep->str() << std::endl;
					//else std::cout << "no response.\n";
				}
			}
			else if (cmd == "test") {
				if (casclient_isconnected()) std::cout << "connected.\n";
				else std::cout << "not connected!\n";
			}
			else if (cmd == "eval") docmd(CASCMD_EVAL);
			else if (cmd == "evalnum") docmd(CASCMD_EVALNUM);
			else if (cmd == "simplify") docmd(CASCMD_SIMPLIFY);
			else if (cmd == "solve") docmd(CASCMD_SOLVE);
			else if (cmd == "expand") docmd(CASCMD_EXPAND);
			else if (cmd == "factor") docmd(CASCMD_FACTOR);
			else if (cmd == "det") docmd(CASCMD_DETERMINANT);
			else if (cmd == "inv") docmd(CASCMD_INVERSE);
			else if (cmd == "rank") docmd(CASCMD_RANK);
			else if (cmd == "nullspace") docmd(CASCMD_NULLSPACE);
		}
	}

	if ((e = casclient_disconnect()) != CASREP_OK) {
		std::cerr << "caclient: error: " << casclient_replystring(e) << std::endl;
	}

	return 0;
}
