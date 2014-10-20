#include "MathRecognizer.h"
#include "../caserver/expr.h"
#include "error.h"
#include <cstring>

namespace scg {

int
ExpressionTree::serialize(char *buf, size_t *n) const {
	if (buf && !n) return E_INVALID;

	exprtree *tree = mkexprtree(this);
	if (!tree) {
		return E_INVALID;
	}
	std::vector<char> vbuf;
	int e = tree->writestream(vbuf);
	delete tree;
	if (e != 0) {
		return E_INTERNAL;
	}

	
	if (buf) {
		if (vbuf.size() <= *n) {
			memcpy(buf, &vbuf[0], vbuf.size());
		}
		else {
			e = E_NOTREADY;
		}
	}
	if (n) {
		*n = vbuf.size();
	}
	return e;
}

ExpressionTree *
CreateExpressionTree(const char *buf, size_t n) {
	if (!buf) {
		ERR(E_INVALID, "");
		return 0;
	}

	exprtree *tree = mkexprtree(&buf, &n);
	if (!tree || n != 0) {
		ERR(E_INVALID, "");
		return 0;
	}
	scg::ExpressionTree *expr = tree->mkscgtree();
	if (!expr) {
		ERR(E_INTERNAL, "");
	}
	delete tree;
	return expr;
}

}
