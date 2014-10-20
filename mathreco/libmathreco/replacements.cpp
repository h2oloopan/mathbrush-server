#include "replacements.h"
#include "expr-node.h"

namespace scg {

basic_tree *
trailing_subscript_removal(const basic_tree *tree, basic_tree *&var) {
	if (tree->type() == SUBSCR_EXPR) {
		const basic_tree *head = tree->child(SUBSCR_BASE);
		const basic_tree *sub = tree->child(SUBSCR_SUB);
		if (head->type() == VAR_EXPR && sub->type() == VAR_EXPR && head->child(VAR_NAME)->ccstr() == "d") {
			var = sub->mkcopy();
		}
		return 0;
	}
	else if (tree->type() == MULT_EXPR) {
		const basic_tree *lhs = tree->child(EXPR_LHS);
		const basic_tree *rhs = tree->child(EXPR_RHS);
		if (lhs->type() == VAR_EXPR && rhs->type() == VAR_EXPR && lhs->child(VAR_NAME)->ccstr() == "d") {
			var = rhs->mkcopy();
		}
		else {
			basic_tree *rhs_repl = trailing_subscript_removal(tree->child(EXPR_RHS), var);
			if (rhs_repl || var) {
				if (rhs_repl) {
					basic_tree *repl = tree->mkcopy();
					repl->replace_child(EXPR_RHS, rhs_repl);
					return repl;
				}
				else {
					return tree->child(EXPR_LHS)->mkcopy();
				}
			}
		}
	}
	return 0;
}

static basic_tree *
integral_trailing_d_repl(const basic_tree *tree) {
	if (tree->type() == INTEGRAL_EXPR) {
		const basic_tree *var = tree->child(INTEGRAL_VAR);
		if (var->type() == BLANK_EXPR || var->type() == PLACEHOLDER_EXPR) {
			const basic_tree *integrand = tree->child(RANGE_OPERAND);
			basic_tree *var = 0;
			basic_tree *repl_integrand = trailing_subscript_removal(integrand, var);
			if (repl_integrand || var) {
				basic_tree *repl = tree->mkcopy();
				if (repl_integrand) repl->replace_child(RANGE_OPERAND, repl_integrand);
				else repl->replace_child(RANGE_OPERAND, mkblank());
				repl->replace_child(INTEGRAL_VAR, var);
				return repl;
			}
		}
	}
	return 0;
}

basic_tree *
do_replacements(const basic_tree *tree) {
	basic_tree *repl;

	repl = integral_trailing_d_repl(tree);
	if (repl) return repl;

	return 0;
}

basic_tree *
replace_tree(const basic_tree *tree) {
	basic_tree *repl = 0;

	for (size_t i = 0; i < tree->nchildren(); ++i) {
		const basic_tree *child = tree->child(i);
		basic_tree *child_repl = replace_tree(child);
		if (child_repl) {
			if (!repl) repl = tree->mkcopy();
			repl->replace_child(i, child_repl);
		}
	}

	basic_tree *root_repl;
	root_repl = do_replacements(repl ? repl : tree);
	if (root_repl) {
		if (repl) repl->release();
		return root_repl;
	}

	return repl;
}


}