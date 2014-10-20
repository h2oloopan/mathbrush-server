#ifndef EXTERN_ITER_H_
#define EXTERN_ITER_H_

#include "MathRecoTypes.h"
#include <vector>

namespace scg {

class interpreter;
class math_recognizer_base;

class external_iterator : public ExpressionIterator {
private:
	struct replacement_entry {
		const basic_tree *replacee;
		const basic_tree *replacement;
		replacement_entry() : replacee(0), replacement(0) { }
		replacement_entry(const basic_tree *replacee_, const basic_tree *replacement_) : replacee(replacee_), replacement(replacement_) { }
	};

public:
	external_iterator(interpreter *src_, bool rm_, bool wrap_mathml_ = true);
	~external_iterator();
	void release();
	const ExpressionTree *next();
	const ExpressionTree *nth(size_t i);
	size_t count() const;

	interpreter *getsrc();

private:
	const basic_tree *tree_from_interpretation(interpretation *intrp);
	const basic_tree *postprocess_tree(const basic_tree *tree);

private:
	interpreter *src;
	size_t at;
	std::vector<replacement_entry> replacement_cache;
	bool rm;
	bool wrap_mathml;
};

}

#endif
