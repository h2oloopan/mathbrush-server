#include "derivation.h"


namespace scg
{


bool
derivation::less(const derivation *rhs) const
{
	if (P < rhs->P) {
		return true;
	}
	else if (P == rhs->P) {
		if (children.size() < rhs->children.size()) {
			return true;
		}
		else if (children.size() == rhs->children.size()) {
			std::vector<derivation *>::const_iterator i, j;
			i = children.begin();
			j = rhs->children.begin();
			for (; i < children.end(); ++i, ++j) {
				if (!(*i)->less(*j)) {
					return false;
				}
			}
			return true;
		}
	}
	return false;
}

bool
derivation::equal(const derivation *rhs) const
{
	if (P != rhs->P) {
		return false;
	}
	else if (P == rhs->P) {
		if (children.size() != rhs->children.size()) {
			return false;
		}
		std::vector<derivation *>::const_iterator i, j;
		i = children.begin();
		j = rhs->children.begin();
		for (; i < children.end(); ++i, ++j) {
			if (!(*i)->less(*j)) {
				return false;
			}
		}
		return true;
	}
	return false;
}


class derivation_subtree_accessor : public subtree_accessor {
public:
	explicit derivation_subtree_accessor(derivation *d_) : d(d_) { }
	expression_node operator[](unsigned i) const { return d->children[i]->expr; }
	bool verify(unsigned, const expression_node &) const { return true; }

private:
	derivation *d;
};


derivation *
derivation::copy() const
{
	derivation *d = new derivation(P);
	d->expr = expression_node(expr);
	for (std::vector<derivation *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		d->add_child((*i)->copy());
	}
	return d;
}

void
derivation::build()
{
	for (std::vector<derivation *>::iterator i = children.begin(); i != children.end(); ++i) {
		if (!(*i)->expr) {
			(*i)->build();
		}
	}

	if (!expr) {
		if (P->tbuild) {
			derivation_subtree_accessor acc(this);
			expr = P->tbuild->build(acc);
			if (P->sbuild) {
				expr->set_long_string(P->sbuild->build(acc));
			}
		}
	}
}


void
derivation::write(std::ostream &os, unsigned indent) const
{
}


std::ostream &
operator<<(std::ostream &os, const derivation &d)
{
	d.write(os);
	return os;
}


/*
static const CNFProduction *
find_production(const CNFGrammar &cnf, const CNFNonTerminal *nt, CNFGrammar::gid_t left, derivation *right)
{
	/*std::cout << "looking for production in " << nt->name << " matching " << left << ", ";
	if (right) {
		if (right->rootgid() == CNFGrammar::INVALID_GID) {
			std::cout << "(unknown)\n";
		}
		std::cout << right->rootgid() << std::endl;
	}
	else {
		std::cout << CNFGrammar::INVALID_GID << std::endl;
	}
	* /

	for (std::list<CNFProduction *>::const_iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
		const CNFProduction *P = *i;
		if (P->first == left) {
			if (right) {
				if (P->second == right->rootgid()) {
					return P;
				}
				else if ((right->rootgid() == CNFGrammar::INVALID_GID) && CNFGrammar::is_nonterminal(P->second)) {
					// FIXME
					nt_derivation *ntd = dynamic_cast<nt_derivation *>(right);
					const CNFNonTerminal *right_nt = cnf.nonterminal(P->second);
					const CNFProduction *childP = find_production(cnf, right_nt, ntd->left->rootgid(), ntd->right);
					if (childP) {
						ntd->P = childP;
						ntd->nt = right_nt;
						//std::cout << " found " << *P << std::endl;
						return P;
					}
				}
			}
			else if (P->second == CNFGrammar::INVALID_GID) {
				//std::cout << " found " << *P << std::endl;
				return P;
			}
		}
	}
	return 0;
}

derivation *
read_derivation(std::istream &is, const CNFGrammar &g)
{
	std::string name;
	is >> name;
	if (!is) {
		return 0;
	}

	//std::cout << "read name " << name << std::endl;

	if (name[0] == '[') {
		std::string ntname = name.substr(1, name.length() - 2);
		const CNFNonTerminal *nt = 0;
		if (ntname != "_anon") {
			nt = g.nonterminal(ntname);
			if (!nt) {
				return 0;
			}
		}

		unsigned nchildren;
		is >> nchildren;
		if (!is) {
			return 0;
		}

		const CNFProduction *P = 0;
		derivation *left = 0, *right = 0;
		if (nchildren > 0) {
			left = read_derivation(is, g);
			if (!left) {
				return 0;
			}
			if (nchildren == 2) {
				right  = read_derivation(is, g);
				if (!right) {
					delete left;
					return 0;
				}
			}
		}

		//std::cout << "seeking production with gids " << left->rootgid() << " and " << (right ? right->rootgid() : CNFGrammar::INVALID_GID) << std::endl;
		if (nt) {
			P = find_production(g, nt, left->rootgid(), right);
			if (!P) {
				delete left;
				delete right;
				return 0;
			}
		}

		derivation *d = new nt_derivation(P, nt, left, right, derivation_create_flags::defer_build);
		//std::cout << "created new non-terminal derivation with gid " << d->rootgid() << ":\n" << *d << std::endl;
		return d;
	}
	else {
		const CNFTerminal *term = g.terminal(name);
		if (!term) {
			return 0;
		}
		derivation *d = new term_derivation(term, derivation_create_flags::defer_build);
		//std::cout << "created new terminal derivation with gid " << d->rootgid() << ":\n" << *d << std::endl;
		return d;
	}
}*/


}
