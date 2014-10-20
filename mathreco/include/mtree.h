#ifndef MTREE_H_
#define MTREE_H_

#include "error.h"
#include "iohelp.h"
#include "verb.h"
#include <map>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <cmath>

namespace scg {

const unsigned QUERY_THRESHOLD = 5;
const static double MAX_VOLUME = 0.5;

double splitweight(double x);

template <unsigned NLABELDIMS>
struct label {
	enum {
		N = NLABELDIMS
	};
	int p[NLABELDIMS];
	bool operator<(const label &rhs) const {
		for (size_t i = 0; i < NLABELDIMS; ++i) {
			if (p[i] < rhs.p[i]) return true;
			if (rhs.p[i] < p[i]) return false;
		}
		return false;
	}
};

template <unsigned NDIMS, unsigned NLABELDIMS>
struct pt {
	enum {
		N = NDIMS
	};
	float x[NDIMS];

	label<NLABELDIMS> *labels;
	size_t nlabels;

	struct cmp {
	private:
		size_t i;
	public:
		explicit cmp(size_t i_) : i(i_) { }
		bool operator()(const pt &lhs, const pt &rhs) const { return lhs.x[i] < rhs.x[i]; }
	};
};

template <unsigned NDIMS, unsigned NLABELDIMS>
float
ptdist(const pt<NDIMS, NLABELDIMS> &x1, const pt<NDIMS, NLABELDIMS> &x2) {
	float D = 0.f;
	for (size_t i = 0; i < NDIMS; ++i) {
		float d = x1.x[i] - x2.x[i];
		D += d*d;
	}
	return std::sqrt(D);
}

template <unsigned NDIMS, unsigned NLABELDIMS>
struct margintree {
/*private:
	void update_margin(pt &p, size_t i) {
		if (margin[i][0] < p.x[i] && p.x[i] < margin[i][1]) {
			if (p.x[i] - margin[i][0] < margin[i][1] - p.x[i]) {
				margin[i][0] = p.x[i];
			}
			else {
				margin[i][1] = p.x[i];
			}
		}
	}*/
	typedef label<NLABELDIMS> label_t;
	typedef pt<NDIMS, NLABELDIMS> pt_t;


public:
	explicit margintree(const pt_t &data_) : low(0), high(0), dsplit(~0),  data(data_) {
		for (size_t i = 0; i < NDIMS; ++i) {
			dmin[i] = std::numeric_limits<float>::infinity();
			dmax[i] = -std::numeric_limits<float>::infinity();
		}
	}
	~margintree() { delete low; delete high; }

private:
	margintree() : low(0), high(0), dsplit(~0) { }

public:
	bool write(writer &wr) const {
		wr.write((long)NDIMS);
		wr.write((long)NLABELDIMS);
		if (!wr.isok()) return false;
		return write_internal(wr);
	}
private:
	bool write_internal(writer &wr) const {
		for (size_t i = 0; i < NDIMS; ++i) {
			wr.write(dmin[i]);
			wr.write(dmax[i]);
		}
		if (low) {
			wr.write(dsplit);
			wr.write(splitval);
			/*wr.write(n.size());
			for (typename std::map<label_t, unsigned>::const_iterator i = n.begin(); i != n.end(); ++i) {
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					wr.write((long)i->first.p[j]);
				}
				wr.write((long)i->second);
			}*/
			return low->write(wr) && high->write(wr);
		}
		else {
			long dum = NDIMS;
			wr.write(dum);
			for (size_t j = 0; j < NDIMS; ++j) {
				wr.write(data.x[j]);
			}
			wr.write(n.size());
			for (typename std::map<label_t, unsigned>::const_iterator i = n.begin(); i != n.end(); ++i) {
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					wr.write((long)i->first.p[j]);
				}
				wr.write((long)i->second);
			}
			return wr.isok();
		}
	}

	int read(reader &re);

public:
	margintree *findcell(const pt_t &p, size_t i, unsigned thres) {
		if (n[p.labels[i]] < thres) return 0;
		if (!low && !high) return this;
		margintree *qchild = (p.x[dsplit] < splitval) ? low : high;
		margintree *cell = qchild->findcell(p, i, thres);
		return cell ? cell : this;
	}

	margintree *nearest(const pt_t &p, size_t i) {
		/*std::cout << "label: [";
		for (size_t j = 0; j < NLABELDIMS; ++j) {
			std::cout << ' ' << p.labels[i].p[j];
		}
		std::cout << "]\n";

		std::cout << "at " << this << ", n = " << n[p.labels[i]] << " (";
		if (low) std::cout << low->n[p.labels[i]];
		else std::cout << "?";
		std::cout << "/";
		if (high) std::cout << high->n[p.labels[i]];
		else std::cout << "?";*/

		if (n[p.labels[i]] == 0) return 0;

		if (!low && !high) {
			//std::cout << " -> " << this << "\n";
			return this;
		}
		if (low && low->n[p.labels[i]] > 0 && (p.x[dsplit] < splitval || (high && high->n[p.labels[i]] == 0))) {
			//std::cout << " -> low\n";
			return low->nearest(p, i);
		}
		else if (high) {
			//std::cout << " -> high\n";
			return high->nearest(p, i);
		}
		//std::cout << " -> none\n";
		return 0;
	}
	
	float distsum(const pt_t &p, size_t i) {
		float d = 0.f;
		if (n[p.labels[i]] == 0) return 0.f;
		if (!low && !high) {
			d = ptdist(p, this->data);
		}
		if (low) d += low->distsum(p, i);
		if (high) d += high->distsum(p, i);
		return d;
	}

	bool query(const pt_t &p, double *output, size_t agg, size_t allagg, unsigned thres = QUERY_THRESHOLD, double C = 1.0) {
		VERBOSE(
			*verb_out << "mtree query pt [";
			for (size_t i = 0; i < NDIMS; ++i) {
				*verb_out << ' ' << p.x[i];
			}
			*verb_out << " ] with labels {";
			for (size_t i = 0; i < p.nlabels; ++i) {
				*verb_out << " [";
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << ' ' << p.labels[i].p[j];
				}
				*verb_out << "]";
			}
			*verb_out << " }\n";
		);
		margintree *dcell = findcell(p, allagg, thres);
		if (!dcell) return false;
		unsigned ndpts = dcell->n[p.labels[allagg]];
		double dist = dcell->distsum(p, allagg) / ndpts;
		double pnil = dist / (dist+3);
		//double invdensity = dcell->volume(p) / ndpts;
		//double pnil = invdensity / (invdensity + 1);

		if (!intquery(p, output, agg, allagg, thres, C)) {
			return false;
		}
		VERBOSE(
			*verb_out << "got output [";
			for (size_t i = 0; i < p.nlabels; ++i) {
				*verb_out << ' ' << output[i];
			}
			*verb_out << " ]\n";
		);
		VERBOSE2(*verb_out << " pnil = avg dist " << dist << " to " << ndpts << " points -> " << pnil << std::endl);
		//VERBOSE(*verb_out << " pnil = inv density " << invdensity << " with " << ndpts << " points -> " << pnil << std::endl);
		/*for (size_t i = 0; i < p.nlabels; ++i) {
			output[i] = output[i]/pnil - output[i];
			VERBOSE2(
				*verb_out << "  [";
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << ' ' << p.labels[i].p[j];
				}
				*verb_out << " ] -> " << output[i] << std::endl;
			);
		}*/
		return true;
	}

	bool intquery(const pt_t &p, double *output, size_t agg, size_t allagg, unsigned thres = QUERY_THRESHOLD, double C = 1.0) {
		const static double DECAY = 0.25;
		
		unsigned nagg = n[p.labels[agg]];

		VERBOSE2(
			*verb_out << "in cell [";
			for (size_t i = 0; i < NDIMS; ++i) {
				*verb_out << ' ' << dmin[i] << "->" << dmax[i];
			}
			*verb_out << " ] with " << nagg << " agg pts\n";
		);
		if (nagg < thres) {
			double v = volume(p);
			if (v > MAX_VOLUME) {
				size_t i;
				for (i = 1; i < NLABELDIMS; ++i) {
					if (p.labels[0].p[i] != 0) {
						break;
					}
				}
				if (i != NLABELDIMS) {
					VERBOSE(*verb_out << "volume " << v << " is too big...\n");
					return false;
				}
			}

			//C = (1.0 - DECAY) / (DECAY - C);
			//VERBOSE(*verb_out << " only " << n[p.labels[agg]] << " pts left; scaling by const " << C << std::endl);
			for (size_t i = 0; i < p.nlabels; ++i) {
				//output[i] *= C;
				VERBOSE2(
					*verb_out << "  [";
					for (size_t j = 0; j < NLABELDIMS; ++j) {
						*verb_out << ' ' << p.labels[i].p[j];
					}
					*verb_out << " ] -> " << output[i] << std::endl;
				);
			}
			return true;
		}

		VERBOSE2(
			*verb_out << "descending to split " << splitval << " on dim " << dsplit << std::endl;
			/**verb_out << "Freq. table: [ ";
			for (typename std::map<label_t, unsigned>::const_iterator i = n.begin(); i != n.end(); ++i) {
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << i->first.p[j] << '/';
				}
				*verb_out << " -> " << i->second << ' ';
			}
			*verb_out << "]\n";*/
		);
		
		for (size_t i = 0; i < p.nlabels; ++i) {
			double phere = (double)n[p.labels[i]] / nagg;
			if (output[i] == 0.0) output[i] = phere;
			else output[i] = DECAY*output[i] + (1.0 - DECAY)*phere;
			VERBOSE2(
				*verb_out << "  [";
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << ' ' << p.labels[i].p[j];
				}
				*verb_out << " ] adding p " << n[p.labels[i]] << "/" << nagg << " = " << phere << " -> " << output[i] << std::endl;
			);
		}
		margintree *qchild = (p.x[dsplit] < splitval) ? low : high;
		return qchild ? qchild->intquery(p, output, agg, allagg, thres, C * DECAY) : false;
	}

	/*bool query(const pt_t &p, double *output, size_t agg, size_t allagg, unsigned thres = QUERY_THRESHOLD) {
		margintree *aggcell = findcell(p, allagg, thres);
		margintree *cell = findcell(p, agg, thres);
		//margintree *dcell = nearest(p, 0);
		margintree *dcell = findcell(p, allagg, 2);

		if (!aggcell || !cell || !dcell) return false;

		VERBOSE(
			*verb_out << "Query pt [ ";
			for (size_t i = 0; i < NDIMS; ++i) {
				*verb_out << p.x[i] << ' ';
			}
			*verb_out << "]\n";
			*verb_out << "Query labels [ ";
			for (size_t i = 0; i < p.nlabels; ++i) {
				*verb_out << "[ ";
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << p.labels[i].p[j] << ' ';
				}
				*verb_out << "] ";
			}
			*verb_out << "]\n";
			*verb_out << "Freq. table: [ ";
			for (typename std::map<label_t, unsigned>::const_iterator i = cell->n.begin(); i != cell->n.end(); ++i) {
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << i->first.p[j] << '/';
				}
				*verb_out << " -> " << i->second << ' ';
			}
			*verb_out << "]\n";
		);
		
		unsigned ndpts = dcell->n[p.labels[allagg]];
		double dist = dcell->distsum(p, allagg) / ndpts;
		//double dist = ptdist(dcell->data, p);
		/*VERBOSE(*verb_out << "nearest pt is [ ");
		for (size_t i = 0; i < NDIMS; ++i) {
			double d = dcell->data.x[i] - p.x[i];
			VERBOSE(*verb_out << dcell->data.x[i] << ' ');
			dist += d*d;
		}
		dist = std::sqrt(dist);* /
		VERBOSE(*verb_out << "]\n");

		//unsigned nallagg = aggcell->n[p.labels[allagg]];
		unsigned nagg = cell->n[p.labels[agg]];
		//double vol = aggcell->volume(p);
		//double D = std::max(3.0, dcell->diag(p));
		//double D = dcell->diag(p);
		//double pnil = dist / std::pow(vol, 1.0 / NDIMS);
		//double pnil = std::min(1.0, dist / D);
		//pnil *= vol / nallagg;
		//pnil = pnil / (pnil + 1.0);
		//double invdensity = vol;
		//invdensity = std::pow(invdensity, 1.0/NDIMS)/nallagg;
		//invdensity = vol / nallagg;
		double pnil = dist / (dist+3);
		//double pnil = dist / cell->diag(p);
		//pnil *= invdensity / (invdensity + 1);
		//pnil = std::sqrt(pnil);
		//pnil = std::max(pnil, 0.1);
		VERBOSE(
			*verb_out << "nagg: " << nagg
			          //<< "; nallagg: " << nallagg
			          //<< "; volume: " << vol
					  //<< "; diag: " << D
			          //<< "; invdensity: " << invdensity
			          << "; dist: " << dist << " on " << ndpts << " points"
			          << "; pnil: " << pnil
			          << std::endl
		);
		pnil = std::min(1.0, pnil);
		//VERBOSE(*verb_out << "volume = " << vol << std::endl);
		for (size_t i = 0; i < p.nlabels; ++i) {
			unsigned N = cell->n[p.labels[i]];
			//output[i] = N*N / (nagg*vol);
			double plab = (double)N/nagg;
			output[i] = plab/pnil - plab;
			//if (N == nagg) return false;
			VERBOSE(
				*verb_out << " [";
				for (size_t j = 0; j < NLABELDIMS; ++j) {
					*verb_out << ' ' << p.labels[i].p[j];
				}
				*verb_out << " ] -> " << plab << " -> " << output[i] << std::endl;
			);
		}
		return true;

		/*unsigned nallagg = n[p.labels[allagg]];
		unsigned nagg = n[p.labels[agg]];
		if (nagg < thres) {
			return false;
		}
		margintree *qchild = (p.x[dsplit] < splitval) ? low : high;
		bool q = qchild->query(p, output, agg, allagg, thres);
		if (!q) {
			VERBOSE(
				*verb_out << "Query pt [ ";
				for (size_t i = 0; i < NDIMS; ++i) {
					*verb_out << p.x[i] << ' ';
				}
				*verb_out << "]\n";
				*verb_out << "Query labels [ ";
				for (size_t i = 0; i < p.nlabels; ++i) {
					*verb_out << "[ ";
					for (size_t j = 0; j < NLABELDIMS; ++j) {
						*verb_out << p.labels[i].p[j] << ' ';
					}
					*verb_out << "] ";
				}
				*verb_out << "]\n";
				*verb_out << "Freq. table: [ ";
				for (typename std::map<label_t, unsigned>::const_iterator i = n.begin(); i != n.end(); ++i) {
					for (size_t j = 0; j < NLABELDIMS; ++j) {
						*verb_out << i->first.p[j] << '/';
					}
					*verb_out << " -> " << i->second << ' ';
				}
				*verb_out << "]\n";
			);
			double vol = volume(NONE);
			double invdensity = vol;
			double pnil = invdensity / (invdensity + 1);
			VERBOSE(*verb_out << "nagg: " << nagg << "; nallagg: " << nallagg << "; volume: " << vol << "; invdensity: " << invdensity << "; pnil: " << pnil << std::endl);
			//VERBOSE(*verb_out << "volume = " << vol << std::endl);
			for (size_t i = 0; i < p.nlabels; ++i) {
				unsigned N = n[p.labels[i]];
				//output[i] = N*N / (nagg*vol);
				output[i] = ((double)N/nagg)/pnil;///nagg;
				//if (N == nagg) return false;
				VERBOSE(
					*verb_out << " [";
					for (size_t j = 0; j < NLABELDIMS; ++j) {
						*verb_out << ' ' << p.labels[i].p[j];
					}
					*verb_out << " ] -> " << output[i] << std::endl;
				);
			}
		}
		return true;* /
	}*/

	enum {
		NONE, LO, HI
	};

	double volume(const pt_t &p) const {
		VERBOSE2(*verb_out << " calc volume:\n");
		double v = 1;
		for (size_t i = 0; i < NDIMS; ++i) {
			float lo, hi;
			lo = std::min(dmin[i], p.x[i]);
			hi = std::max(dmax[i], p.x[i]);
			v *= std::max(0.01f, hi-lo);
			VERBOSE2(*verb_out << "   dim " << i << ": " << hi << " - " << lo << " = " << hi-lo << std::endl);
		}
		VERBOSE2(*verb_out << "  = " << v << std::endl);
		return v;
	}

	double diag(const pt_t &p) const {
		VERBOSE2(*verb_out << " calc diag:\n");
		double D = 0;
		for (size_t i = 0; i < NDIMS; ++i) {
			float lo, hi;
			lo = std::min(dmin[i], p.x[i]);
			hi = std::max(dmax[i], p.x[i]);
			float d = std::max(0.01f, hi-lo);
			D += d*d;
			VERBOSE2(*verb_out << "   dim " << i << ": " << hi << " - " << lo << " = " << hi-lo << std::endl);
		}
		D = std::sqrt(D);
		VERBOSE2(*verb_out << "  = " << D << std::endl);
		return D;
	}

	/*unsigned density(const pt_t &p, unsigned *labn, unsigned *total, double *vol) {
		std::cout << "In mtree " << this << " there are " << n[p.labels[0]] << " things\n";
		std::cout << " low is " << low << " and high is " << high << std::endl;
		*vol = 0;
		if (n[p.labels[0]] == 0) {
			return 0;
		}
		margintree *qchild = (p.x[dsplit] < splitval) ? low : high;
		*labn = qchild ? qchild->density(p, labn, total, vol) : 0;
		if (*labn == 0) {
			*labn = n[p.labels[0]];
		}
		else if (*vol == 0) {
			*vol = volume(p.x[dsplit] < splitval ? LO : HI);
		}
		*total = n[p.labels[0]];
		return *labn;
	}*/

	/*
	void insert(float *p) {
		if (!low) { // ==> !high
			float maxmargin = 0.f;
			for (size_t i = 0; i < NDIMS; ++i) {
				if (p[i] < data[i]) {
					margin[i][0] = p[i];
					margin[i][1] = data[i];
				}
				else {
					margin[i][0] = data[i];
					margin[i][1] = p[i];
				}
				if (margin[i][1] - margin[i][0] > maxmargin) {
					maxmargin = margin[i][1] - margin[i][0];
					dsplit = i;
				}
			}
			if (p[dsplit] < data[dsplit]) {
				low = new margintree(p);
				high = new margintree(data);
			}
			else {
				low = new margintree(data);
				high = new margintree(p);
			}
		}
		else {
			size_t origsplit = dsplit;
			update_margin(p, dsplit);
			float splitmargin = margin[dsplit][1] - margin[dsplit][0];
			for (size_t i = 0; i < NDIMS; ++i) {
				update_margin(p, i);
				if (margin[i][1] - margin[i][0] > splitmargin) {
					splitmargin = margin[i][1] - margin[i][0];
					dsplit = i;
				}
			}
			if (dsplit != origsplit) {
			}
			else {
				float split = 0.5f * (margin[dsplit][0] + margin[dsplit][1]);
				if (p[dsplit] < split) {
					low->insert(p);
				}
				else {
					high->insert(p);
				}
			}
		}
	}*/

	static margintree *read(reader &, int *);

private:
	margintree<NDIMS,NLABELDIMS> *low;
	margintree<NDIMS,NLABELDIMS> *high;
	pt_t data;
	size_t dsplit;
	float splitval;
	std::map<label_t, unsigned> n;
	float dmin[NDIMS];
	float dmax[NDIMS];

private:
	template <unsigned ND, unsigned NLD>
	friend margintree<ND, NLD> *mkmargintree(pt<ND,NLD> *, size_t, pt<ND,NLD> **);

};


template <unsigned NDIMS, unsigned NLABELDIMS>
margintree<NDIMS, NLABELDIMS> *
mkmargintree(pt<NDIMS,NLABELDIMS> *p, size_t n, pt<NDIMS,NLABELDIMS> **buf) {
	assert(n > 0);
	margintree<NDIMS,NLABELDIMS> *tree;
	if (n == 1) {
		tree = new margintree<NDIMS,NLABELDIMS>(p[0]);
		for (size_t j = 0; j < p[0].nlabels; ++j) {
			++tree->n[p[0].labels[j]];
		}
		for (size_t j = 0; j < NDIMS; ++j) {
			tree->dmin[j] = tree->dmax[j] = p[0].x[j];
		}
		return tree;
	}

	tree = new margintree<NDIMS,NLABELDIMS>;

	for (size_t i = 0; i < n; ++i) {
		const pt<NDIMS,NLABELDIMS> &pp = p[i];
		for (size_t j = 0; j < pp.nlabels; ++j) {
			++tree->n[pp.labels[j]];
		}
		for (size_t j = 0; j < NDIMS; ++j) {
			tree->dmin[j] = std::min(tree->dmin[j], pp.x[j]);
			tree->dmax[j] = std::max(tree->dmax[j], pp.x[j]);
		}
	}

	for (size_t j = 0; j < NDIMS; ++j) {
		assert(tree->dmax[j] >= tree->dmin[j]);
		assert(tree->dmax[j] != std::numeric_limits<double>::infinity() && tree->dmax[j] != -std::numeric_limits<double>::infinity());
		assert(tree->dmin[j] != std::numeric_limits<double>::infinity() && tree->dmin[j] != -std::numeric_limits<double>::infinity());
	}

	/*
	std::cout << "my counters:\n";
	for (typename std::map<label<NLABELDIMS>, unsigned>::const_iterator i = tree->n.begin(); i != tree->n.end(); ++i) {
		std::cout << " ";
		for (size_t j = 0; j < NLABELDIMS; ++j) {
			std::cout << " " << i->first.p[j];
		}
		std::cout << ": " << i->second << std::endl;
	}*/

	float maxmargin = 0.f;
	tree->dsplit = 0;
	size_t splitpt = 0;
	for (size_t i = 0; i < NDIMS; ++i) {
		pt<NDIMS,NLABELDIMS> *dimbuf = buf[i];
		memcpy(buf[i], p, sizeof(*p)*n);
		std::sort(dimbuf, dimbuf + n, (typename pt<NDIMS,NLABELDIMS>::cmp)(i));

		for (size_t j = 0; j < n-1; ++j) {
			float marg = dimbuf[j+1].x[i] - dimbuf[j].x[i];
			marg *= splitweight((double)j/n);
			if (marg > maxmargin) {
				maxmargin = marg;
				tree->dsplit = i;
				tree->splitval = dimbuf[j].x[i] + marg/2.f;
				splitpt = j + 1;
			}
		}
	}

	//tree->dsplit = splitdim + 1;
	//if (tree->dsplit == NDIMS) tree->dsplit = 0;
	//splitpt = n/2;
	if (splitpt == 0) splitpt = n/2; // handle case where all points are identical

	//std::cout << "split on " << tree->dsplit << " @ pos " << splitpt << std::endl;
	tree->low = mkmargintree(buf[tree->dsplit], splitpt, buf);
	pt<NDIMS,NLABELDIMS> *hibuf[NDIMS];
	for (size_t i = 0; i < NDIMS; ++i) {
		hibuf[i] = buf[i] + splitpt;
	}
	tree->high = mkmargintree(buf[tree->dsplit] + splitpt, n - splitpt, hibuf);

	return tree;
}


template <unsigned NDIMS, unsigned NLABELDIMS>
margintree<NDIMS, NLABELDIMS> *
margintree<NDIMS, NLABELDIMS>::read(reader &re, int *e) {
	long dum;
	re.read(dum);
	if (dum != NDIMS) {
		if (e) *e = E_INVALID;
		return 0;
	}
	re.read(dum);
	if (dum != NLABELDIMS) {
		if (e) *e = E_INVALID;
		return 0;
	}
	margintree<NDIMS,NLABELDIMS> *mt = new margintree<NDIMS,NLABELDIMS>;
	int err = mt->read(re);
	if (e) *e = err;
	return mt;
}

template <unsigned NDIMS, unsigned NLABELDIMS>
int
margintree<NDIMS, NLABELDIMS>::read(reader &re) {
	for (size_t i = 0; i < NDIMS; ++i) {
		double dum;
		re.read(dum);
		dmin[i] = dum;
		re.read(dum);
		dmax[i] = dum;
	}
	re.read(dsplit);
	if (dsplit == NDIMS) {
		for (size_t j = 0; j < NDIMS; ++j) {
			double dum;
			re.read(dum);
			data.x[j] = (float)dum;
		}
		size_t nlabels;
		re.read(nlabels);
		while (nlabels--) {
			label_t L;
			for (size_t j = 0; j < NLABELDIMS; ++j) {
				long dum;
				re.read(dum);
				L.p[j] = dum;
			}
			long count;
			re.read(count);
			n[L] = count;
		}
		return re.isok() ? 0 : E_IO;
	}
	else {
		double dum;
		re.read(dum);
		splitval = (float)dum;
		int e;
		low = margintree<NDIMS,NLABELDIMS>::read(re, &e);
		if (FAILURE(e)) return e;
		high = margintree<NDIMS,NLABELDIMS>::read(re, &e);
		if (FAILURE(e)) return e;
		n = low->n;
		for (typename std::map<label_t, unsigned>::const_iterator j = high->n.begin(); j != high->n.end(); ++j) {
			n[j->first] += j->second;
		}
		return e;
	}
}


}

#endif
