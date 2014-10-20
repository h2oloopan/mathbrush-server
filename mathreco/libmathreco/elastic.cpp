#include "elastic.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstddef>

#include "algorithm.h"
#include "error.h"
#include "group.h"
#include "parms.h"
#include "stroke.h"
#include "stkutils.h"
#include "verb.h"
#include "ink-io.h"


namespace scg
{


static double NormScale = RegisterParameterDouble("NormalizationScale", &NormScale);
static double PositionWeight = RegisterParameterDouble("PositionWeight", &PositionWeight);

void calc_pos_weight() {
	PositionWeight /= NormScale;
}

static int _e = RegisterParameterCallback(&calc_pos_weight);

static unsigned model_sz = 32;
static unsigned input_sz = 32;
static double *model_angles = new double[model_sz];
static double *input_angles = new double[input_sz];

void
elastic_shutdown()
{
	delete[] model_angles;
	delete[] input_angles;
}


struct stroke_cursor {
	double *x;
	double *y;
	double *t;

	stroke_cursor(double *x_, double *y_, double *t_) : x(x_), y(y_), t(t_) { }

	inline stroke_cursor &operator++()
	{
		++x, ++y, ++t;
		return *this;
	}

	inline stroke_cursor operator++(int)
	{
		stroke_cursor tmp = *this;
		++x, ++y, ++t;
		return tmp;
	}

	inline stroke_cursor &operator--()
	{
		--x, --y, --t;
		return *this;
	}

	inline stroke_cursor operator--(int)
	{
		stroke_cursor tmp = *this;
		--x, --y, --t;
		return tmp;
	}

	inline stroke_cursor &operator=(const stroke_cursor &rhs)
		{ x = rhs.x, y = rhs.y, t = rhs.t; return *this; }

	inline stroke_cursor &operator+=(size_t n)
		{ x += n, y += n, t += n; return *this; }
	inline stroke_cursor &operator-=(size_t n)
		{ x -= n, y -= n, t -= n; return *this; }

	inline stroke_cursor operator+(size_t n) const
		{ return stroke_cursor(x + n, y + n, t + n); }
	inline stroke_cursor operator-(size_t n) const
		{ return stroke_cursor(x - n, y - n, t - n); }

	inline ptrdiff_t operator-(const stroke_cursor &rhs) const
		{ return x - rhs.x; }

	inline bool operator==(const stroke_cursor &rhs) const
		{ return x == rhs.x && y == rhs.y && t == rhs.t; }
	inline bool operator!=(const stroke_cursor &rhs) const
		{ return !(*this == rhs); }
	inline bool operator<(const stroke_cursor &rhs) const
		{ return x < rhs.x || y < rhs.y || t < rhs.t; }
	inline bool operator>(const stroke_cursor &rhs) const
		{ return x > rhs.x || y > rhs.y || t > rhs.t; }
	inline bool operator<=(const stroke_cursor &rhs) const
		{ return x <= rhs.x || y <= rhs.y || t <= rhs.t; }
	inline bool operator>=(const stroke_cursor &rhs) const
		{ return x >= rhs.x || y >= rhs.y || t >= rhs.t; }
};


static double
anglediff(double a, double b) {
	a -= b;
	while (a > M_PI) {
		a -= 2*M_PI;
	}
	while (a < -M_PI) {
		a += 2*M_PI;
	}
	return a;
}

static double
elastic_dist(const stroke_cursor &lhs, const stroke_cursor &rhs) {
	double angle_cost = std::abs(anglediff(*lhs.t, *rhs.t));
	//double angle_cost = std::abs(*lhs.t - *rhs.t) * (180/M_PI);
	double dx_cost = std::abs(*lhs.x - *rhs.x);
	double dy_cost = std::abs(*lhs.y - *rhs.y);
	//return 0.005 * (std::min(angle_cost, 360-angle_cost) + 0.5 * 180 * (dx_cost + dy_cost));
	double d = 1/M_PI*angle_cost + 0.5*(dx_cost + dy_cost);
	//double d = dx_cost + dy_cost;
	//assert(d <= 2);
	return d;
}


static int
match_stroke(const NormalizedStroke &model_stroke, const NormalizedStroke &input_stroke, const Rect<double> &model_bounds, const Rect<double> &input_bounds, double &stroke_cost/*,
             const double *weights, unsigned **choices*/)
{
	while (model_stroke.npoints >= model_sz) {
		delete[] model_angles;
		model_sz *= 2;
		model_angles = new double[model_sz];
	}
	while (input_stroke.npoints >= input_sz) {
		delete[] input_angles;
		input_sz *= 2;
		input_angles = new double[input_sz];
	}

	tangents_along_stroke(model_stroke, model_angles);
	tangents_along_stroke(input_stroke, input_angles);
	
		stroke_cursor model_front(model_stroke.x, model_stroke.y, model_angles);
		stroke_cursor input_start(input_stroke.x, input_stroke.y, input_angles);
		stroke_cursor model_back = model_front + model_stroke.npoints - 1;
		stroke_cursor input_end = input_start + input_stroke.npoints - 1;

		stroke_cursor input_front = input_start;
		stroke_cursor input_back = input_end;
		if (stroke_is_closed(model_stroke, model_bounds) && stroke_is_closed(input_stroke, input_bounds)) {
			double best_cost = elastic_dist(input_front, model_front);
			stroke_cursor best_start = input_front;
			++input_front;
			while (input_front <= input_end) {
				double cost = elastic_dist(input_front, model_front);
				if (cost < best_cost) {
					best_cost = cost;
					best_start = input_front;
				}
				++input_front;
			}
			input_front = best_start;
			if (input_front == input_start) {
				input_back = input_end;
			}
			else {
				input_back = best_start - 1;
			}
		}

		/*VERBOSE(
			if (input_front != input_start) {
				*verb_out << "  BEGinning match at offset " << input_front - input_start << std::endl;
			}
		);*/

		double fcost0 = elastic_dist(input_front, model_front);
		double bcost0 = elastic_dist(input_back, model_back);

		stroke_cost = fcost0 + bcost0;

		++input_front;
		if (input_front > input_end) {
			input_front = input_start;
		}
		--input_back;
		if (input_back < input_start) {
			input_back = input_end;
		}

		fcost0 = elastic_dist(input_front, model_front);
		bcost0 = elastic_dist(input_back, model_back);

		while (input_front != input_back) {
			ptrdiff_t model_rem = model_back - model_front;

			if (model_rem >= 1) {
				double cost1 = elastic_dist(input_front, model_front + 1);
				if (cost1 < fcost0) {
					if (model_rem >= 2) {
						double cost2 = elastic_dist(input_front, model_front + 2);
						if (cost2 < cost1) {
							model_front += 2;
							stroke_cost += cost2;
							fcost0 = cost2;
						}
						else {
							++model_front;
							stroke_cost += cost1;
							fcost0 = cost1;
						}
					}
					else {
						++model_front;
						stroke_cost += cost1;
						fcost0 = cost1;
					}
				}
				else {
					++model_front;
					stroke_cost += cost1;
					fcost0 = cost1;
				}
			}
			else {
				++model_front;
				stroke_cost += fcost0;
				//fcost0 = cost1;
			}

			if (model_front > model_back) {
				model_front = model_back;
			}

			++input_front;
			if (input_front > input_end) {
				input_front = input_start;
			}

			if (input_front != input_back) {
				if (model_rem >= 1) {
					double cost1 = elastic_dist(input_back, model_back - 1);
					if (cost1 < bcost0) {
						if (model_rem >= 2) {
							double cost2 = elastic_dist(input_back, model_back - 2);
							if (cost2 < cost1) {
								model_back -= 2;
								stroke_cost += cost2;
								bcost0 = cost2;
							}
							else {
								--model_back;
								stroke_cost += cost1;
								bcost0 = cost1;
							}
						}
						else {
							--model_back;
							stroke_cost += cost1;
							bcost0 = cost1;
						}
					}
					else {
						--model_back;
						stroke_cost += cost1;
						bcost0 = cost1;
					}
				}
				else {
					--model_back;
					stroke_cost += bcost0;
					//bcost0 = cost1;
				}

				--input_back;
				if (input_back < input_start) {
					input_back = input_end;
				}
			}

			if (model_front > model_back) {
				model_back = model_front;
			}
		}

		while (model_front < model_back) {
			//stroke_cost += 1.5 * elastic_dist(input_front, model_front);
			stroke_cost += elastic_dist(input_front, model_front);
			//++model_front;
			model_front += 2;
		}
    
    //stroke_cost /= num_points(input_stroke);
    
    return 0;
}

static double
do_elastic_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder) {
	int e;
	double total_cost = 1.0;

	Rect<double> model_bounds = model.bounds();
	Rect<double> input_bounds = input.bounds();

	for (size_t i = 0; i < model.nstrokes; ++i) {
		double stroke_cost = 0.0;
		const NormalizedStroke &mod_stroke = model.strokes[i];
		const NormalizedStroke &in_stroke = input.strokes[inputorder[i]];
		if (!(stroke_is_dot(mod_stroke, model_bounds) && stroke_is_dot(in_stroke, input_bounds))) {
			e = match_stroke(mod_stroke, in_stroke, model_bounds, input_bounds, stroke_cost);
			stroke_cost /= in_stroke.npoints;//std::max(in_stroke.npoints, mod_stroke.npoints);
			//assert(stroke_cost <= 2);
			if (FAILURE(e)) {
				return e;
			}
		}
		else {
			stroke_cost = std::abs(mod_stroke.x_center() - in_stroke.x_center()) + std::abs(mod_stroke.y_center() - in_stroke.y_center());
		}
        
		total_cost *= stroke_cost;
	}

	return std::pow(total_cost, 1.0/input.nstrokes);//total_cost / num_points(input);
}


double
elasticdist(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder) {
	double forward_score = do_elastic_match(model, input, inputorder);
	std::vector<size_t> modelorder(inputorder.size());
	for (size_t i = 0; i < inputorder.size(); ++i) {
		modelorder[inputorder[i]] = i;
	}
	double backward_score = do_elastic_match(input, model, modelorder);
	/*static bool ref = false;
	if (!ref) {
		std::cout << "model: " << model << std::endl;
		std::cout << "input: " << input << std::endl;
		std::cout << "score: " << std::sqrt(forward_score * backward_score) << std::endl;
	}*/
	return std::sqrt(forward_score * backward_score);
}


}
