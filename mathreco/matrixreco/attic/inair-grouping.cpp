#include "grouping-strat.h"

#include <algorithm>
#include "rect.h"
#include "interval.h"
#include "error.h"

namespace scg
{



const double InAirGrouping::ROW_BREAK_TOL = 0.5;
const double InAirGrouping::COL_BREAK_TOL = 0.5;
const double InAirGrouping::ALPHA_ROW = 1.0;
const double InAirGrouping::ALPHA_COL = 0.5;

// hybrid constants
const double InAirGrouping::MIN_OVERLAPPING_PROPORTION = 0.25;
const double InAirGrouping::MIN_CHAR_HEIGHT_INCHES = 0.06;
long InAirGrouping::MIN_CHAR_HEIGHT = 0;
const double InAirGrouping::SPAN_OVERLAP_ON_ELEM_UPDATE_TOL = 0.5;
const double InAirGrouping::ALPHA_DIST = 1.0;

const double InAirGrouping::Span::BUFFER_PCT = 0.3;

InAirGrouping::InAirGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices) : ElemGroupingStrat(parser, matrices)
{
}

int InAirGrouping::addStroke(const RawStroke* stroke)
{
	strokes.push_back(stroke);
	analyse(strokes);
	return NO_ERROR;
}

int InAirGrouping::addStrokes(const RawStrokeGroup& strokeGroup)
{
	for (RawStrokeGroup::const_iterator s = strokeGroup.begin(); s != strokeGroup.end(); ++s) {
		const RawStroke* stroke = &*s;
		strokes.push_back(stroke);
	}

	analyse(strokes);
	return NO_ERROR;
}

int InAirGrouping::removeStroke(const RawStroke* stroke)
{
	//TODO
	return NO_ERROR;
}

void InAirGrouping::forceGrouping(const RawStroke** strokes, size_t nstrokes)
{
	//TODO
}

int InAirGrouping::getConfidence(const Matrix* m, double& confidence) const
{
	//TODO how to calculate?
	confidence = 1.0;
	return NO_ERROR;
}

void InAirGrouping::analyse(const std::vector<const RawStroke*> &strokes)
{
	if (strokes.size() == 0) return;

	Matrix* m = matrices->front();
	
	// TODO make more iterative
	m->clear();

	Rect<long> matrixDim = bbox(**strokes.begin());
	for (std::vector<const RawStroke*>::const_iterator it = strokes.begin(); it != strokes.end(); it++) {
		matrixDim = merge(matrixDim, bbox(**it));
	}
	long matrixWidth = width(matrixDim);
	long matrixHeight = height(matrixDim);

	// Calculate pen-down and in-air stroke widths and heights
	std::vector<long> penDownWidths = getStrokeWidths(strokes);
	std::vector<long> penDownHeights = getStrokeHeights(strokes);
	std::vector<long> inAirWidths = getInAirStrokeWidths(strokes);
	std::vector<long> inAirHeights = getInAirStrokeHeights(strokes);

	// Determine if it is row-major, col-major, or hybrid
	std::vector<BreakPoint> rowBreakIndices = findMajorRowBreaks(inAirWidths, matrixWidth);
	std::vector<BreakPoint> colBreakIndices = findMajorColBreaks(inAirHeights, matrixHeight);

	std::vector<BreakPoint> majorBreaks;
	std::vector<BreakPoint> minorBreaks;
	if (rowBreakIndices.size() > 0 && colBreakIndices.size() == 0) {
		// row-major
		majorBreaks = rowBreakIndices;
		minorBreaks = findMinorColBreaks(inAirWidths, penDownWidths, rowBreakIndices);

		size_t majIdx = 0;
		size_t minIdx = 0;
		size_t row = 0;
		size_t col = 0;
		for (size_t i = 0; i < strokes.size(); /*increased below*/) {
			size_t nextBreak = 0;
			bool isRowBreak;
			if (majIdx == majorBreaks.size() && minIdx == minorBreaks.size()) {
				nextBreak = strokes.size() - 1;
				isRowBreak = false;
			} else if (majIdx == majorBreaks.size() && minIdx < minorBreaks.size()) {
				nextBreak = minorBreaks[minIdx++].inAirIndex;
				isRowBreak = false;
			} else if (majIdx < majorBreaks.size() && minIdx == minorBreaks.size()) {
				nextBreak = majorBreaks[majIdx++].inAirIndex;
				isRowBreak = true;
			} else {
				if (majorBreaks[majIdx].inAirIndex < minorBreaks[minIdx].inAirIndex) {
					nextBreak = majorBreaks[majIdx++].inAirIndex;
					isRowBreak = true;
				} else {
					nextBreak = minorBreaks[minIdx++].inAirIndex;
					isRowBreak = false;
				}
			}

			MatrixElement* elem = new MatrixElement();
			for (; i <= nextBreak; i++) {
				elem->addStroke(strokes[i]);
			}
			m->insert(*elem, row, col);

			if (isRowBreak) {
				row++;
				col = 0;
			} else {
				col++;
			}
		}
	} else if (rowBreakIndices.size() == 0 && colBreakIndices.size() > 0) {
		// col-major
		majorBreaks = colBreakIndices;
		minorBreaks = findMinorRowBreaks(inAirHeights, penDownHeights, colBreakIndices);
		
		size_t majIdx = 0;
		size_t minIdx = 0;
		size_t row = 0;
		size_t col = 0;
		for (size_t i = 0; i < strokes.size(); /*increased below*/) {
			size_t nextBreak = 0;
			bool isColBreak;
			if (majIdx == majorBreaks.size() && minIdx == minorBreaks.size()) {
				nextBreak = strokes.size() - 1;
				isColBreak = false;
			} else if (majIdx == majorBreaks.size() && minIdx < minorBreaks.size()) {
				nextBreak = minorBreaks[minIdx++].inAirIndex;
				isColBreak = false;
			} else if (majIdx < majorBreaks.size() && minIdx == minorBreaks.size()) {
				nextBreak = majorBreaks[majIdx++].inAirIndex;
				isColBreak = true;
			} else {
				if (majorBreaks[majIdx].inAirIndex < minorBreaks[minIdx].inAirIndex) {
					nextBreak = majorBreaks[majIdx++].inAirIndex;
					isColBreak = true;
				} else {
					nextBreak = minorBreaks[minIdx++].inAirIndex;
					isColBreak = false;
				}
			}

			MatrixElement* elem = new MatrixElement();
			for (; i <= nextBreak; i++) {
				elem->addStroke(strokes[i]);
			}
			m->insert(*elem, row, col);

			if (isColBreak) {
				col++;
				row = 0;
			} else {
				row++;
			}
		}
	} else {
		// hybrid - use euclidean distances
		rowSpan.clear();
		colSpan.clear();

		if (inAirWidths.size() == 0) return;

		std::vector<long> inAirLengths, penDownLengths;
		for (size_t i = 0; i < inAirWidths.size(); i++) {
			long delta_x = inAirWidths[i];
			long delta_y = inAirHeights[i];
			inAirLengths.push_back((long)std::sqrt((double)(delta_x * delta_x) + (delta_y * delta_y)));
		}
		for (size_t i = 0; i < penDownWidths.size(); i++) {
			long delta_x = penDownWidths[i];
			long delta_y = penDownHeights[i];
			penDownLengths.push_back((long)std::sqrt((double)(delta_x * delta_x) + (delta_y * delta_y)));
		}

		minorBreaks = findMinorHybridBreaks(inAirLengths, penDownLengths);

		// for each element, build it from the strokes and insert it into the matrix according to a spatial segmentation algorithm (taken from InAirGrouping)
		size_t minIdx = 0;
		for (size_t i = 0; i < strokes.size(); /*increased below*/) {
			size_t nextBreak = (minIdx < minorBreaks.size()) ? minorBreaks[minIdx++].inAirIndex : strokes.size() - 1;

			MatrixElement* elem = new MatrixElement();
			for (; i <= nextBreak; i++) {
				elem->addStroke(strokes[i]);
			}
			

			// If this is the first element, default to it being the top-left
			if (m->isEmpty()) {
				m->insert(*elem, 0, 0);
				
				Span rSpan(elem, Span::HEIGHT), cSpan(elem, Span::WIDTH);
				rowSpan.push_back(rSpan);
				colSpan.push_back(cSpan);
			} else {
				size_t row, col;
				bool newRow, newCol;
				determineRowCol(*elem, row, col, newRow, newCol);
				m->insert(*elem, row, col, newRow, newCol);
			}
		}
	}
}

std::vector<long> InAirGrouping::getInAirStrokeWidths(const std::vector<const RawStroke*> &strokes) const
{
	std::vector<long> widths;

	if (strokes.size() == 0) return widths;

	const RawStroke* prev = *strokes.begin();
	for (std::vector<const RawStroke*>::const_iterator it = strokes.begin() + 1; it != strokes.end(); it++) {
		const RawStroke* cur = *it;
		widths.push_back(cur->x[0] - prev->x[prev->npoints - 1]);
		prev = cur;
	}

	return widths;
}

std::vector<long> InAirGrouping::getInAirStrokeHeights(const std::vector<const RawStroke*> &strokes) const
{
	std::vector<long> heights;

	if (strokes.size() == 0) return heights;

	const RawStroke* prev = *strokes.begin();
	for (std::vector<const RawStroke*>::const_iterator it = strokes.begin() + 1; it != strokes.end(); it++) {
		const RawStroke* cur = *it;
		heights.push_back(cur->y[0] - prev->y[prev->npoints - 1]);
		prev = cur;
	}

	return heights;
}

std::vector<long> InAirGrouping::getStrokeWidths(const std::vector<const RawStroke*> &strokes) const
{
	std::vector<long> widths;

	if (strokes.size() == 0) return widths;

	for (std::vector<const RawStroke*>::const_iterator it = strokes.begin(); it != strokes.end(); it++) {
		const RawStroke* stk = *it;
		widths.push_back(stk->x[0] - stk->x[stk->npoints - 1]);
	}

	return widths;
}

std::vector<long> InAirGrouping::getStrokeHeights(const std::vector<const RawStroke*> &strokes) const
{
	std::vector<long> heights;

	if (strokes.size() == 0) return heights;

	for (std::vector<const RawStroke*>::const_iterator it = strokes.begin(); it != strokes.end(); it++) {
		const RawStroke* stk = *it;
		heights.push_back(stk->y[0] - stk->y[stk->npoints - 1]);
	}

	return heights;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMajorRowBreaks(const std::vector<long> &inAirStrokeWidths, long matrixWidth) const
{
	std::vector<BreakPoint> breaks;
	const double THRESHOLD = ROW_BREAK_TOL * matrixWidth;

	for (size_t i = 0; i < inAirStrokeWidths.size(); i++) {
		long inAirWidth = inAirStrokeWidths[i];
		if (inAirWidth < 0 && std::abs(inAirWidth) > THRESHOLD) {
			breaks.push_back(BreakPoint(i));
		}
	}

	return breaks;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMajorColBreaks(const std::vector<long> &inAirStrokeHeights, long matrixHeight) const
{
	std::vector<BreakPoint> breaks;
	const double THRESHOLD = COL_BREAK_TOL * matrixHeight;

	for (size_t i = 0; i < inAirStrokeHeights.size(); i++) {
		long inAirHeight = inAirStrokeHeights[i];
		if (inAirHeight < 0 && std::abs(inAirHeight) > THRESHOLD) {
			breaks.push_back(BreakPoint(i));
		}
	}

	return breaks;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMinorColBreaks(const std::vector<long> &inAirStrokeWidths, const std::vector<long> &penDownStrokeWidths, const std::vector<InAirGrouping::BreakPoint> &rowBreakPoints) const
{
	return findMinorBreaks(inAirStrokeWidths, penDownStrokeWidths, rowBreakPoints, ROW); 
	//std::vector<BreakPoint> minorBreaks;

	//std::vector<int> allStrokeWidths(inAirStrokeWidths);
	//allStrokeWidths.insert(allStrokeWidths.end(), penDownStrokeWidths.begin(), penDownStrokeWidths.end());

	//int numMajorSegments = rowBreakPoints.size() + 1;
	//int* numberOfSegments = new int[numMajorSegments];
	//for (int i = 0; i < numMajorSegments; i++) {
	//	numberOfSegments[i] = 1; //start by assuming that each major segment consists of only 1 minor segment
	//}

	//// Construct majorSegmentsInAirStrokes, which contains each major segment's in-air strokes
	//std::vector<std::vector<int> > majorSegmentsInAirStrokes;
	//for (int i = 0; i < numMajorSegments; i++) {
	//	std::vector<int>::iterator start, end;
	//	if (i == 0) {
	//		// first stroke to the first row break (exclusive)
	//		start = inAirStrokeWidths.begin();
	//		end = inAirStrokeWidths.begin() + rowBreakPoints[0].inAirIndex;
	//	} else if (i == numMajorSegments - 1) {
	//		// previous row break to next row break (exclusive)
	//		start = inAirStrokeWidths.begin() + rowBreakPoints[i - 1].inAirIndex + 1;
	//		end = inAirStrokeWidths.end();
	//	} else {
	//		// previous row break to end
	//		start = inAirStrokeWidths.begin() + rowBreakPoints[i - 1].inAirIndex + 1;
	//		end = inAirStrokeWidths.begin() + rowBreakPoints[i].inAirIndex;
	//	}

	//	majorSegmentsInAirStrokes.push_back(std::vector<int>(start, end));
	//}
	//
	//// Grubb's test for outliers
	//std::map<StatTypes, double> stats = calcStats(allStrokeWidths);
	//for (int i = 0; i < numMajorSegments; i++) {
	//	std::vector<int> rowStrokes = majorSegmentsInAirStrokes[i];

	//	for (std::vector<int>::iterator it2 = rowStrokes.begin(); it2 != rowStrokes.end(); it2++) {
	//		int strokeWidth = *it2;
	//		double alpha = std::abs(std::abs(strokeWidth) - stats[MEAN]) / stats[STANDARD_DEVIATION];
	//		if (strokeWidth > 0 && ALPHA_ROW < alpha) {
	//			numberOfSegments[i]++;
	//		}
	//	}
	//}

	//// Get the median num segments
	//std::sort(numberOfSegments, numberOfSegments + numMajorSegments);
	//int n = numberOfSegments[numMajorSegments / 2];

	//int startIdx = 0;
	//int numBreaks = n-1; //#breaks = (#minor segments) - 1
	//if (numBreaks > 0) {
	//	for (int i = 0; i < numMajorSegments; i++) {
	//		std::vector<int> segmentInAirStrokes = majorSegmentsInAirStrokes[i];
	//		std::vector<int> outlierIndices = getNLargestOutlierIndices(numBreaks, segmentInAirStrokes);
	//		for (std::vector<int>::iterator it = outlierIndices.begin(); it != outlierIndices.end(); it++) {
	//			minorBreaks.push_back(BreakPoint(startIdx + *it));
	//		}
	//		startIdx += segmentInAirStrokes.size() + 1; //+1 for the major break
	//	}
	//}

	//delete [] numberOfSegments;
	//std::sort(minorBreaks.begin(), minorBreaks.end());
	//return minorBreaks;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMinorRowBreaks(const std::vector<long> &inAirStrokeHeights, const std::vector<long> &penDownStrokeHeights, const std::vector<InAirGrouping::BreakPoint> &colBreakPoints) const
{
	return findMinorBreaks(inAirStrokeHeights, penDownStrokeHeights, colBreakPoints, COL); 
	//std::vector<BreakPoint> minorBreaks;

	//std::vector<int> allStrokeHeights(inAirStrokeHeights);
	//allStrokeHeights.insert(allStrokeHeights.end(), penDownStrokeHeights.begin(), penDownStrokeHeights.end());

	//int numMajorSegments = colBreakPoints.size() + 1;
	//int* numberOfSegments = new int[numMajorSegments];
	//for (int i = 0; i < numMajorSegments; i++) {
	//	numberOfSegments[i] = 1; //start by assuming that each major segment consists of only 1 minor segment
	//}

	//// Construct majorSegmentsInAirStrokes, which contains each row's in-air strokes
	//std::vector<std::vector<int> > majorSegmentsInAirStrokes;
	//for (int i = 0; i < numMajorSegments; i++) {
	//	std::vector<int>::iterator start, end;
	//	if (i == 0) {
	//		// first stroke to the first row break (exclusive)
	//		start = inAirStrokeHeights.begin();
	//		end = inAirStrokeHeights.begin() + colBreakPoints[0].inAirIndex;
	//	} else if (i == numMajorSegments - 1) {
	//		// previous row break to next row break (exclusive)
	//		start = inAirStrokeHeights.begin() + colBreakPoints[i - 1].inAirIndex + 1;
	//		end = inAirStrokeHeights.end();
	//	} else {
	//		// previous row break to end
	//		start = inAirStrokeHeights.begin() + colBreakPoints[i - 1].inAirIndex + 1;
	//		end = inAirStrokeHeights.begin() + colBreakPoints[i].inAirIndex;
	//	}

	//	majorSegmentsInAirStrokes.push_back(std::vector<int>(start, end));
	//}
	//
	//// Grubb's test for outliers
	//std::map<StatTypes, double> stats = calcStats(allStrokeHeights);
	//for (int i = 0; i < numMajorSegments; i++) {
	//	std::vector<int> colStrokes = majorSegmentsInAirStrokes[i];

	//	for (std::vector<int>::iterator it2 = colStrokes.begin(); it2 != colStrokes.end(); it2++) {
	//		int strokeHeight = *it2;
	//		double alpha = std::abs(std::abs(strokeHeight) - stats[MEAN]) / stats[STANDARD_DEVIATION];
	//		if (strokeHeight < 0 && ALPHA_COL < alpha) {
	//			numberOfSegments[i]++;
	//		}
	//	}
	//}

	//// Get the median num segments
	//std::sort(numberOfSegments, numberOfSegments + numMajorSegments);
	//int n = numberOfSegments[numMajorSegments / 2];

	//int startIdx = 0;
	//int numBreaks = n-1; //#breaks = (#minor segments) - 1
	//if (numBreaks > 0) {
	//	for (int i = 0; i < numMajorSegments; i++) {
	//		std::vector<int> segmentInAirStrokes = majorSegmentsInAirStrokes[i];
	//		std::vector<int> outlierIndices = getNLargestOutlierIndices(numBreaks, segmentInAirStrokes);
	//		for (std::vector<int>::iterator it = outlierIndices.begin(); it != outlierIndices.end(); it++) {
	//			minorBreaks.push_back(BreakPoint(startIdx + *it));
	//		}
	//		startIdx += segmentInAirStrokes.size() + 1; //+1 for the major break
	//	}
	//}

	//delete [] numberOfSegments;
	//std::sort(minorBreaks.begin(), minorBreaks.end());
	//return minorBreaks;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMinorBreaks(const std::vector<long> &inAirStrokeLengths, const std::vector<long> &penDownStrokeLengths, const std::vector<InAirGrouping::BreakPoint> &majorBreakPoints, InAirGrouping::MajorBreakType majorType) const
{
	std::vector<BreakPoint> minorBreaks;

	std::vector<long> allStrokeLengths(inAirStrokeLengths);
	allStrokeLengths.insert(allStrokeLengths.end(), penDownStrokeLengths.begin(), penDownStrokeLengths.end());

	size_t numMajorSegments = majorBreakPoints.size() + 1;
	size_t* numberOfSegments = new size_t[numMajorSegments];
	for (size_t i = 0; i < numMajorSegments; i++) {
		numberOfSegments[i] = 1; //start by assuming that each major segment consists of only 1 minor segment
	}

	// Construct majorSegmentsInAirStrokes, which contains each major segment's in-air strokes
	std::vector<std::vector<long> > majorSegmentsInAirStrokes;

	if (majorType == HYBRID) {
		// if it's hybrid, we treat the whole matrix as the major axis
		majorSegmentsInAirStrokes.push_back(inAirStrokeLengths);
	} else {
		for (size_t i = 0; i < numMajorSegments; i++) {
			std::vector<long>::const_iterator start, end;
			if (i == 0) {
				// first stroke to the first major break (exclusive)
				start = inAirStrokeLengths.begin();
				end = inAirStrokeLengths.begin() + majorBreakPoints[0].inAirIndex;
			} else if (i == numMajorSegments - 1) {
				// previous major break to next major break (exclusive)
				start = inAirStrokeLengths.begin() + majorBreakPoints[i - 1].inAirIndex + 1;
				end = inAirStrokeLengths.end();
			} else {
				// previous major break to end
				start = inAirStrokeLengths.begin() + majorBreakPoints[i - 1].inAirIndex + 1;
				end = inAirStrokeLengths.begin() + majorBreakPoints[i].inAirIndex;
			}

			majorSegmentsInAirStrokes.push_back(std::vector<long>(start, end));
		}
	}
	
	// Grubb's test for outliers
	std::map<StatTypes, double> stats = calcStats(allStrokeLengths);
	for (size_t i = 0; i < numMajorSegments; i++) {
		std::vector<long> inAirStrokes = majorSegmentsInAirStrokes[i];

		for (std::vector<long>::const_iterator it2 = inAirStrokes.begin(); it2 != inAirStrokes.end(); it2++) {
			long strokeLength = *it2;
			double alpha = std::abs(std::abs(strokeLength) - stats[MEAN]) / stats[STANDARD_DEVIATION];
			if (majorType == COL) {
				if (strokeLength < 0 && ALPHA_COL < alpha) {
					numberOfSegments[i]++;
				}
			} else if (majorType == ROW) {
				if (strokeLength > 0 && ALPHA_ROW < alpha) {
					numberOfSegments[i]++;
				}
			} else if (majorType == HYBRID) {
				if (ALPHA_DIST < alpha) {
					numberOfSegments[i]++;
				}
			}
		}
	}

	// Get the median num segments
	std::sort(numberOfSegments, numberOfSegments + numMajorSegments);
	size_t n = numberOfSegments[numMajorSegments / 2];

	size_t startIdx = 0;
	size_t numBreaks = n-1; //#breaks = (#minor segments) - 1
	if (numBreaks > 0) {
		for (size_t i = 0; i < numMajorSegments; i++) {
			std::vector<long> segmentInAirStrokes = majorSegmentsInAirStrokes[i];
			std::vector<size_t> outlierIndices = getNLargestOutlierIndices(numBreaks, segmentInAirStrokes);
			for (std::vector<size_t>::iterator it = outlierIndices.begin(); it != outlierIndices.end(); it++) {
				minorBreaks.push_back(BreakPoint(startIdx + *it));
			}
			startIdx += segmentInAirStrokes.size() + 1; //+1 for the major break
		}
	}

	delete [] numberOfSegments;
	std::sort(minorBreaks.begin(), minorBreaks.end());
	return minorBreaks;
}

std::vector<InAirGrouping::BreakPoint> InAirGrouping::findMinorHybridBreaks(const std::vector<long> &inAirStrokeLengths, const std::vector<long> &penDownStrokeLengths) const
{
	std::vector<BreakPoint> noMajorAxis;
	return findMinorBreaks(inAirStrokeLengths, penDownStrokeLengths, noMajorAxis, HYBRID);
}

std::map<InAirGrouping::StatTypes, double> InAirGrouping::calcStats(const std::vector<long> &data) const
{
	std::map<StatTypes, double> stats;

	// calculate mean
	double mean = 0.0;
	for (std::vector<long>::const_iterator it = data.begin(); it != data.end(); it++) {
		mean += std::abs(*it);
	}
	mean /= data.size();
	stats[MEAN] = mean;

	// calculate (sample) standard dev
	double stddev = 0.0;
	for (std::vector<long>::const_iterator it = data.begin(); it != data.end(); it++) {
		stddev += std::pow(std::abs(*it) - mean, 2);
	}
	stddev = std::sqrt(stddev / (data.size() - 1));
	stats[STANDARD_DEVIATION] = stddev;

	return stats;
}

std::vector<size_t> InAirGrouping::getNLargestOutlierIndices(size_t n, const std::vector<long> &data) const
{
	std::vector<size_t> outlierIndices;

	if (n == 0) return outlierIndices;

	for (size_t i = 0; i < data.size(); i++) {
		long val = data[i];
		if (outlierIndices.size() == 0) {
			outlierIndices.push_back(i);
		} else {
			bool inserted = false;
			for (size_t j = outlierIndices.size(); j > 0; j--) {
				size_t idx = outlierIndices[j-1];
				if (val >= data[idx]) {
					if (outlierIndices.size() == n) {
						outlierIndices.erase(outlierIndices.begin());
						outlierIndices.insert(outlierIndices.begin() + j - 1, i);
					} else {
						outlierIndices.insert(outlierIndices.begin() + j, i);
					}
					inserted = true;
					break;
				}
			}

			if (!inserted && outlierIndices.size() < n) {
				outlierIndices.insert(outlierIndices.begin(), i);
			}
		}
	}

	return outlierIndices;
}

bool InAirGrouping::determineRowCol(MatrixElement& elem, size_t& row, size_t& col, bool& newIntermRow, bool& newIntermCol)
{
	bool ret = false;
	newIntermRow = false;
	newIntermCol = false;

	const Rect<long>& rect = elem.getBoundingRect();
	Interval heightInterval(rect.top, rect.bottom);
	Interval widthInterval(rect.left, rect.right);

	// Find row
	for (size_t i = 0; i < rowSpan.size(); i++) {
		Interval rowInterval = rowSpan[i].getSpanInterval();
		double overlapProportion = intervalOverlapProportion(rowInterval, heightInterval);
		if (overlapProportion > 0) {
			if (overlapProportion > MIN_OVERLAPPING_PROPORTION) {
				// Belongs to this row
				row = i;
				break;
			} else if (heightInterval.start < rowInterval.start) {
				// Need to insert before this row
				row = i;
				newIntermRow = true;
				break;
			}
		} else if (heightInterval.end < rowInterval.start) {
			// Needs to insert before this row
			row = i;
			newIntermRow = true;
			break;
		}
		
		row = i + 1;
	}

	if (newIntermRow || row == rowSpan.size()) {
		Span rSpan(&elem, Span::HEIGHT);
		sorted_insert(rSpan, rowSpan);
		ret |= true;
	} else {
		rowSpan[row].addElement(&elem, Span::HEIGHT);
		ret |= false;
	}

	// Find column
	for (size_t i = 0; i < colSpan.size(); i++) {
		Interval colInterval = colSpan[i].getSpanInterval();
		double overlapProportion = intervalOverlapProportion(colInterval, widthInterval);
		if (overlapProportion > 0) {
			if (overlapProportion > MIN_OVERLAPPING_PROPORTION) {
				// Belongs to this column
				col = i;
				break;
			} else if (widthInterval.start < colInterval.start) {
				// Need to insert before this column
				col = i;
				newIntermCol = true;
				break;
			}
		} else if (widthInterval.end < colInterval.start) {
			// Needs to insert before this column
			col = i;
			newIntermCol = true;
			break;
		}
		
		col = i + 1;
	}

	if (newIntermCol || col == colSpan.size()) {
		Span cSpan(&elem, Span::WIDTH);
		sorted_insert(cSpan, colSpan);
		ret |= true;
	} else {
		colSpan[col].addElement(&elem, Span::WIDTH);
		ret |= false;
	}

	return ret;
}


// Span, taken from InAirGrouping
InAirGrouping::Span::Span() {}
InAirGrouping::Span::Span(MatrixElement* elem, SpanType spanType)
{
	addElement(elem, spanType);
}

void InAirGrouping::Span::addElement(MatrixElement* elem, SpanType spanType)
{
	elems.push_back(elem);
	spanTypes.push_back(spanType);
}

void InAirGrouping::Span::removeElement(MatrixElement* elem)
{
	size_t i;
	for (i = 0; i < elems.size() && elems[i] != elem; i++) {}
	if (i != elems.size()) {
		elems.erase(elems.begin() + i);
		spanTypes.erase(spanTypes.begin() + i);
	}
}

bool InAirGrouping::Span::isEmpty() const
{
	return (elems.size() == 0);
}

Interval InAirGrouping::Span::getSpanInterval() const
{
	// We must compute the interval every time this is called, in case an element has changed (added/removed strokes)
	Interval interval; //	TODO if no elements?
	bool firstStroke = true;
	double avgHeight = 0.0;
	size_t elemsInAvgHeight = 0;

	// Go over each element in the span and merge all of its strokes' intervals
	for (size_t i = 0; i < elems.size(); i++) {
		std::vector<const RawStroke*> strokes = elems[i]->getStrokes();
		long elemHeight = height(elems[i]->getBoundingRect());
		if (elemHeight > MIN_CHAR_HEIGHT) {
			avgHeight += elemHeight;
			elemsInAvgHeight++;
		}

		// TODO make better? (observer pattern on MatrixElement to observe changes?)
		
		for (std::vector<const RawStroke*>::iterator strokeIter = strokes.begin(); strokeIter != strokes.end(); strokeIter++) {
			Rect<long> rect = bbox(**strokeIter);
			
			// Construct the stroke interval
			Interval strokeInterval;
			if (spanTypes[i] == HEIGHT) {
				strokeInterval.start = rect.top;
				strokeInterval.end = rect.bottom;
			} else {
				strokeInterval.start = rect.left;
				strokeInterval.end = rect.right;
			}

			interval = (firstStroke) ? strokeInterval : mergeIntervals(interval, strokeInterval);
			firstStroke = false;
		}
	}

	if (elemsInAvgHeight > 0 && interval.length() < avgHeight) {
		avgHeight /= elemsInAvgHeight;
		interval.start -= static_cast<int>(avgHeight * BUFFER_PCT);
		interval.end += static_cast<int>(avgHeight * BUFFER_PCT);
	}

	return interval;
}


}