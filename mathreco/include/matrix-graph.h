#ifndef SCG_MATRIX_GRAPH_H
#define SCG_MATRIX_GRAPH_H

#include <vector>
#include <map>
#include "matrix.h"
#include "RangeExpander.h"

namespace scg
{


/* Defines a graph whose nodes/vertices are containers of MatrixElement instances and edges connect two vertices
   such that the weight of the edge is the distance between the vertices (i.e. the number of elements between the two).
*/
class Graph
{
public:
	struct Edge;

	struct Node {
		Node();
		Node(MatrixElement* elem);

		MatrixElement* elem;

		bool operator==(const Node& n) const;
		bool operator!=(const Node& n) const;
	};

	struct Edge {
		enum Type {HORIZONTAL, VERTICAL, DIAGONAL, ANTIDIAGONAL };

		Edge() {}
		Edge(Node n1, Node n2, Type type, RangeExpander expander);

		Node n1, n2;
		Type type;
		
		Node& otherNode(Node refNode);
		size_t getWeight() const;
		void setWeight(size_t weight);

		const ExpressionTree* expand(size_t i);

		bool operator==(const Edge& e2) const;

	private:
		size_t weight; //used if we have a 0-weighted edge (e.g "1 ... 1") and set externally
		RangeExpander expander;
	};

private:
	std::vector<Node> nodes;
	std::map<Node, std::vector<Edge> > edges;

public:
	// Have to manually copy vector contents or else it fails..
	Graph& operator=(const Graph& other) {
		if (this != &other) {
			nodes.clear();
			for (std::vector<Node>::const_iterator it = other.nodes.begin(); it != other.nodes.end(); it++) {
				nodes.push_back(*it);
			}

			edges = other.edges;
		}
		return *this;
	}

	void addNode(Node node);
	void removeNode(Node node);
	std::vector<Node> getNodes() const;

	void addEdge(Node& sNode, Node& eNode, Edge::Type type, RangeExpander expander);
	std::vector<Edge> getEdges(Node& n) const;
	Edge getEdge(Node n1, Node n2) const;
	void setEdgeWeight(Node s, Node e, size_t weight);

	bool isEmpty() const;

	Node* elemToNode(MatrixElement* elem);
	MatrixElement* nodeToElem(Node& node);

	bool operator==(const Graph& other);
	bool operator!=(const Graph& other) { return !operator==(other); }
}; //Graph

//for mapping
bool operator<(const Graph::Node &n1, const Graph::Node &n2);
bool operator<(const Graph::Edge &e1, const Graph::Edge &e2);


}

#endif