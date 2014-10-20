#include "matrix-graph.h"
#include <cstring>
#include <cstddef>

namespace scg
{
	
	
Graph::Node::Node() : elem(NULL) {}
Graph::Node::Node(MatrixElement* elem) : elem(elem) {}
bool Graph::Node::operator==(const Node& n) const { return elem == n.elem; }
bool Graph::Node::operator!=(const Node& n) const { return !(elem == n.elem); }


Graph::Edge::Edge(Node n1, Node n2, Type type, RangeExpander expander) : n1(n1), n2(n2), type(type), weight(0), expander(expander) {}
Graph::Node& Graph::Edge::otherNode(Node refNode) { return (refNode == n1) ? n2 : n1; }
bool Graph::Edge::operator==(const Edge& e2) const { return type == e2.type && ((n1 == e2.n1 && n2 == e2.n2) || (n1 == e2.n2 && n2 == e2.n1)); }
void Graph::Edge::setWeight(size_t weight) { this->weight = weight; }
size_t Graph::Edge::getWeight() const
{
	std::vector<const ExpressionTree *> startExprs = expander.getStartExpressions();

	if (expander.size() > 0) return expander.size() - (startExprs.size() + 1); //account for starting and ending elements

	if (std::strcmp(startExprs.back()->long_str(), expander.getEndExpression()->long_str()) == 0) //TODO long_str ok?
		return weight;
	else
		return 0;
}

const ExpressionTree* Graph::Edge::expand(size_t i)
{
	if (i >= getWeight()) return NULL;

	if (expander.size() == 0) return expander.at(0);
	else return expander.at(i + expander.getStartExpressions().size()); //+expander.size() to skip first elements (only expand the edge, not its bounding elements)
}

void Graph::addNode(Node node)
{
	nodes.push_back(node);
}

void Graph::removeNode(Node node)
{
	for (size_t i = 0; i < nodes.size(); i++) {
		if (nodes[i] == node) {
			nodes.erase(nodes.begin() + i);
			break;
		}
	}
}

std::vector<Graph::Node> Graph::getNodes() const
{
	return nodes;
}

void Graph::addEdge(Node& sNode, Node& eNode, Edge::Type type, RangeExpander expander)
{
	Edge e(sNode, eNode, type, expander);
	edges[sNode].push_back(e);
	edges[eNode].push_back(e);
}

std::vector<Graph::Edge> Graph::getEdges(Graph::Node& n) const
{
	std::map<Node, std::vector<Edge> >::const_iterator edgesIter = edges.find(n);
	if (edgesIter != edges.end()) return edgesIter->second;
	else return std::vector<Edge>();
}

Graph::Edge Graph::getEdge(Node n1, Node n2) const
{
	std::vector<Edge> n1_edges = edges.find(n1)->second;
	for (std::vector<Edge>::const_iterator it = n1_edges.begin(); it != n1_edges.end(); it++) {
		Edge edge = *it;
		if (edge.otherNode(n1) == n2) return edge;
	}
	return Edge();
}

void Graph::setEdgeWeight(Node s, Node e, size_t weight)
{
	std::vector<Edge>& s_edges = (edges.find(s))->second;
	for (std::vector<Edge>::iterator it = s_edges.begin(); it != s_edges.end(); it++) {
		Edge& edge = *it;
		if (edge.otherNode(s) == e) {
			edge.setWeight(weight);
			break;
		}
	}
}

Graph::Node* Graph::elemToNode(MatrixElement* elem)
{
	for (size_t i = 0; i < nodes.size(); i++) {
		Node& n = nodes[i];
		if (n.elem == elem) return &n;
	}
	return NULL;
}

MatrixElement* Graph::nodeToElem(Node& node)
{
	for (size_t i = 0; i < nodes.size(); i++) {
		if (nodes[i] == node) return nodes[i].elem;
	}
	return NULL;
}

bool Graph::isEmpty() const
{
	return nodes.size() == 0;
}

bool Graph::operator==(const Graph& other)
{
	return edges == other.edges;
}

bool operator<(const Graph::Node &n1, const Graph::Node &n2) { return n1.elem->getStrokes() < n2.elem->getStrokes(); }
bool operator<(const Graph::Edge &e1, const Graph::Edge &e2) { return ((ptrdiff_t)e1.n1.elem + (ptrdiff_t)e1.n2.elem) < ((ptrdiff_t)e2.n1.elem + (ptrdiff_t)e2.n2.elem); }


}
