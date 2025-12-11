#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "utility/utility.h"
#include "configuration/types.h"
#include "configuration/config.h"

class Graph
{
private:
	std::string graph_id;	// graph id
	ui 			n;			// number of vertices
	ept 		m;			// number of edges
	ept 		*pstart;	// pstart[u]: starting position of vertex u in edges array
	ui 		*edges;		// edges array
	LabelID 	*vlabels;	// vertex labels array
// #ifdef ENABLE_EDGE_LABEL
	LabelID 	*elabels;	// edge labels array
// #endif

	ui			_max_degree; // maximum degree of the graph
	std::unordered_map<LabelID, ui> vertex_labels_frequency; // frequency of each vertex label
public:
	Graph()
	{
		graph_id = "";
		n = 0;
		m = 0;
		pstart = nullptr;
		edges = nullptr;
		vlabels = nullptr;
#ifdef ENABLE_EDGE_LABEL
		elabels = nullptr;
#endif
		_max_degree = 0;
	}

	~Graph()
	{
		if(pstart != nullptr) { delete[] pstart; pstart = nullptr; }
		if(edges != nullptr) { delete[] edges; edges = nullptr; }
		if(vlabels != nullptr) { delete[] vlabels; vlabels = nullptr; }
#ifdef ENABLE_EDGE_LABEL
		if(elabels != nullptr) { delete[] elabels; elabels = nullptr; }
#endif
	}

	void build_graph(const std::string &id,
					 std::vector<std::pair<ui, LabelID> > &vertices, 
					 std::vector<std::pair<std::pair<ui, ui>, LabelID> > &edges_list);
	
	void print_graph() const;

	ui getVerticesCount() const {
		return n;
	}

	ept getEdgesCount() const {
		return m;
	}

	ui getVertexDegree(ui u) const {
		assert(u >= 0 && u < n);
		return pstart[u + 1] - pstart[u];
	}

	ui getMaxDegree() const {
		return _max_degree;
	}

	LabelID getVertexLabel(ui u) const {
		assert(u >= 0 && u < n);
		return vlabels[u];
	}

	ui getVertexLabelsFrequency(LabelID label) {
		return vertex_labels_frequency.find(label) == vertex_labels_frequency.end() ? 0 : vertex_labels_frequency.at(label);
	}

    const ui * getVertexNeighbors(const ui id, ui& count) const {
        count = pstart[id + 1] - pstart[id];
        return edges + pstart[id];
    }

	// Check if there is an edge u -> v
	bool hasEdge(ui u, ui v) const {
		assert(u < n && v < n);
		ept l = pstart[u], r = pstart[u + 1];
		// binary search in sorted adjacency list
		return std::binary_search(edges + l, edges + r, v);
	}
};

#endif