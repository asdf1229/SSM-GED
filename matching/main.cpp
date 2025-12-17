#include "utility/utility.h"
#include "graph/graph.h"
#include "utility/popl.hpp"
#include "matching/algo.h"

using namespace popl;
using namespace std;

// Usage: [0]exe [1]input_graph [2]k
int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** SSM-GED (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** SSM-GED (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

    int threshold = 0;

    OptionParser op("Allowed options");
    auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
    auto data_option = op.add<Value<string>>("d", "data graph", "\'data graph file name\'");
	auto query_option = op.add<Value<string>>("q", "query graph", "\'query graph file name\'");
    auto threshold_option = op.add<Value<int>>("t", "threshold", "\'threshold", 0, &threshold);
    op.parse(argc, argv);

    if(help_option->is_set() || argc == 1) cout << op << endl;
	if(!data_option->is_set()||!query_option->is_set()) {
		printf("!!! Data graph file name or query graph file name is not provided! Exit !!!\n");
		return 0;
	}
    if(!threshold_option->is_set()) {
        printf("!!! Threshold is not provided! Exit !!!\n");
        return 0;
    }

    string data_graph_file = data_option->value();
	string query_graph_file = query_option->value();

    // query graph q, data graph G, 
    Graph *query_graph = new Graph();
    Graph *data_graph = new Graph();

    map<string, int> vM, eM;

    load_graph(query_graph_file, query_graph, vM, eM);
    load_graph(data_graph_file, data_graph, vM, eM);

    // // print graphs
    // query_graph->print_graph();
    // data_graph->print_graph();

    // candidate filtering
    vector<vector<ui> > candidates;
    bool res = calVerticesFilter(query_graph, data_graph, candidates);

    // // print candidates
    // printf("Number of query vertices: %lu %u\n", candidates.size(), query_graph->getVerticesCount());
    // for(ui i = 0; i < candidates.size(); i++) {
    //     printf("Query vertex %d: ", i);
    //     for(ui j = 0; j < candidates[i].size(); j++) {
    //         printf("%d ", candidates[i][j]);
    //     }
    //     printf("\n");
    // }

    if(res == false) {
        printf("!!! No candidates for some query vertex. No possible mapping. Exit !!!\n");
        assert(false);
        return 1;
    }

    // generate order
    vector<ui> order;
    generateQueryOrder(query_graph, data_graph, candidates, order);

    // print order
    printf("Query vertex order: ");
    for(ui i = 0; i < order.size(); i++) {
        printf("%d ", order[i]);
    }
    printf("\n");

    vector<vector<pair<ui, ui> > > M_ANS;
    M_ANS.clear();

#ifdef APPROXIMATE_MATCHING_V2
    Approximate_Matching_v2(query_graph, data_graph, candidates, order, M_ANS, threshold);
#else
    Approximate_Matching(query_graph, data_graph, candidates, order, M_ANS, threshold);
#endif
    ui count = 0;
    printf("count: %lu\n", M_ANS.size());
    for(auto &m: M_ANS) {
        printf("Mapping %d:\n", ++count);
        for(auto v: m) {
            printf("%d %d\n", v.first, v.second);
        }
        printf("\n");
    }

    return 0;
}