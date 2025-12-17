#include "utility/utility.h"
// #include "graph/graph.h"

using namespace std;

LabelID label2int(const std::string str, std::map<std::string, LabelID> &M) {
	if(M.find(str) == M.end()) M[str] = M.size();
	return M[str];
}

void load_graph(const std::string &input_graph, Graph *graph,
                std::map<std::string, LabelID> &vM, std::map<std::string, LabelID> &eM) {

    std::ifstream fin(input_graph);
    std::string line;

    while (std::getline(fin, line)) if (!line.empty() && line[0] == 't') break;

    if(!fin) {
        printf("!!! Cannot open graph file %s !!!\n", input_graph.c_str());
        assert(false);
        return;
    }

    std::istringstream head(line);
    char tchar; 
    std::string sharp, id;
    head >> tchar >> sharp >> id;
    assert(tchar == 't');

    std::vector<std::pair<ui, LabelID> > vertices;
    std::vector<std::pair<std::pair<ui, ui>, LabelID> > undirected_edges;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        char type = line[0];
        if (type == 't') break;

        std::istringstream iss(line);

        if(type == 'v') {
            char c;
            ui vid;
            std::string vlab;
            iss >> c >> vid >> vlab;
            vertices.emplace_back(vid, label2int(vlab, vM));
        }
        else if (type == 'e') {
            char c;
            ui u, v;
            std::string elab;
            iss >> c >> u >> v >> elab;
            if(u == v) continue;
            LabelID L = label2int(elab, eM);
            undirected_edges.emplace_back(std::make_pair(std::min(u, v), std::max(u, v)), L);
        }
    }

    std::sort(vertices.begin(), vertices.end());
    std::sort(undirected_edges.begin(), undirected_edges.end());

    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
    undirected_edges.erase(unique(undirected_edges.begin(), undirected_edges.end()), undirected_edges.end());

    std::vector<std::pair<std::pair<ui, ui>, LabelID> > edges;

    for(auto &e: undirected_edges) {
        edges.emplace_back(std::make_pair(e.first.first, e.first.second), e.second);
        edges.emplace_back(std::make_pair(e.first.second, e.first.first), e.second);
    }

    std::sort(edges.begin(), edges.end());

    graph->build_graph(id, vertices, edges);

// #ifndef NDEBUG
//     graph.print_graph();
// #endif

	return;
}

bool calVerticesFilter(const Graph *query_graph, const Graph *data_graph,
                       vector<vector<ui> > &candidates) {
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    candidates.clear();
    candidates.resize(qn);

    for (ui i = 0; i < qn; ++i) {
        int label_i = query_graph->getVertexLabel(i);

        for (ui j = 0; j < gn; ++j) {
            int label_j = data_graph->getVertexLabel(j);
            if (label_i == label_j) candidates[i].push_back(j);
        }

        if (candidates[i].empty())  return false;
    }

    return true;
}

void generateQueryOrder(const Graph *query_graph, const Graph *data_graph,
                        const vector<vector<ui> > &candidates,
                        vector<ui> &order) {
    ui qn = query_graph->getVerticesCount();

    order.clear();
    order.resize(qn, 0);

    vector<pair<ui, ui> > degree_order; // (degree, vertex id)
    for (ui i = 0; i < qn; ++i) {
        degree_order.emplace_back(query_graph->getVertexDegree(i), i);
    }

    // ascending order
    sort(degree_order.begin(), degree_order.end());

    for (ui i = 0; i < qn; ++i) {
        order[i] = degree_order[i].second;
    }

    return;
}

void Approximate_Matching(const Graph *query_graph, const Graph *data_graph,
                          vector<vector<ui> > &candidates, 
                          vector<ui> &order,
                          vector<vector<pair<ui, ui> > > &M_ANS,
                          ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    vector<ui> embedding(qn, (ui)(-1));  // embedding[depth] → data vertex
    vector<ui> vis(gn, 0);

    vector<ui> idx(qn, 0);
    vector<ui> idx_count(qn, 0);

    vector<ui> missing_edges(qn, 0);

    M_ANS.clear();
    if (qn == 0) return;

    ui cur_depth = 0;
    ui max_depth = qn;

    {   // initial state
        ui u0 = order[0];
        idx_count[0] = (ui)candidates[u0].size();
        missing_edges[0] = 0;
    }

    // -------------------------
    // Helper lambdas
    // -------------------------
    auto compute_missing_edges_final = [&](const vector<ui>& emb) {
        ui missing = 0;
        for (ui i = 0; i < qn; ++i) {
            ui uq = order[i];
            ui vq = emb[i];
            for (ui j = i + 1; j < qn; ++j) {
                ui uw = order[j];
                ui vw = emb[j];
                if (query_graph->hasEdge(uq, uw)) {
                    if (!data_graph->hasEdge(vq, vw)) {
                        missing++;
                    }
                }
            }
        }
        return missing;
    };

    auto check_connectivity_final = [&](const vector<ui>& emb) {
        // 用 query 的边，构造映射后的 “有效边”
        vector<vector<ui>> adj(qn);
        for (ui i = 0; i < qn; ++i) {
            ui uq = order[i];
            ui vq = emb[i];
            for (ui j = 0; j < qn; ++j) {
                ui uw = order[j];
                ui vw = emb[j];
                if (query_graph->hasEdge(uq, uw)) {
                    if (data_graph->hasEdge(vq, vw)) {
                        adj[i].push_back(j);
                    }
                }
            }
        }

        // BFS/DFS 判断是否连通
        vector<int> visited(qn, 0);
        vector<ui> stack{0};
        visited[0] = 1;

        while (!stack.empty()) {
            ui x = stack.back();
            stack.pop_back();
            for (ui y : adj[x]) {
                if (!visited[y]) {
                    visited[y] = 1;
                    stack.push_back(y);
                }
            }
        }

        for (ui i = 0; i < qn; ++i)
            if (!visited[i]) return false;

        return true;
    };

    // -------------------------
    // Backtracking
    // -------------------------
    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {

            ui u = order[cur_depth];
            ui v = candidates[u][idx[cur_depth]];
            idx[cur_depth]++;

            if (vis[v]) continue;

            // compute new missing edges incrementally (based only on query edges)
            ui new_missing = (cur_depth == 0) ? 0 : missing_edges[cur_depth - 1];

            for (ui i = 0; i < cur_depth; ++i) {
                ui u_prev = order[i];
                ui v_prev = embedding[i];

                if (query_graph->hasEdge(u, u_prev)) {
                    if (!data_graph->hasEdge(v, v_prev)) {
                        new_missing++;
                    }
                }
            }

            if (new_missing > threshold) continue;

            embedding[cur_depth] = v;
            vis[v] = 1;
            missing_edges[cur_depth] = new_missing;

            // -------------------------
            // Reached full mapping
            // -------------------------
            if (cur_depth == max_depth - 1) {

                // 1. 重新精确计算最终缺失边（只用 query 的边）
                ui final_missing = compute_missing_edges_final(embedding);

                // 2. 再检查是否连通（只用 query 的边映射后形成的有效边）
                if (final_missing <= threshold &&
                    check_connectivity_final(embedding))
                {
                    vector<pair<ui, ui>> mapping;
                    for (ui i = 0; i < qn; ++i) {
                        mapping.push_back(make_pair(order[i], embedding[i]));
                    }
                    sort(mapping.begin(), mapping.end());
                    M_ANS.push_back(mapping);
                }

                // 回溯
                vis[v] = 0;
                embedding[cur_depth] = (ui)(-1);
                missing_edges[cur_depth] = 0;
            }
            else {
                // deeper search
                cur_depth++;
                ui u_next = order[cur_depth];
                idx[cur_depth] = 0;
                idx_count[cur_depth] = (ui)candidates[u_next].size();
                missing_edges[cur_depth] = missing_edges[cur_depth - 1];
            }
        }

        if (cur_depth == 0) break;

        // backtrack
        cur_depth--;
        ui v_back = embedding[cur_depth];
        vis[v_back] = 0;
        embedding[cur_depth] = (ui)(-1);
    }
}

void computeFrontierAndMcand(const Graph *query_graph, const Graph *data_graph,
                     const vector<int> &mapped_q, const vector<int> &mapped_g,
                     const vector<pair<ui, ui>> &part_M, const vector<vector<ui> > &candidates,
                     const unordered_set<pair<ui, ui>, PairHash> X, vector<pair<ui,ui> > &Mcand)
{
    Mcand.clear();

    // part_M = {}
    if(part_M.empty()) {
        // select 0 as starting vertex
        ui u = 0;
        assert(mapped_q[u] == -1);
        for(ui v : candidates[u]) {
            assert(mapped_g[v] == -1);
            if(X.count({u, v})) continue;
            Mcand.emplace_back(u, v);
        }
        return;
    }

    unordered_set<uint64_t> seen;
    for(auto p : part_M) {
        ui u = p.first;
        ui v = p.second;

        ui q_count;
        const ui* q_neighbors = query_graph->getVertexNeighbors(u, q_count);
        for(ui i = 0; i < q_count; ++i) {
            ui u2 = q_neighbors[i];
            if(mapped_q[u2] != -1) continue;

            ui g_count;
            const ui* g_neighbors = data_graph->getVertexNeighbors(v, g_count);
            for(ui j = 0; j < g_count; ++j) {
                ui v2 = g_neighbors[j];

                if(mapped_g[v2] != -1) continue;

                if(find(candidates[u2].begin(), candidates[u2].end(), v2) == candidates[u2].end())
                    continue;
                
                if(X.count({u2, v2})) continue;

                uint64_t key = (uint64_t(u2) << 32) | v2;
                if(seen.insert(key).second) {
                    Mcand.emplace_back(u2, v2);
                }
            }
        }

    }
}


ui calPartialMissing(const Graph *query_graph, const Graph *data_graph,
                     const vector<pair<ui, ui>> &part_M)
{
    ui missing = 0;
    ui mn = part_M.size();
    for (ui i = 0; i < mn; ++i) {
        ui u1 = part_M[i].first;
        ui v1 = part_M[i].second;
        for (ui j = i + 1; j < mn; ++j) {
            ui u2 = part_M[j].first;
            ui v2 = part_M[j].second;
            if (query_graph->hasEdge(u1, u2)) {
                if (!data_graph->hasEdge(v1, v2))
                    missing++;
            }
        }
    }
    return missing;
}

void DFS_Approximate(const Graph *query_graph, const Graph *data_graph,
                     const vector<vector<ui>> &candidates, unordered_set<pair<ui, ui>, PairHash> X,
                     vector<int> &mapped_q, vector<int> &mapped_g,
                     vector<pair<ui, ui>> &part_M, vector<vector<pair<ui,ui>>> &M_ANS,
                     ui threshold)
{
#ifndef NDEBUG
    auto msize = part_M.size();
    // printf("Current mapping:\n");
    for(ui i = 0; i < msize; ++i) {
        ui u = part_M[i].first;
        ui v = part_M[i].second;
        assert(u < query_graph->getVerticesCount());
        assert(v < data_graph->getVerticesCount());
        assert(mapped_q[u] == static_cast<int>(v));
        assert(mapped_g[v] == static_cast<int>(u));
        // printf("%d %d\n", u, v);
    }
    // printf("----\n");
#endif

    ui qn = query_graph->getVerticesCount();

    if (part_M.size() == qn) {
        ui missing = calPartialMissing(query_graph, data_graph, part_M);
        assert(missing <= threshold);
        M_ANS.push_back(part_M);
        return;
    }

    vector<pair<ui,ui> > Mcand;
    computeFrontierAndMcand(query_graph, data_graph, mapped_q, mapped_g, part_M, candidates, X, Mcand);

    if(Mcand.empty()) {
        return;
    }

    for (auto p : Mcand) {
        ui u = p.first;
        ui v = p.second;

        assert(X.count({u, v}) == 0);
        assert(mapped_q[u] == -1);
        assert(mapped_g[v] == -1);
        mapped_q[u] = v;
        mapped_g[v] = u;
        part_M.emplace_back(u, v);

        ui missing = calPartialMissing(query_graph, data_graph, part_M);
        if (missing <= threshold) {
            DFS_Approximate(query_graph, data_graph, candidates, X,
                            mapped_q, mapped_g, part_M, M_ANS, threshold);
        }

        mapped_q[u] = -1;
        mapped_g[v] = -1;
        part_M.pop_back();

        X.insert({u, v});
    }
}

// ============================================================
// Top-level function: Approximate_Matching_v2
// ============================================================
void Approximate_Matching_v2(const Graph *query_graph, const Graph *data_graph,
                             vector<vector<ui> > &candidates, 
                             vector<ui> &order, // unused in v2
                             vector<vector<pair<ui, ui> > > &M_ANS,
                             ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();

    M_ANS.clear();
    assert(qn && gn);

    vector<int> mapped_q(qn, -1);
    vector<int> mapped_g(gn, -1);
    unordered_set<pair<ui, ui>, PairHash> X;
    vector<pair<ui, ui>> part_M; // (Mq, Mg)

    DFS_Approximate(query_graph, data_graph, candidates, X,
                    mapped_q, mapped_g, part_M, M_ANS, threshold);
}