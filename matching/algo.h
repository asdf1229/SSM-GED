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

bool calVerticesFilter(const Graph *data_graph, const Graph *query_graph,
                       vector<vector<ui> > &candidates) {
    ui qn = query_graph->getVerticesCount();
    ui dn = data_graph->getVerticesCount();

    candidates.clear();
    candidates.resize(qn);

    for (ui i = 0; i < qn; ++i) {
        int label_i = query_graph->getVertexLabel(i);

        for (ui j = 0; j < dn; ++j) {
            int label_j = data_graph->getVertexLabel(j);
            if (label_i == label_j) candidates[i].push_back(j);
        }

        if (candidates[i].empty())  return false;
    }

    return true;
}

void generateQueryOrder(const Graph *data_graph, const Graph *query_graph,
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

void Approximate_Matching(const Graph *data_graph, const Graph *query_graph,
                          vector<vector<ui> > &candidates, 
                          vector<ui> &order,
                          vector<vector<pair<ui, ui> > > &M_ANS,
                          ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui dn = data_graph->getVerticesCount();

    vector<ui> embedding(qn, (ui)(-1));  // embedding[depth] → data vertex
    vector<ui> vis(dn, 0);

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

// ============================================================
// Helper: Compute frontier in the query graph
// ============================================================
vector<ui> ComputeFrontierQuery(const Graph* query_graph,
                                const vector<ui>& Mq,
                                const vector<int>& mapped_q)
{
    ui qn = query_graph->getVerticesCount();
    vector<char> mark(qn, 0);
    for (ui u : Mq) mark[u] = 1;

    vector<ui> Nq;

    for (ui u : Mq) {
        ui count;
        const ui* nbrs = query_graph->getVertexNeighbors(u, count);
        for (ui i = 0; i < count; ++i) {
            ui w = nbrs[i];
            if (!mark[w] && mapped_q[w] == -1) {
                Nq.push_back(w);
            }
        }
    }
    return Nq;
}


// ============================================================
// Helper: Compute frontier in the data graph
// ============================================================
vector<ui> ComputeFrontierData(const Graph* data_graph,
                               const vector<ui>& Mq,
                               const vector<ui>& Mg,
                               const vector<int>& used_v)
{
    ui dn = data_graph->getVerticesCount();
    vector<char> mark(dn, 0);
    for (ui v : Mg) mark[v] = 1;

    vector<ui> Ng;

    for (ui mapped_v : Mg) {
        ui count;
        const ui* nbrs = data_graph->getVertexNeighbors(mapped_v, count);
        for (ui i = 0; i < count; ++i) {
            ui w = nbrs[i];
            if (!mark[w] && !used_v[w]) {
                Ng.push_back(w);
            }
        }
    }
    return Ng;
}


// ============================================================
// Helper: Compute current missing edges based on Mq→Mg mapping
// ============================================================
ui CalPartialMissing(const Graph* query_graph,
                     const Graph* data_graph,
                     const vector<ui>& Mq,
                     const vector<ui>& Mg)
{
    ui missing = 0;
    for (ui i = 0; i < Mq.size(); ++i) {
        ui u = Mq[i];
        ui vu = Mg[i];
        for (ui j = i + 1; j < Mq.size(); ++j) {
            ui w = Mq[j];
            ui vw = Mg[j];
            if (query_graph->hasEdge(u, w)) {
                if (!data_graph->hasEdge(vu, vw))
                    missing++;
            }
        }
    }
    return missing;
}

// ============================================================
// Main DFS for Approximate Matching via Bidirectional Adjacency Expansion
// ============================================================
void DFS_Approximate(const Graph* data_graph,
                     const Graph* query_graph,
                     vector<vector<ui>>& candidates,
                     vector<vector<char>>& X,             // exclusion set
                     vector<int>& mapped_q,
                     vector<int>& used_v,
                     vector<ui>& Mq,
                     vector<ui>& Mg,
                     vector<vector<pair<ui,ui>>>& M_ANS,
                     ui threshold)
{
    ui qn = query_graph->getVerticesCount();

    // ---------------------------
    // Full mapping reached
    // ---------------------------
    if (Mq.size() == qn) {
        ui missing = CalPartialMissing(query_graph, data_graph, Mq, Mg);
        if (missing <= threshold) {
            vector<pair<ui,ui>> res;
            for (ui i = 0; i < Mq.size(); ++i)
                res.emplace_back(Mq[i], Mg[i]);
            M_ANS.push_back(res);
        }
        return;
    }

    // ---------------------------
    // Compute frontiers
    // ---------------------------
    vector<ui> Nq = ComputeFrontierQuery(query_graph, Mq, mapped_q);
    vector<ui> Ng = ComputeFrontierData(data_graph, Mq, Mg, used_v);

    // 初始时没有 frontier，则选任意未匹配的点作为起点
    if (Nq.empty()) {
        for (ui u = 0; u < qn; ++u)
            if (mapped_q[u] == -1) { Nq.push_back(u); break; }
    }

    // ---------------------------
    // Construct expandable mapping set Mcand
    // ---------------------------
    vector<pair<ui,ui>> Mcand;

    for (ui u : Nq) {
        for (ui v : candidates[u]) {
            if (mapped_q[u] == -1 && !used_v[v] && !X[u][v]) {

                // 如果 Ng 非空，则 v 必须在 Ng 中
                if (!Ng.empty() && find(Ng.begin(), Ng.end(), v) == Ng.end())
                    continue;

                Mcand.emplace_back(u, v);
            }
        }
    }

    // ---------------------------
    // Enumerate branches
    // ---------------------------
    for (auto &p : Mcand) {
        ui u = p.first;
        ui v = p.second;

        // extend mapping
        mapped_q[u] = v;
        used_v[v] = 1;
        Mq.push_back(u);
        Mg.push_back(v);

        ui missing = CalPartialMissing(query_graph, data_graph, Mq, Mg);
        if (missing <= threshold) {
            DFS_Approximate(data_graph, query_graph, candidates, X,
                            mapped_q, used_v, Mq, Mg, M_ANS, threshold);
        }

        // backtrack
        mapped_q[u] = -1;
        used_v[v] = 0;
        Mq.pop_back();
        Mg.pop_back();

        // exclusion set update
        X[u][v] = 1;
    }
}

// ============================================================
// Top-level function: Approximate_Matching_v2
// ============================================================
void Approximate_Matching_v2(const Graph *data_graph,
                             const Graph *query_graph,
                             vector<vector<ui> > &candidates, 
                             vector<ui> &order, // unused in v2
                             vector<vector<pair<ui, ui> > > &M_ANS,
                             ui threshold)
{
    ui qn = query_graph->getVerticesCount();
    ui dn = data_graph->getVerticesCount();

    M_ANS.clear();
    if (qn == 0) return;

    vector<int> mapped_q(qn, -1);
    vector<int> used_v(dn, 0);
    vector<vector<char>> X(qn, vector<char>(dn, 0));  // exclusion set

    vector<ui> Mq, Mg;  // matched sets (order varies dynamically)

    DFS_Approximate(data_graph, query_graph, candidates, X,
                    mapped_q, used_v, Mq, Mg, M_ANS, threshold);
}