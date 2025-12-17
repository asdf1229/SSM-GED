void generateQueryOrder(const Graph *data_graph, const Graph *query_graph,
                        const vector<vector<ui>> &candidates,
                        vector<ui> &order) {
    ui qn = query_graph->getVerticesCount();

    order.clear();
    order.reserve(qn);

    // 预计算所有点的度数
    vector<ui> degrees(qn);
    for (ui i = 0; i < qn; ++i) {
        degrees[i] = query_graph->getVertexDegree(i);
    }

    vector<bool> selected(qn, false);
    ui selected_count = 0;

    // -----------------------------
    // 1. 第一个点：度数最大的点
    // -----------------------------
    ui first = 0;
    for (ui i = 1; i < qn; ++i) {
        if (degrees[i] > degrees[first]) {
            first = i;
        }
    }
    order.push_back(first);
    selected[first] = true;
    selected_count++;

    // -----------------------------
    // 2. 后续点：与已选集合连通的点中度数最大的
    // -----------------------------
    while (selected_count < qn) {
        ui best = UINT32_MAX;
        ui best_deg = 0;
        bool found_connected = false;

        // 在所有未选择的点里找与当前已选集合连通的点
        for (ui v = 0; v < qn; ++v) {
            if (selected[v]) continue;

            // 检查是否与已选集合连通
            bool connected = false;
            for (ui u : order) {
                if (query_graph->hasEdge(u, v)) {
                    connected = true;
                    break;
                }
            }

            if (connected) {
                if (!found_connected || degrees[v] > best_deg) {
                    best = v;
                    best_deg = degrees[v];
                    found_connected = true;
                }
            }
        }

        if (!found_connected) {
            // 图可能不连通 → 从未选择中选度数最大的点
            for (ui v = 0; v < qn; ++v) {
                if (!selected[v] && degrees[v] > best_deg) {
                    best = v;
                    best_deg = degrees[v];
                }
            }
        }

        // 选入 order
        order.push_back(best);
        selected[best] = true;
        selected_count++;
    }
}

// vector<ui> ComputeFrontierQuery(const Graph* query_graph,
//                                 const vector<ui>& Mq,
//                                 const vector<int>& mapped_q)
// {
//     vector<ui> Nq;

//     ui qn = query_graph->getVerticesCount();

//     if(Mq.size() == 0) {
//         Nq = vector<ui>{0};
//         return Nq;
//     }

//     vector<char> mark(qn, 0);
//     for (ui u : Mq) {
//         assert(u < qn);
//         mark[u] = 1;
//     }

//     for (ui u : Mq) {
//         ui count;
//         const ui* neighbors = query_graph->getVertexNeighbors(u, count);
//         for (ui i = 0; i < count; ++i) {
//             ui w = neighbors[i];
//             if (!mark[w] && mapped_q[w] == -1) {
//                 Nq.push_back(w);
//             }
//         }
//     }
//     return Nq;
// }

// vector<ui> ComputeFrontierData(const Graph *query_graph, const Graph *data_graph,
//                                const vector<ui>& Mq,
//                                const vector<ui>& Mg,
//                                const vector<int>& visited_v,
//                                vector<vector<ui> > &candidates)
// {
//     vector<ui> Ng;

//     ui gn = data_graph->getVerticesCount();
//     ui qn = query_graph->getVerticesCount();

//     vector<char> mark_g(gn, 0);
//     for (ui v : Mg) {
//         assert(v < gn);
//         mark_g[v] = 1;
//     }

//     vector<char> mark_q(qn, 0);
//     for (ui u : Mq) {
//         assert(u < qn);
//         mark_q[u] = 1;
//     }

//     vector<char> in_Ng(gn, 0);

//     vector<char> is_neighbor_of_v(gn, 0);

//     for (size_t idx = 0; idx < Mq.size(); ++idx) {
//         ui u = Mq[idx];
//         ui v = Mg[idx];

//         ui deg_v;
//         const ui* neigh_v = data_graph->getVertexNeighbors(v, deg_v);
//         for (ui i = 0; i < deg_v; ++i) {
//             ui w = neigh_v[i];
//             is_neighbor_of_v[w] = 1;
//         }

//         ui deg_u;
//         const ui* neigh_u = query_graph->getVertexNeighbors(u, deg_u);
//         for (ui i = 0; i < deg_u; ++i) {
//             ui u2 = neigh_u[i];
//             if (mark_q[u2]) continue;

//             const auto& cand_u2 = candidates[u2];
//             for (ui w : cand_u2) {
//                 assert(w < gn);
//                 if (!mark_g[w] && !visited_v[w] && is_neighbor_of_v[w] && !in_Ng[w]) {
//                     Ng.push_back(w);
//                     in_Ng[w] = 1;
//                 }
//             }
//         }

//         for (ui i = 0; i < deg_v; ++i) {
//             is_neighbor_of_v[neigh_v[i]] = 0;
//         }
//     }

//     return Ng;
// }


void computeFrontier(const Graph *query_graph, const Graph *data_graph,
                     const vector<int> &mapped_q, const vector<int> &mapped_g,
                     const vector<pair<ui, ui>> &part_M, const vector<vector<ui> > &candidates,
                     vector<ui> &Nq, vector<ui> &Ng)
{
    ui qn = query_graph->getVerticesCount();
    ui gn = data_graph->getVerticesCount();
    ui mn = part_M.size();

    // part_M = {}
    if(mn == 0) {
        // select 0 as starting vertex
        ui u = 0;
        Nq = vector<ui>{u};
        Ng = candidates[u];
        return;
    }



    // calculate Nq
    // find all vertices in query graph adjacent to Mq but not in Mq
    unordered_set<ui> Nq_set;
    for (auto p : part_M) {
        ui u = p.first;
        ui count;
        const ui* neighbors = query_graph->getVertexNeighbors(u, count);
        for (ui i = 0; i < count; ++i) {
            ui w = neighbors[i];
            if (mapped_q[w] == -1) {
                Nq_set.insert(w);
            }
        }
    }
    Nq.assign(Nq_set.begin(), Nq_set.end());

    // calculate Ng
    // find all vertices in data graph adjacent to Mg (not in Mg) and are candidates of some vertex in Nq
    unordered_set<ui> Ng_set;
    
    for(auto p : part_M) {
        ui v = p.second;

        ui count;
        const ui* neighbors = data_graph->getVertexNeighbors(v, count);
        for (ui i = 0; i < count; ++i) {
            ui w = neighbors[i];
            if (mapped_g[w] == -1) {
                Ng_set.insert(w);
            }
        }
    }

    unordered_set<ui> Ng_filtered;
    for(ui u : Nq) {
        for(ui v : candidates[u]) {
            if(Ng_set.count(v)) {
                Ng_filtered.insert(v);
            }
        }
    }

    Ng.assign(Ng_filtered.begin(), Ng_filtered.end());
}