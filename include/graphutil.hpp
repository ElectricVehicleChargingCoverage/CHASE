#pragma once

#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
// #include <scotch.h>
#include <metis.h>

#include "stringutil.hpp"

using namespace std;
using namespace RoutingKit;

struct Graph {
   private:
    int count;

   public:
    Graph(int _count, string _name = "") : count{_count}, name{_name} {}
    struct Side {
        vector<unsigned> first_out;
        vector<unsigned> head;
        vector<unsigned> weight;
        vector<long long> original_arc;
    };
    string name;
    Side forward, backward;
    Side forward_up, forward_down, backward_up, backward_down;
    vector<float> latitude, longitude;
    vector<int> order;
    unordered_map<int, int> new_id;
    unordered_map<int, int> old_id;
    vector<bool> boundary_nodes;
    int node_count() { return count; };
    void compute_boundary_node(vector<int> partition) {
        boundary_nodes.clear();
        boundary_nodes.resize(count);
        assert(partition.size() == count);
        for (int x = 0; x < count; ++x) {
            for (int arc = forward.first_out[x]; arc < forward.first_out[x + 1]; ++arc) {
                int y = forward.head[arc];
                if (partition[x] != partition[y]) {
                    boundary_nodes[x] = true;
                    boundary_nodes[y] = true;
                }
            }
        }
    }
};

using namespace RoutingKit;
using namespace std;

vector<vector<int>> ccList;

using QueueElement = pair<int, double>;
auto cmp = [](QueueElement& a, QueueElement& b) {
    return a.second > b.second;
};

Graph loadGraph(string path) {
    ifstream file(path);
    string line;
    unordered_map<long long, pair<double, double>> nodes;
    unordered_map<long long, vector<pair<long long, int>>> adj;
    unordered_map<long long, set<long long>> temp_adj_set;
    getline(file, line);  // skip first line
    while (getline(file, line)) {
        vector<string> csv = split(line, ",", false);
        int cost = stod(csv[5]) * 3600.0 * 1000;
        int reverse_cost = stod(csv[6]) * 3600.0 * 1000;
        long long source_id = stoll(csv[1]);
        long long target_id = stoll(csv[2]);
        double from_lng = stod(csv[7]);
        double from_lat = stod(csv[8]);
        double to_lng = stod(csv[9]);
        double to_lat = stod(csv[10]);
        nodes[source_id] = make_pair(from_lat, from_lng);
        nodes[target_id] = make_pair(to_lat, to_lng);
        if (source_id != target_id) {
            if (temp_adj_set[source_id].find(target_id) == temp_adj_set[source_id].end()) {
                adj[source_id].push_back(make_pair(target_id, cost));
                temp_adj_set[source_id].insert(target_id);
            }
            if (stod(csv[6]) < 1000000.0 && temp_adj_set[target_id].find(source_id) == temp_adj_set[target_id].end()) {
                adj[target_id].push_back(make_pair(source_id, reverse_cost));
                temp_adj_set[target_id].insert(source_id);
            }
        }
    }
    vector<long long> n;
    unordered_map<long long, int> new_id;
    for (auto& [from, latlng] : nodes) {
        new_id[from] = n.size();
        n.push_back(from);
    }
    Graph g(nodes.size(), path);
    for (int i = 0; i < n.size(); ++i) {
        long long from = n[i];
        auto latlng = nodes[from];
        g.forward.first_out.push_back(g.forward.head.size());
        g.latitude.push_back(latlng.first);
        g.longitude.push_back(latlng.second);
        for (auto& [to, cost] : adj[from]) {
            g.forward.head.push_back(new_id[to]);
            g.forward.weight.push_back(cost);
        }
    }
    g.forward.first_out.push_back(g.forward.head.size());
    return g;
}

Graph removeSmallCCs(Graph& graph, int minSize) {
    cout << "Graph size before cleanup: ";
    cout << graph.node_count() << " nodes, " << graph.forward.head.size() << " edges" << endl;
    vector<bool> deleted(graph.node_count(), false);
    for (auto cc : ccList)
        if (cc.size() < minSize)
            for (int node : cc) deleted[node] = true;

    vector<vector<pair<int, int>>> adj(graph.node_count());
    for (int x = 0; x < graph.node_count(); ++x) {
        if (deleted[x]) continue;
        for (int arc = graph.forward.first_out[x]; arc < graph.forward.first_out[x + 1]; ++arc) {
            int y = graph.forward.head[arc];
            if (deleted[y]) continue;
            adj[x].push_back(make_pair(y, graph.forward.weight[arc]));
        }
    }
    unordered_map<int, int> id_map;
    vector<int> nodes;
    for (int x = 0, count = 0; x < graph.node_count(); ++x) {
        if (deleted[x]) continue;
        id_map[x] = count;
        count++;
        nodes.push_back(x);
    }
    Graph cleaned = Graph(nodes.size(), graph.name);
    for (int x : nodes) {
        cleaned.forward.first_out.push_back(cleaned.forward.head.size());
        cleaned.latitude.push_back(graph.latitude[x]);
        cleaned.longitude.push_back(graph.longitude[x]);
        for (auto pair : adj[x]) {
            cleaned.forward.head.push_back(id_map[pair.first]);
            cleaned.forward.weight.push_back(pair.second);
        }
    }
    cleaned.forward.first_out.push_back(cleaned.forward.head.size());

    cout << "Graph size after cleanup: ";
    cout << cleaned.node_count() << " nodes, " << cleaned.forward.head.size() << " edges" << endl;
    return cleaned;
}

void bfs(Graph& graph, int src, vector<bool>& visited, int comp) {
    visited[src] = true;
    queue<int> q{{{src}}};
    while (!q.empty()) {
        auto v = q.front();
        ccList[comp].push_back(v);
        q.pop();
        for (int arc = graph.forward.first_out[v]; arc < graph.forward.first_out[v + 1]; ++arc) {
            int y = graph.forward.head[arc];
            if (visited[y]) continue;
            visited[y] = true;
            q.push(y);
        }
    }
}

int findCCs(Graph& graph) {
    vector<bool> visited(graph.node_count(), false);
    int scc = 0;
    cout << "start cc search" << endl;
    for (int x = 0; x < graph.node_count(); ++x) {
        if (visited[x]) continue;
        ccList.push_back(vector<int>());
        bfs(graph, x, visited, scc);
        scc++;
    }
    cout << "cc search finished" << endl;
    cout << "cc size: " << ccList.size() << endl;
    vector<int> sizes(ccList.size());
    for (auto scc : ccList) sizes.push_back(scc.size());
    sort(sizes.begin(), sizes.end());
    reverse(sizes.begin(), sizes.end());
    cout << "The sizes of the first 50 connected components are: ";
    for (int i = 0; i < 50; i++) cout << sizes[i] << ",";
    cout << endl;
    cout << sizes[0] << endl;
    return sizes[0];
}

Graph build_ch_complete_graph(string name, Graph& graph, ContractionHierarchy& ch, float best_percentage) {
    struct Info {
        Info(unsigned _y, long long _arc, unsigned _weight) : y{_y}, arc{_arc}, weight{_weight} {}
        unsigned y;
        long long arc;
        unsigned weight;
    };
    vector<vector<Info>> adj(ch.node_count());
    vector<vector<Info>> back_adj(ch.node_count());
    assert(ch.forward.first_out.size() == ch.backward.first_out.size());

    for (int x = 0; x < ch.node_count(); ++x) {
        if (x < ch.node_count() * (1 - best_percentage)) continue;
        for (int arc = ch.forward.first_out[x]; arc < ch.forward.first_out[x + 1]; ++arc) {
            int y = ch.forward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[x].push_back(Info(y, arc, ch.forward.weight[arc]));
            back_adj[y].push_back(Info(x, arc, ch.forward.weight[arc]));
        }
        for (int arc = ch.backward.first_out[x]; arc < ch.backward.first_out[x + 1]; ++arc) {
            int y = ch.backward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[y].push_back(Info(x, -arc, ch.backward.weight[arc]));
            back_adj[x].push_back(Info(y, -arc, ch.backward.weight[arc]));
        }
    }
    unordered_map<int, int> id_map;
    vector<int> nodes;
    for (int x = 0, count = 0; x < graph.node_count(); ++x) {
        if (x < ch.node_count() * (1 - best_percentage)) continue;
        id_map[x] = count;
        count++;
        nodes.push_back(x);
    }
    Graph g = Graph(nodes.size(), name);
    for (int x : nodes) {
        g.new_id[x] = g.forward.first_out.size();
        g.old_id[g.forward.first_out.size()] = x;

        g.forward.first_out.push_back(g.forward.head.size());
        g.backward.first_out.push_back(g.backward.head.size());

        g.latitude.push_back(graph.latitude[ch.order[x]]);
        g.longitude.push_back(graph.longitude[ch.order[x]]);
        g.order.push_back(ch.order[x]);

        for (auto& info : adj[x]) {
            g.forward.head.push_back(id_map[info.y]);
            g.forward.weight.push_back(info.weight);
            g.forward.original_arc.push_back(info.arc);
        }
        for (auto& info : back_adj[x]) {
            g.backward.head.push_back(id_map[info.y]);
            g.backward.weight.push_back(info.weight);
            g.backward.original_arc.push_back(info.arc);
        }
    }
    g.forward.first_out.push_back(g.forward.head.size());
    g.backward.first_out.push_back(g.backward.head.size());
    return g;
}

void build_up_down_cores(Graph& g) {
    for (int x = 0; x < g.node_count(); ++x) {
        g.forward_up.first_out.push_back(g.forward_up.head.size());
        g.forward_down.first_out.push_back(g.forward_down.head.size());

        g.backward_up.first_out.push_back(g.backward_up.head.size());
        g.backward_down.first_out.push_back(g.backward_down.head.size());

        for (int arc = g.forward.first_out[x]; arc < g.forward.first_out[x + 1]; ++arc) {
            int y = g.forward.head[arc];
            if (g.old_id[x] < g.old_id[y]) {
                g.forward_up.head.push_back(y);
                g.forward_up.weight.push_back(g.forward.weight[arc]);
                g.forward_up.original_arc.push_back(g.forward.original_arc[arc]);
            } else {
                g.forward_down.head.push_back(y);
                g.forward_down.weight.push_back(g.forward.weight[arc]);
                g.forward_down.original_arc.push_back(g.forward.original_arc[arc]);
            }
        }
        for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x + 1]; ++arc) {
            int y = g.backward.head[arc];
            if (g.old_id[x] < g.old_id[y]) {
                g.backward_up.head.push_back(y);
                g.backward_up.weight.push_back(g.backward.weight[arc]);
                g.backward_up.original_arc.push_back(g.backward.original_arc[arc]);
            } else {
                g.backward_down.head.push_back(y);
                g.backward_down.weight.push_back(g.backward.weight[arc]);
                g.backward_down.original_arc.push_back(g.backward.original_arc[arc]);
            }
        }
    }
    g.forward_up.first_out.push_back(g.forward_up.head.size());
    g.backward_up.first_out.push_back(g.backward_up.head.size());
    g.forward_down.first_out.push_back(g.forward_down.head.size());
    g.backward_down.first_out.push_back(g.backward_down.head.size());
}

vector<int> partition_graph(Graph& g, int nparts, float best_percentage) {
    int n = g.node_count();

    vector<int> xadj;
    vector<int> adjncy;
    for (int x = 0; x < g.node_count(); ++x) {
        xadj.push_back(adjncy.size());
        vector<bool> is_set(g.node_count(), false);
        // shortcut edges
        for (int arc = g.forward.first_out[x]; arc < g.forward.first_out[x + 1]; ++arc) {
            int y = g.forward.head[arc];
            if (is_set[y] || y == x) continue;
            adjncy.push_back(y);
            is_set[y] = true;
        }
        for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x + 1]; ++arc) {
            int y = g.backward.head[arc];
            if (is_set[y] || y == x) continue;
            adjncy.push_back(y);
            is_set[y] = true;
        }
    }
    xadj.push_back(adjncy.size());

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    int ncon = 1, edgecut = 0;
    vector<int> part(n);
    METIS_PartGraphRecursive(&n, &ncon, &(xadj[0]), &(adjncy[0]), NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, &(part[0]));

    // SCOTCH_Graph grafdat;
    // SCOTCH_Strat stradat;
    // SCOTCH_Num baseval;
    // SCOTCH_Num vertnbr;
    // int o;

    // SCOTCH_graphInit(&grafdat);
    // baseval = 0;
    // vertnbr = n;
    // o = 1;
    // if (SCOTCH_graphBuild(&grafdat, baseval, vertnbr, &(xadj[0]), &(xadj[1]), NULL, NULL, adjncy.size(), &(adjncy[0]), NULL) == 0) {
    //     SCOTCH_stratInit(&stradat);
    //     if (SCOTCH_graphCheck(&grafdat) == 0) {
    //         o = SCOTCH_graphPart(&grafdat, nparts, &stradat, &(part[0]));
    //     }
    //     SCOTCH_stratExit(&stradat);
    // }
    // SCOTCH_graphExit(&grafdat);


    vector<int> partition(n);
    for (int i = 0; i < n; ++i) partition[i] = part[i];
    return partition;
}

vector<int> read_partition(string path) {
    vector<int> partition;
    ifstream file(path);
    string line;
    while (getline(file, line)) {
        partition.push_back(stoi(line));
    }
    return partition;
}

void export_partition(vector<int> partition, string location) {
    ofstream out(location);
    for (int cell_id : partition) {
        out << cell_id << endl;
    }
    out.close();
}