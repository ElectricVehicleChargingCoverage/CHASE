#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <scotch.h>

#include <arcflags.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
// #include <metis.h>
// #include "kaHIP_interface.h"

using namespace RoutingKit;
using namespace std;

vector<vector<int>> ccList;

using QueueElement = pair<int, double>;
auto cmp = [](QueueElement& a, QueueElement& b) {
    return a.second > b.second;
};

SimpleOSMCarRoutingGraph removeSmallCCs(SimpleOSMCarRoutingGraph& graph, int minSize) {
    cout << "Graph size before cleanup: ";
    cout << graph.node_count() << " nodes, " << graph.arc_count() << " edges" << endl;
    vector<bool> deleted(graph.node_count(), false);
    for (auto cc : ccList)
        if (cc.size() < minSize)
            for (int node : cc) deleted[node] = true;

    vector<vector<pair<int, int>>> adj(graph.node_count());
    for (int x = 0; x < graph.node_count(); ++x) {
        if (deleted[x]) continue;
        for (int arc = graph.first_out[x]; arc < graph.first_out[x + 1]; ++arc) {
            int y = graph.head[arc];
            if (deleted[y]) continue;
            adj[x].push_back(make_pair(y, graph.travel_time[arc]));
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
    SimpleOSMCarRoutingGraph cleaned = SimpleOSMCarRoutingGraph();
    for (int x : nodes) {
        cleaned.first_out.push_back(cleaned.head.size());
        cleaned.latitude.push_back(graph.latitude[x]);
        cleaned.longitude.push_back(graph.longitude[x]);
        for (auto pair : adj[x]) {
            cleaned.head.push_back(id_map[pair.first]);
            cleaned.travel_time.push_back(pair.second);
        }
    }
    cleaned.first_out.push_back(cleaned.head.size());

    cout << "Graph size after cleanup: ";
    cout << cleaned.node_count() << " nodes, " << cleaned.arc_count() << " edges" << endl;
    return cleaned;
}

void bfs(SimpleOSMCarRoutingGraph& graph, int src, vector<bool>& visited, int comp) {
    visited[src] = true;
    queue<int> q{{{src}}};
    while (!q.empty()) {
        auto v = q.front();
        ccList[comp].push_back(v);
        q.pop();
        for (int arc = graph.first_out[v]; arc < graph.first_out[v + 1]; ++arc) {
            int y = graph.head[arc];
            if (visited[y]) continue;
            visited[y] = true;
            q.push(y);
        }
    }
}

int findCCs(SimpleOSMCarRoutingGraph& graph) {
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

Graph build_ch_complete_graph(string name, SimpleOSMCarRoutingGraph& graph, ContractionHierarchy& ch, float best_percentage) {
    vector<vector<pair<int, int>>> adj(ch.node_count());
    vector<vector<pair<int, int>>> back_adj(ch.node_count());
    assert(ch.forward.first_out.size() == ch.backward.first_out.size());

    for (int x = 0; x < ch.node_count(); ++x) {
        if (x < ch.node_count() * (1 - best_percentage)) continue;
        for (int arc = ch.forward.first_out[x]; arc < ch.forward.first_out[x + 1]; ++arc) {
            int y = ch.forward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[x].push_back(make_pair(y, ch.forward.weight[arc]));
            back_adj[y].push_back(make_pair(x, ch.forward.weight[arc]));
        }
        for (int arc = ch.backward.first_out[x]; arc < ch.backward.first_out[x + 1]; ++arc) {
            int y = ch.backward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[y].push_back(make_pair(x, ch.backward.weight[arc]));
            back_adj[x].push_back(make_pair(y, ch.backward.weight[arc]));
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
        g.forward.first_out.push_back(g.forward.head.size());
        g.backward.first_out.push_back(g.backward.head.size());
        g.latitude.push_back(graph.latitude[ch.order[x]]);
        g.longitude.push_back(graph.longitude[ch.order[x]]);
        g.order.push_back(ch.order[x]);
        for (auto& pair : adj[x]) {
            int y = pair.first;
            int w = pair.second;
            g.forward.head.push_back(id_map[y]);
            g.forward.weight.push_back(w);
        }
        for (auto& pair : back_adj[x]) {
            int y = pair.first;
            int w = pair.second;
            g.backward.head.push_back(id_map[y]);
            g.backward.weight.push_back(w);
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
            if (x < y) {
                g.forward_up.head.push_back(y);
                g.forward_up.weight.push_back(g.forward.weight[arc]);
            } else {
                g.forward_down.head.push_back(y);
                g.forward_down.weight.push_back(g.forward.weight[arc]);
            }
        }
        for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x + 1]; ++arc) {
            int y = g.backward.head[arc];
            if (x < y) {
                g.backward_up.head.push_back(y);
                g.backward_up.weight.push_back(g.forward.weight[arc]);
            } else {
                g.backward_down.head.push_back(y);
                g.backward_down.weight.push_back(g.backward.weight[arc]);
            }
        }
    }
    g.forward_up.first_out.push_back(g.forward_up.head.size());
    g.backward_up.first_out.push_back(g.backward_up.head.size());
    g.forward_down.first_out.push_back(g.forward_down.head.size());
    g.backward_down.first_out.push_back(g.backward_down.head.size());
}

vector<int> partition_graph(SimpleOSMCarRoutingGraph& graph, ContractionHierarchy& ch, Graph& g, int nparts, float best_percentage) {
    int n = g.node_count();

    vector<SCOTCH_Num> xadj;
    vector<SCOTCH_Num> adjncy;
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

    SCOTCH_Graph grafdat;
    SCOTCH_Strat stradat;
    SCOTCH_Num baseval;
    SCOTCH_Num vertnbr;
    int o;

    SCOTCH_graphInit(&grafdat);
    baseval = 0;
    vertnbr = n;
    o = 1;
    vector<SCOTCH_Num> part(n);
    if (SCOTCH_graphBuild(&grafdat, baseval, vertnbr, &(xadj[0]), &(xadj[1]), NULL, NULL, adjncy.size(), &(adjncy[0]), NULL) == 0) {
        SCOTCH_stratInit(&stradat);
        if (SCOTCH_graphCheck(&grafdat) == 0) {
            o = SCOTCH_graphPart(&grafdat, nparts, &stradat, &(part[0]));
        }
        SCOTCH_stratExit(&stradat);
    }
    SCOTCH_graphExit(&grafdat);
    vector<int> partition(n);
    for (int i = 0; i < n; ++i) partition[i] = part[i];
    return partition;
}

int main(int argc, const char* argv[]) {
    try {
        namespace po = boost::program_options;
        po::options_description description("ParkPlacement");
        description.add_options()("help,h", "Display this help message")("graph,g", po::value<string>(), "Graph pbf file")("core,c", po::value<float>()->default_value(1.0), "Number in [0,1] describing the core size")("partition,p", po::value<int>()->default_value(100), "Partition size");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            std::cout << description;
            return 0;
        }
        if (!vm.count("graph")) {
            cout << "please provide an input graph!" << endl;
            return 0;
        }
        string pbf_file = vm["graph"].as<string>();

        cout << "start loading graph" << endl;
        auto graph = simple_load_osm_car_routing_graph_from_pbf(pbf_file);
        graph = removeSmallCCs(graph, findCCs(graph));
        cout << "graph loaded!" << endl;
        cout << "start contracting graph" << endl;

        auto tail = invert_inverse_vector(graph.first_out);

        // Build the shortest path index
        cout << "start building contraction hierarchy" << endl;
        string ch_save = pbf_file + ".ch";
        auto ch = boost::filesystem::exists(ch_save) ? ContractionHierarchy::load_file(ch_save) : ContractionHierarchy::build(graph.node_count(), tail, graph.head, graph.travel_time, [](string msg) { std::cout << msg << endl; });
        if (!boost::filesystem::exists(ch_save)) ch.save_file(ch_save);
        cout << "contraction hierarchy finished!" << endl;

        float core = max(0.0f, min(1.0f, vm["core"].as<float>()));
        Graph g = build_ch_complete_graph("ch_core", graph, ch, core);
        build_up_down_cores(g);

        cout << "ch search graph build" << endl;

        cout << "start partition: " << g.node_count() << " nodes, " << g.forward.head.size() << " edges" << endl;
        auto partition = partition_graph(graph, ch, g, vm["partition"].as<int>(), core);
        g.compute_boundary_node(partition);
        ArcFlags arcflags = ArcFlags(g, partition, vm["partition"].as<int>());
        arcflags.precompute(0, vm["partition"].as<int>());

    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}