#include <bits/stdc++.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <routingkit/osm_simple.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <arcflags.hpp>
#include <metis.h>
#include "kaHIP_interface.h"

using namespace RoutingKit;
using namespace std;

vector<vector<int>> ccList;


SimpleOSMCarRoutingGraph removeSmallCCs(SimpleOSMCarRoutingGraph& graph, int minSize) {
    cout << "Graph size before cleanup: ";
    cout << graph.node_count() << " nodes, " << graph.arc_count() << " edges" << endl;
    vector<bool> deleted(graph.node_count(), false);
    for (auto cc : ccList)
        if (cc.size() < minSize) for (int node : cc) deleted[node] = true;

    vector<vector<pair<int, int>>> adj(graph.node_count());
    for (int x = 0; x < graph.node_count(); ++x) {
        if (deleted[x]) continue;
        for (int arc = graph.first_out[x]; arc < graph.first_out[x+1]; ++arc) {
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
 
    ofstream out("graph2.csv");
    out << "from_lat,from_lng,to_lat,to_lng" << endl;
    for (int x = 0; x < cleaned.node_count(); ++x) {
        int d = 0;
        for (int arc = cleaned.first_out[x]; arc < cleaned.first_out[x+1]; ++arc) {
            int y = cleaned.head[arc];
            out << cleaned.latitude[x] << "," << cleaned.longitude[x] << "," << cleaned.latitude[y] << "," << cleaned.longitude[y] << endl;
            d++;
        }
        assert(d <= 7);
    }
    out.close();

    cout << "Graph size after cleanup: ";
    cout << cleaned.node_count() << " nodes, " << cleaned.arc_count() << " edges" << endl;
    return cleaned;
}

void bfs(SimpleOSMCarRoutingGraph& graph, int src, vector<bool>& visited, int comp) {
    visited[src] = true;
    queue<int> q{ {{src}} };
    while (!q.empty()) {
        auto v = q.front();
        ccList[comp].push_back(v);
        q.pop();
        for (int arc = graph.first_out[v]; arc < graph.first_out[v+1]; ++arc) {
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

Graph build_ch_complete_graph(string name, SimpleOSMCarRoutingGraph& graph, ContractionHierarchy& ch, float best_percentage){
    vector<vector<pair<int, int>>> adj(ch.node_count());
    vector<vector<pair<int, int>>> back_adj(ch.node_count());
    assert(ch.forward.first_out.size() == ch.backward.first_out.size());

    for (int x = 0; x < ch.node_count(); ++x ){
        if (x < ch.node_count() * (1-best_percentage)) continue;
        for (int arc = ch.forward.first_out[x]; arc < ch.forward.first_out[x+1]; ++arc) {
            int y = ch.forward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[x].push_back(make_pair(y, ch.forward.weight[arc]));
            back_adj[y].push_back(make_pair(x, ch.forward.weight[arc]));
        }
        for (int arc = ch.backward.first_out[x]; arc < ch.backward.first_out[x+1]; ++arc) {
            int y = ch.backward.head[arc];
            if (y < ch.node_count() * (1 - best_percentage)) continue;
            adj[y].push_back(make_pair(x, ch.backward.weight[arc]));
            back_adj[x].push_back(make_pair(y, ch.backward.weight[arc]));
        }
    }
    unordered_map<int, int> id_map;
    vector<int> nodes;
    for (int x = 0, count = 0; x < graph.node_count(); ++x) {
        if (x < ch.node_count() * (1-best_percentage)) continue;
        id_map[x] = count;
        count++;
        nodes.push_back(x);
    }
    Graph g = Graph(nodes.size(), name);
    for (int x : nodes) {
        g.forward.first_out.push_back(g.forward.head.size());
        g.backward.first_out.push_back(g.backward.head.size());
        g.latitude.push_back(graph.latitude[x]);
        g.longitude.push_back(graph.longitude[x]);
        g.order.push_back(ch.order[x]);
        for (auto& pair : adj[x]) {
            int y = pair.first; int w = pair.second;
            g.forward.head.push_back(id_map[y]);
            g.forward.weight.push_back(w);
        }
        for (auto& pair : back_adj[x]) {
            int y = pair.first; int w = pair.second;
            g.backward.head.push_back(id_map[y]);
            g.backward.weight.push_back(w);
        }
    }
    assert(g.order.size() == g.node_count());
    cout << g.node_count() << ", " << ch.node_count() << " -> " << (double)g.node_count() / (double)ch.node_count() << endl;
    g.forward.first_out.push_back(g.forward.head.size());
    g.backward.first_out.push_back(g.backward.head.size());
    return g;
}

vector<idx_t> paritition_graph(SimpleOSMCarRoutingGraph& graph, ContractionHierarchy& ch, Graph& g, idx_t nparts) {
    vector<int> xadj = g.forward.first_out;
    vector<int> adjncy = g.forward.head;

    cout << xadj.size() << endl;
    cout <<adjncy.size() << endl;

    // idx_t options[METIS_NOPTIONS];
    // METIS_SetDefaultOptions(options);
    // options[METIS_OPTION_NUMBERING] = 0;

    int edgecut = 0;
    int ncon = 1;
    int n = xadj.size()-1;
    double imbalance = 0.03;
    vector<int> partitions;
    partitions.resize(n);

    cout << "start partitioning" << endl;
    kaffpa(&n, NULL, &(xadj[0]), NULL, &(adjncy[0]), &nparts, &imbalance, false, 0, ECO, &edgecut, &(partitions[0]));
    // auto ios = METIS_PartGraphRecursive(&n, &ncon, &(xadj[0]), &(adjncy[0]), NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, &(partitions[0]));
    // auto ios = METIS_PartGraphKway(&n, &ncon, &(xadj[0]), &(adjncy[0]), NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, &(partitions[0]));
    // cout << "finished partitioning into "<< nparts << " parts" << endl;
    // if (ios != METIS_OK) {
    //     cout << "METIS failed with error: " << ios << endl;
    // }else {
    //     cout << "METIS succeeded" << endl;
    // }

    return partitions;
}

int main(int argc, const char* argv[]) {
    try {
        namespace po = boost::program_options;
        po::options_description description("ParkPlacement");
        description.add_options()
            ("help,h", "Display this help message")
            ("graph,g", po::value<string>(), "Graph pbf file")
            ("core,c", po::value<float>()->default_value(1.0), "Number in [0,1] describing the core size")
            ("partition,p", po::value<int>()->default_value(100), "Partition size");
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
        auto tail = invert_inverse_vector(graph.first_out);

        // Build the shortest path index
        cout << "start building contraction hierarchy" << endl;
        string ch_save = pbf_file + ".ch";
        auto ch = boost::filesystem::exists(ch_save)? ContractionHierarchy::load_file(ch_save) : ContractionHierarchy::build(
            graph.node_count(), 
            tail, graph.head, 
            graph.travel_time,
            [](string msg) { std::cout << msg << endl; }
        );
        if (!boost::filesystem::exists(ch_save)) ch.save_file(ch_save);
        cout << "contraction hierarchy finished!" << endl;

        float core = max(0.0f, min(1.0f, vm["core"].as<float>()));
        Graph g = build_ch_complete_graph("test", graph, ch, core);
        cout << "ch search graph build" << endl;

        auto partition = paritition_graph(graph, ch, g, vm["partition"].as<int>());
        ofstream out("nodes.csv");
        out << "lat,lng,cell" << endl;
        for (int x = 0; x < ch.node_count(); ++x) {
            out << graph.latitude[x] << "," << graph.longitude[x] << "," << partition[x] << endl;
        }
        out.close();
        // // ArcFlags arcflags = ArcFlags(g);

        // // ContractionHierarchyQuery ch_query(ch);
        
    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}