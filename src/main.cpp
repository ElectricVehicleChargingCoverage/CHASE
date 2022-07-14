#include <bits/stdc++.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <routingkit/osm_simple.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <arcflags.hpp>
#include <metis.h>

using namespace RoutingKit;
using namespace std;

Graph build_ch_complete_graph(string name, SimpleOSMCarRoutingGraph& graph, ContractionHierarchy& ch, float best_percentage){
    Graph g = Graph(graph.node_count(), name);
    vector<vector<pair<int, int>>> adj(ch.node_count());
    vector<vector<pair<int, int>>> back_adj(ch.node_count());
    for (int x = 0; x < ch.node_count(); ++x ){
        if (x < ch.node_count() * best_percentage) continue;
        for (int arc = ch.forward.first_out[x]; arc < ch.forward.first_out[x+1]; ++arc) {
            int y = ch.forward.head[arc];
            adj[x].push_back(make_pair(y, ch.forward.weight[arc]));
            back_adj[y].push_back(make_pair(x, ch.forward.weight[arc]));
        }
        for (int arc = ch.backward.first_out[x]; arc < ch.backward.first_out[x+1]; ++arc) {
            int y = ch.backward.head[y];
            adj[y].push_back(make_pair(x, ch.forward.weight[arc]));
            back_adj[x].push_back(make_pair(y, ch.forward.weight[arc]));
        }
    }
    for (int x = 0; x < ch.node_count(); ++x) {
        g.forward.first_out.push_back(g.forward.head.size());
        g.backward.first_out.push_back(g.backward.head.size());
        for (auto& pair : adj[x]) {
            int y = pair.first; int w = pair.second;
            g.forward.head.push_back(y);
            g.forward.weight.push_back(w);
        }
        for (auto& pair : back_adj[x]) {
            int y = pair.first; int w = pair.second;
            g.backward.head.push_back(y);
            g.backward.weight.push_back(w);
        }
    }
    g.forward.first_out.push_back(g.forward.head.size());
    g.backward.first_out.push_back(g.backward.head.size());
    return g;
}

void paritition_graph(Graph& g, int nparts = 20) {

    // vector<int> first_out(g.node_count()+1);
    // vector<int> head;
    // vector<int> cell(graph.node_count());
    // vector<int> weight;

    // first_out[0] = 1;
    // int pos = 0;
    // for (int x = 0; x < graph.node_count(); ++x) {
    //     for (auto& pair : adj[x]) {
    //         head.push_back(pair.first);
    //         weight.push_back(pair.second);
    //         ++pos;
    //     }
    //     first_out[x+1] = pos+1;
    // }

    vector<int> xadj = g.forward.first_out;//{0, 2, 5, 8, 11, 13 ,16, 20, 24, 28, 31, 33, 36, 39, 42, 44};
    vector<int> adjncy =  g.forward.head; //{1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14 ,5 ,11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};
    vector<int> w(adjncy.size(), 1);
    assert(xadj.back() == adjncy.size());

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    int edgecut = 0;
    int ncon = 1;
    int n = xadj.size()-1;

    vector<int> partitions;
    partitions.resize(n);

    cout << "start partitioning" << endl;
    METIS_PartGraphRecursive(&n, &ncon, &(xadj[0]), &(adjncy[0]), NULL, NULL,NULL, &nparts, NULL, NULL, options, &edgecut, &(partitions[0]));
    cout << "finished partitioning into "<< nparts << " parts" << endl;
}

int main(int argc, const char* argv[]) {
    try {
        namespace po = boost::program_options;
        po::options_description description("ParkPlacement");
        description.add_options()
            ("help,h", "Display this help message")
            ("graph,g", po::value<string>(), "graph pbf file");
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

        cout << graph.node_count() << "," << graph.head.size() << endl;
        cout
        // Graph g = build_ch_complete_graph("test", graph, ch, 0.05);
        // cout << g.forward.head.size() << endl;
        // paritition_graph(g);
        // ArcFlags arcflags = ArcFlags(g);
        
    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}