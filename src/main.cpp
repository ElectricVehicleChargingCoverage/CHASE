#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <scotch.h>
#include <arcflags.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

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
    getline(file, line); // skip first line
    while (getline(file, line)) {
        vector<string> csv = split(line, ",", false);
        int cost = stod(csv[5]) * 3600.0 * 1000;
        int reverse_cost = stod(csv[6]) * 3600.0 * 1000;
        long long source_id = stoll(csv[1]);
        long long target_id = stoll(csv[2]);
        double from_lng = stod(csv[7]); double from_lat = stod(csv[8]);
        double to_lng = stod(csv[9]); double to_lat = stod(csv[10]);
        nodes[source_id] = make_pair(from_lat, from_lng);
        nodes[target_id] = make_pair(to_lat, to_lng);
        if (source_id != target_id) {
            adj[source_id].push_back(make_pair(target_id, cost));
            if(stod(csv[6]) < 1000000.0)
                adj[target_id].push_back(make_pair(source_id, reverse_cost));
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
        Info(unsigned _y, unsigned _arc, unsigned _weight): y{_y}, arc{_arc}, weight{_weight} {}
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
            if (x < y) {
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
            if (x < y) {
                g.backward_up.head.push_back(y);
                g.backward_up.weight.push_back(g.forward.weight[arc]);
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
        auto graph = loadGraph(pbf_file + "/edges.csv"); //simple_load_osm_car_routing_graph_from_pbf(pbf_file);

        graph = removeSmallCCs(graph, findCCs(graph));
        cout << "graph loaded!" << endl;
        cout << "start contracting graph" << endl;

        auto tail = invert_inverse_vector(graph.forward.first_out);

        // Build the shortest path index
        cout << "start building contraction hierarchy" << endl;
        string ch_save = pbf_file + ".ch";
        auto ch = boost::filesystem::exists(ch_save) ? ContractionHierarchy::load_file(ch_save) : ContractionHierarchy::build(
            graph.node_count(), 
            tail, 
            graph.forward.head, 
            graph.forward.weight, 
            [](string msg) { std::cout << msg << endl; 
        });
        if (!boost::filesystem::exists(ch_save)) ch.save_file(ch_save);
        cout << "contraction hierarchy finished!" << endl;

        float core = max(0.0f, min(1.0f, vm["core"].as<float>()));

        vector<string> csv = split(pbf_file, "/", false);
        string name = csv.back() + "_core_" + to_string(core);
        Graph g = build_ch_complete_graph(name, graph, ch, core);
        build_up_down_cores(g);

        cout << "ch search graph build" << endl;
        
        cout << "start partition: " << g.node_count() << " nodes, " << g.forward.head.size() << " edges" << endl;
        auto partition = partition_graph(g, vm["partition"].as<int>(), core);
        g.compute_boundary_node(partition);
        ofstream out("partition.txt");
        for (int i = 0; i < partition.size(); ++i) {
            out << partition[i] << endl;
        }
        out.close();
        ArcFlags flags(g, partition, vm["partition"].as<int>());
        flags.precompute(0, vm["partition"].as<int>());

        // ContractionHierarchyQuery query(ch);
        // query.edge_hashes = flags.label_hashes;
        // query.hash_flags = flags.labels;

    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}