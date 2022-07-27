#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <arcflags.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace RoutingKit;

void find_shortest_path_backwards(Graph& g, unsigned src, unsigned trgt, int cell_idx, ArcFlags& flags, ContractionHierarchy& ch) {
    vector<unsigned> tentative_dist(g.node_count(), inf_weight);
    vector<bool> down(g.node_count(), false);
    vector<bool> pop(g.node_count(), false);
    vector<long long> pred(g.node_count(), -1);
    vector<int> pred_n(g.node_count(), -1);

    MinIDQueue queue(g.node_count());
    TimestampFlags was_pushed(g.node_count());
    queue.push({trgt, 0});

    while (!queue.empty()) {
        auto popped = queue.pop();
        int v = popped.id;
        pop[v] = true;
        auto distance_to_popped_node = popped.key;
        if (v == src) {
            break;
        }
        if (v == trgt || !down[v]){
            for (int arc = g.backward_up.first_out[v]; arc < g.backward_up.first_out[v+1]; ++arc) {
                unsigned u = g.backward_up.head[arc];
                if (pop[u]) continue;
                unsigned d = distance_to_popped_node + g.backward_up.weight[arc];

                long long orig_arc = g.backward_up.original_arc[arc];
                unsigned w = (orig_arc >= 0)? ch.forward.weight[(unsigned)orig_arc] : ch.backward.weight[-(unsigned)orig_arc];
                assert(g.backward_up.weight[arc] == w);

                if(was_pushed.is_set(u)){
                    if (d < tentative_dist[u]) {
                        pred[u] = g.backward_up.original_arc[arc];
                        pred_n[u] = v;
                        down[u] = false;
                        queue.decrease_key({u, d});
                        tentative_dist[u] = d;
                    }
                }else if(d < inf_weight){
                    pred[u] = g.backward_up.original_arc[arc];
                    pred_n[u] = v;
                    down[u] = false;
                    was_pushed.set(u);
                    queue.push({u, d});
                    tentative_dist[u] = d;
                }
            }
        }
        for (int arc = g.backward_down.first_out[v]; arc < g.backward_down.first_out[v+1]; ++arc) {
            unsigned u = g.backward_down.head[arc];
            if (pop[u]) continue;
            unsigned d = distance_to_popped_node + g.backward_down.weight[arc];

            long long orig_arc = g.backward_down.original_arc[arc];
            unsigned w = (orig_arc >= 0)? ch.forward.weight[(unsigned)orig_arc] : ch.backward.weight[-(unsigned)orig_arc];
            assert(g.backward_down.weight[arc] == w);

            if(was_pushed.is_set(u)){
                if (d < tentative_dist[u]) {
                    pred[u] = g.backward_down.original_arc[arc];
                    pred_n[u] = v;
                    down[u] = true;
                    queue.decrease_key({u, d});
                    tentative_dist[u] = d;
                }
            }else if(d < inf_weight){
                pred[u] = g.backward_down.original_arc[arc];
                pred_n[u] = v;
                down[u] = true;
                was_pushed.set(u);
                tentative_dist[u] = d;
                queue.push({u, d});
            }
        }
    }
    cout << "dist: " << tentative_dist[src] << endl;
    list<long long> path;
    int x = src;
    while (pred_n[x] != -1) {
        long long arc = pred[x];
        // if (flags.partition[x] == cell_idx) cout << "bn: " << x << " " << arc << " bn:" << g.boundary_nodes[x] << endl;
        cout << "x: " << x << " " << flags.partition[x] << " " << g.boundary_nodes[x] << endl;
        path.push_back(pred[x]);
        x = pred_n[x];
    }
    cout << "----- arc path ---------" << endl;
    auto tail_f = invert_inverse_vector(ch.forward.first_out);
    auto tail_b = invert_inverse_vector(ch.backward.first_out);
    for (long long arc : path) {
        unsigned x,y,w;
        if (arc >= 0){
            w = ch.forward.weight[(unsigned)arc];
            x = g.new_id[tail_f[(unsigned)arc]];
            y = g.new_id[ch.forward.head[(unsigned)arc]];
        }else {
            w = ch.backward.weight[(unsigned)(-arc)];
            y = g.new_id[tail_b[(unsigned)(-arc)]];
            x = g.new_id[ch.backward.head[(unsigned)(-arc)]];
        }
        // unsigned w = (arc >= 0)? ch.forward.weight[(unsigned)arc] : ch.backward.weight[-(unsigned)arc];
        cout << arc << " " << w << ": " << x << "[" << flags.partition[x] << "] " << y << "[" << flags.partition[y] << "]" << endl;
    }
}

void find_shortest_path(Graph& g, unsigned src, unsigned trgt, int cell_idx, ArcFlags& flags) {
    vector<unsigned> tentative_dist(g.node_count(), inf_weight);
    vector<bool> down(g.node_count(), false);
    vector<bool> pop(g.node_count(), false);
    vector<long long> pred(g.node_count(), -1);
    vector<unsigned> pred_n(g.node_count(), -1);

    MinIDQueue queue(g.node_count());
    TimestampFlags was_pushed(g.node_count());
    queue.push({src, 0});
    
    while (!queue.empty()) {
        auto popped = queue.pop();
        int v = popped.id;
        pop[v] = true;
        auto distance_to_popped_node = popped.key;
        if (v == trgt) {
            break;
        }
        if (v == src || !down[v]){
            for (int arc = g.forward_up.first_out[v]; arc < g.forward_up.first_out[v+1]; ++arc) {
                // if (flags.labels[flags.label_hashes[g.forward_up.original_arc[arc]]][cell_idx] == 0) continue;
                unsigned u = g.forward_up.head[arc]; 
                if (pop[u]) continue;
                unsigned d = distance_to_popped_node + g.forward_up.weight[arc];
                if(was_pushed.is_set(u)){
                    if (d < tentative_dist[u]) {
                        pred[u] = g.forward_up.original_arc[arc];
                        pred_n[u] = v;
                        down[u] = false;
                        queue.decrease_key({u, d});
                        tentative_dist[u] = d;
                    }
                }else if(d < inf_weight){
                    pred[u] = g.forward_up.original_arc[arc];
                    pred_n[u] = v;
                    down[u] = false;
                    was_pushed.set(u);
                    queue.push({u, d});
                    tentative_dist[u] = d;
                }
            }
        }
        for (int arc = g.forward_down.first_out[v]; arc < g.forward_down.first_out[v+1]; ++arc) {
            unsigned u = g.forward_down.head[arc];
            if (pop[u]) continue;
            unsigned d = distance_to_popped_node + g.forward_down.weight[arc];
            if(was_pushed.is_set(u)){
                if (d < tentative_dist[u]) {
                    pred[u] = g.forward_down.original_arc[arc];
                    pred_n[u] = v;
                    down[u] = true;
                    queue.decrease_key({u, d});
                    tentative_dist[u] = d;
                }
            }else if(d < inf_weight){
                pred[u] = g.forward_down.original_arc[arc];
                pred_n[u] = v;
                down[u] = true;
                was_pushed.set(u);
                queue.push({u, d});
                tentative_dist[u] = d;
            }
        }
    }
    cout << tentative_dist[trgt] << endl;
    list<long long> path;
    int x = trgt;
    while (pred_n[x] != -1) {
        if (flags.partition[x] == cell_idx && g.boundary_nodes[x]) cout << "bn: " << x << endl;
        cout << "x: " << x << endl;
        long long arc = pred[x];
        path.push_front(pred[x]);
        x = pred_n[x];
    }
    cout << "----- arc path ---------" << endl;
    for (long long arc : path) {
        cout << arc << endl;
    }
}

void start_experiments(int tries, ContractionHierarchy& ch, Graph& g, ArcFlags& flags, unsigned min_rank) {
    ContractionHierarchyQuery query(ch);
    vector<int> parts(ch.node_count(), -1);
    for (int i = 0; i < ch.node_count(); ++i) {
        if (g.new_id.find(i) != g.new_id.end()) {
            parts[i] = flags.partition[g.new_id[i]];
        }
    }
    query.partition = parts;
    query.partition_size = flags.partition_size;
    query.edge_hashes = flags.label_hashes;
    query.hash_flags = flags.labels;


    // run test queries
    string seed_str = "very_random_seed";
    seed_seq seed (seed_str.begin(), seed_str.end());
    mt19937 gen(seed);
    uniform_int_distribution<> dist(0, g.node_count()-1);

    vector<pair<unsigned, unsigned>> routes(tries);
    vector<unsigned> ch_distances(tries);
    vector<unsigned> chase_distances(tries);
    vector<int64_t> ch_time(tries);
    vector<int64_t> chase_time(tries);
    vector<pair<unsigned, unsigned>> metric_ch(tries);
    vector<pair<unsigned, unsigned>> metric_chase(tries);
    
    auto f_tail = invert_inverse_vector(ch.forward.first_out);
    auto b_tail = invert_inverse_vector(ch.backward.first_out);
    for (int i = 0; i < tries; ++i) {
        cout << "ch i:" << i+1 << "/" << tries << endl; 
        int x = g.order[dist(gen)];
        int y = g.order[dist(gen)];
        routes[i] = make_pair(x,y);

        query.reset().add_source(x).add_target(y);
        auto start = chrono::high_resolution_clock::now();
        query.run();
        auto finish = chrono::high_resolution_clock::now();

        ch_time[i] = chrono::duration_cast<chrono::nanoseconds>(finish-start).count();
        ch_distances[i] = query.get_distance();
        metric_ch[i] = make_pair(query.relaxed, query.visited);
    }
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    auto time = oss.str();
    ofstream out("../logs/log_" + time + ".txt");

    int f = 0;
    for (int i = 0; i < tries; ++i) {
        auto [x, y] = routes[i];

        query.reset().add_source(x).add_target(y);
        auto start = chrono::high_resolution_clock::now();
        query.run_chase(min_rank);
        auto finish = chrono::high_resolution_clock::now();
        
        chase_time[i] = chrono::duration_cast<chrono::nanoseconds>(finish-start).count();
        chase_distances[i] = query.get_distance();
        metric_chase[i] = make_pair(query.relaxed, query.visited);

        if (chase_distances[i] != ch_distances[i]) f++;
        cout << "[" << (chase_distances[i] == ch_distances[i]) << "] i:" << i+1 << "/" << tries << endl; 
        out << (chase_distances[i] == ch_distances[i]) << "," << chase_distances[i] << "," << ch_distances[i] << "," << metric_ch[i].first << "," << query.relaxed << "," << metric_ch[i].second << "," << query.visited << "," << ch_time[i] << "," << chase_time[i] << endl;
    }
    cout << "failed: " << f << endl;
    out.close();
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

        string input_path = "../data/" + vm["graph"].as<string>();

        cout << "start loading graph" << endl;
        auto graph = loadGraph(input_path + "/edges.csv");

        graph = removeSmallCCs(graph, findCCs(graph));
        cout << "graph loaded!" << endl;
        cout << "start contracting graph" << endl;

        auto tail = invert_inverse_vector(graph.forward.first_out);
        float core = max(0.0f, min(1.0f, vm["core"].as<float>()));

        vector<string> csv = split(input_path, "/", false);
        string name = csv.back() + "_core_" + to_string(core);

        // Build the shortest path index
        cout << "start building contraction hierarchy" << endl;
        string ch_save = input_path + "/" + name + ".ch";
        auto ch = boost::filesystem::exists(ch_save) ? ContractionHierarchy::load_file(ch_save) : ContractionHierarchy::build(
            graph.node_count(), 
            tail, 
            graph.forward.head, 
            graph.forward.weight, 
            [](string msg) { std::cout << msg << endl; 
        });
        if (!boost::filesystem::exists(ch_save)) ch.save_file(ch_save);
        cout << "contraction hierarchy finished!" << endl;

        Graph g = build_ch_complete_graph(name, graph, ch, core);
        build_up_down_cores(g);
        unsigned min_rank = ch.node_count() * (1 - core);
        cout << "ch search graph build" << endl;
        
        cout << "start partition: " << g.node_count() << " nodes, " << g.forward.head.size() << " edges" << endl;
        int partition_size = vm["partition"].as<int>();
        string partition_location = input_path + "/" + name + "_" + to_string(partition_size) + ".part";
        auto partition = boost::filesystem::exists(partition_location)? read_partition(partition_location) : partition_graph(g, vm["partition"].as<int>(), core);
        if (!boost::filesystem::exists(partition_location)) 
            export_partition(partition, partition_location);

        g.compute_boundary_node(partition);

        ArcFlags flags(g, partition, partition_size);
        string import_file = "../flags/" + g.name + "_" + to_string(partition_size);
        flags.importFlags(import_file + ".csv", import_file + ".bin");
        // flags.precompute(0, vm["partition"].as<int>());
        // flags.mergeFlags();

        start_experiments(10000, ch, g, flags, min_rank);
    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}