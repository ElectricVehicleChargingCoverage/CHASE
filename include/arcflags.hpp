#pragma once
#include <bits/stdc++.h>
#include <thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <routingkit/id_queue.h>
#include <routingkit/timestamp_flag.h>
#include "stringutil.hpp"

using namespace std;
using namespace RoutingKit;

struct Graph {
private:
    int count;
public:
    Graph(int _count, string _name=""): count{_count}, name{_name} {}
    struct Side{
        vector<unsigned>first_out;
        vector<unsigned>head;
        vector<unsigned>weight;
        vector<long long>original_arc;
    };
    string name;
	Side forward, backward;
    Side forward_up, forward_down, backward_up, backward_down;
    vector<float> latitude, longitude;
    vector<int> order;
    unordered_map<int, int> new_id;
    vector<bool> boundary_nodes;
    int node_count() { return count; };
    void compute_boundary_node(vector<int> partition){
        boundary_nodes.clear(); boundary_nodes.resize(count);
        assert(partition.size() == count);
        for (int x = 0; x < count; ++x) {
            for (int arc = forward.first_out[x]; arc < forward.first_out[x+1]; ++arc) {
                int y = forward.head[arc];
                if (partition[x] != partition[y])
                    boundary_nodes[x] = true;
            }
        }
    }
};

class ArcFlags {
public:
    Graph& g;
    vector<int> partition;
    int partition_size;
    ArcFlags(Graph& _g, vector<int>& _partition, int _partition_size): g{_g}, partition{_partition}, partition_size{_partition_size} {};
    void precompute(int start, int end);
    void mergeFlags(string folder);
    void compress(unordered_map<long long, boost::dynamic_bitset<>>& labels);
    void exportFlags(string folder);
    void importFlags(string edges_path, string flags_path);
    unordered_map<size_t, boost::dynamic_bitset<>> preprocessing_labels;
    unordered_map<long long, size_t> preprocessing_label_hash;
    unordered_map<long long, size_t> label_hashes;
    unordered_map<size_t, boost::dynamic_bitset<>> labels;
    bool precomputed = false;
};

void ArcFlags::precompute(int start, int end) {
    thread_pool pool; // amount of threads
    synced_stream sync_out;
    unordered_map<long long, bool> cell_maps_arc_flags[2 * partition_size];

    using QueueElement = pair<int, double>;
    auto cmp = [](QueueElement &a, QueueElement &b) {
        return a.second > b.second;
    };
    auto markEdgesOnSptTo = [&](unsigned src, int cell_idx) {
        vector<unsigned> tentative_dist(g.node_count(), inf_weight);
        vector<bool> down(g.node_count(), false);
        vector<long long> pred(g.node_count(), -1);

        MinIDQueue queue(g.node_count());
        TimestampFlags was_pushed(g.node_count());
        queue.push({src, 0});

        while (!queue.empty()) {
            auto popped = queue.pop();
            int v = popped.id;
            auto distance_to_popped_node = popped.key;
            if (pred[v] != -1) {
                cell_maps_arc_flags[cell_idx][pred[v]] = true; 
            }
            if (v == src || !down[v]){
                for (int arc = g.backward_up.first_out[v]; arc < g.backward_up.first_out[v+1]; ++arc) {
                    unsigned u = g.backward_up.head[arc]; 
                    unsigned d = distance_to_popped_node + g.backward_up.weight[arc];
                    if(was_pushed.is_set(u)){
                        if (d < tentative_dist[u]) {
                            pred[u] = g.backward_up.original_arc[arc];
                            down[u] = false;
                            queue.decrease_key({u, d});
                            tentative_dist[u] = d;
                        }
                    }else if(d < inf_weight){
                        pred[u] = g.backward_up.original_arc[arc];
                        down[u] = false;
                        was_pushed.set(u);
                        queue.push({u, d});
                        tentative_dist[u] = d;
                    }
                }
            }
            for (int arc = g.backward_down.first_out[v]; arc < g.backward_down.first_out[v+1]; ++arc) {
                unsigned u = g.backward_down.head[arc];
                unsigned d = distance_to_popped_node + g.backward_down.weight[arc];
                if(was_pushed.is_set(u)){
                    if (d < tentative_dist[u]) {
                        pred[u] = g.backward_down.original_arc[arc];
                        down[u] = true;
                        queue.decrease_key({u, d});
                        tentative_dist[u] = d;
                    }
                }else if(d < inf_weight){
                    pred[u] = g.backward_down.original_arc[arc];
                    down[u] = true;
                    was_pushed.set(u);
                    tentative_dist[u] = d;
                    queue.push({u, d});
                }
            }
        }
    };

    auto markEdgesOnSptFrom = [&](unsigned src, int cell_idx) {
        vector<unsigned> tentative_dist(g.node_count(), inf_weight);
        vector<bool> down(g.node_count(), false);
        vector<long long> pred(g.node_count(), -1);

        MinIDQueue queue(g.node_count());
        TimestampFlags was_pushed(g.node_count());
        queue.push({src, 0});
        
        while (!queue.empty()) {
            auto popped = queue.pop();
            int v = popped.id;
            auto distance_to_popped_node = popped.key;
            if (pred[v] != -1) {
                cell_maps_arc_flags[cell_idx + partition_size][pred[v]] = true; 
            }
            if (v == src || !down[v]){
                for (int arc = g.forward_up.first_out[v]; arc < g.forward_up.first_out[v+1]; ++arc) {
                    unsigned u = g.forward_up.head[arc]; 
                    unsigned d = distance_to_popped_node + g.forward_up.weight[arc];
                    if(was_pushed.is_set(u)){
                        if (d < tentative_dist[u]) {
                            pred[u] = g.forward_up.original_arc[arc];
                            down[u] = false;
                            queue.decrease_key({u, d});
                            tentative_dist[u] = d;
                        }
                    }else if(d < inf_weight){
                        pred[u] = g.forward_up.original_arc[arc];
                        down[u] = false;
                        was_pushed.set(u);
                        queue.push({u, d});
                        tentative_dist[u] = d;
                    }
                }
            }
            for (int arc = g.forward_down.first_out[v]; arc < g.forward_down.first_out[v+1]; ++arc) {
                unsigned u = g.forward_down.head[arc];
                unsigned d = distance_to_popped_node + g.forward_down.weight[arc];
                if(was_pushed.is_set(u)){
                    if (d < tentative_dist[u]) {
                        pred[u] = g.forward_down.original_arc[arc];
                        down[u] = true;
                        queue.decrease_key({u, d});
                        tentative_dist[u] = d;
                    }
                }else if(d < inf_weight){
                    pred[u] = g.forward_down.original_arc[arc];
                    down[u] = true;
                    was_pushed.set(u);
                    queue.push({u, d});
                    tentative_dist[u] = d;
                }
            }
        }
    };

    auto saveFlags = [&](int cell_idx, const unordered_map<long long, bool> &map) {
        ofstream out("../flags/" + g.name + "_" + to_string(partition_size) + "_" + to_string(cell_idx));
        for (auto [e, _] : map)
            if (_) out << e << endl;
        out.close();
    };

    auto precomputeCell = [&](int cell_idx) {
        sync_out.println("Start computation for cell ", cell_idx);
        for (int x = 0; x < g.node_count(); ++x) {
            if (g.boundary_nodes[x]){
                markEdgesOnSptTo(x, cell_idx);
                markEdgesOnSptFrom(x, cell_idx);
            }
            for (int arc = g.forward.first_out[x]; arc < g.forward.first_out[x+1]; ++arc) {
                int y = g.forward.head[arc];
                if (partition[x] == partition[y])
                    cell_maps_arc_flags[cell_idx + partition_size][g.forward.original_arc[arc]] = true;
            }
            for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x+1]; ++arc) {
                int y = g.backward.head[arc];
                if (partition[x] == partition[y])
                    cell_maps_arc_flags[cell_idx][g.backward.original_arc[arc]] = true;
            }
        }
        saveFlags(cell_idx, cell_maps_arc_flags[cell_idx]);
        saveFlags(cell_idx + partition_size, cell_maps_arc_flags[cell_idx + partition_size]);
        cell_maps_arc_flags[cell_idx].clear();
        cell_maps_arc_flags[cell_idx + partition_size].clear();
        sync_out.println("Finished cell ", cell_idx);
    };

    for (int i = start; i < end; i++) {
        pool.push_task(precomputeCell, i);
    }

    while (pool.get_tasks_total() > 0) {
        sync_out.println(pool.get_tasks_total(),
            " tasks total, ",
            pool.get_tasks_running(),
            " tasks running, ",
            pool.get_tasks_queued(),
            " tasks queued.");
        this_thread::sleep_for(std::chrono::milliseconds(5000));
    }
    pool.wait_for_tasks();
    precomputed = true;
}

void ArcFlags::mergeFlags(string folder) {
    unordered_map<long long, boost::dynamic_bitset<>> labels;
    for (int i = 0; i < 2*partition_size; i++) {
        cout << "merge cell " << i << endl;
        string file_name = "../flags/" + g.name + "_" + to_string(partition_size) + "_" + to_string(i);
        ifstream file(file_name);
        string line;
        while (getline(file, line)) {
            long long arc = stoll(line);
            if (labels.find(arc) == labels.end()) labels[arc].resize(2*partition_size, 0);
            labels[arc][i] = 1;
        }
    }
    cout << "Start compressing" << endl;
    compress(labels);
    cout << "Start exporting" << endl;
    exportFlags(folder);
    cout << "Finished exporting" << endl;
    precomputed = true;
}

void ArcFlags::exportFlags(string folder) {
    ofstream file("../flags/" + g.name + "_" + to_string(partition_size) + ".csv");
    file << "edge_id,key\n";
    for (auto& [edge, key] : preprocessing_label_hash) {
        file << edge << "," << key << "\n";
    }
    file.close();

    vector<byte> bytes;
    for (auto& [_key, label] : preprocessing_labels) {
        long long key = _key;
        boost::dynamic_bitset<> flag = label;
        bitset<64> key_bits(key_bits);
        vector<byte> keyBytes;
        for (int i = 0; i < sizeof(size_t); i++) {
            byte nextByte{ static_cast<int>(key & ((1 << 8) - 1)) };
            keyBytes.push_back(nextByte);
            key >>= 8;
        }
        for (int i = 0; i < keyBytes.size(); i++)
            bytes.push_back(keyBytes[keyBytes.size() - 1 - i]);

        vector<byte> flagBytes;
        for (int i = 0; i < 2*partition_size; i += 8) {
            boost::dynamic_bitset<> copy = flag;
            copy &= boost::dynamic_bitset<>(2*partition_size, 255);
            byte nextByte{ static_cast<int>(copy.to_ulong()) };
            flagBytes.push_back(nextByte);
            flag >>= 8;
        }
        for (int i = 0; i < flagBytes.size(); i++)
            bytes.push_back(flagBytes[flagBytes.size() - 1 - i]);
    }
    ofstream outfile("../flags/" + g.name + "_" + to_string(partition_size) + ".bin", ios::out | ios::binary);
    outfile.write((const char*)&bytes[0], bytes.size());
}

void ArcFlags::importFlags(string edges_path, string flags_path) {
    ifstream input(flags_path, ios::binary);
    vector<unsigned char> buffer(istreambuf_iterator<char>(input), {});
    int rowByteSize = sizeof(size_t) + ceil(2*partition_size / 8.0);

    bitset<64> keyBits;
    boost::dynamic_bitset<> flag(2*partition_size);
    size_t key;
    for (int i = 0; i < buffer.size(); i++) {
        int pos = i % rowByteSize;
        if (pos == 0) {
            keyBits = bitset<64>{ 0 };
            flag = boost::dynamic_bitset<>(2*partition_size, 0);
        }
        if (pos < 8) { // read key
            bitset<64> mask{ buffer[i] };
            keyBits |= mask;
            if (pos == 7) {
                key = keyBits.to_ulong();
            } else keyBits <<= 8;
        } else { // read flag
            flag |= boost::dynamic_bitset<>(2*partition_size, buffer[i]);
            if (pos < rowByteSize - 1) {
                flag <<= 8;
            } else labels[key] = flag;
        }
    }
    ifstream file(edges_path);
    string line;
    getline(file, line); // skip first line
    vector<string> lines;
    while (getline(file, line)) {
        vector<string> csv = split(line, ",", false);
        int edge_id = stoi(csv[0]);
        size_t key = stoul(csv[1]);
        label_hashes[edge_id] = key;
    }
    precomputed = true;
}

void ArcFlags::compress(unordered_map<long long, boost::dynamic_bitset<>>& labels) {
    hash<boost::dynamic_bitset<>> hash_f;
    for (auto& [edge, label] : labels) {
        size_t key = hash_f(label);
        if (preprocessing_labels.find(key) != preprocessing_labels.end()) {
            if (label != preprocessing_labels[key])
                cout << "Hash collision!" << endl;
        } else preprocessing_labels[key] = label;
        preprocessing_label_hash[edge] = key;
    }
}