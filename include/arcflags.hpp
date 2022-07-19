#pragma once
#include <bits/stdc++.h>
#include <thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace std;

struct Graph {
private:
    int count;
public:
    Graph(int _count, string _name=""): count{_count}, name{_name} {}
    struct Side{
        vector<unsigned>first_out;
        vector<unsigned>head;
        vector<unsigned>weight;
        vector<unsigned>original_arc;
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
    void compress(unordered_map<int, boost::dynamic_bitset<>>& labels);
    unordered_map<size_t, boost::dynamic_bitset<>> preprocessing_labels;
    unordered_map<int, size_t> preprocessing_label_hash;
    bool precomputed = false;
};

void ArcFlags::precompute(int start, int end) {
    thread_pool pool; // amount of threads
    synced_stream sync_out;
    unordered_map<int, bool> cell_maps_arc_flags[2 * partition_size];

    using QueueElement = pair<int, double>;
    auto cmp = [](QueueElement &a, QueueElement &b) {
        return a.second > b.second;
    };
    auto markEdgesOnSptTo = [&](int src, int cell_idx) {
        vector<double> d(g.node_count(), numeric_limits<double>::max());
        vector<bool> down(g.node_count(), false);
        vector<int> p(g.node_count(), -1);
        vector<bool> visited(g.node_count(), false);
        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({src, 0});
        d[src] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v])
                continue;
            visited[v] = true;
            if (p[v] != -1) {
                cell_maps_arc_flags[cell_idx][p[v]] = true; 
            }
            if (v == src || !down[v]){
                for (int arc = g.backward_up.first_out[v]; arc < g.backward_up.first_out[v+1]; ++arc) {
                    int u = g.backward_up.head[arc];
                    if (visited[u]) continue;
                    if (d[v] + g.backward_up.weight[arc] < d[u]) {
                        d[u] = d[v] + g.backward_up.weight[arc];
                        p[u] = g.backward_up.original_arc[arc];
                        down[u] = false;
                        q.push({u, d[u]});
                    }
                }
            }
            for (int arc = g.backward_down.first_out[v]; arc < g.backward_down.first_out[v+1]; ++arc) {
                int u = g.backward_down.head[arc];
                if (visited[u]) continue;
                if (d[v] + g.backward_down.weight[arc] < d[u]) {
                    d[u] = d[v] + g.backward_down.weight[arc];
                    p[u] = g.backward_down.original_arc[arc];
                    down[u] = true;
                    q.push({u, d[u]});
                }
            }
        }
    };

    auto markEdgesOnSptFrom = [&](int src, int cell_idx) {
        vector<double> d(g.node_count(), numeric_limits<double>::max());
        vector<int> p(g.node_count(), -1);
        vector<bool> down(g.node_count(), false);
        vector<bool> visited(g.node_count(), false);
        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({src, 0});
        d[src] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v])
                continue;
            visited[v] = true;
            if (p[v] != -1) {
                cell_maps_arc_flags[cell_idx + partition_size][down[v]? -p[v] : p[v]] = true; 
            }
            if (v == src || !down[v]){
                for (int arc = g.forward_up.first_out[v]; arc < g.forward_up.first_out[v+1]; ++arc) {
                    int u = g.forward_up.head[arc];
                    if (visited[u]) continue;
                    if (d[v] + g.forward_up.weight[arc] < d[u]) {
                        d[u] = d[v] + g.forward_up.weight[arc];
                        p[u] = g.forward_up.original_arc[arc];
                        down[u] = false;
                        q.push({u, d[u]});
                    }
                }
            }
            for (int arc = g.forward_down.first_out[v]; arc < g.forward_down.first_out[v+1]; ++arc) {
                int u = g.forward_down.head[arc];
                if (visited[u]) continue;
                if (d[v] + g.forward_down.weight[arc] < d[u]) {
                    d[u] = d[v] + g.forward_down.weight[arc];
                    p[u] = g.forward_down.original_arc[arc];
                    down[u] = true;
                    q.push({u, d[u]});
                }
            }
        }
    };

    auto saveFlags = [&](int cell_idx, const unordered_map<int, bool> &map) {
        ofstream out("../flags/" + g.name + "_" + to_string(partition_size) + "_" + to_string(cell_idx));
        for (auto [e, _] : map)
            if (_) out << e << endl;
        out.close();
    };

    auto precomputeCell = [&](int cell_idx) {
        sync_out.println("Start computation for cell ", cell_idx);
        for (int x = 0; x < g.node_count(); ++x) {
            cout << "cell: " << cell_idx << " x: " << x << "/" << g.node_count()-1 << endl;
            if (g.boundary_nodes[x]){
                markEdgesOnSptTo(x, cell_idx);
                markEdgesOnSptFrom(x, cell_idx);
            }
            for (int arc = g.forward.first_out[x]; arc < g.forward.first_out[x+1]; ++arc) {
                cell_maps_arc_flags[cell_idx + partition_size][arc] = true;
            }
            for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x+1]; ++arc) {
                cell_maps_arc_flags[cell_idx][arc] = true;
            }
        }
        saveFlags(cell_idx, cell_maps_arc_flags[cell_idx]);
        saveFlags(cell_idx + partition_size, cell_maps_arc_flags[cell_idx + partition_size]);
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
    unordered_map<int, boost::dynamic_bitset<>> labels;
    for (int i = 0; i < 2*partition_size; i++) {
        cout << "merge cell " << i << endl;
        string file_name = "../flags/" + g.name + "_" + to_string(partition_size) + "_" + to_string(i);
        ifstream file(file_name);
        string line;
        while (getline(file, line)) {
            int arc = stoi(line);
            if (labels.find(arc) == labels.end()) labels[arc].resize(2*partition_size, 0);
            labels[arc][i] = 1;
        }
    }
    cout << "Start compressing" << endl;
    compress(labels);
    cout << "Start exporting" << endl;
    // exportFlags(folder);
    cout << "Finished exporting" << endl;
    precomputed = true;
}

void ArcFlags::compress(unordered_map<int, boost::dynamic_bitset<>>& labels) {
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