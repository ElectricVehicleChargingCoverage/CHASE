#pragma once
#include <bits/stdc++.h>
#include <thread_pool.hpp>

#define PARTITION_SIZE 500

using namespace std;
using Label = bitset<2*PARTITION_SIZE>;

struct Graph {
private:
    int count;
public:
    Graph(int _count, string _name): count{_count}, name{_name} {}
    struct Side{
        vector<int>first_out;
        vector<int>head;
        vector<int>weight;
    };
    string name;
	Side forward, backward;
    vector<bool> boundary_nodes;
    int node_count() { return count; };
};

class ArcFlags {
    Graph& g;
public:
    ArcFlags(Graph& _g): g{_g} {};
    void precompute(int start=0, int end=PARTITION_SIZE);
};

void ArcFlags::precompute(int start, int end) {
    thread_pool pool(15); // amount of threads
    synced_stream sync_out;
    unordered_map<int, bool> cell_maps_arc_flags[2 * PARTITION_SIZE];

    using QueueElement = pair<int, double>;
    auto cmp = [](QueueElement &a, QueueElement &b) {
        return a.second > b.second;
    };
    auto markEdgesOnSptTo = [&](int src, int cell_idx) {
        vector<int> d(g.node_count(), numeric_limits<int>::max());
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
            if (p[v]) {
                cell_maps_arc_flags[cell_idx][p[v]] = true; 
            }
            for (int arc_index = g.backward.first_out[v]; arc_index < g.backward.first_out[v + 1]; ++arc_index) {
                int u = g.backward.head[arc_index];
                if (visited[u])
                    continue;
                if (d[v] + g.backward.weight[arc_index] < d[u]) {
                    d[u] = d[v] + g.backward.weight[arc_index];
                    p[u] = arc_index;
                    q.push({u, d[u]});
                }
            }
        }
    };

    auto markEdgesOnSptFrom = [&](int src, int cell_idx) {
        vector<int> d(g.node_count(), numeric_limits<int>::max());
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
            if (p[v]) {
                cell_maps_arc_flags[cell_idx + PARTITION_SIZE][p[v]] = true; 
            }
            for (int arc_index = g.forward.first_out[v]; arc_index < g.forward.first_out[v + 1]; ++arc_index) {
                int u = g.forward.head[arc_index];
                if (visited[u])
                    continue;
                if (d[v] + g.forward.weight[arc_index] < d[u]) {
                    d[u] = d[v] + g.forward.weight[arc_index];
                    p[u] = arc_index;
                    q.push({u, d[u]});
                }
            }
        }
    };

    auto saveFlags = [&](int cell_idx, const unordered_map<int, bool> &map) {
        ofstream out("../flags/" + g.name + "_" + to_string(PARTITION_SIZE) + "_" + to_string(cell_idx));
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
                cell_maps_arc_flags[cell_idx + PARTITION_SIZE][arc] = true;
            }
            for (int arc = g.backward.first_out[x]; arc < g.backward.first_out[x+1]; ++arc) {
                cell_maps_arc_flags[cell_idx][arc] = true;
            }
        }
        saveFlags(cell_idx, cell_maps_arc_flags[cell_idx]);
        saveFlags(cell_idx + PARTITION_SIZE, cell_maps_arc_flags[cell_idx + PARTITION_SIZE]);
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
}