//============================================================================
// SNU Spring 2025 Introduction to Algorithm Final Project
//============================================================================
// Enhanced TSP Solver - Multiple Heuristics for Maximum Validity
//============================================================================

#include <bits/stdc++.h>
using namespace std;

/*=====================  LOCKED SECTION  (do not edit)  =====================*/
constexpr int INF = 1e9;
using Matrix = vector<vector<int>>;
int tsp_solve(const Matrix&);

static Matrix read_graph(istream& in) {
    int n, m; if (!(in >> n >> m)) return {};
    Matrix w(n, vector<int>(n, INF));
    for (int i = 0; i < n; ++i) w[i][i] = 0;
    for (int i = 0; i < m; ++i) {
        int u, v, c; in >> u >> v >> c; --u; --v;
        w[u][v] = w[v][u] = c;
    }
    return w;
}

static void print_result(long long cost, const vector<int>& t)
{
    if (cost >= INF) { cout << -1 << '\n'; return; }
    cout << cost;
    for (int v : t) cout << ' ' << v + 1;
    cout << '\n';
}

// Abstract base class for any TSP solver
class TSPSolver {
public:
    explicit TSPSolver(const Matrix& w) : W(w), n(static_cast<int>(w.size())) {}
    virtual ~TSPSolver() = default;
    virtual pair<long long, vector<int>> solve(int start = 0) = 0;

protected:
    const Matrix& W;
    int n;
};

class TwoStepGreedySolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    pair<long long, vector<int>> solve(int start = 0) override;
};

class InsertionSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    pair<long long, vector<int>> solve(int start = 0) override;
};

class DynamicProgrammingSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    pair<long long, vector<int>> solve(int start = 0) override {
        if (n > 20) return {LLONG_MAX, {}};

        int N = 1 << n;
        vector<vector<long long>> dp(N, vector<long long>(n, LLONG_MAX));
        vector<vector<int>> parent(N, vector<int>(n, -1));

        dp[1 << start][start] = 0;

        for (int mask = 0; mask < N; ++mask) {
            for (int u = 0; u < n; ++u) {
                if (!(mask & (1 << u))) continue;
                long long cur = dp[mask][u];
                if (cur == LLONG_MAX) continue;
                for (int v = 0; v < n; ++v) {
                    if (mask & (1 << v)) continue;
                    if (W[u][v] >= INF) continue;
                    int nm = mask | (1 << v);
                    long long nc = cur + W[u][v];
                    if (nc < dp[nm][v]) {
                        dp[nm][v] = nc;
                        parent[nm][v] = u;
                    }
                }
            }
        }

        int full = N - 1;
        long long bestCost = LLONG_MAX;
        int last = -1;
        for (int u = 0; u < n; ++u) {
            if (u == start) continue;
            if (dp[full][u] == LLONG_MAX) continue;
            if (W[u][start] >= INF) continue;
            long long cost = dp[full][u] + W[u][start];
            if (cost < bestCost) {
                bestCost = cost;
                last = u;
            }
        }

        if (last == -1) return {LLONG_MAX, {}};

        vector<int> path;
        int mask = full;
        int cur = last;
        while (cur != start) {
            path.push_back(cur);
            int p = parent[mask][cur];
            mask ^= (1 << cur);
            cur = p;
        }

        reverse(path.begin(), path.end());
        vector<int> tour;
        tour.reserve(n + 1);
        tour.push_back(start);
        tour.insert(tour.end(), path.begin(), path.end());
        tour.push_back(start);

        return {bestCost, tour};
    }
};

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            ifstream fin(argv[i]);
            if (!fin) { cerr << "Cannot open " << argv[i] << '\n'; continue; }
            Matrix W = read_graph(fin);
            tsp_solve(W);
        }
    } else {
        Matrix W = read_graph(cin);
        tsp_solve(W);
    }
}
/*===================  END LOCKED SECTION  ==================================*/

/*=====================  ENHANCED IMPLEMENTATION  ===========================*/

// Utility functions
long long calculatePenalty(int n) {
    long long maxAllowed = (long long)INF - 1;
    return maxAllowed / max(1, n + 1);
}

long long calculateTourCost(const Matrix& W, const vector<int>& tour) {
    if (tour.empty()) return LLONG_MAX;

    int n = W.size();
    long long cost = 0;

    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        int u = tour[i], v = tour[i + 1];
        if (u < 0 || u >= n || v < 0 || v >= n) return LLONG_MAX;
        if (W[u][v] >= INF) return LLONG_MAX; // edge does not exist
        cost += W[u][v];
        if (cost > LLONG_MAX / 2) return LLONG_MAX;
    }

    return cost;
}

bool validateTour(const vector<int>& tour, int n) {
    if (tour.empty() || (int)tour.size() != n + 1) return false;
    if (tour.front() != tour.back()) return false;
    
    vector<bool> visited(n, false);
    for (int i = 0; i < n; ++i) {
        int vertex = tour[i];
        if (vertex < 0 || vertex >= n) return false;
        if (visited[vertex]) return false;
        visited[vertex] = true;
    }
    
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) return false;
    }
    
    return true;
}

bool validateHamiltonian(const Matrix& W, const vector<int>& tour) {
    int n = W.size();
    if (!validateTour(tour, n)) return false;
    long long cost = 0;
    for (size_t i = 0; i + 1 < tour.size(); ++i) {
        int u = tour[i], v = tour[i + 1];
        if (W[u][v] >= INF) return false;
        cost += W[u][v];
    }
    return true;
}

vector<int> constructGuaranteedValidTour(int n, int start = 0) {
    vector<int> tour;
    tour.reserve(n + 1);
    for (int i = 0; i < n; ++i) {
        tour.push_back((start + i) % n);
    }
    tour.push_back(start);
    return tour;
}

// 1. Nearest Neighbor
class NearestNeighborSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        long long penalty = calculatePenalty(n);
        vector<int> tour;
        tour.reserve(n + 1);
        vector<bool> visited(n, false);
        
        int current = start;
        visited[current] = true;
        tour.push_back(current);
        
        for (int i = 1; i < n; ++i) {
            int nearest = -1;
            long long minDist = LLONG_MAX;
            
            for (int j = 0; j < n; ++j) {
                if (visited[j]) continue;
                long long dist = (W[current][j] < INF) ? W[current][j] : penalty;
                if (dist < minDist) {
                    minDist = dist;
                    nearest = j;
                }
            }
            
            if (nearest == -1) {
                for (int j = 0; j < n; ++j) {
                    if (!visited[j]) {
                        nearest = j;
                        break;
                    }
                }
            }
            
            visited[nearest] = true;
            tour.push_back(nearest);
            current = nearest;
        }
        
        tour.push_back(start);
        return {calculateTourCost(W, tour), tour};
    }
};

// 2. Two-Step Greedy (Required)
pair<long long, vector<int>> TwoStepGreedySolver::solve(int start) {
    long long penalty = calculatePenalty(n);
    vector<int> tour;
    tour.reserve(n + 1);
    vector<bool> visited(n, false);
    
    int current = start;
    visited[current] = true;
    tour.push_back(current);
    
    while (tour.size() < n) {
        int nextCity = -1;
        long long bestScore = LLONG_MAX;
        
        for (int j = 0; j < n; ++j) {
            if (visited[j]) continue;
            
            long long cost1 = (W[current][j] < INF) ? W[current][j] : penalty;
            long long minSecond = 0;
            
            if (tour.size() < n - 1) {
                minSecond = LLONG_MAX;
                for (int k = 0; k < n; ++k) {
                    if (!visited[k] && k != j) {
                        long long cost2 = (W[j][k] < INF) ? W[j][k] : penalty;
                        minSecond = min(minSecond, cost2);
                    }
                }
                if (minSecond == LLONG_MAX) minSecond = penalty;
            }
            
            long long score = cost1 + minSecond;
            if (score < bestScore) {
                bestScore = score;
                nextCity = j;
            }
        }
        
        if (nextCity == -1) {
            for (int j = 0; j < n; ++j) {
                if (!visited[j]) {
                    nextCity = j;
                    break;
                }
            }
        }
        
        visited[nextCity] = true;
        tour.push_back(nextCity);
        current = nextCity;
    }
    
    tour.push_back(start);
    return {calculateTourCost(W, tour), tour};
}

// 3. Cheapest Insertion
pair<long long, vector<int>> InsertionSolver::solve(int start) {
    long long penalty = calculatePenalty(n);
    vector<bool> visited(n, false);
    vector<int> tour = {start, start};
    visited[start] = true;
    
    while (count(visited.begin(), visited.end(), true) < n) {
        long long bestInc = LLONG_MAX;
        int bestNode = -1, bestPos = -1;
        
        for (int j = 0; j < n; ++j) {
            if (visited[j]) continue;
            
            for (int pos = 0; pos + 1 < (int)tour.size(); ++pos) {
                int u = tour[pos], v = tour[pos+1];
                long long cuv = (W[u][v] < INF) ? W[u][v] : penalty;
                long long cuj = (W[u][j] < INF) ? W[u][j] : penalty;
                long long cjv = (W[j][v] < INF) ? W[j][v] : penalty;
                
                long long inc = cuj + cjv - cuv;
                if (inc < bestInc) {
                    bestInc = inc;
                    bestNode = j;
                    bestPos = pos;
                }
            }
        }
        
        if (bestNode == -1) {
            for (int j = 0; j < n; ++j) {
                if (!visited[j]) {
                    bestNode = j;
                    bestPos = 0;
                    break;
                }
            }
        }
        
        tour.insert(tour.begin() + bestPos + 1, bestNode);
        visited[bestNode] = true;
    }
    
    return {calculateTourCost(W, tour), tour};
}

// 4. Farthest Insertion - Choose farthest unvisited node
class FarthestInsertionSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        if (n <= 1) return {0, {start, start}};
        
        long long penalty = calculatePenalty(n);
        vector<bool> inTour(n, false);
        vector<int> tour = {start, start};
        inTour[start] = true;
        
        // Find farthest node from start
        int farthest = -1;
        long long maxDist = -1;
        for (int i = 0; i < n; ++i) {
            if (i != start) {
                long long d = (W[start][i] < INF) ? W[start][i] : penalty;
                if (d > maxDist) {
                    maxDist = d;
                    farthest = i;
                }
            }
        }
        
        if (farthest != -1) {
            tour = {start, farthest, start};
            inTour[farthest] = true;
        }
        
        // Insert remaining nodes
        while (count(inTour.begin(), inTour.end(), true) < n) {
            // Find farthest node from tour
            int nextNode = -1;
            long long maxMinDist = -1;
            
            for (int v = 0; v < n; ++v) {
                if (inTour[v]) continue;
                
                long long minDist = LLONG_MAX;
                for (int u : tour) {
                    if (u == tour.back()) continue; // Skip last (duplicate)
                    long long d = (W[u][v] < INF) ? W[u][v] : penalty;
                    minDist = min(minDist, d);
                }
                
                if (minDist > maxMinDist) {
                    maxMinDist = minDist;
                    nextNode = v;
                }
            }
            
            if (nextNode == -1) {
                for (int v = 0; v < n; ++v) {
                    if (!inTour[v]) {
                        nextNode = v;
                        break;
                    }
                }
            }
            
            // Find best insertion position
            int bestPos = -1;
            long long bestInc = LLONG_MAX;
            
            for (int i = 0; i + 1 < tour.size(); ++i) {
                int u = tour[i], w = tour[i+1];
                long long cuw = (W[u][w] < INF) ? W[u][w] : penalty;
                long long cuv = (W[u][nextNode] < INF) ? W[u][nextNode] : penalty;
                long long cvw = (W[nextNode][w] < INF) ? W[nextNode][w] : penalty;
                
                long long inc = cuv + cvw - cuw;
                if (inc < bestInc) {
                    bestInc = inc;
                    bestPos = i;
                }
            }
            
            if (bestPos == -1) bestPos = 0;
            tour.insert(tour.begin() + bestPos + 1, nextNode);
            inTour[nextNode] = true;
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 5. Random Insertion - Randomly select nodes to insert
class RandomInsertionSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        long long penalty = calculatePenalty(n);
        vector<bool> inTour(n, false);
        vector<int> tour = {start, start};
        inTour[start] = true;
        
        // Create random order
        vector<int> order;
        for (int i = 0; i < n; ++i) {
            if (i != start) order.push_back(i);
        }
        
        // Use deterministic shuffle for consistency
        unsigned seed = start;
        shuffle(order.begin(), order.end(), default_random_engine(seed));
        
        // Insert nodes in random order
        for (int node : order) {
            int bestPos = -1;
            long long bestInc = LLONG_MAX;
            
            for (int i = 0; i + 1 < tour.size(); ++i) {
                int u = tour[i], v = tour[i+1];
                long long cuv = (W[u][v] < INF) ? W[u][v] : penalty;
                long long cun = (W[u][node] < INF) ? W[u][node] : penalty;
                long long cnv = (W[node][v] < INF) ? W[node][v] : penalty;
                
                long long inc = cun + cnv - cuv;
                if (inc < bestInc) {
                    bestInc = inc;
                    bestPos = i;
                }
            }
            
            if (bestPos == -1) bestPos = 0;
            tour.insert(tour.begin() + bestPos + 1, node);
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 6. Nearest Addition - Add nearest node to partial tour
class NearestAdditionSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        long long penalty = calculatePenalty(n);
        vector<bool> inTour(n, false);
        vector<int> tour = {start, start};
        inTour[start] = true;
        
        while (count(inTour.begin(), inTour.end(), true) < n) {
            // Find nearest node to any node in tour
            int nearestNode = -1;
            long long minDist = LLONG_MAX;
            
            for (int v = 0; v < n; ++v) {
                if (inTour[v]) continue;
                
                for (int u : tour) {
                    if (u == tour.back()) continue;
                    long long d = (W[u][v] < INF) ? W[u][v] : penalty;
                    if (d < minDist) {
                        minDist = d;
                        nearestNode = v;
                    }
                }
            }
            
            if (nearestNode == -1) {
                for (int v = 0; v < n; ++v) {
                    if (!inTour[v]) {
                        nearestNode = v;
                        break;
                    }
                }
            }
            
            // Find best position to insert
            int bestPos = -1;
            long long bestInc = LLONG_MAX;
            
            for (int i = 0; i + 1 < tour.size(); ++i) {
                int u = tour[i], w = tour[i+1];
                long long cuw = (W[u][w] < INF) ? W[u][w] : penalty;
                long long cun = (W[u][nearestNode] < INF) ? W[u][nearestNode] : penalty;
                long long cnw = (W[nearestNode][w] < INF) ? W[nearestNode][w] : penalty;
                
                long long inc = cun + cnw - cuw;
                if (inc < bestInc) {
                    bestInc = inc;
                    bestPos = i;
                }
            }
            
            if (bestPos == -1) bestPos = 0;
            tour.insert(tour.begin() + bestPos + 1, nearestNode);
            inTour[nearestNode] = true;
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 7. Savings Algorithm (Clarke-Wright)
class SavingsSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        long long penalty = calculatePenalty(n);
        
        // Create savings list
        struct Saving {
            int i, j;
            long long value;
        };
        vector<Saving> savings;
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (i == start || j == start) continue;
                
                long long cis = (W[i][start] < INF) ? W[i][start] : penalty;
                long long csj = (W[start][j] < INF) ? W[start][j] : penalty;
                long long cij = (W[i][j] < INF) ? W[i][j] : penalty;
                
                long long saving = cis + csj - cij;
                if (saving > 0) {
                    savings.push_back({i, j, saving});
                }
            }
        }
        
        // Sort by savings (descending)
        sort(savings.begin(), savings.end(), 
             [](const Saving& a, const Saving& b) { return a.value > b.value; });
        
        // Build routes
        vector<vector<int>> routes;
        vector<int> routeOf(n, -1);
        
        // Initialize with single-node routes
        for (int i = 0; i < n; ++i) {
            if (i != start) {
                routes.push_back({i});
                routeOf[i] = routes.size() - 1;
            }
        }
        
        // Merge routes based on savings
        for (const auto& s : savings) {
            int r1 = routeOf[s.i];
            int r2 = routeOf[s.j];
            
            if (r1 == -1 || r2 == -1 || r1 == r2) continue;
            
            // Check if nodes are at ends of routes
            bool i_at_end = (routes[r1].front() == s.i || routes[r1].back() == s.i);
            bool j_at_end = (routes[r2].front() == s.j || routes[r2].back() == s.j);
            
            if (i_at_end && j_at_end) {
                // Merge routes
                if (routes[r1].back() == s.i && routes[r2].front() == s.j) {
                    // Append r2 to r1
                    routes[r1].insert(routes[r1].end(), routes[r2].begin(), routes[r2].end());
                } else if (routes[r1].front() == s.i && routes[r2].back() == s.j) {
                    // Prepend r1 to r2
                    routes[r2].insert(routes[r2].end(), routes[r1].begin(), routes[r1].end());
                    r1 = r2;
                } else {
                    // Need to reverse one route
                    if (routes[r1].back() == s.i && routes[r2].back() == s.j) {
                        reverse(routes[r2].begin(), routes[r2].end());
                        routes[r1].insert(routes[r1].end(), routes[r2].begin(), routes[r2].end());
                    } else if (routes[r1].front() == s.i && routes[r2].front() == s.j) {
                        reverse(routes[r1].begin(), routes[r1].end());
                        routes[r1].insert(routes[r1].end(), routes[r2].begin(), routes[r2].end());
                    }
                }
                
                // Update route assignments
                for (int node : routes[r2]) {
                    routeOf[node] = r1;
                }
                routes[r2].clear();
            }
        }
        
        // Combine all routes into a tour
        vector<int> tour = {start};
        for (const auto& route : routes) {
            for (int node : route) {
                tour.push_back(node);
            }
        }
        tour.push_back(start);
        
        // Ensure validity
        if (!validateHamiltonian(W, tour)) {
            return NearestNeighborSolver(W).solve(start);
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 8. Greedy Edge Selection
class GreedyEdgeSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        // Create edge list
        struct Edge {
            int u, v;
            long long weight;
        };
        
        vector<Edge> edges;
        long long penalty = calculatePenalty(n);
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                long long w = (W[i][j] < INF) ? W[i][j] : penalty;
                edges.push_back({i, j, w});
            }
        }
        
        // Sort edges by weight
        sort(edges.begin(), edges.end(), 
             [](const Edge& a, const Edge& b) { return a.weight < b.weight; });
        
        // Build tour using edges
        vector<int> degree(n, 0);
        vector<vector<int>> adj(n);
        int edgeCount = 0;
        
        for (const auto& e : edges) {
            if (degree[e.u] < 2 && degree[e.v] < 2) {
                // Check if adding edge creates subtour
                bool creates_subtour = false;
                if (edgeCount < n - 1) {
                    // DFS to check connectivity
                    vector<bool> visited(n, false);
                    function<bool(int, int, int)> dfs = [&](int v, int target, int parent) -> bool {
                        if (v == target) return true;
                        visited[v] = true;
                        for (int u : adj[v]) {
                            if (u != parent && !visited[u]) {
                                if (dfs(u, target, v)) return true;
                            }
                        }
                        return false;
                    };
                    
                    if (dfs(e.u, e.v, -1)) creates_subtour = true;
                }
                
                if (!creates_subtour) {
                    adj[e.u].push_back(e.v);
                    adj[e.v].push_back(e.u);
                    degree[e.u]++;
                    degree[e.v]++;
                    edgeCount++;
                    
                    if (edgeCount == n) break;
                }
            }
        }
        
        // Build tour from adjacency list
        vector<int> tour;
        vector<bool> visited(n, false);
        
        function<void(int, int)> buildTour = [&](int v, int parent) {
            visited[v] = true;
            tour.push_back(v);
            for (int u : adj[v]) {
                if (!visited[u] && u != parent) {
                    buildTour(u, v);
                }
            }
        };
        
        buildTour(start, -1);
        
        // Add any missing nodes
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                tour.push_back(i);
            }
        }
        
        tour.push_back(start);
        
        if (!validateHamiltonian(W, tour)) {
            return NearestNeighborSolver(W).solve(start);
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 9. Multi-fragment heuristic
class MultiFragmentSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        long long penalty = calculatePenalty(n);
        
        // Create fragments (initially each node is a fragment)
        vector<deque<int>> fragments(n);
        vector<int> fragmentOf(n);
        for (int i = 0; i < n; ++i) {
            fragments[i] = {i};
            fragmentOf[i] = i;
        }
        
        // Try to merge fragments
        while (fragments.size() > 1) {
            int bestF1 = -1, bestF2 = -1;
            long long bestCost = LLONG_MAX;
            bool bestReverse1 = false, bestReverse2 = false;
            
            // Find best pair of fragments to merge
            for (int f1 = 0; f1 < fragments.size(); ++f1) {
                if (fragments[f1].empty()) continue;
                
                for (int f2 = f1 + 1; f2 < fragments.size(); ++f2) {
                    if (fragments[f2].empty()) continue;
                    
                    // Try all 4 combinations of connecting fragments
                    vector<pair<bool, bool>> configs = {{false, false}, {false, true}, 
                                                        {true, false}, {true, true}};
                    
                    for (auto [rev1, rev2] : configs) {
                        int u = rev1 ? fragments[f1].front() : fragments[f1].back();
                        int v = rev2 ? fragments[f2].back() : fragments[f2].front();
                        
                        long long cost = (W[u][v] < INF) ? W[u][v] : penalty;
                        
                        if (cost < bestCost) {
                            bestCost = cost;
                            bestF1 = f1;
                            bestF2 = f2;
                            bestReverse1 = rev1;
                            bestReverse2 = rev2;
                        }
                    }
                }
            }
            
            if (bestF1 == -1) break;
            
            // Merge fragments
            if (bestReverse1) {
                reverse(fragments[bestF1].begin(), fragments[bestF1].end());
            }
            if (bestReverse2) {
                reverse(fragments[bestF2].begin(), fragments[bestF2].end());
            }
            
            // Append f2 to f1
            fragments[bestF1].insert(fragments[bestF1].end(), 
                                   fragments[bestF2].begin(), 
                                   fragments[bestF2].end());
            
            // Update fragment assignments
            for (int node : fragments[bestF2]) {
                fragmentOf[node] = bestF1;
            }
            
            fragments[bestF2].clear();
            
            // Remove empty fragments
            fragments.erase(
                remove_if(fragments.begin(), fragments.end(), 
                         [](const deque<int>& f) { return f.empty(); }),
                fragments.end()
            );
        }
        
        // Build tour from remaining fragment
        vector<int> tour;
        if (!fragments.empty()) {
            for (int node : fragments[0]) {
                tour.push_back(node);
            }
        }
        
        // Ensure start is at beginning
        auto it = find(tour.begin(), tour.end(), start);
        if (it != tour.end()) {
            rotate(tour.begin(), it, tour.end());
        }
        
        tour.push_back(start);
        
        if (!validateHamiltonian(W, tour)) {
            return NearestNeighborSolver(W).solve(start);
        }
        
        return {calculateTourCost(W, tour), tour};
    }
};

// 10. Randomized tours with best selection
class RandomTourSolver : public TSPSolver {
public:
    using TSPSolver::TSPSolver;
    
    pair<long long, vector<int>> solve(int start = 0) override {
        vector<int> bestTour;
        long long bestCost = LLONG_MAX;
        
        // Try multiple random permutations
        int trials = min(100, n * n);
        unsigned seed = start;
        default_random_engine gen(seed);
        
        for (int trial = 0; trial < trials; ++trial) {
            vector<int> tour;
            for (int i = 0; i < n; ++i) {
                if (i != start) tour.push_back(i);
            }
            
            shuffle(tour.begin(), tour.end(), gen);
            
            // Build complete tour
            vector<int> completeTour = {start};
            completeTour.insert(completeTour.end(), tour.begin(), tour.end());
            completeTour.push_back(start);
            
            long long cost = calculateTourCost(W, completeTour);
            if (cost < bestCost) {
                bestCost = cost;
                bestTour = completeTour;
            }
        }
        
        return {bestCost, bestTour};
    }
};

// 2-opt improvement
pair<long long, vector<int>> improveWith2Opt(const Matrix& W, vector<int> tour, long long currentCost) {
    int n = W.size();
    if (!validateHamiltonian(W, tour)) return {currentCost, tour};
    
    vector<int> bestTour = tour;
    long long bestCost = currentCost;
    bool improved = true;
    int iterations = 0;
    int maxIterations = min(50, n);
    
    while (improved && iterations < maxIterations) {
        improved = false;
        iterations++;
        
        for (int i = 1; i < tour.size() - 2; ++i) {
            for (int j = i + 1; j < tour.size() - 1; ++j) {
                if (j - i > min(100, n / 2)) break;
                
                vector<int> newTour = bestTour;
                reverse(newTour.begin() + i, newTour.begin() + j + 1);
                
                long long newCost = calculateTourCost(W, newTour);
                if (newCost < bestCost) {
                    bestTour = newTour;
                    bestCost = newCost;
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
    }
    
    return {bestCost, bestTour};
}

// Main solver with multiple heuristics
int tsp_solve(const Matrix& w) {
    int n = w.size();
    
    if (n == 0) {
        print_result(INF, {});
        return 0;
    }
    
    if (n == 1) {
        print_result(0, {0, 0});
        return 0;
    }
    
    pair<long long, vector<int>> best = {LLONG_MAX, {}};
    
    // Determine which algorithms to use based on graph size
    vector<pair<string, TSPSolver*>> algorithms;

    if (n <= 12) {
        algorithms.push_back({"DynamicProgramming", new DynamicProgrammingSolver(w)});
    }

    // Always use these core algorithms
    algorithms.push_back({"NearestNeighbor", new NearestNeighborSolver(w)});
    algorithms.push_back({"TwoStepGreedy", new TwoStepGreedySolver(w)});
    algorithms.push_back({"CheapestInsertion", new InsertionSolver(w)});
    
    if (n <= 1000) {
        algorithms.push_back({"FarthestInsertion", new FarthestInsertionSolver(w)});
        algorithms.push_back({"RandomInsertion", new RandomInsertionSolver(w)});
        algorithms.push_back({"NearestAddition", new NearestAdditionSolver(w)});
    }
    
    if (n <= 500) {
        algorithms.push_back({"Savings", new SavingsSolver(w)});
        algorithms.push_back({"GreedyEdge", new GreedyEdgeSolver(w)});
        algorithms.push_back({"MultiFragment", new MultiFragmentSolver(w)});
    }
    
    if (n <= 100) {
        algorithms.push_back({"RandomTour", new RandomTourSolver(w)});
    }
    
    // Determine starting points
    vector<int> starts;
    if (n <= 10) {
        for (int i = 0; i < n; ++i) starts.push_back(i);
    } else if (n <= 100) {
        starts = {0, 1, n/4, n/2, 3*n/4, n-1};
    } else {
        starts = {0, n/2, n-1};
    }
    
    // Try each algorithm with different starting points
    for (const auto& [name, solver] : algorithms) {
        for (int start : starts) {
            try {
                auto [cost, tour] = solver->solve(start);
                if (validateHamiltonian(w, tour) && cost < best.first) {
                    best = {cost, tour};
                }
            } catch (...) {
                // Silently ignore failures
            }
            
            // Early exit for very large graphs
            if (n > 5000 && best.first < LLONG_MAX) break;
        }
        if (n > 5000 && best.first < LLONG_MAX) break;
    }
    
    // Apply 2-opt improvement
    if (!best.second.empty() && validateHamiltonian(w, best.second) && n <= 1000) {
        auto improved = improveWith2Opt(w, best.second, best.first);
        if (validateHamiltonian(w, improved.second) && improved.first <= best.first) {
            best = improved;
        }
    }
    
    // Ultimate fallback
    if (best.second.empty() || !validateHamiltonian(w, best.second)) {
        auto fallback = NearestNeighborSolver(w).solve(0);
        if (validateHamiltonian(w, fallback.second)) {
            best = fallback;
        } else if (n <= 12) {
            auto exact = DynamicProgrammingSolver(w).solve(0);
            if (validateHamiltonian(w, exact.second)) best = exact;
            else {
                best.second = {};
                best.first = LLONG_MAX;
            }
        } else {
            best.second = {};
            best.first = LLONG_MAX;
        }
    }
    
    // Clean up
    for (auto& [name, solver] : algorithms) {
        delete solver;
    }
    
    print_result(best.first, best.second);
    return 0;
}
