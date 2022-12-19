#include <bits/stdc++.h>
//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <random>
//#include <cmath>

using namespace std;

struct Activity {
    vector<int> rows;
    int removed = -1; // уйдет в следующем слое
    int added = -1; // пришел в этом слое
};

struct Data {
    double sum;
    int parent;
    int depth;
};

int n, k;
vector<vector<bool>> g;
vector<vector<bool>> gSpan;
vector<Activity> activity_rows;
vector<pair<int, int>> graph;

vector<bool> encode(const vector<bool>& v) {
    vector<bool> res;

    for (int j = 0; j < n; j++) {
        res.push_back(false);
        for (int i = 0; i < k; i++) {
            res[j] = res[j] ^ (v[i] & g[i][j]);
        }
    }

    return res;
}

void sum_matrix_rows(int row_from, int row_to) {
    for (int j = 0; j < n; j++) {
        gSpan[row_to][j] = gSpan[row_to][j] ^ gSpan[row_from][j];
    }
}

void make_triangle(int curr_row) {
    for (int i = curr_row + 1; i < k; i++) {
        if (gSpan[i][curr_row]) {
            if (!gSpan[curr_row][curr_row]) {
                sum_matrix_rows(i, curr_row);
            }
            sum_matrix_rows(curr_row, i);
        }
    }
}

void make_unique_starts(int curr_column, bool used[]) {
    int latest_true_row = -1;
    for (int i = 0; i < k; i++) {
        if (!used[i] && gSpan[i][curr_column]) {
            if (latest_true_row == -1) {
                latest_true_row = i;
                continue;
            }

            for (int j = n - 1; j >= 0; j--) {
                if (gSpan[i][j]) {
                    sum_matrix_rows(latest_true_row, i);
                    break;
                }
                if (gSpan[latest_true_row][j]) {
                    sum_matrix_rows(i, latest_true_row);
                    latest_true_row = i;
                    break;
                }
            }

        }
    }

    if (latest_true_row != -1) {
        used[latest_true_row] = true;
    }
}

void make_unique_ends(int curr_column, bool used[]) {
    int latest_true_row = -1;
    for (int i = k - 1; i >= 0; i--) {
        if (!used[i] && gSpan[i][curr_column]) {
            if (latest_true_row == -1) {
                latest_true_row = i;
                used[i] = true;
                continue;
            }

            sum_matrix_rows(latest_true_row, i);
        }
    }
}

void construct_activities() {
    vector<pair<int, int>> activity_range;
    for (int i = 0; i < k; i++) {
        int j_start = 0;
        while (j_start < n && !gSpan[i][j_start]) {
            j_start++;
        }

        int j_end = n - 1;
        while (j_end >= 0 && !gSpan[i][j_end]) {
            j_end--;
        }

        activity_range.emplace_back(j_start, j_end);
    }

    Activity empty_activity;
    activity_rows.push_back(empty_activity);

    for (int j = 0; j < n; j++) {
        vector<int> active_rows;
        int removed = -1;
        int added = -1;

        for (int i = 0; i < k; i++) {
            if (activity_range[i].first <= j && j < activity_range[i].second) {
                active_rows.push_back(i);
            }
            if (activity_range[i].first == j) {
                added = i;
            }
            if (activity_range[i].second == j + 1) {
                removed = i;
            }
        }

        activity_rows.emplace_back(Activity{active_rows, removed, added});
    }
}

void construct_span_matrix() {
    for (int i = 0; i < k; i++) {
        gSpan.push_back(g[i]);
    }

    bool used[k];

    for (int i = 0; i < k; i++) {
        make_triangle(i);
        used[i] = false;
    }

    for (int j = n - 1; j >= 0; j--) {
        make_unique_ends(j, used);
    }

    for (int i = 0; i < k; i++) {
        used[i] = false;
    }

    for (int j = 0; j < n; j++) {
        make_unique_starts(j, used);
    }

    construct_activities();
}

vector<bool> get_mask(int num, int power) {
    vector<bool> res;

    while (num > 0) {
        res.push_back(num % 2 == 1);
        num /= 2;
    }

    while (res.size() < power) {
        res.push_back(false);
    }

    return res;
}

vector<bool> gen_mask(int layer_num, const vector<bool>& mask, bool added_val) {
    vector<bool> res;
    size_t must_be = activity_rows[layer_num + 1].rows.size();

    int i = 0;
    int j = 0;
    while (i < activity_rows[layer_num].rows.size() && j < activity_rows[layer_num + 1].rows.size()) {
        if (activity_rows[layer_num].rows[i] == activity_rows[layer_num + 1].rows[j]) {
            res.push_back(mask[i]);
            i++; j++;
        } else {
            if (activity_rows[layer_num].removed == activity_rows[layer_num].rows[i]) {
                if (activity_rows[layer_num + 1].added == activity_rows[layer_num + 1].rows[j]) {
                    res.push_back(added_val);
                    i++; j++;
                } else {
                    i++;
                }
            } else {
                if (activity_rows[layer_num + 1].added == activity_rows[layer_num + 1].rows[j]) {
                    res.push_back(added_val);
                    j++;
                }
            }
        }
    }

    while (res.size() < must_be && i < activity_rows[layer_num].rows.size()) {
        res.push_back(mask[i]);
        i++;
    }

    if (j < activity_rows[layer_num + 1].rows.size()) {
        res.push_back(added_val);
    }

    return res;
}

int parse_from_mask(const vector<bool>& mask, int power) {
    int res = 0;
    for (int i = 0; i < mask.size(); i++) {
        if (mask[i]) {
            res += (1 << i);
        }
    }

    return res;
}

bool mul_two_vectors(const vector<bool>& a, const vector<bool>& b) {
    bool res = false;

    for (int i = 0; i < a.size(); i++) {
        res = res ^ (a[i] & b[i]);
    }

    return res;
}

bool get_edge_val(int layer_num, const vector<bool>& now_mask, bool added_val) {
    vector<pair<int, bool>> row_and_value;
    int added_row = activity_rows[layer_num + 1].added;

    if (added_row != -1) {
        row_and_value.emplace_back(added_row, added_val);
    }

    for (int i = 0; i < activity_rows[layer_num].rows.size(); i++) {
        row_and_value.emplace_back(activity_rows[layer_num].rows[i], now_mask[i]);
    }

    vector<bool> from_gSpan;
    vector<bool> from_vertexes;

    for (int i = 0; i < row_and_value.size(); i++) {
        int row = row_and_value[i].first;

        from_gSpan.push_back(gSpan[row][layer_num]);
        from_vertexes.push_back(row_and_value[i].second);
    }

    return mul_two_vectors(from_vertexes, from_gSpan);
}

void construct_graph() {
    int vertexes_was = 0;

    for (int i = 0; i < activity_rows.size() - 1; i++) {
        int now_v_cnt = (int) activity_rows[i].rows.size();
        int next_v_cnt = (int) activity_rows[i + 1].rows.size();

        for (int j = 0; j < (1 << now_v_cnt); j++) {
            pair<int, int> next_v = {-1, -1};
            vector<bool> mask = get_mask(j, now_v_cnt);

            vector<bool> new_mask = gen_mask(i, mask, false);
            bool edge_val = get_edge_val(i, mask, false);
            int parsed_num = parse_from_mask(new_mask, next_v_cnt);
            if (edge_val) {
                next_v.second = parsed_num + vertexes_was + (1 << now_v_cnt);
            } else {
                next_v.first = parsed_num + vertexes_was + (1 << now_v_cnt);
            }

            if (activity_rows[i + 1].added != -1) {
                new_mask = gen_mask(i, mask, true);
                edge_val = get_edge_val(i, mask, true);
                parsed_num = parse_from_mask(new_mask, next_v_cnt);
                if (edge_val) {
                    next_v.second = parsed_num + vertexes_was + (1 << now_v_cnt);
                } else {
                    next_v.first = parsed_num + vertexes_was + (1 << now_v_cnt);
                }
            }

            graph.push_back(next_v);
        }

        vertexes_was += (1 << now_v_cnt);
    }

    graph.emplace_back(-1, -1);
}

vector<bool> decode(const vector<double>& v) {
    Data dp[graph.size()];
    for (int i = 0; i < graph.size(); i++) {
        dp[i] = Data{-1e9, -1, -1};
    }

    dp[0] = Data{0, -1, 0};
    for (int i = 0; i < graph.size(); i++) {
        double now_sum = dp[i].sum;
        int now_depth = dp[i].depth;
        auto data = graph[i];

        if (data.first != -1 && dp[data.first].sum < now_sum + v[now_depth]) {
            dp[data.first] = {now_sum + v[now_depth], i, now_depth + 1};
        }

        if (data.second != -1 && dp[data.second].sum < now_sum - v[now_depth]) {
            dp[data.second] = {now_sum - v[now_depth], i, now_depth + 1};
        }
    }

    int curr = (int) graph.size() - 1;
    vector<bool> res;
    while (dp[curr].parent != -1) {
        int before = dp[curr].parent;
        res.push_back(graph[before].second == curr);
        curr = before;
    }

    reverse(res.begin(), res.end());
    return res;
}

vector<bool> gen_word() {
    vector<bool> res;

    for (int i = 0; i < k; i++) {
        int rand_gen = round(((double) rand() / (RAND_MAX)));
        res.push_back(rand_gen == 1);
    }

    return res;
}

double simulate(double noise_lvl, int num_of_operations, int max_errors) {
    int curr_operations = 0;
    int curr_errors = 0;

    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<> nd{0, sqrt(0.5 * pow(10,  -noise_lvl / 10.0) * ((double) n / k))};

    while (curr_operations < num_of_operations && curr_errors < max_errors) {
        curr_operations++;

        vector<bool> word = gen_word();
        vector<bool> encoded_word = encode(word);

        vector<double> converted_word;
        for (int i = 0; i < encoded_word.size(); i++) {
            converted_word.push_back(encoded_word[i] ? -1.0 : 1.0);
            converted_word[i] += nd(gen);
        }

        vector<bool> decoded_word = decode(converted_word);

        for (int i = 0; i < encoded_word.size(); i++) {
            if (encoded_word[i] != decoded_word[i]) {
                curr_errors++;
                break;
            }
        }

    }

    return (double) curr_errors / curr_operations;
}

int main() {
//    ifstream in("/Users/st1580/CLionProjects/professional/input.txt");
//    ofstream out("/Users/st1580/CLionProjects/professional/output.txt");
    ifstream in("input.txt");
    ofstream out("output.txt");

    ios_base::sync_with_stdio(false);
    in.tie(0);
    out.tie(0);

    in >> n >> k;

    vector<bool> empty;
    for (int i = 0; i < k; i++) {
        g.push_back(empty);
        for (int j = 0; j < n; j++) {
            int x;
            in >> x;
            g[i].push_back(x == 1);
        }
    }

    construct_span_matrix();
    construct_graph();

    for (int i = 0; i < activity_rows.size(); i++) {
        out << (1 << activity_rows[i].rows.size()) << " ";
    }

    string s;
    while (in >> s) {
        out << "\n";

        if (s[0] == 'E') {
            vector<bool> v;
            for (int i = 0; i < k; i++) {
                int x;
                in >> x;
                v.push_back(x == 1);
            }

            for (bool val : encode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else if (s[0] == 'D') {
            vector<double> v;
            for (int i = 0; i < n; i++) {
                double x;
                in >> x;
                v.push_back(x);
            }

            for (bool val : decode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else {
            double noise_lvl;
            int num_of_operations, max_errors;

            in >> noise_lvl >> num_of_operations >> max_errors;
            out << setprecision(15) << simulate(noise_lvl, num_of_operations, max_errors);
        }
    }

    in.close();
    out.close();

    return 0;
}
