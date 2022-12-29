#include <bits/stdc++.h>

using namespace std;

// структура для информации о текущем слое в решетке
// rows - активные строки на текущем слое
// removed - слой, который уйдет в следующем слое (-1 если никакой)
// added - слой, который пришел в этом слое (-1 если никакой)
struct Activity {
    vector<int> rows;
    int removed = -1;
    int added = -1;
};

// структура для подсчета динамики в решетке
// sum - текущая насчитаная сумма
// parent - вершина из которой пришли
// depth - глубина от стартовой вершины
struct Data {
    double sum;
    int parent;
    int depth;
};

int n, k;
vector<vector<bool>> g;         // порождающая матрица
vector<vector<bool>> gSpan;     // порождающая матрица в МСФ
vector<Activity> activity_rows; // вектор слоев
vector<pair<int, int>> graph;   // решетка, first - ребро по 0; second - ребро по 1

// добавить к строке row_to строку row_from из gSpan
void sum_matrix_rows(int row_from, int row_to) {
    for (int j = 0; j < n; j++) {
        gSpan[row_to][j] = gSpan[row_to][j] ^ gSpan[row_from][j];
    }
}

// делаем все концы спэнов уникальными
// curr_column - текущая колонка для проверки уникальности концов спэнов
// used - массив уже проверенных строк
void make_unique_ends(int curr_column, bool used[]) {
    // строка, в которой начало спэна в этом столбце
    int latest_true_row = -1;
    for (int i = k - 1; i >= 0; i--) {
        if (!used[i] && gSpan[i][curr_column]) {
            if (latest_true_row == -1) {
                latest_true_row = i;
                // эта строка будет в конце итерации единственной с концом спена в curr_column
                used[i] = true;
                continue;
            }

            // будем в этом месте, когда в текущем столбце начинаются два спена
            // это строки i и latest_true_row

            // ищем строку у которой конец спена наибольший
            // прибавляем к одной строке другую
            sum_matrix_rows(latest_true_row, i);
        }
    }
}

// делаем все начала спэнов уникальными
// curr_column - текущая колонка для проверки уникальности начал спэнов
// used - массив уже проверенных строк
void make_unique_starts(int curr_column, bool used[]) {
    // строка, в которой начало спэна в этом столбце
    int latest_true_row = -1;
    for (int i = 0; i < k; i++) {
        if (!used[i] && gSpan[i][curr_column]) {
            if (latest_true_row == -1) {
                latest_true_row = i;
                continue;
            }

            // будем в этом месте, когда в текущем столбце начинаются два спена
            // это строки i и latest_true_row

            // ищем строку у которой конец спена наибольший
            // добавляем к такой строке другую
            // необходимо, чтобы не сломать уникальность концов спэнов
            for (int j = n - 1; j >= 0; j--) {
                if (gSpan[i][j]) {
                    sum_matrix_rows(latest_true_row, i);
                    break;
                }
                if (gSpan[latest_true_row][j]) {
                    sum_matrix_rows(i, latest_true_row);
                    // необходимо сохранить инвариант last_true_row
                    latest_true_row = i;
                    break;
                }
            }
        }
    }

    // помечаем строку, в которой текущей столбец - начало спэна
    if (latest_true_row != -1) {
        used[latest_true_row] = true;
    }
}

// строим слои решетки
void construct_activities() {
    // массив, которой хранит начало и конец спэна для каждой строки
    vector<pair<int, int>> activity_range;
    // подсчет начала и конца спжна для каждой строки
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
    // первый слой состоит из одной вершины, поэтому заполним его пустыми данными
    activity_rows.push_back(empty_activity);

    // заполняем слои
    for (int j = 0; j < n; j++) {
        vector<int> active_rows;
        int removed = -1;
        int added = -1;

        for (int i = 0; i < k; i++) {
            // если текущая строка активна в j-том столбце, то добавляем ее
            if (activity_range[i].first <= j && j < activity_range[i].second) {
                active_rows.push_back(i);
            }

            // если начало спэна строки в текущем столбце, то помечаем это
            if (activity_range[i].first == j) {
                added = i;
            }

            // если конец спэна строки в следующем столбце, то помечаем это
            if (activity_range[i].second == j + 1) {
                removed = i;
            }
        }

        // добавляем данные по текущему слою
        activity_rows.emplace_back(Activity{active_rows, removed, added});
    }
}

// строим порождающую матрицу в МСФ
void construct_span_matrix() {
    // изначально она равна порождающей
    for (int i = 0; i < k; i++) {
        gSpan.push_back(g[i]);
    }

    // массив показывающий, что конец спэна строки i - уникален
    bool used[k];

    // делаем уникальными концы спэнов
    for (int j = n - 1; j >= 0; j--) {
        make_unique_ends(j, used);
    }

    // массив показывающий, что начало спэна строки i - уникално
    for (int i = 0; i < k; i++) {
        used[i] = false;
    }

    // делаем уникальными начала спэнов
    for (int j = 0; j < n; j++) {
        make_unique_starts(j, used);
    }

    // на текущем этапе получена порождающая матрица в МСФ
    // подсчитаем начало и концов спэнов и для каждого слоя найдем какие строки активны в нем
    construct_activities();
}

// переводим число в битовую маску
// power - размер маски
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

// строим битовую маску значений информационных символов для текущей вершины
// порядок информационных символов в маске соответсвует упорядоченным по возрастанию номерам строк
// layer_num - номер текщуго слоя
// mask - маска значений информационных символов на предыдущем слое
// added_val - значение, которое будет у нового для текущего слоя информационного символа
vector<bool> gen_mask(int layer_num, const vector<bool>& mask, bool added_val) {
    // итоговая битовая маска значений информационных символов
    vector<bool> res;
    // итоговое число информационных символов
    size_t must_be = activity_rows[layer_num + 1].rows.size();

    int i = 0;
    int j = 0;
    // при помощи методов двух указателей посортируем номера строк
    // поскольку на новом слое может стать активной одна строка (из-за посторения МСФ)
    // то достаточно обратить внимание на появившуюся строку,
    // ее и будем пытаться вставить в правильное место (упорядоченное по возрастанию номеров строк)
    while (i < activity_rows[layer_num].rows.size() && j < activity_rows[layer_num + 1].rows.size()) {
        if (activity_rows[layer_num].rows[i] == activity_rows[layer_num + 1].rows[j]) {
            // если у нас совпали номера строк на предыдущем и текущем слое => эти строки не добавились и не удалились
            // берем уже существующее значение
            res.push_back(mask[i]);
            i++; j++;
        } else {
            if (activity_rows[layer_num].removed == activity_rows[layer_num].rows[i]) {
                if (activity_rows[layer_num + 1].added == activity_rows[layer_num + 1].rows[j]) {
                    // случай когда указатели стоят на удалившейся и добавившейся строке
                    // значит необходимо добавить новое значение и сдвинуть оба указателя
                    res.push_back(added_val);
                    i++; j++;
                } else {
                    // случай когда указатель стоит на удалившейся строке
                    // значит необходимо сдвинуть соответствующий указатель не записывая значение
                    i++;
                }
            } else {
                if (activity_rows[layer_num + 1].added == activity_rows[layer_num + 1].rows[j]) {
                    // случай когда указатель стоит на добавившейся строке
                    // значит необходимо добавить новое значение и сдвинуть соответствующий указатель
                    res.push_back(added_val);
                    j++;
                }
            }
        }
    }

    // заполняем оставшиеся информационные символы
    while (res.size() < must_be && i < activity_rows[layer_num].rows.size()) {
        res.push_back(mask[i]);
        i++;
    }

    // добавляем новый информационный символ, если он еще не добавлен
    if (j < activity_rows[layer_num + 1].rows.size()) {
        res.push_back(added_val);
    }

    return res;
}

// получаем число из маски
int parse_from_mask(const vector<bool>& mask) {
    int res = 0;
    for (int i = 0; i < mask.size(); i++) {
        if (mask[i]) {
            res += (1 << i);
        }
    }

    return res;
}

// умножение двух векторов
bool mul_two_vectors(const vector<bool>& a, const vector<bool>& b) {
    bool res = false;

    for (int i = 0; i < a.size(); i++) {
        res = res ^ (a[i] & b[i]);
    }

    return res;
}

// получаем значение на ребре между текущим слоем и будущим
// это вектор a - значений информационных символах на вершинах
// перемноженный с вектором b - значений из gSpan
// (строки - активные строки на текущих слоях; столбец - текущий слой)
// layer_num - номер текщуго слоя
// mask - маска значений информационных символов на текущем слое
// added_val - значение, которое будет у нового информационного символа
bool get_edge_val(int layer_num, const vector<bool>& now_mask, bool added_val) {
    // вектор значений для каждой строки
    vector<pair<int, bool>> row_and_value;
    // добавленная на будущем слое строка
    int added_row = activity_rows[layer_num + 1].added;

    if (added_row != -1) {
        row_and_value.emplace_back(added_row, added_val);
    }

    for (int i = 0; i < activity_rows[layer_num].rows.size(); i++) {
        row_and_value.emplace_back(activity_rows[layer_num].rows[i], now_mask[i]);
    }

    // значения из порождающей матрицы в МСФ
    vector<bool> from_gSpan;
    // значения информационных символов
    vector<bool> from_vertexes;

    // для каждой активной строки находим соответствующее значение в gSpan
    for (int i = 0; i < row_and_value.size(); i++) {
        int row = row_and_value[i].first;

        from_gSpan.push_back(gSpan[row][layer_num]);
        from_vertexes.push_back(row_and_value[i].second);
    }

    // возвращаем результат умножения вектров
    return mul_two_vectors(from_vertexes, from_gSpan);
}

// строим решетку
void construct_graph() {
    // сколько вершин уже в решетке
    int vertexes_was = 0;

    // итерируемся по слоям
    for (int i = 0; i < activity_rows.size() - 1; i++) {
        // число активных строк на текущем слое
        int now_v_cnt = (int) activity_rows[i].rows.size();

        // итерируемся по вершинам текущего слоя
        for (int j = 0; j < (1 << now_v_cnt); j++) {
            // следующие вершины
            pair<int, int> next_v = {-1, -1};
            // генерим информационные символы текущей вершины
            vector<bool> mask = get_mask(j, now_v_cnt);

            // генерим информационные символы следующей вершины
            vector<bool> new_mask = gen_mask(i, mask, false);
            // получаем значение на ребре
            bool edge_val = get_edge_val(i, mask, false);
            // получаем индекс следующей вершины
            int parsed_num = parse_from_mask(new_mask);

            // добавляем ребро
            if (edge_val) {
                next_v.second = parsed_num + vertexes_was + (1 << now_v_cnt);
            } else {
                next_v.first = parsed_num + vertexes_was + (1 << now_v_cnt);
            }

            // если на следующем слое стала активной новая строка
            // то мы сможем перейти в две новых вершины
            // ту, у которой новый информационный символ 0
            // и в ту, у которой новый информационный симвлол 1
            if (activity_rows[i + 1].added != -1) {
                // повторяем уже известные шаги
                new_mask = gen_mask(i, mask, true);
                edge_val = get_edge_val(i, mask, true);
                parsed_num = parse_from_mask(new_mask);
                if (edge_val) {
                    next_v.second = parsed_num + vertexes_was + (1 << now_v_cnt);
                } else {
                    next_v.first = parsed_num + vertexes_was + (1 << now_v_cnt);
                }
            }

            // добавляем в решетку текущую вершину
            graph.push_back(next_v);
        }

        vertexes_was += (1 << now_v_cnt);
    }

    // добавляем финишную вершину
    graph.emplace_back(-1, -1);
}

vector<bool> encode(const vector<bool>& v) {
    vector<bool> res;

    // простое перемножение пришедшего вектора с порождающей матрицей
    for (int j = 0; j < n; j++) {
        res.push_back(false);
        for (int i = 0; i < k; i++) {
            res[j] = res[j] ^ (v[i] & g[i][j]);
        }
    }

    return res;
}

vector<bool> decode(const vector<double>& v) {
    // инициализируем динамику
    Data dp[graph.size()];
    for (int i = 0; i < graph.size(); i++) {
        dp[i] = Data{-1e9, -1, -1};
    }

    dp[0] = Data{0, -1, 0};
    // простая динамика
    // в качестве метрики используется скалярное произведение
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

    // восстанавливаем ответ (путь по которому прошли через решетку)
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

// генерация вектора для encode
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

        // перевод 1 - 2х и зашумлление
        vector<double> converted_word;
        for (int i = 0; i < encoded_word.size(); i++) {
            converted_word.push_back(encoded_word[i] ? -1.0 : 1.0);
            converted_word[i] += nd(gen);
        }

        vector<bool> decoded_word = decode(converted_word);

        // проверка
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

    // вывод размера слоев
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