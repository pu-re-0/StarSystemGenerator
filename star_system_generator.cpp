#include <iostream>
#include <vector>
#include <string>
#include <cmath>     // For std::acos, std::pow, std::log10, std::cbrt, std::asin
#include <random>    // For std::random_device, std::mt19937, std::uniform_real_distribution, std::discrete_distribution
#include <fstream>   // For std::ofstream
#include <iomanip>   // For std::fixed, std::setprecision
#include <map>       // For std::map
#include <limits>    // For std::numeric_limits
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::clamp
#include "cluster_generator.h"
#include "star_map_drawer.h"
#include <iterator>  // For std::make_move_iterator
#include <omp.h>     // OpenMPヘッダー
#include <optional>  // For std::optional
#include <cctype>    // for std::isdigit
#include <sstream>   // For std::stringstream

#ifdef _WIN32
// windows.h が min/max マクロを定義するのを防ぐ
// これにより std::numeric_limits<T>::max() などとの衝突を回避できる
#define NOMINMAX
#include <windows.h> // SetConsoleOutputCPのため
#endif

// 3Dビューア専用の軽量CSVを書き出す関数
void output_3d_viewer_csv(const std::vector<StarSystem>& all_stars, const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "エラー: 3Dビューア用CSVファイル " << filename << " を書き込み用に開けませんでした。" << std::endl;
        return;
    }

    // UTF-8 BOM
    ofs << "\xEF\xBB\xBF";

    // ヘッダー: main.jsのparseCSVが必要とする列のみ
    ofs << "x,y,z,スペクトル,星系視等級\n";

    // --- 高速化のためのバッファリング書き込み ---
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4);
    const size_t buffer_flush_size = 4096; // この行数ごとにファイルに書き出す
    size_t line_count = 0;

    for (const auto& s : all_stars) {
        // 伴星（等級-1）はビューアに不要なのでスキップする
        if (s.system_apparent_magnitude == -1.0f) {
            continue;
        }

        ss << s.x << ","
           << s.y << ","
           << s.z << ","
           << spectral_type_to_string(s.spectral_type) << ","
           << s.system_apparent_magnitude << "\n";
        line_count++;

        if (line_count >= buffer_flush_size) {
            ofs << ss.rdbuf(); // バッファの内容をファイルに書き出す
            ss.str("");      // バッファをクリア
            ss.clear();      // ストリームの状態をリセット
            line_count = 0;
        }
    }

    // ループ終了後、バッファに残っているデータを書き出す
    if (line_count > 0) {
        ofs << ss.rdbuf();
    }

    ofs.close();
    std::cout << "3Dビューア用のデータが " << filename << " に正常に書き込まれました。" << std::endl;
}

// Helper to split a string by a delimiter
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// CSVから星と星団のデータを読み込む関数。等級によるフィルタリング機能付き。
std::optional<std::pair<std::vector<StarSystem>, std::vector<ClusterProperties>>>
load_data_from_csv(const std::string& filename, std::optional<float> max_magnitude_filter = std::nullopt) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "エラー: ファイル " << filename << " を開けませんでした。" << std::endl;
        return std::nullopt;
    }

    // --- メモリ確保の効率化 ---
    // ファイルサイズから読み込む星のおおよその数を推測し、メモリを事前に確保する。
    // これにより、ループ中の複数回のメモリ再確保を防ぎ、断片化を抑制し、bad_allocのリスクを低減する。
    ifs.seekg(0, std::ios::end);
    std::streampos file_size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    // --- UTF-8 BOMの読み飛ばし処理 ---
    // ファイルの先頭3バイトを読み込み、BOM（Byte Order Mark）かどうかを判定します。
    // BOMが存在する場合、それ以降からファイルの読み込みを開始します。
    char bom[3];
    ifs.read(bom, 3);
    if (!(ifs.gcount() == 3 && bom[0] == (char)0xEF && bom[1] == (char)0xBB && bom[2] == (char)0xBF)) {
        // BOMでなかった場合（またはファイルが3バイト未満の場合）、ファイルポインタを先頭に戻します。
        ifs.seekg(0);
    }

    std::vector<StarSystem> stars;
    std::vector<ClusterProperties> clusters;

    if (file_size > 0) {
        // 1行あたり平均120バイトと仮定して、必要な要素数を見積もる
        size_t estimated_stars = static_cast<size_t>(file_size) / 120;
        stars.reserve(estimated_stars);
    }

    std::string line;
    bool header_skipped = false;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;

        // コメント行の処理
        if (line[0] == '#') {
            // 星団プロパティの行（データまたはヘッダー）を処理
            if (line.rfind("# CLUSTER_PROP,", 0) == 0) { // C++20のstarts_withの代わり
                std::string prop_data = line.substr(std::string("# CLUSTER_PROP,").length());
                // 最初の文字が数字である場合のみ、データ行として解析を試みる
                // これにより、"id"で始まるヘッダー行を安全にスキップできる
                if (!prop_data.empty() && std::isdigit(static_cast<unsigned char>(prop_data[0]))) {
                    std::vector<std::string> tokens = split(prop_data, ',');
                    if (tokens.size() == 11) {
                        try {
                            ClusterProperties prop;
                            prop.id = std::stoi(tokens[0]);
                            prop.center_x = std::stof(tokens[1]);
                            prop.center_y = std::stof(tokens[2]);
                            prop.center_z = std::stof(tokens[3]);
                            prop.axis_a = std::stof(tokens[4]);
                            prop.axis_b = std::stof(tokens[5]);
                            prop.axis_c = std::stof(tokens[6]);
                            prop.rot_alpha = std::stof(tokens[7]);
                            prop.rot_beta = std::stof(tokens[8]);
                            prop.rot_gamma = std::stof(tokens[9]);
                            prop.is_old = std::stoi(tokens[10]) != 0;
                            clusters.push_back(prop);
                        } catch (const std::exception& e) {
                            std::cerr << "警告: 星団プロパティ行の解析に失敗しました: " << line << " (" << e.what() << ")" << std::endl;
                        }
                    } else {
                        std::cerr << "警告: 星団プロパティ行の列数が不正です（期待値: 11, 実際: " << tokens.size() << "）。行をスキップします: " << line << std::endl;
                    }
                }
            }
            continue; // 他のコメント行はスキップ
        }

        // ヘッダー行をスキップ
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }

        // データ行の処理
        std::vector<std::string> tokens = split(line, ',');
        if (tokens.size() == 15) {
            try {
                // フィルタリングのために、先に視等級をチェックする
                float system_apparent_magnitude = std::stof(tokens[13]);

                // 伴星は描画しないため、常にスキップ
                if (system_apparent_magnitude == -1.0f) {
                    continue;
                }

                // 等級フィルタが指定されていれば、適用する
                if (max_magnitude_filter && system_apparent_magnitude > *max_magnitude_filter) {
                    continue; // この星は暗すぎるのでスキップ
                }

                StarSystem s;
                s.id = std::stoll(tokens[0]);
                s.binary_system_id = std::stoll(tokens[1]);
                s.apparent_magnitude = std::stof(tokens[2]);
                s.spectral_type = string_to_spectral_type(tokens[3]);
                s.absolute_magnitude = std::stof(tokens[4]);
                s.distance_r = std::stof(tokens[5]); // 距離r
                s.longitude_theta = std::stof(tokens[6]); // 経度θ
                s.latitude_phi = std::stof(tokens[7]); // 緯度φ
                s.x = std::stof(tokens[8]); // x
                s.y = std::stof(tokens[9]); // y
                s.z = std::stof(tokens[10]); // z
                s.cluster_id = std::stoi(tokens[11]); // 星団ID
                s.system_absolute_magnitude = std::stof(tokens[12]); // 星系絶対等級
                s.system_apparent_magnitude = system_apparent_magnitude; // 星系視等級（フィルタリングでパース済み）
                if (!tokens[14].empty()) {
                    s.name = std::make_shared<std::string>(tokens[14]);
                }
                stars.push_back(s);
            } catch (const std::exception& e) {
                std::cerr << "警告: データ行の解析に失敗しました: " << line << " (" << e.what() << ")" << std::endl;
            }
        } else if (!tokens.empty()) { // 空行でない場合のみ警告
            std::cerr << "警告: データ行の列数が不正です（期待値: 15, 実際: " << tokens.size() << "）。行をスキップします: " << line << std::endl;
        }
    }

    if (stars.empty()) {
        std::cerr << "エラー: " << filename << " から星データを読み込めませんでした。ファイルが空か、形式が正しくない可能性があります。" << std::endl;
        return std::nullopt;
    }

    std::cout << filename << " から " << stars.size() << " 個の星と " << clusters.size() << " 個の星団データを読み込みました。" << std::endl;
    return std::make_pair(stars, clusters);
}

// 描画パラメータをユーザーから取得するヘルパー関数
void get_drawing_parameters_from_user(
    bool& show_celestial_lines,
    float& max_magnitude_for_svg,
    bool& show_cluster_regions,
    bool& show_star_labels,
    float& label_magnitude_threshold
);

// 星生成のメインロジックをこの関数にカプセル化
int run_simulation(unsigned int seed) {
    std::cout << "\n--- シミュレーション開始 (シード値: " << seed << ") ---\n";

    // 1. Setup random number generator with the given seed
    std::mt19937 gen(seed);

    // --- 星の名前リストをファイルから読み込む ---
    std::vector<std::string> star_names;
    std::ifstream name_file("starnamelist.txt");
    if (name_file.is_open()) {
        std::string name;
        while (std::getline(name_file, name)) {
            if (!name.empty()) {
                star_names.push_back(name);
            }
        }
        name_file.close();
        if (star_names.empty()) {
            std::cout << "starnamelist.txt は空です。星の名前は割り当てられません。" << std::endl;
        } else {
            std::cout << "starnamelist.txt から " << star_names.size() << " 個の星の名前を読み込みました。" << std::endl;
        }
    } else {
        std::cout << "警告: starnamelist.txt が見つかりません。星の名前は割り当てられません。" << std::endl;
    }

    // Helper function to normalize a vector of weights to sum to 1.0f
    auto normalize_weights = [](std::vector<float>& weights) {
        float total_weight = std::accumulate(weights.begin(), weights.end(), 0.0f);
        if (total_weight > std::numeric_limits<float>::epsilon()) {
            for (float& weight : weights) {
                weight /= total_weight;
            }
        }
    };

    // 2. Define spectral type properties (probabilities and absolute magnitude ranges)
    std::vector<std::string> spectral_types_vec = {"M", "K", "G", "F", "A", "B", "O", "WD", "ES"};
    std::vector<float> spectral_type_weights = {0.725, 0.129, 0.041, 0.029, 0.006, 0.00039, 0.0000005, 0.059, 0.011}; // Approximate weights, will be normalized
    normalize_weights(spectral_type_weights);

    // Absolute magnitude ranges for each spectral type {min_M, max_M}
    std::map<std::string, std::pair<float, float>> abs_mag_ranges;
    abs_mag_ranges["O"] = {-6.0f, -4.0f};
    abs_mag_ranges["B"] = {-4.5f, 2.0f};
    abs_mag_ranges["A"] = {1.5f, 4.0f};
    abs_mag_ranges["F"] = {3.5f, 5.5f};
    abs_mag_ranges["G"] = {4.5f, 6.2f};
    abs_mag_ranges["K"] = {5.8f, 7.6f};
    abs_mag_ranges["M"] = {7.4f, 12.0f};
    abs_mag_ranges["WD"] = {10.0f, 16.0f}; // White Dwarf
    abs_mag_ranges["ES"] = {-8.0f, 3.0f}; // non-Red Giant Evolved Stars

    // Multiplicity (平均伴星数) for each spectral type
    std::map<std::string, float> multiplicity_map;
    multiplicity_map["O"] = 1.3f;  // O型星は多重連星系を形成することが多い
    multiplicity_map["B"] = 1.0f;
    multiplicity_map["A"] = 1.0f;
    multiplicity_map["F"] = 0.62f;
    multiplicity_map["G"] = 0.62f; // 太陽のようなG型星
    multiplicity_map["K"] = 0.62f;
    multiplicity_map["M"] = 0.33f; // M型（赤色矮星）は連星率が低い
    multiplicity_map["WD"] = 0.1f; // 白色矮星
    multiplicity_map["ES"] = 0.3f; // 進化した星

    // --- 連星間距離のパラメータ設定 ---
    // 2つのピーク（近接連星と遠距離連星）を持つ分布を定義
    // { close_log_mean, close_log_stddev, wide_log_mean, wide_log_stddev, close_binary_prob }
    std::map<std::string, BinarySeparationParams> binary_separation_params;
    // O型: 近い伴星も遠い伴星も持ちうる (中央値: 8 AU, 200 AU)
    binary_separation_params["O"]  = {std::log(8.0f), 1.2f, std::log(200.0f), 1.5f, 0.5f};
    // B型 (中央値: 6 AU, 150 AU)
    binary_separation_params["B"]  = {std::log(6.0f), 1.2f, std::log(150.0f), 1.4f, 0.55f};
    // A型 (中央値: 5 AU, 80 AU)
    binary_separation_params["A"]  = {std::log(5.0f), 1.1f, std::log(80.0f),  1.3f, 0.6f};
    // F型 (中央値: 4 AU, 60 AU)
    binary_separation_params["F"]  = {std::log(4.0f), 1.1f, std::log(60.0f),  1.3f, 0.65f};
    // G型 (太陽) (中央値: 3 AU, 50 AU)
    binary_separation_params["G"]  = {std::log(3.0f), 1.0f, std::log(50.0f),  1.2f, 0.7f};
    // K型 (中央値: 2.5 AU, 40 AU)
    binary_separation_params["K"]  = {std::log(2.5f), 1.0f, std::log(40.0f),  1.2f, 0.75f};
    // M型 (赤色矮星): ほとんどが近接連星 (中央値: 2 AU, 20 AU)
    binary_separation_params["M"]  = {std::log(2.0f), 0.9f, std::log(20.0f),  1.1f, 0.85f};
    // 白色矮星 (中央値: 1 AU, 30 AU)
    binary_separation_params["WD"] = {std::log(1.0f), 0.9f, std::log(30.0f),  1.2f, 0.8f};
    // 進化した星: 軌道が広がっていることが多い (中央値: 10 AU, 300 AU)
    binary_separation_params["ES"] = {std::log(10.0f), 1.3f, std::log(300.0f), 1.6f, 0.4f};


    // For open clusters (younger, hotter stars are more common)
    // Weights for:                                M,    K,    G,    F,    A,    B,      O,        WD,    ES
    std::vector<float> cluster_spectral_weights = {0.60f, 0.20f, 0.10f, 0.05f, 0.02f, 0.004f, 0.0005f, 0.015f, 0.0105f}; // Approximate weights, will be normalized
    normalize_weights(cluster_spectral_weights);

    // 3. Define the radius of the spherical space
    float R_parsecs = 100.0f; // Default radius in parsecs (pc)
    std::cout << "半径Rパーセクの球形状空間の半径を入力してください (デフォルト: 100.0f): ";
    std::string line;
    std::getline(std::cin, line);
    if (!line.empty()) {
        try {
            float input_r = std::stof(line);
            if (input_r > 0) {
                R_parsecs = input_r;
            } else {
                std::cout << "半径は正の値である必要があります。デフォルト値 " << R_parsecs << " を使用します。" << std::endl;
            }
        } catch (const std::exception&) {
            std::cout << "無効な入力です。デフォルト値 " << R_parsecs << " を使用します。" << std::endl;
        }
    }

    // 描画パラメータをユーザーから取得
    bool show_celestial_lines = true;
    float max_magnitude_for_svg = 5.0f;
    bool show_cluster_regions = true;
    bool show_star_labels = true;
    float label_magnitude_threshold = 2.5f;
    get_drawing_parameters_from_user(
        show_celestial_lines,
        max_magnitude_for_svg,
        show_cluster_regions,
        show_star_labels,
        label_magnitude_threshold
    );

    // 4. Calculate the number of stars
    // Correct volume of a sphere: (4/3) * PI * R^3
    const float PI = std::acos(-1.0f); // Portable way to get PI
    float volume = (4.0f / 3.0f) * PI * std::pow(R_parsecs, 3);
    float star_density = 0.13f; // Star system density in pc^-3

    long long num_stars = static_cast<long long>(volume * star_density);
    if (num_stars == 0) { // Ensure at least one star if the calculated number is zero
        num_stars = 1;
    }

    // 安全装置: 生成する星の数が多すぎるとメモリを大量に消費しクラッシュするため、上限を設ける
    const long long MAX_STARS_LIMIT = 50000000; // 5千万個の主星系を上限とする (メモリ使用量 約7.5GBに相当)
    if (num_stars > MAX_STARS_LIMIT) {
        std::cout << "\n警告: 計算された星の数 (" << num_stars << ") が上限 (" << MAX_STARS_LIMIT << ") を超えています。" << std::endl;
        std::cout << "このまま実行すると、メモリ不足でプログラムがクラッシュする可能性があります。" << std::endl;
        std::cout << "星の数を上限である " << MAX_STARS_LIMIT << " 個に制限して続行しますか？ (yes/no, デフォルト: no): ";
        std::string user_choice;
        std::getline(std::cin, user_choice);
        if (user_choice == "yes" || user_choice == "y") {
            num_stars = MAX_STARS_LIMIT;
            std::cout << "星の数を " << num_stars << " 個に制限して続行します。" << std::endl;
        } else {
            std::cout << "シミュレーションをキャンセルしました。メニューに戻ります。" << std::endl;
            return 1; // キャンセルを示すコードを返す
        }
    }

    std::cout << R_parsecs << "パーセクの半径を持つ球状空間内に " << num_stars << " 個の星系を計算します。" << std::endl;
    std::cout << "半径が非常に大きい場合、時間がかかることがあります。" << std::endl;

    // 6. Discrete distribution for spectral types
    std::discrete_distribution<> spectral_dist(spectral_type_weights.begin(), spectral_type_weights.end());

    // 7. Vector to store star systems
    std::vector<StarSystem> star_systems; // Initially empty

    // 連星系IDを採番するためのカウンター
    long long next_binary_system_id = 1;

    // 8. Generate stars (clusters and field stars)

    // First, generate open clusters
    std::cout << "\n--- 散開星団の生成 ---" << std::endl;
    auto [cluster_stars, cluster_properties] = generate_open_clusters(
        num_stars,
        R_parsecs,
        gen,
        spectral_types_vec,
        cluster_spectral_weights,
        abs_mag_ranges,
        multiplicity_map,
        binary_separation_params,
        next_binary_system_id // 星団生成後にIDカウンターが更新される
    );

    // Count the number of primary systems generated in clusters to accurately calculate remaining field stars.
    // A primary system is identified by its system_apparent_magnitude not being -1.0f.
    long long num_cluster_primary_stars = 0;
    for(const auto& star : cluster_stars) {
        if (star.system_apparent_magnitude != -1.0f) {
            num_cluster_primary_stars++;
        }
    }

    // Move the generated cluster stars into the main vector to avoid copying.
    star_systems.reserve(cluster_stars.size());
    star_systems.insert(
        star_systems.end(),
        std::make_move_iterator(cluster_stars.begin()),
        std::make_move_iterator(cluster_stars.end())
    );
    cluster_stars.clear(); // Free up memory immediately
    cluster_stars.shrink_to_fit();

    // Then, generate the remaining field stars
    long long num_field_stars = num_stars - num_cluster_primary_stars;
    if (num_field_stars < 0) num_field_stars = 0; // Should not happen, but as a safeguard
    std::cout << "\n--- 現場星の生成 ---" << std::endl;
    std::cout << num_field_stars << " 個の主星系を並列で生成します..." << std::endl;

    std::vector<StarSystem> generated_field_stars; // This will collect stars from all threads
    // Reserve memory for field stars, also considering companions.
    // A factor of 1.5 is a reasonable guess.
    generated_field_stars.reserve(static_cast<size_t>(num_field_stars * 1.5));

    // --- 再現性を保証する決定論的な並列処理 ---
    // 1. 乱数生成のためのマスターシードを取得
    unsigned int master_seed = gen();
    // 2. 星団生成後に使われた最後の連星IDを取得
    const long long last_cluster_binary_id = next_binary_system_id;

    #pragma omp parallel
    {
        std::vector<StarSystem> local_stars; // 各スレッドが生成した星を一時的に格納するベクター
        // スレッドごとに一度だけ乱数生成器と分布オブジェクトを生成する
        std::mt19937 thread_gen;
        std::discrete_distribution<> thread_spectral_dist(spectral_type_weights.begin(), spectral_type_weights.end());

        // 5. Distributions for position (uniform distribution within a sphere)
        // スレッドセーフ性を確保するため、各スレッドが自分自身の分布オブジェクトを持つように、
        // parallelブロック内で宣言します。これにより、複数のスレッドが同時に分布オブジェクトの
        // 内部状態を変更しようとするデータ競合を防ぎます。
        std::uniform_real_distribution<float> dist_r_cube(std::numeric_limits<float>::epsilon(), 1.0f);
        std::uniform_real_distribution<float> dist_theta(0.0f, 2.0f * PI); // Longitude [0, 2*PI)
        std::uniform_real_distribution<float> dist_phi_sin(-1.0f, 1.0f); // For sin(phi) [-1, 1]


        // ループの反復をスレッドに分割
        #pragma omp for schedule(dynamic, 100) // チャンクサイズを指定してスケジューリングのオーバーヘッドを削減
        for (long long i = 0; i < num_field_stars; ++i) {
            // 3. 各星系ごとに決定論的な乱数生成器を初期化する
            // iがunsigned intの最大値を超えてもユニークなシードを生成するための堅牢な方法
            unsigned int seed_high = static_cast<unsigned int>(i >> 32); // long longの上位32ビット
            unsigned int seed_low = static_cast<unsigned int>(i);      // long longの下位32ビット
            std::seed_seq seed_sequence = {master_seed, seed_high, seed_low};
            thread_gen.seed(seed_sequence); // 構築済みの生成器にシードを設定する

            StarSystem s;
            s.cluster_id = 0; // Field star
            // 5. 連星IDを決定論的に割り当て
            s.binary_system_id = last_cluster_binary_id + i;

            // Generate spherical position
            s.distance_r = R_parsecs * std::cbrt(dist_r_cube(thread_gen));
            s.longitude_theta = dist_theta(thread_gen);
            s.latitude_phi = std::asin(dist_phi_sin(thread_gen));

            // Convert spherical to Cartesian coordinates
            float cos_phi = std::cos(s.latitude_phi);
            s.x = s.distance_r * cos_phi * std::cos(s.longitude_theta);
            s.y = s.distance_r * cos_phi * std::sin(s.longitude_theta);
            s.z = s.distance_r * std::sin(s.latitude_phi);

            // 共通関数を使って星の核となる特性を生成
            generate_star_core_properties(s, thread_gen, spectral_types_vec, thread_spectral_dist, abs_mag_ranges);

            // --- 連星生成ロジック ---
            StarGenerationParameters field_star_params = {spectral_types_vec, abs_mag_ranges, thread_spectral_dist, binary_separation_params};
            std::vector<StarSystem> companions_to_add = generate_companions_for_primary(
                s,
                multiplicity_map,
                thread_gen,
                field_star_params
            );

            // Calculate apparent magnitude
            s.apparent_magnitude = s.absolute_magnitude + 5.0f * (std::log10(s.distance_r) - 1.0f);

            // 星系全体の等級を計算して設定
            set_system_magnitudes(s, companions_to_add);

            // 主星をローカルリストに追加
            local_stars.push_back(s);
            // 生成された伴星をローカルリストに追加
            for (const auto& comp : companions_to_add) {
                local_stars.push_back(comp);
            }
        }

        // 各スレッドが生成した星を、排他制御しながらグローバルなベクターにマージする
        #pragma omp critical
        {
            // Move elements from local_stars to generated_field_stars to avoid expensive copies.
            // This is a critical optimization for memory usage in a parallel context.
            if (!local_stars.empty()) {
                generated_field_stars.insert(
                    generated_field_stars.end(),
                    std::make_move_iterator(local_stars.begin()),
                    std::make_move_iterator(local_stars.end())
                );
                local_stars.clear(); // Explicitly clear to release memory.
            }
        }
    }
    // Reserve enough space in the main vector to avoid reallocations during the final merge.
    star_systems.reserve(star_systems.size() + generated_field_stars.size());
    // Move the generated field stars into the main vector to avoid copying
    star_systems.insert(
        star_systems.end(),
        std::make_move_iterator(generated_field_stars.begin()),
        std::make_move_iterator(generated_field_stars.end())
    );
    // Assign unique IDs to all generated stars (clusters + field)
    std::cout << "\n--- 最終処理 ---" << std::endl;
    std::cout << "全星系にユニークIDを割り当てています..." << std::endl;
    for(long long i = 0; i < star_systems.size(); ++i) {
        star_systems[i].id = i + 1;
    }

    // Group stars by binary system ID to assign names correctly
    std::cout << "連星系ごとに名前を割り当てています..." << std::endl;
    std::map<long long, std::vector<StarSystem*>> binary_systems_map;
    for (auto& star : star_systems) {
        binary_systems_map[star.binary_system_id].push_back(&star);
    }

    // Prepare for random name assignment from the list
    std::uniform_int_distribution<size_t> name_dist(0, star_names.empty() ? 0 : star_names.size() - 1);

    // Assign names to each system
    for (auto& pair : binary_systems_map) {
        auto& stars_in_system = pair.second;

        if (stars_in_system.empty()) {
            continue;
        }

        // Find the primary star in the system.
        // The primary is the one whose system_apparent_magnitude is not -1.0f.
        StarSystem* primary = nullptr;
        for (auto* star_ptr : stars_in_system) {
            if (star_ptr->system_apparent_magnitude != -1.0f) {
                primary = star_ptr;
                break;
            }
        }
        // If no primary is found by the rule (e.g., single star system), take the first one.
        if (!primary) {
            primary = stars_in_system[0];
        }

        // Assign a random name to the primary star.
        if (!star_names.empty()) {
            primary->name = std::make_shared<std::string>("rd_" + star_names[name_dist(gen)]);
        } else {
            primary->name = nullptr;
        }

        // Collect and sort companions by absolute magnitude to assign names consistently.
        std::vector<StarSystem*> companions;
        for (auto* star_ptr : stars_in_system) {
            if (star_ptr != primary) {
                companions.push_back(star_ptr);
            }
        }

        // Sort companions by absolute magnitude (brighter first).
        std::sort(companions.begin(), companions.end(), [](const StarSystem* a, const StarSystem* b) {
            return a->absolute_magnitude < b->absolute_magnitude;
        });

        // Assign names to companions like "PrimaryName B", "PrimaryName C", etc.
        char companion_letter = 'B';
        for (auto* comp_ptr : companions) {
            if (primary->name && !primary->name->empty()) {
                comp_ptr->name = std::make_shared<std::string>(*primary->name + " " + companion_letter);
                companion_letter++;
            } else {
                comp_ptr->name = nullptr;
            }
        }
    }

    // Count stars per spectral type for reporting
    std::map<std::string, long long> total_spectral_counts;
    std::map<int, std::map<std::string, long long>> cluster_spectral_counts;
    std::map<int, long long> stars_per_cluster_count;

    for (const auto& s : star_systems) {
        const std::string spectral_type_str = spectral_type_to_string(s.spectral_type);
        total_spectral_counts[spectral_type_str]++;
        if (s.cluster_id > 0) {
            cluster_spectral_counts[s.cluster_id][spectral_type_str]++;
            stars_per_cluster_count[s.cluster_id]++;
        }
    }

    // Display counts on console
    std::cout << "\n--- スペクトル型ごとの星の数 (全体) ---" << std::endl;
    // spectral_types_vec を使って、表示順を固定する
    for (const auto& type : spectral_types_vec) {
        std::cout << type << ": " << total_spectral_counts[type] << " 個" << std::endl;
    }

    std::cout << "\n--- 星団ごとの詳細 ---" << std::endl;
    for(const auto& pair : stars_per_cluster_count) {
        int cluster_id = pair.first;
        long long total_stars = pair.second;
        std::cout << "星団ID: " << cluster_id << " (合計 " << total_stars << " 個)" << std::endl;
        if (cluster_spectral_counts.count(cluster_id)) {
            // spectral_types_vec を使って、表示順を固定する
            for (const auto& type : spectral_types_vec) {
                if (cluster_spectral_counts[cluster_id].count(type)) {
                    std::cout << "  " << type << ": " << cluster_spectral_counts[cluster_id][type] << " 個" << std::endl;
                }
            }
        }
    }

    // 9. Output to CSV file
    const std::string output_filename = "star_systems.csv";
    std::ofstream ofs(output_filename);
    if (!ofs.is_open()) {
        std::cerr << "エラー: " << output_filename << " を書き込み用に開けませんでした。" << std::endl;
        return 1;
    }

    // UTF-8 BOM (Byte Order Mark) を書き込み、ファイルがUTF-8として正しく解釈されるようにする
    ofs << "\xEF\xBB\xBF";

    // Write spectral counts summary to the top of the CSV
    ofs << "# スペクトル型ごとの星の数 (全体)\n";
    // spectral_types_vec を使って、書き込み順を固定する
    for (const auto& type : spectral_types_vec) {
        ofs << "# " << type << "," << total_spectral_counts[type] << "\n";
    }
    ofs << "#\n";

    // Write per-cluster details to the CSV
    ofs << "# 星団ごとの詳細\n";
    for(const auto& pair : stars_per_cluster_count) {
        int cluster_id = pair.first;
        long long total_stars = pair.second;
        ofs << "# 星団ID," << cluster_id << ",合計," << total_stars << "\n";
        if (cluster_spectral_counts.count(cluster_id)) {
            for (const auto& type : spectral_types_vec) {
                if (cluster_spectral_counts[cluster_id].count(type)) {
                    ofs << "#  " << type << "," << cluster_spectral_counts[cluster_id][type] << "\n";
                }
            }
        }
    }
    ofs << "#\n";

    // Write cluster 3D properties to the CSV for re-drawing
    ofs << "# 星団の物理的プロパティ (再描画用)\n";
    ofs << "# CLUSTER_PROP,id,center_x,center_y,center_z,axis_a,axis_b,axis_c,rot_alpha,rot_beta,rot_gamma,is_old\n";
    ofs << std::fixed << std::setprecision(6); // Use higher precision for properties
    for (const auto& prop : cluster_properties) {
        ofs << "# CLUSTER_PROP,"
            << prop.id << ","
            << prop.center_x << ","
            << prop.center_y << ","
            << prop.center_z << ","
            << prop.axis_a << ","
            << prop.axis_b << ","
            << prop.axis_c << ","
            << prop.rot_alpha << ","
            << prop.rot_beta << ","
            << prop.rot_gamma << ","
            << prop.is_old << "\n";
    }
    ofs << "#\n";

    // CSV Header (Japanese as requested)
    ofs << "通し番号,連星系ID,視等級,スペクトル,絶対等級,距離r,経度θ,緯度φ,x,y,z,星団ID,星系絶対等級,星系視等級,星の名前\n";

    // --- 高速化のためのバッファリング書き込み ---
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4); // Set precision for floating-point numbers
    const size_t buffer_flush_size = 1000; // この行数ごとにファイルに書き出す
    size_t line_count = 0;

    for (const auto& s : star_systems) {
        ss << s.id << ","
           << s.binary_system_id << ","
           << s.apparent_magnitude << ","
           << spectral_type_to_string(s.spectral_type) << ","
           << s.absolute_magnitude << ","
           << s.distance_r << ","
           << s.longitude_theta << ","
           << s.latitude_phi << ","
           << s.x << ","
           << s.y << ","
           << s.z << ","
           << s.cluster_id << ","
           << s.system_absolute_magnitude << ","
           << s.system_apparent_magnitude << ","
           << (s.name ? *s.name : "") << "\n";
        line_count++;

        if (line_count >= buffer_flush_size) {
            ofs << ss.rdbuf(); // バッファの内容をファイルに書き出す
            ss.str("");      // バッファをクリア
            ss.clear();      // ストリームの状態をリセット
            line_count = 0;
        }
    }
    // ループ終了後、バッファに残っているデータを書き出す
    if (line_count > 0) {
        ofs << ss.rdbuf();
    }

    ofs.close();
    std::cout << "星系データは " << output_filename << " に正常に書き込まれました。" << std::endl;

    // 3Dビューア用の軽量CSVを書き出す
    output_3d_viewer_csv(star_systems, "3d_viewer/3d_viewer_data.csv");

    // 10. Generate Star Maps
    std::cout << "\n--- 星図の生成 ---" << std::endl;
    // ステレオ投影（平射図法）の星図を生成
    generate_star_map_svg(star_systems, "star_map_stereographic_north.svg", cluster_properties, ProjectionType::Stereographic, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);
    // 正距方位図法の星図を生成
    generate_star_map_svg(star_systems, "star_map_azimuthal_equidistant_north.svg", cluster_properties, ProjectionType::AzimuthalEquidistant, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);

    std::cout << "\n--- シミュレーション完了 ---\n";
    return 0;
}

// 描画パラメータをユーザーから取得するヘルパー関数の実装
void get_drawing_parameters_from_user(
    bool& show_celestial_lines,
    float& max_magnitude_for_svg,
    bool& show_cluster_regions,
    bool& show_star_labels,
    float& label_magnitude_threshold
) {
    std::string line;

    // 赤道・黄道の描画有無をユーザーに確認
    std::cout << "星図に赤道と黄道を描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_celestial_lines = false;
    } else {
        show_celestial_lines = true;
    }

    // SVGに描画する最大等級をユーザーに確認
    std::cout << "星図に描画する星の最大等級を入力してください (デフォルト: 4.0f): ";
    std::getline(std::cin, line);
    if (!line.empty()) {
        try {
            max_magnitude_for_svg = std::stof(line);
        } catch (const std::exception&) {
            std::cout << "無効な入力です。デフォルト値 " << max_magnitude_for_svg << " を使用します。" << std::endl;
        }
    }

    // 星団領域の描画有無をユーザーに確認
    std::cout << "星図に星団の領域を描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_cluster_regions = false;
    } else {
        show_cluster_regions = true;
    }

    // 星のラベル描画有無をユーザーに確認
    std::cout << "星図に明るい星の名前ラベルを描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_star_labels = false;
    } else {
        show_star_labels = true;
    }

    // ラベルを描画する等級のしきい値をユーザーに確認
    if (show_star_labels) {
        std::cout << "名前ラベルを描画する星の最大等級を入力してください (デフォルト: 2.5): ";
        std::getline(std::cin, line);
        if (!line.empty()) {
            try {
                label_magnitude_threshold = std::stof(line);
            } catch (const std::exception&) {
                std::cout << "無効な入力です。デフォルト値 " << label_magnitude_threshold << " を使用します。" << std::endl;
            }
        }
    }
}

// CSVからデータを読み込んで星図を再描画する関数
void redraw_from_csv() {
    std::cout << "\n--- CSVから星図を再生成 ---" << std::endl;
    std::cout << "読み込むCSVファイル名を入力してください (デフォルト: star_systems.csv): ";
    std::string filename;
    std::getline(std::cin, filename);
    if (filename.empty()) {
        filename = "star_systems.csv";
    }

    // 先に描画パラメータをユーザーから取得
    std::cout << "\n新しい描画パラメータを入力してください:" << std::endl;
    bool show_celestial_lines = true;
    float max_magnitude_for_svg = 4.0f;
    bool show_cluster_regions = true;
    bool show_star_labels = true;
    float label_magnitude_threshold = 2.5f;
    get_drawing_parameters_from_user(
        show_celestial_lines,
        max_magnitude_for_svg,
        show_cluster_regions,
        show_star_labels,
        label_magnitude_threshold
    );

    // 指定された等級でフィルタリングしながらCSVからデータをロード
    std::cout << "\n" << filename << " から " << max_magnitude_for_svg << " 等級より明るい星を読み込んでいます..." << std::endl;
    auto loaded_data = load_data_from_csv(filename, max_magnitude_for_svg);
    if (!loaded_data) {
        // エラーメッセージはload_data_from_csv内で表示される
        std::cout << "メニューに戻ります。" << std::endl;
        return;
    }

    // 読み込んだデータを展開 (コピーを避けるためconst参照を使用)
    const auto& star_systems = loaded_data->first;
    const auto& cluster_properties = loaded_data->second;

    // 新しい設定で星図を再生成 (ファイル名を変えて上書きを防ぐ)
    std::cout << "\n--- 星図の再生成 ---" << std::endl;
    generate_star_map_svg(star_systems, "redraw_star_map_stereographic_north.svg", cluster_properties, ProjectionType::Stereographic, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);
    generate_star_map_svg(star_systems, "redraw_star_map_azimuthal_equidistant_north.svg", cluster_properties, ProjectionType::AzimuthalEquidistant, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);

    std::cout << "\n--- 再描画完了 ---\n";
}

// --- インタラクティブモード用ヘルパー関数 ---

void display_main_menu() {
    std::cout << "\n========================================" << std::endl;
    std::cout << " Star System Generator" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "[1] 新しい星空をランダムに生成" << std::endl;
    std::cout << "[2] シード値を指定して星空を生成" << std::endl;
    std::cout << "[3] CSVファイルから星図を再生成" << std::endl;
    std::cout << "[4] 3Dビューアの起動方法を表示" << std::endl;
    std::cout << "[0] 終了" << std::endl;
    std::cout << "選択してください: ";
}

void handle_random_generation() {
    std::random_device rd;
    unsigned int seed = rd();
    if (run_simulation(seed) != 0) {
        std::cout << "エラーが発生しました。メニューに戻ります。" << std::endl;
    }
}

void handle_seeded_generation() {
    std::cout << "使用するシード値を入力してください: ";
    std::string seed_str;
    std::getline(std::cin, seed_str);
    try {
        unsigned long long seed_val = std::stoull(seed_str); // 1. より大きな型で受け取る
        if (seed_val > std::numeric_limits<unsigned int>::max()) { // 2. 範囲チェック
            std::cout << "エラー: シード値が大きすぎます。unsigned intの最大値 ("
                      << std::numeric_limits<unsigned int>::max() << ") 以下で入力してください。" << std::endl;
            return; // メニューに戻る
        }
        unsigned int seed = static_cast<unsigned int>(seed_val); // 3. 安全にキャスト
        if (run_simulation(seed) != 0) {
            std::cout << "エラーが発生しました。メニューに戻ります。" << std::endl;
        }
    } catch (const std::exception& e) {
        // stoullが失敗した場合 (数字でない、またはullの範囲外)
        std::cout << "エラー: 無効なシード値です。数字で入力してください。" << std::endl;
    }
}

void show_3d_viewer_instructions() {
    std::cout << "\n--- 3Dビューアの起動方法 ---" << std::endl;
    std::cout << "生成された star_systems.csv ファイルをインタラクティブな3Dビューアで表示できます。" << std::endl;
    std::cout << "手順:" << std::endl;
    std::cout << "1. Python3.xがPCにインストールされていることを確認してください。" << std::endl;
    std::cout << "2. このプログラムの実行ファイルがあるディレクトリでターミナル（コマンドプロンプトやPowerShell）を開きます。" << std::endl;
    std::cout << "3. 以下のコマンドを実行して、簡易Webサーバーを起動します:" << std::endl;
    std::cout << "   py -m http.server" << std::endl;
    std::cout << "   (もし 'python' コマンドで動作しない場合は 'python -m http.server' を試してください)" << std::endl;
    std::cout << "4. Webブラウザを開き、アドレスバーに以下のURLを入力してアクセスします:" << std::endl;
    std::cout << "   http://localhost:8000/3d_viewer/" << std::endl;
    std::cout << "5. サーバーを停止するには、ターミナルで Ctrl+C を押してください。" << std::endl;
    std::cout << "\n--- !! セキュリティに関する注意 !! ---" << std::endl;
    std::cout << "この簡易サーバーは開発・確認用です。認証機能はなく、同じネットワーク上の誰でもアクセスできる可能性があります。" << std::endl;
    std::cout << "・信頼できるネットワーク（自宅のWi-Fiなど）でのみ使用してください。" << std::endl;
    std::cout << "・公共のWi-Fiなど、信頼できないネットワークでは絶対に使用しないでください。" << std::endl;
    std::cout << "・使い終わったら、必ずターミナルで Ctrl+C を押してサーバーを停止してください。" << std::endl;
    std::cout << "\n何かキーを押してメニューに戻ってください...";
    std::string dummy;
    std::getline(std::cin, dummy);
}

int main(int argc, char* argv[]) {
    #ifdef _WIN32
        // WindowsのコンソールでUTF-8を正しく表示するための設定
        SetConsoleOutputCP(CP_UTF8);
    #endif

    // コマンドライン引数でシードが指定された場合は、一度だけ実行して終了する（従来通りの動作）
    if (argc > 1) {
        try {
            unsigned long long seed_val = std::stoull(argv[1]);
            if (seed_val > std::numeric_limits<unsigned int>::max()) {
                std::cerr << "エラー: シード値が大きすぎます。unsigned intの最大値 ("
                          << std::numeric_limits<unsigned int>::max() << ") 以下で指定してください。" << std::endl;
                return 1;
            }
            unsigned int seed = static_cast<unsigned int>(seed_val);
            std::cout << "コマンドラインで指定されたシード値で実行します: " << seed << std::endl;
            if (run_simulation(seed) != 0) {
                std::cerr << "シミュレーション中にエラーが発生しました。" << std::endl;
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "エラー: 無効なシード値です。数字で指定してください。 (" << e.what() << ")" << std::endl;
            return 1;
        }
        return 0;
    }

    // インタラクティブモード
    while (true) {
        display_main_menu();

        std::string choice_str;
        std::getline(std::cin, choice_str);

        if (choice_str == "1") {
            handle_random_generation();
        } else if (choice_str == "2") {
            handle_seeded_generation();
        } else if (choice_str == "3") {
            redraw_from_csv();
        } else if (choice_str == "4") {
            show_3d_viewer_instructions();
        } else if (choice_str == "0" || choice_str.empty()) { // Enterキーだけでも終了
            std::cout << "プログラムを終了します。" << std::endl;
            break;
        } else {
            std::cout << "無効な選択です。もう一度入力してください。" << std::endl;
        }
    }

    return 0;
}