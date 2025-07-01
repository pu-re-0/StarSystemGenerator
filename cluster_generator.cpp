#include "cluster_generator.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm> // For std::shuffle, std::clamp
#include <optional>  // For std::optional
#include <limits>

// 複数の等級を合算して一つの等級にするヘルパー関数
// 等級は対数スケールのため、一度光度（フラックス）に変換し、合算してから再度等級に戻す
double combine_magnitudes(const std::vector<double>& magnitudes) {
    if (magnitudes.empty()) {
        // 空の場合は非常に暗い等級を返す（エラー値として）
        return 99.0;
    }
    if (magnitudes.size() == 1) {
        return magnitudes[0];
    }

    double total_flux = 0.0;
    for (double mag : magnitudes) {
        // 等級(m)から相対的な光度(F)への変換: F = 10^(-0.4 * m)
        total_flux += std::pow(10.0, -0.4 * mag);
    }

    // 合計光度(F_total)から合成等級(m_combined)への変換: m = -2.5 * log10(F)
    // total_fluxが0以下になることは理論上ないが、念のためチェック
    if (total_flux <= 0) return 99.0;
    return -2.5 * std::log10(total_flux);
}

// 主星と伴星のリストから星系全体の等級を計算し、各星に設定する関数
void set_system_magnitudes(StarSystem& primary, std::vector<StarSystem>& companions) {
    if (companions.empty()) {
        // 単独星の場合、星系全体の等級は主星（単独星）の等級と同じ
        primary.system_absolute_magnitude = primary.absolute_magnitude;
        primary.system_apparent_magnitude = primary.apparent_magnitude;
    } else {
        // 連星系の場合
        std::vector<double> abs_mags_to_combine;
        std::vector<double> app_mags_to_combine;

        // 主星の等級を追加
        abs_mags_to_combine.push_back(primary.absolute_magnitude);
        app_mags_to_combine.push_back(primary.apparent_magnitude);

        // 伴星の等級を追加
        for (const auto& comp : companions) {
            abs_mags_to_combine.push_back(comp.absolute_magnitude);
            app_mags_to_combine.push_back(comp.apparent_magnitude);
        }

        // 主星に合算等級を設定
        primary.system_absolute_magnitude = combine_magnitudes(abs_mags_to_combine);
        primary.system_apparent_magnitude = combine_magnitudes(app_mags_to_combine);

        // 伴星に-1を設定
        for (auto& comp : companions) {
            comp.system_absolute_magnitude = -1.0;
            comp.system_apparent_magnitude = -1.0;
        }
    }
}

// 星の核となる特性（スペクトル型、絶対等級）を生成する共通関数
void generate_star_core_properties(
    StarSystem& s,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    std::discrete_distribution<>& spectral_dist,
    const std::map<std::string, std::pair<double, double>>& abs_mag_ranges
) {
    // スペクトル型を生成
    int spectral_idx = spectral_dist(gen);
    s.spectral_type = spectral_types_vec[spectral_idx];

    // スペクトル型に基づいて絶対等級を生成
    auto range_it = abs_mag_ranges.find(s.spectral_type);
    if (range_it != abs_mag_ranges.end()) {
        double min_mag = range_it->second.first;
        double max_mag = range_it->second.second;

        if (s.spectral_type == "ES") {
            std::lognormal_distribution<double> log_dist(0, 0.55);
            double log_val;
            do { log_val = log_dist(gen); } while (log_val >= max_mag - min_mag);
            s.absolute_magnitude = max_mag - log_val;
        } else if (s.spectral_type == "B") {
            std::lognormal_distribution<double> log_dist(0.3, 0.75);
            double log_val;
            do { log_val = log_dist(gen); } while (log_val >= std::abs(max_mag - min_mag));
            s.absolute_magnitude = max_mag - log_val;
        } else if (s.spectral_type == "A") {
            std::lognormal_distribution<double> log_dist(0.75, 1.0);
            double log_val;
            do { log_val = log_dist(gen); } while (log_val >= std::abs(max_mag - min_mag));
            s.absolute_magnitude = max_mag - log_val;
        } else {
            std::uniform_real_distribution<double> dist_abs_mag(min_mag, max_mag);
            s.absolute_magnitude = dist_abs_mag(gen);
        }
    } else {
        // このケースは現在のロジックでは到達不能ですが、安全策として記述します。
        s.absolute_magnitude = 5.0; // デフォルト値（太陽相当）
        std::cerr << "警告: スペクトル型 " << s.spectral_type << " の絶対等級範囲が見つかりません。" << std::endl;
    }
}

// 天文単位 (AU) からパーセクへの変換定数
const double AU_TO_PARSEC = 4.848e-6; // 1 AU in parsecs

// 伴星を生成するヘルパー関数 (このファイル内でのみ使用)
static std::optional<std::pair<StarSystem, double>> generate_companion( // `std::optional` を使用
    const StarSystem& primary,
    const StarGenerationParameters& params,
    std::mt19937& gen,
    double min_separation_au = 0.0,
    const std::optional<std::vector<std::string>>& spectral_type_whitelist = std::nullopt
) {
    StarSystem companion;
    const double PI = std::acos(-1.0);

    companion.cluster_id = primary.cluster_id;
    companion.binary_system_id = primary.binary_system_id; // 主星の連星系IDを継承

    // --- Step 1: 伴星の核となる特性（スペクトル型と絶対等級）の決定 ---
    // 使用するスペクトル型と分布を決定する
    const std::vector<std::string>* types_to_use = &params.spectral_types_vec;
    std::discrete_distribution<>* dist_to_use = &params.spectral_dist;

    // ホワイトリストが指定されている場合のカスタム分布準備
    std::vector<std::string> temp_allowed_types;
    std::vector<double> temp_allowed_weights;
    std::discrete_distribution<> temp_whitelist_dist;

    if (spectral_type_whitelist && !spectral_type_whitelist->empty()) {
        auto original_probs = params.spectral_dist.probabilities();
        for (size_t i = 0; i < params.spectral_types_vec.size(); ++i) {
            const auto& type = params.spectral_types_vec[i];
            // ホワイトリストに型が含まれているかチェック
            if (std::find(spectral_type_whitelist->begin(), spectral_type_whitelist->end(), type) != spectral_type_whitelist->end()) {
                temp_allowed_types.push_back(type);
                temp_allowed_weights.push_back(original_probs[i]);
            }
        }

        if (!temp_allowed_types.empty()) {
            // 新しい重みを正規化して合計を1.0にする
            double total_weight = std::accumulate(temp_allowed_weights.begin(), temp_allowed_weights.end(), 0.0);
            if (total_weight > std::numeric_limits<double>::epsilon()) {
                for (double& weight : temp_allowed_weights) {
                    weight /= total_weight;
                }
            }
            temp_whitelist_dist = std::discrete_distribution<>(temp_allowed_weights.begin(), temp_allowed_weights.end());
            
            types_to_use = &temp_allowed_types;
            dist_to_use = &temp_whitelist_dist;
        } else {
             std::cerr << "警告: 指定されたスペクトル型ホワイトリストに一致する型がありません。通常の分布を使用します。" << std::endl;
        }
    }

    // 決定された分布を使用して星を生成
    int attempts = 0;
    do {
        generate_star_core_properties(companion, gen, *types_to_use, *dist_to_use, params.abs_mag_ranges);
        attempts++;
    } while (companion.absolute_magnitude < primary.absolute_magnitude && attempts < 50);

    if (companion.absolute_magnitude < primary.absolute_magnitude) {
        companion.absolute_magnitude = primary.absolute_magnitude + std::uniform_real_distribution<double>(0.0, 0.1)(gen);
    }

    // --- Step 2: 主星からの距離（公転軌道の大きさ）の決定 ---
    // 主星のスペクトル型に基づいて分布パラメータを取得
    BinarySeparationParams sep_params;
    auto it = params.binary_separation_params.find(primary.spectral_type);
    if (it != params.binary_separation_params.end()) {
        sep_params = it->second;
    } else {
        // フォールバック用のデフォルト値 (G型星相当の2ピーク分布)
        sep_params = {std::log(3.0), 1.0, std::log(50.0), 1.2, 0.7}; 
    }

    // 混合モデルからのサンプリングは、まずどちらの分布からサンプリングするかを確率的に決定する。
    // 近接連星か遠距離連星かを確率で決定
    std::bernoulli_distribution choose_close_binary(sep_params.close_binary_probability);

    std::lognormal_distribution<double> separation_dist;
    if (choose_close_binary(gen)) {
        // 近接連星の分布を使用
        separation_dist = std::lognormal_distribution<double>(sep_params.close_log_mean, sep_params.close_log_stddev);
    } else {
        // 遠距離連星の分布を使用
        separation_dist = std::lognormal_distribution<double>(sep_params.wide_log_mean, sep_params.wide_log_stddev);
    }

    double separation_au;
    int sep_attempts = 0;
    do {
        separation_au = separation_dist(gen);
        sep_attempts++;
    } while (separation_au < min_separation_au && sep_attempts < 50);

    if (separation_au < min_separation_au) {
        separation_au = min_separation_au * std::uniform_real_distribution<double>(1.05, 1.2)(gen);
    }

    double separation_pc = separation_au * AU_TO_PARSEC;

    std::uniform_real_distribution<double> dist_u(0.0, 1.0);
    double u = dist_u(gen);
    double v = dist_u(gen);
    double theta = 2 * PI * u;
    double phi = std::acos(2 * v - 1);

    double dx = separation_pc * std::sin(phi) * std::cos(theta);
    double dy = separation_pc * std::sin(phi) * std::sin(theta);
    double dz = separation_pc * std::cos(phi);

    companion.x = primary.x + dx;
    companion.y = primary.y + dy;
    companion.z = primary.z + dz;

    companion.distance_r = std::sqrt(companion.x * companion.x + companion.y * companion.y + companion.z * companion.z);
    if (companion.distance_r < 0.001) {
        companion.distance_r = 0.001;
    }
    companion.latitude_phi = std::asin(std::clamp(companion.z / companion.distance_r, -1.0, 1.0));
    companion.longitude_theta = std::atan2(companion.y, companion.x);
    companion.apparent_magnitude = companion.absolute_magnitude + 5.0 * (std::log10(companion.distance_r) - 1.0);

    return std::make_pair(companion, separation_au);
}

// 連星系を生成するための共通関数
std::vector<StarSystem> generate_companions_for_primary(
    const StarSystem& primary,
    const std::map<std::string, double>& multiplicity_map,
    std::mt19937& gen,
    const StarGenerationParameters& params
) {
    std::vector<StarSystem> companions_to_add;
    double multiplicity = multiplicity_map.count(primary.spectral_type) ? multiplicity_map.at(primary.spectral_type) : 0.0;

    if (multiplicity > 0) {
        std::uniform_real_distribution<double> dist_binary_check(0.0, 1.0);
        double prob_comp1 = std::min(multiplicity, 0.8);

        if (dist_binary_check(gen) < prob_comp1) {
            // --- 主星の型に基づいて、伴星のスペクトル型ホワイトリストを作成 ---
            std::optional<std::vector<std::string>> companion_whitelist;
            const std::vector<std::string> spectral_hierarchy = {"O", "B", "A", "F", "G", "K", "M"};
            
            auto it = std::find(spectral_hierarchy.begin(), spectral_hierarchy.end(), primary.spectral_type);

            if (it != spectral_hierarchy.end()) {
                // 主系列星の場合: 自分と同じかそれ以下の質量の主系列星、および白色矮星を許可
                std::vector<std::string> allowed_types(it, spectral_hierarchy.end());
                allowed_types.push_back("WD");
                companion_whitelist = allowed_types;
            } else if (primary.spectral_type == "WD") {
                // 主星が白色矮星の場合: 低質量の主系列星か、別の白色矮星を許可
                companion_whitelist = std::vector<std::string>{"K", "M", "WD"};
            } else if (primary.spectral_type == "ES") {
                // 主星が進化した星の場合: 中～低質量の主系列星か、白色矮星を許可
                companion_whitelist = std::vector<std::string>{"F", "G", "K", "M", "WD"};
            }
            // 上記のいずれでもない場合 (例: 未知の型)、ホワイトリストは std::nullopt のままとなり、全ての型が許可される。

            auto companion1_opt = generate_companion(primary, params, gen, 0.0, companion_whitelist);
            if (companion1_opt) {
                companions_to_add.push_back(companion1_opt->first);
                double separation1_au = companion1_opt->second;

                double prob_comp2 = std::max(0.0, multiplicity - 0.8);
                const std::string& st = primary.spectral_type;
                if (st == "F" || st == "G" || st == "K" || st == "A" || st == "B" || st == "O") {
                    prob_comp2 = std::max(prob_comp2, 0.10);
                }

                if (dist_binary_check(gen) < prob_comp2) {
                    // 2番目の伴星にも同じホワイトリストを適用
                    auto companion2_opt = generate_companion(primary, params, gen, separation1_au * 2.5, companion_whitelist);
                    if (companion2_opt) {
                        companions_to_add.push_back(companion2_opt->first);
                    }
                }
            }
        }
    }
    return companions_to_add;
}

// 散開星団を生成する関数の実装
std::pair<std::vector<StarSystem>, std::vector<ClusterProperties>> generate_open_clusters(
    long long total_num_stars,
    double R_parsecs,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    const std::vector<double>& cluster_spectral_weights,
    const std::map<std::string, std::pair<double, double>>& abs_mag_ranges,
    const std::map<std::string, double>& multiplicity_map,
    const std::map<std::string, BinarySeparationParams>& binary_separation_params,
    long long& next_binary_system_id
) {
    const double PI = std::acos(-1.0);
    std::vector<StarSystem> cluster_stars_vec;
    std::vector<ClusterProperties> cluster_properties_vec;

    // この関数内で使用する乱数分布
    std::uniform_real_distribution<double> dist_uniform_0_1(0.0, 1.0);
    std::uniform_real_distribution<double> dist_cluster_star_fraction(0.005, 0.01);
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> dist_phi_sin(-1.0, 1.0);

    // --- 星団の数と所属する星の数を決定 ---
    // (星の総数1,300個/星団に含まれる割合)あたり1個の星団を基準に、±30%の範囲で数を決定する
    double base_num_clusters = static_cast<double>(total_num_stars) * dist_cluster_star_fraction(gen) / 1300.0;
    int min_num_clusters = static_cast<int>(std::round(base_num_clusters * 0.7));
    int max_num_clusters = static_cast<int>(std::round(base_num_clusters * 1.3));

    // 最小値が最大値を上回らないように調整
    if (max_num_clusters < min_num_clusters) {
        max_num_clusters = min_num_clusters;
    }

    int num_clusters = 0;
    if (max_num_clusters > 0) {
        std::uniform_int_distribution<int> dist_num_clusters(min_num_clusters, max_num_clusters);
        num_clusters = dist_num_clusters(gen);
    }

    if (num_clusters == 0) {
        std::cout << "星の総数が少ないため、散開星団は生成されませんでした。" << std::endl;
        return {cluster_stars_vec, cluster_properties_vec}; // 空のペアを返す
    }

    // --- B型・O型を持たない「古い」星団を決定 ---
    // 50% ± 20% (30%～70%)の星団を「古い」星団とする
    std::uniform_real_distribution<double> dist_old_fraction(0.3, 0.7);
    int num_old_clusters = static_cast<int>(std::round(num_clusters * dist_old_fraction(gen)));

    // どの星団が「古い」かランダムに割り当てる
    std::vector<int> cluster_indices(num_clusters);
    std::iota(cluster_indices.begin(), cluster_indices.end(), 0); // 0, 1, 2, ... と連番で埋める
    std::shuffle(cluster_indices.begin(), cluster_indices.end(), gen);

    std::vector<bool> is_old_cluster(num_clusters, false);
    for (int i = 0; i < num_old_clusters; ++i) {
        is_old_cluster[cluster_indices[i]] = true;
    }

    // --- 「古い」星団用のスペクトル分布を作成 ---
    std::vector<double> old_cluster_weights = cluster_spectral_weights;
    // B型とO型の重みを0にする
    for (size_t i = 0; i < spectral_types_vec.size(); ++i) {
        if (spectral_types_vec[i] == "B" || spectral_types_vec[i] == "O") {
            old_cluster_weights[i] = 0.0;
        }
        else if (spectral_types_vec[i] == "WD") {
            old_cluster_weights[i] = 0.05;
        }
    }
    // 新しい重みを正規化して合計を1.0にする
    double total_old_weight = std::accumulate(old_cluster_weights.begin(), old_cluster_weights.end(), 0.0);
    if (total_old_weight > std::numeric_limits<double>::epsilon()) {
        for (double& weight : old_cluster_weights) {
            weight /= total_old_weight;
        }
    }

    // 2種類のスペクトル分布を用意
    std::discrete_distribution<> young_cluster_spectral_dist(cluster_spectral_weights.begin(), cluster_spectral_weights.end());
    std::discrete_distribution<> old_cluster_spectral_dist(old_cluster_weights.begin(), old_cluster_weights.end());

    const long long total_cluster_stars = static_cast<long long>(total_num_stars * dist_cluster_star_fraction(gen));
    long long stars_in_clusters_generated = 0;

    std::cout << num_clusters << " 個の散開星団を生成し、合計 " << total_cluster_stars << " 個の星を配置します。" << std::endl;

    // --- 各星団の星の数を配分 ---
    std::vector<double> cluster_star_dist_weights;
    for (int i = 0; i < num_clusters; ++i) {
        cluster_star_dist_weights.push_back(pow(dist_uniform_0_1(gen), 1.5));
    }
    double total_weight = std::accumulate(cluster_star_dist_weights.begin(), cluster_star_dist_weights.end(), 0.0);

    std::vector<long long> stars_per_cluster(num_clusters);
    if (total_weight > 0) {
        for (int i = 0; i < num_clusters - 1; ++i) {
            stars_per_cluster[i] = static_cast<long long>(total_cluster_stars * (cluster_star_dist_weights[i] / total_weight));
            stars_in_clusters_generated += stars_per_cluster[i];
        }
        stars_per_cluster[num_clusters - 1] = total_cluster_stars - stars_in_clusters_generated;
    } else if (num_clusters > 0) {
        // Fallback: もし全ての重みが0なら均等に配分
        long long stars_each = total_cluster_stars / num_clusters;
        for(int i=0; i<num_clusters; ++i) stars_per_cluster[i] = stars_each;
        stars_per_cluster[num_clusters-1] += total_cluster_stars % num_clusters;
    }

    // --- 各星団を生成し、星を配置 ---
    for (int i = 0; i < num_clusters; ++i) {
        long long current_cluster_star_count = stars_per_cluster[i];
        if (current_cluster_star_count <= 0) continue;

        // 星団のプロパティを定義
        std::uniform_real_distribution<double> dist_major_axis(1.5, 4.5);
        double major_axis_a = dist_major_axis(gen);
        std::uniform_real_distribution<double> dist_aspect_ratio(0.4, 0.6);
        double minor_axis_b = major_axis_a * dist_aspect_ratio(gen);

        // 回転楕円体のランダムな向きを決定 (オイラー角)
        double rot_alpha = dist_theta(gen); // Z軸周り
        double rot_beta = dist_theta(gen);  // Y軸周り
        double rot_gamma = dist_theta(gen); // X軸周り

        // 星団全体がシミュレーション空間に収まるように中心位置を決定
        double safe_radius = R_parsecs - major_axis_a;
        if (safe_radius < 0) safe_radius = 0;
        std::uniform_real_distribution<double> dist_center_pos(-safe_radius, safe_radius);
        double center_x = dist_center_pos(gen);
        double center_y = dist_center_pos(gen);
        double center_z = dist_center_pos(gen);

        // この星団のプロパティを保存
        cluster_properties_vec.push_back({
            i + 1,
            center_x, center_y, center_z,
            major_axis_a,
            minor_axis_b,
            rot_alpha,
            rot_beta,
            rot_gamma,
            is_old_cluster[i] // 星団の年齢フラグを設定
        });

        for (long long j = 0; j < current_cluster_star_count; ++j) {
            StarSystem s;
            s.cluster_id = i + 1; // 星団IDを設定
            s.binary_system_id = next_binary_system_id++; // 新しい連星系IDを割り当て

            // --- Step 1: 共通関数を使って星の特性（種類と明るさ）を決定 ---
            auto& current_spectral_dist = is_old_cluster[i] ? old_cluster_spectral_dist : young_cluster_spectral_dist;
            generate_star_core_properties(s, gen, spectral_types_vec, current_spectral_dist, abs_mag_ranges);

            // --- 連星生成ロジック ---
            StarGenerationParameters cluster_params = {spectral_types_vec, abs_mag_ranges, current_spectral_dist, binary_separation_params};
            std::vector<StarSystem> companions_to_add = generate_companions_for_primary(
                s,
                multiplicity_map,
                gen,
                cluster_params
            );

            // --- Step 2: 星の特性に基づいて中心への集中度を決定 ---
            double concentration_exponent = 1.5; // 通常の星の集中度
            // 明るく重い星（O, B, A型）はより中心に集める
            if (s.spectral_type == "O" || s.spectral_type == "B" || s.spectral_type == "A") {
                concentration_exponent = 3.0;
            }

            // --- Step 3: 決定した集中度で星の位置を生成 ---
            // 中心が濃くなるように、単位球内に点を生成
            double r_norm = std::pow(dist_uniform_0_1(gen), concentration_exponent);
            double theta_norm = dist_theta(gen);
            double phi_norm_sin = dist_phi_sin(gen);
            double phi_norm_cos = std::sqrt(1.0 - phi_norm_sin * phi_norm_sin);

            double x_unit = r_norm * phi_norm_cos * std::cos(theta_norm);
            double y_unit = r_norm * phi_norm_cos * std::sin(theta_norm);
            double z_unit = r_norm * phi_norm_sin;

            // 回転楕円体の形状にスケーリング
            double x_scaled = minor_axis_b * x_unit;
            double y_scaled = minor_axis_b * y_unit;
            double z_scaled = major_axis_a * z_unit;

            // オイラー角でランダムに回転
            double x1 = x_scaled * std::cos(rot_alpha) - y_scaled * std::sin(rot_alpha);
            double y1 = x_scaled * std::sin(rot_alpha) + y_scaled * std::cos(rot_alpha);
            double z1 = z_scaled;
            double x2 = x1 * std::cos(rot_beta) + z1 * std::sin(rot_beta);
            double y2 = y1;
            double z2 = -x1 * std::sin(rot_beta) + z1 * std::cos(rot_beta);
            s.x = x2;
            s.y = y2 * std::cos(rot_gamma) - z2 * std::sin(rot_gamma);
            s.z = y2 * std::sin(rot_gamma) + z2 * std::cos(rot_gamma);

            // 星団の中心位置へ移動
            s.x += center_x;
            s.y += center_y;
            s.z += center_z;

            // --- Step 4: 残りのプロパティを計算 ---
            s.distance_r = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
            // 距離が0に近すぎる場合、log10(0)などの計算エラーを防ぐために最小値を設定
            if (s.distance_r < 0.001) {
                s.distance_r = 0.001;
            }
            // 球面座標を計算 (distance_r > 0 が保証されているため安全)
            s.latitude_phi = std::asin(s.z / s.distance_r);
            s.longitude_theta = std::atan2(s.y, s.x);

            // 視等級
            s.apparent_magnitude = s.absolute_magnitude + 5.0 * (std::log10(s.distance_r) - 1.0);

            // 星系全体の等級を計算して設定
            set_system_magnitudes(s, companions_to_add);

            // 主星をリストに追加
            cluster_stars_vec.push_back(s);

            // 生成された伴星をリストに追加し、ループカウンタを進める
            for (const auto& comp : companions_to_add) {
                // ループの上限を超えないように最終チェック
                if (j < current_cluster_star_count - 1) {
                    cluster_stars_vec.push_back(comp);
                    j++; // 伴星を生成したので、この星団で生成した星の総数を1増やす
                }
            }
        }
    }
    return {cluster_stars_vec, cluster_properties_vec};
}