#include "cluster_generator.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm> // For std::shuffle, std::clamp
#include <optional>  // For std::optional
#include <limits>
#include <unordered_set> // std::unordered_set を使用するために追加
#include <memory>    // For std::shared_ptr

// Enumと文字列を相互変換するヘルパー関数の実装
std::string spectral_type_to_string(SpectralType type) {
    switch (type) {
        case SpectralType::O: return "O";
        case SpectralType::B: return "B";
        case SpectralType::A: return "A";
        case SpectralType::F: return "F";
        case SpectralType::G: return "G";
        case SpectralType::K: return "K";
        case SpectralType::M: return "M";
        case SpectralType::WD: return "WD";
        case SpectralType::ES: return "ES";
        default: return "UNKNOWN";
    }
}

SpectralType string_to_spectral_type(const std::string& s) {
    if (s == "O") return SpectralType::O;
    if (s == "B") return SpectralType::B;
    if (s == "A") return SpectralType::A;
    if (s == "F") return SpectralType::F;
    if (s == "G") return SpectralType::G;
    if (s == "K") return SpectralType::K;
    if (s == "M") return SpectralType::M;
    if (s == "WD") return SpectralType::WD;
    if (s == "ES") return SpectralType::ES;
    return SpectralType::UNKNOWN;
}

// 複数の等級を合算して一つの等級にするヘルパー関数
// 等級は対数スケールのため、一度光度（フラックス）に変換し、合算してから再度等級に戻す
float combine_magnitudes(const std::vector<float>& magnitudes) {
    if (magnitudes.empty()) {
        // 空の場合は非常に暗い等級を返す（エラー値として）
        return 99.0f;
    }
    if (magnitudes.size() == 1) {
        return magnitudes[0];
    }

    float total_flux = 0.0f;
    for (float mag : magnitudes) {
        // 等級(m)から相対的な光度(F)への変換: F = 10^(-0.4 * m)
        total_flux += std::powf(10.0f, -0.4f * mag);
    }

    // 合計光度(F_total)から合成等級(m_combined)への変換: m = -2.5 * log10(F)
    // total_fluxが0以下になることは理論上ないが、念のためチェック
    if (total_flux <= 0.0f) return 99.0f;
    return -2.5f * std::log10f(total_flux);
}

// 主星と伴星のリストから星系全体の等級を計算し、各星に設定する関数
void set_system_magnitudes(StarSystem& primary, std::vector<StarSystem>& companions) {
    if (companions.empty()) {
        // 単独星の場合、星系全体の等級は主星（単独星）の等級と同じ
        primary.system_absolute_magnitude = primary.absolute_magnitude;
        primary.system_apparent_magnitude = primary.apparent_magnitude;
    } else {
        // 連星系の場合
        std::vector<float> abs_mags_to_combine;
        std::vector<float> app_mags_to_combine;

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
            comp.system_absolute_magnitude = -1.0f;
            comp.system_apparent_magnitude = -1.0f;
        }
    }
}

// 星の核となる特性（スペクトル型、絶対等級）を生成する共通関数
void generate_star_core_properties(
    StarSystem& s,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    std::discrete_distribution<>& spectral_dist,
    const std::map<std::string, std::pair<float, float>>& abs_mag_ranges
) {
    // スペクトル型を生成
    int spectral_idx = spectral_dist(gen);
    const std::string& spectral_type_str = spectral_types_vec[spectral_idx];
    s.spectral_type = string_to_spectral_type(spectral_type_str);

    // スペクトル型に基づいて絶対等級を生成
    auto range_it = abs_mag_ranges.find(spectral_type_str);
    if (range_it != abs_mag_ranges.end()) {
        float min_mag = range_it->second.first;
        float max_mag = range_it->second.second;

        if (spectral_type_str == "ES") {
            std::lognormal_distribution<float> log_dist(0.3f, 0.75f);
            float log_val;
            do { log_val = log_dist(gen); } while (log_val >= std::abs(max_mag - min_mag));
            s.absolute_magnitude = max_mag - log_val;
        } else if (s.spectral_type == SpectralType::A) {
            std::lognormal_distribution<float> log_dist(0.75f, 1.0f);
            float log_val;
            do { log_val = log_dist(gen); } while (log_val >= std::abs(max_mag - min_mag));
            s.absolute_magnitude = max_mag - log_val;
        } else {
            std::uniform_real_distribution<float> dist_abs_mag(min_mag, max_mag);
            s.absolute_magnitude = dist_abs_mag(gen);
        }
    } else {
        // このケースは現在のロジックでは到達不能ですが、安全策として記述します。
        s.absolute_magnitude = 5.0f; // デフォルト値（太陽相当）
        std::cerr << "警告: スペクトル型 " << spectral_type_to_string(s.spectral_type) << " の絶対等級範囲が見つかりません。" << std::endl;
    }
}

// 天文単位 (AU) からパーセクへの変換定数
const float AU_TO_PARSEC = 4.848e-6f; // 1 AU in parsecs

// 伴星を生成するヘルパー関数 (このファイル内でのみ使用)
static std::optional<std::pair<StarSystem, float>> generate_companion( // `std::optional` を使用
    const StarSystem& primary,
    const StarGenerationParameters& params,
    std::mt19937& gen,
    float min_separation_au = 0.0f,
    const std::optional<std::vector<std::string>>& spectral_type_whitelist = std::nullopt
) {
    StarSystem companion;
    const float PI = std::acos(-1.0f);

    companion.cluster_id = primary.cluster_id;
    companion.binary_system_id = primary.binary_system_id; // 主星の連星系IDを継承

    // --- Step 1: 伴星の核となる特性（スペクトル型と絶対等級）の決定 ---
    // 使用するスペクトル型と分布を決定する
    const std::vector<std::string>* types_to_use = &params.spectral_types_vec;
    std::discrete_distribution<>* dist_to_use = &params.spectral_dist;

    // ホワイトリストが指定されている場合のカスタム分布準備
    std::vector<std::string> temp_allowed_types;
    std::vector<float> temp_allowed_weights;
    std::discrete_distribution<> temp_whitelist_dist;

    if (spectral_type_whitelist && !spectral_type_whitelist->empty()) {
        // 高速なルックアップのためにホワイトリストをunordered_setに変換
        const std::unordered_set<std::string> whitelist_set(spectral_type_whitelist->begin(), spectral_type_whitelist->end());
        auto original_probs = params.spectral_dist.probabilities();
        for (size_t i = 0; i < params.spectral_types_vec.size(); ++i) {
            const auto& type = params.spectral_types_vec[i];
            // 高速化されたハッシュセットで検索
            if (whitelist_set.count(type)) {
                temp_allowed_types.push_back(type);
                temp_allowed_weights.push_back(original_probs[i]);
            }
        }

        if (!temp_allowed_types.empty()) {
            temp_whitelist_dist = std::discrete_distribution<>(temp_allowed_weights.begin(), temp_allowed_weights.end());
            types_to_use = &temp_allowed_types;
            dist_to_use = &temp_whitelist_dist;
        } else {
            std::cerr << "警告: 指定されたスペクトル型ホワイトリストに一致する型がありません。通常の分布を使用します。" << std::endl;
        }
    }

    // 決定された分布を使用して星を生成
    // リジェクションサンプリング（do-whileループ）を廃止し、一度だけ生成する
    generate_star_core_properties(companion, gen, *types_to_use, *dist_to_use, params.abs_mag_ranges);

    // もし伴星が主星より明るくなってしまった場合、物理的に妥当な範囲に調整（クランプ）する
    // これにより、無限ループや多数の試行を防ぎ、パフォーマンスを安定させる
    if (companion.absolute_magnitude < primary.absolute_magnitude) {
        // 主星の等級にわずかな値を加えることで、伴星が少しだけ暗くなるように調整
        companion.absolute_magnitude = primary.absolute_magnitude + std::uniform_real_distribution<float>(0.0f, 0.1f)(gen);
    }

    // --- Step 2: 主星からの距離（公転軌道の大きさ）の決定 ---
    // 主星のスペクトル型に基づいて分布パラメータを取得
    BinarySeparationParams sep_params;
    auto it = params.binary_separation_params.find(spectral_type_to_string(primary.spectral_type));
    if (it != params.binary_separation_params.end()) {
        sep_params = it->second;
    } else {
        // フォールバック用のデフォルト値 (G型星相当の2ピーク分布)
        sep_params = {std::logf(3.0f), 1.0f, std::logf(50.0f), 1.2f, 0.7f}; 
    }

    // 混合モデルからのサンプリングは、まずどちらの分布からサンプリングするかを確率的に決定する。
    // 近接連星か遠距離連星かを確率で決定
    std::bernoulli_distribution choose_close_binary(sep_params.close_binary_probability);

    std::lognormal_distribution<float> separation_dist;
    if (choose_close_binary(gen)) {
        // 近接連星の分布を使用
        separation_dist = std::lognormal_distribution<float>(sep_params.close_log_mean, sep_params.close_log_stddev);
    } else {
        // 遠距離連星の分布を使用
        separation_dist = std::lognormal_distribution<float>(sep_params.wide_log_mean, sep_params.wide_log_stddev);
    }

    float separation_au;
    // リジェクションサンプリングを廃止
    separation_au = separation_dist(gen);

    // 階層的連星系で、内側の軌道より小さくならないように調整（クランプ）
    if (separation_au < min_separation_au) {
        // 最小距離にわずかなマージンを加えて設定
        separation_au = min_separation_au * std::uniform_real_distribution<float>(1.05f, 1.2f)(gen);
    }

    float separation_pc = separation_au * AU_TO_PARSEC;

    std::uniform_real_distribution<float> dist_u(0.0f, 1.0f);
    float u = dist_u(gen);
    float v = dist_u(gen);
    float theta = 2.0f * PI * u;
    float phi = std::acosf(2.0f * v - 1.0f);

    float dx = separation_pc * std::sinf(phi) * std::cosf(theta);
    float dy = separation_pc * std::sinf(phi) * std::sinf(theta);
    float dz = separation_pc * std::cosf(phi);

    companion.x = primary.x + dx;
    companion.y = primary.y + dy;
    companion.z = primary.z + dz;

    companion.distance_r = std::sqrtf(companion.x * companion.x + companion.y * companion.y + companion.z * companion.z);
    if (companion.distance_r < 0.001f) {
        companion.distance_r = 0.001f;
    }
    companion.latitude_phi = std::asinf(std::clamp(companion.z / companion.distance_r, -1.0f, 1.0f));
    companion.longitude_theta = std::atan2f(companion.y, companion.x);
    companion.apparent_magnitude = companion.absolute_magnitude + 5.0f * (std::log10f(companion.distance_r) - 1.0f);

    return std::make_pair(companion, separation_au);
}

// 連星系を生成するための共通関数
std::vector<StarSystem> generate_companions_for_primary(
    const StarSystem& primary,
    const std::map<std::string, float>& multiplicity_map,
    std::mt19937& gen,
    const StarGenerationParameters& params
) {
    std::vector<StarSystem> companions_to_add;
    const std::string primary_spectral_type_str = spectral_type_to_string(primary.spectral_type);

    float multiplicity = multiplicity_map.count(primary_spectral_type_str) ? multiplicity_map.at(primary_spectral_type_str) : 0.0f;

    if (multiplicity > 0) {
        std::uniform_real_distribution<float> dist_binary_check(0.0f, 1.0f);
        float prob_comp1 = std::min(multiplicity, 0.8f);

        if (dist_binary_check(gen) < prob_comp1) {
            // --- 主星の型に基づいて、伴星のスペクトル型ホワイトリストを作成 ---
            std::optional<std::vector<std::string>> companion_whitelist;
            const std::vector<std::string> spectral_hierarchy = {"O", "B", "A", "F", "G", "K", "M"};
            
            auto it = std::find(spectral_hierarchy.begin(), spectral_hierarchy.end(), primary_spectral_type_str);

            if (it != spectral_hierarchy.end()) {
                // 主系列星の場合: 自分と同じかそれ以下の質量の主系列星、および白色矮星を許可
                std::vector<std::string> allowed_types(it, spectral_hierarchy.end());
                allowed_types.push_back("WD");
                companion_whitelist = allowed_types;
            } else if (primary.spectral_type == SpectralType::WD) {
                // 主星が白色矮星の場合: 低質量の主系列星か、別の白色矮星を許可
                companion_whitelist = std::vector<std::string>{"K", "M", "WD"};
            } else if (primary.spectral_type == SpectralType::ES) {
                // 主星が進化した星の場合: 中～低質量の主系列星か、白色矮星を許可
                companion_whitelist = std::vector<std::string>{"F", "G", "K", "M", "WD"};
            }
            // 上記のいずれでもない場合 (例: 未知の型)、ホワイトリストは std::nullopt のままとなり、全ての型が許可される。

            auto companion1_opt = generate_companion(primary, params, gen, 0.0f, companion_whitelist);
            if (companion1_opt) {
                companions_to_add.push_back(companion1_opt->first);
                float separation1_au = companion1_opt->second;

                float prob_comp2 = std::max(0.0f, multiplicity - 1.0f);
                const std::string& st = primary_spectral_type_str;
                if (st == "F" || st == "G" || st == "K" || st == "A" || st == "B" || st == "O") {
                    prob_comp2 = std::max(prob_comp2, 0.10f);
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
    float R_parsecs,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    const std::vector<float>& cluster_spectral_weights,
    const std::map<std::string, std::pair<float, float>>& abs_mag_ranges,
    const std::map<std::string, float>& multiplicity_map,
    const std::map<std::string, BinarySeparationParams>& binary_separation_params,
    long long& next_binary_system_id
) {
    const float PI = std::acosf(-1.0f);
    std::vector<StarSystem> cluster_stars_vec;
    std::vector<ClusterProperties> cluster_properties_vec;

    // 生成される星の総数をおおよそで見積もり、メモリを事前に確保する。
    // これにより、ループ内での push_back によるメモリ再確保のオーバーヘッドを大幅に削減する。
    // 星団の星の割合の最大値(0.01)と、伴星の生成(多めに見積もって1.5倍)を考慮する。
    const long long total_cluster_stars_estimate = static_cast<long long>(total_num_stars * 0.01f * 1.5f);
    cluster_stars_vec.reserve(static_cast<size_t>(total_cluster_stars_estimate));

    // この関数内で使用する乱数分布
    std::uniform_real_distribution<float> dist_uniform_0_1(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_cluster_star_fraction(0.005f, 0.01f);
    std::uniform_real_distribution<float> dist_theta(0.0f, 2.0f * PI);
    std::uniform_real_distribution<float> dist_phi_sin(-1.0f, 1.0f);
    std::normal_distribution<float> dist_normal(0.0f, 1.0f); // 正規分布 N(0, 1)

    // --- 星団の数と所属する星の数を決定 ---
    // (星の総数1,300個/星団に含まれる割合)あたり1個の星団を基準に、±30%の範囲で数を決定する
    float base_num_clusters = static_cast<float>(total_num_stars) * dist_cluster_star_fraction(gen) / 1300.0f;
    int min_num_clusters = static_cast<int>(std::roundf(base_num_clusters * 0.7f));
    int max_num_clusters = static_cast<int>(std::roundf(base_num_clusters * 1.3f));

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
    std::uniform_real_distribution<float> dist_old_fraction(0.3f, 0.7f);
    int num_old_clusters = static_cast<int>(std::roundf(num_clusters * dist_old_fraction(gen)));

    // どの星団が「古い」かランダムに割り当てる
    std::vector<int> cluster_indices(num_clusters);
    std::iota(cluster_indices.begin(), cluster_indices.end(), 0); // 0, 1, 2, ... と連番で埋める
    std::shuffle(cluster_indices.begin(), cluster_indices.end(), gen);

    std::vector<bool> is_old_cluster(num_clusters, false);
    for (int i = 0; i < num_old_clusters; ++i) {
        is_old_cluster[cluster_indices[i]] = true;
    }

    // --- 「古い」星団用のスペクトル分布を作成 ---
    std::vector<float> old_cluster_weights = cluster_spectral_weights;
    // B型とO型の重みを0にする
    for (size_t i = 0; i < spectral_types_vec.size(); ++i) {
        if (spectral_types_vec[i] == "B" || spectral_types_vec[i] == "O") {
            old_cluster_weights[i] = 0.0f;
        }
        else if (spectral_types_vec[i] == "WD") {
            old_cluster_weights[i] = 0.05f;
        }
    }
    // 新しい重みを正規化して合計を1.0fにする
    float total_old_weight = std::accumulate(old_cluster_weights.begin(), old_cluster_weights.end(), 0.0f);
    if (total_old_weight > std::numeric_limits<float>::epsilon()) {
        for (float& weight : old_cluster_weights) {
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
    std::vector<float> cluster_star_dist_weights;
    for (int i = 0; i < num_clusters; ++i) {
        cluster_star_dist_weights.push_back(pow(dist_uniform_0_1(gen), 1.5));
    }
    float total_weight = std::accumulate(cluster_star_dist_weights.begin(), cluster_star_dist_weights.end(), 0.0f);

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
        // 三軸不等楕円体の3つの軸長を生成
        std::uniform_real_distribution<float> dist_axis(2.0f, 5.0f);
        float axis_a = dist_axis(gen); // 最も長い軸
        std::uniform_real_distribution<float> dist_aspect_ratio1(0.5f, 0.9f);
        float axis_b = axis_a * dist_aspect_ratio1(gen); // 2番目に長い軸
        std::uniform_real_distribution<float> dist_aspect_ratio2(0.4f, 0.8f);
        float axis_c = axis_b * dist_aspect_ratio2(gen); // 最も短い軸

        // 回転楕円体のランダムな向きを決定 (オイラー角)
        float rot_alpha = dist_theta(gen); // Z軸周り
        float rot_beta = dist_theta(gen);  // Y軸周り
        float rot_gamma = dist_theta(gen); // X軸周り

        // 楕円体の最大半径を計算
        float max_axis = std::max({axis_a, axis_b, axis_c});
        // 星団全体がシミュレーション空間に収まるように中心位置を決定する。
        // 星団は長半径 'major_axis_a' を持つ回転楕円体である。
        // 星団全体が半径 R_parsecs のシミュレーション空間内に収まるためには、
        // 星団の中心は、原点から半径 (R_parsecs - major_axis_a) の球内に配置される必要がある。
        float safe_radius = R_parsecs - max_axis;
        if (safe_radius < 0) safe_radius = 0;

        // 'safe_radius' の球内に一様にランダムな点を生成して、星団の中心位置とする。
        // (以前の方法では立方体内に点を生成していたため、中心が遠すぎる可能性があった)
        float center_r = safe_radius * std::cbrt(dist_uniform_0_1(gen)); // 体積で一様分布させるため三乗根を使用
        float center_theta = dist_theta(gen);
        float center_phi_sin = dist_phi_sin(gen);
        float center_phi_cos = std::sqrt(1.0f - center_phi_sin * center_phi_sin);
        float center_x = center_r * center_phi_cos * std::cos(center_theta);
        float center_y = center_r * center_phi_cos * std::sin(center_theta);
        float center_z = center_r * center_phi_sin;

        // この星団のプロパティを保存
        cluster_properties_vec.push_back({
            i + 1,
            static_cast<float>(center_x), static_cast<float>(center_y), static_cast<float>(center_z),
            static_cast<float>(axis_a),
            static_cast<float>(axis_b),
            static_cast<float>(axis_c),
            static_cast<float>(rot_alpha),
            static_cast<float>(rot_beta),
            static_cast<float>(rot_gamma),
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

            // --- Step 2: 星の特性に基づいて中心への集中度（正規分布の広がり）を決定 ---
            // sigma_scaleが小さいほど、星は中心に強く集まる。
            float sigma_scale = 1.0f; // 通常の星の分布の広がり
            // 明るく重い星（O, B, A型）はより中心に集める (分布を狭める)
            if (s.spectral_type == SpectralType::O || s.spectral_type == SpectralType::B || s.spectral_type == SpectralType::A) {
                sigma_scale = 0.5f; // 分布の広がりを半分にする
            }

            // --- Step 3: 3次元正規分布に従って星の位置を生成 ---
            // 3つの独立した標準正規分布 N(0,1) の乱数を生成
            float nx = dist_normal(gen);
            float ny = dist_normal(gen);
            float nz = dist_normal(gen);

            // --- Step 4: 回転楕円体の形状と質量分離に応じてスケーリング ---
            // 標準偏差を (軸長 / 3) に設定することで、星の大部分が楕円体内に収まるようにする。
            // ローカル座標系の各軸に沿って、対応する軸長でスケーリングする。
            float sigma_x = (axis_a / 3.0f) * sigma_scale;
            float sigma_y = (axis_b / 3.0f) * sigma_scale;
            float sigma_z = (axis_c / 3.0f) * sigma_scale;

            float x_scaled = nx * sigma_x;
            float y_scaled = ny * sigma_y;
            float z_scaled = nz * sigma_z; // ここはローカル座標系での位置

            // オイラー角でランダムに回転
            float x1 = x_scaled * std::cosf(rot_alpha) - y_scaled * std::sinf(rot_alpha);
            float y1 = x_scaled * std::sinf(rot_alpha) + y_scaled * std::cosf(rot_alpha);
            float z1 = z_scaled;
            float x2 = x1 * std::cosf(rot_beta) + z1 * std::sinf(rot_beta);
            float y2 = y1;
            float z2 = -x1 * std::sinf(rot_beta) + z1 * std::cosf(rot_beta);
            s.x = x2;
            s.y = y2 * std::cosf(rot_gamma) - z2 * std::sinf(rot_gamma);
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
            s.apparent_magnitude = s.absolute_magnitude + 5.0f * (std::log10(s.distance_r) - 1.0f);

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