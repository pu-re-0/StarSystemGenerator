#ifndef CLUSTER_GENERATOR_H
#define CLUSTER_GENERATOR_H

#include <vector>
#include <string>
#include <random>
#include <map>
#include <memory> // For std::shared_ptr
#include <cstdint> // For uint8_t

// スペクトル型を効率的に管理するためのenum class
enum class SpectralType : uint8_t {
    M, K, G, F, A, B, O, WD, ES, UNKNOWN
};

// Enumと文字列を相互変換するヘルパー関数
std::string spectral_type_to_string(SpectralType type);
SpectralType string_to_spectral_type(const std::string& s);

// 両方のファイルで共有する星系の構造体
struct StarSystem {
    long long id;
    long long binary_system_id; // 連星系ID
    float apparent_magnitude;
    SpectralType spectral_type; // メモリ効率の良いenumに変更
    float absolute_magnitude;
    float distance_r;
    float longitude_theta; // 経度 (radians)
    float latitude_phi;    // 緯度 (radians)
    float x, y, z;         // デカルト座標 (Cartesian coordinates)
    int cluster_id;         // 0: 現場星, >0: 星団ID
    float system_absolute_magnitude; // 星系全体の絶対等級
    float system_apparent_magnitude; // 星系全体の視等級
    std::shared_ptr<std::string> name; // 名前を持つ星だけメモリを確保するように変更
};

// 星団の物理的特性を保持する構造体
struct ClusterProperties {
    int id;
    float center_x, center_y, center_z;
    float axis_a;       // ローカルX軸方向の半径
    float axis_b;       // ローカルY軸方向の半径
    float axis_c;       // ローカルZ軸方向の半径
    float rot_alpha;    // Z軸周りの回転 (Yaw)
    float rot_beta;     // Y軸周りの回転
    float rot_gamma;    // X軸周りの回転
    bool is_old;         // 古い星団かどうか
};

// 複数の等級を合算して一つの等級にするヘルパー関数
float combine_magnitudes(const std::vector<float>& magnitudes);

// 主星と伴星のリストから星系全体の等級を計算し、各星に設定する関数
void set_system_magnitudes(StarSystem& primary, std::vector<StarSystem>& companions);

// 星の核となる特性（スペクトル型、絶対等級）を生成する共通関数
void generate_star_core_properties(
    StarSystem& s,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    std::discrete_distribution<>& spectral_dist,
    const std::map<std::string, std::pair<float, float>>& abs_mag_ranges
);

// 連星間距離の分布パラメータ
struct BinarySeparationParams {
    // 近接連星用
    float close_log_mean;
    float close_log_stddev;
    // 遠距離連星用
    float wide_log_mean;
    float wide_log_stddev;
    // 近接連星になる確率
    float close_binary_probability;
};

// 星生成に必要なパラメータをまとめる構造体
struct StarGenerationParameters {
    const std::vector<std::string>& spectral_types_vec;
    const std::map<std::string, std::pair<float, float>>& abs_mag_ranges;
    std::discrete_distribution<>& spectral_dist; // 分布は状態を持つので非const参照
    const std::map<std::string, BinarySeparationParams>& binary_separation_params;
};

// 連星系を生成するための共通関数
std::vector<StarSystem> generate_companions_for_primary(
    const StarSystem& primary,
    const std::map<std::string, float>& multiplicity_map,
    std::mt19937& gen,
    const StarGenerationParameters& params
);

// 散開星団を生成する関数の宣言
std::pair<std::vector<StarSystem>, std::vector<ClusterProperties>> generate_open_clusters(
    long long total_num_stars,
    float R_parsecs,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    const std::vector<float>& cluster_spectral_weights,
    const std::map<std::string, std::pair<float, float>>& abs_mag_ranges,
    const std::map<std::string, float>& multiplicity_map,
    const std::map<std::string, struct BinarySeparationParams>& binary_separation_params,
    long long& next_binary_system_id
);

#endif // CLUSTER_GENERATOR_H