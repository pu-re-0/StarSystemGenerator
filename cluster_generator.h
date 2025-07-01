#ifndef CLUSTER_GENERATOR_H
#define CLUSTER_GENERATOR_H

#include <vector>
#include <string>
#include <random>
#include <map>

// 両方のファイルで共有する星系の構造体
struct StarSystem {
    long long id;
    long long binary_system_id; // 連星系ID
    double apparent_magnitude;
    std::string spectral_type;
    double absolute_magnitude;
    double distance_r;
    double longitude_theta; // 経度 (radians)
    double latitude_phi;    // 緯度 (radians)
    double x, y, z;         // デカルト座標 (Cartesian coordinates)
    int cluster_id;         // 0: 現場星, >0: 星団ID
    double system_absolute_magnitude; // 星系全体の絶対等級
    double system_apparent_magnitude; // 星系全体の視等級
};

// 星団の物理的特性を保持する構造体
struct ClusterProperties {
    int id;
    double center_x, center_y, center_z;
    double major_axis_a; // 長軸
    double minor_axis_b; // 短軸
    double rot_alpha;    // Z軸周りの回転
    double rot_beta;     // Y軸周りの回転
    double rot_gamma;    // X軸周りの回転
    bool is_old;         // 古い星団かどうか
};

// 複数の等級を合算して一つの等級にするヘルパー関数
double combine_magnitudes(const std::vector<double>& magnitudes);

// 主星と伴星のリストから星系全体の等級を計算し、各星に設定する関数
void set_system_magnitudes(StarSystem& primary, std::vector<StarSystem>& companions);

// 星の核となる特性（スペクトル型、絶対等級）を生成する共通関数
void generate_star_core_properties(
    StarSystem& s,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    std::discrete_distribution<>& spectral_dist,
    const std::map<std::string, std::pair<double, double>>& abs_mag_ranges
);

// 連星間距離の分布パラメータ
struct BinarySeparationParams {
    // 近接連星用
    double close_log_mean;
    double close_log_stddev;
    // 遠距離連星用
    double wide_log_mean;
    double wide_log_stddev;
    // 近接連星になる確率
    double close_binary_probability;
};

// 星生成に必要なパラメータをまとめる構造体
struct StarGenerationParameters {
    const std::vector<std::string>& spectral_types_vec;
    const std::map<std::string, std::pair<double, double>>& abs_mag_ranges;
    std::discrete_distribution<>& spectral_dist; // 分布は状態を持つので非const参照
    const std::map<std::string, BinarySeparationParams>& binary_separation_params;
};

// 連星系を生成するための共通関数
std::vector<StarSystem> generate_companions_for_primary(
    const StarSystem& primary,
    const std::map<std::string, double>& multiplicity_map,
    std::mt19937& gen,
    const StarGenerationParameters& params
);

// 散開星団を生成する関数の宣言
std::pair<std::vector<StarSystem>, std::vector<ClusterProperties>> generate_open_clusters(
    long long total_num_starts,
    double R_parsecs,
    std::mt19937& gen,
    const std::vector<std::string>& spectral_types_vec,
    const std::vector<double>& cluster_spectral_weights,
    const std::map<std::string, std::pair<double, double>>& abs_mag_ranges,
    const std::map<std::string, double>& multiplicity_map,
    const std::map<std::string, struct BinarySeparationParams>& binary_separation_params,
    long long& next_binary_system_id
);

#endif // CLUSTER_GENERATOR_H