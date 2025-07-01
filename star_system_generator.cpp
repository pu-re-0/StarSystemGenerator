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
#include <omp.h>     // OpenMPヘッダー
#include <optional>  // For std::optional


// ヘルパー関数：スペクトル型に基づいてSVGで使う色を返す
std::string get_star_color_for_plot(const std::string& spectral_type) {
    if (spectral_type == "O") return "#9bb0ff";  // 青白
    if (spectral_type == "B") return "#aabfff";  // 青白
    if (spectral_type == "A") return "#cad8ff";  // 白
    if (spectral_type == "F") return "#f8f7ff";  // 黄白
    if (spectral_type == "G") return "#fff4ea";  // 黄
    if (spectral_type == "K") return "#ffd2a1";  // 橙
    if (spectral_type == "M") return "#ffb56c";  // 赤
    if (spectral_type == "WD") return "#f0f0f0"; // 白
    if (spectral_type == "ES") return "#ffeea8"; // 明るい黄
    return "white";
}

// 投影法の種類を定義するenum
enum class ProjectionType {
    Stereographic,          // 平射図法（ステレオ投影）
    AzimuthalEquidistant    // 正距方位図法
};

// 指定された投影法で星図をSVGファイルとして生成する汎用関数
void generate_star_map_svg(const std::vector<StarSystem>& all_stars, const std::string& filename, const std::vector<ClusterProperties>& clusters, ProjectionType proj_type, bool show_celestial_lines = true, double max_magnitude_to_plot = 5.0, bool show_cluster_regions = true, bool show_star_labels = true, double label_magnitude_threshold = 2.5) {
    const double PI = std::acos(-1.0);
    const int svg_size = 1000;
    const double center = svg_size / 2.0;
    const double plot_scale = svg_size / 2.0;
    const double max_plot_radius = 8.0; // 最も明るい星の半径
    const double min_plot_radius = 1.0; // 6.0等星の半径

    // 投影法ごとの設定
    double limit_latitude_deg;
    std::string projection_name;

    if (proj_type == ProjectionType::Stereographic) {
        limit_latitude_deg = 0.0; // 北極から赤道まで
        projection_name = "Stereographic Projection (North Pole to Equator)";
    } else { // AzimuthalEquidistant
        limit_latitude_deg = -55.0; // 北極から南緯55度まで
        projection_name = "Azimuthal Equidistant Projection (North Pole to -55 deg)";
    }
    const double limit_latitude_rad = limit_latitude_deg * PI / 180.0;

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "エラー: SVGファイル " << filename << " を書き込み用に開けませんでした。" << std::endl;
        return;
    }

    // SVGヘッダー
    ofs << "<svg width=\"" << svg_size << "\" height=\"" << svg_size << "\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background-color:black;\">\n";
    ofs << "<!-- Generated Star Map (" << projection_name << ") -->\n";

    // クリッピングパスの定義
    ofs << "<defs>\n";
    ofs << "  <clipPath id=\"map-circle\">\n";
    ofs << "    <circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << plot_scale << "\"/>\n";
    ofs << "  </clipPath>\n";
    ofs << "</defs>\n";

    // 描画範囲の外周円
    ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << plot_scale << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"1\"/>\n";

    // 緯度線の描画
    if (proj_type == ProjectionType::Stereographic) {
        // 赤道が外周になる。緯度30度、60度線を描画
        double lat60_colatitude = PI / 2.0 - (60.0 * PI / 180.0);
        double r60 = plot_scale * std::tan(lat60_colatitude / 2.0);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << r60 << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
        double lat30_colatitude = PI / 2.0 - (30.0 * PI / 180.0);
        double r30 = plot_scale * std::tan(lat30_colatitude / 2.0);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << r30 << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
    } else { // AzimuthalEquidistant
        // 南緯55度が外周になる。赤道と南緯30度線を描画
        const double max_angular_distance = (PI / 2.0) - limit_latitude_rad;
        double equator_angular_distance = PI / 2.0 - (0.0 * PI / 180.0);
        double equator_radius = plot_scale * (equator_angular_distance / max_angular_distance);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << equator_radius << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
        double lat30s_angular_distance = PI / 2.0 - (-30.0 * PI / 180.0);
        double lat30s_radius = plot_scale * (lat30s_angular_distance / max_angular_distance);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << lat30s_radius << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
    }

    // 座標投影ヘルパー (ラムダ式)
    auto project_point = [&](double lon, double lat) -> std::pair<double, double> {
        double px, py;
        // 中心からの角距離 c = 90度 - 緯度
        double colatitude_c = (PI / 2.0) - lat;

        if (proj_type == ProjectionType::Stereographic) {
            double r_plot = plot_scale * std::tan(colatitude_c / 2.0);
            px = center + r_plot * std::cos(lon);
            py = center - r_plot * std::sin(lon);
        } else { // AzimuthalEquidistant
            const double max_angular_distance = (PI / 2.0) - limit_latitude_rad;
            double r_plot = plot_scale * (colatitude_c / max_angular_distance);
            px = center + r_plot * std::cos(lon);
            py = center - r_plot * std::sin(lon);
        }
        return {px, py};
    };

    // クリッピングを適用するグループを開始
    ofs << "<g clip-path=\"url(#map-circle)\">\n";

    if (show_celestial_lines) {
        // --- 天の赤道の描画 (赤色) ---
        ofs << "  <polyline points=\"";
        for (int i = 0; i <= 360; ++i) {
            auto p = project_point(static_cast<double>(i) * PI / 180.0, 0.0);
            ofs << p.first << "," << p.second << " ";
        }
        ofs << "\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1.5\"/>\n";

        // --- 黄道の描画 (黄色) ---
        const double obliquity_rad = 23.44 * PI / 180.0; // 黄道傾斜角
        ofs << "  <polyline points=\"";
        for (int i = 0; i <= 360; ++i) {
            double ecliptic_lon = static_cast<double>(i) * PI / 180.0;
            double eq_lon = std::atan2(std::sin(ecliptic_lon) * std::cos(obliquity_rad), std::cos(ecliptic_lon));
            double eq_lat = std::asin(std::sin(obliquity_rad) * std::sin(ecliptic_lon));
            auto p = project_point(eq_lon, eq_lat);
            ofs << p.first << "," << p.second << " ";
        }
        ofs << "\" fill=\"none\" stroke=\"#ffff00\" stroke-width=\"1.5\" stroke-dasharray=\"5 5\"/>\n";
    }

    if (show_cluster_regions) {
        // --- 星団の領域を描画 ---
        ofs << "  <!-- Cluster Regions -->\n";

        // 星団の年齢に応じた色を定義
        const std::string young_cluster_fill = "rgba(100, 100, 255, 0.1)";
        const std::string young_cluster_stroke = "rgba(150, 150, 255, 0.5)";
        const std::string old_cluster_fill = "rgba(255, 200, 100, 0.1)";
        const std::string old_cluster_stroke = "rgba(255, 220, 150, 0.5)";

        for (const auto& cluster : clusters) {
            const std::string& fill_color = cluster.is_old ? old_cluster_fill : young_cluster_fill;
            const std::string& stroke_color = cluster.is_old ? old_cluster_stroke : young_cluster_stroke;
            ofs << "  <polygon points=\"";
            const int num_points = 36; // 輪郭を描くための点の数
            bool first_point = true;
            for (int i = 0; i <= num_points; ++i) {
                double angle = 2.0 * PI * i / num_points;

                // 楕円体の赤道上の点をローカル座標で生成
                double x_local = cluster.minor_axis_b * std::cos(angle);
                double y_local = cluster.minor_axis_b * std::sin(angle);
                double z_local = 0;

                // オイラー角で回転 (Z-Y-X順)
                // Z軸周り (alpha)
                double x1 = x_local * std::cos(cluster.rot_alpha) - y_local * std::sin(cluster.rot_alpha);
                double y1 = x_local * std::sin(cluster.rot_alpha) + y_local * std::cos(cluster.rot_alpha);
                double z1 = z_local;
                // Y軸周り (beta)
                double x2 = x1 * std::cos(cluster.rot_beta) + z1 * std::sin(cluster.rot_beta);
                double y2 = y1;
                double z2 = -x1 * std::sin(cluster.rot_beta) + z1 * std::cos(cluster.rot_beta);
                // X軸周り (gamma)
                double x_rot = x2;
                double y_rot = y2 * std::cos(cluster.rot_gamma) - z2 * std::sin(cluster.rot_gamma);
                double z_rot = y2 * std::sin(cluster.rot_gamma) + z2 * std::cos(cluster.rot_gamma);

                // 星団の中心位置へ移動
                double world_x = x_rot + cluster.center_x;
                double world_y = y_rot + cluster.center_y;
                double world_z = z_rot + cluster.center_z;

                // ワールド座標から球面座標へ変換
                double dist_r = std::sqrt(world_x * world_x + world_y * world_y + world_z * world_z);
                if (dist_r < 1e-6) continue; // 中心点はスキップ

                double lat = std::asin(std::clamp(world_z / dist_r, -1.0, 1.0));
                double lon = std::atan2(world_y, world_x);

                // 2Dに投影
                auto p = project_point(lon, lat);

                if (!first_point) {
                    ofs << " ";
                }
                ofs << p.first << "," << p.second;
                first_point = false;
            }
            // ポリゴンを閉じる
            ofs << "\" fill=\"" << fill_color << "\" stroke=\"" << stroke_color << "\" stroke-width=\"1\"/>\n";
        }
        ofs << "  <!-- End Cluster Regions -->\n";
    }

    // 1. 描画対象の星をフィルタリングし、その中での最小等級（最も明るい）を見つける
    std::vector<StarSystem> drawable_stars;
    double min_plot_magnitude = max_magnitude_to_plot;

    for (const auto& star : all_stars) {
        // 緯度でフィルタリング
        if (star.latitude_phi < limit_latitude_rad) {
            continue;
        }

        // 伴星 (-1) は描画しない
        if (star.system_apparent_magnitude == -1.0) {
            continue;
        }

        // 描画に使う等級は、星系全体の等級
        double magnitude_to_check = star.system_apparent_magnitude;

        if (magnitude_to_check <= max_magnitude_to_plot) {
            drawable_stars.push_back(star);
            if (magnitude_to_check < min_plot_magnitude) {
                min_plot_magnitude = magnitude_to_check;
            }
        }
    }

    if (drawable_stars.empty()) {
        std::cout << "警告: 描画対象の星（" << max_magnitude_to_plot << "等級より明るく、緯度" << limit_latitude_deg << "度より北）が見つかりませんでした。" << std::endl;
        ofs << "</svg>\n";
        ofs.close();
        return;
    }

    // 2. フィルタリングされた星を描画
    double magnitude_range = max_magnitude_to_plot - min_plot_magnitude;

    for (const auto& star : drawable_stars) {
        auto p = project_point(star.longitude_theta, star.latitude_phi);

        // 描画に使う等級は、星系全体の等級
        double magnitude_for_size = star.system_apparent_magnitude;

        // 等級に基づいて星の描画半径を計算
        double plot_radius;
        if (magnitude_range < 1e-6) { // 等級の範囲がほぼ0の場合（星が1つ、または全て同じ明るさ）
            plot_radius = max_plot_radius;
        } else {
            plot_radius = min_plot_radius + (max_plot_radius - min_plot_radius) * (max_magnitude_to_plot - magnitude_for_size) / magnitude_range;
        }

        ofs << "  <circle cx=\"" << p.first << "\" cy=\"" << p.second << "\" r=\"" << plot_radius << "\" fill=\"" << get_star_color_for_plot(star.spectral_type) << "\"/>\n";
    }

    // 3. 星のラベルを描画
    if (show_star_labels) {
        ofs << "  <!-- Star Labels -->\n";
        for (const auto& star : drawable_stars) {
            // ラベル表示対象の等級を決定
            double magnitude_for_label = star.system_apparent_magnitude;

            if (magnitude_for_label <= label_magnitude_threshold) {
                auto p = project_point(star.longitude_theta, star.latitude_phi);
                double text_x = p.first + 5;  // x座標を少し右にずらす
                double text_y = p.second + 5; // y座標を少し下にずらす

                // ラベルがクリッピング範囲外にはみ出さないように簡易チェック
                if (text_x > svg_size - 30) text_x = p.first - 10;
                if (text_y > svg_size - 5) text_y = p.second - 5;

                // 連星系IDをラベルとして表示
                ofs << "  <text x=\"" << text_x << "\" y=\"" << text_y
                    << "\" font-family=\"sans-serif\" font-size=\"8\" fill=\"#aaa\">"
                    << star.binary_system_id << "</text>\n";
            }
        }
        ofs << "  <!-- End Star Labels -->\n";

        // --- 星団のラベルを描画 ---
        ofs << "  <!-- Cluster Labels -->\n";
        for (const auto& cluster : clusters) {
            // 星団の中心座標をワールド座標から球面座標へ変換
            double dist_r = std::sqrt(cluster.center_x * cluster.center_x + cluster.center_y * cluster.center_y + cluster.center_z * cluster.center_z);
            if (dist_r < 1e-6) continue;

            double lat = std::asin(std::clamp(cluster.center_z / dist_r, -1.0, 1.0));
            double lon = std::atan2(cluster.center_y, cluster.center_x);

            // 描画範囲外の星団ラベルはスキップ
            if (lat < limit_latitude_rad) {
                continue;
            }

            auto p = project_point(lon, lat);

            ofs << "  <text x=\"" << p.first << "\" y=\"" << p.second
                << "\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"bold\" fill=\"#87ceeb\" text-anchor=\"middle\" dy=\".3em\">"
                << "C" << cluster.id << "</text>\n";
        }
        ofs << "  <!-- End Star Labels -->\n";
    }

    // グループを閉じる
    ofs << "</g>\n";
    
    ofs << "</svg>\n";
    ofs.close();
    std::cout << "星図が " << filename << " に正常に書き込まれました。" << std::endl;
}

int main() {
    // 1. Setup random number generator
    // Using std::random_device for non-deterministic seeding, or a fixed seed for reproducibility
    std::random_device rd;
    //std::mt19937 gen(rd()); // Seed with a random device 
    std::mt19937 gen(2025062809); // For reproducible results, use a fixed seed like 42  数値はペルド・ペルソーネより

    // Helper function to normalize a vector of weights to sum to 1.0
    auto normalize_weights = [](std::vector<double>& weights) {
        double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
        if (total_weight > std::numeric_limits<double>::epsilon()) {
            for (double& weight : weights) {
                weight /= total_weight;
            }
        }
    };

    // 2. Define spectral type properties (probabilities and absolute magnitude ranges)
    std::vector<std::string> spectral_types_vec = {"M", "K", "G", "F", "A", "B", "O", "WD", "ES"};
    std::vector<double> spectral_type_weights = {0.725, 0.129, 0.041, 0.029, 0.006, 0.00039, 0.0000005, 0.059, 0.011}; // Approximate weights, will be normalized
    normalize_weights(spectral_type_weights);

    // Absolute magnitude ranges for each spectral type {min_M, max_M}
    std::map<std::string, std::pair<double, double>> abs_mag_ranges;
    abs_mag_ranges["O"] = {-6.0, -4.0};
    abs_mag_ranges["B"] = {-4.5, 2.0};
    abs_mag_ranges["A"] = {1.5, 4.0};
    abs_mag_ranges["F"] = {3.5, 5.5};
    abs_mag_ranges["G"] = {4.5, 6.2};
    abs_mag_ranges["K"] = {5.8, 7.6};
    abs_mag_ranges["M"] = {7.4, 12.0};
    abs_mag_ranges["WD"] = {10.0, 16.0}; // White Dwarf
    abs_mag_ranges["ES"] = {-8.0, 3.0}; // non-Red Giant Evolved Stars

    // Multiplicity (平均伴星数) for each spectral type
    std::map<std::string, double> multiplicity_map;
    multiplicity_map["O"] = 1.3;  // O型星は多重連星系を形成することが多い
    multiplicity_map["B"] = 1.0;
    multiplicity_map["A"] = 1.0;
    multiplicity_map["F"] = 0.62;
    multiplicity_map["G"] = 0.62; // 太陽のようなG型星
    multiplicity_map["K"] = 0.62;
    multiplicity_map["M"] = 0.33; // M型（赤色矮星）は連星率が低い
    multiplicity_map["WD"] = 0.1; // 白色矮星
    multiplicity_map["ES"] = 0.3; // 進化した星

    // --- 連星間距離のパラメータ設定 ---
    // 2つのピーク（近接連星と遠距離連星）を持つ分布を定義
    // { close_log_mean, close_log_stddev, wide_log_mean, wide_log_stddev, close_binary_prob }
    std::map<std::string, BinarySeparationParams> binary_separation_params;
    // O型: 近い伴星も遠い伴星も持ちうる (中央値: 8 AU, 200 AU)
    binary_separation_params["O"]  = {std::log(8.0), 1.2, std::log(200.0), 1.5, 0.5};
    // B型 (中央値: 6 AU, 150 AU)
    binary_separation_params["B"]  = {std::log(6.0), 1.2, std::log(150.0), 1.4, 0.55};
    // A型 (中央値: 5 AU, 80 AU)
    binary_separation_params["A"]  = {std::log(5.0), 1.1, std::log(80.0),  1.3, 0.6};
    // F型 (中央値: 4 AU, 60 AU)
    binary_separation_params["F"]  = {std::log(4.0), 1.1, std::log(60.0),  1.3, 0.65};
    // G型 (太陽) (中央値: 3 AU, 50 AU)
    binary_separation_params["G"]  = {std::log(3.0), 1.0, std::log(50.0),  1.2, 0.7};
    // K型 (中央値: 2.5 AU, 40 AU)
    binary_separation_params["K"]  = {std::log(2.5), 1.0, std::log(40.0),  1.2, 0.75};
    // M型 (赤色矮星): ほとんどが近接連星 (中央値: 2 AU, 20 AU)
    binary_separation_params["M"]  = {std::log(2.0), 0.9, std::log(20.0),  1.1, 0.85};
    // 白色矮星 (中央値: 1 AU, 30 AU)
    binary_separation_params["WD"] = {std::log(1.0), 0.9, std::log(30.0),  1.2, 0.8};
    // 進化した星: 軌道が広がっていることが多い (中央値: 10 AU, 300 AU)
    binary_separation_params["ES"] = {std::log(10.0), 1.3, std::log(300.0), 1.6, 0.4};


    // For open clusters (younger, hotter stars are more common)
    // Weights for:                                M,    K,    G,    F,    A,    B,      O,        WD,    ES
    std::vector<double> cluster_spectral_weights = {0.60, 0.20, 0.10, 0.05, 0.02, 0.004, 0.0005, 0.015, 0.0105}; // Approximate weights, will be normalized
    normalize_weights(cluster_spectral_weights);

    // 3. Define the radius of the spherical space
    double R_parsecs = 100.0; // Default radius in parsecs (pc)
    std::cout << "半径Rパーセクの球形状空間の半径を入力してください (デフォルト: 100.0): ";
    std::string line;
    std::getline(std::cin, line);
    if (!line.empty()) {
        try {
            double input_r = std::stod(line);
            if (input_r > 0) {
                R_parsecs = input_r;
            } else {
                std::cout << "半径は正の値である必要があります。デフォルト値 " << R_parsecs << " を使用します。" << std::endl;
            }
        } catch (const std::exception&) {
            std::cout << "無効な入力です。デフォルト値 " << R_parsecs << " を使用します。" << std::endl;
        }
    }

    // 赤道・黄道の描画有無をユーザーに確認
    bool show_celestial_lines = true;
    std::cout << "星図に赤道と黄道を描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_celestial_lines = false;
    }

    // SVGに描画する最大等級をユーザーに確認
    double max_magnitude_for_svg = 5.0;
    std::cout << "星図に描画する星の最大等級を入力してください (デフォルト: 5.0): ";
    std::getline(std::cin, line);
    if (!line.empty()) {
        try {
            max_magnitude_for_svg = std::stod(line);
        } catch (const std::exception&) {
            std::cout << "無効な入力です。デフォルト値 " << max_magnitude_for_svg << " を使用します。" << std::endl;
        }
    }

    // 星団領域の描画有無をユーザーに確認
    bool show_cluster_regions = true;
    std::cout << "星図に星団の領域を描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_cluster_regions = false;
    }

    // 星のラベル描画有無をユーザーに確認
    bool show_star_labels = true;
    std::cout << "星図に明るい星のIDラベルを描画しますか？ (T/F, デフォルト: T): ";
    std::getline(std::cin, line);
    if (!line.empty() && (line[0] == 'F' || line[0] == 'f')) {
        show_star_labels = false;
    }

    // ラベルを描画する等級のしきい値をユーザーに確認
    double label_magnitude_threshold = 2.5;
    if (show_star_labels) {
        std::cout << "IDラベルを描画する星の最大等級を入力してください (デフォルト: 2.5): ";
        std::getline(std::cin, line);
        if (!line.empty()) {
            try {
                label_magnitude_threshold = std::stod(line);
            } catch (const std::exception&) {
                std::cout << "無効な入力です。デフォルト値 " << label_magnitude_threshold << " を使用します。" << std::endl;
            }
        }
    }

    // 4. Calculate the number of stars
    // Correct volume of a sphere: (4/3) * PI * R^3
    const double PI = std::acos(-1.0); // Portable way to get PI
    double volume = (4.0 / 3.0) * PI * std::pow(R_parsecs, 3);
    double star_density = 0.13; // Star system density in pc^-3

    long long num_stars = static_cast<long long>(volume * star_density);
    if (num_stars == 0) { // Ensure at least one star if the calculated number is zero
        num_stars = 1;
    }

    std::cout << R_parsecs << "パーセクの半径を持つ球状空間内に " << num_stars << " 個の星系を計算します。" << std::endl;
    std::cout << "半径が非常に大きい場合、時間がかかることがあります。" << std::endl;

    // 5. Distributions for position (uniform distribution within a sphere)
    // For uniform distribution in a sphere, r^3 should be uniform.
    // To avoid log10(0) for distance, ensure the lower bound is slightly above zero.
    std::uniform_real_distribution<double> dist_r_cube(std::numeric_limits<double>::epsilon(), 1.0);
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * PI); // Longitude [0, 2*PI)
    // For uniform distribution over a sphere's surface, sin(phi) should be uniform.
    std::uniform_real_distribution<double> dist_phi_sin(-1.0, 1.0); // For sin(phi) [-1, 1]

    // 6. Discrete distribution for spectral types
    std::discrete_distribution<> spectral_dist(spectral_type_weights.begin(), spectral_type_weights.end());

    // 7. Vector to store star systems
    std::vector<StarSystem> star_systems;
    star_systems.reserve(num_stars); // Pre-allocate memory for efficiency

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
    star_systems.insert(star_systems.end(), cluster_stars.begin(), cluster_stars.end());

    // Then, generate the remaining field stars
    long long num_field_stars = num_stars - star_systems.size();
    if (num_field_stars < 0) num_field_stars = 0; // Should not happen, but as a safeguard
    std::cout << "\n--- 現場星の生成 ---" << std::endl;
    std::cout << num_field_stars << " 個の主星系を並列で生成します..." << std::endl;

    std::vector<StarSystem> generated_field_stars;
    generated_field_stars.reserve(num_field_stars * 1.2); // 伴星の分も考慮して多めに確保

    // --- 再現性を保証する決定論的な並列処理 ---
    // 1. 乱数生成のためのマスターシードを取得
    unsigned int master_seed = gen();
    // 2. 星団生成後に使われた最後の連星IDを取得
    const long long last_cluster_binary_id = next_binary_system_id;

    #pragma omp parallel
    {
        std::vector<StarSystem> local_stars; // 各スレッドが生成した星を一時的に格納するベクター

        // ループの反復をスレッドに分割
        #pragma omp for schedule(dynamic)
        for (long long i = 0; i < num_field_stars; ++i) {
            // 3. 各星系ごとに決定論的な乱数生成器を初期化
            std::mt19937 thread_gen(master_seed + static_cast<unsigned int>(i));
            // 4. 状態を持つ分布オブジェクトもループ内で生成
            std::discrete_distribution<> thread_spectral_dist(spectral_type_weights.begin(), spectral_type_weights.end());

            StarSystem s;
            s.cluster_id = 0; // Field star
            // 5. 連星IDを決定論的に割り当て
            s.binary_system_id = last_cluster_binary_id + i;

            // Generate spherical position
            s.distance_r = R_parsecs * std::cbrt(dist_r_cube(thread_gen));
            s.longitude_theta = dist_theta(thread_gen);
            s.latitude_phi = std::asin(dist_phi_sin(thread_gen));

            // Convert spherical to Cartesian coordinates
            double cos_phi = std::cos(s.latitude_phi);
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
            s.apparent_magnitude = s.absolute_magnitude + 5.0 * (std::log10(s.distance_r) - 1.0);

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
            generated_field_stars.insert(generated_field_stars.end(), local_stars.begin(), local_stars.end());
        }
    }
    star_systems.insert(star_systems.end(), generated_field_stars.begin(), generated_field_stars.end());

    // Assign unique IDs to all generated stars (clusters + field)
    std::cout << "\n--- 最終処理 ---" << std::endl;
    std::cout << "全星系にユニークIDを割り当てています..." << std::endl;
    for(long long i = 0; i < star_systems.size(); ++i) {
        star_systems[i].id = i + 1;
    }

    // Count stars per spectral type for reporting
    std::map<std::string, long long> total_spectral_counts;
    std::map<int, std::map<std::string, long long>> cluster_spectral_counts;
    std::map<int, long long> stars_per_cluster_count;

    for (const auto& s : star_systems) {
        total_spectral_counts[s.spectral_type]++;
        if (s.cluster_id > 0) {
            cluster_spectral_counts[s.cluster_id][s.spectral_type]++;
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

    // CSV Header (Japanese as requested)
    ofs << "通し番号,連星系ID,視等級,スペクトル,絶対等級,距離r,経度θ,緯度φ,x,y,z,星団ID,星系絶対等級,星系視等級\n";

    // CSV Data
    ofs << std::fixed << std::setprecision(4); // Set precision for floating-point numbers
    for (const auto& s : star_systems) {
        ofs << s.id << ","
            << s.binary_system_id << ","
            << s.apparent_magnitude << ","
            << s.spectral_type << ","
            << s.absolute_magnitude << ","
            << s.distance_r << ","
            << s.longitude_theta << ","
            << s.latitude_phi << ","
            << s.x << ","
            << s.y << ","
            << s.z << ","
            << s.cluster_id << ","
            << s.system_absolute_magnitude << ","
            << s.system_apparent_magnitude << "\n";
    }

    ofs.close();
    std::cout << "星系データは " << output_filename << " に正常に書き込まれました。" << std::endl;

    // 10. Generate Star Maps
    std::cout << "\n--- 星図の生成 ---" << std::endl;
    // ステレオ投影（平射図法）の星図を生成
    generate_star_map_svg(star_systems, "star_map_stereographic_north.svg", cluster_properties, ProjectionType::Stereographic, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);
    // 正距方位図法の星図を生成
    generate_star_map_svg(star_systems, "star_map_azimuthal_equidistant_north.svg", cluster_properties, ProjectionType::AzimuthalEquidistant, show_celestial_lines, max_magnitude_for_svg, show_cluster_regions, show_star_labels, label_magnitude_threshold);

    return 0;
}