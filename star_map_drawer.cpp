#include "star_map_drawer.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> // For std::clamp
#include <utility>   // For std::pair

// ヘルパー関数：スペクトル型に基づいてSVGで使う色を返す
// このファイル内でのみ使用するため、static宣言します。

static std::string get_star_color_for_plot(SpectralType spectral_type) {
    switch (spectral_type) {
        case SpectralType::O:  return "#9bb0ff";  // 青白
        case SpectralType::B:  return "#aabfff";  // 青白
        case SpectralType::A:  return "#cad8ff";  // 白
        case SpectralType::F:  return "#f8f7ff";  // 黄白
        case SpectralType::G:  return "#fff4ea";  // 黄
        case SpectralType::K:  return "#ffd2a1";  // 橙
        case SpectralType::M:  return "#ffb56c";  // 赤
        case SpectralType::WD: return "#f0f0f0"; // 白
        case SpectralType::ES: return "#ffeea8"; // 明るい黄
        default:               return "white";
    }
}

// 指定された投影法で星図をSVGファイルとして生成する汎用関数
void generate_star_map_svg(const std::vector<StarSystem>& all_stars, const std::string& filename, const std::vector<ClusterProperties>& clusters, ProjectionType proj_type, bool show_celestial_lines = true, float max_magnitude_to_plot = 4.0f, bool show_cluster_regions = true, bool show_star_labels = true, float label_magnitude_threshold = 2.5) {
    const float PI = std::acos(-1.0f);
    const int svg_size = 1000;
    const float center = svg_size / 2.0f;
    const float plot_scale = svg_size / 2.0f;
    const float max_plot_radius = 8.0f; // 最も明るい星の半径
    const float min_plot_radius = 1.0f; // 6.0f等星の半径

    // 投影法ごとの設定
    float limit_latitude_deg;
    std::string projection_name;

    if (proj_type == ProjectionType::Stereographic) {
        limit_latitude_deg = 0.0f; // 北極から赤道まで
        projection_name = "Stereographic Projection (North Pole to Equator)";
    } else { // AzimuthalEquidistant
        limit_latitude_deg = -55.0f; // 北極から南緯55度まで
        projection_name = "Azimuthal Equidistant Projection (North Pole to -55 deg)";
    }
    const float limit_latitude_rad = limit_latitude_deg * PI / 180.0f;

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
        float lat60_colatitude = PI / 2.0f - (60.0f * PI / 180.0f);
        float r60 = plot_scale * std::tan(lat60_colatitude / 2.0f);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << r60 << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
        float lat30_colatitude = PI / 2.0f - (30.0f * PI / 180.0f);
        float r30 = plot_scale * std::tan(lat30_colatitude / 2.0f);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << r30 << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
    } else { // AzimuthalEquidistant
        // 南緯55度が外周になる。赤道と南緯30度線を描画
        const float max_angular_distance = (PI / 2.0f) - limit_latitude_rad;
        float equator_angular_distance = PI / 2.0f - (0.0f * PI / 180.0f);
        float equator_radius = plot_scale * (equator_angular_distance / max_angular_distance);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << equator_radius << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
        float lat30s_angular_distance = PI / 2.0f - (-30.0f * PI / 180.0f);
        float lat30s_radius = plot_scale * (lat30s_angular_distance / max_angular_distance);
        ofs << "<circle cx=\"" << center << "\" cy=\"" << center << "\" r=\"" << lat30s_radius << "\" fill=\"none\" stroke=\"#333\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>\n";
    }

    // 座標投影ヘルパー (ラムダ式)
    auto project_point = [&](float lon, float lat) -> std::pair<float, float> {
        float px, py;
        // 中心からの角距離 c = 90度 - 緯度
        float colatitude_c = (PI / 2.0f) - lat;

        if (proj_type == ProjectionType::Stereographic) {
            float r_plot = plot_scale * std::tan(colatitude_c / 2.0f);
            px = center + r_plot * std::cos(lon);
            py = center - r_plot * std::sin(lon);
        } else { // AzimuthalEquidistant
            const float max_angular_distance = (PI / 2.0f) - limit_latitude_rad;
            float r_plot = plot_scale * (colatitude_c / max_angular_distance);
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
            auto p = project_point(static_cast<float>(i) * PI / 180.0f, 0.0f);
            ofs << p.first << "," << p.second << " ";
        }
        ofs << "\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1.5\"/>\n";

        // --- 黄道の描画 (黄色) ---
        const float obliquity_rad = 23.44 * PI / 180.0f; // 黄道傾斜角
        ofs << "  <polyline points=\"";
        for (int i = 0; i <= 360; ++i) {
            float ecliptic_lon = static_cast<float>(i) * PI / 180.0f;
            float eq_lon = std::atan2(std::sin(ecliptic_lon) * std::cos(obliquity_rad), std::cos(ecliptic_lon));
            float eq_lat = std::asin(std::sin(obliquity_rad) * std::sin(ecliptic_lon));
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
                float angle = 2.0f * PI * i / num_points;

                // 三軸不等楕円体のXY平面上の輪郭（楕円）をローカル座標で生成
                // これにより、星図上での星団の見た目がより多様になる
                float x_local = cluster.axis_a * std::cos(angle);
                float y_local = cluster.axis_b * std::sin(angle);
                float z_local = 0;

                // オイラー角で回転 (Z-Y-X順)
                // Z軸周り (alpha)
                float x1 = x_local * std::cos(cluster.rot_alpha) - y_local * std::sin(cluster.rot_alpha);
                float y1 = x_local * std::sin(cluster.rot_alpha) + y_local * std::cos(cluster.rot_alpha);
                float z1 = z_local;
                // Y軸周り (beta)
                float x2 = x1 * std::cos(cluster.rot_beta) + z1 * std::sin(cluster.rot_beta);
                float y2 = y1;
                float z2 = -x1 * std::sin(cluster.rot_beta) + z1 * std::cos(cluster.rot_beta);
                // X軸周り (gamma)
                float x_rot = x2;
                float y_rot = y2 * std::cos(cluster.rot_gamma) - z2 * std::sin(cluster.rot_gamma);
                float z_rot = y2 * std::sin(cluster.rot_gamma) + z2 * std::cos(cluster.rot_gamma);

                // 星団の中心位置へ移動
                float world_x = x_rot + cluster.center_x;
                float world_y = y_rot + cluster.center_y;
                float world_z = z_rot + cluster.center_z;

                // ワールド座標から球面座標へ変換
                float dist_r = std::sqrt(world_x * world_x + world_y * world_y + world_z * world_z);
                if (dist_r < 1e-6) continue; // 中心点はスキップ

                float lat = std::asinf(std::clamp(world_z / dist_r, -1.0f, 1.0f));
                float lon = std::atan2f(world_y, world_x);

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
    // オブジェクトのコピーを避け、ポインタを格納するように変更
    std::vector<const StarSystem*> drawable_stars;
    float min_plot_magnitude = max_magnitude_to_plot;

    for (const auto& star : all_stars) {
        // 緯度でフィルタリング
        if (star.latitude_phi < limit_latitude_rad) {
            continue;
        }

        // 伴星 (-1) は描画しない
        if (star.system_apparent_magnitude == -1.0f) {
            continue;
        }

        // 描画に使う等級は、星系全体の等級
        float magnitude_to_check = star.system_apparent_magnitude;

        if (magnitude_to_check <= max_magnitude_to_plot) {
            drawable_stars.push_back(&star); // ポインタを追加
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

    std::cout << "星図 (" << projection_name << "): " << drawable_stars.size() << " 個の星を描画しました。" << std::endl;

    // 2. フィルタリングされた星を描画
    float magnitude_range = max_magnitude_to_plot - min_plot_magnitude;

    for (const auto* star_ptr : drawable_stars) {
        auto p = project_point(star_ptr->longitude_theta, star_ptr->latitude_phi);

        // 描画に使う等級は、星系全体の等級
        float magnitude_for_size = star_ptr->system_apparent_magnitude;

        // 等級に基づいて星の描画半径を計算
        float plot_radius;
        if (magnitude_range < 1e-6) { // 等級の範囲がほぼ0の場合（星が1つ、または全て同じ明るさ）
            plot_radius = max_plot_radius;
        } else {
            plot_radius = min_plot_radius + (max_plot_radius - min_plot_radius) * (max_magnitude_to_plot - magnitude_for_size) / magnitude_range;
        }

        ofs << "  <circle cx=\"" << p.first << "\" cy=\"" << p.second << "\" r=\"" << plot_radius << "\" fill=\"" << get_star_color_for_plot(star_ptr->spectral_type) << "\"/>\n";
    }

    // 3. 星のラベルを描画
    if (show_star_labels) {
        long long labeled_star_count = 0;
        ofs << "  <!-- Star Labels -->\n";
        for (const auto* star_ptr : drawable_stars) {
            // ラベル表示対象の等級を決定
            float magnitude_for_label = star_ptr->system_apparent_magnitude;

            // 等級がしきい値以下で、かつ名前が割り当てられている場合のみラベルを描画
            if (magnitude_for_label <= label_magnitude_threshold && star_ptr->name && !star_ptr->name->empty()) {
                auto p = project_point(star_ptr->longitude_theta, star_ptr->latitude_phi);
                float text_x = p.first + 5;  // x座標を少し右にずらす
                float text_y = p.second + 5; // y座標を少し下にずらす

                std::string text_anchor = "start"; // デフォルトは左揃え

                // ラベルがクリッピング範囲外にはみ出さないように簡易チェック
                // 右端に近い場合は、星の左側に右揃えで表示する
                if (text_x > svg_size - 60) {
                    text_x = p.first - 5;
                    text_anchor = "end";
                }
                if (text_y > svg_size - 5) text_y = p.second - 5;

                // 星の名前をラベルとして表示
                ofs << "  <text x=\"" << text_x << "\" y=\"" << text_y << "\" font-family=\"sans-serif\" font-size=\"8\" fill=\"#aaa\" text-anchor=\"" << text_anchor << "\">" << *star_ptr->name << "</text>\n";
                labeled_star_count++;
            }
        }
        ofs << "  <!-- End Star Labels -->\n";
        std::cout << "星図 (" << projection_name << "): " << labeled_star_count << " 個の星にラベルを付けました。" << std::endl;

        // --- 星団のラベルを描画 ---
        ofs << "  <!-- Cluster Labels -->\n";
        for (const auto& cluster : clusters) {
            // 星団の中心座標をワールド座標から球面座標へ変換
            float dist_r = std::sqrt(cluster.center_x * cluster.center_x + cluster.center_y * cluster.center_y + cluster.center_z * cluster.center_z);
            if (dist_r < 1e-6) continue;

            float lat = std::asinf(std::clamp(cluster.center_z / dist_r, -1.0f, 1.0f));
            float lon = std::atan2f(cluster.center_y, cluster.center_x);

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