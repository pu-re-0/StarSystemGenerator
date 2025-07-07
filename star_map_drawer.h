#ifndef STAR_MAP_DRAWER_H
#define STAR_MAP_DRAWER_H

#include <vector>
#include <string>
#include "cluster_generator.h" // StarSystem, ClusterProperties を使うため

// 投影法の種類を定義するenum
enum class ProjectionType {
    Stereographic,          // 平射図法（ステレオ投影）
    AzimuthalEquidistant    // 正距方位図法
};

// 指定された投影法で星図をSVGファイルとして生成する汎用関数
void generate_star_map_svg(
    const std::vector<StarSystem>& all_stars,
    const std::string& filename,
    const std::vector<ClusterProperties>& clusters,
    ProjectionType proj_type,
    bool show_celestial_lines = true,
    float max_magnitude_to_plot = 4.0f,
    bool show_cluster_regions = true,
    bool show_star_labels = true,
    float label_magnitude_threshold = 2.5
);

#endif // STAR_MAP_DRAWER_H