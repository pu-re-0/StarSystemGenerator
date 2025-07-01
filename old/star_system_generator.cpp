#include <iostream>
#include <vector>
#include <string>
#include <cmath>     // For std::acos, std::pow, std::log10, std::cbrt, std::asin
#include <random>    // For std::random_device, std::mt19937, std::uniform_real_distribution, std::discrete_distribution
#include <fstream>   // For std::ofstream
#include <iomanip>   // For std::fixed, std::setprecision
#include <map>       // For std::map
#include <limits>    // For std::numeric_limits

// Struct to hold star system properties
struct StarSystem {
    long long id;
    double apparent_magnitude;
    std::string spectral_type;
    double absolute_magnitude;
    double distance_r;
    double longitude_theta; // 経度 (radians)
    double latitude_phi;    // 緯度 (radians)
};

int main() {
    // 1. Setup random number generator
    // Using std::random_device for non-deterministic seeding, or a fixed seed for reproducibility
    std::random_device rd;
    // std::mt19937 gen(rd()); // Seed with a random device 
    std::mt19937 gen(2025062809); // For reproducible results, use a fixed seed like 42  数値はペルド・ペルソーネより

    // 2. Define spectral type properties (probabilities and absolute magnitude ranges)
    // Probabilities are approximate and sum to 1.0
    std::vector<std::string> spectral_types_vec = {"M", "K", "G", "F", "A", "B", "O", "WD", "RG", "ES"};
    std::vector<double> spectral_type_weights = {0.70, 0.15, 0.07, 0.03, 0.01, 0.005, 0.0005, 0.02, 0.01, 0.0045}; // Sums to 1.0

    // Absolute magnitude ranges for each spectral type {min_M, max_M}
    std::map<std::string, std::pair<double, double>> abs_mag_ranges;
    abs_mag_ranges["O"] = {-6.0, -4.0};
    abs_mag_ranges["B"] = {-4.0, -1.0};
    abs_mag_ranges["A"] = {-1.0, 2.0};
    abs_mag_ranges["F"] = {2.0, 3.5};
    abs_mag_ranges["G"] = {3.5, 5.0};
    abs_mag_ranges["K"] = {5.0, 7.5};
    abs_mag_ranges["M"] = {7.5, 12.0};
    abs_mag_ranges["WD"] = {10.0, 16.0}; // White Dwarf
    abs_mag_ranges["RG"] = {-3.0, 0.0};  // Red Giant
    abs_mag_ranges["ES"] = {-8.0, -4.0}; // Exotic/Supergiant (e.g., Supergiants)

    // 3. Define the radius of the spherical space
    double R_parsecs = 100.0; // Default radius in parsecs (pc)
    std::cout << "半径Rパーセクの球形状空間の半径を入力してください (例: 100.0): ";
    std::cin >> R_parsecs;

    if (R_parsecs <= 0) {
        std::cerr << "エラー: 半径は正の値である必要があります。" << std::endl;
        return 1;
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

    // 8. Generate stars
    for (long long i = 0; i < num_stars; ++i) {
        StarSystem s;
        s.id = i + 1;

        // Generate position
        s.distance_r = R_parsecs * std::cbrt(dist_r_cube(gen));
        s.longitude_theta = dist_theta(gen);
        s.latitude_phi = std::asin(dist_phi_sin(gen));

        // Generate spectral type
        int spectral_idx = spectral_dist(gen);
        s.spectral_type = spectral_types_vec[spectral_idx];

        // Generate absolute magnitude based on spectral type
        auto range_it = abs_mag_ranges.find(s.spectral_type);
        if (range_it != abs_mag_ranges.end()) {
            std::uniform_real_distribution<double> dist_abs_mag(range_it->second.first, range_it->second.second);
            s.absolute_magnitude = dist_abs_mag(gen);
        } else {
            // Fallback if spectral type not found (shouldn't happen with current setup)
            s.absolute_magnitude = 5.0; // Default to Sun-like if unknown
        }

        // Calculate apparent magnitude: m = M + 5 * (log10(d) - 1)
        // This formula is equivalent to m = M + 5 * log10(d/10pc)
        s.apparent_magnitude = s.absolute_magnitude + 5.0 * (std::log10(s.distance_r) - 1.0);

        star_systems.push_back(s);
    }

    // 9. Output to CSV file
    const std::string output_filename = "star_systems.csv";
    std::ofstream ofs(output_filename);
    if (!ofs.is_open()) {
        std::cerr << "エラー: " << output_filename << " を書き込み用に開けませんでした。" << std::endl;
        return 1;
    }

    // CSV Header (Japanese as requested)
    ofs << "通し番号,視等級,スペクトル,絶対等級,距離r,経度θ,緯度φ\n";

    // CSV Data
    ofs << std::fixed << std::setprecision(4); // Set precision for floating-point numbers
    for (const auto& s : star_systems) {
        ofs << s.id << ","
            << s.apparent_magnitude << ","
            << s.spectral_type << ","
            << s.absolute_magnitude << ","
            << s.distance_r << ","
            << s.longitude_theta << ","
            << s.latitude_phi << "\n";
    }

    ofs.close();
    std::cout << "星系データは " << output_filename << " に正常に書き込まれました。" << std::endl;

    return 0;
}