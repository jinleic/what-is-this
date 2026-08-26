#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int N = 20;
constexpr int F = 5;

struct Domain {
  std::array<int8_t, N> x{};
  std::array<int8_t, N> w{};
};

struct Search {
  std::array<uint8_t, N> masks{};
  std::array<uint8_t, N> internal{};
  std::array<std::vector<Domain>, N> domains;
  std::array<std::array<int, N>, N> target{};
  std::array<int, N> selected{};
  std::array<bool, N> assigned{};
  uint64_t nodes = 0;
  uint64_t checks = 0;
  std::vector<std::string> solutions;

  bool compatible(int left, int left_domain, int right, int right_domain) {
    ++checks;
    const auto& a = domains[left][left_domain];
    const auto& b = domains[right][right_domain];
    if (a.x[right] != b.x[left]) return false;
    int dot = 0;
    for (int k = 0; k < N; ++k) dot += int(a.w[k]) * int(b.w[k]);
    return dot == target[left][right];
  }

  bool visit(std::vector<std::vector<int>> candidates, int depth) {
    ++nodes;
    if (depth == N) {
      std::string witness;
      for (int left = 0; left < N; ++left) {
        for (int right = left + 1; right < N; ++right) {
          const int a = domains[left][selected[left]].x[right];
          const int b = domains[right][selected[right]].x[left];
          if (a != b) throw std::runtime_error("W solution is not symmetric");
          witness.push_back(char('1' + a));
        }
      }
      solutions.push_back(std::move(witness));
      return false;
    }
    int row = -1;
    size_t best = std::numeric_limits<size_t>::max();
    for (int i = 0; i < N; ++i) {
      if (!assigned[i] && candidates[i].size() < best) {
        row = i;
        best = candidates[i].size();
      }
    }
    if (row < 0 || best == 0) return false;
    for (int candidate : candidates[row]) {
      bool agrees = true;
      for (int other = 0; other < N && agrees; ++other) {
        if (assigned[other]) {
          agrees = compatible(row, candidate, other, selected[other]);
        }
      }
      if (!agrees) continue;
      auto next = candidates;
      next[row].clear();
      bool viable = true;
      for (int other = 0; other < N && viable; ++other) {
        if (assigned[other] || other == row) continue;
        std::vector<int> filtered;
        filtered.reserve(next[other].size());
        for (int option : next[other]) {
          if (compatible(row, candidate, other, option)) filtered.push_back(option);
        }
        if (filtered.empty()) {
          viable = false;
        } else {
          next[other] = std::move(filtered);
        }
      }
      if (!viable) continue;
      assigned[row] = true;
      selected[row] = candidate;
      visit(std::move(next), depth + 1);
      assigned[row] = false;
    }
    return false;
  }
};

uint32_t read_u32(std::ifstream& stream) {
  unsigned char bytes[4];
  stream.read(reinterpret_cast<char*>(bytes), 4);
  if (!stream) throw std::runtime_error("truncated W domain dump");
  return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) |
         (uint32_t(bytes[2]) << 8) | uint32_t(bytes[3]);
}

Search load(const std::string& path, uint32_t& source_index) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open W domain dump");
  char magic[5];
  stream.read(magic, 5);
  if (std::string(magic, 5) != "WDOM1") {
    throw std::runtime_error("bad W domain magic");
  }
  source_index = read_u32(stream);
  Search search;
  stream.read(reinterpret_cast<char*>(search.masks.data()), N);
  stream.read(reinterpret_cast<char*>(search.internal.data()), N);
  if (!stream) throw std::runtime_error("truncated W support data");
  for (int row = 0; row < N; ++row) {
    const uint32_t count = read_u32(stream);
    search.domains[row].reserve(count);
    for (uint32_t item = 0; item < count; ++item) {
      uint32_t code = read_u32(stream);
      Domain domain;
      for (int index = N - 1; index >= 0; --index) {
        domain.x[index] = int8_t(code % 3) - 1;
        code /= 3;
      }
      if (code != 0 || domain.x[row] != 0) {
        throw std::runtime_error("bad W domain code");
      }
      for (int index = 0; index < N; ++index) {
        domain.w[index] = 2 * domain.x[index];
      }
      domain.w[row] = 1 - 2 * search.internal[row];
      search.domains[row].push_back(domain);
    }
  }
  char extra;
  if (stream.read(&extra, 1)) {
    throw std::runtime_error("trailing W domain bytes");
  }
  for (int left = 0; left < N; ++left) {
    for (int right = 0; right < N; ++right) {
      int gram = 0;
      for (int vertex = 0; vertex < F; ++vertex) {
        const int a = 1 - 2 * ((search.masks[left] >> vertex) & 1);
        const int b = 1 - 2 * ((search.masks[right] >> vertex) & 1);
        gram += a * b;
      }
      search.target[left][right] = 45 * (left == right) - 2 - 2 * gram;
    }
  }
  return search;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: involution_f5_w_enumerate DOMAIN_DUMP OUTPUT\n";
    return 2;
  }
  try {
    uint32_t source_index = 0;
    Search search = load(argv[1], source_index);
    std::vector<std::vector<int>> candidates(N);
    for (int row = 0; row < N; ++row) {
      candidates[row].reserve(search.domains[row].size());
      for (size_t domain = 0; domain < search.domains[row].size(); ++domain) {
        candidates[row].push_back(int(domain));
      }
    }
    const auto started = std::chrono::steady_clock::now();
    search.visit(std::move(candidates), 0);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::ofstream output(argv[2], std::ios::binary);
    if (!output) throw std::runtime_error("cannot open W solution output");
    auto put_u32 = [&](uint32_t value) {
      const unsigned char bytes[4] = {
          static_cast<unsigned char>(value >> 24),
          static_cast<unsigned char>(value >> 16),
          static_cast<unsigned char>(value >> 8),
          static_cast<unsigned char>(value)};
      output.write(reinterpret_cast<const char*>(bytes), 4);
    };
    output.write("WLIST1", 6);
    put_u32(source_index);
    put_u32(search.solutions.size());
    for (const std::string& witness : search.solutions) output << witness;
    if (!output) throw std::runtime_error("failed while writing W solutions");
    std::cout << "{\"source_index\":" << source_index
              << ",\"w_solutions\":" << search.solutions.size()
              << ",\"nodes\":" << search.nodes
              << ",\"compatibility_checks\":" << search.checks
              << ",\"wall_seconds\":" << seconds << "}\n";
    return search.solutions.empty() ? 20 : 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << "\n";
    return 1;
  }
}
