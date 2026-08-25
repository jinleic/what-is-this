#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int N = 20;

struct Domain {
  std::array<int8_t, N> half{};
  std::array<int8_t, N> row{};
};

struct Search {
  std::array<uint8_t, N> internal{};
  std::array<std::array<int8_t, N>, N> w_half{};
  std::array<std::vector<Domain>, N> domains;
  std::array<std::array<bool, N>, N> tree_edge{};
  int tree_edges = 0;
  std::array<int, N> selected{};
  std::array<bool, N> assigned{};
  uint64_t nodes = 0;
  uint64_t checks = 0;

  bool compatible(int left, int left_domain, int right, int right_domain) {
    ++checks;
    const auto& a = domains[left][left_domain];
    const auto& b = domains[right][right_domain];
    if (a.half[right] != b.half[left]) return false;
    int dot = 0;
    for (int k = 0; k < N; ++k) dot += int(a.row[k]) * int(b.row[k]);
    return dot == 0;
  }

  bool visit(std::vector<std::vector<int>> candidates, int depth) {
    ++nodes;
    if (depth == N) return true;
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
        if (filtered.empty()) viable = false;
        else next[other] = std::move(filtered);
      }
      if (!viable) continue;
      assigned[row] = true;
      selected[row] = candidate;
      if (visit(std::move(next), depth + 1)) return true;
      assigned[row] = false;
    }
    return false;
  }

  void build_domains() {
    std::array<bool, N> seen{};
    std::queue<int> frontier;
    seen[0] = true;
    frontier.push(0);
    while (!frontier.empty()) {
      const int vertex = frontier.front();
      frontier.pop();
      for (int other = 0; other < N; ++other) {
        if (other == vertex || w_half[vertex][other] != 0 || seen[other]) continue;
        seen[other] = true;
        frontier.push(other);
        tree_edge[vertex][other] = tree_edge[other][vertex] = true;
        ++tree_edges;
      }
    }
    if (tree_edges != N - 1) {
      throw std::runtime_error("T support graph is unexpectedly disconnected");
    }
    for (int row = 0; row < N; ++row) {
      std::array<int, N> free_support{};
      int degree = 0;
      int free_degree = 0;
      Domain base;
      base.row[row] = 2 * internal[row] - 1;
      for (int other = 0; other < N; ++other) {
        if (other == row || w_half[row][other] != 0) continue;
        ++degree;
        if (tree_edge[row][other]) {
          base.half[other] = 1;
          base.row[other] = 2;
        } else {
          free_support[free_degree++] = other;
        }
      }
      if (degree != 11) {
        throw std::runtime_error("W complement does not have T degree 11");
      }
      domains[row].reserve(1u << free_degree);
      for (uint32_t signs = 0; signs < (1u << free_degree); ++signs) {
        Domain domain = base;
        for (int position = 0; position < free_degree; ++position) {
          const int other = free_support[position];
          const int8_t sign = signs >> position & 1 ? 1 : -1;
          domain.half[other] = sign;
          domain.row[other] = 2 * sign;
        }
        domains[row].push_back(domain);
      }
    }
  }
};

uint32_t read_u32(std::ifstream& stream) {
  unsigned char bytes[4];
  stream.read(reinterpret_cast<char*>(bytes), 4);
  if (!stream) throw std::runtime_error("truncated T input");
  return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) |
         (uint32_t(bytes[2]) << 8) | uint32_t(bytes[3]);
}

Search load(const std::string& path, uint32_t& source_index) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open T input");
  char magic[5];
  stream.read(magic, 5);
  if (std::string(magic, 5) != "TDOM1") throw std::runtime_error("bad T input magic");
  source_index = read_u32(stream);
  Search search;
  stream.read(reinterpret_cast<char*>(search.internal.data()), N);
  if (!stream) throw std::runtime_error("truncated T internal bits");
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right) {
      char encoded;
      stream.get(encoded);
      if (!stream || encoded < '0' || encoded > '2') {
        throw std::runtime_error("bad W half-trit input");
      }
      const int8_t value = int8_t(encoded - '1');
      search.w_half[left][right] = search.w_half[right][left] = value;
    }
  }
  char extra;
  if (stream.read(&extra, 1)) throw std::runtime_error("trailing T input bytes");
  search.build_domains();
  return search;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: involution_f5_t_square T_INPUT\n";
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
    const bool found = search.visit(std::move(candidates), 0);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"source_index\":" << source_index
              << ",\"t_square_satisfiable\":" << (found ? "true" : "false")
              << ",\"nodes\":" << search.nodes
              << ",\"compatibility_checks\":" << search.checks
              << ",\"switching_gauge_tree_edges\":" << search.tree_edges
              << ",\"wall_seconds\":" << seconds;
    if (found) {
      std::string witness;
      for (int left = 0; left < N; ++left) {
        for (int right = left + 1; right < N; ++right) {
          const int a = search.domains[left][search.selected[left]].half[right];
          const int b = search.domains[right][search.selected[right]].half[left];
          if (a != b) throw std::runtime_error("final T witness is not symmetric");
          witness.push_back(char('1' + a));
        }
      }
      std::cout << ",\"t_half_trits\":\"" << witness << "\"";
    }
    std::cout << "}\n";
    return found ? 0 : 20;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
