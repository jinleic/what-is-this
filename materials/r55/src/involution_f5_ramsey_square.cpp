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
constexpr int F = 5;
bool has_clique(const std::array<uint64_t, 45>& adjacency,
                uint64_t candidates, int need) {
  if (need == 0) return true;
  if (__builtin_popcountll(candidates) < need) return false;
  while (candidates) {
    const uint64_t bit = candidates & -candidates;
    candidates ^= bit;
    const int vertex = __builtin_ctzll(bit);
    if (has_clique(adjacency, candidates & adjacency[vertex], need - 1)) {
      return true;
    }
  }
  return false;
}



struct TDomain {
  std::array<int8_t, N> half{};
  std::array<int8_t, N> row{};
};

struct TSearch {
  std::array<uint8_t, N> masks{};
  std::array<uint8_t, N> internal{};
  std::array<std::array<int8_t, N>, N> w_half{};
  std::array<std::vector<TDomain>, N> domains;
  std::array<int, N> selected{};
  std::array<bool, N> assigned{};
  uint64_t nodes = 0;
  uint64_t checks = 0;
  uint64_t signed_completions = 0;
  uint64_t with_k5 = 0;
  uint64_t with_i5 = 0;

  void build_domains() {
    std::array<std::array<bool, N>, N> tree_edge{};
    std::array<bool, N> seen{};
    std::queue<int> frontier;
    int tree_edges = 0;
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
    if (tree_edges != N - 1) throw std::runtime_error("T support is disconnected");
    for (int row = 0; row < N; ++row) {
      std::array<int, N> free_support{};
      int degree = 0;
      int free_degree = 0;
      TDomain base;
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
      if (degree != 11) throw std::runtime_error("T support degree is not 11");
      domains[row].reserve(1u << free_degree);
      for (uint32_t signs = 0; signs < (1u << free_degree); ++signs) {
        TDomain domain = base;
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

  bool compatible(int left, int left_domain, int right, int right_domain) {
    ++checks;
    const auto& a = domains[left][left_domain];
    const auto& b = domains[right][right_domain];
    if (a.half[right] != b.half[left]) return false;
    int dot = 0;
    for (int k = 0; k < N; ++k) dot += int(a.row[k]) * int(b.row[k]);
    return dot == 0;
  }
  bool ramsey_good() {
    ++signed_completions;
    std::array<uint64_t, 45> adjacency{};
    auto add_edge = [&](int left, int right) {
      adjacency[left] |= uint64_t{1} << right;
      adjacency[right] |= uint64_t{1} << left;
    };
    for (int left = 0; left < F; ++left) {
      for (int right = left + 1; right < F; ++right) {
        if ((left - right + F) % F == 1 || (right - left + F) % F == 1) {
          add_edge(left, right);
        }
      }
    }
    for (int orbit = 0; orbit < N; ++orbit) {
      const int first = F + 2 * orbit;
      const int second = first + 1;
      for (int fixed = 0; fixed < F; ++fixed) {
        if (masks[orbit] >> fixed & 1) {
          add_edge(fixed, first);
          add_edge(fixed, second);
        }
      }
      if (internal[orbit]) add_edge(first, second);
    }
    for (int left = 0; left < N; ++left) {
      for (int right = left + 1; right < N; ++right) {
        const int w = w_half[left][right];
        const int t = domains[left][selected[left]].half[right];
        const int parallel = (1 - w - t) / 2;
        const int crossed = (1 - w + t) / 2;
        const int left_first = F + 2 * left;
        const int right_first = F + 2 * right;
        if (parallel) {
          add_edge(left_first, right_first);
          add_edge(left_first + 1, right_first + 1);
        }
        if (crossed) {
          add_edge(left_first, right_first + 1);
          add_edge(left_first + 1, right_first);
        }
      }
    }
    const uint64_t full = (uint64_t{1} << 45) - 1;
    std::array<uint64_t, 45> complement{};
    for (int vertex = 0; vertex < 45; ++vertex) {
      complement[vertex] = full ^ (uint64_t{1} << vertex) ^ adjacency[vertex];
    }
    const bool k5 = has_clique(adjacency, full, 5);
    const bool i5 = has_clique(complement, full, 5);
    with_k5 += k5;
    with_i5 += i5;
    return !k5 && !i5;
  }



  bool visit(std::vector<std::vector<int>> candidates, int depth) {
    ++nodes;
    if (depth == N) return ramsey_good();
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

  bool solve() {
    build_domains();
    std::vector<std::vector<int>> candidates(N);
    for (int row = 0; row < N; ++row) {
      for (size_t domain = 0; domain < domains[row].size(); ++domain) {
        candidates[row].push_back(int(domain));
      }
    }
    return visit(std::move(candidates), 0);
  }
};

struct WDomain {
  std::array<int8_t, N> half{};
  std::array<int8_t, N> row{};
};

struct SignedSearch {
  std::array<uint8_t, N> masks{};
  std::array<uint8_t, N> internal{};
  std::array<std::vector<WDomain>, N> domains;
  std::array<std::array<int, N>, N> target{};
  std::array<int, N> selected{};
  std::array<bool, N> assigned{};
  std::array<std::array<int8_t, N>, N> t_witness{};
  uint64_t w_nodes = 0;
  uint64_t w_checks = 0;
  uint64_t w_completions = 0;
  uint64_t t_nodes = 0;
  uint64_t t_checks = 0;
  uint64_t signed_completions = 0;
  uint64_t signed_with_k5 = 0;
  uint64_t signed_with_i5 = 0;

  bool compatible(int left, int left_domain, int right, int right_domain) {
    ++w_checks;
    const auto& a = domains[left][left_domain];
    const auto& b = domains[right][right_domain];
    if (a.half[right] != b.half[left]) return false;
    int dot = 0;
    for (int k = 0; k < N; ++k) dot += int(a.row[k]) * int(b.row[k]);
    return dot == target[left][right];
  }

  bool complete_t() {
    ++w_completions;
    TSearch search;
    search.internal = internal;
    search.masks = masks;
    for (int left = 0; left < N; ++left) {
      for (int right = 0; right < N; ++right) {
        search.w_half[left][right] = domains[left][selected[left]].half[right];
      }
    }
    const bool found = search.solve();
    t_nodes += search.nodes;
    t_checks += search.checks;
    signed_completions += search.signed_completions;
    signed_with_k5 += search.with_k5;
    signed_with_i5 += search.with_i5;
    if (found) {
      for (int left = 0; left < N; ++left) {
        for (int right = 0; right < N; ++right) {
          t_witness[left][right] = search.domains[left][search.selected[left]].half[right];
        }
      }
    }
    return found;
  }

  bool visit(std::vector<std::vector<int>> candidates, int depth) {
    ++w_nodes;
    if (depth == N) return complete_t();
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
};

uint32_t read_u32(std::ifstream& stream) {
  unsigned char bytes[4];
  stream.read(reinterpret_cast<char*>(bytes), 4);
  if (!stream) throw std::runtime_error("truncated W domain dump");
  return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) |
         (uint32_t(bytes[2]) << 8) | uint32_t(bytes[3]);
}

SignedSearch load(const std::string& path, uint32_t& source_index) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open W domain dump");
  char magic[5];
  stream.read(magic, 5);
  if (std::string(magic, 5) != "WDOM1") throw std::runtime_error("bad W domain magic");
  source_index = read_u32(stream);
  SignedSearch search;
  stream.read(reinterpret_cast<char*>(search.masks.data()), N);
  stream.read(reinterpret_cast<char*>(search.internal.data()), N);
  if (!stream) throw std::runtime_error("truncated W support data");
  for (int row = 0; row < N; ++row) {
    const uint32_t count = read_u32(stream);
    search.domains[row].reserve(count);
    for (uint32_t item = 0; item < count; ++item) {
      uint32_t code = read_u32(stream);
      WDomain domain;
      for (int index = N - 1; index >= 0; --index) {
        domain.half[index] = int8_t(code % 3) - 1;
        code /= 3;
      }
      if (code != 0 || domain.half[row] != 0) throw std::runtime_error("bad W code");
      for (int index = 0; index < N; ++index) domain.row[index] = 2 * domain.half[index];
      domain.row[row] = 1 - 2 * search.internal[row];
      search.domains[row].push_back(domain);
    }
  }
  char extra;
  if (stream.read(&extra, 1)) throw std::runtime_error("trailing W domain bytes");
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
  if (argc != 2) {
    std::cerr << "usage: involution_f5_ramsey_square W_DOMAIN_DUMP\n";
    return 2;
  }
  try {
    uint32_t source_index = 0;
    SignedSearch search = load(argv[1], source_index);
    std::vector<std::vector<int>> candidates(N);
    for (int row = 0; row < N; ++row) {
      for (size_t domain = 0; domain < search.domains[row].size(); ++domain) {
        candidates[row].push_back(int(domain));
      }
    }
    const auto started = std::chrono::steady_clock::now();
    const bool found = search.visit(std::move(candidates), 0);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"source_index\":" << source_index
              << ",\"ramsey_square_satisfiable\":" << (found ? "true" : "false")
              << ",\"w_nodes\":" << search.w_nodes
              << ",\"w_compatibility_checks\":" << search.w_checks
              << ",\"w_completions_tested\":" << search.w_completions
              << ",\"t_nodes\":" << search.t_nodes
              << ",\"t_compatibility_checks\":" << search.t_checks
              << ",\"signed_completions_tested\":" << search.signed_completions
              << ",\"signed_completions_with_k5\":" << search.signed_with_k5
              << ",\"signed_completions_with_i5\":" << search.signed_with_i5
              << ",\"wall_seconds\":" << seconds;
    if (found) {
      std::string w_witness, t_witness;
      for (int left = 0; left < N; ++left) {
        for (int right = left + 1; right < N; ++right) {
          const int w = search.domains[left][search.selected[left]].half[right];
          const int t = search.t_witness[left][right];
          w_witness.push_back(char('1' + w));
          t_witness.push_back(char('1' + t));
        }
      }
      std::cout << ",\"w_half_trits\":\"" << w_witness
                << "\",\"t_half_trits\":\"" << t_witness << "\"";
    }
    std::cout << "}\n";
    return found ? 0 : 20;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
