#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
constexpr int N = 20;

struct Domain {
  std::array<int8_t, N> half{};
  std::array<int8_t, N> row{};
};

struct Instance {
  std::array<uint8_t, N> internal{};
  std::array<std::array<int8_t, N>, N> w_half{};
  std::array<std::vector<Domain>, N> domains;
  uint32_t source_index = 0;
};

struct ProofNode {
  std::vector<std::pair<int, int>> partial;
  std::vector<size_t> children;
  int next_row = -1;
};

using SupportKey = std::tuple<int, int, int>;

uint32_t read_u32(std::ifstream& stream) {
  unsigned char bytes[4];
  stream.read(reinterpret_cast<char*>(bytes), 4);
  if (!stream) throw std::runtime_error("truncated T input");
  return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) |
         (uint32_t(bytes[2]) << 8) | uint32_t(bytes[3]);
}

Instance load(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open T input");
  char magic[5];
  stream.read(magic, 5);
  if (std::string(magic, 5) != "TDOM1") {
    throw std::runtime_error("bad T input magic");
  }
  Instance instance;
  instance.source_index = read_u32(stream);
  stream.read(reinterpret_cast<char*>(instance.internal.data()), N);
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right) {
      char encoded;
      stream.get(encoded);
      if (!stream || encoded < '0' || encoded > '2') {
        throw std::runtime_error("bad W half-trit input");
      }
      const int8_t value = int8_t(encoded - '1');
      instance.w_half[left][right] = instance.w_half[right][left] = value;
    }
  }
  char extra;
  if (stream.read(&extra, 1)) throw std::runtime_error("trailing T input");

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
      if (other == vertex || instance.w_half[vertex][other] != 0
          || seen[other]) {
        continue;
      }
      seen[other] = true;
      frontier.push(other);
      tree_edge[vertex][other] = tree_edge[other][vertex] = true;
      ++tree_edges;
    }
  }
  if (tree_edges != N - 1) {
    throw std::runtime_error("T support is disconnected");
  }
  for (int row = 0; row < N; ++row) {
    std::array<int, N> free_support{};
    int degree = 0;
    int free_degree = 0;
    Domain base;
    base.row[row] = 2 * instance.internal[row] - 1;
    for (int other = 0; other < N; ++other) {
      if (other == row || instance.w_half[row][other] != 0) continue;
      ++degree;
      if (tree_edge[row][other]) {
        base.half[other] = 1;
        base.row[other] = 2;
      } else {
        free_support[free_degree++] = other;
      }
    }
    if (degree != 11) throw std::runtime_error("T support degree is not 11");
    instance.domains[row].reserve(1u << free_degree);
    for (uint32_t signs = 0; signs < (1u << free_degree); ++signs) {
      Domain domain = base;
      for (int position = 0; position < free_degree; ++position) {
        const int other = free_support[position];
        const int8_t sign = signs >> position & 1 ? 1 : -1;
        domain.half[other] = sign;
        domain.row[other] = 2 * sign;
      }
      instance.domains[row].push_back(domain);
    }
  }
  return instance;
}

std::vector<int> parse_rows(const std::string& text) {
  std::vector<int> rows;
  std::stringstream stream(text);
  std::string token;
  std::array<bool, N> seen{};
  while (std::getline(stream, token, ',')) {
    const int row = std::stoi(token);
    if (row < 0 || row >= N || seen[row]) throw std::runtime_error("invalid rows");
    seen[row] = true;
    rows.push_back(row);
  }
  if (rows.size() < 2) throw std::runtime_error("at least two rows are required");
  return rows;
}

bool compatible(const Instance& instance, int left, int left_domain,
                int right, int right_domain) {
  const auto& a = instance.domains[left][left_domain];
  const auto& b = instance.domains[right][right_domain];
  if (a.half[right] != b.half[left]) return false;
  int dot = 0;
  for (int k = 0; k < N; ++k) dot += int(a.row[k]) * int(b.row[k]);
  return dot == 0;
}

void verify_full_rows(const std::vector<int>& rows) {
  if (rows.size() != N) {
    throw std::runtime_error("T certificate requires all twenty rows");
  }
  std::array<bool, N> seen{};
  for (int row : rows) seen[row] = true;
  if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
    throw std::runtime_error("T certificate row set is incomplete");
  }
}

struct TraceBuilder {
  const Instance& instance;
  const std::vector<int>& rows;
  std::array<bool, N> assigned{};
  std::array<int, N> selected{};
  std::map<SupportKey, uint64_t> support_ids;
  std::vector<ProofNode> nodes;
  TraceBuilder(const Instance& input, const std::vector<int>& active_rows)
      : instance(input), rows(active_rows) {}


  size_t add_node(int next_row, std::vector<size_t> children = {}) {
    ProofNode node;
    node.next_row = next_row;
    node.children = std::move(children);
    for (int row : rows) {
      if (assigned[row]) node.partial.push_back({row, selected[row]});
    }
    if (node.partial.empty() && nodes.empty()) {
      throw std::runtime_error("empty root cannot be a leaf");
    }
    for (const auto [row, domain] : node.partial) {
      support_ids.try_emplace({row, domain, next_row}, 0);
    }
    nodes.push_back(std::move(node));
    return nodes.size() - 1;
  }

  size_t visit(std::vector<std::vector<int>> candidates, size_t depth) {
    if (depth == rows.size()) {
      throw std::runtime_error("row core is satisfiable");
    }
    int row = -1;
    size_t best = std::numeric_limits<size_t>::max();
    for (int candidate_row : rows) {
      if (!assigned[candidate_row]
          && candidates[candidate_row].size() < best) {
        row = candidate_row;
        best = candidates[candidate_row].size();
      }
    }
    if (row < 0 || best == 0) {
      throw std::runtime_error("invalid DFS state");
    }
    std::vector<size_t> children;
    for (int domain : candidates[row]) {
      auto next = candidates;
      next[row].clear();
      bool viable = true;
      int blocking_row = -1;
      for (int other : rows) {
        if (!viable || assigned[other] || other == row) continue;
        std::vector<int> filtered;
        filtered.reserve(next[other].size());
        for (int option : next[other]) {
          if (compatible(instance, row, domain, other, option)) {
            filtered.push_back(option);
          }
        }
        if (filtered.empty()) {
          viable = false;
          blocking_row = other;
        } else {
          next[other] = std::move(filtered);
        }
      }
      assigned[row] = true;
      selected[row] = domain;
      if (viable) {
        children.push_back(visit(std::move(next), depth + 1));
      } else {
        children.push_back(add_node(blocking_row));
      }
      assigned[row] = false;
    }
    return add_node(row, std::move(children));
  }

  void build() {
    std::vector<std::vector<int>> candidates(N);
    for (int row : rows) {
      for (size_t domain = 0; domain < instance.domains[row].size();
           ++domain) {
        candidates[row].push_back(int(domain));
      }
    }
    const size_t root = visit(std::move(candidates), 0);
    if (root + 1 != nodes.size() || !nodes.back().partial.empty()) {
      throw std::runtime_error("last proof node is not the root");
    }
  }
};
}  // namespace

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: involution_f5_t_dfs_certificate T_INPUT ROWS FORMULA PROOF\n";
    return 2;
  }
  try {
    const Instance instance = load(argv[1]);
    const std::vector<int> rows = parse_rows(argv[2]);
    verify_full_rows(rows);
    TraceBuilder trace(instance, rows);
    trace.build();

    std::array<uint32_t, N> offset{};
    uint32_t variables = 0;
    for (int row : rows) {
      offset[row] = variables;
      variables += instance.domains[row].size();
    }
    const uint64_t formula_constraints = 2 * rows.size() + trace.support_ids.size();
    std::ofstream formula(argv[3]);
    if (!formula) throw std::runtime_error("cannot open formula output");
    formula << "* #variable= " << variables
            << " #constraint= " << formula_constraints
            << " #equal= 0 intsize= 64\n";
    formula << "* exact T row-domain DFS certificate, source_index="
            << instance.source_index << "\n";
    std::array<uint64_t, N> atleast_id{};
    std::array<uint64_t, N> atmost_id{};
    uint64_t formula_id = 0;
    for (int row : rows) {
      atleast_id[row] = ++formula_id;
      for (size_t domain = 0; domain < instance.domains[row].size(); ++domain) {
        formula << "+1 x" << (offset[row] + domain + 1) << ' ';
      }
      formula << ">= 1;\n";
      atmost_id[row] = ++formula_id;
      for (size_t domain = 0; domain < instance.domains[row].size(); ++domain) {
        formula << "-1 x" << (offset[row] + domain + 1) << ' ';
      }
      formula << ">= -1;\n";
    }
    for (auto& [key, id] : trace.support_ids) {
      id = ++formula_id;
      const auto [row, domain, other] = key;
      formula << "-1 x" << (offset[row] + domain + 1) << ' ';
      for (size_t option = 0; option < instance.domains[other].size(); ++option) {
        if (compatible(instance, row, domain, other, option)) {
          formula << "+1 x" << (offset[other] + option + 1) << ' ';
        }
      }
      formula << ">= 0;\n";
    }
    if (formula_id != formula_constraints || !formula) {
      throw std::runtime_error("formula constraint count mismatch");
    }

    std::ofstream proof(argv[4]);
    if (!proof) throw std::runtime_error("cannot open proof output");
    proof << "pseudo-Boolean proof version 2.0\n";
    proof << "f " << formula_constraints << "\n";
    uint64_t proof_id = formula_constraints;
    std::vector<uint64_t> node_clause_ids(trace.nodes.size(), 0);
    for (size_t node_index = 0; node_index < trace.nodes.size();
         ++node_index) {
      const ProofNode& node = trace.nodes[node_index];
      uint64_t aggregate_id = 0;
      if (!node.partial.empty()) {
        proof << "pol ";
        size_t terms = 0;
        for (const auto [row, domain] : node.partial) {
          const uint64_t id =
              trace.support_ids.at({row, domain, node.next_row});
          proof << id << ' ';
          if (++terms > 1) proof << "+ ";
        }
        if (node.partial.size() > 1) {
          proof << atmost_id[node.next_row] << ' '
                << (node.partial.size() - 1) << " * + ";
        }
        proof << "\n";
        aggregate_id = ++proof_id;
      }
      proof << "rup ";
      for (const auto [row, domain] : node.partial) {
        proof << "+1 ~x" << (offset[row] + domain + 1) << ' ';
      }
      proof << ">= 1 ;";
      proof << ' ' << (node.partial.empty()
                           ? atleast_id[node.next_row]
                           : aggregate_id);
      for (size_t child : node.children) {
        if (child >= node_clause_ids.size() || node_clause_ids[child] == 0) {
          throw std::runtime_error("child proof ID is unavailable");
        }
        proof << ' ' << node_clause_ids[child];
      }
      proof << "\n";
      node_clause_ids[node_index] = ++proof_id;
    }
    proof << "output NONE\nconclusion UNSAT : " << proof_id
          << "\nend pseudo-Boolean proof\n";
    if (!proof) throw std::runtime_error("failed while writing proof");
    std::cout << "{\"source_index\":" << instance.source_index
              << ",\"rows\":" << rows.size()
              << ",\"variables\":" << variables
              << ",\"formula_constraints\":" << formula_constraints
              << ",\"support_constraints\":" << trace.support_ids.size()
              << ",\"proof_nodes\":" << trace.nodes.size()
              << ",\"proof_lines\":" << (proof_id - formula_constraints)
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
