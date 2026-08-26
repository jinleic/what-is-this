#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
constexpr int N = 20;
constexpr int F = 5;

struct Domain {
  std::array<int8_t, N> half{};
  std::array<int8_t, N> row{};
};

struct Instance {
  std::array<uint8_t, N> masks{};
  std::array<uint8_t, N> internal{};
  std::array<std::vector<Domain>, N> domains;
  std::array<std::array<int, N>, N> target{};
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
  if (!stream) throw std::runtime_error("truncated W domain dump");
  return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) |
         (uint32_t(bytes[2]) << 8) | uint32_t(bytes[3]);
}

Instance load(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open W domain dump");
  char magic[5];
  stream.read(magic, 5);
  if (std::string(magic, 5) != "WDOM1") throw std::runtime_error("bad W domain magic");
  Instance instance;
  instance.source_index = read_u32(stream);
  stream.read(reinterpret_cast<char*>(instance.masks.data()), N);
  stream.read(reinterpret_cast<char*>(instance.internal.data()), N);
  if (!stream) throw std::runtime_error("truncated W support data");
  for (int row = 0; row < N; ++row) {
    const uint32_t count = read_u32(stream);
    instance.domains[row].reserve(count);
    for (uint32_t item = 0; item < count; ++item) {
      uint32_t code = read_u32(stream);
      Domain domain;
      for (int index = N - 1; index >= 0; --index) {
        domain.half[index] = int8_t(code % 3) - 1;
        code /= 3;
      }
      if (code != 0 || domain.half[row] != 0) throw std::runtime_error("bad W code");
      for (int index = 0; index < N; ++index) domain.row[index] = 2 * domain.half[index];
      domain.row[row] = 1 - 2 * instance.internal[row];
      instance.domains[row].push_back(domain);
    }
  }
  char extra;
  if (stream.read(&extra, 1)) throw std::runtime_error("trailing W domain bytes");
  for (int left = 0; left < N; ++left) {
    for (int right = 0; right < N; ++right) {
      int gram = 0;
      for (int vertex = 0; vertex < F; ++vertex) {
        const int a = 1 - 2 * ((instance.masks[left] >> vertex) & 1);
        const int b = 1 - 2 * ((instance.masks[right] >> vertex) & 1);
        gram += a * b;
      }
      instance.target[left][right] = 45 * (left == right) - 2 - 2 * gram;
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
  return dot == instance.target[left][right];
}

void verify_residual_closure(const Instance& instance,
                             const std::vector<int>& rows) {
  std::array<bool, N> active{};
  for (int row : rows) active[row] = true;
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right) {
      if (instance.masks[left] == instance.masks[right]
          && instance.internal[left] == instance.internal[right]
          && active[left] != active[right]) {
        throw std::runtime_error("rows are not residual-symmetry closed");
      }
    }
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
    std::cerr << "usage: involution_f5_w_dfs_certificate DUMP ROWS FORMULA PROOF\n";
    return 2;
  }
  try {
    const Instance instance = load(argv[1]);
    const std::vector<int> rows = parse_rows(argv[2]);
    verify_residual_closure(instance, rows);
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
    formula << "* exact W row-domain DFS certificate, source_index="
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
