#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
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
    if (row < 0 || row >= N || seen[row]) throw std::runtime_error("invalid core rows");
    seen[row] = true;
    rows.push_back(row);
  }
  if (rows.empty()) throw std::runtime_error("empty core rows");
  return rows;
}

bool compatible(const Instance& instance, int left, size_t left_domain,
                int right, size_t right_domain) {
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
        throw std::runtime_error("core rows are not residual-symmetry closed");
      }
    }
  }
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: involution_f5_w_core_opb DOMAIN_DUMP CORE_ROWS OUTPUT\n";
    return 2;
  }
  try {
    const Instance instance = load(argv[1]);
    const std::vector<int> rows = parse_rows(argv[2]);
    verify_residual_closure(instance, rows);
    std::array<uint32_t, N> offset{};
    uint32_t variables = 0;
    for (int row : rows) {
      offset[row] = variables;
      variables += instance.domains[row].size();
    }
    uint64_t constraints = rows.size();
    struct Pair { int antecedent; int consequent; };
    std::vector<Pair> pairs;
    for (size_t i = 0; i < rows.size(); ++i) {
      for (size_t j = i + 1; j < rows.size(); ++j) {
        int left = rows[i], right = rows[j];
        if (instance.domains[left].size() > instance.domains[right].size()) {
          std::swap(left, right);
        }
        pairs.push_back({left, right});
        constraints += instance.domains[left].size();
      }
    }
    std::ofstream output(argv[3]);
    if (!output) throw std::runtime_error("cannot open OPB output");
    output << "* #variable= " << variables << " #constraint= " << constraints
           << " #equal= " << rows.size() << " intsize= 64\n";
    output << "* exact W row-domain core for source_index="
           << instance.source_index << "\n";
    output << "* no symmetry breaking and no projected counting\n";
    for (int row : rows) {
      output << "* exactly one W domain for row " << row << "\n";
      for (size_t domain = 0; domain < instance.domains[row].size(); ++domain) {
        output << "+1 x" << (offset[row] + domain + 1) << ' ';
      }
      output << "= 1;\n";
    }
    for (const Pair pair : pairs) {
      const auto& left_domains = instance.domains[pair.antecedent];
      const auto& right_domains = instance.domains[pair.consequent];
      for (size_t left = 0; left < left_domains.size(); ++left) {
        output << "-1 x" << (offset[pair.antecedent] + left + 1) << ' ';
        for (size_t right = 0; right < right_domains.size(); ++right) {
          if (compatible(instance, pair.antecedent, left,
                         pair.consequent, right)) {
            output << "+1 x" << (offset[pair.consequent] + right + 1) << ' ';
          }
        }
        output << ">= 0;\n";
      }
    }
    if (!output) throw std::runtime_error("failed while writing OPB output");
    std::cout << "{\"source_index\":" << instance.source_index
              << ",\"variables\":" << variables
              << ",\"constraints\":" << constraints
              << ",\"core_rows\":" << rows.size() << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
