// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <iostream>
#include <stdexcept>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

static struct option options[] = {
  {"annotate", no_argument, nullptr, 'a'},
  {"kmer-len", required_argument, nullptr, 'k'},
  {"window-len", required_argument, nullptr, 'w'},
  {"frequency", required_argument, nullptr, 'f'},
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[merlion::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: merlion [options ...] <sequences> [<sequences> ...]\n"
      "\n"
      "  # default output is to stdout in JSON format\n"
      "  <sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -a, --annotate\n"
      "      use heuristics to find contained and chimeric sequences\n"
      "    -k, --kmer-len <int>\n"
      "      default: 15\n"
      "      length of minimizers used to find overlaps\n"
      "    -w, --window-len <int>\n"
      "      default: 5\n"
      "      length of sliding window from which minimizers are sampled\n"
      "    -f, --frequency <double>\n"
      "      default: 0.001\n"
      "      threshold for ignoring most frequent minimizers\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  bool annotate = false;

  std::uint8_t kmer_len = 15;
  std::uint8_t window_len = 5;
  double freq = 0.001;

  std::uint32_t num_threads = 1;

  std::string optstr = "ak:w:f:t:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'a': annotate = true; break;
      case 'k': kmer_len = std::atoi(optarg); break;
      case 'w': window_len = std::atoi(optarg); break;
      case 'f': freq = std::atof(optarg); break;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc) {
    std::cerr << "[merlion::] error: missing input file(s)!" << std::endl;
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  for (int i = optind; i < argc; ++i) {
    auto sparser = CreateParser(argv[i]);
    if (sparser == nullptr) {
      return 1;
    }

    decltype(sequences) chunk;
    try {
      chunk = sparser->Parse(-1);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << " (" << argv[i] << ")" << std::endl;
      return 1;
    }

    if (chunk.empty()) {
      std::cerr << "[merlion::] warning: file " << argv[i] << " is empty"
                << std::endl;
      continue;
    }

    sequences.insert(
        sequences.end(),
        std::make_move_iterator(chunk.begin()),
        std::make_move_iterator(chunk.end()));
  }
  if (sequences.empty()) {
    std::cerr << "[merlion::] error: empty sequences set!" << std::endl;
    return 1;
  }

  std::cerr << "[merlion::] loaded " << sequences.size() << " sequences "
            << std::fixed << timer.Stop() << "s"
            << std::endl;

  return 0;
}
