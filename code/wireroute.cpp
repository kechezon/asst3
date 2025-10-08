/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <limits.h>

#include <unistd.h>
#include <omp.h>

void print_stats(const std::vector<std::vector<int>>& occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto& row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

void write_output(const std::vector<Wire>& wires, const int num_wires, const std::vector<std::vector<int>>& occupancy, const int dim_x, const int dim_y, const int num_threads, std::string input_filename) {
  if (std::size(input_filename) >= 4 && input_filename.substr(std::size(input_filename) - 4) == ".txt") {
    input_filename.resize(std::size(input_filename) - 4);
  }

  const std::string occupancy_filename = input_filename + "_occupancy_" + std::to_string(num_threads) + ".txt";
  const std::string wires_filename = input_filename + "_wires_" + std::to_string(num_threads) + ".txt";

  std::ofstream out_occupancy(occupancy_filename, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_occupancy << dim_x << ' ' << dim_y << '\n';
  for (const auto& row : occupancy) {
    for (const int count : row) {
      out_occupancy << count << ' ';
    }
    out_occupancy << '\n';
  }

  out_occupancy.close();

  std::ofstream out_wires(wires_filename, std::fstream:: out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n' << num_wires << '\n';

  for (const auto& [start_x, start_y, end_x, end_y, bend1_x, bend1_y] : wires) {
    out_wires << start_x << ' ' << start_y << ' ' << bend1_x << ' ' << bend1_y << ' ';

    if (start_y == bend1_y) {
    // first bend was horizontal

      if (end_x != bend1_x) {
        // two bends

        out_wires << bend1_x << ' ' << end_y << ' ';
      }
    } else if (start_x == bend1_x) {
      // first bend was vertical

      if(end_y != bend1_y) {
        // two bends

        out_wires << end_x << ' ' << bend1_y << ' ';
      }
    }
    out_wires << end_x << ' ' << end_y << '\n';
  }

  out_wires.close();
}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'n':
        num_threads = atoi(optarg);
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'm':
        parallel_mode = *optarg;
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 || (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);

  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x));

  for (auto& wire : wires) {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
  }

  /* Initialize any additional data structures needed in the algorithm */

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  // INCLUSIVE
  auto horizontal_line = [&](int x1, int x2, int y, int val) {
    if (x1 < 0 || x2 >= dim_x || x2 < 0 || x2 >= dim_x || y < 0 || y >= dim_y) {
      printf("Tried to draw horizontal line out of bounds! (%i->%i, %i)\n", x1, x2, y);
      abort();
    }

    if (x1 <= x2) {
      for (int i = x1; i <= x2; i++) occupancy[y][i] += val;
    }
    else {
      for (int i = x1; i >= x2; i--) occupancy[y][i] += val;
    }

    return;
  };

  // INCLUSIVE
  auto vertical_line = [&](int y1, int y2, int x, int val) {
    if (y1 < 0 || y2 >= dim_y || y2 < 0 || y2 >= dim_x || x < 0 || x >= dim_x) {
        printf("Tried to draw horizontal line out of bounds! (%i, %i->%i)\n", x, y1, y2);
        abort();
    }

    if (y1 <= y2) {
      for (int i = y1; i <= y2; i++) occupancy[i][x] += val;
    }
    else {
      for (int i = y1; i >= y2; i--) occupancy[i][x] += val;
    }

    return;
  };

  auto overall_cost = [&]() {
    int cost = 0;
    for (const auto& row : occupancy) {
      for (const int count : row) {
        cost += count * count;
      }
    }
    return cost;
  };

  /********************************************************
   * the cost of adding `wire` to the occupancy matrix
   * with a potential bend location
   * (assumes the route is not already present)
   *
   * Travel from start->bend1->bend2->end,
   * calculating how much you *would* add to cost
   * by reading the occupancy matrix
   * (if o.m. position has value v, then you would add
   * ((v+1)^2 - v^2 = 2v+1)
   ********************************************************/
  auto route_cost = [&](Wire wire, int bend1_x, int bend1_y) {
    int cost = 0;

    // first move is horizontal
    if (wire.start_y == bend1_y) {
      if (wire.start_x <= bend1_x)
        for (int i = wire.start_x; i <= bend1_x; i++) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }
      else
        for (int i = wire.start_x; i >= bend1_x; i--) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }

      // First bend
      if (wire.bend1_y <= wire.end_y)
        for (int i = bend1_y + 1; i <= wire.end_y; i++) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }
      else
        for (int i = bend1_y - 1; i >= wire.end_y; i--) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }

      // Second bend
      if (wire.bend1_x <= wire.end_x)
        for (int i = bend1_x + 1; i <= wire.end_x; i++) {
          cost += 2*(occupancy[wire.end_y][i]) + 1;
        }
      else
        for (int i = bend1_x - 1; i >= wire.end_x; i--) {
          cost += 2*(occupancy[wire.end_y][i]) + 1;
        }

      return cost;
    }

    // first move is vertical
    else {
      if (wire.start_y <= bend1_y)
        for (int i = wire.start_y; i <= bend1_y; i++) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }
      else
        for (int i = wire.start_y; i >= bend1_y; i--) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }

      // First bend
      if (bend1_x <= wire.end_x)
        for (int i = bend1_x + 1; i <= wire.end_x; i++) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }
      else
        for (int i = bend1_x - 1; i >= wire.end_x; i--) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }

      // Second bend
      if (wire.bend1_y <= wire.end_y)
        for (int i = bend1_y + 1; i <= wire.end_y; i++) {
          cost += 2*(occupancy[i][wire.end_x]) + 1;
        }
      else
        for (int i = bend1_y - 1; i >= wire.end_y; i--) {
          cost += 2*(occupancy[i][wire.end_x]) + 1;
        }

      return cost;
    }
  };

  /**
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm.
   * TODO
   */
  if (parallel_mode == 'W') {

    // First pass:
    // 1) Add random route for each wire and draw it
    for (auto& wire : wires) { // TODO: OpenMP this to fork based on wire
      // Chose random bend axis and location
      if (rand() % 2 == 1) { // horizontal first
        wire.bend1_y = wire.start_y;
        wire.bend1_x = wire.start_x + (rand() % (abs(wire.end_x - wire.start_x) + 1));

        horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, 1);
        vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
        horizontal_line(wire.bend1_x + 1, wire.end_x, wire.end_y, 1);
      }
      else { // vertical first
        wire.bend1_x = wire.start_x;
        wire.bend1_y = wire.start_y + (rand() % (abs(wire.end_y - wire.start_y) + 1));

        vertical_line(wire.start_y, wire.bend1_y, wire.start_x, 1);
        horizontal_line(wire.bend1_x + 1, wire.end_x, wire.bend1_y, 1);
        vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
      }
    }
    // TODO: OpenMP Join


    // 2) Second pass onwards:
    for (int iter = 0; iter < SA_iters; iter++) {

      // For each wire:
      for (auto &wire : wires) {
        int min_path_bend1_x = 0;
        int min_path_bend1_y = 0;
        // a) Compute overall cost in current occupancy matrix
        //    with and without wire
        int cost_with = overall_cost();

        // Remove the line:
        if (wire.start_y == wire.bend1_y) { // horizontal first
          horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, -1);
          vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, -1);
          horizontal_line(wire.bend1_x + 1, wire.end_x, wire.end_y, -1);

        }
        else { // vertical first
          vertical_line(wire.start_y, wire.bend1_y, wire.start_x, -1);
          horizontal_line(wire.bend1_x + 1, wire.end_x, wire.bend1_y, -1);
          vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, -1);
        }
        int cost_without = overall_cost();

        int min_cost = INT_MAX;

        // OpenMP "Fork" to schedule
        // dx + dy threadto calculate minimum path
        int dx = wire.end_x - wire.start_x;
        int dy = wire.end_y - wire.start_y;

        omp_lock_t *cost_lock;
        omp_init_lock(cost_lock);
        #pragma omp parallel num_threads(abs(dx)+abs(dy))
        {
          unsigned int thread_idx = omp_get_thread_num();
          int my_cost = 0;
          int bend1_x = wire.start_x;
          int bend1_y = wire.start_y;
          // b) Determine minimum horizontal-first path
          if (thread_idx < abs(dx)) { // Determine minimum horizontal-first path
            int shift = thread_idx + 1;
            int dir = dx >= 0 ? 1 : -1;
            bend1_x += shift * dir;

            my_cost = route_cost(wire, bend1_x, wire.start_y);
          }
          else { // c) Determine minimum vertical-first path
            int shift = thread_idx - abs(dx) + 1;
            int dir = dy >= 0 ? 1 : -1;
            bend1_y += shift * dir;

            my_cost = route_cost(wire, wire.start_x, bend1_y);
          }

          omp_set_lock(cost_lock);
          if (my_cost < min_cost) {
              min_cost = my_cost;
              wire.bend1_x = bend1_x;
              wire.bend1_y = bend1_y;
          }
          omp_unset_lock(cost_lock);
        }
        // OpenMP "Join"

        // d) With probability (1-P) (or straight line is min), choose minimum path.
        //    Otherwise, choose path randomly.
        if ((rand() / (float)(RAND_MAX)) <= 1 - SA_prob) {
          // Add to occupancy matrix
          if (wire.start_y == wire.bend1_y) { // horizontal first
            horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, 1);
            vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
            horizontal_line(wire.bend1_x + 1, wire.end_x, wire.end_y, 1);
          }
          else { // vertical first
            vertical_line(wire.start_y, wire.bend1_y, wire.start_x, 1);
            horizontal_line(wire.bend1_x + 1, wire.end_x, wire.bend1_y, 1);
            vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
          }
        }
        else {
          // TODO: Random Path
          if (rand() % 2 == 1) { // horizontal first
            wire.bend1_y = wire.start_y;
            wire.bend1_x = wire.start_x + (rand() % (abs(wire.end_x - wire.start_x) + 1));

            horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, 1);
            vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
            horizontal_line(wire.bend1_x + 1, wire.end_x, wire.end_y, 1);
          }
          else { // vertical first
            wire.bend1_x = wire.start_x;
            wire.bend1_y = wire.start_y + (rand() % (abs(wire.end_y - wire.start_y) + 1));

            vertical_line(wire.start_y, wire.bend1_y, wire.start_x, 1);
            horizontal_line(wire.bend1_x + 1, wire.end_x, wire.bend1_y, 1);
            vertical_line(wire.bend1_y + 1, wire.end_y, wire.bend1_x, 1);
          }
        }
      }
    }
  }
  else { assert(parallel_mode == 'A'); // parallel_mode == 'A'
    // TODO: What can be done independently betwen wires?
  }

  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* Write wires and occupancy matrix to files */

  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
}

validate_wire_t Wire::to_validate_format(void) const {
  /* TODO(student): Implement this if you want to use the wr_checker. */
  /* See wireroute.h for details on validate_wire_t. */
  throw std::logic_error("to_validate_format not implemented.");
}
