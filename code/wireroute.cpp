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

  // INCLUSIVE
  auto horizontal_line = [&](int x1, int x2, int y, int val) {
    if (x1 < 0 || x2 >= dim_x || x2 < 0 || x2 >= dim_x || y < 0 || y >= dim_y) {
      printf("Tried to draw horizontal line out of bounds! (%i->%i, %i) in (%i, %i)\n", x1, x2, y, dim_x, dim_y);
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
    if (y1 < 0 || y2 >= dim_y || y2 < 0 || y2 >= dim_y || x < 0 || x >= dim_x) {
        printf("Tried to draw vertical line out of bounds! (%i, %i->%i) in (%i, %i)\n", x, y1, y2, dim_x, dim_y);
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

  // To draw wires that have at least one bend
  auto draw_wire = [&](Wire wire, int val) {
    if (wire.start_y == wire.bend1_y) { // horizontal first
      if (wire.start_y == wire.end_y) {
          printf("Attempted to draw horizontal line... Please only use the draw_wire function for wires with at least one bend!");
          abort();
      }
      horizontal_line(wire.start_x, wire.bend1_x, wire.start_y, val);

      int line_start = wire.bend1_y + (wire.end_y > wire.start_y ? 1 : -1);
      vertical_line(line_start, wire.end_y, wire.bend1_x, val);

      if (wire.bend1_x != wire.end_x) {
        line_start = wire.bend1_x + (wire.end_x > wire.start_x ? 1 : -1);
        horizontal_line(line_start, wire.end_x, wire.end_y, val);
      }
    }

    else { assert(wire.start_x == wire.bend1_x); // vertical first
      if (wire.start_x == wire.end_x) {
          printf("Attempted to draw vertical line... Please only use the draw_wire function for wires with at least one bend!");
          abort();
      }
      vertical_line(wire.start_y, wire.bend1_y, wire.start_x, val);

      int line_start = wire.bend1_x + (wire.end_x > wire.start_x ? 1 : -1);
      horizontal_line(line_start, wire.end_x, wire.bend1_y, val);

      if (wire.bend1_y != wire.end_y) {
        line_start = wire.bend1_y + (wire.end_y > wire.start_y ? 1 : -1);
        vertical_line(line_start, wire.end_y, wire.end_x, val);
      }
    }
  };

  /********************************************************
   * the cost of adding `wire` to the occupancy matrix
   * with a potential bend location
   * (assumes the route is not already present,
   *  and has at least one bend)
   *
   * Travel from start->bend1->bend2->end,
   * calculating how much you *would* add to cost
   * by reading the occupancy matrix
   * (if o.m. position has value v, then you would add
   * ((v+1)^2 - v^2 = 2v+1)
   ********************************************************/
  auto route_cost_help = [&](Wire wire, int bend1_x, int bend1_y) {
    int cost = 0;

    // first move is horizontal
    if (wire.start_y == bend1_y) {
      if (wire.start_x < wire.end_x)
        for (int i = wire.start_x; i <= bend1_x; i++) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }
      else
        for (int i = wire.start_x; i >= bend1_x; i--) {
          cost += 2*(occupancy[wire.start_y][i]) + 1;
        }

      // First bend
      if (wire.start_y < wire.end_y)
        for (int i = bend1_y + 1; i <= wire.end_y; i++) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }
      else
        for (int i = bend1_y - 1; i >= wire.end_y; i--) {
          cost += 2*(occupancy[i][bend1_x]) + 1;
        }

      if (bend1_x == wire.end_x) return cost;

      // Second bend
      if (wire.start_x < wire.end_x)
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
      if (wire.start_y < wire.end_y)
        for (int i = wire.start_y; i <= bend1_y; i++) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }
      else
        for (int i = wire.start_y; i >= bend1_y; i--) {
          cost += 2*(occupancy[i][wire.start_x]) + 1;
        }

      // First bend
      if (wire.start_x < wire.end_x)
        for (int i = bend1_x + 1; i <= wire.end_x; i++) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }
      else
        for (int i = bend1_x - 1; i >= wire.end_x; i--) {
          cost += 2*(occupancy[bend1_y][i]) + 1;
        }

      if (bend1_y == wire.end_y) return cost;

      // Second bend
      if (wire.start_y < wire.end_y)
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

  auto route_cost_iteration = [&](Wire wire, int route_idx, int dx, int dy, int num_routes) {
    int bend1_x = wire.start_x;
    int bend1_y = wire.start_y;

    // Determine minimum horizontal-first path
    if (route_idx < abs(dx)) { // Determine minimum horizontal-first path
      int shift = route_idx + 1;
      int dir = dx >= 0 ? 1 : -1;
      bend1_x += shift * dir;

      return route_cost_help(wire, bend1_x, wire.start_y);
    }
    else if (route_idx < num_routes) { // 3) Determine minimum vertical-first path
      int shift = route_idx - abs(dx) + 1;
      int dir = dy >= 0 ? 1 : -1;
      bend1_y += shift * dir;

      return route_cost_help(wire, wire.start_x, bend1_y);
    }

    printf("Attempting to find cost of invalid route %i on wire with %i routes\n", route_idx, num_routes);
    abort();
  };

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
  auto compute_start = std::chrono::steady_clock::now();

  /**
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm.
   */
  if (parallel_mode == 'W') {
    compute_start = std::chrono::steady_clock::now();
    for (int iter = 0; iter < SA_iters; iter++) {
      // For each wire:

      auto iter_start = std::chrono::steady_clock::now();
      for (int w = 0; w < num_wires; w++) {
        auto wire = wires[w];

        if (iter > 0) { // remove route
          if (wire.start_x != wire.end_x && wire.start_y != wire.end_y) {
            // Remove the line:
            auto remove_start = std::chrono::steady_clock::now();
            draw_wire(wire, -1);
            const double remove_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - remove_start).count();
            //std::cout << "Removing Wire " << w << " time (sec):     " << remove_time << '\n';
          }

        }

        // OpenMP "Fork" to schedule
        // dx + dy threads to calculate minimum path
        int dx = wire.end_x - wire.start_x;
        int dy = wire.end_y - wire.start_y;

        if (dx == 0) {
          if (iter == 0) vertical_line(wire.start_y, wire.end_y, wire.start_x, 1);
          continue;
        }

        if (dy == 0) {
          if (iter == 0) horizontal_line(wire.start_x, wire.end_x, wire.start_y, 1);
          continue;
        }

        assert(dx != 0 && dy != 0 && "Either dx or dy is 0?!");
        assert(abs(dx) > 0 && "What? abs(dx) is <= 0?!");
        assert(abs(dy) > 0 && "What? abs(dy) is <= 0?!");

        unsigned int min_idx = 0;
        int min_cost = INT_MAX;

        int num_routes = abs(dx) + abs(dy);
        auto wire_start = std::chrono::steady_clock::now();
        //#pragma omp parallel num_threads(abs(dx)+abs(dy))
        #pragma omp parallel num_threads(num_threads)
        {
          //auto route_start = std::chrono::steady_clock::now();

          unsigned int thread_idx = omp_get_thread_num();
          unsigned int my_min_idx = thread_idx;
          int my_min_cost = INT_MAX;

          //
          for (unsigned int route_idx = thread_idx; route_idx < num_routes; route_idx += num_threads) {
            // 2) and 3) Determine all paths that travel horizontal and vertical
            int my_route_cost = INT_MAX;
            my_route_cost = route_cost_iteration(wire, route_idx, dx, dy, num_routes);

            if (my_route_cost < my_min_cost) {
              my_min_cost = my_route_cost;
              my_min_idx = route_idx;
            }
          }

          #pragma omp barrier

          #pragma omp critical
          {
            if (my_min_cost < min_cost) {
              min_cost = my_min_cost;
              min_idx = my_min_idx;
            }
          }
          /*
          const double route_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - route_start).count();*/
        }

        auto assign_start = std::chrono::steady_clock::now();
        if (min_idx < abs(dx)) {
          wire.bend1_x = wire.start_x + ((min_idx + 1) * (dx >= 0 ? 1 : -1));
          wire.bend1_y = wire.start_y;
        }
        else {
          wire.bend1_x = wire.start_x;
          wire.bend1_y = wire.start_y + (((min_idx - abs(dx)) + 1) * (dy >= 0 ? 1 : -1));
        }
        const double assign_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - assign_start).count();
        //std::cout << "Wire " << w << "/" << num_wires << " Assign time (sec):       " << assign_time << '\n';

        const double wire_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - wire_start).count();
        //std::cout << "Wire " << w << "/" << num_wires << " (" << abs(dx) + abs(dy) << " routes) time (sec): " << wire_time << '\n';
        // OpenMP "Join"

        // 4) With probability (1-P) (or straight line is min), choose minimum path.
        //    Otherwise, choose path randomly.
        auto add_start = std::chrono::steady_clock::now();
        if ((rand() / (float)(RAND_MAX)) <= SA_prob) {// Random Path
          if (rand() % 2 == 1) { // horizontal first
            wire.bend1_y = wire.start_y;
            wire.bend1_x = wire.start_x + ((rand() % (abs(wire.end_x - wire.start_x) + 1)) * (dx > 0 ? 1 : -1));
          }
          else { // vertical first
            wire.bend1_x = wire.start_x;
            wire.bend1_y = wire.start_y + ((rand() % (abs(wire.end_y - wire.start_y) + 1)) * (dy > 0 ? 1 : -1));
          }
        }
        draw_wire(wire, 1);
        const double add_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - add_start).count();
        //std::cout << "Adding wire " << w << "/" << num_wires << " time (sec):       " << add_time << '\n';
      }
      const double iter_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - iter_start).count();
      //std::cout << "Iteration " << iter << " time (sec): " << iter_time << '\n';
    }
  }
  else { assert(parallel_mode == 'A'); // parallel_mode == 'A'
    compute_start = std::chrono::steady_clock::now();
    for (int iter = 0; iter < SA_iters; iter++) {
      // Determine the batches: grab the wires in blocks and just assign them sequentially
      int B = 0; // wires already handled
      while (B < num_wires) {
        int best_route_idces[batch_size];
        int my_min_costs[batch_size];

        #pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int i = 0; i < batch_size; i++) { // Grab a batch and prepare the matrix
          my_min_costs[i] = INT_MAX;
          best_route_idces[i] = -1;
          if (i + B >= num_wires) continue;

          Wire wire;
          int dx, dy, num_routes;
          wire = wires[i + B];
          dx = wire.end_x - wire.start_x;
          dy = wire.end_y - wire.start_y;
          num_routes = abs(dx) + abs(dy);

          if (iter > 0) {
            if (dx != 0 && dy != 0) draw_wire(wire, -1);
            else if (dx == 0) vertical_line(wire.start_y, wire.end_y, wire.start_x, -1);
            else if (dy == 0) horizontal_line(wire.start_x, wire.end_x, wire.start_y, -1);
          }
        }

        #pragma omp parallel for schedule(semi-static) num_threads(num_threads)
        for (int i = 0; i < batch_size; i++) { // find the best route
          if (i + B >= num_wires) continue;

          auto wire = wires[i + B];
          int dx = wire.end_x - wire.start_x;
          int dy = wire.end_y - wire.start_y;
          int num_routes = abs(dx) + abs(dy);
          if (dx == 0) { // best route is horizontal
            best_route_idces[i] = num_routes;
            continue;
          }

          if (dy == 0) { // best route is vertical
            best_route_idces[i] = num_routes + 1;
            continue;
          }

          for (int route_idx = 0; route_idx < num_routes; route_idx++) {
            // determine new route for wire using another loop
            int my_route_cost = route_cost_iteration(wire, route_idx, dx, dy, num_routes);

            if (my_route_cost < my_min_costs[i]) {
              best_route_idces[i] = route_idx; // remember it!
              my_min_costs[i] = my_route_cost;
            }
          }
        }

        #pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int i = 0; i < batch_size; i++) {
          // update occupancy matrix
          // reader-writer locks are unnecessary
          if (i + B >= num_wires) continue;

          auto wire = wires[i + B];
          int dx = wire.end_x - wire.start_x;
          int dy = wire.end_y - wire.start_y;
          int num_routes = abs(dx) + abs(dy);

          if (best_route_idces[i] == num_routes) { // horizontal
            horizontal_line(wire.start_x, wire.end_x, wire.start_y, 1);
            continue;
          }

          if (best_route_idces[i] == num_routes + 1) { // vertical
            vertical_line(wire.start_y, wire.end_y, wire.start_x, 1);
            continue;
          }

          if (best_route_idces[i] < abs(dx)) {
            wire.bend1_x = wire.start_x + ((best_route_idces[i] + 1) * (dx >= 0 ? 1 : -1));
            wire.bend1_y = wire.start_y;
          }
          else if (best_route_idces[i] < num_routes) {
            wire.bend1_x = wire.start_x;
            wire.bend1_y = wire.start_y + ((best_route_idces[i] - abs(dx) + 1) * (dy >= 0 ? 1 : -1));
          }

          if ((rand() / (float)(RAND_MAX)) <= SA_prob) { // Random Path
            if (rand() % 2 == 1) { // horizontal first
              wire.bend1_y = wire.start_y;
              wire.bend1_x = wire.start_x + ((rand() % (abs(wire.end_x - wire.start_x) + 1)) * (dx > 0 ? 1 : -1));
            }
            else { // vertical first
              wire.bend1_x = wire.start_x;
              wire.bend1_y = wire.start_y + ((rand() % (abs(wire.end_y - wire.start_y) + 1)) * (dy > 0 ? 1 : -1));
            }
          }
          draw_wire(wire, 1);
        }

        B += batch_size;
      }
    }
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
