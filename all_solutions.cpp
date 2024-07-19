#include <map>
#include <iostream>
#include <climits>
#include <stdlib.h>
#include <math.h>

#define min(X, Y)(((X) < (Y)) ? (X) : (Y))

#define N 5

class Graph {
public: int matrix[N][N];
      bool is_symmetric();

      std::map < int,
          int > storage;

      Graph();

      std::pair < int, int > bruteforce(int position = 0, int mask = 0);

      int nearest_neighbour(int position = 0, int mask = 0);
};

bool Graph::is_symmetric() {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

Graph::Graph() {}

std::pair < int, int > Graph::bruteforce(int position, int mask) {
    if ((mask | (1 << position)) == (1 << N) - 1)
        return std::make_pair(matrix[position][0], position + 1);

    auto min_route = std::make_pair(INT_MAX, 0);

    for (int i = 0; i < N; i++) {
        bool passed = (mask & (1 << i)) != 0;

        if (!passed && i != position) {

            auto route = bruteforce(i, mask | (1 << position));
            if (route.first + matrix[position][i] < min_route.first) {
                route.first += matrix[position][i];
                route.second *= 10;
                route.second += position + 1;
                min_route = route;
            }

        }
    }

    return min_route;
}

int Graph::nearest_neighbour(int position, int mask) {
    if ((mask | (1 << position)) == (1 << N) - 1) {
        return position + 1;
    }
    int nearest = INT_MAX, distance = INT_MAX;
    for (int i = 0; i < N; i++) {
        bool passed = (mask & (1 << i)) != 0;
        if (!passed && position != i && matrix[position][i] < distance) {
            nearest = i;
            distance = matrix[position][i];
        }
    }
    int route = nearest_neighbour(nearest, mask | (1 << position));
    route *= 10;
    route += position + 1;
    return route;
}

#define ALPHA 0
#define BETA 0
#define EVAPORATION_INDEX 0.5
#define Q 100

#define ANT_COUNT 10
#define ITERATIONS 100

class AntGraph : public Graph {
private: int citiesCount;
       float** distances;
       float** pheromones;
       int alphaMod;
       int betaMod;
       int alphaBest;
       int betaBest;

       void initialize_pheromones();

       float ant_product(int cityX, int cityY);

       float select_next_city(int currentCity, int* visitedCities);

       int* generate_solution();

       int calculate_distance(int* solution);

       void update_pheromones();

       int* solveTSP();
public: int* ant_alhorithm();
};

void AntGraph::initialize_pheromones() {
    pheromones = (float**)malloc(N * sizeof(float*));
    if (pheromones == nullptr) {
        throw "No memory left";
    }
    for (int i = 0; i < N; i++) {
        pheromones[i] = (float*)malloc(N * sizeof(float));
        if (pheromones[i] == nullptr) {
            throw "No memory left";
        }
        for (int j = 0; j < N; j++) {
            pheromones[i][j] = 1.0;
        }
    }
}

float AntGraph::ant_product(int city1, int city2) {
    return 1.0 / matrix[city1][city2];
}

float AntGraph::select_next_city(int currentCity, int* visitedCities) {
    float probability = 0.0;
    for (int i = 0; i < N; i++) {
        if (visitedCities[i] == 0) {
            probability += pow(pheromones[currentCity][i], (ALPHA + alphaMod)) *
                pow(ant_product(currentCity, i), (BETA + betaMod));
        }
    }
    float random = ((float)rand() / RAND_MAX);
    float cumulativeProbability = 0.0;
    for (int i = 0; i < N; i++) {
        if (visitedCities[i] == 0) {
            cumulativeProbability += pow(pheromones[currentCity][i], (ALPHA + alphaMod)) *
                pow(ant_product(currentCity, i), (BETA + betaMod));
            if ((cumulativeProbability / probability) >= random) {
                return i;
            }
        }
    }
    throw "Something went wrong, cannot be here!";
}

int* AntGraph::generate_solution() {
    int* visitedCities = (int*)malloc(N * sizeof(int));
    if (visitedCities == nullptr) {
        throw "No memory left";
    }
    int* solution = (int*)malloc(N * sizeof(int));
    if (solution == nullptr) {
        throw "No memory left";
    }
    int i;
    for (i = 0; i < N; i++) {
        visitedCities[i] = 0;
    }
    solution[0] = 0;
    visitedCities[0] = 1;
    for (i = 1; i < N; i++) {
        solution[i] = select_next_city(solution[i - 1], visitedCities);
        visitedCities[solution[i]] = 1;
    }
    free(visitedCities);
    return solution;
}

int AntGraph::calculate_distance(int* solution) {
    int distance = 0.0;
    for (int i = 0; i < N - 1; i++) {
        distance += matrix[solution[i]][solution[i + 1]];
    }
    distance += matrix[solution[N - 1]][solution[0]];
    return distance;
}

void AntGraph::update_pheromones() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            pheromones[i][j] *= (1.0 - EVAPORATION_INDEX);
        }
    }
    for (i = 0; i < ANT_COUNT; i++) {
        int* solution = generate_solution();
        float distance = calculate_distance(solution);
        for (j = 0; j < N - 1; j++) {
            pheromones[solution[j]][solution[j + 1]] += (Q / distance);
        }
        pheromones[solution[N - 1]][solution[0]] += (Q / distance);
        free(solution);
    }
}

int* AntGraph::solveTSP() {
    initialize_pheromones();
    int i;
    for (i = 0; i < ITERATIONS; i++) {
        update_pheromones();
    }
    int* bestSolution = generate_solution();
    return bestSolution;
}

int* AntGraph::ant_alhorithm() {
    alphaMod = 9;
    betaMod = 0;
    int* solution;
    float best = 999999999.9;
    int current;
    int* solutionBest = solveTSP();
    for (int i = 0; i < 10; i++) {
        solution = solveTSP();
        current = calculate_distance(solution);
        if (current < best) {
            best = current;
            alphaBest = alphaMod;
            betaBest = betaMod;
            for (int j = 0; j < N; j++) {
                solutionBest[j] = solution[j] + 1;
            }
        }
        alphaMod--;
        betaMod++;
    }
    return solutionBest;
}

int main() {
    setlocale(LC_ALL, "Russian");

    int matrix[N][N] = {
      {0, 7, 8, 5, 5},
      {7, 0, 6, 2, 1},
      {8, 6, 0, 12, 8},
      {5, 2, 12, 0, 4},
      {5, 1, 8, 4, 0},
    };

    AntGraph graph = AntGraph();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            graph.matrix[i][j] = matrix[i][j];
    if (!graph.is_symmetric()) {
        throw "Матрица несимметричная";
        return 1;
    }

    auto routeb = graph.bruteforce();

    std::cout << "Результат перебора:" << std::endl << "Расстояние: " << routeb.first << std::endl << "Маршрут: " << routeb.second << std::endl << std::endl;

    int* route = graph.ant_alhorithm();
    int distance = 0;

    for (int i = 0; i < N - 1; i++) {
        distance += matrix[route[i] - 1][route[i + 1] - 1];
    }
    distance += matrix[route[N - 1] - 1][0];

    std::cout << "Муравьиный алгоритм:" << std::endl << "Расстояние: " << distance << std::endl << "Маршрут: ";
    for (int i = 0; i < N; i++) {
        std::cout << route[i];
    }
    std::cout << std::endl << std::endl;

    distance = graph.nearest_neighbour();
    for (int i = 0; i < N; i++) {
        route[i] = distance % 10;
        distance /= 10;
    }
    distance = 0;
    for (int i = 0; i < N - 1; i++) {
        distance += matrix[route[i] - 1][route[i + 1] - 1];
    }
    distance += matrix[route[N - 1] - 1][0];

    std::cout << "Алгоритм ближайшего соседа:" << std::endl << "Расстояние: " << distance << std::endl << "Маршрут: ";
    for (int i = 0; i < N; i++) {
        std::cout << route[i];
    }
    std::cout << std::endl << std::endl;
    return 0;
}
