# TSP_Solvers
Решение задачи коммивояжера несколькими способами. Полный перебор, Алгоритм ближайшего соседа, Муравьиный алгоритм.

- **Алгоритм полного перебора**: исследует все возможные маршруты для нахождения оптимального решения.
- **Алгоритм ближайшего соседа**: алгоритм, выбирающий следующий город с минимальным расстоянием от текущего.
- **Муравьиный алгоритм**: алгоритм оптимизации, вдохновленный поведением муравьев при поиске пути от места к месту. В алгоритме создается множество виртуальных муравьев, которые последовательно принимают решения о выборе следующего города для посещения на основе локальной информации о расстояниях и феромонах

## Сборка программы

g++ [PATH_TO_FILE] -o [PATH_TO_EXECUTABLE]
### Пример ввода

Входные данные представлены в виде матрицы смежности, в которую записаны расстояния между городами. 
Matrix[i][j] - расстояние между городами i и j.
```
int matrix[N][N] = {
      {0, 7, 8, 5, 5},
      {7, 0, 6, 2, 1},
      {8, 6, 0, 12, 8},
      {5, 2, 12, 0, 4},
      {5, 1, 8, 4, 0},
    }
```

### Пример вывода

```
Результат перебора:
Расстояние: 26
Маршрут: 12345

Муравьиный алгоритм:
Расстояние: 26
Маршрут: 12435

Алгоритм ближайшего соседа:
Расстояние: 26
Маршрут: 13524

```

