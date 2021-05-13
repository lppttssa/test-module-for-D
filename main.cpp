#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <time.h>
#include <set>

using namespace std;

struct coordStruct {
    int x;
    int y;
};

struct part {
    int vert;
    double weight;
};

struct newVertPlace {
    int vert;
    int x;
    int y;
};

struct newEdgesStruct {
    int start;
    int end;
    double weight;
};

vector<int> parents(0);
vector<int> path(0);

int getMinVert(vector <vector <part>> list_, vector <coordStruct> coord_, int from, int to) {
    double minPath = INT16_MAX;
    int nextVert = 1;
    for (auto j: list_[from]) {
        if (minPath > j.weight + sqrt(pow(coord_[j.vert].x - coord_[to].x, 2) + pow(coord_[j.vert].y - coord_[to].y, 2))) {
            minPath = j.weight + sqrt(pow(coord_[j.vert].x - coord_[to].x, 2) + pow(coord_[j.vert].y - coord_[to].y, 2));
            nextVert = j.vert;
        }
    }
    return nextVert;
}

void getPath(int from, int to) {
    path.clear();
    path.resize(0);
    int j = to;
    while (j != from) {
        path.push_back(parents[j]);
        j = parents[j];
    }
    vector<int> tmp(0);
    for (int i = path.size() - 1; i >= 0; i--) {
        tmp.push_back(path[i]);
    }
    path = tmp;
}

double Dijkstra(vector <vector <part>> list_, int fromVert, int endVert, int vertNum) {

    parents.clear();
    parents.resize(vertNum);
    vector <double> distDijkstra(vertNum+1);
    for (int i = 0; i <= vertNum; i++) {
        distDijkstra[i] = INT_MAX;
    }
    for (int i = 0; i < vertNum; i++) {
        parents[i] = -1;
    }

    distDijkstra[fromVert] = 0;
    set <pair <int, double> > s;
    s.insert(pair<int, double>(fromVert, 0));

    while (!s.empty()) {
        int from = s.begin()->first;
        s.erase(s.begin());

        for (int i = 0; i < list_[from].size(); i++) {
            part to = list_[from][i];

            if (distDijkstra[from] + to.weight < distDijkstra[to.vert]) {
                s.erase(pair<int, double>(to.vert, to.weight));
                parents[to.vert] = from;
                distDijkstra[to.vert] = distDijkstra[from] + to.weight;
                s.insert(pair<int, double>(to.vert, distDijkstra[to.vert]));
            }
        }
    }
    if (distDijkstra[endVert] == INT_MAX) {
        return(-1);
    }
    else {
        return(distDijkstra[endVert]);
    }
}

//input : vert number, edge number, from, to
//vert number lines with x, y coordinates
//edge number lines with a vert, b vert and edge value

//if changed
//Y
//changed vert number, changed edge number
//changed vert number with vert, x, y
//changed edge number with start, end, changed edge value


int main() {
    int minN = 500;
    int maxN = 500;
    int max = 0;
    int num = 0;


    srand(9);
    for (int i = 0; i < 100; ++i) {

        int testNum = i;
        string testName = "input" + to_string(testNum) + ".txt";

        ofstream out("D:\\6sem\\algorithm\\testDirective\\tests100\\" + testName);


//    Получаем minN <= nVertex <= maxN
        int nVertex = minN + rand() % (maxN - minN + 1);


        int minM = 3;
        int maxM = (nVertex * (nVertex - 1)) / 2;
//    Получаем minM <= nEdges <= maxM
        int nEdges = minM + rand() % (maxM - minM + 1);


        vector <vector <part>> list(nVertex);
        vector <coordStruct> coord(0);


        int minC = 0;
        int maxC = 50;
        for (int i = 0; i < nVertex; ++i) {
            bool isContinued = false;
            //Добавляем координаты вершие в массив координат
            int x = minC + rand() % (maxC - minC + 1);
            int y = minC + rand() % (maxC - minC + 1);
            for (auto i: coord) {
                if ((i.x == x) && (i.y == y)) {
                    isContinued = true;
                }
            }
            if (isContinued) {
                i--;
                continue;
            }
            coordStruct mem{x, y};
            coord.push_back(mem);
        }

        for (int i = 0; i < nEdges; ++i) {
            bool isContinued = false;
//        Генерируем случайные рёбра для графа и высчитываем их веса, используя массив координат
            int from = rand() % nVertex;
            int to = rand() % nVertex;

//        Тут можешь проверять, можно ли такое ребро добавлять, не сломает ли оно граф.
//        Я вот проверил, чтобы не было петель:
            if (from == to) {
//            Херовый тест получился, потому что ребро ведёт само в себя.
                isContinued = true;
            }
            for (auto j: list[from]) {
                if (j.vert == to) {
                    isContinued = true;
                }
            }
            for (auto j: list[to]) {
                if (j.vert == from) {
                    isContinued = true;
                }
            }
            if (isContinued) {
                i--;
                continue;
            }

            part mem{};
            mem.vert = to;
            mem.weight = sqrt(pow(coord[from].x - coord[to].x, 2) + pow(coord[from].y - coord[to].y, 2));
            list[from].push_back(mem);
        }



        if (Dijkstra(list, 0, nVertex-1, nVertex) == -1) {

            i--;
            continue;
        }



        out << nVertex << " " << nEdges << " " << 0 << " " << nVertex-1 << "\n";
        for (auto i : coord) {
            out << i.x << " " << i.y << "\n";
        }
        for (int i = 0; i < nVertex; i++) {
            for (auto j : list[i]) {
                out << i << " " << j.vert << " " << j.weight << "\n";
            }
        }


        getPath(0, nVertex-1);



        //Генерируем часть с изменением графа///
        int changes = path.size();
        for (int i = 1; i < path.size(); ++i) {
            vector<newVertPlace> changedVerts(0);
            vector<newEdgesStruct> newEdges(0);
            bool isContinued = false;
            int changedVert = path[i];
            int newX = minC + rand() % (maxC - minC + 1);
            int newY = minC + rand() % (maxC - minC + 1);
            for (auto i: coord) {
                if ((i.x == newX) && (i.y == newY)) {
                    isContinued = true;
                }
            }
            if (isContinued) {
                i--;
                continue;
            }
            //coord[changedVert].x = newX;
            //coord[changedVert].y = newY;
            changedVerts.push_back({changedVert, coord[changedVert].x, coord[changedVert].y});


            for (auto k: list[changedVert]) {
                double minWeight = sqrt(pow(coord[changedVert].x - coord[k.vert].x, 2) + pow(coord[changedVert].y - coord[k.vert].y, 2));
                double maxWeight = minWeight * 2;
                double f = (double)rand() / RAND_MAX;
                double changedWeight = minWeight + f * (maxWeight - minWeight);
                k.weight = changedWeight;
                newEdges.push_back({changedVert, k.vert, changedWeight});
            }

            /*for (int j = 0; j < list.size(); j++) {
                for (auto k: list[j]) {
                    if (j == changedVert) {
                        double changedWeight = sqrt(pow(coord[j].x - coord[k.vert].x, 2) + pow(coord[j].y - coord[k.vert].y, 2));
                        k.weight = changedWeight;
                        newEdges.push_back({j, k.vert, changedWeight});
                    }
                    if (k.vert == changedVert) {
                        double changedWeight = sqrt(pow(coord[j].x - coord[k.vert].x, 2) + pow(coord[j].y - coord[k.vert].y, 2));
                        k.weight = changedWeight;
                        newEdges.push_back({j, k.vert, changedWeight});
                    }
                }
            }*/

            out << "Y\n" << changedVerts.size() << " " << newEdges.size() << "\n";
            for (auto j: changedVerts) {
                out << j.vert << " " << j.x << " " << j.y << "\n";
            }
            for (auto j: newEdges) {
                out << j.start << " " << j.end << " " << j.weight << "\n";
            }

            //testNum++;

            double res = Dijkstra(list, 0, nVertex - 1, nVertex);
            getPath(0, nVertex - 1);

        }

        cout << testNum << " " << path.size() << "\n";


        /*int minChangedVert = 1;
        int maxChangedVert = 1;
//    Находим количество измененных вершин
        int changedVertNum = minChangedVert + rand() % (maxChangedVert - minChangedVert + 1);
        //Генерируем номера измененных вершин и их координаты
        vector<newVertPlace> changedVerts(0);
        vector<newEdgesStruct> newEdges(0);
        for (int i = 0; i < changedVertNum; ++i) {
            bool isContinued = false;
            int minVert = 0;
            int maxVert = nVertex;
            //Генерируем номер изменяемой вершины
            //int changedVert = minVert + rand() % (maxVert - minVert + 1);
            //int changedVert = getMinVert(list, coord, 0, nVertex - 1);
            int changedVert = path[1];
            //Генерируем новые координаты для измененной вершины
            int newX = minC + rand() % (maxC - minC + 1);
            int newY = minC + rand() % (maxC - minC + 1);
            for (auto i: coord) {
                if ((i.x == newX) && (i.y == newY)) {
                    isContinued = true;
                }
            }
            if (isContinued) {
                i--;
                continue;
            }
            coord[changedVert].x = newX;
            coord[changedVert].y = newY;
            changedVerts.push_back({changedVert, newX, newY});
            //Высчитываем новые веса для всех ребер, входящих в и выходящих из измененной вершины

            for (int j = 0; j < list.size(); j++) {
                for (auto k: list[j]) {
                    if (j == changedVert) {
                        double changedWeight = sqrt(pow(coord[j].x - coord[k.vert].x, 2) + pow(coord[j].y - coord[k.vert].y, 2));
                        newEdges.push_back({j, k.vert, changedWeight});
                    }
                    if (k.vert == changedVert) {
                        double changedWeight = sqrt(pow(coord[j].x - coord[k.vert].x, 2) + pow(coord[j].y - coord[k.vert].y, 2));
                        newEdges.push_back({j, k.vert, changedWeight});
                    }
                }
            }
        }
        out << "Y\n" << changedVerts.size() << " " << newEdges.size() << "\n";
        for (auto j: changedVerts) {
            out << j.vert << " " << j.x << " " << j.y << "\n";
        }
        for (auto j: newEdges) {
            out << j.start << " " << j.end << " " << j.weight << "\n";
        }

        testNum++;*/
        //cout << testNum << " \n";
        if (max < path.size()) {
            max = path.size();
            num = testNum;
        }
    }

    cout << max << " " << num;
}
