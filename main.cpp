#include <QCoreApplication>
#include <QDebug>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include "kdtree.h"
#include "point.h"

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);
    //return a.exec();
#if 0
    PointV<2> pv3;
    pv3[0] = 2.5f;
    pv3[1] = 1.8f;

    qDebug() << pv3[2];

    KDTree<2> *myTree = new KDTree<2>();
    myTree->insert(pv3, 0);

    pv3[0] = 4; pv3[1] = 3;
    myTree->insert(pv3, 0);

    pv3[0] = 2; pv3[1] = 3;
    myTree->insert(pv3, 0);

    pv3[0] = 6; pv3[1] = 1;
    myTree->insert(pv3, 0);

    pv3[0] = 6; pv3[1] = 1;
    qDebug() << myTree->contains(pv3) << myTree->at(pv3);
#endif
    ifstream fs;
    fs.open("dump.txt");

    if(!fs.is_open())
        return -1;

    std::vector<PointV<2>> _positions;
    _positions.reserve(10000000);

    std::string line;
    size_t counter = 10;

    while(!fs.eof())// && counter-- > 0)
    {
        std::getline(fs, line);

        if(line.empty())
            continue;

        std::istringstream sstr(line);
        std::vector<std::string> tokens;

        std::copy(std::istream_iterator<std::string>(sstr),
                  std::istream_iterator<std::string>(),
                  back_inserter(tokens));

        PointV<2> newPoint;
        newPoint[0] = (float)atof(tokens[1].c_str());
        newPoint[1] = (float)atof(tokens[2].c_str());
        _positions.push_back(newPoint);
    }

    //Sort _positions based on x-coords using Lamda function
    std::sort(_positions.begin(), _positions.end(),
              [ ](const PointV<2>&a, const PointV<2>&b) -> bool
              { return (a[0] < b[0]); } );

    KDTree<2> *myTree = new KDTree<2>(_positions.begin(), _positions.end());
    qDebug() << myTree->size();

    PointV<2> min, max;
    min[0] = -125.;     max[0] = -105.0;
    min[1] =  75.0;     max[1] =  95.0;
    myTree->searchRange(min, max);

    _positions.clear();
    std::vector<PointV<2>>(_positions).swap(_positions);
    delete myTree;

    return 1;
}
