#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <vector>

#define PI 3.14
#define EARTHRADIUSINKM 6356.752
#define MAXRANGE 320 // range in KM
#define CARSPEED 105.f //spedd in KM
#define CHARGESPEED 800 //speed in KM/hr

//Helpers
double ToDeg(double Rad)
{
    return Rad * (180.f / PI);
}

double ToRad(double Deg)
{
    return Deg * (PI / 180.f);
}

double RoundTwoDecimals(double val)
{
    double res = static_cast<double>(val * 100 + 0.5);
    return res / 100;
}

struct GeoCoordinates
{
    double Latitude = 0.0;
    double Longitude = 0.0;
};

double GetDistanceFromCoordinates(const GeoCoordinates& i_Coords1, const GeoCoordinates& i_Coords2)
{

    long double lat1 = ToRad(i_Coords1.Latitude);
    long double lat2 = ToRad(i_Coords2.Latitude);
    long double long1 = ToRad(i_Coords1.Longitude);
    long double long2 = ToRad(i_Coords2.Longitude);

    long double dLat =  lat2 - lat1;
    long double dLong = long2 - long1;

    long double res = std::pow(std::sin(dLat / 2), 2) + std::cos(lat1) * std::cos(lat2) * std::pow(sin(dLong / 2), 2);

    res = 2 * asin(std::sqrt(res));
    
    res = res * EARTHRADIUSINKM;
    return res;
}

//Graph
struct Edge
{
    int start = 0;
    int sink = 0;
    float cost = 0;
    Edge(int begin, int end, float timeCost)
        :start(begin)
        ,sink(end)
        ,cost(timeCost)
    {}
};

struct  Node
{
    std::vector<Edge*> Connections;
    GeoCoordinates Coord;
    std::string Name;
    int Id;
    Node() = default;
    Node(int id)
        : Id(id)
    {}
};

//Graph
class NetworkGraph
{
public:
    NetworkGraph(const char* NetworkFilePath)
    {
        if (IsFileFormatSupported(NetworkFilePath))
        {
            std::fstream file;
            file.open(NetworkFilePath, std::ios::in);
            if (!file.is_open())
            {
                std::cout << "Error: File Could not be opened" << std::endl;
                std::cout << "FilePath:" << NetworkFilePath << std::endl;
            }
            else
            {
                std::string line, col;
                while (std::getline(file, line))
                {
                    std::stringstream ss(line);
                    int i = 0;
                    std::string NodeName;
                    GeoCoordinates NodeCoord;
                    while (std::getline(ss, col, ','))
                    {
                        if (i == 0)
                        {
                            NodeName = col;
                        }
                        else if (i == 1)
                        {
                            NodeCoord.Latitude = std::stod(col);
                        }
                        else if (i == 2)
                        {
                            NodeCoord.Longitude = std::stod(col);
                        }
                        ++i;
                    }
                    auto NodeItr = GraphMap.emplace(++ChargersCount, ChargersCount);
                    NodeItr.first->second.Coord = NodeCoord;
                    NodeItr.first->second.Name = NodeName;
                    GraphNameMap.emplace(NodeName, NodeItr.first->second.Id);
                }
                file.close();
            }
        }
        else
        {
            std::cout << "Error: FileFormat Not Supported" << std::endl;
        }
    }

    bool IsValid() const { return GraphMap.size() > 0; }

    void SetupConnections()
    {
        if (IsValid())
        {
            for (auto it = GraphMap.begin(); it != GraphMap.end(); ++it)
            {
                for (auto it2 = GraphMap.begin(); it2 != GraphMap.end(); ++it2)
                {
                    if (it->first == it2->first)
                        continue;
                    long double dis = GetDistanceFromCoordinates(it->second.Coord, it2->second.Coord);
                    if (dis < MAXRANGE)
                    {
                        float timeCost = dis / CARSPEED;
                        it->second.Connections.emplace_back(new Edge(it->first, it2->first, timeCost));
                    }

                }
            }
        }
    }
    std::unordered_map<int, Node> GraphMap;
    std::unordered_map<std::string,int> GraphNameMap;
private:
    bool IsFileFormatSupported(const char* FilePath)
    {
        std::string File(FilePath);
        size_t pos = File.find(".csv");
        return pos != std::string::npos;
    }
    static int ChargersCount;
};

int NetworkGraph::ChargersCount = 0;

//PathFinding
struct  NodeRecord
{
    int Node;
    Edge* IncomingEdge;
    float CostSoFar;
};

inline bool operator<(const NodeRecord& left, const NodeRecord& right)
{
    return left.CostSoFar < right.CostSoFar;
}

inline bool operator>(const NodeRecord& left, const NodeRecord& right)
{
    return left.CostSoFar > right.CostSoFar;
}

inline bool operator==(const NodeRecord& left, const NodeRecord& right)
{
    return left.Node == right.Node;
}

inline  bool operator!=(const NodeRecord& left, const NodeRecord& right)
{
    return !(left == right);
}

class PathFinder
{
public:
    PathFinder(NetworkGraph* i_GS)
        :GraphSpace(i_GS)
    {}

    std::vector<Edge*> FindPath(const std::string& source, const std::string& sink)
    {
        if (GraphSpace)
        {
            auto& Map = GraphSpace->GraphNameMap;
            if (Map.find(source) != Map.end() && Map.find(sink) != Map.end())
            {
                return FindPath(Map[source], Map[sink]);
            }
        }
        return std::vector<Edge*>();
    }

    std::vector<Edge*> FindPath(int source, int sink)
    {
        ClearFringe();
        ClosedList.clear();
        OpenList.clear();

        NodeRecord startNode{ source, nullptr, 0 };
        Fringe.push(startNode);
        OpenList.push_back(startNode);
        NodeRecord current = startNode;
        while (!Fringe.empty())
        {
            current = Fringe.top();
            Fringe.pop();

            auto it = std::find(OpenList.begin(), OpenList.end(), current);
            if (it != OpenList.end())
            {
                OpenList.erase(it);
            }

            if (current.Node == sink)
                break;
            std::vector<Edge*> connections = GraphSpace->GraphMap[current.Node].Connections;
            for (auto edge : connections)
            {
                int toNodeFromEdge = edge->sink;
                float toNodeCost = current.CostSoFar + edge->cost;
                NodeRecord endNode{ toNodeFromEdge,edge,toNodeCost };
                if (std::find(ClosedList.begin(), ClosedList.end(), endNode) != ClosedList.end())
                {
                    continue;
                }
                else if (std::find(OpenList.begin(), OpenList.end(), endNode) != OpenList.end())
                {
                    NodeRecord node = *std::find(OpenList.begin(), OpenList.end(), endNode);
                    if (node.CostSoFar <= endNode.CostSoFar)
                        continue;
                }
                UpdateOpenList(endNode);
            }
            ClosedList.push_back(current);
        }
        std::vector<Edge*> path;
        if (current.Node != sink)
        {
            return path;
        }
        while (current.Node != source)
        {
            path.push_back(current.IncomingEdge);
            NodeRecord node{ current.IncomingEdge->start,nullptr,0 };
            auto it = std::find(ClosedList.begin(), ClosedList.end(), node);
            current = *it;
        }
        std::reverse(path.begin(), path.end());
        return path;
    }

private:
    void UpdateOpenList(NodeRecord nodeToInsert)
    {
        std::vector<NodeRecord> recordsInFringe;
        while (Fringe.size() > 0)
        {
            NodeRecord record = Fringe.top();
            Fringe.pop();
            if (record == nodeToInsert)
                continue;
            recordsInFringe.push_back(record);
        }
        for (auto& record : recordsInFringe)
        {
            Fringe.push(record);
        }
        Fringe.push(nodeToInsert);
        auto it = std::find(OpenList.begin(), OpenList.end(), nodeToInsert);
        if (it != OpenList.end())
        {
            OpenList.erase(it);
        }
        OpenList.push_back(nodeToInsert);
    }
    void ClearFringe()
    {
        while (Fringe.size() > 0)
        {
            Fringe.pop();
        }
    }
    NetworkGraph* GraphSpace = nullptr;
    std::priority_queue<NodeRecord, std::vector<NodeRecord>, std::greater<NodeRecord>> Fringe;
    std::vector<NodeRecord> ClosedList;
    std::vector<NodeRecord> OpenList;
};
int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cout << "Error: requires csv path, initial supercharger names, and final supercharger names" << std::endl;
        return -1;
    }

    std::string charger_csv_path = argv[1];
    std::string initial_charger_name = argv[2];
    std::string goal_charger_name = argv[3];


    NetworkGraph NG(charger_csv_path.c_str());
    NG.SetupConnections();

    PathFinder P(&NG);
    auto Vec = P.FindPath(initial_charger_name, goal_charger_name);

    double CarRange = MAXRANGE;
    float TotalTimeCost = 0.f;
    for (auto CS : Vec)
    {
        double dist = GetDistanceFromCoordinates(NG.GraphMap[CS->start].Coord, NG.GraphMap[CS->sink].Coord);
        double TimeToCharge = 0.f;
        if (CarRange < dist)
        {
            double diffInRange = (dist - CarRange);
            CarRange += diffInRange;
            TimeToCharge = diffInRange / CHARGESPEED;
        }
        std::cout << NG.GraphMap[CS->start].Name<<",";
        if (TimeToCharge > 0.0) 
        {
            printf("%.3f,", RoundTwoDecimals(TimeToCharge));
        }
        CarRange -= dist;
        TotalTimeCost += (CS->cost + TimeToCharge);
    }
    std::cout << NG.GraphMap[Vec[Vec.size() - 1]->sink].Name << std::endl;
    std::cout << "TotalTime:" << TotalTimeCost << std::endl;
    return 0;
}