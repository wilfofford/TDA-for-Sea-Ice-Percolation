#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>

class Edge
{
    public:
        int start;
        int end;
        float weight;

        Edge(int start, int end, float weight) 
            : start(start), end(end), weight(weight) {}
};

class Graph
{
    public:
        std::vector<Edge> edges;
        std::vector< std::pair<int,int> > vertices;
        std::vector<float> births;

        Graph(std::vector< std::pair<int,int> > vertices, 
              std::vector<Edge> edges, std::vector<float> births)
            : vertices(vertices), edges(edges), births(births) {}
};

class Union_Find
{
    std::vector<int> parent;
    std::vector<int> rank;
    std::vector<float> birth;
    std::vector<int> birthplace;

    public:
        Union_Find(int n): parent(n), rank(n,1), birth(n,0), birthplace(n)
        {
            for (int i = 0; i<n; ++i)
            {
                parent[i]=i;
                birthplace[i]=i;
            }
        }

        void set_birth(int i, float val) {birth[i]=val;}

        float get_birth(int i) {return birth[i];}

        void set_birthplace(int i, int val) {birthplace[i]=val;}

        int get_birthplace(int i) {return birthplace[i];}

        int find(int x)
        {
            int y = x, z = parent[y];
            while (z != y) {
                y = z;
                z = parent[y];
            }
            y = parent[x];
            while (z != y) {
                parent[x] = z;
                x = y;
                y = parent[x];
            }
            return z;
        }

        void link(int x, int y)
        {
            x = find(x);
            y = find(y);
            if (x == y)
                return;
            if (rank[x] > rank[y]) {
                parent[y] = x;
                if (birth[x]>birth[y]) {
                    birth[x]=birth[y];
                    birthplace[x]=birthplace[y];
                }
                
            } else {
                parent[x] = y;
                if (birth[y]>birth[x]) {
                    birth[y]=birth[x];
                    birthplace[y]=birthplace[x];
                }
                if (rank[x] == rank[y])
                    ++rank[y];
            }
        }
};

class Barcode
{
    public:
    std::vector< std::pair<float,float> > dgm;
    std::vector< std::pair< std::pair<int,int>, std::pair<int,int> > > locations;

    Barcode(std::vector< std::pair<float,float> > dgm, std::vector< std::pair< std::pair<int,int>, std::pair<int,int> > > locations)
        : dgm(dgm), locations(locations) {}
};

struct compare_edges
{
    bool operator()(const Edge& e1, const Edge& e2) const
    {
        return (e1.weight > e2.weight);
    }    
};

Barcode compute_dim0_pairs(Graph graph)
{
    std::vector<std::pair<float,float>> dgm;
    std::vector< std::pair< std::pair<int,int>, std::pair<int,int> > > locations;
    int n = graph.vertices.size();
    Union_Find dset(n);
    for (int i=0; i<n; i++)
    {
        dset.set_birth(i, graph.births[i]);
    }
    std::vector<Edge> edges = graph.edges;
    std::sort(edges.rbegin(),edges.rend(), compare_edges());
    for (auto e : edges)
    {
        int u = dset.find(e.start), v = dset.find(e.end);
        if (u!=v)
        {
            float birth;
            int birthplace;
            int deathplace;
            float b_u=dset.get_birth(u), b_v=dset.get_birth(v);
            if (b_u>b_v)
            {
                birth = b_u;
                birthplace = dset.get_birthplace(u);
                deathplace = e.end;
            } else {
                birth = b_v;
                birthplace = dset.get_birthplace(v);
                deathplace=e.start;
            }
            float death = e.weight;
            if (death>birth)
            {
                dgm.push_back(std::make_pair(birth,death));
                locations.push_back(std::make_pair(graph.vertices[birthplace],graph.vertices[deathplace]));
            }
            dset.link(u,v);
        }
        
    }
    for (int i=0; i<n; ++i)
        if (dset.find(i)==i)
        {
            dgm.push_back(std::make_pair(dset.get_birth(i),std::numeric_limits<float>::infinity()));
            locations.push_back(std::make_pair(graph.vertices[dset.get_birthplace(i)],std::make_pair(-1,-1)));
        }
    Barcode bcode(dgm,locations);
    return bcode;
}


Graph lower_star_graph(std::vector<float> img, int m, int n)
{

    std::vector<std::vector<int>> locations(m,std::vector<int>(n));
    std::vector<Edge> edges;
    std::vector<std::pair<int,int>> flatcoords(m*n);
    std::vector<float> flatimg(m*n);


    for(int i=0; i<m; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            locations[i][j]=n*i+j;
            flatcoords[n*i+j]=std::make_pair(i,j);            
        }
    }

    int added[4][2] = { {-1,-1},{-1,0},{-1,1},{0,1} };


    for(int i=0; i<m; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            for (int t=0; t<4; ++t)
            {
                if (i+added[t][0]>=0 && i+added[t][0]<m
                    &&j+added[t][1]>=0 && j+added[t][1]<n)
                    {
                        edges.push_back(Edge(  locations[i][j],  locations[i+added[t][0]][j+added[t][1]],  std::max(img[n*i+j],img[n*(i+added[t][0])+j+added[t][1]]) ));
                    }
            }
        }
    }
    return Graph(flatcoords,edges,img);
}



int main()
{   
    std::ifstream in;
    in.open("storing_img.txt");
    
    std::string entry;
    int m,n;
    
    in >> entry;
    m = std::stoi(entry);
    in >> entry;
    n = std::stoi(entry);
    std::vector< float > vect(m*n);
    int i=0;
    while (in >> entry)
    {
        vect[i]=std::stof(entry);
        i++;
    }

    in.close();

    Graph graph = lower_star_graph(vect,m,n);

    

    Barcode barcode = compute_dim0_pairs(graph);

   

    std::ofstream temp;   
    temp.open("storing_barcode.txt", std::ofstream::out | std::ofstream::trunc);
    
    int l=barcode.dgm.size();
    int j=barcode.locations.size();
    
    for (int i=0;i<l;++i)
    {
        temp<<barcode.dgm[i].first<<" "<<barcode.dgm[i].second<<" "<<barcode.locations[i].first.first<<" "<<barcode.locations[i].first.second<<" "
        <<barcode.locations[i].second.first<<" "<<barcode.locations[i].second.second<<"\n";
    }
 
    
    temp.close();


    

}
