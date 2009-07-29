
#include <exception>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <list>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>

using namespace std;

class Graph {
protected:
    int * cumul_nodes;
    int cumul_max_size;
public:
    virtual bool is_connected(int i, int j) = 0;
    virtual void connect(int i, int j) = 0;
    virtual int addnode() = 0;
    virtual int degree(int i) = 0;
    virtual int nodenum() = 0;
    virtual int edgenum() = 0;
    virtual int totaldegree() = 0;

    void init_cumulnodes(int n, int m){
        cumul_max_size = ((m * (m-1)) / 2 + m * (n - m)) * 2;
        cumul_nodes = new int[cumul_max_size];
    }
    void free_cumulnodes(){
        delete[](cumul_nodes);
    }

    void choose_for_small(int m, int * ret, int * degrees/* cache */, int totaldegree){
        int nn = nodenum();
        memset(ret, -1, sizeof(int) * m);
        
        if (m > nn){
            throw logic_error("error : choose : m > nn");
            for(int i = 0;i < nn;i++){
                ret[i] = i;
            }
            return;
        }
        
        for(int t = 0;t < m;t++){
            int degreestock = rand() % totaldegree;
            int chosen_node = -1;
            
            for(int i = 0;i < nn;i++){
                bool not_chosen = true;
                for(int j = 0;j < t;j++){
                    if (i == ret[j]){ not_chosen = false; }
                }
                if (not_chosen) {
                    degreestock -= degrees[i];
                    if (degreestock < 0){
                        chosen_node = i;
                        break;
                    }
                }
            }

            if (chosen_node < 0){
                throw logic_error("error : choose");
            }

            ret[t] = chosen_node;
            totaldegree -= degree(chosen_node);
        }
    }

    void choose(int m, int * ret, int * degrees/* cache */, int totaldegree){
        int nn = nodenum();
        memset(ret, -1, sizeof(int) * m);

        if (nn < 1000){
            choose_for_small(m, ret, degrees, totaldegree);
        } else {
            for(int t = 0;t < m;t++){
                int nodeid;
                int nodeidx;
                bool chosen;
                
                do {
                    nodeidx = rand() % totaldegree;
                    nodeid = cumul_nodes[nodeidx];
                    chosen = false;
                    for(int j = 0;j < t;j++){
                        if (nodeid == ret[j]){
                            chosen = true;
                            // break;
                        }
                    }
                    if (!chosen){
                        break;
                    }
                }while(true);
                
                ret[t] = nodeid;
            }
        }
    }
};

class AdjMatrixGraph : public Graph {
private:
    int ** adj_matrix;
    int totalnodesize;
    int nodenum_;
    int edgenum_;
    int * degree_;
    int totaldegree_;

    bool valid_range(int i){
        return 0 <= i && i < nodenum_;
    }

    void setadj(int i, int j, int val){
        if (i == j) { throw invalid_argument("setadj"); }
        if (!valid_range(i) || !valid_range(j)){
            throw out_of_range("setadj");
        }
        if(i > j){
            adj_matrix[i][j] = val;
        } else {
            adj_matrix[j][i] = val;
        }
    }

public:
    bool is_connected(int i, int j){
        if (!valid_range(i) || !valid_range(j)){
            throw out_of_range("is_connected");
        }
        if (i == j) return false;
        return (i > j ? adj_matrix[i][j] : adj_matrix[j][i]) > 0; 
    }

    void connect(int i, int j){
        if (!is_connected(i, j)){
            setadj(i, j, 1);
            degree_[i]++;
            degree_[j]++;
            cumul_nodes[totaldegree_] = i;
            cumul_nodes[totaldegree_ + 1] = j;
            totaldegree_ += 2;
            edgenum_++;
        }
    }

    int addnode(){
        if (nodenum_ >= totalnodesize){
            throw overflow_error("addnode");
        }
        return nodenum_++;
    }

    AdjMatrixGraph(int nodeN): totalnodesize(nodeN){
        // parameter check
        if (nodeN < 0){
            throw invalid_argument("Invalid parameters for BA model");
        }

        // initialize variables
        nodenum_ = 0;
        edgenum_ = 0;
        totaldegree_ = 0;
        adj_matrix = new int*[nodeN];
        degree_ = new int[nodeN];
        for(int i = 1;i < nodeN;i++){
            adj_matrix[i] = new int[i];
            degree_[i] = 0;
            for(int j = 0;j < i;j++){
                adj_matrix[i][j] = 0;
            }
        }
    }

    ~AdjMatrixGraph(){
        for(int i = 0;i < totalnodesize;i++){
            delete[](adj_matrix[i]);
        }
        delete[](adj_matrix);
        delete[](degree_);
        free_cumulnodes();
    }

    int degree(int i){
        if(!valid_range(i)){
            throw out_of_range("degree");
        }
        return degree_[i];
    };
    
    int nodenum(){ return nodenum_;};
    int edgenum(){ return edgenum_;};
    int totaldegree(){
        return totaldegree_;
    }
};

class AdjListGraph : public Graph {
    set<int> ** adjlists;
    int * degree_;
    int totalnodesize;
    int nodenum_;
    int edgenum_;
    int totaldegree_;

    bool valid_range(int i){
        return 0 <= i && i < nodenum_;
    }

public:
    AdjListGraph(int nodeN): totalnodesize(nodeN){
        nodenum_ = 0;
        edgenum_ = 0;
        totaldegree_ = 0;
        degree_ = new int[nodeN];
        adjlists = new set<int>*[nodeN];
    }

    ~AdjListGraph(){
        free_cumulnodes();
        delete[](degree_);
        delete[](adjlists);
    }
    
    bool is_connected(int i, int j){
        if (!valid_range(i) || !valid_range(j)){
            throw out_of_range("is_connected");
        }
        return adjlists[i]->find(j) != adjlists[i]->end();
    }
    
    void connect(int i, int j){
        if (!is_connected(i, j)){
            adjlists[i]->insert(j);
            adjlists[j]->insert(i);
            edgenum_++;
            degree_[j]++;
            degree_[i]++;
            cumul_nodes[totaldegree_] = i;
            cumul_nodes[totaldegree_ + 1] = j;
            totaldegree_ += 2;
        }
    }
    
    int addnode(){
        if (nodenum_ >= totalnodesize) { throw overflow_error(""); }
        adjlists[nodenum_] = new set<int>();
        return nodenum_++;
    }
    
    inline int degree(int i){
        if (!valid_range(i)){ throw out_of_range("degree"); }
        return degree_[i];
    }
    
    inline int nodenum(){ return nodenum_; }
    inline int edgenum(){ return edgenum_; }
    inline int totaldegree() { return totaldegree_; }
};

class BAGraph {
private:
    // parameters
    int n; // total number of nodes
    int m;

    // AdjMatrixGraph g;
    AdjListGraph g;

public:

    // primitive interfaces
    
    inline int getnodenum(){ return g.nodenum(); }
    inline int getedgenum(){ return g.edgenum(); }
    inline int getdegree(int node){ return g.degree(node); }

    BAGraph(int nodeN, int edgeM): n(nodeN), m(edgeM), g(nodeN){
        g.init_cumulnodes(nodeN, edgeM);
        
        // make m-node complete graph
        for(int i = 0;i < m;i++){
            g.addnode();
            for(int j = 0;j < i;j++){
                g.connect(i, j);
            }
        }

        // add a node iteratively
        int chosen[m];
        int *degrees = new int[n];
        for(int i = 0;i < g.nodenum();i++){ degrees[i] = g.degree(i); }
        while(g.nodenum() < n){
            g.choose(m, chosen, degrees, g.totaldegree());
            
            int newnode = g.addnode();
            for(int i = 0;i < m;i++){
                g.connect(chosen[i], newnode);
                degrees[chosen[i]]++;
            }
            degrees[newnode] = m;
        }
        delete[](degrees);

        if (g.nodenum() != nodeN) {
            cerr << "edge: " << g.nodenum() << endl;
            throw logic_error("node num error");
        }
        if (g.edgenum() != (edgeM * (edgeM - 1)) / 2 + (nodeN - edgeM) * edgeM) {
            cerr << "edge: " << g.edgenum() << endl;
            throw logic_error("edge num error");
        }
    }


    double degreecorrelation(int degree){
        int nodenum = 0;
        int sum = 0;
        for(int i = 0;i < getnodenum();i++){
            if (degree == getdegree(i)){
                for(int j = 0;j < getnodenum();j++){
                    if (i == j) continue;
                    if(g.is_connected(i, j)){
                        nodenum++;
                        sum += getdegree(j);
                    }
                }
            }
        }

        return (nodenum == 0 ? 0 : (double)sum / nodenum);
    }

    double degreermsdev(int degree){
        int nodenum = 0;
        int sum = 0;
        double v = 0;
        for(int i = 0;i < getnodenum();i++){
            if (degree == getdegree(i)){
                for(int j = 0;j < getnodenum();j++){
                    if (i == j) continue;
                    if(g.is_connected(i, j)){
                        nodenum++;
                        int deg = getdegree(j);
                        sum += deg;
                        v += deg * deg;
                    }
                }
            }
        }

        if (nodenum == 0){
            return 0;
        }
        double avg = (double)sum / nodenum;
        v -= avg * avg;
        v *= nodenum;
        v /= nodenum - 1;
        v = sqrt(v);

        return v;
    }
};

void print_degree_distribution(BAGraph & graph){
    int count[graph.getnodenum()];
    for(int i = 0;i < graph.getnodenum();i++){
        count[i] = 0;
    }
    for(int i = 0;i < graph.getnodenum();i++){
        count[graph.getdegree(i)]++;
    }
    for(int i = 0;i < graph.getnodenum();i++){
        double p = count[i];
        p /= graph.getnodenum();
        cout << i << "\t" << p << endl;
    }
}

void print_degree_correlation(BAGraph & graph){
    for(int i = 0; i < graph.getnodenum();i++){
        cout << i << "\t" << graph.degreecorrelation(i) << endl;
    }
}

void print_degree_rmsdev(BAGraph & graph){
    for(int i = 0; i < graph.getnodenum();i++){
        cout << i << "\t" << graph.degreermsdev(i) << endl;
    }
}

void show_help(char * program){
    cerr << "** Barabasi-Albert model **" << endl << endl;
    cerr << "Usage: " << program << " nodenum initsize command" << endl << endl;
    cerr << "nodenum : the number of nodes" << endl;
    cerr << "initsize: initial size of the (complete) graph in BA model" << endl;
    cerr << "command : dist ... Degree distribution" << endl
         << "          corr ... Degree correlation" << endl;
    exit(1);
}

typedef enum {
    CMD_DIST, CMD_CORR,
} command_t;

int main(int argc, char ** argv){
    // srand(0);
    srand(time(NULL));

    if (argc < 4){
        show_help(argv[0]);
    }

    int n = strtol(argv[1], NULL, 10);
    int m = strtol(argv[2], NULL, 10);
    string cmdstr(argv[3]);

    command_t command = CMD_DIST; // pseudo initialization
    if(cmdstr.compare("dist") == 0){
        command = CMD_DIST;
    } else if (cmdstr.compare("corr") == 0){
        command = CMD_CORR;
    } else {
        show_help(argv[0]);
    }

    cerr << "Building graph(n=" << n << ") ... ";
    cerr.flush();
    BAGraph graph(n, m);
    cerr << "done" << endl;
    cerr << "nodenum=" << graph.getnodenum()
         << ", edgenum=" << graph.getedgenum() << endl;

    switch(command){
    case CMD_DIST:
        print_degree_distribution(graph);
        break;
    case CMD_CORR:
        print_degree_correlation(graph);
        break;
    }
    
    // print_degree_rmsdev(graph);

    return 0;
}
