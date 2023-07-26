#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>
#include <iterator>
#include "read_data.cpp"
#include "motif_id.cpp"
#include <string>
#include <omp.h>
#include <atomic>

using namespace std;

inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

void print(std::unordered_set<int> const &s)
{
    std::copy(s.begin(),
            s.end(),
            std::ostream_iterator<int>(std::cout, ","));
}

int main(int argc, char *argv[])
{
	chrono::system_clock::time_point start;
	chrono::system_clock::time_point run_start;
	chrono::duration<double> dur;
	int progress;

	int num_threads = stoi(argv[1]);
	// cout << argv[2] << endl;
	// cout << argv[3] << endl;

	omp_set_num_threads(num_threads);

	// string graphFile = "/home/isis/HgNonLinear/data/Sociopatterns/hgs/aggr_1min_cliques_thr2_LH10.hg";
	// /home/isis/HgNonLinear/data/Sociopatterns/hgs/aggr_1min_cliques_thr2_LH10.hg
	// string dir = "/home/isis/HgNonLinear/data/Sociopatterns/Opinions_Data/";
	// string fname = "LH10_1.random";
	// /Users/ddevin/Desktop/MoCHy/reddit_so/reddit_trims/trimestre1/
	// /Users/ddevin/Documents/vscode/DevCommunities/so_data/2023-07-14_11-53-29/so_dueanni.hgf
	// trimestre1.txt
	// string dir = argv[2];
	// string fname = argv[3];
	string graphFile = "/Users/ddevin/Documents/vscode/DevCommunities/so_data/2023-07-14_11-53-29/so_dueanni.hgf";
	// string graphFile = "/Users/ddevin/Desktop/MoCHy/test_hgs.txt";
	// print graphFile
	cout << graphFile << endl;
	// Read data
	start = chrono::system_clock::now();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);
	cout << "OKOKOKOK" << endl;

	vector< vector<int> > motifs;

	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "Vsize" << V << endl;

	// cout << "# of nodes: " << V << '\n';
	// cout << "# of hyperedges: " << E << '\n';
	// dur = std::chrono::system_clock::now() - start;
	// cout << "Reading data done: "
	// 	<< dur.count() << " sec" << endl;
	// cout << "------------------------------------------" << endl << endl;


	// Make adjacency list
	start = chrono::system_clock::now();
	run_start = chrono::system_clock::now();
	hyperedge2node.resize(E); hyperedge2node_set.resize(E);
	vector< vector<pair<int, int> > > hyperedge_adj;
	vector< unordered_map<int, int> > hyperedge_inter;
	hyperedge_adj.resize(E);
	hyperedge_inter.resize(E);
	vector< vector<long long> > upd_time(num_threads);
	for (int i = 0; i < num_threads; i++)
		upd_time[i].resize(E, -1LL);
	
	cout << "build adj list start" << endl;


	#pragma omp parallel for
	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		int tid = omp_get_thread_num();
		int deg_a = 0;
		long long l_hyperedge_a = (long long)hyperedge_a;
		for (const int &node: hyperedge2node[hyperedge_a]){
			for (const int &hyperedge_b: node2hyperedge[node]){
				if (hyperedge_b == hyperedge_a) continue;
				if ((upd_time[tid][hyperedge_b] >> 31) ^ hyperedge_a){
					upd_time[tid][hyperedge_b] = (l_hyperedge_a << 31) + deg_a; deg_a++;
					hyperedge_adj[hyperedge_a].push_back({hyperedge_b, 0});
				}
				hyperedge_adj[hyperedge_a][(int)(upd_time[tid][hyperedge_b] & 0x7FFFFFFFLL)].second++;
			}
		}
	}

	#pragma omp parallel for 
	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		sort(hyperedge_adj[hyperedge_a].begin(), hyperedge_adj[hyperedge_a].end());
		int deg_a = hyperedge_adj[hyperedge_a].size();
		hyperedge_inter[hyperedge_a].rehash(deg_a);
		for (int i = 0; i < deg_a; i++){
			int hyperedge_b = hyperedge_adj[hyperedge_a][i].first;
			int C_ab = hyperedge_adj[hyperedge_a][i].second;
			hyperedge_inter[hyperedge_a].insert({hyperedge_b, C_ab});
		}
	}

	// dur = std::chrono::system_clock::now() - start;
	cout << "Adjacency list construction done: " << endl;
	// 	<< dur.count() << " sec" << endl;
	// cout << "------------------------------------------" << endl << endl;


	// Exact hypergraph motif counting
	start = chrono::system_clock::now();
	vector< vector<long long> > h_motif(num_threads);
	vector< vector< vector<int> > > data_motifs(num_threads);
	vector< vector<int> > intersection(num_threads);

	for (int i = 0; i < num_threads; i++){
		std::fill(upd_time[i].begin(), upd_time[i].end(), -1LL);
		intersection[i].resize(V);
		h_motif[i].resize(30, 0);
		
	}

	vector< std::atomic_flag > mutex(E);

	#pragma omp parallel for
	for(int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		int tid = omp_get_thread_num();

		long long l_hyperedge_a = (long long)hyperedge_a;
		int size_a = (int)hyperedge2node[hyperedge_a].size();
		int deg_a = (int)hyperedge_adj[hyperedge_a].size();

		for (int i = 0; i < deg_a; i++){
			int hyperedge_b = hyperedge_adj[hyperedge_a][i].first, C_ab = hyperedge_adj[hyperedge_a][i].second;
			int size_b = (int)hyperedge2node[hyperedge_b].size();
			int deg_b = (int)hyperedge_adj[hyperedge_b].size();

			const auto &nodes = hyperedge2node_set[hyperedge_b]; auto it_end = nodes.end(); int cnt = 0;
			for (const int &node: hyperedge2node[hyperedge_a]){ if(nodes.find(node) != it_end) intersection[tid][cnt++] = node;}

			for (int j = i + 1; j < deg_a; j++){
				int hyperedge_c = hyperedge_adj[hyperedge_a][j].first, C_ca = hyperedge_adj[hyperedge_a][j].second;
				int size_c = (int)hyperedge2node[hyperedge_c].size();
				int deg_c = (int)hyperedge_adj[hyperedge_c].size();

				int C_bc = 0;
				while (mutex[hyperedge_b].test_and_set(std::memory_order_acquire));
				C_bc = hyperedge_inter[hyperedge_b][hyperedge_c];
				mutex[hyperedge_b].clear(std::memory_order_release);
				//auto it = hyperedge_inter[hyperedge_b].find(hyperedge_c);
				//if (it == hyperedge_inter[hyperedge_b].end()) C_bc = 0;
				//else C_bc = (*it).second;
				if (C_bc){
					if (hyperedge_a < hyperedge_b){
						int g_abc = 0;
						const auto &nodes = hyperedge2node_set[hyperedge_c]; auto it_end = nodes.end();

						//cout << "set " << hyperedge2node_set[hyperedge_a] << endl;
						
						for (int k = 0; k < C_ab; k++){ if(nodes.find(intersection[tid][k]) != it_end) g_abc++;}
						int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);

						data_motifs[tid].push_back({motif_index,hyperedge_a,hyperedge_b,hyperedge_c});
						h_motif[tid][motif_index]++;
					}
				} else {
					int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, 0, C_ca, 0);
						data_motifs[tid].push_back({motif_index,hyperedge_a,hyperedge_b,hyperedge_c});

					h_motif[tid][motif_index]++;
				}
			}
			
		}	
	}

	cout << "KOKOKOKOKOKO" << endl;


	int index = 0;	
	vector<long long> h_motif_final(30, 0);
	for (int i = 0; i < 30; i++){	
		for (int j = 0; j < num_threads; j++) h_motif_final[i] += h_motif[j][i];
		if (i == 0 || i == 1 || i == 4 || i == 6) continue;
		cout << fixed << "motif " << ++index << ": " << fixed << h_motif_final[i] << endl;
	}

	dur = std::chrono::system_clock::now() - run_start;
	double runtime = (double)dur.count();

	dur = std::chrono::system_clock::now() - start;
	cout << "\nMotif counting: "
		<< dur.count() << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;
	cout << "------------------------------------------" << endl << endl;

	fstream fout;
	
	int indici[30] = {-1, -1, 1, 2, -1, 3, -1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	// print length of
	string output = "/Users/ddevin/Documents/vscode/DevCommunities/so_data/2023-07-14_11-53-29/so_dueanni.motifs400random";
	cout << output << endl;
	fout.open(output, ios::app);
	int limit = 400;
	// set an array of 30 elements to 0
	int h_motif_size[30] = {0};
	// set a random seed
	srand(time(NULL));
	for(int i = 0; i < num_threads; i++){
		for(int j = 0; j < data_motifs[i].size(); j++){

			if (data_motifs[i][j][0] == 0 || data_motifs[i][j][0] == 1 || data_motifs[i][j][0] == 4 || data_motifs[i][j][0] == 6) continue;
			
			if (h_motif_size[data_motifs[i][j][0]] > limit) continue;
			// if that motif has been found more than 3000 times, pick at random
			if (h_motif_final[data_motifs[i][j][0]] > 3000) {
				if (rand() % 2 == 0) 
					continue;
			}

			h_motif_size[data_motifs[i][j][0]]++;

			auto a = hyperedge2node_set[data_motifs[i][j][1]];
			auto b = hyperedge2node_set[data_motifs[i][j][2]];
			auto c = hyperedge2node_set[data_motifs[i][j][3]];
			auto index = indici[data_motifs[i][j][0]];

			fout << index << ",";
			fout << data_motifs[i][j][1] << ",";
			fout << data_motifs[i][j][2] << ",";
			fout << data_motifs[i][j][3];
			// for (auto it = a.begin(); it != a.end(); it++)
			// 	// build a string with *it
			// 	fout << to_string(*it) << " ";
			// 	// sa.append(to_string(*it) + " ");
			// fout << ",";

			// for (auto it = b.begin(); it != b.end(); it++) 
			// 	// build a string with *it
			// 	fout << to_string(*it) << " ";
			// fout << ",";

			// for (auto it = c.begin(); it != c.end(); it++) 
			// 	// build a string with *it
			// 	fout << to_string(*it) << " ";

			fout << endl;
		}
	}
	fout.close();

	node2hyperedge.clear();
	hyperedge2node.clear();
	hyperedge_adj.clear();
	hyperedge_inter.clear();
	data_motifs.clear();
	return 0;
}
