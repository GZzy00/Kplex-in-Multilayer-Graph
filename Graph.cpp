#include "Graph.h"
#include "KPlex_BB_matrix.h"

using namespace std;

Graph::Graph(const char *_dir, const int _K, const int _Sigma) /*定义图的构造函数，接收目录字符串和整数_K作为参数。*/
{
	dir = string(_dir); /*将传入的目录字符串转换为string类型并赋值给成员变量dir。*/
	K = _K;				/*将传入的_K赋值给成员变量K。*/
	Sigma = _Sigma;
	n = m = l = 0;			
	pstart = nullptr;	/*指向邻接表相关数组的指针初始化为nullptr。*/
	pend = pend_buf = nullptr;
	edges = nullptr;
	lstart = nullptr;
	layers = nullptr;
	kplex.clear();					 /*清空kplex集合，为后面的优化做准备。*/
	s_degree = s_edges = NULL;		 /*初始化用于排序的数组指针为NULL*/
	s_pstart = s_pend = NULL;		 /*初始化用于排序的起始和结束位置指针为NULL*/
	s_peel_sequence = s_core = NULL; /* 初始化剥离序列和核心度数数组指针为NULL*/
	s_vis = NULL;					 /*初始化访问标记数组指针为NULL*/
	s_heap = NULL;					 /*初始化线性堆指针为NULL*/
	s_edgelist_pointer = NULL;		 /*初始化边列表指针数组为NULL*/
	s_tri_cnt = s_edge_list = NULL;	 /*初始化三角形计数和边列表数组指针为NULL*/
	s_active_edgelist = NULL;		 /*初始化活跃边列表指针为NULL*/
	s_deleted = NULL;				 /*初始化已删除标记数组指针为NULL*/
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != NULL)
	{
		delete[] pend_buf;
		pend_buf = NULL;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (layers != nullptr)
	{
		delete[] layers;
		layers = nullptr;
	}
	if (lstart != nullptr)
	{
		delete[] lstart;
		lstart = nullptr;
	}
	if (s_degree != NULL)
	{
		delete[] s_degree;
		s_degree = NULL;
	}
	if (s_pstart != NULL)
	{
		delete[] s_pstart;
		s_pstart = NULL;
	}
	if (s_pend != NULL)
	{
		delete[] s_pend;
		s_pend = NULL;
	}
	if (s_edges != NULL)
	{
		delete[] s_edges;
		s_edges = NULL;
	}
	if (s_peel_sequence != NULL)
	{
		delete[] s_peel_sequence;
		s_peel_sequence = NULL;
	}
	if (s_core != NULL)
	{
		delete[] s_core;
		s_core = NULL;
	}
	if (s_vis != NULL)
	{
		delete[] s_vis;
		s_vis = NULL;
	}
	if (s_heap != NULL)
	{
		delete s_heap;
		s_heap = NULL;
	}
	if (s_edgelist_pointer != NULL)
	{
		delete[] s_edgelist_pointer;
		s_edgelist_pointer = NULL;
	}
	if (s_active_edgelist != NULL)
	{
		delete[] s_active_edgelist;
		s_active_edgelist = NULL;
	}
	if (s_deleted != NULL)
	{
		delete[] s_deleted;
		s_deleted = NULL;
	}
}

void Graph::read_graph() 
{
	printf("# Reading graph from file\n");
	FILE *f = Utility::open_file((dir).c_str(), "r");

	fscanf(f, "%u%u%lu", &l, &n, &m);
	m *= 2;
	printf("*\tl = %s; n = %s; m = %s (undirected)\n", Utility::integer_to_string(l).c_str(), Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	vector<tuple<ui,ui,ui>> vtp;
	for(ui i = 0; i < m/2; i ++) 
	{
		ui a, b, c;
		fscanf(f, "%u%u%u", &a, &b, &c);
		a--;
		b--;
		c--;
		if(a >= l || b >= n || c >= n) {
			printf("!!! Vertex IDs must be between 0 and n-1 and Layer IDS must be between 0 and l. Exit !!!\n");
			return ;
		}
		vtp.pb(mtp(a, b, c));
		vtp.pb(mtp(a, c, b));
	}
	sort(vtp.begin(), vtp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ept[l*n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];
	if(layers != nullptr) delete[] layers;
	layers = new ui[m];
	if(lstart != nullptr) delete[] lstart;
	lstart = new ept[l+1];

	pstart[0] = 0;
	lstart[0] = 0;
	ui idx = 0;
	for (ui i = 0; i < l; i ++)
	{
		lstart[i + 1] = lstart[i];
		for (ui j = 0; j < n; j ++)
		{
			pstart[n*i + j + 1] = pstart[n*i + j];
			while (idx < vtp.size() && get<0>(vtp[idx]) == i && get<1>(vtp[idx]) == j)
			{
				layers[pstart[n*i + j + 1]] = get<0>(vtp[idx]);
				edges[pstart[n*i + j + 1] ++] = get<2>(vtp[idx ++]);
				lstart[i + 1] ++;
			}
		}
	}
	// printf("%lu  %lu  %lu  %lu\n", pstart[0*n + 1], pstart[0*n + 2], lstart[2], lstart[1]);
	// for(ui i = 0;i < n;i ++) 
	// {
	// 	pstart[i+1] = pstart[i];
	// 	while(idx < vtp.size() && get<0>(vtp[idx]) == i)
	// 	{
	// 		layers[pstart[i+1]] = get<1>(vtp[idx]);
	// 		edges[pstart[i+1] ++] = get<2>(vtp[idx ++]);
	// 	}
	// }

	fclose(f);

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::kPlex_exact() {
	Timer t;

	assert(K > 0&&K < n);

	kplex.clear();
	
	// TODO
	// **启发式搜索一个较大的解**
	heuristic_kplex_max_degree(10);

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *multilayer_core = new ui[n*l];
	ui *degree = new ui[n*l];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);
	
	
	// TODO
	// **寻找上界，并更新启发式解**
	ui UB = degen(n, peel_sequence, core, multilayer_core, pstart, edges, degree, vis, heap, true);

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] multilayer_core;

	assert(kplex.size() >= K);

	if(kplex.size() < UB) {
		ui old_size = kplex.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];

		// **对图进行约减**
		shrink_graph(n, m, peel_sequence, core, out_mapping, NULL, rid, pstart, edges, true);
		delete[] core; core = NULL;

		ui *degree = new ui[n*l];
		for(ui i = 0;i < n*l;i ++) degree[i] = pstart[i+1] - pstart[i];

		ui *total_degree = new ui[n];
		memset(total_degree, 0, sizeof(ui)*n);
		for (ui i = 0; i < l; i ++)
		{
			for (ui j = 0; j < n; j ++)
			{
				total_degree[j] += degree[n*i + j];
			}
		}
		ListLinearHeap *linear_heap = new ListLinearHeap(n, (n-1) * l);
		linear_heap->init(n, (n-1)*l, peel_sequence, total_degree);
		assert(pend == nullptr);
		pend = new ept[n*l];

		ui *edgelist_pointer = new ui[m];
		oriented_triangle_counting(n, m, peel_sequence, pstart, pend, edges, edgelist_pointer, rid); // edgelist_pointer currently stores triangle_counts
		
		delete[] peel_sequence; peel_sequence = NULL;

		pend_buf = new ept[n*l];
		ui *edge_list = new ui[m];
		ui *tri_cnt = new ui[m/2];
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer, rid);

		for(ui i = 0;i < n*l;i ++) pend[i] = pstart[i+1];
		ui *active_edgelist = new ui[m>>1];
		ui active_edgelist_n = m>>1;
		for(ui i = 0;i < (m>>1);i ++) active_edgelist[i] = i;

		ui *Qe = new ui[m>>1];
		char *deleted = new char[m>>1];
		memset(deleted, 0, sizeof(char)*(m>>1));
		char *deleted_vertex = new char[n];
		memset(deleted_vertex, 0, sizeof(char)*n);
		char *exists = new char[n];
		memset(exists, 0, sizeof(char)*n);

		ui *Qv = new ui[n];
		ui Qv_n = 0;
		if(kplex.size()+1 > 2*K) {
			m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, true, kplex.size()+1-2*K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
			printf("*** After core-truss co-pruning: n = %s, m = %s\n", Utility::integer_to_string(n-Qv_n).c_str(), Utility::integer_to_string(m/2).c_str());
		}

		Timer tt;

		int max_n = n - Qv_n;
		s_degree = new ui[max_n*l];
		s_pstart = new ept[max_n*l+1];
		s_pend = new ept[max_n*l];
		s_edges = new ui[m];
		s_peel_sequence = new ui[max_n];
		s_core = new ui[max_n];
		s_vis = new char[max_n];
		s_heap = new ListLinearHeap(max_n,max_n-1);
		s_edgelist_pointer = new ui[m];
		s_tri_cnt = new ui[m/2];
		s_edge_list = new ui[m];
		s_active_edgelist = new ui[m/2];
		s_deleted = new char[m/2];

		KPLEX_BB_MATRIX *kplex_solver = new KPLEX_BB_MATRIX();
		kplex_solver->allocateMemory(max_n, m/2, l);

		vector<tuple<int,int,int> > vtp; vtp.reserve(m/2);
		ui *t_degree = new ui[n];

// 		ui max_n_prune = 0, max_n_search = 0, prune_cnt = 0, search_cnt = 0;
// 		double min_density_prune = 1, min_density_search = 1, total_density_prune = 0, total_density_search = 0;

		for(int i = 0;i < n&&m&&kplex.size() < UB;i ++) { 
			ui u, key;
			ui degree_u[l];
			ui total_degree_u = 0;
			bool ret_tmp = linear_heap->pop_min(u, key);
			assert(ret_tmp);
			for (ui t = 0; t < l; t ++) {
				degree_u[t] = degree[n*t + u];
				total_degree_u += degree[n*t + u];
			}
			for (ui t = 0; t < l - 1; t ++) {
				for (ui k = 0; k < l - t - 1; k ++) {
					if (degree_u[k] < degree_u[k+1]) exchange_value(degree_u[k], degree_u[k+1]);
				}
			}
			if(degree_u[Sigma-1] < kplex.size()+1-K) {
				if(total_degree_u != 0) { // degree[u] == 0 means u is deleted. it could be the case that degree[u] == 0, but key[u] > 0, as key[u] is not fully updated in linear_heap
					Qv[0] = u; Qv_n = 1;
					if(kplex.size()+1>2*K) m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, false, kplex.size()+1-2*K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
					else m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, false, 0, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
				}
				continue;
			}
			if(m == 0) break;
#ifndef NDEBUG
			if(total_degree_u != key) printf("u=%u, total_degree_u=%u, key=%u, kplex.size=%lu\n", u, degree[u], key, kplex.size());
#endif
			assert(total_degree_u == key);

			ui *ids = Qv;
			ui ids_n = 0;
#ifndef NDEBUG
			for(ui i = 0;i < n;i ++) assert(!exists[i]);
#endif
			if(kplex.size()+1 >= 2*K) {
				extract_subgraph_and_prune(u, ids, ids_n, rid, vtp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
				// if(ids_n) {
				// 	double density = (double(vp.size()*2))/ids_n/(ids_n-1);
				// 	total_density_prune += density; ++ prune_cnt;
				// 	if(density < min_density_prune) min_density_prune = density;
				// 	if(ids_n > max_n_prune) max_n_prune = ids_n;
				// }
			}
			else {
				extract_subgraph(u, ids, ids_n, rid, vtp, exists, pstart, pend, edges, deleted, edgelist_pointer);
				// double density = (double(vp.size()*2))/ids_n/(ids_n-1);
				// total_density_prune += density; ++ prune_cnt;
				// if(density < min_density_prune) min_density_prune = density;
				// if(ids_n > max_n_prune) max_n_prune = ids_n;
			}
			ui pre_size = kplex.size();
			if(ids_n > kplex.size()) {
				// double density = (double(vp.size()*2))/ids_n/(ids_n-1);
				// total_density_search += density; ++ search_cnt;
				// if(density < min_density_search) min_density_search = density;
				// if(ids_n > max_n_search) max_n_search = ids_n;
				kplex_solver->load_graph(ids_n, l, vtp);
				kplex_solver->kPlex(K, Sigma, kplex, true);
			}
			Qv[0] = u; Qv_n = 1;
			if(kplex.size() != pre_size&&kplex.size()+1 > 2*K) {
				for(ui j = 0;j < kplex.size();j ++) kplex[j] = ids[kplex[j]];
				//output_one_kplex(); break;
				m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, true, kplex.size()+1-2*K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
			}
			else if(kplex.size()+1>2*K) m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, false, kplex.size()+1-2*K, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
			else m -= 2*peeling(n, linear_heap, Qv, Qv_n, kplex.size()+1-K, Qe, false, 0, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deleted_vertex, degree, pstart, pend, edges, exists);
#ifndef NDEBUG
			printf("Number of remaining undirected edges: %s\n", Utility::integer_to_string(m/2).c_str());
#endif
		}

		// if(prune_cnt == 0) ++ prune_cnt;
		// if(search_cnt == 0) ++ search_cnt;
		// printf("prune_cnt: %u, max_n: %u, min_density: %.4lf, avg_density: %.4lf\n", prune_cnt, max_n_prune, min_density_prune, total_density_prune/prune_cnt);
		// printf("search_cnt: %u, max_n: %u, min_density: %.4lf, avg_density: %.4lf\n", search_cnt, max_n_search, min_density_search, total_density_search/search_cnt);

		if(kplex.size() > old_size) {
			for(ui i = 0;i < kplex.size();i ++) {
				assert(kplex[i] < n);
				kplex[i] = out_mapping[kplex[i]];
			}
		}
		printf("Solution size: %lu\n", kplex.size());
		printf("Solution: ");
		for(ui i = 0;i < kplex.size();i ++) {
			printf("%u ", kplex[i] + 1);
		}
		printf("\n");

		delete kplex_solver;
		delete linear_heap;
		delete[] t_degree;
		delete[] exists;
		delete[] out_mapping;
		delete[] rid;
		delete[] degree;
		delete[] total_degree;
		delete[] edgelist_pointer;
		delete[] tri_cnt;
		delete[] active_edgelist;
		delete[] Qe;
		delete[] Qv;
		delete[] deleted;
		delete[] deleted_vertex;
		printf("*** Search time: %s (microseconds)\n", Utility::integer_to_string(tt.elapsed()).c_str());
	}
	else {
		delete[] core;
		delete[] peel_sequence;
	}

// 	printf("\tMaximum kPlex Size: %lu, Total Time: %s (microseconds)\n", kplex.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// 启发式算法求一个较大的解
void Graph::heuristic_kplex_max_degree(ui processed_threshold)
{
	for(ui i = 0; i < 2*K; i++)
	{
		kplex.pb(i);
	}
}

// 寻找上界
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ui *multilayer_core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output) 
{
	Timer t;

	ui threshold = (kplex.size()+1 > K? kplex.size()+1-K: 0);
	// printf("%u\n", threshold);
	
	for(ui i = 0;i < l*n;i ++) 
	{
		degree[i] = pstart[i+1] - pstart[i]; 
	}


	ui queue_n = 0, new_size = 0;
	memset(vis, 0, sizeof(char)*n);

	// 将第Sigma大的度数都没有阈值大的节点加入剥离序列
	for (ui i = 0; i < n; i ++)
	{
		ui cnt = 0;
		for (ui j = 0; j < l; j ++)
		{
			if (degree[n*j + i] >= threshold) cnt++;
		}
		if (cnt < Sigma)
		{
			peel_sequence[queue_n ++] = i;
			vis[i] = 1;
			// printf("peel vertex %u\n", i);
		}
	}

	// 将不符合条件的点删除，更新度数，并继续检查是否符合条件
	for (ui i = 0; i < queue_n; i ++)
	{
		ui u = peel_sequence[i];
		for (ui j = 0; j < l; j ++)
		{
			for (ept k = pstart[n*j + u]; k < pstart[n*j + u + 1]; k ++)
			{
				if (!vis[edges[k]])
				{
					if ((degree[n*j + edges[k]] --) == threshold)
					{
						ui cnt = 0;
						for (ui j = 0; j < l; j ++)
						{
							if (degree[n*j + edges[k]] >= threshold) cnt++;
						}
						if (cnt < Sigma)
						{
							peel_sequence[queue_n ++] = edges[k];
							vis[edges[k]] = 1;
						}
					}
				}
			}
		}
	}
	printf("*** Degen delete vertexs number %u\n", queue_n);
	ui UB = n;
	ui multilayer_UB[l];
	if(queue_n == n)
	{
		UB = kplex.size();
		return UB;
	}
	
	// 寻找每层的上界UB，并寻找每层每个节点的核数
	for (ui i = 0; i < l; i ++)
	{
		memset(vis, 0, sizeof(char)*n);
		new_size = 0;
		for (ui j = 0; j < queue_n; j ++)
		{
			vis[peel_sequence[j]] = 1;
			multilayer_core[n*i + peel_sequence[j]] = 0;
		}
		for (ui j = 0; j < n; j ++)
		{
			if (!vis[j]) peel_sequence[queue_n + (new_size ++)] = j;
		}
		assert(queue_n + new_size == n);
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree + n*i);
		ui max_core = 0;
		UB = 0;
		for (ui j = 0; j < new_size; j ++)
		{
			ui u, key;
			heap->pop_min(u, key);
			if (key > max_core) max_core = key;
			multilayer_core[n*i + u] = max_core;
			peel_sequence[queue_n + j] = u;
			ui t_UB = multilayer_core[n*i + u] + K;
			if (new_size - j < t_UB) t_UB = new_size - j;
			if (t_UB > UB) UB = t_UB;
			for (ept k = pstart[n*i + u]; k < pstart[n*i + u + 1]; k ++) if (vis[edges[k]] == 0) {
				heap->decrement(edges[k], 1);
			}
		}
		multilayer_UB[i] = UB;
	}

	// 提取第Sigma大的上界为总上界，和Sigma大的核数为总核数
	for (ui i = 0; i < l - 1; i ++)
	{
		for (ui j = 0; j < l - i - 1; j ++)
		{
			if (multilayer_UB[j] < multilayer_UB[j + 1]) exchange_value(multilayer_UB[j], multilayer_UB[j + 1]);
			for (ui k = 0; k < n; k ++)
			{
				if (multilayer_core[n*j + k] < multilayer_core[n*(j+1) + k]) exchange_value(multilayer_core[n*j + k], multilayer_core[n*(j+1) + k]);
			}
		}
	}

	for (ui i = 0; i < n; i ++)
	{
		core[i] = multilayer_core[n*(Sigma - 1) + i];
	}

	heap->init(new_size, new_size - 1, peel_sequence + queue_n, core);
	for (ui j = 0; j < new_size; j ++)
	{
		ui u, key;
		heap->pop_min(u, key);
		peel_sequence[queue_n + j] = u;
	}

	UB = multilayer_UB[Sigma - 1];

	if(output) printf("*** UB: %u, Time: %s (microseconds)\n", UB, Utility::integer_to_string(t.elapsed()).c_str());

	return UB;
}

void Graph::shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *pstart, ui *edges, bool output) {
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] + K > kplex.size()) {
		rid[i] = cnt;
		if(in_mapping == NULL) out_mapping[cnt] = i;
		else out_mapping[cnt] = in_mapping[i];
		++ cnt;
	}

	ui t_n = cnt;

	if (cnt != n)
	{
		ept pos = 0;
		for (ui i = 0; i < l; i ++)
		{
			cnt = 0;
			for (ui j = 0; j < n; j ++) if (core[j] + K > kplex.size()) {
				ept t_start = pstart[n*i + j];
				pstart[t_n*i + cnt] = pos;
				for (ept k = t_start; k < pstart[n*i + j + 1]; k ++) if(core[edges[k]] + K > kplex.size()) {
					edges[pos ++] = rid[edges[k]];
				}
				++ cnt;
			}
			pstart[t_n*i + cnt] = pos;
		}

		// printf("%u %u %u %u %u\n", n, cnt, t_n, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n-cnt-1]] == 0||core[peel_sequence[n-cnt-1]] + K <= kplex.size());
		assert(cnt == 0||core[peel_sequence[n-cnt]] + K > kplex.size());
		for(ui i = 0;i < cnt;i ++) {
			peel_sequence[i] = rid[peel_sequence[n-cnt+i]];
			//core[i] = core[out_mapping[i]];
		}

		n = cnt;
		m = pos;
	}

	if(output) printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());
}

// orient graph and triangle counting
void Graph::oriented_triangle_counting(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj) {
	ui *rid = adj;
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for (ui j = 0; j < l; j ++)
	{
		for (ui i = 0; i < n; i ++)
		{
			ept &end = pend[n*j + i] = pstart[n*j + i];
			for (ept k = pstart[n*j + i]; k < pstart[n*j + i + 1]; k ++) if (rid[edges[k]] > rid[i]) edges[end ++] = edges[k];
		}
	}

#ifndef NDEBUG
	long long sum = 0;
	for(int i = 0;i < n*l;i ++) sum += pend[i] - pstart[i];
	// printf("%lld %lld\n", sum, m);
	assert(sum*2 == m);
#endif

	memset(adj, 0, sizeof(ui)*n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui)*m);
	for (ui i = 0; i < l; i ++)
	{
		for (ui u = 0; u < n; u ++)
		{
			for (ept j = pstart[n*i + u]; j < pend[n*i + u]; j ++) adj[edges[j]] = j + 1;
			for (ept j = pstart[n*i + u]; j < pend[n*i + u]; j ++)
			{
				ui v = edges[j];
				for (ept k = pstart[n*i + v]; k < pend[n*i + v]; k ++) if (adj[edges[k]]) {
					++ tri_cnt[j];
					++ tri_cnt[k];
					++ tri_cnt[adj[edges[k]] - 1];
					++ cnt;
				}
			}
			for (ept j = pstart[n*i + u]; j < pend[n*i + u]; j ++) adj[edges[j]] = 0;
		}
	}

#ifndef NDEBUG
	printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf) {
	for(ui i = 0;i < n*l;i ++) pend2[i] = pend[i];
	ept pos = 0;
	for (ui t = 0; t < l; t ++)
	{
		for (ui i = 0; i < n; i ++)
		{
			for (ept j = pstart[n*t + i]; j < pend[n*t + i]; j ++)
			{
				tri_cnt[pos>>1] = edgelist_pointer[j]; edge_list[pos++] = n*t + i; edge_list[pos++] = n*t + edges[j];

				ept &k = pend2[n*t + edges[j]];
				edgelist_pointer[k] = edgelist_pointer[j] = (pos>>1) - 1;
				edges[k ++] = i;
			}
		}
	}
	
#ifndef NDEBUG
	for(ui i = 0;i < n*l;i ++) assert(pend2[i] == pstart[i+1]);
#endif

	for(ui i = 0;i < n*l;i ++) {
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui t = 0; t < l; t ++)
	{
		for (ui i = 0; i < n; i ++)
		{
			for (ept j = pend2[n*t + i]; j < pstart[n*t + i + 1]; j ++)
			{
				ept &k = pend[n*t + edges[j]];
				edgelist_pointer[k] = edgelist_pointer[j];
				edges[k ++] = i;
			}
		}
	}

	ept *ids = pend2;
	for(ui i = 0;i < n*l;i ++) {
		if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while(j < pend[i]&&k < pstart[i+1]) {
			if(edges[j] < edges[k]) {
				ids[pos] = edges[j];
				buf[pos ++] = edgelist_pointer[j ++];
			}
			else {
				ids[pos] = edges[k];
				buf[pos ++] = edgelist_pointer[k ++];
			}
		}
		while(j < pend[i]) {
			ids[pos] = edges[j];
			buf[pos ++] = edgelist_pointer[j ++];
		}
		while(k < pstart[i+1]) {
			ids[pos] = edges[k];
			buf[pos ++] = edgelist_pointer[k ++];
		}
		for(ept j = 0;j < pos;j ++) {
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
#ifndef NDEBUG
		assert(pos == pstart[i+1] - pstart[i]);
#endif
	}
// #ifndef NDEBUG
// 	for (ui i = 0; i <n*l; i ++){
// 		if (pstart[i+1] - pstart[i] > 1) for (ui j = pstart[i]; j < pstart[i+1] - 1; j ++){
// 			assert(edges[j] < edges[j+1]);
// 		}
// 	}
// 	printf("check order!\n");
// #endif
}

char Graph::find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges) {
	if(b >= e) return 0;

	while(b+1 < e) {
		idx = b + (e-b)/2;
		if(edges[idx] > w) e = idx;
		else b = idx;
	}

	if(edges[b] == w) {
		idx = edgelist_pointer[b];
		if(!deleted[idx]) return 1;
	}

	return 0;
}

// return the number of peeled edges
ept Graph::peeling(ui critical_vertex, ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, char *deleted_vertex, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists) {
	ept Qe_n = 0;
	if(initialize_Qe) {
		ept active_edgelist_newn = 0;
		for(ept j = 0;j < active_edgelist_n;j ++) if(!deleted[active_edgelist[j]]) {
			if(tri_cnt[active_edgelist[j]] < t_threshold) Qe[Qe_n++] = active_edgelist[j];
			else active_edgelist[active_edgelist_newn ++] = active_edgelist[j];
		}
		active_edgelist_n = active_edgelist_newn;
	}
	ept deleted_edges_n = 0;
	ui Qv_idx = 0;
	while(Qv_idx < Qv_n || Qe_n) {
		if(Qe_n == 0) {
			//printf("hit\n");
			ui u = Qv[Qv_idx ++]; // delete u from the graph due to have a degree < d_threshold
			
			for (ui t = 0; t < l; t ++)
			{
				ept u_n = pstart[n*t + u];
				for (ept k = pstart[n*t + u]; k < pend[n*t + u]; k ++) if (!deleted[edgelist_pointer[k]]) {
					edges[u_n] = edges[k]; edgelist_pointer[u_n++] = edgelist_pointer[k];
					exists[edges[k]] = 1;
				}
				pend[n*t + u] = u_n;
				for (ept k = pstart[n*t + u]; k < pend[n*t + u]; k ++) deleted[edgelist_pointer[k]] = 1;
				deleted_edges_n += pend[n*t + u] - pstart[n*t + u];
				degree[n*t + u] = 0;

				for (ept k = pstart[n*t + u]; k < pend[n*t + u]; k ++)
				{
					ui v = edges[k];
					ept v_n = pstart[n*t + v];
					for (ept x = pstart[n*t + v]; x < pend[n*t + v]; x ++) if (!deleted[edgelist_pointer[x]]) {
						edges[v_n] = edges[x]; edgelist_pointer[v_n ++] = edgelist_pointer[x];
						if (edges[x] > v && exists[edges[x]]) {
							if( (tri_cnt[edgelist_pointer[x]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[x];
						}
					}
					pend[n*t + v] = v_n;

					if ( (degree[n*t + v]--) == d_threshold) if (!deleted_vertex[v]) {
						ui cnt = 0;
						for (ui i = 0; i < l; i ++)
						{
							if (degree[n*i + v] >= d_threshold) cnt++;
						}
						if (cnt < Sigma)
						{
							Qv[Qv_n++] = v;
							deleted_vertex[v] = 1;
							if (v == critical_vertex) {
								for(ept k = pstart[n*t + u]; k < pend[n*t + u]; k ++) exists[edges[k]] = 0;
								return 0;
							}
						}
					}
					if(linear_heap != NULL) linear_heap->decrement(v, 1);
				}
				for(ept k = pstart[n*t + u]; k < pend[n*t + u]; k ++) exists[edges[k]] = 0;
			}
		}
		for(ept j = 0;j < Qe_n;j ++) {
// #ifndef NDEBUG
// 			for (ui a = 0; a < m / 2; a ++) if (!deleted[a])
// 			{
// 				ui u = edge_list[a<<1] % n, v = edge_list[(a<<1)+1] % n, layer = edge_list[a<<1] / n;
// 				ui exi[n];
// 				ui cnt = 0;
// 				memset(exi, 0, sizeof(ui)*n);
// 				for (ui i = pstart[n*layer + u]; i < pend[n*layer + u]; i ++) if(!deleted[edgelist_pointer[i]]) exi[edges[i]] = 1;
// 				for (ui i = pstart[n*layer + v]; i < pend[n*layer + v]; i ++) if(!deleted[edgelist_pointer[i]]) if(exi[edges[i]]) cnt++;
// 				assert(cnt == tri_cnt[a]);
// 			}
// 			printf("A: Triangle number is correct!\n");
// #endif
			ept idx = Qe[j];
			ui u = edge_list[idx<<1] % n, v = edge_list[(idx<<1)+1] % n, layer = edge_list[idx<<1] / n;
			ui tri_n = tri_cnt[idx];
#ifndef NDEBUG
			assert(edge_list[idx<<1] / n == edge_list[(idx<<1)+1] / n);
#endif
			// printf("remove %u %u %u %u %lu\n", layer, u, v, tri_n, idx);
			deleted[idx] = 1;
			if( (degree[n*layer + u] --) == d_threshold) if (!deleted_vertex[u]){
				ui cnt = 0;
				for (ui i = 0; i < l; i ++)
				{
					if (degree[n*i + u] >= d_threshold) cnt++;
				}
				if (cnt < Sigma)
				{
					Qv[Qv_n++] = u;
					deleted_vertex[u] = 1;
					if (u == critical_vertex) return 0;
				}
			}
			if( (degree[n*layer + v] --) == d_threshold) if (!deleted_vertex[u]){
				ui cnt = 0;
				for (ui i = 0; i < l; i ++)
				{
					if (degree[n*i + v] >= d_threshold) cnt++;
				}
				if (cnt < Sigma)
				{
					Qv[Qv_n++] = v;
					deleted_vertex[u] = 1;
					if (v == critical_vertex) return 0;
				}
			}
			//printf("before\n");
			if(linear_heap != NULL) {
				linear_heap->decrement(u, 1);
				linear_heap->decrement(v, 1);
			}
			//printf("after\n");
			deleted_edges_n ++;
			
			if(degree[n*layer + u] < degree[n*layer + v]) swap(u,v);
			//printf("here\n");

			if(degree[n*layer + u] > degree[n*layer + v]*2) { // binary search
			//if(false) {
				ept v_n = pstart[n*layer + v], start = pstart[n*layer + u];
				for(ept k = pstart[n*layer + v];k < pend[n*layer + v];k ++) if(!deleted[edgelist_pointer[k]]) {
					edges[v_n] = edges[k]; edgelist_pointer[v_n++] = edgelist_pointer[k];
					if(tri_n && find(u, edges[k], start, pend[n*layer + u], deleted, idx, edgelist_pointer, edges)) {
						-- tri_n;
						if( (tri_cnt[idx]--) == t_threshold) Qe[Qe_n++] = idx;
						if( (tri_cnt[edgelist_pointer[k]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[k];
					}
				}
				pend[n*layer + v] = v_n;
// #ifndef NDEBUG
// 				for (ui a = 0; a < m / 2; a ++) if (!deleted[a])
// 				{
// 					ui u = edge_list[a<<1] % n, v = edge_list[(a<<1)+1] % n, layer = edge_list[a<<1] / n;
// 					ui exi[n];
// 					ui cnt = 0;
// 					memset(exi, 0, sizeof(ui)*n);
// 					for (ui i = pstart[n*layer + u]; i < pend[n*layer + u]; i ++) if(!deleted[edgelist_pointer[i]]) exi[edges[i]] = 1;
// 					for (ui i = pstart[n*layer + v]; i < pend[n*layer + v]; i ++) if(!deleted[edgelist_pointer[i]]) if(exi[edges[i]]) cnt++;
// 					// printf("%u %u %u %u %u %u\n",a, u, v, layer, cnt, tri_cnt[a] );
// 					assert(cnt == tri_cnt[a]);
// 				}
// 				printf("B: Triangle number is correct!\n");
// #endif
				assert(tri_n == 0);
			}
			else { // sorted_merge
				ept ii = pstart[n*layer + u], jj = pstart[n*layer + v];
				ept u_n = pstart[n*layer + u], v_n = pstart[n*layer + v];
				while(true) {
					while(ii < pend[n*layer + u]&&deleted[edgelist_pointer[ii]]) ++ ii;
					while(jj < pend[n*layer + v]&&deleted[edgelist_pointer[jj]]) ++ jj;
					if(ii >= pend[n*layer + u]||jj >= pend[n*layer + v]) break;
					if(edges[ii] == edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if( (tri_cnt[edgelist_pointer[ii]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[ii];
						if( (tri_cnt[edgelist_pointer[jj]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[jj];
						++ ii;
						++ jj;
					}
					else if(edges[ii] < edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ ii;
					}
					else {
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++ jj;
					}
				}
				while(ii < pend[n*layer + u]) {
					if(!deleted[edgelist_pointer[ii]]) {
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ ii;
				}
				while(jj < pend[n*layer + v]) {
					if(!deleted[edgelist_pointer[jj]]) {
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++ jj;
				}
				pend[n*layer + u] = u_n; pend[n*layer + v] = v_n;
// #ifndef NDEBUG
// 				for (ui a = 0; a < m / 2; a ++) if (!deleted[a])
// 				{
// 					ui u = edge_list[a<<1] % n, v = edge_list[(a<<1)+1] % n, layer = edge_list[a<<1] / n;
// 					ui exi[n];
// 					ui cnt = 0;
// 					memset(exi, 0, sizeof(ui)*n);
// 					for (ui i = pstart[n*layer + u]; i < pend[n*layer + u]; i ++) if(!deleted[edgelist_pointer[i]]) exi[edges[i]] = 1;
// 					for (ui i = pstart[n*layer + v]; i < pend[n*layer + v]; i ++) if(!deleted[edgelist_pointer[i]]) if(exi[edges[i]]) cnt++;
// 					// printf("%u %u %u %u %u %u\n", a, u, v, layer, cnt, tri_cnt[a]);
// 					assert(cnt == tri_cnt[a]);
// 				}
// 				printf("C: Triangle number is correct!\n");
// #endif
			}
			//printf("finish %u %u\n", u, v);
		}
		Qe_n = 0;
	}
#ifndef NDEBUG
	printf("*** Truss removed %s undirected edges\n", Utility::integer_to_string(deleted_edges_n).c_str());
#endif
	return deleted_edges_n;
}

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, vector<tuple<int,int,int> > &vtp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer) {
	ids_n = 0; vtp.clear();
	ids[ids_n++] = u; exists[u] = 1; rid[u] = 0;
	char t_exists[n];
	ui t_ids[n];
	ui t_ids_n;
	ui old_u = u;
	
	for (ui t = 0; t < l; t ++)
	{
		t_ids_n = 0;
		memset(t_exists, 0, sizeof(char)*n);
		u = old_u;
		t_ids[t_ids_n++] = u; t_exists[u] = 1;
		ui u_n = pstart[n*t + u];
		for (ept i = pstart[n*t + u]; i < pend[n*t + u]; i ++) if (!deleted[edgelist_pointer[i]]) {
			edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			if(t_exists[v]) continue;
			t_ids[t_ids_n++] = v; t_exists[v] = 1;
			if(exists[v]) continue;
			rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
		}
		pend[n*t + u] = u_n;
		ui old_size = t_ids_n;
		for (ui i = 1; i < old_size; i ++)
		{
			u = t_ids[i];
			u_n = pstart[n*t + u];
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if (!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				ui v = edges[j];
				if(t_exists[v]) continue;
				t_ids[t_ids_n++] = v; t_exists[v] = 1;
				if(exists[v]) continue;
				rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
			}
			pend[n*t + u] = u_n;
		}

		for (ui i = 0; i < old_size; i ++)
		{
			u = t_ids[i];
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if(edges[j] > u) {
				vtp.pb(mtp(t, rid[u], rid[edges[j]]));
			}
		}
		for (ui i = old_size; i < t_ids_n; i ++)
		{
			u = t_ids[i];
			u_n = pstart[n*t + u];
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if(!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				if(edges[j] > u && t_exists[edges[j]]) vtp.pb(mtp(t, rid[u], rid[edges[j]]));
			}
			pend[n*t + u] = u_n;
		}
	}
	for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<tuple<int,int,int> > &vtp, ui *Q, ui* degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer) {
	vtp.clear();
	ids_n = 0; ids[ids_n++] = u; exists[u] = 1; rid[u] = 0;
	char t_exists[n];
	ui t_ids[n];
	ui t_ids_n;
	ui old_u = u;

	for (ui t = 0; t < l; t ++)
	{
		t_ids_n = 0;
		memset(t_exists, 0, sizeof(char)*n);
		u = old_u;
		t_ids[t_ids_n++] = u; t_exists[u] = 1;
		ui u_n = pstart[n*t + u];
		for (ept i = pstart[n*t + u]; i < pend[n*t + u]; i ++) if(!deleted[edgelist_pointer[i]]) {
			edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			t_ids[t_ids_n++] = v; t_exists[v] = 2;
		}
		pend[n*t + u] = u_n;

		ui Q_n = 0;
		for (ui i = 1; i <t_ids_n; i ++)
		{
			u = t_ids[i];
			u_n = pstart[n*t + u];
			degree[u] = 0;
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if(!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (t_exists[edges[j]] == 2) ++ degree[u];
			}
			pend[n*t + u] = u_n;
			if (degree[u] + 2*K <= kplex.size()) Q[Q_n++] = u;
		}
		for (ui i = 0; i < Q_n; i ++)
		{
			u = Q[i];
			t_exists[u] = 10;
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if (t_exists[edges[j]] == 2) {
				if ((degree[edges[j]]--) + 2*K == kplex.size()+1) {
					assert(Q_n < m/2);
					Q[Q_n++] = edges[j];
				}
			}
		}
		assert(Q_n <= t_ids_n);
		if (t_ids_n - 1 - Q_n + K <= kplex.size()) continue;
		ui nr_size = t_ids_n;
		for (ui i = 1; i < nr_size; i ++) if (t_exists[t_ids[i]] == 2) {
			u = t_ids[i];
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++)
			{
				// printf("A %u\n", t_exists[edges[j]]);
				if (!t_exists[edges[j]]) {
					t_ids[t_ids_n++] = edges[j];
					t_exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (t_exists[edges[j]] == 3) ++ degree[edges[j]];
			}
		}
		ui new_size = 1;
		for (ui i = 1; i < nr_size; i ++)
		{
			if (t_exists[t_ids[i]] == 10) t_exists[t_ids[i]] = 0;
			else t_ids[new_size++] = t_ids[i];
		}
		assert(new_size + Q_n == nr_size);
		ui old_nr_size = nr_size;
		nr_size = new_size;
		for (ui i = old_nr_size; i <t_ids_n; i ++) {
			if(degree[t_ids[i]] + 2*K <= kplex.size() + 2) t_exists[t_ids[i]] = 0;
			else t_ids[new_size++] = t_ids[i];
		}
		t_ids_n = new_size;
#ifndef NDEBUG
		assert(t_exists[t_ids[0]] == 1);
		for(ui i = 1;i < nr_size;i ++) assert(t_exists[t_ids[i]] == 2);
		for(ui i = nr_size;i < t_ids_n;i ++) assert(t_exists[t_ids[i]] == 3);
#endif
		for (ui i = 1; i < t_ids_n; i ++)
		{
			assert(t_exists[t_ids[i]]);
			ui v = t_ids[i];
			if (!exists[v])
			{
				rid[v] = ids_n; ids[ids_n ++] = v; exists[v] = 1;
			}
		}
		for (ui i = 0; i < nr_size; i ++)
		{
			u = t_ids[i];
			for (ept j = pstart[n*t + u]; j < pend[n*t + u]; j ++) if(t_exists[edges[j]] && edges[j] > u) {
				assert(!deleted[edgelist_pointer[j]]);
				vtp.pb(mtp(t, rid[u], rid[edges[j]]));
			}
		}
		for (ui i = nr_size; i <t_ids_n; i ++)
		{
			u = t_ids[i];
			u_n = pstart[n*t + u];
			for (ept j = pstart[n*t +u]; j < pend[n*t + u]; j ++) if(!deleted[edgelist_pointer[j]]) {
				edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && t_exists[edges[j]]) vtp.pb(mtp(t, rid[u], rid[edges[j]]));
			}
			pend[n*t + u] = u_n;
		}
	}
	for(ui i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(exists[i] == 0);
#endif
}

void Graph::exchange_value(ui &a, ui &b)
{
	ui temp;
	temp = a;
	a = b;
	b = temp;
}

