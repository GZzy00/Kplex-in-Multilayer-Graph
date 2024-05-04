#ifndef _KPLEX_BB_MATRIX_
#define _KPLEX_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"

// #define _SECOND_ORDER_PRUNING_

//#define set_bit(array, pos) (((array)[pos]) = 1)
//#define reverse_bit(array, pos) (((array)[pos]) = 1- ((array)[pos]))
//#define test_bit(array, pos) ((array)[pos])

class KPLEX_BB_MATRIX
{
private:
  long long n;
	ui l;

  char *matrix;
  long long matrix_size;

#ifdef _SECOND_ORDER_PRUNING_
  ui *cn;
  std::queue<std::tuple<ui,ui,ui> > Qe;
  std::vector<std::tuple<ui,ui,ui> > removed_edges;
  long long removed_edges_n;
#endif

  ui *degree;
  ui *degree_in_S;

  char *efficient_layer;

  ui K;
	ui Sigma;
  ui *best_solution;
  ui best_solution_size;

  ui *neighbors;
  ui *nonneighbors;

  ui *SR; // union of S and R, where S is at the front
  ui *SR_rid; // reverse ID for SR
  std::queue<ui> Qv;
  ui *level_id;

  std::vector<ui> must_include_vertices;

public:
  KPLEX_BB_MATRIX() 
	{
    n = 0;
		l = 0;
    matrix = NULL;
    matrix_size = 0;

#ifdef _SECOND_ORDER_PRUNING_
    cn = NULL;
    removed_edges_n = 0;
#endif

    degree = degree_in_S = NULL;

    efficient_layer = NULL;
        
    best_solution = NULL;
    Sigma = K = best_solution_size = 0;

    neighbors = nonneighbors = NULL;

    SR = SR_rid = NULL;
    level_id = NULL;
  }

  ~KPLEX_BB_MATRIX()
	{
    if(matrix != NULL) {
      delete[] matrix;
      matrix = NULL;
    }
#ifdef _SECOND_ORDER_PRUNING_
    if(cn != NULL) {
      delete[] cn;
      cn = NULL;
    }
#endif
    if(degree != NULL) {
      delete[] degree;
      degree = NULL;
    }
    if(degree_in_S != NULL) {
      delete[] degree_in_S;
      degree_in_S = NULL;
    }
    if(efficient_layer != NULL) {
      delete[] efficient_layer;
      efficient_layer = NULL;
    }
    if(best_solution != NULL) {
      delete[] best_solution;
      best_solution = NULL;
    }
    if(SR != NULL) {
      delete[] SR;
      SR = NULL;
    }
    if(SR_rid != NULL) {
      delete[] SR_rid;
      SR_rid = NULL;
    }
    if(neighbors != NULL) {
      delete[] neighbors;
      neighbors = NULL;
    }
    if(nonneighbors != NULL) {
      delete[] nonneighbors;
      nonneighbors = NULL;
    }
    if(level_id != NULL) {
    delete[] level_id;
    level_id = NULL;
    }
  }

  void allocateMemory(ui n, ui m, ui l)
	{
	  if(n <= 0||m <= 0) return ;

    matrix = new char[m*2];
#ifdef _SECOND_ORDER_PRUNING_
    cn = new ui[m*2];
#endif
    matrix_size = m*2;
        
    degree = new ui[n*l];
    degree_in_S = new ui[n*l];
    efficient_layer = new char[l];
    best_solution = new ui[n];
    SR = new ui[n];
    SR_rid = new ui[n];
    neighbors = new ui[n*l];
    nonneighbors = new ui[n*l];
    level_id = new ui[n];
  }

  void load_graph(ui _n, ui _l, const std::vector<std::tuple<int,int,int> > &vtp) 
	{
    n = _n;
		l = _l;
    if(((long long)n)*n*l > matrix_size) 
		{
      do {
        matrix_size *= 2;
      } while(((long long)n)*n*l > matrix_size);
      delete[] matrix; matrix = new char[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
      delete[] cn; cn = new ui[matrix_size];
#endif
    }

#ifdef _SECOND_ORDER_PRUNING_
    memset(cn, 0, sizeof(ui)*((long long)n)*n*l);
#endif
    memset(matrix, 0, sizeof(char)*((long long)n)*n*l);
    for(ui i = 0; i < n*l; i++) degree[i] = 0;
    for(ui i = 0;i < vtp.size();i ++)
		{
			ui a = std::get<0>(vtp[i]), b = std::get<1>(vtp[i]), c = std::get<2>(vtp[i]);
      assert(a >= 0 && a < l && b >= 0 && b < n && c >= 0 && c < n);
      degree[n*a + b] ++;
      degree[n*a + c] ++;
      matrix[n*n*a + n*b + c] = matrix[n*n*a +n*c + b] = 1;
      // printf("edge: %u %u %u\n", a, b, c);
    }

#ifndef NDEBUG
    printf("load graph of size l = %u n=%lld, m=%lu\n", l, n, vtp.size());
    //for(ui i = 0;i < vp.size();i ++) printf("%d %d\n", vp[i].first, vp[i].second);
#endif
  }

  void kPlex(ui K_, ui _Sigma, std::vector<ui> &kplex, bool must_include_0)
	{
    K = K_;
		Sigma = _Sigma;
    if(K == 1) {
      printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
      return ;
    }
    best_solution_size = kplex.size();
    ui R_end;
    initialization(R_end, must_include_0);
    if(R_end) BB_search(0, R_end, 1, must_include_0, 0);
    if(best_solution_size > kplex.size()) {
      kplex.clear();
      for(int i = 0;i < best_solution_size;i ++) kplex.push_back(best_solution[i]);
    }
  }

private:
  void initialization(ui &R_end, bool must_include_0)
	{
    // the following computes a degeneracy ordering and a heuristic solution
  	ui *core = nonneighbors;
  	ui *vis = SR;
		ui multilayer_core[n*l];
		ui multilayer_UB[l];
		memset(multilayer_core, 0, sizeof(ui)*n*l);
		memset(multilayer_UB, 0, sizeof(ui)*n);
		for (ui t = 0; t < l; t ++)
		{
			ui max_core = 0, UB = 0;
			memset(vis, 0, sizeof(ui)*n);
			for (ui i = 0; i < n; i ++)
			{
				ui u, min_degree = n;
				for(ui j = 0; j < n; j ++) if(!vis[j] && degree[n*t + j] < min_degree) {
					u = j;
					min_degree = degree[n*t + j];
				}
				if (min_degree > max_core) max_core = min_degree;
				multilayer_core[n*t + u] = max_core;
				vis[u] = 1;
				ui t_UB = multilayer_core[n*t + u] + K;
				if (n - i < t_UB) t_UB = n - i;
				if (t_UB > UB) UB = t_UB;
				
				for (ui j = 0; j < n; j ++) if(!vis[j]&&matrix[n*n*t + n*u + j]) --degree[n*t + j];
			}
			multilayer_UB[t] = UB;
		}

		for (ui i = 0; i < l - 1; i ++) {
			for (ui j = 0; j < l - i - 1; j ++) {
				if (multilayer_UB[j] < multilayer_UB[j + 1]) std::swap(multilayer_UB[j], multilayer_UB[j + 1]);
				for (ui k = 0; k < n; k ++) {
					if (multilayer_core[n*j + k] < multilayer_core[n*(j+1) + k]) std::swap(multilayer_core[n*j + k], multilayer_core[n*(j+1) + k]);
				}
			}
		}

		for (ui i = 0; i < n; i ++) {
			core[i] = multilayer_core[n*(Sigma - 1) + i];
		}

    memset(degree_in_S, 0, sizeof(ui)*n*l);
    R_end = 0;
    for(ui i = 0;i < n;i ++) SR_rid[i] = n;
    for(ui i = 0;i < n;i ++) if(core[i] + K > best_solution_size) {
      SR[R_end] = i; SR_rid[i] = R_end;
      ++ R_end;
    }
    if(must_include_0&&SR_rid[0] == n) {
      R_end = 0;
      return ;
    }
		for (ui t = 0; t < l; t ++) {
			for (ui i = 0; i < R_end; i ++) {
				ui u = SR[i];
				degree[n*t + u] = 0;
				for (ui j = 0; j < R_end; j ++) if(matrix[n*n*t + n*u + SR[j]]) ++ degree[n*t + u];
			}
		}

    ui cnt = 0;
    for (ui t = 0; t < l; t ++) if(degree[n*t] + K > best_solution_size) cnt ++;
    // Core 符合要求，但是可能邻居点被删了，导致初始点0不符合要求
    if (cnt < Sigma) {
      R_end = 0;
      return ;
    }

    memset(level_id, 0, sizeof(ui)*n);
    for(ui i = 0;i < R_end;i ++) level_id[SR[i]] = n;

    assert(Qv.empty());

// #ifdef _SECOND_ORDER_PRUNING_
// 		for (ui t = 0; t < l; t ++)
// 		{
// 			for (ui i = 0; i < R_end; i ++)
// 			{
// 				ui neighbors_n = 0;
// 				char *t_matrix = matrix + n*n*t + n*SR[i];
// 				for (ui j = 0; j < R_end; j ++) if (t_matrix[SR[j]]) neighbors[neighbors_n ++] = SR[j];
// 				for (ui j = 0; j < neighbors_n; j ++) for(ui k = j+1; k <neighbors_n; k ++)
// 				{
// 					++ cn[n*n*t + neighbors[j]*n + neighbors[k]];
// 					++ cn[n*n*t + neighbors[k]*n + neighbors[j]];
// 				}
// 			}
// 		}
//     while(!Qe.empty()) Qe.pop();
// 		for (ui t = 0; t < l; t ++)
// 		{
// 			for (ui i = 0; i < R_end; i ++) for (ui j = i+1; j < R_end; j ++)
// 			{
// 				if(matrix[n*n*t + SR[i]*n + SR[j]] && upper_bound_based_prune(0, t, SR[i], SR[j])) {
// 					Qe.push(std::mtp(t, SR[i], SR[j]));
// 				}
// 			}
// 		}
//     removed_edges_n = 0;
// #endif
    for (ui t = 0; t < l; t ++) efficient_layer[t] = 1;

    must_include_vertices.resize(n);
    if(remove_vertices_and_edges_with_prune(0, R_end, 0)) R_end = 0;
  }

	bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level)
	{
#ifdef _SECOND_ORDER_PRUNING_
    while(!Qv.empty()||!Qe.empty()) {
#else
    while(!Qv.empty()) {
#endif
      while(!Qv.empty())
			{
				ui u = Qv.front(); Qv.pop(); // remove u
				assert(SR[SR_rid[u]] == u);
				assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
				-- R_end;
				swap_pos(SR_rid[u], R_end);

				bool terminate = false;
				ui neighbors_n[l];
				char *t_matrix;
        for (ui t = 0; t < l; t ++) {
          t_matrix = matrix + n*n*t + n*u;
          neighbors_n[t] = 0;
          for (ui i = 0; i < R_end; i ++) if(t_matrix[SR[i]]) {
            ui w = SR[i];
            neighbors[n*t + neighbors_n[t]++] = w;
            -- degree[n*t + w];
            if (degree[n*t + w] + K == best_solution_size) if(efficient_layer[t]){
              ui cnt = 0;
              if (i < S_end) {
                for (ui j = 0; j < l; j ++) if(efficient_layer[j]) if(degree[n*j + w] + K > best_solution_size) {
                  cnt ++;
                }
                if (cnt < Sigma) terminate = true;
              } else {
                for (ui j = 0; j < l; j ++) if(efficient_layer[j]) if(degree[n*j + w] + K > best_solution_size && degree_in_S[n*j + w] + K > S_end) {
                  cnt ++;
                }
                if (cnt < Sigma) if (level_id[w] > level) {
                  level_id[w] = level;
                  Qv.push(w);
                }
              }
            }
          }
        }
				// UB1
        if (terminate) {
          for (ui t = 0; t < l; t ++) {
            for (ui i = 0; i < neighbors_n[t]; i++)
              ++degree[n*t + neighbors[n*t + i]];
          }
          level_id[u] = n;
          ++ R_end;
          return true;
        }
      }
    }

    return false;
  }

  void BB_search(ui S_end, ui R_end, ui level, bool choose_zero, ui must_include_vertices_n)
  {
    assert(!choose_zero||!must_include_vertices_n);
#ifndef NDEBUG
    for (ui t = 0; t < l; t ++)
    {
      efficient_layer[t] = 1;
      for (ui i = 0; i < S_end; i ++) if(degree_in_S[n*t + SR[i]] + K < S_end || degree[n*t + SR[i]] + K <= best_solution_size) {
        efficient_layer[t] = 0;
        break;
      }
    }
#endif
    if(S_end > best_solution_size)
    { // find a larger solution
      best_solution_size = S_end;
      for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
      ui cnt = 0;
      bool is_kplex = true;
      for(ui t = 0; t < l; t ++)
      {
        is_kplex = true;
        for (ui i = 0; i < best_solution_size; i ++) if(degree_in_S[n*t + best_solution[i]] + K < best_solution_size) is_kplex = false;
        if (is_kplex) cnt ++;
      }
      assert(cnt >= Sigma);
      printf("Find a k-plex of size: %u\n", best_solution_size);
#endif
    }
    if(R_end <= best_solution_size) return ;

#ifndef NDEBUG
    ui cnt_promise = 0;
    bool have_promise = true;
    for (ui t = 0; t < l; t ++)
    {
      have_promise = true;
      for (ui i = 0; i < S_end; i ++)
      {
        if (degree_in_S[n*t + SR[i]] + K < S_end) {
          have_promise = false;
          break;
        } else if (degree[n*t + SR[i]] + K <= best_solution_size) {
          have_promise = false;
          break;
        }
      }
      if (have_promise) cnt_promise ++;
    }
    assert(cnt_promise >= Sigma);
#endif
    bool is_kplex = true;
    ui cnt = 0;
    for (ui t = 0; t < l; t ++)
    {
      is_kplex = true;
      for (ui i = 0; i < R_end; i ++) if(degree[n*t + SR[i]] + K < R_end) is_kplex = false;
      if (is_kplex) cnt ++;
    }
    if(cnt >= Sigma) {
      best_solution_size = R_end;
      for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
      printf("Greedy find a k-plex of size: %u\n", R_end);
#endif
      return ;
    }
    else if(R_end == best_solution_size + 1) return ;

#ifndef NDEBUG
    for(ui t = 0; t < l; t ++) for(ui i = 0;i < R_end;i ++) {
      ui d1 = 0, d2 = 0;
      for(ui j = 0;j < S_end;j ++) if(matrix[n*n*t + n*SR[i] + SR[j]]) ++ d1;
      d2 = d1;
      for(ui j = S_end;j < R_end;j ++) if(matrix[n*n*t + n*SR[i] + SR[j]]) ++ d2;
      assert(d1 == degree_in_S[n*t + SR[i]]);
      assert(d2 == degree[n*t + SR[i]]);
    }
#endif
    bool must_include = false;
    ui u = n; // u is the branching vertex
    if(choose_zero) {
      assert(S_end == 0&&must_include_vertices_n == 0);
      if(SR_rid[0] >= R_end) return ;
      u = 0;
      must_include = true;
    }
    else if(must_include_vertices_n) {
      must_include = true;
      u = must_include_vertices[-- must_include_vertices_n];
      assert(SR_rid[u] >= S_end);
      if(SR_rid[u] >= R_end) return ;
    }
    else {
      // for(ui i = S_end;i < R_end;i ++) if(degree[SR[i]] + 2 >= R_end) {
      //   u = SR[i];
      //   must_include = true;
      //   break;
      // }
      if(u == n) {
        // u = choose_branch_vertex(S_end, R_end);
        u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end);
      }
    }
#ifndef NDEBUG
    cnt = 0;
    for (ui t = 0; t < l; t ++) if(efficient_layer[t]) if(degree[n*t + u] + K > best_solution_size && degree_in_S[n*t + u] + K > S_end) cnt ++;
    assert(cnt >= Sigma);
#endif
        // the first branch includes u into S
    assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
    swap_pos(S_end, SR_rid[u]);
    ++ S_end;

    ui pre_best_solution_size = best_solution_size, old_R_end = R_end;
    ui old_removed_edges_n = 0;
#ifdef  _SECOND_ORDER_PRUNING_
    old_removed_edges_n = removed_edges_n;
#endif
    if(!move_u_to_S_with_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1, false, must_include_vertices_n);
    restore_SR_and_edges(S_end, R_end, old_R_end, level, old_removed_edges_n);

    if(must_include) {
      move_u_to_R_wo_prune(S_end, R_end, level);
      return ;
    }
    // the second branch exclude u from S
    assert(Qv.empty());
    bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
    if(!pruned)
    {
      collect_removable_vertices_and_edges(S_end, R_end, level);
      if(!remove_vertices_and_edges_with_prune(S_end, R_end, level))
      {
        BB_search(S_end, R_end, level+1, false, must_include_vertices_n);
      }
    }
    restore_SR_and_edges(S_end, R_end, old_R_end, level, old_removed_edges_n);
  }

  bool move_u_to_S_with_prune(ui S_end, ui &R_end, ui level)
  {
    assert(S_end > 0);
    ui u = SR[S_end-1];
    char *t_matrix;
    ui neighbors_n[l], nonneighbors_n[l];
    ui deleted[n*l];
    memset(deleted, 0, sizeof(ui)*n*l);
    for (ui t = 0; t < l; t ++)
    {
      t_matrix = matrix + n*n*t + n*u;
      neighbors_n[t] = 0;
      nonneighbors_n[t] = 0;
      for (ui i = 0; i < R_end; i ++) if(i != S_end-1) {
        if(t_matrix[SR[i]]) neighbors[n*t + neighbors_n[t]++] = SR[i];
        else nonneighbors[n*t + nonneighbors_n[t]++] = SR[i];
      }

      for(ui i = 0; i < neighbors_n[t]; i ++) ++ degree_in_S[n*t + neighbors[n*t + i]];
    }
    assert(Qv.empty());
    for (ui t = 0; t < l; t ++)
    {
      efficient_layer[t] = 1;
      for (ui i = 0; i < S_end; i ++) if(degree_in_S[n*t + SR[i]] + K < S_end || degree[n*t + SR[i]] + K <= best_solution_size) {
        efficient_layer[t] = 0;
        break;
      }
    }
    for (ui t = 0; t < l; t ++) if (efficient_layer[t]) {
      // 删除导致新增点u不符合条件的点
      if (degree_in_S[n*t + u] + K == S_end) {
        ui i = 0;
        while(i < nonneighbors_n[t]&&SR_rid[nonneighbors[n*t + i]] < S_end) ++ i;
        for(;i < nonneighbors_n[t];i ++) { // remove non-neighbors from R
          deleted[n*t + nonneighbors[n*t + i]] = 1;
        }
      }

      // 删除本身不符合条件的候选集R中的点
      for (ui i = S_end; i < R_end; i ++) if (degree_in_S[n*t + SR[i]] + K <= S_end || degree[n*t + SR[i]] + K <= best_solution_size) {
        deleted[n*t + SR[i]] = 1;
      }

      // S中点处于临界，删除其非邻居点
      for(ui i = 0;i < S_end-1;i ++) if(degree_in_S[n*t + SR[i]] + K == S_end) {
        char *tt_matrix = matrix +n*n*t + n*SR[i];
        for(ui j = S_end;j < R_end;j ++) if(!tt_matrix[SR[j]]) {
          deleted[n*t + SR[j]] = 1;
        }
      }
    }
    for (ui i = S_end; i < R_end; i ++)
    {
      ui cnt = 0;
      for (ui t = 0; t < l; t ++) if (efficient_layer[t]) {
        if(!deleted[n*t + SR[i]]) cnt ++;
      }
      if (cnt < Sigma) {
        assert(level_id[SR[i]] > level);
        level_id[SR[i]] = level;
        Qv.push(SR[i]);
      }
    }
    return remove_vertices_and_edges_with_prune(S_end, R_end, level);
  }

  void restore_SR_and_edges(ui S_end, ui &R_end, ui old_R_end, ui level, ui old_removed_edges_n)
  {
    while(!Qv.empty()) {
      ui u = Qv.front(); Qv.pop();
      assert(level_id[u] == level&&SR_rid[u] < R_end);
      level_id[u] = n;
    }
    while(R_end < old_R_end) { // insert u back into R
      ui u = SR[R_end];
      assert(level_id[u] == level&&SR_rid[u] == R_end);
      level_id[u] = n;

      char *t_matrix;
      for (ui t = 0; t < l; t ++)
      {
        t_matrix = matrix + n*n*t + n*u;
        for (ui i = 0; i < R_end; i ++) if (t_matrix[SR[i]]) {
          ui w = SR[i];
          ++ degree[n*t + w];
        }
      }
      ++ R_end;
    }
  }

  void move_u_to_R_wo_prune(ui &S_end, ui &R_end, ui level)
  {
    assert(S_end);
    ui u = SR[-- S_end];
    ui neighbors_n = 0;
    char *t_matrix;
    for (ui t = 0; t < l; t ++) {
      t_matrix = matrix + n*n*t + n*u;
      for (ui i = 0; i < R_end; i ++) if (t_matrix[SR[i]]) -- degree_in_S[n*t + SR[i]];
    }
  }

  bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
  {
    assert(S_end);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool ret = false;
    ui neighbors_n[l];
    memset(neighbors_n, 0, sizeof(ui)*l);
    char *t_matrix;

    for (ui t = 0; t < l;t ++){
      t_matrix = matrix + n*n*t + n*u;
      for (ui i = 0; i < R_end; i ++) if (t_matrix[SR[i]]) neighbors[n*t + neighbors_n[t]++] = SR[i];
      for (ui i = 0; i < neighbors_n[t]; i ++) {
        -- degree_in_S[n*t +neighbors[n*t + i]];
        -- degree[n*t + neighbors[n*t + i]];
      }
    }

    for (ui t = 0; t < l; t ++)
    {
      efficient_layer[t] = 1;
      for (ui i = 0; i < S_end; i ++) if(degree_in_S[n*t + SR[i]] + K < S_end || degree[n*t + SR[i]] + K <= best_solution_size) {
        efficient_layer[t] = 0;
        break;
      }
    }

    ui cnt = 0;
    for (ui t = 0; t < l; t ++) if (efficient_layer[t]) cnt ++;
    if (cnt < Sigma) ret = true;

    if(ret) return true;
		return false;
	}

  void collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level)
  {
    ui deleted[n*l];
    memset(deleted, 0, sizeof(ui)*n*l);

    for (ui t = 0; t < l; t ++) if (efficient_layer[t]) {
      // 删除本身不符合条件的候选集R中的点
      for (ui i = S_end; i < R_end; i ++) if (degree_in_S[n*t + SR[i]] + K <= S_end || degree[n*t + SR[i]] + K <= best_solution_size) {
        deleted[n*t + SR[i]] = 1;
      }

      // S中点处于临界，删除其非邻居点
      for(ui i = 0;i < S_end;i ++) if(degree_in_S[n*t + SR[i]] + K == S_end) {
        char *tt_matrix = matrix +n*n*t + n*SR[i];
        for(ui j = S_end;j < R_end;j ++) if(!tt_matrix[SR[j]]) {
          deleted[n*t + SR[j]] = 1;
        }
      }
    }
    for (ui i = S_end; i < R_end; i ++)
    {
      ui cnt = 0;
      for (ui t = 0; t < l; t ++) if (efficient_layer[t]) {
        if(!deleted[n*t + SR[i]]) cnt ++;
      }
      if (cnt < Sigma) {
        assert(level_id[SR[i]] > level);
        level_id[SR[i]] = level;
        Qv.push(SR[i]);
      }
    }
  }

  ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end)
  {
    ui u = n, min_degree_in_S = n*l;
    for(ui i = S_end;i < R_end;i ++)
    {
    	ui v = SR[i];
      ui total_degree_in_S = 0;
      for (ui t = 0; t < l; t ++) total_degree_in_S += degree_in_S[n*t + v];
      if (total_degree_in_S < min_degree_in_S) {
        u = v;
        min_degree_in_S = total_degree_in_S;
      }
    }
    assert(u != n);
    return u;
  }


#ifdef _SECOND_ORDER_PRUNING_
  bool upper_bound_based_prune(ui S_end, ui t,ui u, ui v)
	{
    // ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
    ui ub = 2*K + degree_in_S[n*t + u] - S_end + degree_in_S[n*t + v] + cn[n*n*t + u*n + v];
    if(SR_rid[u] >= S_end) {
    	-- ub; // S_end ++
    	if(matrix[n*n*t+n*u+v]) ++ ub; // degree_in_S[v] ++
    }
    if(SR_rid[v] >= S_end) {
    	-- ub;
    	if(matrix[n*n*t+n*v+u]) ++ ub;
    }
    return ub <= best_solution_size;
  }
#endif

#ifdef _SECOND_ORDER_PRUNING_
  bool total_upper_bound_based_prune(ui S_end, ui u, ui v)
	{
    ui cnt = 0;
    for (ui t = 0; t < l; t ++)
    {
      if(!upper_bound_based_prune(S_end, t, u, v)) cnt++;
    }
    // ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
    return cnt < Sigma;
  }
#endif

  void swap_pos(ui i, ui j)
	{
    std::swap(SR[i], SR[j]);
    SR_rid[SR[i]] = i;
    SR_rid[SR[j]] = j;
  }
};

#endif
